/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <sys/times.h>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/ArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/Event.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactEventHandler.h>

using namespace Ravelin;
using namespace Moby;
using std::list;
using boost::shared_ptr;
using std::vector;
using std::map;
using std::endl;
using std::cerr;
using std::pair;
using std::min_element;
using boost::dynamic_pointer_cast;

/// Solves the quadratic program (potentially solves two QPs, actually)
void ImpactEventHandler::solve_qp(const VectorNd& zf, EventProblemData& q, double poisson_eps, double max_time)
{
  const double TOL = poisson_eps;

  // set z to frictionless solution to start
  const unsigned N_TOTAL = q.N_VARS + q.N_CONTACTS + q.N_LIMITS + q.N_K_TOTAL + 1;
  _z.set_zero(N_TOTAL);
  _z.set_sub_vec(0, zf);

  // mark starting time
  tms cstart;
  times(&cstart);

  // keep solving until we run out of time or all contact points are active
  while (true)
  {
    FILE_LOG(LOG_EVENT) << "Running QP solve iteration with " << (q.N_ACT_CONTACTS) << " active contacts" << std::endl;

    // solve the QP
    solve_qp_work(q, _z);

    // check whether we can mark any more contacts as active
    tms cstop;
    times(&cstop);
    if ((double) (cstop.tms_stime-cstart.tms_stime)/CLOCKS_PER_SEC > max_time ||
        q.N_ACT_CONTACTS == q.N_CONTACTS)
      break;

    // we can; mark next contact for solving
    q.N_ACT_K += q.contact_events[q.N_ACT_CONTACTS]->contact_NK/2;
    q.active_contacts[q.N_ACT_CONTACTS++] = true;
  }

  // apply (Poisson) restitution to contacts
  for (unsigned i=0, j=q.CN_IDX; i< q.N_ACT_CONTACTS; i++, j++)
    _z[j] *= ((double) 1.0 + q.contact_events[i]->contact_epsilon);

  // apply (Poisson) restitution to limits
  for (unsigned i=0, j=q.L_IDX; i< q.N_LIMITS; i++, j++)
    _z[j] *= ((double) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  q.update_from_stacked_qp(_z);

  // update Cn_v
  q.Cn_v += q.Cn_iM_CnT.mult(q.cn, _a);
  q.Cn_v += q.Cn_iM_CsT.mult(q.cs, _a);
  q.Cn_v += q.Cn_iM_CtT.mult(q.ct, _a);
  q.Cn_v += q.Cn_iM_LT.mult(q.l, _a);
  q.Cn_v += q.Cn_iM_JxT.mult(q.alpha_x, _a);

  // update Cs_v
  q.Cs_v += q.Cn_iM_CsT.transpose_mult(q.cn, _a);
  q.Cs_v += q.Cs_iM_CsT.mult(q.cs, _a);
  q.Cs_v += q.Cs_iM_CtT.mult(q.ct, _a);
  q.Cs_v += q.Cs_iM_LT.mult(q.l, _a);
  q.Cs_v += q.Cs_iM_JxT.mult(q.alpha_x, _a);

  // update Ct_v
  q.Ct_v += q.Cn_iM_CtT.transpose_mult(q.cn, _a);
  q.Ct_v += q.Cs_iM_CtT.transpose_mult(q.cs, _a);
  q.Ct_v += q.Ct_iM_CtT.mult(q.ct, _a);
  q.Ct_v += q.Ct_iM_LT.mult(q.l, _a);
  q.Ct_v += q.Ct_iM_JxT.mult(q.alpha_x, _a);

  // update L_v
  q.L_v += q.Cn_iM_LT.transpose_mult(q.cn, _a);
  q.L_v += q.Cs_iM_LT.transpose_mult(q.cs, _a);
  q.L_v += q.Ct_iM_LT.transpose_mult(q.ct, _a);
  q.L_v += q.L_iM_LT.mult(q.l, _a);
  q.L_v += q.L_iM_JxT.mult(q.alpha_x, _a);

  // update Jx_v
  q.Jx_v += q.Cn_iM_JxT.transpose_mult(q.cn, _a);
  q.Jx_v += q.Cs_iM_JxT.transpose_mult(q.cs, _a);
  q.Jx_v += q.Ct_iM_JxT.transpose_mult(q.ct, _a);
  q.Jx_v += q.L_iM_JxT.transpose_mult(q.l, _a);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, _a);

  // output results
  FILE_LOG(LOG_EVENT) << "results: " << std::endl;
  FILE_LOG(LOG_EVENT) << "cn: " << q.cn << std::endl;
  FILE_LOG(LOG_EVENT) << "cs: " << q.cs << std::endl;
  FILE_LOG(LOG_EVENT) << "ct: " << q.ct << std::endl;
  FILE_LOG(LOG_EVENT) << "l: " << q.l << std::endl;
  FILE_LOG(LOG_EVENT) << "alpha_x: " << q.alpha_x << std::endl;
  FILE_LOG(LOG_EVENT) << "new Cn_v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Cs_v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Ct_v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new L_v: " << q.L_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Jx_v: " << q.Jx_v << std::endl;

  // see whether another QP must be solved
  if (q.Cn_v.size() > 0 && *min_element(q.Cn_v.column_iterator_begin(), q.Cn_v.column_iterator_end()) < -TOL)
  {
    FILE_LOG(LOG_EVENT) << "minimum Cn*v: " << *min_element(q.Cn_v.column_iterator_begin(), q.Cn_v.column_iterator_end()) << std::endl;
    FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
    solve_qp_work(q, _z);
    q.update_from_stacked_qp(_z);
  }
  else if (q.L_v.size() > 0 && *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) < -TOL)
    {
      FILE_LOG(LOG_EVENT) << "minimum L*v: " << *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, _z);
      q.update_from_stacked_qp(_z);
    }
  else
  {
    pair<ColumnIteratord, ColumnIteratord> mm = boost::minmax_element(q.Jx_v.column_iterator_begin(), q.Jx_v.column_iterator_end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_EVENT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_EVENT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, _z);
      q.update_from_stacked_qp(_z);
    }
  }

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save contact impulses
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_events[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = q.contact_events[i]->contact_normal * q.cn[i];
    j += q.contact_events[i]->contact_tan1 * q.cs[i];
    j += q.contact_events[i]->contact_tan2 * q.ct[i];

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);    

    // transform the impulse to the global frame
    q.contact_events[i]->contact_impulse = Pose3d::transform(GLOBAL, jx);
  }

  // save limit impulses
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    q.limit_events[i]->limit_impulse = q.l[i]; 
    if (q.limit_events[i]->limit_upper)
      q.limit_events[i]->limit_impulse = -q.limit_events[i]->limit_impulse;
  }
}

/// Updates the solution by adding a contact
/**
 * TODO: prepare to delete this
void ImpactEventHandler::update_solution(const EventProblemData& q, const VectorNd& x, const vector<bool>& old_working_set, unsigned jidx, VectorNd& z)
{
  unsigned start, end;

  // setup NEW constants
  const unsigned N_VARS = q.N_LIMITS + q.L_IDX;
  const unsigned NP_IDX = N_VARS;
  const unsigned MU_IDX = NP_IDX + q.N_CONTACTS + q.N_LIMITS; 

  // setup old constants
  const unsigned OLD_N_CONTACTS = q.N_CONTACTS - 1;
  const unsigned OLD_CN_IDX = 0;
  const unsigned OLD_CS_IDX = OLD_N_CONTACTS;
  const unsigned OLD_CT_IDX = OLD_N_CONTACTS + OLD_CS_IDX;
  const unsigned OLD_NCS_IDX = OLD_N_CONTACTS + OLD_CT_IDX;
  const unsigned OLD_NCT_IDX = OLD_N_CONTACTS + OLD_NCS_IDX;
  const unsigned OLD_L_IDX = OLD_N_CONTACTS + OLD_NCT_IDX;
  const unsigned OLD_N_VARS = q.N_LIMITS + OLD_L_IDX;
  const unsigned OLD_NP_IDX = OLD_N_VARS;
  const unsigned OLD_MU_IDX = OLD_NP_IDX + OLD_N_CONTACTS + q.N_LIMITS; 

  // determine NK indices 
  unsigned OLD_MU_JIDX = OLD_MU_IDX;
  unsigned J_NK = 0;
  for (unsigned i=0, j=0; i< q.N_CONTACTS; i++)
  {
    if (!old_working_set[i])
      continue;
    if (j++ == jidx)
    {
      J_NK += q.contact_events[i]->contact_NK/2; 
      break;
    }
    OLD_MU_JIDX += q.contact_events[i]->contact_NK/2;
  }

  // initialize z -- x is the old solution to the LCP
  // NOTE: each contact adds six variables plus friction directions
  z.set_zero(x.size() + 6 + J_NK);

  // setup normal contact variables
  z.segment(q.CN_IDX, q.CN_IDX+jidx) = x.segment(OLD_CN_IDX, OLD_CN_IDX+jidx);
  start = q.CN_IDX + jidx + 1;
  end = start + OLD_CS_IDX - OLD_CN_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_CN_IDX+jidx, OLD_CS_IDX); 

  // setup positive first direction tangent contact variables
  z.segment(q.CS_IDX, q.CS_IDX+jidx) = x.segment(OLD_CS_IDX, OLD_CS_IDX+jidx);
  start = q.CS_IDX+jidx+1;
  end = start + OLD_CT_IDX - OLD_CS_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_CS_IDX+jidx, OLD_CT_IDX);

  // setup positive second direction tangent contact variables
  z.segment(q.CT_IDX, q.CT_IDX+jidx) = x.segment(OLD_CT_IDX, OLD_CT_IDX+jidx);
  start = q.CT_IDX+jidx+1;
  end = start + OLD_NCS_IDX - OLD_CT_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_CT_IDX+jidx, OLD_NCS_IDX);

  // setup negative first direction tangent contact variables
  z.segment(q.NCS_IDX, q.NCS_IDX+jidx) = x.segment(OLD_NCS_IDX, OLD_NCS_IDX+jidx);
  start = q.NCS_IDX+jidx+1;
  end = start + OLD_NCT_IDX - OLD_NCS_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_NCS_IDX+jidx, OLD_NCT_IDX);

  // setup negative second direction tangent contact variables AND
  // all other variables
  z.segment(q.NCT_IDX, q.NCT_IDX+jidx) = x.segment(OLD_NCT_IDX, OLD_NCT_IDX+jidx);
  start = q.NCT_IDX+jidx+1;
  end = start + OLD_NP_IDX - OLD_NCT_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_NCT_IDX+jidx, OLD_N_VARS);

  // setup constraint equation variables
  // constraint ordering: noninterpenetration, joint limits, 
  // contact friction
  // 1. noninterpenetration and joint limit constraints
  z.segment(NP_IDX, NP_IDX+jidx) = x.segment(OLD_NP_IDX, OLD_NP_IDX+jidx);
  start = NP_IDX+jidx+1;
  end = start + OLD_MU_IDX - OLD_NP_IDX - jidx;
  z.segment(start, end) = x.segment(OLD_NP_IDX+jidx, OLD_MU_IDX);

  // 2. contact friction constraints
  z.segment(MU_IDX, MU_IDX+OLD_MU_JIDX-OLD_MU_IDX) = x.segment(OLD_MU_IDX, OLD_MU_JIDX);  
  start = MU_IDX+OLD_MU_JIDX-OLD_MU_IDX+J_NK;
  end = start + x.size() - OLD_MU_JIDX;
  z.segment(start, end) = x.segment(OLD_MU_JIDX, x.size());
}
*/

/// Computes the kinetic energy of the system using the current impulse set 
double ImpactEventHandler::calc_ke(EventProblemData& q, const VectorNd& z)
{
  static VectorNd cn, cs, ct, l, alpha_x;

  // save the current impulses
  cn = q.cn;
  cs = q.cs;
  ct = q.ct;
  l = q.l;
  alpha_x = q.alpha_x;

  // update impulses
  q.update_from_stacked_qp(z);

  // calculate KE
  double KE = (double) 0.0;
  apply_impulses(q);
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    KE += q.super_bodies[i]->calc_kinetic_energy();

  // reset impulses
  q.cn = cn;
  q.cs = cs;
  q.ct = ct;
  q.l = l;
  q.alpha_x = alpha_x;

  return KE;
}

/// Solves the quadratic program (does all of the work) 
/**
 * \note this is the version without joint friction forces
 * \param z the solution is returned here; zeros are returned at appropriate
 *        places for inactive contacts
 */
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorNd& z)
{
  // implicit constraints not handled at the moment
  assert(q.N_CONSTRAINT_EQNS_IMP == 0);

  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cs': " << std::endl << q.Cn_iM_CsT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Ct': " << std::endl << q.Cn_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * Cs': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * Ct': " << std::endl << q.Cs_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * L': " << std::endl << q.Cs_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Ct * inv(M) * Ct': " << std::endl << q.Ct_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Ct * inv(M) * L': " << std::endl << q.Ct_iM_LT;
  FILE_LOG(LOG_EVENT) << "  L * inv(M) * L': " << std::endl << q.L_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Cs * v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Ct * v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  L * v: " << q.L_v << std::endl;

  // determine number of active variables
  const unsigned N_ACT_CONTACTS = q.N_ACT_CONTACTS;

  // setup new indices
  const unsigned xCN_IDX = 0;
  const unsigned xCS_IDX = xCN_IDX + N_ACT_CONTACTS;
  const unsigned xCT_IDX = xCS_IDX + N_ACT_CONTACTS;
  const unsigned xNCS_IDX = xCT_IDX + N_ACT_CONTACTS;
  const unsigned xNCT_IDX = xNCS_IDX + N_ACT_CONTACTS;
  const unsigned xL_IDX = xNCT_IDX + N_ACT_CONTACTS;
  const unsigned N_ACT_VARS = xL_IDX + q.N_LIMITS;

  // init the QP matrix and vector
  const unsigned N_INEQUAL = q.N_ACT_CONTACTS + q.N_ACT_K + q.N_LIMITS + 1;
  _MM.set_zero(N_ACT_VARS + N_INEQUAL, N_ACT_VARS + N_INEQUAL);
  _qq.resize(_MM.rows());

  // get useful blocks of _MM and segments of _qq
  SharedMatrixNd A = _MM.block(N_ACT_VARS, _MM.rows(), 0, N_ACT_VARS);
  SharedMatrixNd H = _MM.block(0, N_ACT_VARS, 0, N_ACT_VARS);
  SharedVectorNd c = _qq.segment(0, N_ACT_VARS);
  SharedVectorNd nb = _qq.segment(N_ACT_VARS, _qq.size()).set_zero();

  // setup row (block) 1 -- Cn * iM * [Cn' Cs Ct' -Cs' -Ct' L' ]
  unsigned col_start = 0, col_end = q.N_ACT_CONTACTS;
  unsigned row_start = 0, row_end = q.N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cn_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_ACT_CONTACTS;
  col_start = 0; col_end = q.N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_ACT_CONTACTS;
  col_start = 0; col_end = q.N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Ct_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 4 -- -Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_ACT_CONTACTS;
  SharedMatrixNd NCs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 5 -- -Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_ACT_CONTACTS;
  SharedMatrixNd NCt_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block 6) -- L * iM *  [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_LIMITS;
  col_start = 0; col_end = q.N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd L_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_ACT_CONTACTS;
  SharedMatrixNd L_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // copy to row block 1 (contact normals)
  q.Cn_iM_CnT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CnT);
  q.Cn_iM_CsT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CsT);
  q.Cn_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CtT);
  (Cn_iM_NCsT = Cn_iM_CsT).negate();
  (Cn_iM_NCtT = Cn_iM_CtT).negate();
  q.Cn_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, q.N_LIMITS, Cn_iM_LT);
  
  // copy to row block 2 (first contact tangents)
  MatrixNd::transpose(Cn_iM_CsT, Cs_iM_CnT);
  q.Cs_iM_CsT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cs_iM_CsT);
  q.Cs_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cs_iM_CtT);
  (Cs_iM_NCsT = Cs_iM_CsT).negate();
  (Cs_iM_NCtT = Cs_iM_CtT).negate();
  q.Cs_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, q.N_LIMITS, Cs_iM_LT);

  // copy to row block 3 (second contact tangents)
  MatrixNd::transpose(Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(Cs_iM_CtT, Ct_iM_CsT);
  q.Ct_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Ct_iM_CtT);
  (Ct_iM_NCsT = Ct_iM_CsT).negate();
  (Ct_iM_NCtT = Ct_iM_CtT).negate();
  q.Ct_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, q.N_LIMITS, Ct_iM_LT);

  // copy to row block 4 (negative first contact tangents)
  (NCs_block = Cs_block).negate();

  // copy to row block 5 (negative second contact tangents)
  (NCt_block = Ct_block).negate();

  // copy to row block 6 (limits)
  MatrixNd::transpose(Cn_iM_LT, L_iM_CnT);
  MatrixNd::transpose(Cs_iM_LT, L_iM_CsT);
  MatrixNd::transpose(Ct_iM_LT, L_iM_CtT);
  (L_iM_NCsT = L_iM_CsT).negate();
  (L_iM_NCtT = L_iM_CtT).negate();
  q.L_iM_LT.get_sub_mat(0, q.N_LIMITS, 0, q.N_LIMITS, L_iM_LT);

  // get shared vectors to components of c
  SharedVectorNd Cn_v = c.segment(xCN_IDX, xCS_IDX);
  SharedVectorNd Cs_v = c.segment(xCS_IDX, xCT_IDX);
  SharedVectorNd Ct_v = c.segment(xCT_IDX, xNCS_IDX);
  SharedVectorNd nCs_v = c.segment(xNCS_IDX, xNCT_IDX);
  SharedVectorNd nCt_v = c.segment(xNCT_IDX, xL_IDX);
  SharedVectorNd L_v = c.segment(xL_IDX, N_ACT_VARS);

  // setup c 
  q.Cn_v.get_sub_vec(0, N_ACT_CONTACTS, Cn_v);
  q.Cs_v.get_sub_vec(0, N_ACT_CONTACTS, Cs_v);
  q.Ct_v.get_sub_vec(0, N_ACT_CONTACTS, Ct_v);
  (nCs_v = Cs_v).negate();
  (nCt_v = Ct_v).negate();
  L_v = q.L_v;         

  // ------- setup A/-b -------
  // setup the Cn*v+ >= 0 constraint
  // Cn*(inv(M)*impulses + v) >= 0, Cn*inv(M)*impulses >= -Cn*v
  FILE_LOG(LOG_EVENT) << "Cn block: " << std::endl << Cn_block;
  row_start = 0; row_end = q.N_ACT_CONTACTS;
  A.block(row_start, row_end, xCN_IDX, N_ACT_VARS) = Cn_block;
  FILE_LOG(LOG_EVENT) << "A: " << std::endl << A;
  nb.set_sub_vec(row_start, Cn_v);
  row_start = row_end; row_end += q.N_LIMITS;  

  // setup the L*v+ >= 0 constraint
  A.block(row_start, row_end, xCN_IDX, N_ACT_VARS) = L_block;
  nb.set_sub_vec(row_start, q.L_v);
  row_start = row_end; row_end += q.N_ACT_CONTACTS;  

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0; i< q.N_ACT_CONTACTS; i++)
  {
    // initialize the contact velocity
    double vel = std::sqrt(sqr(q.Cs_v[i]) + sqr(q.Ct_v[i]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK/2; j++)
    {
      double theta = (double) j/(q.contact_events[i]->contact_NK/2-1) * M_PI_2;
      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      A(row_start, xCN_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row_start, xCS_IDX+i) = -ct;
      A(row_start, xNCS_IDX+i) = -ct;
      A(row_start, xCT_IDX+i) = -st;
      A(row_start, xNCT_IDX+i) = -st;

      // setup the viscous friction component
      nb[row_start] = q.contact_events[i]->contact_mu_viscous * vel;
      row_start++;
    }
  }

  // setup the normal velocity constraint
  // 1'N*v+ <= kappa  (in our form) -1'N*v+ >= -kappa
  for (unsigned i=0; i< Cn_block.columns(); i++)
  {
    SharedConstVectorNd Cn_col = Cn_block.column(i);
    A(row_start, i) = -std::accumulate(Cn_col.begin(), Cn_col.end(), 0.0);
  }
  nb[row_start] = -std::accumulate(q.Cn_v.row_iterator_begin(), q.Cn_v.row_iterator_end(), 0.0) + q.kappa;

  // set A = -A'
  SharedMatrixNd AT = _MM.block(0, N_ACT_VARS, N_ACT_VARS, _MM.rows());
  MatrixNd::transpose(A, AT);
  AT.negate();

  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "-b vector: " << nb << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << _MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << _qq << std::endl; 

  // init z to zero
  z.set_zero(_qq.rows());

  // use warm starting, if possible; last solution should be for one fewer
  // contact. New variables (primal and dual) will have values of zero.
  const unsigned N_LACT_CONTACTS = N_ACT_CONTACTS - 1;
  const unsigned N_LACT_VARS = N_LACT_CONTACTS*5 + q.N_LIMITS;
  if (_zlast.rows() == N_LACT_VARS)
  {
    // setup last indices
    const unsigned yCN_IDX = 0;
    const unsigned yCS_IDX = N_LACT_CONTACTS;
    const unsigned yCT_IDX = N_LACT_CONTACTS*2;
    const unsigned yNCS_IDX = N_LACT_CONTACTS*3;
    const unsigned yNCT_IDX = N_LACT_CONTACTS*4;
    const unsigned yL_IDX = N_LACT_CONTACTS*5;

    // populate normal contact forces
    SharedVectorNd cn_last = _zlast.segment(yCN_IDX, yCS_IDX);
    z.set_sub_vec(xCN_IDX, cn_last);

    // populate first direction positive tangent contact forces
    SharedVectorNd cs_last = _zlast.segment(yCS_IDX, yCT_IDX);
    z.set_sub_vec(xCS_IDX, cs_last);

    // populate second direction positive tangent contact forces
    SharedVectorNd ct_last = _zlast.segment(yCT_IDX, yNCS_IDX);
    z.set_sub_vec(xCT_IDX, ct_last);
  
    // populate first direction negative tangent contact forces
    SharedVectorNd ncs_last = _zlast.segment(yNCS_IDX, yNCT_IDX);
    z.set_sub_vec(xNCS_IDX, ncs_last);

    // populate second direction negative tangent contact forces
    SharedVectorNd nct_last = _zlast.segment(yNCT_IDX, yL_IDX);
    z.set_sub_vec(xNCT_IDX, nct_last);

    // populate limit forces
    SharedVectorNd l_last = _zlast.segment(yL_IDX,yL_IDX+q.N_LIMITS);
    z.set_sub_vec(xL_IDX, l_last);

    // populate Cn*v+ >= 0 constraint multipliers
    SharedVectorNd c1_last = _zlast.segment(N_LACT_VARS, N_LACT_VARS+N_LACT_CONTACTS);
    z.set_sub_vec(N_ACT_VARS, c1_last);

    // populate L*v+ >= 0 constraint multipliers
    SharedVectorNd c2_last = _zlast.segment(N_LACT_VARS+N_LACT_CONTACTS, N_LACT_VARS+N_LACT_CONTACTS+q.N_LIMITS);
    z.set_sub_vec(N_ACT_VARS+N_ACT_CONTACTS, c2_last);

    // populate friction coefficient constraint multipliers
    SharedVectorNd c3_last = _zlast.segment(N_LACT_VARS+N_LACT_CONTACTS+q.N_LIMITS, _zlast.rows());
    z.set_sub_vec(N_ACT_VARS+N_ACT_CONTACTS+q.N_LIMITS, c3_last);
  }

  // solve the LCP using Lemke's algorithm
  if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
    throw std::runtime_error("Unable to solve event QP!");

  // output reported LCP solution
  FILE_LOG(LOG_EVENT) << "LCP solution: " << z << std::endl;

  // store zlast
  _zlast = z;

  // get relevant forces
  SharedVectorNd cn = _zlast.segment(0, N_ACT_CONTACTS);
  SharedVectorNd cs = _zlast.segment(N_ACT_CONTACTS, N_ACT_CONTACTS*2);
  SharedVectorNd ct = _zlast.segment(N_ACT_CONTACTS*2, N_ACT_CONTACTS*3);
  SharedVectorNd ncs = _zlast.segment(N_ACT_CONTACTS*3, N_ACT_CONTACTS*4);
  SharedVectorNd nct = _zlast.segment(N_ACT_CONTACTS*4, N_ACT_CONTACTS*5);
  SharedVectorNd l = _zlast.segment(N_ACT_CONTACTS*5, N_ACT_CONTACTS*5+q.N_LIMITS);

  // put z in the expected format (full contact forces)
  z.set_zero(q.N_VARS);
  z.set_sub_vec(q.CN_IDX, cn);
  z.set_sub_vec(q.CS_IDX, cs);
  z.set_sub_vec(q.CT_IDX, ct);
  z.set_sub_vec(q.NCS_IDX, ncs);
  z.set_sub_vec(q.NCT_IDX, nct);
  z.set_sub_vec(q.L_IDX, l);

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  if (LOGGING(LOG_EVENT))
  {
    VectorNd workv;
    SharedVectorNd zsub = z.segment(0, c.rows());
    H.mult(zsub, workv) *= 0.5;
    workv += c;
    FILE_LOG(LOG_EVENT) << "(signed) computed energy dissipation: " << zsub.dot(workv) << std::endl;
  }    
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work() exited" << std::endl;
}

/// Updates problem matrices and vectors based on the working set
/**
 * TODO: prepare to remove this function
 */
/*
void ImpactEventHandler::update_problem(const EventProblemData& q, EventProblemData& qworking)
{
  static vector<unsigned> indices;
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // determine how many contacts we have
  const vector<bool>& working_set = qworking.contact_working_set;
  qworking.N_CONTACTS = std::count(working_set.begin(), working_set.end(), true);

  // fix vector of contact events
  qworking.contact_events.clear();
  for (unsigned i=0; i< q.contact_events.size(); i++)
    if (working_set[i])
      qworking.contact_events.push_back(q.contact_events[i]);
  assert(qworking.contact_events.size() == qworking.N_CONTACTS);

  // determine which contact constraints use a true friction cone
  qworking.N_K_TOTAL = 0;
  qworking.N_LIN_CONE = 0;
  for (unsigned i=0; i< qworking.contact_events.size(); i++)
    if (qworking.contact_events[i]->contact_NK < UINF)
    {
      qworking.N_K_TOTAL += qworking.contact_events[i]->contact_NK/2;
      qworking.N_LIN_CONE++;
    }

  // setup indices set
  indices.clear();
  for (unsigned i=0, j=0; i< working_set.size(); i++, j+=2)
    if (working_set[i])
      indices.push_back(i);

  // resize impulse vectors
  qworking.cn.set_zero(qworking.N_CONTACTS);
  qworking.cs.set_zero(qworking.N_CONTACTS);
  qworking.ct.set_zero(qworking.N_CONTACTS);

  // update indices
  qworking.CN_IDX = 0;
  qworking.CS_IDX = qworking.CN_IDX + qworking.N_CONTACTS;
  qworking.CT_IDX = qworking.CS_IDX + qworking.N_CONTACTS;
  qworking.NCS_IDX = qworking.CT_IDX + qworking.N_CONTACTS;
  qworking.NCT_IDX = qworking.NCS_IDX + qworking.N_LIN_CONE;
  qworking.L_IDX = qworking.NCT_IDX + qworking.N_LIN_CONE;
  qworking.ALPHA_X_IDX = qworking.L_IDX + qworking.N_LIMITS;
  qworking.N_VARS = qworking.ALPHA_X_IDX + qworking.N_CONSTRAINT_EQNS_IMP;

  // select appropriate parts of vectors 
  q.Cn_v.select(indices.begin(), indices.end(), qworking.Cn_v);
  q.Cs_v.select(indices.begin(), indices.end(), qworking.Cs_v);
  q.Ct_v.select(indices.begin(), indices.end(), qworking.Ct_v);

  // select appropriate parts of Cn matrices 
  q.Cn_iM_CnT.select_square(indices.begin(), indices.end(), qworking.Cn_iM_CnT);
  q.Cn_iM_CsT.select_square(indices.begin(), indices.end(), qworking.Cn_iM_CsT);
  q.Cn_iM_CtT.select_square(indices.begin(), indices.end(), qworking.Cn_iM_CtT);
  q.Cn_iM_LT.select_rows(indices.begin(), indices.end(), qworking.Cn_iM_LT);
  q.Cn_iM_JxT.select_rows(indices.begin(), indices.end(), qworking.Cn_iM_JxT);

  // select appropriate parts of Cs matrices
  q.Cs_iM_CsT.select_square(indices.begin(), indices.end(), qworking.Cs_iM_CsT);
  q.Cs_iM_CtT.select_square(indices.begin(), indices.end(), qworking.Cs_iM_CtT);
  q.Cs_iM_LT.select_rows(indices.begin(), indices.end(), qworking.Cs_iM_LT);
  q.Cs_iM_JxT.select_rows(indices.begin(), indices.end(), qworking.Cs_iM_JxT);

  // select appropriate parts of Ct matrices
  q.Ct_iM_CtT.select_square(indices.begin(), indices.end(), qworking.Ct_iM_CtT);
  q.Ct_iM_LT.select_rows(indices.begin(), indices.end(), qworking.Ct_iM_LT);
  q.Ct_iM_JxT.select_rows(indices.begin(), indices.end(), qworking.Ct_iM_JxT);
}
*/

