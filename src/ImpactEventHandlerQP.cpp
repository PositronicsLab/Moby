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
#include <Moby/LCPSolverException.h>
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

  // setup last successful solution
  _zsuccess = _z;

  // mark starting time
  tms cstart;
  clock_t start = times(&cstart);

  // reset warm starting variables
  _last_contacts = _last_contact_constraints = _last_contact_nk = 0;
  _last_limits = 0;
 
  // keep solving until we run out of time or all contact points are active
  while (true)
  {
    FILE_LOG(LOG_EVENT) << "Running QP solve iteration with " << (q.N_ACT_CONTACTS) << " active contacts" << std::endl;

    // solve the QP
    try
    {
      solve_qp_work(q, _z);
    }
    catch (LCPSolverException e)
    {
      FILE_LOG(LOG_EVENT) << "Failed to solve QP: returning best solution so far" << std::endl;
      _z = _zsuccess;
      break;
    }

    // save our successful solve
    _zsuccess = _z;

    // get the elapsed time
    const long TPS = sysconf(_SC_CLK_TCK);
    tms cstop;
    clock_t stop = times(&cstop);
    double elapsed = (double) (stop - start)/TPS;
    FILE_LOG(LOG_EVENT) << "Elapsed time: " << elapsed << std::endl;

    // check whether we can mark any more contacts as active
    if (elapsed > max_time || q.N_ACT_CONTACTS == q.N_CONTACTS)
      break;

    // we can; mark next contact for solving
    q.N_ACT_K += q.contact_events[q.N_ACT_CONTACTS]->contact_NK/2;
    q.N_ACT_CONTACTS++;
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

  // get useful constants
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
  const unsigned N_INEQUAL = q.N_CONTACT_CONSTRAINTS + q.N_ACT_K + q.N_LIMITS + 1;
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

  // get contact constrained Cn blocks
  q.Cn_iM_CnT.select_rows(q.contact_constraints, _Cnstar_Cn);
  q.Cn_iM_CsT.select_rows(q.contact_constraints, _Cnstar_Cs);
  q.Cn_iM_CtT.select_rows(q.contact_constraints, _Cnstar_Ct);
  q.Cn_iM_LT.select_rows(q.contact_constraints, _Cnstar_L);
  q.Cn_v.select(q.contact_constraints, _Cnstar_v);
  SharedConstMatrixNd Cnstar_Cnx = _Cnstar_Cn.block(0, _Cnstar_Cn.rows(), 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd Cnstar_Csx = _Cnstar_Cs.block(0, _Cnstar_Cs.rows(), 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd Cnstar_Ctx = _Cnstar_Ct.block(0, _Cnstar_Ct.rows(), 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd Cnstar_Lx = _Cnstar_L.block(0, _Cnstar_L.rows(), 0, q.N_LIMITS);

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
  row_start = 0; row_end = q.N_CONTACT_CONSTRAINTS;
  A.set_sub_mat(row_start, xCN_IDX, Cnstar_Cnx);
  A.set_sub_mat(row_start, xCS_IDX, Cnstar_Csx);
  A.set_sub_mat(row_start, xCT_IDX, Cnstar_Ctx);
  A.set_sub_mat(row_start, xNCS_IDX, Cnstar_Csx);
  A.set_sub_mat(row_start, xNCT_IDX, Cnstar_Ctx);
  A.set_sub_mat(row_start, xL_IDX, Cnstar_Lx);
  A.block(row_start, row_end, xNCS_IDX, xL_IDX).negate();
  SharedConstMatrixNd full_Cn_block = A.block(row_start, row_end, xCN_IDX, N_ACT_VARS);
  FILE_LOG(LOG_EVENT) << "A: " << std::endl << A;
  nb.set_sub_vec(row_start, _Cnstar_v);
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
  for (unsigned i=0; i< full_Cn_block.columns(); i++)
  {
    SharedConstVectorNd Cn_col = full_Cn_block.column(i);
    A(row_start, i) = -std::accumulate(Cn_col.begin(), Cn_col.end(), 0.0);
  }
  nb[row_start] = -std::accumulate(_Cnstar_v.row_iterator_begin(), _Cnstar_v.row_iterator_end(), 0.0) + q.kappa;

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

  // use warm starting if possible
  // new variables (primal, dual) will have values of zero.
  const unsigned N_LAST_ACT_CONTACTS = _last_contacts;
  const unsigned N_LAST_LIMITS = _last_limits;
  const unsigned N_LAST_CC = _last_contact_constraints;
  const unsigned N_LAST_ACT_K = _last_contact_nk;
  const unsigned N_LAST_VARS = N_LAST_ACT_CONTACTS*5 + N_LAST_LIMITS;
  const unsigned N_LAST_INEQUAL = N_LAST_CC + N_LAST_ACT_K + N_LAST_LIMITS + 1;
  const unsigned N_LAST_TOTAL = N_LAST_VARS + N_LAST_INEQUAL;
  if (N_LAST_ACT_CONTACTS <= q.N_ACT_CONTACTS && 
      N_LAST_CC <= q.N_CONTACT_CONSTRAINTS &&
      N_LAST_LIMITS <= q.N_LIMITS &&
      N_LAST_TOTAL == _zlast.rows())
  {
    // setup last indices
    const unsigned yCN_IDX = 0;
    const unsigned yCS_IDX = N_LAST_ACT_CONTACTS;
    const unsigned yCT_IDX = yCS_IDX + N_LAST_ACT_CONTACTS;
    const unsigned yNCS_IDX = yCT_IDX + N_LAST_ACT_CONTACTS;
    const unsigned yNCT_IDX = yNCS_IDX + N_LAST_ACT_CONTACTS;
    const unsigned yL_IDX = yNCT_IDX + N_LAST_ACT_CONTACTS;
    const unsigned lNPc_IDX = yL_IDX + N_LAST_LIMITS; // np lagrange indices
    const unsigned lNPl_IDX = lNPc_IDX + N_LAST_CC; // limit lagrange indices
    const unsigned lRest_IDX = lNPl_IDX + N_LAST_LIMITS; // remaining lagrange indices

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
    SharedVectorNd c1_last = _zlast.segment(lNPc_IDX, lNPc_IDX+N_LAST_CC);
    z.set_sub_vec(N_ACT_VARS, c1_last);

    // populate L*v+ >= 0 constraint multipliers
    SharedVectorNd c2_last = _zlast.segment(lNPl_IDX, lNPl_IDX+N_LAST_LIMITS);
    z.set_sub_vec(N_ACT_VARS+q.N_CONTACT_CONSTRAINTS, c2_last);

    // populate friction coefficient constraint multipliers and kappa constraint
    SharedVectorNd c3_last = _zlast.segment(lRest_IDX, _zlast.rows());
    z.set_sub_vec(N_ACT_VARS+q.N_CONTACT_CONSTRAINTS+q.N_LIMITS, c3_last);
  }

  // solve the LCP using Lemke's algorithm
  if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
    throw LCPSolverException();

  // output reported LCP solution
  FILE_LOG(LOG_EVENT) << "LCP solution: " << z << std::endl;

  // store zlast
  _zlast = z;

  // save number of variables
  _last_contacts = N_ACT_CONTACTS;
  _last_limits = q.N_LIMITS;
  _last_contact_constraints = q.N_CONTACT_CONSTRAINTS;
  _last_contact_nk = q.N_ACT_K;

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

  // compute full Cn_v solution *if necessary*
  if (q.N_CONTACT_CONSTRAINTS < q.N_CONTACTS) 
  {
    // get all contact blocks
    SharedConstMatrixNd full_Cn_Cn = q.Cn_iM_CnT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);  
    SharedConstMatrixNd full_Cn_Cs = q.Cn_iM_CsT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);  
    SharedConstMatrixNd full_Cn_Ct = q.Cn_iM_CtT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);  
    SharedConstMatrixNd full_Cn_L = q.Cn_iM_LT.block(0, q.N_CONTACTS, 0, q.N_LIMITS);  
    SharedConstVectorNd full_Cn_v = q.Cn_v.segment(0, q.N_CONTACTS);  

    // compute new velocity in the normal direction
    full_Cn_Cn.mult(cn, _new_Cn_v);
    _new_Cn_v += full_Cn_Cs.mult(cs, _workv);
    _new_Cn_v += full_Cn_Ct.mult(ct, _workv);
    _new_Cn_v -= full_Cn_Cs.mult(ncs, _workv);
    _new_Cn_v -= full_Cn_Ct.mult(nct, _workv);
    _new_Cn_v += full_Cn_L.mult(l, _workv);
    _new_Cn_v += full_Cn_v;

    // check whether new contacts need to be added to the set of contact constraints
    bool rerun = false;
    for (unsigned i=0; i< _new_Cn_v.rows(); i++)
      if (_new_Cn_v[i] < -NEAR_ZERO)
      {
        if (!q.contact_constraints[i])
        {
          // make it active 
          q.contact_constraints[i] = true;

          // update the number of contact constraints
          q.N_CONTACT_CONSTRAINTS++;

          // indicate to rerun the optimization
          rerun = true;

          // only do this once
          break;
        }
        else
        {
          // indicate we couldn't solve the problem to the desired tolerance
          throw LCPSolverException();
        }
      }

    // rerun the contact optimization if necessary
    if (rerun)
    {
      FILE_LOG(LOG_EVENT) << "-- constraint violation detected on unincorported constraint(s)" << std::endl;
      FILE_LOG(LOG_EVENT) << "   re-running with " << q.N_CONTACT_CONSTRAINTS " contact constraints << std::endl; 
      solve_qp_work(q, z);
    }
  }

  if (LOGGING(LOG_EVENT))
  {
    VectorNd workv;
    SharedVectorNd zsub = _zlast.segment(0, c.rows());
    H.mult(zsub, workv) *= 0.5;
    workv += c;
    FILE_LOG(LOG_EVENT) << "(signed) computed energy dissipation: " << zsub.dot(workv) << std::endl;
  }    
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work() exited" << std::endl;
}


