/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

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
void ImpactEventHandler::solve_qp(EventProblemData& q, double poisson_eps)
{
  static VectorNd z, tmp, tmp2;
  const double TOL = poisson_eps;

  // solve the QP
  solve_qp_work(q, z);

  // apply (Poisson) restitution to contacts
  for (unsigned i=0; i< q.N_CONTACTS; i++)
    z[i] *= ((double) 1.0 + q.contact_events[i]->contact_epsilon);

  // apply (Poisson) restitution to limits
  for (unsigned i=0; i< q.N_LIMITS; i++)
    z[q.N_CONTACTS*5+i] *= ((double) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  q.update_from_stacked(z);

  // update Cn_v
  q.Cn_v += q.Cn_iM_CnT.mult(q.cn, tmp);
  q.Cn_v += q.Cn_iM_CsT.mult(q.cs, tmp);
  q.Cn_v += q.Cn_iM_CtT.mult(q.ct, tmp);
  q.Cn_v += q.Cn_iM_LT.mult(q.l, tmp);
  q.Cn_v += q.Cn_iM_JxT.mult(q.alpha_x, tmp);

  // update Cs_v
  q.Cs_v += q.Cn_iM_CsT.transpose_mult(q.cn, tmp);
  q.Cs_v += q.Cs_iM_CsT.mult(q.cs, tmp);
  q.Cs_v += q.Cs_iM_CtT.mult(q.ct, tmp);
  q.Cs_v += q.Cs_iM_LT.mult(q.l, tmp);
  q.Cs_v += q.Cs_iM_JxT.mult(q.alpha_x, tmp);

  // update Ct_v
  q.Ct_v += q.Cn_iM_CtT.transpose_mult(q.cn, tmp);
  q.Ct_v += q.Cs_iM_CtT.transpose_mult(q.cs, tmp);
  q.Ct_v += q.Ct_iM_CtT.mult(q.ct, tmp);
  q.Ct_v += q.Ct_iM_LT.mult(q.l, tmp);
  q.Ct_v += q.Ct_iM_JxT.mult(q.alpha_x, tmp);

  // update L_v
  q.L_v += q.Cn_iM_LT.transpose_mult(q.cn, tmp);
  q.L_v += q.Cs_iM_LT.transpose_mult(q.cs, tmp);
  q.L_v += q.Ct_iM_LT.transpose_mult(q.ct, tmp);
  q.L_v += q.L_iM_LT.mult(q.l, tmp);
  q.L_v += q.L_iM_JxT.mult(q.alpha_x, tmp);

  // update Jx_v
  q.Jx_v += q.Cn_iM_JxT.transpose_mult(q.cn, tmp);
  q.Jx_v += q.Cs_iM_JxT.transpose_mult(q.cs, tmp);
  q.Jx_v += q.Ct_iM_JxT.transpose_mult(q.ct, tmp);
  q.Jx_v += q.L_iM_JxT.transpose_mult(q.l, tmp);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, tmp);

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
    solve_qp_work(q, z);
    q.update_from_stacked(z);
  }
  else 
    if (q.L_v.size() > 0 && *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) < -TOL)
    {
      FILE_LOG(LOG_EVENT) << "minimum L*v: " << *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      q.update_from_stacked(z);
    }
  else
  {
    pair<ColumnIteratord, ColumnIteratord> mm = boost::minmax_element(q.Jx_v.column_iterator_begin(), q.Jx_v.column_iterator_end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_EVENT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_EVENT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      q.update_from_stacked(z);
    }
  }

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save normal contact impulses
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

/// Solves the quadratic program (does all of the work) for given active contact set
/**
 * \note this is the version without joint friction forces
 */
/*
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorNd& z)
{
  SAFESTATIC MatrixNd sub, t1, t2, t3, neg1, A, AR, R, RTH;
  SAFESTATIC MatrixNd H, MM;
  SAFESTATIC VectorNd negv, c, qq, nb, tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_IMP = q.N_CONSTRAINT_EQNS_IMP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // setup variable indices
  const unsigned CN_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + L_IDX;
  const unsigned N_VARS = N_CONSTRAINT_EQNS_IMP + ALPHA_X_IDX;

  // first, solve for impulses that satisfy explicit constraint equations
  // and compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_IMP > 0)
  {
    // compute the homogeneous solution
    A. = q.Jx_iM_JxT);
    z. = q.Jx_v).negate();
    try
    {
      LinAlg::solve_LS_fast1(A, z);
    }
    catch (NumericalException e)
    {
      A. = q.Jx_iM_JxT);
      LinAlg::solve_LS_fast2(A, z);
    }

    // compute the nullspace
    A.resize(N_CONSTRAINT_EQNS_IMP, N_VARS);
    MatrixNd::transpose(q.Cn_iM_JxT, t1);
    MatrixNd::transpose(q.Dc_iM_JxT, t2);
    MatrixNd::transpose(q.L_iM_JxT, t3);
    neg1. = t2).negate();
    A.set_sub_row_block(0, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);
    LinAlg::nullspace(A, R);
  }
  else
  {
    R.set_zero(N_VARS,N_VARS);
    for (unsigned i=0; i< N_VARS; i++) 
      R(i,i) = (double) 1.0;
    z.set_zero(N_VARS);
  }

  // get number of qp variables
  const unsigned N_VARS = R.columns();

  // init the QP matrix and vector
  const unsigned N_INEQUAL = N_CONTACTS + N_K_TOTAL + N_LIMITS + 1;
  H.resize(N_VARS, N_VARS);
  c.resize(H.rows());
  A.set_zero(N_INEQUAL, N_VARS);
  nb.set_zero(N_INEQUAL);
  MM.set_zero(N_VARS + N_INEQUAL, N_VARS + N_INEQUAL);
  qq.resize(MM.rows());

  // setup [Q M'; -M 0]
  unsigned col = 0, row = 0;

  // row (block) 1 -- Cn * iM * [Cn' Dc' -Dc' L' Jx']
  neg1. = q.Cn_iM_DcT).negate();
  H.set_sub_row_block(0, &q.Cn_iM_CnT, &q.Cn_iM_DcT, &neg1, &q.Cn_iM_LT, 
                      &q.Cn_iM_JxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Cn' Dc' -Dc' L' Jx']
  MatrixNd::transpose(q.Cn_iM_DcT, t1);
  neg1. = q.Dc_iM_DcT).negate();
  H.set_sub_row_block(row, &t1, &q.Dc_iM_DcT, &neg1, &q.Dc_iM_LT, 
                      &q.Dc_iM_JxT);

  // row (block 3) -- negated block 2
  H.get_sub_mat(row, row+N_CONTACTS*2, 0, H.columns(), sub);
  H.set_sub_mat(row+N_CONTACTS*2, 0, sub.negate());
  row += N_CONTACTS*4;

  // row (block 4) -- L * iM * [Cn' Dc' -Dc' L Jx']
  MatrixNd::transpose(q.Cn_iM_LT, t1);
  MatrixNd::transpose(q.Dc_iM_LT, t2);
  neg1. = t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &q.L_iM_LT, &q.L_iM_JxT);
  row += N_LIMITS;
  
  // row (block 5) -- Jx * iM * [Cn' Dc' -Dc' L Jx']
  MatrixNd::transpose(q.Cn_iM_JxT, t1);
  MatrixNd::transpose(q.Dc_iM_JxT, t2);
  MatrixNd::transpose(q.L_iM_JxT, t3);
  neg1. = t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);

  // setup c 
  c.set_sub_vec(CN_IDX, q.Cn_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  negv. = q.Dc_v).negate();
  c.set_sub_vec(NBETA_C_IDX, negv);           
  c.set_sub_vec(L_IDX, q.L_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // setup the Cn*v+ >= 0 constraint
  // Cn*(inv(M)*impulses + v) >= 0, Cn*inv(M)*impulses >= -Cn*v
  row = 0; col = 0;
  nb.set_sub_vec(row, q.Cn_v);
  H.get_sub_mat(CN_IDX, CN_IDX+N_CONTACTS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the L*v+ >= 0 constraint
  nb.set_sub_vec(row, q.L_v);
  H.get_sub_mat(L_IDX, L_IDX+N_LIMITS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, k=0; i< N_CONTACTS; i++, k+= 2)
  {
    // initialize the contact velocity
    double vel = std::sqrt(sqr(q.Dc_v[k]) + sqr(q.Dc_v[k+1]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK/2; j++)
    {
      double theta = (double) j/(q.contact_events[i]->contact_NK/2-1) * M_PI_2;
      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      A(row, CN_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row, BETA_C_IDX+k) = -ct;
      A(row, NBETA_C_IDX+k) = -ct;
      A(row, BETA_C_IDX+k+1) = -st;
      A(row, NBETA_C_IDX+k+1) = -st;

      // setup the viscous friction component
      nb[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }
  }

  // setup the normal velocity constraint
  // 1'N*v+ <= 1'N*v-  (in our form) -1'N*v+ >= -1'N*v- 
  for (unsigned i=0; i< Cn_block.columns(); i++)
  {
    SharedConstVectorNd Cn_col = Cn_block.column(i);
    A(row, i) = -std::accumulate(Cn_col.begin(), Cn_col.end(), 0.0);
  }
  nb[row] = q.kappa;

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H);
  R.transpose_mult(c, tmpv);
  c. = tmpv);
  RTH.mult(z, tmpv);
  c += tmpv;

  // setup constraints A*x >= b in nullspace, yielding
  // A(R*y + z) >= b, yielding R*y >= b - A*z
  A.mult(z, tmpv);
  nb += tmpv;
  A.mult(R, AR);

  // setup the LCP matrix
  MM.set_sub_mat(0, 0, H);
  MM.set_sub_mat(N_VARS, 0, AR);
  MM.set_sub_mat(0, N_VARS, AR.negate(), true);

  // setup the LCP vector
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N_VARS, nb);

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Dc': " << std::endl << q.Cn_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Jx': " << std::endl << q.Cn_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * L': " << std::endl << q.Dc_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  L * inv(M) * L': " << std::endl << q.L_iM_LT;
  FILE_LOG(LOG_EVENT) << "  L * inv(M) * Jx': " << std::endl << q.L_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!_lcp.lcp_lemke_regularized(MM, qq, tmpv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_EVENT) << "LCP solution: " << tmpv << std::endl; 
  tmpv.get_sub_vec(0, N_VARS, y);
  R.mult(y, tmpv);
  z += tmpv;

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() exited" << std::endl;
}
*/

/// Updates the solution by adding a contact
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
  q.update_from_stacked(z);

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
 */
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorNd& z)
{
  // if there are implicit constraints, cannot use the fast method
  if (q.N_CONSTRAINT_EQNS_IMP > 0)
  {
    solve_qp_work_general(q, z);
    return;
  } 

  // determine the number of QP variables used for contacts
  unsigned N_QP_CONTACT_VARS = q.N_CONTACTS*6 + q.N_K_TOTAL;

  // the algorithm below is inefficient 
  solve_qp_work_ijoints(q, z);

/*
  // setup the working set -- set only first contact to active
  qworking = q;
  vector<bool>& working_set = qworking.contact_working_set;
  working_set.resize(q.N_CONTACTS);
  std::fill(working_set.begin(), working_set.end(), false);
  working_set[0] = true;
  update_problem(q, qworking);

  // solve using warm starting
  solve_qp_work_ijoints(qworking, z);
  double KE = calc_ke(qworking, z);

  // begin looping
  for (unsigned j=1; j< q.N_CONTACTS; j++)
  {
    // see whether contact is already in the working set _or_
    // setting contact forces for y to zero is acceptible
    if (working_set[j] || opt_satisfied(q, working_set, KE, z, j))
      continue;

    // if we're here, optimization required adding contact j to working set
    // and z and KE have been updated 
    working_set[j] = true;
    
    // reset j now (at beginning of next 'for' loop, j = 1)
    j = 0;
  }

  FILE_LOG(LOG_EVENT) << "-- reduced contact events from " << q.N_CONTACTS << " to " << std::count(working_set.begin(), working_set.end(), true) << std::endl;

  // update qworking 
  update_problem(q, qworking);
  q = qworking;
*/
}

/// Updates problem matrices and vectors based on the working set
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
  qworking.CS_U_IDX = qworking.NCT_IDX + qworking.N_LIN_CONE;
  qworking.CT_U_IDX = qworking.CS_U_IDX + qworking.N_TRUE_CONE;
  qworking.L_IDX = qworking.CT_U_IDX + qworking.N_TRUE_CONE;
  qworking.BETA_T_IDX = qworking.L_IDX + qworking.N_LIMITS;
  qworking.ALPHA_X_IDX = qworking.BETA_T_IDX + qworking.N_CONSTRAINT_DOF_EXP;
  qworking.BETA_X_IDX = qworking.ALPHA_X_IDX + qworking.N_CONSTRAINT_EQNS_IMP;
  qworking.N_VARS = qworking.BETA_X_IDX + qworking.N_CONSTRAINT_DOF_IMP;

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
  q.Cn_iM_DxT.select_rows(indices.begin(), indices.end(), qworking.Cn_iM_DxT);

  // select appropriate parts of Cs matrices
  q.Cs_iM_CsT.select_square(indices.begin(), indices.end(), qworking.Cs_iM_CsT);
  q.Cs_iM_CtT.select_square(indices.begin(), indices.end(), qworking.Cs_iM_CtT);
  q.Cs_iM_LT.select_rows(indices.begin(), indices.end(), qworking.Cs_iM_LT);
  q.Cs_iM_JxT.select_rows(indices.begin(), indices.end(), qworking.Cs_iM_JxT);
  q.Cs_iM_DxT.select_rows(indices.begin(), indices.end(), qworking.Cs_iM_DxT);

  // select appropriate parts of Ct matrices
  q.Ct_iM_CtT.select_square(indices.begin(), indices.end(), qworking.Ct_iM_CtT);
  q.Ct_iM_LT.select_rows(indices.begin(), indices.end(), qworking.Ct_iM_LT);
  q.Ct_iM_JxT.select_rows(indices.begin(), indices.end(), qworking.Ct_iM_JxT);
  q.Ct_iM_DxT.select_rows(indices.begin(), indices.end(), qworking.Ct_iM_DxT);
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work_ijoints(EventProblemData& q, VectorNd& z)
{
  MatrixNd MM;
  VectorNd qq;

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

  // init the QP matrix and vector
  const unsigned N_INEQUAL = q.N_CONTACTS + q.N_K_TOTAL + q.N_LIMITS + 1;
  MM.set_zero(q.N_VARS + N_INEQUAL, q.N_VARS + N_INEQUAL);
  qq.resize(MM.rows());

  // get useful blocks of MM and segments of qq
  SharedMatrixNd A = MM.block(q.N_VARS, MM.rows(), 0, q.N_VARS);
  SharedMatrixNd H = MM.block(0, q.N_VARS, 0, q.N_VARS);
  SharedVectorNd c = qq.segment(0, q.N_VARS);
  SharedVectorNd nb = qq.segment(q.N_VARS, qq.size()).set_zero();

  // setup row (block) 1 -- Cn * iM * [Cn' Cs Ct' -Cs' -Ct' L' ]
  unsigned col_start = 0, col_end = q.N_CONTACTS;
  unsigned row_start = 0, row_end = q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cn_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONTACTS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_CONTACTS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Ct_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 4 -- -Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_CONTACTS;
  SharedMatrixNd NCs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 5 -- -Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_CONTACTS;
  SharedMatrixNd NCt_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block 6) -- L * iM *  [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += q.N_LIMITS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // copy to row block 1 (contact normals)
  Cn_iM_CnT = q.Cn_iM_CnT;
  Cn_iM_CsT = q.Cn_iM_CsT;
  Cn_iM_CtT = q.Cn_iM_CtT;
  (Cn_iM_NCsT = Cn_iM_CsT).negate();
  (Cn_iM_NCtT = Cn_iM_CtT).negate();
  Cn_iM_LT = q.Cn_iM_LT;
  
  // copy to row block 2 (first contact tangents)
  MatrixNd::transpose(q.Cn_iM_CsT, Cs_iM_CnT);
  Cs_iM_CsT = q.Cs_iM_CsT;
  Cs_iM_CtT = q.Cs_iM_CtT;
  (Cs_iM_NCsT = Cs_iM_CsT).negate();
  (Cs_iM_NCtT = Cs_iM_CtT).negate();
  Cs_iM_LT = q.Cs_iM_LT;

  // copy to row block 3 (second contact tangents)
  MatrixNd::transpose(q.Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(q.Cs_iM_CtT, Ct_iM_CsT);
  Ct_iM_CtT = q.Ct_iM_CtT;
  (Ct_iM_NCsT = Ct_iM_CsT).negate();
  (Ct_iM_NCtT = Ct_iM_CtT).negate();
  Ct_iM_LT = q.Ct_iM_LT;

  // copy to row block 4 (negative first contact tangents)
  (NCs_block = Cs_block).negate();

  // copy to row block 5 (negative second contact tangents)
  (NCt_block = Ct_block).negate();

  // copy to row block 6 (limits)
  MatrixNd::transpose(q.Cn_iM_LT, L_iM_CnT);
  MatrixNd::transpose(q.Cs_iM_LT, L_iM_CsT);
  MatrixNd::transpose(q.Ct_iM_LT, L_iM_CtT);
  (L_iM_NCsT = L_iM_CsT).negate();
  (L_iM_NCtT = L_iM_CtT).negate();
  L_iM_LT = q.L_iM_LT;

  // setup c 
  c.set_sub_vec(q.CN_IDX, q.Cn_v);         
  c.set_sub_vec(q.CS_IDX, q.Cs_v);         
  c.set_sub_vec(q.CT_IDX, q.Ct_v);         
  (c.segment(q.NCS_IDX, q.NCT_IDX) = q.Cs_v).negate();
  (c.segment(q.NCT_IDX, q.L_IDX) = q.Ct_v).negate();
  c.set_sub_vec(q.L_IDX, q.L_v);         

  // ------- setup A/-b -------
  // setup the Cn*v+ >= 0 constraint
  // Cn*(inv(M)*impulses + v) >= 0, Cn*inv(M)*impulses >= -Cn*v
  FILE_LOG(LOG_EVENT) << "Cn block: " << std::endl << Cn_block;
  row_start = 0; row_end = q.N_CONTACTS;
  A.block(row_start, row_end, q.CN_IDX, q.N_VARS) = Cn_block;
  FILE_LOG(LOG_EVENT) << "A: " << std::endl << A;
  nb.set_sub_vec(row_start, q.Cn_v);
  row_start = row_end; row_end += q.N_LIMITS;  

  // setup the L*v+ >= 0 constraint
  A.block(row_start, row_end, q.CN_IDX, q.N_VARS) = L_block;
  nb.set_sub_vec(row_start, q.L_v);
  row_start = row_end; row_end += q.N_CONTACTS;  

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // initialize the contact velocity
    double vel = std::sqrt(sqr(q.Cs_v[i]) + sqr(q.Ct_v[i]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK/2; j++)
    {
      double theta = (double) j/(q.contact_events[i]->contact_NK/2-1) * M_PI_2;
      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      A(row_start, q.CN_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row_start, q.CS_IDX+i) = -ct;
      A(row_start, q.NCS_IDX+i) = -ct;
      A(row_start, q.CT_IDX+i) = -st;
      A(row_start, q.NCT_IDX+i) = -st;

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
  SharedMatrixNd AT = MM.block(0, q.N_VARS, q.N_VARS, MM.rows());
  MatrixNd::transpose(A, AT);
  AT.negate();

  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "-b vector: " << nb << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  solve_lcp(MM, qq, z);

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  if (LOGGING(LOG_EVENT))
  {
    VectorNd workv;
    SharedVectorNd zsub = z.segment(0, c.size());
    H.mult(zsub, workv) *= 0.5;
    workv += c;
    FILE_LOG(LOG_EVENT) << "(signed) computed energy dissipation: " << zsub.dot(workv) << std::endl;
  }    
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_ijoints() exited" << std::endl;
}

/// Solves the linear complementarity problem intelligently (using sparse matrices, if possible)
void ImpactEventHandler::solve_lcp(const MatrixNd& M, const VectorNd& q, VectorNd& z)
{
  #if 0
  const double SPARSE_PCT = 1.00;  // if 85% of matrix is sparse, triggers sparse arithmetic

  // get the infinity-norm of the matrix
  const double MINF = M.norm_inf();

  // setup a sensible zero tolerance
  const double ZERO_TOL = std::max(1.0, MINF) * std::numeric_limits<double>::epsilon() * M.rows();

  // get the number of zeros in the matrix
  unsigned nz = 0;
  for (ColumnIteratord_const i = M.column_iterator_begin(); i != i.end(); i++)
    if (std::fabs(*i) < ZERO_TOL)
      nz++;

  // use the proper method depending on sparsity 
  FILE_LOG(LOG_EVENT) << "Order of LCP matrix: " << M.rows() << std::endl;
  FILE_LOG(LOG_EVENT) << "Sparsity percentage (higher = sparser): " << ((double) nz/M.size()) << std::endl;

  if ((double) nz/M.size() > SPARSE_PCT)
  {
    // create the new sparse matrix
    SparseMatrixNd sM(SparseMatrixNd::eCSC, M, ZERO_TOL);
    _lcp.lcp_lemke_regularized(sM, q, z);
  }
  else
    if (!_lcp.lcp_lemke_regularized(M, q, z))
      throw std::runtime_error("Unable to solve event QP!");
  #else
  if (!_lcp.lcp_lemke_regularized(M, q, z))
    throw std::runtime_error("Unable to solve event QP!");

  #endif
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work_general(EventProblemData& q, VectorNd& z)
{
  MatrixNd A, AR, R, RTH, H, MM;
  VectorNd c, qq, workv, y;

  // first, solve for impulses that satisfy implicit constraint equations
  // and compute the appropriate nullspace 
  if (q.N_CONSTRAINT_EQNS_IMP > 0)
  {
    // compute the homogeneous solution
    A = q.Jx_iM_JxT;
    (z = q.Jx_v).negate();
    try
    {
      _LA.solve_LS_fast1(A, z);
    }
    catch (NumericalException e)
    {
      A = q.Jx_iM_JxT;
      _LA.solve_LS_fast2(A, z);
    }

    // prepare to compute the nullspace
    A.resize(q.N_CONSTRAINT_EQNS_IMP, q.N_VARS);
    unsigned col_start = 0, col_end = q.N_CONTACTS;
    const unsigned ROW_START = 0, ROW_END = q.N_CONSTRAINT_EQNS_IMP;
    SharedMatrixNd Cn_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_CONTACTS;
    SharedMatrixNd Cs_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_CONTACTS;
    SharedMatrixNd Ct_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_CONTACTS;
    SharedMatrixNd NCs_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_CONTACTS;
    SharedMatrixNd NCt_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_LIMITS;
    SharedMatrixNd L_block = A.block(ROW_START, ROW_END, col_start, col_end);
    col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
    MatrixNd::transpose(q.Cn_iM_JxT, Cn_block);
    MatrixNd::transpose(q.Cs_iM_JxT, Cs_block);
    MatrixNd::transpose(q.Cs_iM_JxT, Ct_block);
    (NCs_block = Cs_block).negate();
    (NCt_block = Ct_block).negate();
    MatrixNd::transpose(q.L_iM_JxT, L_block);

    // compute the nullspace
    _LA.nullspace(A, R);
  }
  else
  {
    R.set_identity(q.N_VARS);
    z.set_zero(q.N_VARS);
  }

  // get number of qp variables
  const unsigned N_VARS = R.columns();

  // init the QP matrix and vector
  const unsigned N_INEQUAL = q.N_CONTACTS + q.N_K_TOTAL + q.N_LIMITS + 1;
  H.resize(q.N_VARS, q.N_VARS);
  c.resize(H.rows());
  A.set_zero(N_INEQUAL, q.N_VARS);
  MM.set_zero(N_VARS + N_INEQUAL, N_VARS + N_INEQUAL);
  qq.resize(MM.rows());
  SharedVectorNd nb = qq.segment(N_VARS, qq.size()).set_zero();
  assert(nb.size() == N_INEQUAL);

  // setup row (block) 1 -- Cn * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  unsigned col_start = 0, col_end = q.N_CONTACTS;
  unsigned row_start = 0, row_end = q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
  SharedMatrixNd Cn_iM_JxT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cn_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONTACTS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
  SharedMatrixNd Cs_iM_JxT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONTACTS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
  SharedMatrixNd Ct_iM_JxT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Ct_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 4 -- -Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONTACTS;
  SharedMatrixNd NCs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 5 -- -Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONTACTS;
  SharedMatrixNd NCt_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block 6) -- L * iM *  [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_LIMITS;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd L_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
  SharedMatrixNd L_iM_JxT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block 5) -- Jx * iM *  [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += q.N_CONSTRAINT_EQNS_IMP;
  col_start = 0; col_end = q.N_CONTACTS;
  SharedMatrixNd Jx_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Jx_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Jx_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Jx_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONTACTS;
  SharedMatrixNd Jx_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_LIMITS;
  SharedMatrixNd Jx_iM_LT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += q.N_CONSTRAINT_EQNS_IMP;
  SharedMatrixNd Jx_iM_JxT = H.block(row_start, row_end, col_start, col_end);

  // copy to row block 1 (contact normals)
  Cn_iM_CnT = q.Cn_iM_CnT;
  Cn_iM_CsT = q.Cn_iM_CsT;
  Cn_iM_CtT = q.Cn_iM_CtT;
  (Cn_iM_NCsT = q.Cn_iM_CsT).negate();
  (Cn_iM_NCtT = q.Cn_iM_CtT).negate();
  Cn_iM_LT = q.Cn_iM_LT;
  Cn_iM_JxT = q.Cn_iM_JxT;
  
  // copy to row block 2 (first contact tangents)
  MatrixNd::transpose(q.Cn_iM_CsT, Cs_iM_CnT);
  Cs_iM_CsT = q.Cs_iM_CsT;
  Cs_iM_CtT = q.Cs_iM_CtT;
  (Cs_iM_NCsT = Cs_iM_CsT).negate();
  (Cs_iM_NCtT = Cs_iM_CtT).negate();
  Cs_iM_LT = q.Cs_iM_LT;
  Cs_iM_JxT = q.Cs_iM_JxT;

  // copy to row block 3 (second contact tangents)
  MatrixNd::transpose(q.Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(q.Cs_iM_CtT, Ct_iM_CsT);
  Ct_iM_CtT = q.Ct_iM_CtT;
  (Ct_iM_NCsT = Ct_iM_CsT).negate();
  (Ct_iM_NCtT = Ct_iM_CtT).negate();
  Ct_iM_LT = q.Ct_iM_LT;
  Ct_iM_JxT = q.Ct_iM_JxT;

  // copy to row block 4 (negative contact tangents)
  (NCs_block = Cs_block).negate();

  // copy to row block 5 (negative contact tangents)
  (NCt_block = Ct_block).negate();

  // copy to row block 6 (limits)
  MatrixNd::transpose(q.Cn_iM_LT, L_iM_CnT);
  MatrixNd::transpose(q.Cs_iM_LT, L_iM_CsT);
  MatrixNd::transpose(q.Ct_iM_LT, L_iM_CtT);
  (L_iM_NCsT = L_iM_CsT).negate();
  (L_iM_NCtT = L_iM_CtT).negate();
  L_iM_LT = q.L_iM_LT;
  L_iM_JxT = q.L_iM_JxT;

  // copy to row block 6 (explicit constraints)
  MatrixNd::transpose(q.Cn_iM_JxT, Jx_iM_CnT);
  MatrixNd::transpose(q.Cs_iM_JxT, Jx_iM_CsT);
  MatrixNd::transpose(q.Ct_iM_JxT, Jx_iM_CtT);
  MatrixNd::transpose(q.L_iM_JxT, Jx_iM_LT);
  (Jx_iM_NCsT = Jx_iM_CsT).negate();
  (Jx_iM_NCtT = Jx_iM_CtT).negate();
  Jx_iM_JxT = q.Jx_iM_JxT;

  // setup c 
  c.set_sub_vec(q.CN_IDX, q.Cn_v);         
  c.set_sub_vec(q.CS_IDX, q.Cs_v);         
  c.set_sub_vec(q.CT_IDX, q.Ct_v);         
  (c.segment(q.NCS_IDX, q.NCT_IDX) = q.Cs_v).negate();
  (c.segment(q.NCT_IDX, q.CS_U_IDX) = q.Ct_v).negate();
  c.set_sub_vec(q.L_IDX, q.L_v);         
  c.set_sub_vec(q.ALPHA_X_IDX, q.Jx_v);         

  // ------- setup A/-b -------
  // setup the Cn*v+ >= 0 constraint
  // Cn*(inv(M)*impulses + v) >= 0, Cn*inv(M)*impulses >= -Cn*v
  row_start = 0; row_end = q.N_CONTACTS;
  A.block(row_start, row_end, q.CN_IDX, q.N_VARS) = Cn_block;
  nb.set_sub_vec(row_start, q.Cn_v);
  row_start = row_end; row_end += q.N_LIMITS;  

  // setup the L*v+ >= 0 constraint
  A.block(row_start, row_end, q.CN_IDX, q.N_VARS) = L_block;
  nb.set_sub_vec(row_start, q.L_v);
  row_start = row_end; row_end += q.N_CONTACTS;  

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // initialize the contact velocity
    double vel = std::sqrt(sqr(q.Cs_v[i]) + sqr(q.Ct_v[i]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK/2; j++)
    {
      double theta = (double) j/(q.contact_events[i]->contact_NK/2-1) * M_PI_2;
      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      A(row_start, q.CN_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row_start, q.CS_IDX+i) = -ct;
      A(row_start, q.NCS_IDX+i) = -ct;
      A(row_start, q.CT_IDX+i) = -st;
      A(row_start, q.NCT_IDX+i) = -st;

      // setup the viscous friction component
      nb[row_start] = q.contact_events[i]->contact_mu_viscous * vel;
      row_start++;
    }
  }

  // setup the normal velocity constraint
  // 1'N*v+ <= 1'N*v-  (in our form) -1'N*v+ >= -1'N*v- 
  for (unsigned i=0; i< Cn_block.columns(); i++)
  {
    SharedConstVectorNd Cn_col = Cn_block.column(i);
    A(row_start, i) = -std::accumulate(Cn_col.begin(), Cn_col.end(), 0.0);
  }
  nb[row_start] = -std::accumulate(q.Cn_v.row_iterator_begin(), q.Cn_v.row_iterator_end(), 0.0) + q.kappa;

  // get useful blocks of MM and segments of qq
  SharedMatrixNd H_block = MM.block(0, N_VARS, 0, N_VARS);
  SharedMatrixNd AR_block = MM.block(N_VARS, MM.rows(), 0, N_VARS);
  SharedMatrixNd ART_block = MM.block(0, N_VARS, N_VARS, MM.rows());
  SharedVectorNd c_seg = qq.segment(0, N_VARS);

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H_block);
  R.transpose_mult(c, c_seg);
  RTH.mult(z, workv);
  c_seg += workv;

  // setup constraints A*x >= b in nullspace, yielding
  // A(R*y + z) >= b, yielding R*y >= b - A*z
  A.mult(z, workv);
  nb += workv;
  A.mult(R, AR_block);
  MatrixNd::transpose(AR_block, ART_block);
  ART_block.negate();

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_general() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cs': " << std::endl << q.Cn_iM_CsT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Ct': " << std::endl << q.Cn_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Jx': " << std::endl << q.Cn_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * Cs': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * Ct': " << std::endl << q.Cs_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * L': " << std::endl << q.Cs_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Cs * inv(M) * Jx': " << std::endl << q.Cs_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Ct * inv(M) * Ct': " << std::endl << q.Ct_iM_CtT;
  FILE_LOG(LOG_EVENT) << "  Ct * inv(M) * L': " << std::endl << q.Ct_iM_LT;
  FILE_LOG(LOG_EVENT) << "  Ct * inv(M) * Jx': " << std::endl << q.Ct_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  L * inv(M) * L': " << std::endl << q.L_iM_LT;
  FILE_LOG(LOG_EVENT) << "  L * inv(M) * Jx': " << std::endl << q.L_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Cs * v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Ct * v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "-b vector: " << nb << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!_lcp.lcp_lemke_regularized(MM, qq, workv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_EVENT) << "LCP solution: " << workv << std::endl; 
  workv.get_sub_vec(0, N_VARS, y);
  R.mult(y, workv);
  z += workv;

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_general() exited" << std::endl;
}

