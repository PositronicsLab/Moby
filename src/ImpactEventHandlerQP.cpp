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
#include <Moby/LinAlg.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/Optimization.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactEventHandler.h>

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
void ImpactEventHandler::solve_qp(EventProblemData& q, Real poisson_eps)
{
  SAFESTATIC VectorN z, tmp, tmp2;
  const Real TOL = poisson_eps;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_CONSTRAINT_DOF_IMP = q.N_CONSTRAINT_DOF_IMP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // solve the QP
  solve_qp_work(q, z);

  // apply (Poisson) restitution to contacts
  for (unsigned i=0; i< N_CONTACTS; i++)
    z[i] *= ((Real) 1.0 + q.contact_events[i]->contact_epsilon);

  // apply (Poisson) restitution to limits
  for (unsigned i=0; i< N_LIMITS; i++)
    z[N_CONTACTS*5+i] *= ((Real) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  update_impulses(q, z);

  // update Jc_v, Dc_v, Jl_v, and Jx_v
  q.Jc_v += q.Jc_iM_JcT.mult(q.alpha_c, tmp);
  q.Jc_v += q.Jc_iM_DcT.mult(q.beta_c, tmp);
  q.Jc_v += q.Jc_iM_JlT.mult(q.alpha_l, tmp);
  q.Jc_v += q.Jc_iM_JxT.mult(q.alpha_x, tmp);
  q.Dc_v += q.Jc_iM_DcT.transpose_mult(q.alpha_c, tmp);
  q.Dc_v += q.Dc_iM_DcT.mult(q.beta_c, tmp);
  q.Dc_v += q.Dc_iM_JlT.mult(q.alpha_l, tmp);
  q.Dc_v += q.Dc_iM_JxT.mult(q.alpha_x, tmp);
  q.Jl_v += q.Jc_iM_JlT.transpose_mult(q.alpha_c, tmp);
  q.Jl_v += q.Dc_iM_JlT.transpose_mult(q.beta_c, tmp);
  q.Jl_v += q.Jl_iM_JlT.mult(q.alpha_l, tmp);
  q.Jl_v += q.Jl_iM_JxT.mult(q.alpha_x, tmp);
  q.Jx_v += q.Jc_iM_JxT.transpose_mult(q.alpha_c, tmp);
  q.Jx_v += q.Dc_iM_JxT.transpose_mult(q.beta_c, tmp);
  q.Jx_v += q.Jl_iM_JxT.transpose_mult(q.alpha_l, tmp);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, tmp);

  // output results
  FILE_LOG(LOG_EVENT) << "results: " << std::endl;
  FILE_LOG(LOG_EVENT) << "new Jc_v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Dc_v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Jl_v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "new Jx_v: " << q.Jx_v << std::endl;

  // see whether another QP must be solved
  if (q.Jc_v.size() > 0 && *min_element(q.Jc_v.begin(), q.Jc_v.end()) < -TOL)
  {
    FILE_LOG(LOG_EVENT) << "minimum Jc*v: " << *min_element(q.Jc_v.begin(), q.Jc_v.end()) << std::endl;
    FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
    solve_qp_work(q, z);
    update_impulses(q, z);
  }
  else 
    if (q.Jl_v.size() > 0 && *min_element(q.Jl_v.begin(), q.Jl_v.end()) < -TOL)
    {
      FILE_LOG(LOG_EVENT) << "minimum Jl*v: " << *min_element(q.Jl_v.begin(), q.Jl_v.end()) << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      update_impulses(q, z);
    }
  else
  {
    pair<Real*, Real*> mm = boost::minmax_element(q.Jx_v.begin(), q.Jx_v.end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_EVENT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_EVENT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      update_impulses(q, z);
    }
  }

  // save normal contact impulses
  for (unsigned i=0; i< N_CONTACTS; i++)
    q.contact_events[i]->contact_impulse = q.contact_events[i]->contact_normal * q.alpha_c[i];

  // save tangent contact impulses
  for (unsigned i=0, j=0; i< N_CONTACTS; i++)
  {
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan1 * q.beta_c[j++];
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan2 * q.beta_c[j++];
  }

  // save limit impulses
  for (unsigned i=0; i< N_LIMITS; i++)
    q.limit_events[i]->limit_impulse = q.alpha_l[i]; 
}

/// Solves the quadratic program (does all of the work) for given active contact set
/**
 * \note this is the version without joint friction forces
 */
/*
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixN sub, t1, t2, t3, neg1, A, AR, R, RTH;
  SAFESTATIC MatrixN H, MM;
  SAFESTATIC VectorN negv, c, qq, nb, tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned ALPHA_L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned NVARS = N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;

  // first, solve for impulses that satisfy explicit constraint equations
  // and compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_EXP > 0)
  {
    // compute the homogeneous solution
    A.copy_from(q.Jx_iM_JxT);
    z.copy_from(q.Jx_v).negate();
    try
    {
      LinAlg::solve_LS_fast1(A, z);
    }
    catch (NumericalException e)
    {
      A.copy_from(q.Jx_iM_JxT);
      LinAlg::solve_LS_fast2(A, z);
    }

    // compute the nullspace
    A.resize(N_CONSTRAINT_EQNS_EXP, NVARS);
    MatrixN::transpose(q.Jc_iM_JxT, t1);
    MatrixN::transpose(q.Dc_iM_JxT, t2);
    MatrixN::transpose(q.Jl_iM_JxT, t3);
    neg1.copy_from(t2).negate();
    A.set_sub_row_block(0, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);
    LinAlg::nullspace(A, R);
  }
  else
  {
    R.set_zero(NVARS,NVARS);
    for (unsigned i=0; i< NVARS; i++) 
      R(i,i) = (Real) 1.0;
    z.set_zero(NVARS);
  }

  // get number of qp variables
  const unsigned N_PRIMAL = R.columns();

  // init the QP matrix and vector
  const unsigned KAPPA = (q.use_kappa) ? 1 : 0;
  const unsigned N_INEQUAL = N_CONTACTS + N_K_TOTAL + N_LIMITS + KAPPA;
  H.resize(N_PRIMAL, N_PRIMAL);
  c.resize(H.rows());
  A.set_zero(N_INEQUAL, N_PRIMAL);
  nb.set_zero(N_INEQUAL);
  MM.set_zero(N_PRIMAL + N_INEQUAL, N_PRIMAL + N_INEQUAL);
  qq.resize(MM.rows());

  // setup [Q M'; -M 0]
  unsigned col = 0, row = 0;

  // row (block) 1 -- Jc * iM * [Jc' Dc' -Dc' Jl' Jx']
  neg1.copy_from(q.Jc_iM_DcT).negate();
  H.set_sub_row_block(0, &q.Jc_iM_JcT, &q.Jc_iM_DcT, &neg1, &q.Jc_iM_JlT, 
                      &q.Jc_iM_JxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Jc' Dc' -Dc' Jl' Jx']
  MatrixN::transpose(q.Jc_iM_DcT, t1);
  neg1.copy_from(q.Dc_iM_DcT).negate();
  H.set_sub_row_block(row, &t1, &q.Dc_iM_DcT, &neg1, &q.Dc_iM_JlT, 
                      &q.Dc_iM_JxT);

  // row (block 3) -- negated block 2
  H.get_sub_mat(row, row+N_CONTACTS*2, 0, H.columns(), sub);
  H.set_sub_mat(row+N_CONTACTS*2, 0, sub.negate());
  row += N_CONTACTS*4;

  // row (block 4) -- Jl * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JlT, t1);
  MatrixN::transpose(q.Dc_iM_JlT, t2);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &q.Jl_iM_JlT, &q.Jl_iM_JxT);
  row += N_LIMITS;
  
  // row (block 5) -- Jx * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JxT, t1);
  MatrixN::transpose(q.Dc_iM_JxT, t2);
  MatrixN::transpose(q.Jl_iM_JxT, t3);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);

  // setup c 
  c.set_sub_vec(ALPHA_C_IDX, q.Jc_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  negv.copy_from(q.Dc_v).negate();
  c.set_sub_vec(NBETA_C_IDX, negv);           
  c.set_sub_vec(ALPHA_L_IDX, q.Jl_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // setup the Jc*v+ >= 0 constraint
  // Jc*(inv(M)*impulses + v) >= 0, Jc*inv(M)*impulses >= -Jc*v
  row = 0; col = 0;
  nb.set_sub_vec(row, q.Jc_v);
  H.get_sub_mat(ALPHA_C_IDX, ALPHA_C_IDX+N_CONTACTS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the Jl*v+ >= 0 constraint
  nb.set_sub_vec(row, q.Jl_v);
  H.get_sub_mat(ALPHA_L_IDX, ALPHA_L_IDX+N_LIMITS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, k=0; i< N_CONTACTS; i++, k+= 2)
  {
    // initialize the contact velocity
    Real vel = std::sqrt(sqr(q.Dc_v[k]) + sqr(q.Dc_v[k+1]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK; j++)
    {
      Real theta = (Real) j/(q.contact_events[i]->contact_NK-1) * M_PI_2;
      const Real ct = std::cos(theta);
      const Real st = std::sin(theta);
      A(row, ALPHA_C_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row, BETA_C_IDX+k) = -ct;
      A(row, NBETA_C_IDX+k) = -ct;
      A(row, BETA_C_IDX+k+1) = -st;
      A(row, NBETA_C_IDX+k+1) = -st;

      // setup the viscous friction component
      nb[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }
  }

  // setup the normal impulse constraint
  // 1'cn <= kappa (equiv. to -1'cn >= -kappa)
  if (q.use_kappa)
  {
    for (unsigned i=0; i< N_CONTACTS; i++)
      A(row, ALPHA_C_IDX+i) = (Real) -1.0;
    nb[row] = q.kappa;
  }

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H);
  R.transpose_mult(c, tmpv);
  c.copy_from(tmpv);
  RTH.mult(z, tmpv);
  c += tmpv;

  // setup constraints A*x >= b in nullspace, yielding
  // A(R*y + z) >= b, yielding R*y >= b - A*z
  A.mult(z, tmpv);
  nb += tmpv;
  A.mult(R, AR);

  // setup the LCP matrix
  MM.set_sub_mat(0, 0, H);
  MM.set_sub_mat(N_PRIMAL, 0, AR);
  MM.set_sub_mat(0, N_PRIMAL, AR.negate(), true);

  // setup the LCP vector
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N_PRIMAL, nb);

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!Optimization::lcp_lemke_regularized(MM, qq, tmpv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_EVENT) << "LCP solution: " << tmpv << std::endl; 
  tmpv.get_sub_vec(0, N_PRIMAL, y);
  R.mult(y, tmpv);
  z += tmpv;

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() exited" << std::endl;
}
*/

void ImpactEventHandler::update_solution(const EventProblemData& q, const VectorN& x, const vector<bool>& working_set, unsigned jidx, VectorN& z)
{
  SAFESTATIC VectorN workv;

  // setup NEW constants
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned ALPHA_L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned N_VARS = N_LIMITS + ALPHA_L_IDX;
  const unsigned NP_IDX = N_VARS;
  const unsigned MU_IDX = NP_IDX + N_LIMITS; 

  // setup old constants
  const unsigned OLD_N_CONTACTS = N_CONTACTS - 1;
  const unsigned OLD_ALPHA_C_IDX = 0;
  const unsigned OLD_BETA_C_IDX = OLD_N_CONTACTS;
  const unsigned OLD_NBETA_C_IDX = OLD_N_CONTACTS*2 + OLD_BETA_C_IDX;
  const unsigned OLD_ALPHA_L_IDX = OLD_N_CONTACTS*2 + OLD_NBETA_C_IDX;
  const unsigned OLD_NVARS = N_LIMITS + OLD_ALPHA_L_IDX;
  const unsigned OLD_NP_IDX = OLD_NVARS;
  const unsigned OLD_MU_IDX = OLD_NP_IDX + N_LIMITS; 

  // determine NK indices 
  unsigned OLD_MU_JIDX = 0;
  unsigned J_NK = 0;
  for (unsigned i=0, j=0; i< N_CONTACTS; i++)
  {
    if (!working_set[i])
      continue;
    if (j++ == jidx)
    {
      J_NK += q.contact_events[i]->contact_NK; 
      break;
    }
    OLD_MU_JIDX += q.contact_events[i]->contact_NK;
  }

  // initialize z
  z.set_zero(x.size() + 6 + J_NK);

  // setup normal contact variables
  x.get_sub_vec(OLD_ALPHA_C_IDX, OLD_ALPHA_C_IDX+jidx, workv);
  z.set_sub_vec(ALPHA_C_IDX, workv);
  x.get_sub_vec(OLD_ALPHA_C_IDX+jidx, OLD_BETA_C_IDX, workv);
  z.set_sub_vec(ALPHA_C_IDX+jidx+1, workv);

  // setup positive tangent contact variables
  x.get_sub_vec(OLD_BETA_C_IDX, OLD_BETA_C_IDX+jidx*2, workv);
  z.set_sub_vec(BETA_C_IDX, workv);
  x.get_sub_vec(OLD_BETA_C_IDX+jidx*2, OLD_NBETA_C_IDX, workv);
  z.set_sub_vec(BETA_C_IDX+(jidx+1)*2, workv);

  // setup negative tangent contact variables and all other primal variables 
  x.get_sub_vec(OLD_NBETA_C_IDX, OLD_NBETA_C_IDX+jidx*2, workv);
  z.set_sub_vec(NBETA_C_IDX, workv);
  x.get_sub_vec(OLD_NBETA_C_IDX+jidx*2, OLD_NVARS, workv);
  z.set_sub_vec(NBETA_C_IDX+(jidx+1)*2, workv);

  // setup constraint equation variables
  // constraint ordering: noninterpenetration, joint limits, 
  // contact friction
  // 1. noninterpenetration and joint limit constraints
  x.get_sub_vec(OLD_NP_IDX, OLD_NP_IDX+jidx, workv);
  z.set_sub_vec(NP_IDX, workv);
  x.get_sub_vec(OLD_NP_IDX+jidx, OLD_MU_IDX, workv);
  z.set_sub_vec(NP_IDX+jidx+1, workv);

  // 2. contact friction constraints
  x.get_sub_vec(OLD_MU_IDX, OLD_MU_JIDX, workv);
  z.set_sub_vec(MU_IDX, workv);
  x.get_sub_vec(OLD_MU_JIDX, x.size(), workv);
  z.set_sub_vec(MU_IDX+J_NK, workv);
}

/// Checks whether the optimization is satisfied *without* adding contact j
bool ImpactEventHandler::opt_satisfied(const EventProblemData& q, const vector<bool>& working_set, Real& KE, VectorN& x, unsigned j)
{
  SAFESTATIC VectorN workv, Jc_v, z;
  SAFESTATIC EventProblemData qcopy;

  // make a copy of q and set it appropriately
  qcopy.copy_from(q);
  vector<bool> ws = working_set;
  ws[j] = true;
  const unsigned JIDX = std::count(ws.begin(), ws.begin()+j, true);
  update_problem(q, ws, qcopy);

  // copy x to appropriate parts of z to facilitate warm-starting
  // setup variable indices
  update_solution(qcopy, x, working_set, JIDX, z);   

  // solve the problem with addition of contact j
  solve_qp_work_ijoints(qcopy, z); 

  // check whether normal force for contact j is zero
  if (z[JIDX] < NEAR_ZERO)
    return true; 

  // check whether there is an appreciable change in K.E.
  Real new_KE = calc_ke(qcopy);
  if (new_KE < KE - NEAR_ZERO)
  {
    x.copy_from(z);
    KE = new_KE;
    return false;
  }
  else
    return true;   
}

/// Computes the kinetic energy of the system using the current impulse set 
Real ImpactEventHandler::calc_ke(const EventProblemData& q)
{
  Real KE = (Real) 0.0;
  set_generalized_velocities(q);
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    KE += q.super_bodies[i]->calc_kinetic_energy();
  return KE;
}

/// Solves the quadratic program (does all of the work) 
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorN& z)
{
  SAFESTATIC vector<bool> working_set;
  SAFESTATIC EventProblemData qworking;

  // if there are explicit constraints, cannot use the fast method
  if (q.N_CONSTRAINT_EQNS_EXP > 0)
  {
    solve_qp_work_general(q, z);
    return;
  } 

  // if we're not dealing with many contacts, exit now 
  if (q.N_CONTACTS < 4)
  {
    solve_qp_work_ijoints(q, z);
    return;
  }

  // setup the working set -- set only first contact to active
  working_set.resize(q.N_CONTACTS);
  std::fill(working_set.begin(), working_set.end(), false);
  working_set[0] = true;
  update_problem(q, working_set, qworking);

  // solve using warm starting
  solve_qp_work_ijoints(qworking, z);
  update_impulses(qworking, z);
  Real KE = calc_ke(qworking);

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

  // copy qworking to q
  q.copy_from(qworking);

/*
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!Optimization::lcp_lemke_regularized(MM, qq, tmpv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_EVENT) << "LCP solution: " << tmpv << std::endl; 
  tmpv.get_sub_vec(0, N_PRIMAL, y);
  R.mult(y, tmpv);
  z += tmpv;
*/

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() exited" << std::endl;
}

/// Updates problem matrices and vectors based on the working set
void ImpactEventHandler::update_problem(const EventProblemData& q, const vector<bool>& working_set, EventProblemData& qworking)
{
  SAFESTATIC vector<unsigned> norm_indices, tan_indices;

  // determine how many contacts we have
  qworking.N_CONTACTS = std::count(working_set.begin(), working_set.end(), true);

  // setup contact indices set
  norm_indices.clear();
  for (unsigned i=0; i< working_set.size(); i++)
    if (working_set[i])
    {
      norm_indices.push_back(i);
      tan_indices.push_back(i*2);
      tan_indices.push_back(i*2+1);
    }

  // select appropriate parts of vectors 
  q.Jc_v.select(norm_indices.begin(), norm_indices.end(), qworking.Jc_v);
  q.Dc_v.select(tan_indices.begin(), tan_indices.end(), qworking.Dc_v);

  // select appropriate parts of Jc matrices 
  q.Jc_iM_JcT.select_square(norm_indices.begin(), norm_indices.end(), qworking.Jc_iM_JcT);
  q.Jc_iM_DcT.select(norm_indices.begin(), norm_indices.end(), tan_indices.begin(), tan_indices.end(), qworking.Jc_iM_DcT);
  q.Jc_iM_JlT.select_rows(norm_indices.begin(), norm_indices.end(), qworking.Jc_iM_JlT);
  q.Jc_iM_JxT.select_rows(norm_indices.begin(), norm_indices.end(), qworking.Jc_iM_JxT);
  q.Jc_iM_DxT.select_rows(norm_indices.begin(), norm_indices.end(), qworking.Jc_iM_DxT);

  // select appropriate parts of Dc matrices
  q.Dc_iM_DcT.select_square(tan_indices.begin(), tan_indices.end(), qworking.Dc_iM_DcT);
  q.Dc_iM_JlT.select_rows(tan_indices.begin(), tan_indices.end(), qworking.Dc_iM_JlT);
  q.Dc_iM_JxT.select_rows(tan_indices.begin(), tan_indices.end(), qworking.Dc_iM_JxT);
  q.Dc_iM_DxT.select_rows(tan_indices.begin(), tan_indices.end(), qworking.Dc_iM_DxT);
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work_ijoints(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixN sub, t1, t2, t3, neg1, A; 
  SAFESTATIC MatrixN H, MM;
  SAFESTATIC VectorN negv, c, qq, nb, tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned ALPHA_L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned NVARS = N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;

  // get number of qp variables
  const unsigned N_PRIMAL = NVARS;

  // init the QP matrix and vector
  const unsigned KAPPA = (q.use_kappa) ? 1 : 0;
  const unsigned N_INEQUAL = N_CONTACTS + N_K_TOTAL + N_LIMITS + KAPPA;
  H.resize(N_PRIMAL, N_PRIMAL);
  c.resize(H.rows());
  A.set_zero(N_INEQUAL, N_PRIMAL);
  nb.set_zero(N_INEQUAL);
  MM.set_zero(N_PRIMAL + N_INEQUAL, N_PRIMAL + N_INEQUAL);
  qq.resize(MM.rows());

  // setup [Q M'; -M 0]
  unsigned col = 0, row = 0;

  // row (block) 1 -- Jc * iM * [Jc' Dc' -Dc' Jl' Jx']
  neg1.copy_from(q.Jc_iM_DcT).negate();
  H.set_sub_row_block(0, &q.Jc_iM_JcT, &q.Jc_iM_DcT, &neg1, &q.Jc_iM_JlT, 
                      &q.Jc_iM_JxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Jc' Dc' -Dc' Jl' Jx']
  MatrixN::transpose(q.Jc_iM_DcT, t1);
  neg1.copy_from(q.Dc_iM_DcT).negate();
  H.set_sub_row_block(row, &t1, &q.Dc_iM_DcT, &neg1, &q.Dc_iM_JlT, 
                      &q.Dc_iM_JxT);

  // row (block 3) -- negated block 2
  H.get_sub_mat(row, row+N_CONTACTS*2, 0, H.columns(), sub);
  H.set_sub_mat(row+N_CONTACTS*2, 0, sub.negate());
  row += N_CONTACTS*4;

  // row (block 4) -- Jl * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JlT, t1);
  MatrixN::transpose(q.Dc_iM_JlT, t2);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &q.Jl_iM_JlT, &q.Jl_iM_JxT);
  row += N_LIMITS;
  
  // row (block 5) -- Jx * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JxT, t1);
  MatrixN::transpose(q.Dc_iM_JxT, t2);
  MatrixN::transpose(q.Jl_iM_JxT, t3);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);

  // setup c 
  c.set_sub_vec(ALPHA_C_IDX, q.Jc_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  negv.copy_from(q.Dc_v).negate();
  c.set_sub_vec(NBETA_C_IDX, negv);           
  c.set_sub_vec(ALPHA_L_IDX, q.Jl_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // setup the Jc*v+ >= 0 constraint
  // Jc*(inv(M)*impulses + v) >= 0, Jc*inv(M)*impulses >= -Jc*v
  row = 0; col = 0;
  nb.set_sub_vec(row, q.Jc_v);
  H.get_sub_mat(ALPHA_C_IDX, ALPHA_C_IDX+N_CONTACTS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the Jl*v+ >= 0 constraint
  nb.set_sub_vec(row, q.Jl_v);
  H.get_sub_mat(ALPHA_L_IDX, ALPHA_L_IDX+N_LIMITS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, k=0; i< N_CONTACTS; i++, k+= 2)
  {
    // initialize the contact velocity
    Real vel = std::sqrt(sqr(q.Dc_v[k]) + sqr(q.Dc_v[k+1]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK; j++)
    {
      Real theta = (Real) j/(q.contact_events[i]->contact_NK-1) * M_PI_2;
      const Real ct = std::cos(theta);
      const Real st = std::sin(theta);
      A(row, ALPHA_C_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row, BETA_C_IDX+k) = -ct;
      A(row, NBETA_C_IDX+k) = -ct;
      A(row, BETA_C_IDX+k+1) = -st;
      A(row, NBETA_C_IDX+k+1) = -st;

      // setup the viscous friction component
      nb[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }
  }

  // setup the normal impulse constraint
  // 1'cn <= kappa (equiv. to -1'cn >= -kappa)
  if (q.use_kappa)
  {
    for (unsigned i=0; i< N_CONTACTS; i++)
      A(row, ALPHA_C_IDX+i) = (Real) -1.0;
    nb[row] = q.kappa;
  }

  // setup the LCP matrix
  MM.set_sub_mat(0, 0, H);
  MM.set_sub_mat(N_PRIMAL, 0, A);
  MM.set_sub_mat(0, N_PRIMAL, A.negate(), true);

  // setup the LCP vector
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N_PRIMAL, nb);

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!Optimization::lcp_lemke_regularized(MM, qq, z))
    throw std::runtime_error("Unable to solve event QP!");

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_ijoints() exited" << std::endl;
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work_general(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixN sub, t1, t2, t3, neg1, A, AR, R, RTH;
  SAFESTATIC MatrixN H, MM;
  SAFESTATIC VectorN negv, c, qq, nb, tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned ALPHA_L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned NVARS = N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;

  // first, solve for impulses that satisfy explicit constraint equations
  // and compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_EXP > 0)
  {
    // compute the homogeneous solution
    A.copy_from(q.Jx_iM_JxT);
    z.copy_from(q.Jx_v).negate();
    try
    {
      LinAlg::solve_LS_fast1(A, z);
    }
    catch (NumericalException e)
    {
      A.copy_from(q.Jx_iM_JxT);
      LinAlg::solve_LS_fast2(A, z);
    }

    // compute the nullspace
    A.resize(N_CONSTRAINT_EQNS_EXP, NVARS);
    MatrixN::transpose(q.Jc_iM_JxT, t1);
    MatrixN::transpose(q.Dc_iM_JxT, t2);
    MatrixN::transpose(q.Jl_iM_JxT, t3);
    neg1.copy_from(t2).negate();
    A.set_sub_row_block(0, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);
    LinAlg::nullspace(A, R);
  }
  else
  {
    R.set_zero(NVARS,NVARS);
    for (unsigned i=0; i< NVARS; i++) 
      R(i,i) = (Real) 1.0;
    z.set_zero(NVARS);
  }

  // get number of qp variables
  const unsigned N_PRIMAL = R.columns();

  // init the QP matrix and vector
  const unsigned KAPPA = (q.use_kappa) ? 1 : 0;
  const unsigned N_INEQUAL = N_CONTACTS + N_K_TOTAL + N_LIMITS + KAPPA;
  H.resize(N_PRIMAL, N_PRIMAL);
  c.resize(H.rows());
  A.set_zero(N_INEQUAL, N_PRIMAL);
  nb.set_zero(N_INEQUAL);
  MM.set_zero(N_PRIMAL + N_INEQUAL, N_PRIMAL + N_INEQUAL);
  qq.resize(MM.rows());

  // setup [Q M'; -M 0]
  unsigned col = 0, row = 0;

  // row (block) 1 -- Jc * iM * [Jc' Dc' -Dc' Jl' Jx']
  neg1.copy_from(q.Jc_iM_DcT).negate();
  H.set_sub_row_block(0, &q.Jc_iM_JcT, &q.Jc_iM_DcT, &neg1, &q.Jc_iM_JlT, 
                      &q.Jc_iM_JxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Jc' Dc' -Dc' Jl' Jx']
  MatrixN::transpose(q.Jc_iM_DcT, t1);
  neg1.copy_from(q.Dc_iM_DcT).negate();
  H.set_sub_row_block(row, &t1, &q.Dc_iM_DcT, &neg1, &q.Dc_iM_JlT, 
                      &q.Dc_iM_JxT);

  // row (block 3) -- negated block 2
  H.get_sub_mat(row, row+N_CONTACTS*2, 0, H.columns(), sub);
  H.set_sub_mat(row+N_CONTACTS*2, 0, sub.negate());
  row += N_CONTACTS*4;

  // row (block 4) -- Jl * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JlT, t1);
  MatrixN::transpose(q.Dc_iM_JlT, t2);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &q.Jl_iM_JlT, &q.Jl_iM_JxT);
  row += N_LIMITS;
  
  // row (block 5) -- Jx * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JxT, t1);
  MatrixN::transpose(q.Dc_iM_JxT, t2);
  MatrixN::transpose(q.Jl_iM_JxT, t3);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);

  // setup c 
  c.set_sub_vec(ALPHA_C_IDX, q.Jc_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  negv.copy_from(q.Dc_v).negate();
  c.set_sub_vec(NBETA_C_IDX, negv);           
  c.set_sub_vec(ALPHA_L_IDX, q.Jl_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // setup the Jc*v+ >= 0 constraint
  // Jc*(inv(M)*impulses + v) >= 0, Jc*inv(M)*impulses >= -Jc*v
  row = 0; col = 0;
  nb.set_sub_vec(row, q.Jc_v);
  H.get_sub_mat(ALPHA_C_IDX, ALPHA_C_IDX+N_CONTACTS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the Jl*v+ >= 0 constraint
  nb.set_sub_vec(row, q.Jl_v);
  H.get_sub_mat(ALPHA_L_IDX, ALPHA_L_IDX+N_LIMITS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, k=0; i< N_CONTACTS; i++, k+= 2)
  {
    // initialize the contact velocity
    Real vel = std::sqrt(sqr(q.Dc_v[k]) + sqr(q.Dc_v[k+1]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK; j++)
    {
      Real theta = (Real) j/(q.contact_events[i]->contact_NK-1) * M_PI_2;
      const Real ct = std::cos(theta);
      const Real st = std::sin(theta);
      A(row, ALPHA_C_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row, BETA_C_IDX+k) = -ct;
      A(row, NBETA_C_IDX+k) = -ct;
      A(row, BETA_C_IDX+k+1) = -st;
      A(row, NBETA_C_IDX+k+1) = -st;

      // setup the viscous friction component
      nb[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }
  }

  // setup the normal impulse constraint
  // 1'cn <= kappa (equiv. to -1'cn >= -kappa)
  if (q.use_kappa)
  {
    for (unsigned i=0; i< N_CONTACTS; i++)
      A(row, ALPHA_C_IDX+i) = (Real) -1.0;
    nb[row] = q.kappa;
  }

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H);
  R.transpose_mult(c, tmpv);
  c.copy_from(tmpv);
  RTH.mult(z, tmpv);
  c += tmpv;

  // setup constraints A*x >= b in nullspace, yielding
  // A(R*y + z) >= b, yielding R*y >= b - A*z
  A.mult(z, tmpv);
  nb += tmpv;
  A.mult(R, AR);

  // setup the LCP matrix
  MM.set_sub_mat(0, 0, H);
  MM.set_sub_mat(N_PRIMAL, 0, AR);
  MM.set_sub_mat(0, N_PRIMAL, AR.negate(), true);

  // setup the LCP vector
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N_PRIMAL, nb);

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_general() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_EVENT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_EVENT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_EVENT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_EVENT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_EVENT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_EVENT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_EVENT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_EVENT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!Optimization::lcp_lemke_regularized(MM, qq, tmpv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_EVENT) << "LCP solution: " << tmpv << std::endl; 
  tmpv.get_sub_vec(0, N_PRIMAL, y);
  R.mult(y, tmpv);
  z += tmpv;

  FILE_LOG(LOG_EVENT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_qp_work_general() exited" << std::endl;
}

