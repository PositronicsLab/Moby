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

/// Solves the nonlinearly constrained quadratic program (potentially solves two nQPs, actually)
void ImpactEventHandler::solve_nqp(const VectorNd& zf, EventProblemData& q, double poisson_eps, double max_time)
{
  const double TOL = poisson_eps;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned CL_IDX = N_CONTACTS*3;

  // TODO: set z to frictionless solution to start

  // mark starting time
  tms cstart;
  times(&cstart);

  // keep solving until we run out of time or all contact points are active
  while (true)
  {
    FILE_LOG(LOG_EVENT) << "Running NQP solve iteration with " << (q.N_ACT_CONTACTS) << " active contacts" << std::endl;

    // solve the nonlinearly constrained QP
    solve_nqp_work(q, _z);

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
  for (unsigned i=0, j=q.L_IDX; i< N_LIMITS; i++, j++)
    _z[j] *= ((double) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  q.update_from_stacked_nqp(_z);

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
    FILE_LOG(LOG_EVENT) << " -- running another interior-point iteration..." << std::endl;
    solve_nqp_work(q, _z);
    q.update_from_stacked_nqp(_z);
  }
  else if (q.L_v.size() > 0 && *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) < -TOL)
    {
      FILE_LOG(LOG_EVENT) << "minimum L*v: " << *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()) << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another interior-point iteration..." << std::endl;
      solve_nqp_work(q, _z);
      q.update_from_stacked_nqp(_z);
    }
  else
  {
    pair<ColumnIteratord, ColumnIteratord> mm = boost::minmax_element(q.Jx_v.column_iterator_begin(), q.Jx_v.column_iterator_end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_EVENT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_EVENT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_EVENT) << " -- running another interior-point iteration..." << std::endl;
      solve_nqp_work(q, _z);
      q.update_from_stacked_nqp(_z);
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

/// Solves the nonlinearly constrained quadratic program (does all of the work)
/**
 * \param x the solution is returned here; zeros will be returned at appropriate indices for inactive contacts
 */
void ImpactEventHandler::solve_nqp_work(EventProblemData& q, VectorNd& x)
{
  const double INF = std::numeric_limits<double>::max();

  // setup constants
  const unsigned N_CONTACTS = q.N_ACT_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_IMP = q.N_CONSTRAINT_EQNS_IMP; 
  const unsigned CN_IDX = 0;
  const unsigned CS_IDX = N_CONTACTS;
  const unsigned CT_IDX = CS_IDX + N_CONTACTS;
  const unsigned CL_IDX = CT_IDX + N_CONTACTS;
  const unsigned NVARS = N_LIMITS + CL_IDX; 

  // setup the optimization data
  _ipsolver->epd = &q;
  _ipsolver->mu_c.resize(N_CONTACTS);
  _ipsolver->mu_visc.resize(N_CONTACTS);

  // setup true friction cone for every contact
  for (unsigned i=0; i< N_CONTACTS; i++)
  {
    _ipsolver->mu_c[i] = sqr(q.contact_events[i]->contact_mu_coulomb);
    _ipsolver->mu_visc[i] = (sqr(q.Cs_v[i]) + sqr(q.Ct_v[i])) *
                       sqr(q.contact_events[i]->contact_mu_viscous);
  }

  // setup matrices
  MatrixNd& R = _ipsolver->R;
  MatrixNd& H = _ipsolver->H;
  VectorNd& c = _ipsolver->c; 
  VectorNd& z = _ipsolver->z; 

  // init z (particular solution) 
  z.set_zero(NVARS);

  // first, compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_IMP > 0)
  {
    // compute the homogeneous solution
    _A = q.Jx_iM_JxT;
    (_workv = q.Jx_v).negate();
    try
    {
      _LA.solve_LS_fast1(_A, _workv);
    }
    catch (NumericalException e)
    {
      _A = q.Jx_iM_JxT;
      _LA.solve_LS_fast2(_A, _workv);
    }
    z.set_sub_vec(q.ALPHA_X_IDX, _workv);

    // setup blocks of A
    _A.resize(N_CONSTRAINT_EQNS_IMP, NVARS);
    SharedMatrixNd b1 = _A.block(0, N_CONSTRAINT_EQNS_IMP, 0, N_CONTACTS);
    SharedMatrixNd b2 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_CONTACTS, N_CONTACTS*2);
    SharedMatrixNd b3 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_CONTACTS*2, N_CONTACTS*3);
    SharedMatrixNd b4 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_CONTACTS*3, N_CONTACTS*3+N_LIMITS);

    // compute the nullspace
    MatrixNd::transpose(q.Cn_iM_JxT, b1);
    MatrixNd::transpose(q.Cs_iM_JxT, b2);
    MatrixNd::transpose(q.Ct_iM_JxT, b3);
    MatrixNd::transpose(q.L_iM_JxT, b4);
    _LA.nullspace(_A, R);
  }
  else
    // clear the nullspace 
    R.resize(0,0);

  // get number of qp variables
  const unsigned N_PRIMAL = (R.columns() > 0) ? R.columns() : NVARS;

  // setup number of nonlinear inequality constraints
  const unsigned NONLIN_INEQUAL = N_CONTACTS;

  // init the QP matrix and vector
  H.resize(N_PRIMAL, N_PRIMAL);
  c.resize(H.rows());

  // setup row (block) 1 -- Cn * iM * [Cn' Cs Ct' L']
  unsigned col_start = 0, col_end = N_CONTACTS;
  unsigned row_start = 0, row_end = N_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cn_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_CONTACTS;
  col_start = 0; col_end = N_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_CONTACTS;
  col_start = 0; col_end = N_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block 4) -- L * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_LIMITS;
  col_start = 0; col_end = N_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // copy to row block 1 (contact normals)
  q.Cn_iM_CnT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CnT);
  q.Cn_iM_CsT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CsT);
  q.Cn_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CtT);
  q.Cn_iM_LT.get_sub_mat(0, N_CONTACTS, 0, N_LIMITS, Cn_iM_LT);
  
  // copy to row block 2 (first contact tangents)
  MatrixNd::transpose(Cn_iM_CsT, Cs_iM_CnT);
  q.Cs_iM_CsT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cs_iM_CsT);
  q.Cs_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cs_iM_CtT);
  q.Cs_iM_LT.get_sub_mat(0, N_CONTACTS, 0, N_LIMITS, Cs_iM_LT);

  // copy to row block 3 (second contact tangents)
  MatrixNd::transpose(Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(Cs_iM_CtT, Ct_iM_CsT);
  q.Ct_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Ct_iM_CtT);
  q.Ct_iM_LT.get_sub_mat(0, N_CONTACTS, 0, N_LIMITS, Ct_iM_LT);

  // copy to row block 6 (limits)
  MatrixNd::transpose(Cn_iM_LT, L_iM_CnT);
  MatrixNd::transpose(Cs_iM_LT, L_iM_CsT);
  MatrixNd::transpose(Ct_iM_LT, L_iM_CtT);
  q.L_iM_LT.get_sub_mat(0, N_LIMITS, 0, N_LIMITS, L_iM_LT);

  // get components of c
  SharedVectorNd Cn_v = c.segment(0, N_CONTACTS); 
  SharedVectorNd Cs_v = c.segment(N_CONTACTS, N_CONTACTS*2); 
  SharedVectorNd Ct_v = c.segment(N_CONTACTS*2, N_CONTACTS*3); 
  SharedVectorNd L_v = c.segment(N_CONTACTS*3, N_CONTACTS*3+N_LIMITS); 

  // setup c 
  q.Cn_v.get_sub_vec(0, N_CONTACTS, Cn_v);
  q.Cn_v.get_sub_vec(0, N_CONTACTS, Cs_v);
  q.Cn_v.get_sub_vec(0, N_CONTACTS, Ct_v);
  L_v = q.L_v;

  // ****** now setup linear inequality constraints ******

  // determine whether to use kappa constraint
  const unsigned KAPPA = 1;

  // determine number of linear inequality constraints
  const unsigned N_INEQUAL = N_CONTACTS + N_LIMITS + KAPPA;

  // verify that Cn_block, L_block can be reset
  _ipsolver->Cn_block.reset();
  _ipsolver->L_block.reset();
  _ipsolver->Cn_block = Cn_block;
  _ipsolver->L_block = L_block;

  // setup optimizations in nullspace (if necessary)
  if (R.columns() > 0)
  {
    R.transpose_mult(H, _RTH);
    _RTH.mult(R, H);
    R.transpose_mult(c, _workv);
    c = _workv;
    _RTH.mult(z, _workv);
    c += _workv;
  }

   FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_nqp_work() entered" << std::endl;
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

  // setup ipopt options
  _app.Options()->SetIntegerValue("print_level", 0);
//  _app.Options()->SetStringValue("derivative_test", "second-order");

  // solve the nonlinear QP using the interior-point algorithm 
  Ipopt::ApplicationReturnStatus status = _app.OptimizeTNLP(_ipsolver);
  if (!(status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level))
    throw std::runtime_error("Could not solve nonlinearly constrained QP");

  // get the final solution out
  x = _ipsolver->z;

  FILE_LOG(LOG_EVENT) << "nonlinear QP solution: " << x << std::endl; 
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::solve_nqp() exited" << std::endl;
}


