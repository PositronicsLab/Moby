/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
#include <Moby/UnilateralConstraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactConstraintHandler.h>

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

/// Special functions that do nothing...
#ifndef HAVE_IPOPT
void ImpactConstraintHandler::solve_nqp(VectorNd& z, UnilateralConstraintProblemData& q, double max_time)
{
  throw std::runtime_error("Build without IPOPT!");
}

void ImpactConstraintHandler::solve_nqp_work(UnilateralConstraintProblemData& q, VectorNd& x)
{
}

#else // #ifndef HAVE_IPOPT
/// Solves the nonlinearly constrained quadratic program (potentially solves two nQPs, actually)
void ImpactConstraintHandler::solve_nqp(VectorNd& z, UnilateralConstraintProblemData& q, double max_time)
{
  // get the number of different types of each constraint
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned CL_IDX = N_CONTACTS*3;

  // TODO: set z to frictionless solution to start

  // mark starting time
  tms cstart;
  clock_t start = times(&cstart);

  // keep solving until we run out of time or all contact points are active
  while (true)
  {
    FILE_LOG(LOG_CONSTRAINT) << "Running NQP solve iteration with " << (q.N_ACT_CONTACTS) << " active contacts" << std::endl;

    // solve the nonlinearly constrained QP
    solve_nqp_work(q, z);

    // get the elapsed time
    const long TPS = sysconf(_SC_CLK_TCK);
    tms cstop;
    clock_t stop = times(&cstop);
    double elapsed = (double) (stop-start)/TPS;
    FILE_LOG(LOG_CONSTRAINT) << "Elapsed time: " << elapsed << std::endl; 

    // check whether we can mark any more contacts as active
    if (elapsed > max_time || q.N_ACT_CONTACTS == q.N_CONTACTS)
      break;

    // we can; mark next contact for solving
    q.N_ACT_K += q.contact_constraints[q.N_ACT_CONTACTS]->contact_NK/2;
    q.N_ACT_CONTACTS++;
  }
}

/// Solves the nonlinearly constrained quadratic program (does all of the work)
/**
 * \param x the solution is returned here; zeros will be returned at appropriate indices for inactive contacts
 */
void ImpactConstraintHandler::solve_nqp_work(UnilateralConstraintProblemData& q, VectorNd& x)
{
  const double INF = std::numeric_limits<double>::max();

  // setup constants
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_ACT_CONTACTS = q.N_ACT_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_IMP = q.N_CONSTRAINT_EQNS_IMP; 
  const unsigned CN_IDX = 0;
  const unsigned CS_IDX = N_ACT_CONTACTS;
  const unsigned CT_IDX = CS_IDX + N_ACT_CONTACTS;
  const unsigned CL_IDX = CT_IDX + N_ACT_CONTACTS;
  const unsigned NVARS = N_LIMITS + CL_IDX; 

  // setup the optimization data
  _ipsolver->epd = &q;
  _ipsolver->mu_c.resize(N_ACT_CONTACTS);
  _ipsolver->mu_visc.resize(N_ACT_CONTACTS);

  // setup true friction cone for every contact
  for (unsigned i=0; i< N_ACT_CONTACTS; i++)
  {
    _ipsolver->mu_c[i] = sqr(q.contact_constraints[i]->contact_mu_coulomb);
    _ipsolver->mu_visc[i] = (sqr(q.Cs_v[i]) + sqr(q.Ct_v[i])) *
                       sqr(q.contact_constraints[i]->contact_mu_viscous);
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
    SharedMatrixNd b1 = _A.block(0, N_CONSTRAINT_EQNS_IMP, 0, N_ACT_CONTACTS);
    SharedMatrixNd b2 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_ACT_CONTACTS, N_ACT_CONTACTS*2);
    SharedMatrixNd b3 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_ACT_CONTACTS*2, N_ACT_CONTACTS*3);
    SharedMatrixNd b4 = _A.block(0, N_CONSTRAINT_EQNS_IMP, N_ACT_CONTACTS*3, N_ACT_CONTACTS*3+N_LIMITS);

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
  const unsigned NONLIN_INEQUAL = N_ACT_CONTACTS;

  // init the QP matrix and vector
  H.resize(N_PRIMAL, N_PRIMAL);
  c.resize(H.rows());

  // setup row (block) 1 -- Cn * iM * [Cn' Cs Ct' L']
  unsigned col_start = 0, col_end = N_ACT_CONTACTS;
  unsigned row_start = 0, row_end = N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_ACT_CONTACTS;
  col_start = 0; col_end = N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_ACT_CONTACTS;
  col_start = 0; col_end = N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block 4) -- L * iM * [Cn' Cs' Ct' L']
  row_start = row_end; row_end += N_LIMITS;
  col_start = 0; col_end = N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_ACT_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // copy to row block 1 (contact normals)
  q.Cn_iM_CnT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CnT);
  q.Cn_iM_CsT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CsT);
  q.Cn_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cn_iM_CtT);
  q.Cn_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_LIMITS, Cn_iM_LT);
  
  // copy to row block 2 (first contact tangents)
  MatrixNd::transpose(Cn_iM_CsT, Cs_iM_CnT);
  q.Cs_iM_CsT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cs_iM_CsT);
  q.Cs_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Cs_iM_CtT);
  q.Cs_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_LIMITS, Cs_iM_LT);

  // copy to row block 3 (second contact tangents)
  MatrixNd::transpose(Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(Cs_iM_CtT, Ct_iM_CsT);
  q.Ct_iM_CtT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_ACT_CONTACTS, Ct_iM_CtT);
  q.Ct_iM_LT.get_sub_mat(0, N_ACT_CONTACTS, 0, N_LIMITS, Ct_iM_LT);

  // copy to row block 6 (limits)
  MatrixNd::transpose(Cn_iM_LT, L_iM_CnT);
  MatrixNd::transpose(Cs_iM_LT, L_iM_CsT);
  MatrixNd::transpose(Ct_iM_LT, L_iM_CtT);
  q.L_iM_LT.get_sub_mat(0, N_LIMITS, 0, N_LIMITS, L_iM_LT);

  // get components of c
  SharedVectorNd Cn_v = c.segment(0, N_ACT_CONTACTS); 
  SharedVectorNd Cs_v = c.segment(N_ACT_CONTACTS, N_ACT_CONTACTS*2); 
  SharedVectorNd Ct_v = c.segment(N_ACT_CONTACTS*2, N_ACT_CONTACTS*3); 
  SharedVectorNd L_v = c.segment(N_ACT_CONTACTS*3, N_ACT_CONTACTS*3+N_LIMITS); 

  // setup c 
  q.Cn_v.get_sub_vec(0, N_ACT_CONTACTS, Cn_v);
  q.Cs_v.get_sub_vec(0, N_ACT_CONTACTS, Cs_v);
  q.Ct_v.get_sub_vec(0, N_ACT_CONTACTS, Ct_v);
  L_v = q.L_v;

  // ****** now setup linear inequality constraints ******

  // determine whether to use kappa constraint
  const unsigned KAPPA = 1;

  // determine number of linear inequality constraints
  const unsigned N_INEQUAL = q.N_CONTACTS + N_LIMITS + KAPPA;

  // get Cn sub blocks
  SharedConstMatrixNd sub_Cn_Cn = q.Cn_iM_CnT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd sub_Cn_Cs = q.Cn_iM_CsT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd sub_Cn_Ct = q.Cn_iM_CtT.block(0, q.N_CONTACTS, 0, q.N_ACT_CONTACTS);
  SharedConstMatrixNd sub_Cn_L = q.Cn_iM_LT.block(0, q.N_CONTACTS, 0, q.N_LIMITS);

  // setup Cn block
  MatrixNd& Cn_block = _ipsolver->Cn_block;
  Cn_block.resize(q.N_CONTACTS, NVARS);
  Cn_block.set_sub_mat(0, 0, sub_Cn_Cn); 
  Cn_block.set_sub_mat(0, N_ACT_CONTACTS, sub_Cn_Cs); 
  Cn_block.set_sub_mat(0, N_ACT_CONTACTS*2, sub_Cn_Ct); 
  Cn_block.set_sub_mat(0, N_ACT_CONTACTS*3, sub_Cn_L); 

  // verify that L_block can be reset
  _ipsolver->L_block.reset();
  _ipsolver->Cn_v.reset();
  _ipsolver->L_block = L_block;
  _ipsolver->Cn_v = q.Cn_v.segment(0, N_CONTACTS);

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

   FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_nqp_work() entered" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cs': " << std::endl << q.Cn_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Ct': " << std::endl << q.Cn_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Jx': " << std::endl << q.Cn_iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Cs': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Ct': " << std::endl << q.Cs_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * L': " << std::endl << q.Cs_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Jx': " << std::endl << q.Cs_iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Ct': " << std::endl << q.Ct_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * L': " << std::endl << q.Ct_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Jx': " << std::endl << q.Ct_iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * L': " << std::endl << q.L_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * Jx': " << std::endl << q.L_iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_CONSTRAINT) << "c vector: " << c << std::endl;

  // setup ipopt options
  _app.Options()->SetIntegerValue("print_level", 0);
  _app.Options()->SetNumericValue("constr_viol_tol", 0.0005);
//  _app.Options()->SetIntegerValue("max_iter", 10000);
//  _app.Options()->SetStringValue("derivative_test", "second-order");

  // set the ipsolver tolerance on the Coulomb friction and kappa constraints
  _ipsolver->_tol = 0.0;

  // solve the nonlinear QP using the interior-point algorithm 
  _app.Initialize();
  Ipopt::ApplicationReturnStatus status = _app.OptimizeTNLP(_ipsolver);

  // look for acceptable solve conditions
  if (!(status == Ipopt::Solve_Succeeded || 
        status == Ipopt::Solved_To_Acceptable_Level)) 
    throw std::runtime_error("Could not solve nonlinearly constrained QP");

  // get the final solution out
  SharedVectorNd cn = _ipsolver->z.segment(0, N_ACT_CONTACTS);
  SharedVectorNd cs = _ipsolver->z.segment(N_ACT_CONTACTS, N_ACT_CONTACTS*2);
  SharedVectorNd ct = _ipsolver->z.segment(N_ACT_CONTACTS*2, N_ACT_CONTACTS*3);
  SharedVectorNd l =  _ipsolver->z.segment(N_ACT_CONTACTS*3, N_ACT_CONTACTS*3+q.N_LIMITS);

  // put x in the expected format
  x.resize(q.N_VARS);
  x.set_sub_vec(q.CN_IDX, cn);
  x.set_sub_vec(q.CS_IDX, cs);
  x.set_sub_vec(q.CT_IDX, ct);
  x.set_sub_vec(q.L_IDX, l);

  FILE_LOG(LOG_CONSTRAINT) << "nonlinear QP solution: " << x << std::endl; 
  if (LOGGING(LOG_CONSTRAINT))
  {
    VectorNd workv;
    SharedVectorNd xsub = _ipsolver->z.segment(0, c.rows());
    H.mult(xsub, workv) *= 0.5;
    workv += c;
    FILE_LOG(LOG_CONSTRAINT) << "(signed) computed energy dissipation: " << xsub.dot(workv) << std::endl;
  }
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_nqp() exited" << std::endl;
}
#endif // #ifndef HAVE_IPOPT

