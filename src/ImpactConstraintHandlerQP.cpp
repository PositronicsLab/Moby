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
#include <Moby/LCPSolverException.h>
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

/// Solves the quadratic program (potentially solves two QPs, actually)
void ImpactConstraintHandler::solve_qp(VectorNd& z, UnilateralConstraintProblemData& q)
{
  const unsigned N_TOTAL = q.N_VARS + q.N_CONTACTS + q.N_LIMITS + q.N_K_TOTAL + 1;

  // mark starting time
  tms cstart;
  clock_t start = times(&cstart);

  // solve the QP
  solve_qp_work(q, z);

  // get the elapsed time
  const long TPS = sysconf(_SC_CLK_TCK);
  tms cstop;
  clock_t stop = times(&cstop);
  double elapsed = (double) (stop - start)/TPS;
  FILE_LOG(LOG_CONSTRAINT) << "Elapsed time: " << elapsed << std::endl;
}

/// Computes the kinetic energy of the system using the current impulse set
double ImpactConstraintHandler::calc_ke(UnilateralConstraintProblemData& q, const VectorNd& z)
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

/// Solves the quadratic program (does all of the work) as an LCP
/**
 * \note this is the version without joint friction forces
 * \param z the solution is returned here; zeros are returned at appropriate
 *        places for inactive contacts
 */
void ImpactConstraintHandler::solve_qp_work(UnilateralConstraintProblemData& epd, VectorNd& z)
{
  // implicit constraints not handled at the moment
  assert(epd.N_CONSTRAINT_EQNS_IMP == 0);

  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << epd.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cs': " << std::endl << epd.Cn_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Ct': " << std::endl << epd.Cn_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * L': " << std::endl << epd.Cn_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Cs': " << std::endl << epd.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Ct': " << std::endl << epd.Cs_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * L': " << std::endl << epd.Cs_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Ct': " << std::endl << epd.Ct_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * L': " << std::endl << epd.Ct_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * L': " << std::endl << epd.L_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << epd.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * v: " << epd.Cs_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * v: " << epd.Ct_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << epd.L_v << std::endl;

  // setup new indices
  const unsigned N_CONTACTS = epd.N_CONTACTS;
  const unsigned xCN_IDX = 0;
  const unsigned xCS_IDX = xCN_IDX + N_CONTACTS;
  const unsigned xCT_IDX = xCS_IDX + N_CONTACTS;
  const unsigned xNCS_IDX = xCT_IDX + N_CONTACTS;
  const unsigned xNCT_IDX = xNCS_IDX + N_CONTACTS;
  const unsigned xL_IDX = xNCT_IDX + N_CONTACTS;
  const unsigned N_VARS = xL_IDX + epd.N_LIMITS;

  // init the QP matrix and vector
  const unsigned N_INEQUAL = epd.N_CONTACTS + epd.N_K_TOTAL + epd.N_LIMITS + 1;
  _MM.set_zero(N_VARS + N_INEQUAL, N_VARS + N_INEQUAL);
  _qq.resize(_MM.rows());

  // get the necessary matrices for setting up the LCP
  SharedMatrixNd M = _MM.block(N_VARS, _MM.rows(), 0, N_VARS);
  SharedMatrixNd H = _MM.block(0, N_VARS, 0, N_VARS);
  SharedVectorNd c = _qq.segment(0, N_VARS);
  SharedVectorNd q = _qq.segment(N_VARS, _qq.size()).set_zero();

  // setup the QP
  MatrixNd J(0, N_VARS);
  VectorNd Jv(0);
  SharedMatrixNd A = J.block(0, J.rows(), 0, J.columns());
  SharedVectorNd b = Jv.segment(0, Jv.rows());
  setup_QP(epd, H, c, M, q, A, b);

  // set M = -M'
  SharedMatrixNd MT = _MM.block(0, N_VARS, N_VARS, _MM.rows());
  MatrixNd::transpose(M, MT);
  MT.negate();

  FILE_LOG(LOG_CONSTRAINT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_CONSTRAINT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "M matrix: " << std::endl << M;
  FILE_LOG(LOG_CONSTRAINT) << "q vector: " << q << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << "LCP vector: " << _qq << std::endl;

  // init z 
  z.resize(_qq.rows());

  // try warmstarting if possible
  if (z.size() == _zlast.size())
    z = _zlast;

  // solve the LCP using Lemke's algorithm
  #ifdef USE_QLCPD
  VectorNd lb(c.size()), ub(c.size());
  lb.set_zero();
  ub.set_one() *= 1e+29;
  if (!_qp.qp_activeset(H, c, lb, ub, M, q, A, b, z))
  {
    FILE_LOG(LOG_CONSTRAINT) << "QLCPD failed to solve; finding closest feasible point" << std::endl;

    // QP solver not successful by default; attempt to find the closest
    // feasible point
    if (!_qp.find_closest_feasible(lb, ub, M, q, A, b, z))
    {
      // QP solver failed completely; use Lemke's Algorithm as backup
      q.negate();
      if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
        throw LCPSolverException();
    }
    else
    {
      FILE_LOG(LOG_CONSTRAINT) << "updating q; q=" << q << std::endl;

      // found closest feasible point; compute M*z - q
      M.mult(z, _workv) -= q;
      for (unsigned i=0; i< _workv.size(); i++)
        if (_workv[i] < 0.0)
          q[i] += _workv[i] - NEAR_ZERO;
      FILE_LOG(LOG_CONSTRAINT) << "            q'=" << q << std::endl;

      // now attempt to solve the QP again
      if (!_qp.qp_activeset(H, c, lb, ub, M, q, A, b, z))
      {
        // QP solver failed on second attempt; use Lemke's Algorithm as backup
        q.negate();
        if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
          throw LCPSolverException();
      }
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "QLCPD solution: " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "M: " << std::endl << M;
  FILE_LOG(LOG_CONSTRAINT) << "q: " << q << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "M*z - q: " << (M.mult(z.get_sub_vec(0,M.columns(),_workv2), _workv) -= q) << std::endl;

  #else
  // negate q (it was in form Mx >= q, needs to be in Mx + q >= 0)
  q.negate();

  // attempt fast lcp solve
  if (!_lcp.lcp_fast_regularized(_MM, _qq, z, -20, 4, -8))
  {
    // Lemke does not like warm starting
    z.set_zero();

    if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
      throw LCPSolverException();
  }

  // output reported LCP solution
  FILE_LOG(LOG_CONSTRAINT) << "LCP solution: " << z << std::endl;
  #endif

  // store zlast
  _zlast = z;

  // get relevant forces
  SharedVectorNd cn = _zlast.segment(0, N_CONTACTS);
  SharedVectorNd cs = _zlast.segment(N_CONTACTS, N_CONTACTS*2);
  SharedVectorNd ct = _zlast.segment(N_CONTACTS*2, N_CONTACTS*3);
  SharedVectorNd ncs = _zlast.segment(N_CONTACTS*3, N_CONTACTS*4);
  SharedVectorNd nct = _zlast.segment(N_CONTACTS*4, N_CONTACTS*5);
  SharedVectorNd l = _zlast.segment(N_CONTACTS*5, N_CONTACTS*5+epd.N_LIMITS);

  // put z in the expected format (full contact forces)
  z.set_zero(epd.N_VARS);
  z.set_sub_vec(epd.CN_IDX, cn);
  z.set_sub_vec(epd.CS_IDX, cs);
  z.set_sub_vec(epd.CT_IDX, ct);
  z.set_sub_vec(epd.NCS_IDX, ncs);
  z.set_sub_vec(epd.NCT_IDX, nct);
  z.set_sub_vec(epd.L_IDX, l);

  FILE_LOG(LOG_CONSTRAINT) << "QP solution: " << z << std::endl;

  if (LOGGING(LOG_CONSTRAINT))
  {
    VectorNd workv;
    SharedVectorNd zsub = _zlast.segment(0, c.rows());
    H.mult(zsub, workv) *= 0.5;
    workv += c;
    FILE_LOG(LOG_CONSTRAINT) << "(signed) computed energy dissipation: " << zsub.dot(workv) << std::endl;
  }
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_qp_work() exited" << std::endl;
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 * \param z the solution is returned here; zeros are returned at appropriate
 *        places for inactive contacts
 */
void ImpactConstraintHandler::setup_QP(UnilateralConstraintProblemData& epd, SharedMatrixNd& H, SharedVectorNd& c, SharedMatrixNd& M, SharedVectorNd& q, SharedMatrixNd& A, SharedVectorNd& b)
{
  // implicit constraints not handled at the moment
  assert(epd.N_CONSTRAINT_EQNS_IMP == 0);

  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << epd.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cs': " << std::endl << epd.Cn_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Ct': " << std::endl << epd.Cn_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * L': " << std::endl << epd.Cn_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Cs': " << std::endl << epd.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Ct': " << std::endl << epd.Cs_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * L': " << std::endl << epd.Cs_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Ct': " << std::endl << epd.Ct_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * L': " << std::endl << epd.Ct_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * L': " << std::endl << epd.L_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << epd.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * v: " << epd.Cs_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * v: " << epd.Ct_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << epd.L_v << std::endl;

  // get useful constants
  const unsigned N_CONTACTS = epd.N_CONTACTS;

  // setup new indices
  const unsigned xCN_IDX = 0;
  const unsigned xCS_IDX = xCN_IDX + N_CONTACTS;
  const unsigned xCT_IDX = xCS_IDX + N_CONTACTS;
  const unsigned xNCS_IDX = xCT_IDX + N_CONTACTS;
  const unsigned xNCT_IDX = xNCS_IDX + N_CONTACTS;
  const unsigned xL_IDX = xNCT_IDX + N_CONTACTS;
  const unsigned N_VARS = xL_IDX + epd.N_LIMITS;

  // init the QP matrix and vector
  const unsigned N_INEQUAL = epd.N_CONTACTS + epd.N_K_TOTAL + epd.N_LIMITS;

  // setup row (block) 1 -- Cn * iM * [Cn' Cs Ct' -Cs' -Ct' L' ]
  unsigned col_start = 0, col_end = epd.N_CONTACTS;
  unsigned row_start = 0, row_end = epd.N_CONTACTS;
  SharedMatrixNd Cn_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cn_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cn_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cn_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_LIMITS;
  SharedMatrixNd Cn_iM_LT = H.block(row_start, row_end, col_start, col_end);

  // setup row (block) 2 -- Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' Jx']
  row_start = row_end; row_end += epd.N_CONTACTS;
  col_start = 0; col_end = epd.N_CONTACTS;
  SharedMatrixNd Cs_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cs_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cs_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Cs_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_LIMITS;
  SharedMatrixNd Cs_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Cs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 3 -- Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += epd.N_CONTACTS;
  col_start = 0; col_end = epd.N_CONTACTS;
  SharedMatrixNd Ct_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Ct_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Ct_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd Ct_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_LIMITS;
  SharedMatrixNd Ct_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd Ct_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 4 -- -Cs * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += epd.N_CONTACTS;
  SharedMatrixNd NCs_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block) 5 -- -Ct * iM * [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += epd.N_CONTACTS;
  SharedMatrixNd NCt_block = H.block(row_start, row_end, 0, col_end);

  // setup row (block 6) -- L * iM *  [Cn' Cs' Ct' -Cs' -Ct' L' ]
  row_start = row_end; row_end += epd.N_LIMITS;
  col_start = 0; col_end = epd.N_CONTACTS;
  SharedMatrixNd L_iM_CnT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd L_iM_CsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd L_iM_CtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd L_iM_NCsT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_CONTACTS;
  SharedMatrixNd L_iM_NCtT = H.block(row_start, row_end, col_start, col_end);
  col_start = col_end; col_end += epd.N_LIMITS;
  SharedMatrixNd L_iM_LT = H.block(row_start, row_end, col_start, col_end);
  SharedMatrixNd L_block = H.block(row_start, row_end, 0, col_end);

  // copy to row block 1 (contact normals)
  epd.Cn_iM_CnT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CnT);
  epd.Cn_iM_CsT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CsT);
  epd.Cn_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cn_iM_CtT);
  (Cn_iM_NCsT = Cn_iM_CsT).negate();
  (Cn_iM_NCtT = Cn_iM_CtT).negate();
  epd.Cn_iM_LT.get_sub_mat(0, N_CONTACTS, 0, epd.N_LIMITS, Cn_iM_LT);

  // copy to row block 2 (first tangents)
  MatrixNd::transpose(Cn_iM_CsT, Cs_iM_CnT);
  epd.Cs_iM_CsT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cs_iM_CsT);
  epd.Cs_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Cs_iM_CtT);
  (Cs_iM_NCsT = Cs_iM_CsT).negate();
  (Cs_iM_NCtT = Cs_iM_CtT).negate();
  epd.Cs_iM_LT.get_sub_mat(0, N_CONTACTS, 0, epd.N_LIMITS, Cs_iM_LT);

  // copy to row block 3 (second tangents)
  MatrixNd::transpose(Cn_iM_CtT, Ct_iM_CnT);
  MatrixNd::transpose(Cs_iM_CtT, Ct_iM_CsT);
  epd.Ct_iM_CtT.get_sub_mat(0, N_CONTACTS, 0, N_CONTACTS, Ct_iM_CtT);
  (Ct_iM_NCsT = Ct_iM_CsT).negate();
  (Ct_iM_NCtT = Ct_iM_CtT).negate();
  epd.Ct_iM_LT.get_sub_mat(0, N_CONTACTS, 0, epd.N_LIMITS, Ct_iM_LT);

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
  epd.L_iM_LT.get_sub_mat(0, epd.N_LIMITS, 0, epd.N_LIMITS, L_iM_LT);

  // get shared vectors to components of c
  SharedVectorNd Cn_v = c.segment(xCN_IDX, xCS_IDX);
  SharedVectorNd Cs_v = c.segment(xCS_IDX, xCT_IDX);
  SharedVectorNd Ct_v = c.segment(xCT_IDX, xNCS_IDX);
  SharedVectorNd nCs_v = c.segment(xNCS_IDX, xNCT_IDX);
  SharedVectorNd nCt_v = c.segment(xNCT_IDX, xL_IDX);
  SharedVectorNd L_v = c.segment(xL_IDX, N_VARS);

  // setup c
  epd.Cn_v.get_sub_vec(0, N_CONTACTS, Cn_v);
  epd.Cs_v.get_sub_vec(0, N_CONTACTS, Cs_v);
  epd.Ct_v.get_sub_vec(0, N_CONTACTS, Ct_v);
  (nCs_v = Cs_v).negate();
  (nCt_v = Ct_v).negate();
  L_v = epd.L_v;

  // ------- setup M/q -------
  // setup the Cn*v+ >= 0 constraint
  // Cn*(inv(M)*impulses + v) >= 0, Cn*inv(M)*impulses >= -Cn*v
  row_start = 0; row_end = N_CONTACTS;
  M.block(row_start, row_end, 0, M.columns()) = H.block(0, N_CONTACTS, 0, H.columns());
/*
  M.set_sub_mat(row_start, xCN_IDX, epd.Cn_iM_CnT, eTranspose);
  M.set_sub_mat(row_start, xCS_IDX, epd.Cn_iM_CsT, eTranspose);
  M.set_sub_mat(row_start, xCT_IDX, epd.Cn_iM_CtT, eTranspose);
  M.set_sub_mat(row_start, xNCS_IDX, epd.Cn_iM_CsT, eTranspose);
  M.set_sub_mat(row_start, xNCT_IDX, epd.Cn_iM_CtT, eTranspose);
  M.set_sub_mat(row_start, xL_IDX, epd.Cn_iM_LT, eTranspose);
  M.block(row_start, row_end, xNCS_IDX, xL_IDX).negate();
*/
  FILE_LOG(LOG_CONSTRAINT) << "M: " << std::endl << M;
  q.set_sub_vec(row_start, epd.Cn_v);
  row_start = row_end; row_end += epd.N_LIMITS;

  // setup the L*v+ >= 0 constraint
  M.block(row_start, row_end, xCN_IDX, N_VARS) = L_block;
  q.set_sub_vec(row_start, epd.L_v);
  row_start = row_end; row_end += epd.N_CONTACTS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
  {
    // initialize the contact velocity
    double vel = std::sqrt(sqr(epd.Cs_v[i]) + sqr(epd.Ct_v[i]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< epd.contact_constraints[i]->contact_NK/2; j++)
    {
      double theta = (double) j/(epd.contact_constraints[i]->contact_NK/2-1) * M_PI_2;
      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      M(row_start, xCN_IDX+i) = epd.contact_constraints[i]->contact_mu_coulomb;
      M(row_start, xCS_IDX+i) = -ct;
      M(row_start, xNCS_IDX+i) = -ct;
      M(row_start, xCT_IDX+i) = -st;
      M(row_start, xNCT_IDX+i) = -st;

      // setup the viscous friction component
      q[row_start] = epd.contact_constraints[i]->contact_mu_viscous * vel;
      row_start++;
    }
  }

  // negate q
  q.negate();
}

