/****************************************************************************
 * Copyright 2016 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/QP.h>

using namespace Ravelin;
using namespace Moby;

/// Converts the QP 0.5*x'*H*x + c'*x, ub >= x >= lb, Mx >= q, Ax = b where H is symmetric, PSD into a generalized linear complementarity problem
void QP::convex_qp_to_glcp(const MatrixNd& H, const VectorNd& c, const MatrixNd& M, const VectorNd& q, const MatrixNd& A, const VectorNd& b, const VectorNd& lb, const VectorNd& ub, MatrixNd& MM, VectorNd& qq, VectorNd& zlo, VectorNd& zhi)
{
  const double INF = (double) std::numeric_limits<double>::max();

  // get number of variables and number of constraints
  const unsigned NVARS = c.size();
  const unsigned NINEQ = q.size();
  const unsigned NEQ = b.size();

  // look for quick exit
  if (NVARS == 0)
  {
    MM.resize(0,0);
    qq.resize(0);
    zlo.resize(0);
    zhi.resize(0);
    return;
  } 

  // initialize MMM and qqq
  MM.resize(NVARS + NINEQ + NEQ, NVARS + NINEQ + NEQ);
  qq.resize(NVARS + NINEQ + NEQ);

  // form MM
  MM.block(0, NVARS, 0, NVARS) = H;
  SharedMatrixNd MT_block = MM.block(0, NVARS, NVARS, NVARS+NINEQ);
  SharedMatrixNd AT_block = MM.block(0, NVARS, NVARS+NINEQ, NVARS+NINEQ+NEQ);
  MatrixNd::transpose(M, MT_block);
  MatrixNd::transpose(A, AT_block);
  MT_block.negate();
  AT_block.negate();
  MM.block(NVARS, NVARS+NINEQ, 0, NVARS) = M;
  MM.block(NVARS+NINEQ, NVARS+NINEQ+NEQ, 0, NVARS) = A;
  MM.block(NVARS, MM.rows(), NVARS, MM.rows()).set_zero();

  // form qq
  qq.segment(0, NVARS) = c;
  (qq.segment(NVARS, NVARS+NINEQ) = q).negate();
  (qq.segment(NVARS+NINEQ, NVARS+NINEQ+NEQ) = b).negate(); 

  // setup the lower and upper bounds for primal variables
  zlo.resize(NVARS + NINEQ + NEQ); 
  zhi.resize(NVARS + NINEQ + NEQ); 
  zlo.set_sub_vec(0, lb);
  zhi.set_sub_vec(0, ub);

  // setup the lower and upper bounds for dual variables
  zlo.segment(NVARS, NVARS+NINEQ).set_zero();
  zlo.segment(NVARS+NINEQ, NVARS+NINEQ+NEQ).set_one() *= -INF;
  zhi.segment(NVARS, zhi.size()).set_one() *= INF;
}

