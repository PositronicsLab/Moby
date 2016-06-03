/****************************************************************************
 * Copyright 2016 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _OPT_QPOASES_H
#define _OPT_QPOASES_H

#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/Constants.h>
#include <qpOASES.hpp>

namespace Moby {

class qpOASES
{
  public:
   qpOASES() {
   }

  template <class Mat1, class Vec1, class Vec2, class Vec3, class Mat2, class Vec4, class Mat3, class Vec5, class Vec6>
  bool qp_activeset(const Mat1& H, const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat2& M, const Vec4& q, const Mat3& A, const Vec5& b, Vec6& z);

  template <class Mat1, class Vec1, class Vec2, class Vec3, class Mat2, class Vec4, class Mat3, class Vec5, class Vec6>
  bool qp_activeset_regularized(const Mat1& H, const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat2& M, const Vec4& q, const Mat3& A, const Vec5& b, Vec6& z, int min_exp=-20, int step_exp=4, int max_exp=10);

private:
   Ravelin::MatrixNd _HH, _X;
   Ravelin::VectorNd _lb, _ub, _lbA, _ubA;

   static bool rel_equal(double x, double y, double tol = NEAR_ZERO) { return std::fabs(x-y) <= tol*std::max(std::fabs(x), std::max(std::fabs(y), (double) 1.0)); }
}; // end class

#include "qpOASES.inl"

} // end namespace 

#endif

