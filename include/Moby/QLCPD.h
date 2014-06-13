/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _OPT_QPACTIVESET_H
#define _OPT_QPACTIVESET_H

#include <signal.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/Constants.h>

namespace Moby {

class QLCPD
{
  public:
    QLCPD();

  template <class Mat1, class Vec1, class Vec2, class Vec3, class Mat2, class Vec4, class Mat3, class Vec5, class Vec6>
  bool qp_activeset(const Mat1& H, const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat2& M, const Vec4& q, const Mat3& A, const Vec5& b, Vec6& z);

  template <class Vec1, class Vec2, class Mat1, class Vec3, class Mat2, class Vec4, class Vec5>
  bool find_closest_feasible(const Vec1& lb, const Vec2& ub, const Mat1& M, const Vec3& q, const Mat2& A, const Vec4& b, Vec5& z);

  template <class Vec1, class Vec2, class Vec3, class Mat1, class Vec4, class Mat2, class Vec5, class Vec6>
  bool lp_activeset(const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat1& M, const Vec4& q, const Mat2& A, const Vec5& b, Vec6& z);

  private:
    static bool rel_equal(double x, double y, double tol = NEAR_ZERO) { return std::fabs(x-y) <= tol*std::max(std::fabs(x), std::max(std::fabs(y), (double) 1.0)); }

    void qlcpd(int* n, int* m, int* k, int* kmax, int* maxg,
             double* a, int* la, double* x, double* bl, double* bu,
             double *f, double* fmin, double *g,
             double* r, double* w, double* e, int* ls, double* alp, int* lp,
             int* mlp, int* peq, double* ws, int* lws, double* v,
             int* nv, int* lin, double* rgtol, int* mode, int* ifail, int* mxgr,
             int* iprint, int* nout);

    // sets values for QLCPD
    void set_values(int n, int m);

    // workspace sizes for active set solver
    int _mxws, _mxlws;

    // temporaries for active set solver
    Ravelin::VectorNd _lb, _ub, _r, _e, _w, _v, _alp, _g;
    Ravelin::MatrixNd _X;
    std::vector<double> _ws;
    std::vector<int> _ls, _lp, _lws;
}; // end class

#include "QLCPD.inl"

} // end namespace 

#endif

