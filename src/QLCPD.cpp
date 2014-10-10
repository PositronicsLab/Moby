/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <numeric>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <iostream>
#include <Ravelin/select>
#include <Ravelin/cblas.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/SingularException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/NumericalException.h>
#include <Moby/QLCPD.h>

using namespace Moby;
using namespace Ravelin;
using std::endl;
using std::vector;
using std::pair;
using std::map;
using std::make_pair;
using boost::shared_ptr;

// setup structures for Fortran
extern "C"
{
  extern struct
  {
    double eps, t0l, emin;
  } epsc_;

  extern struct
  {
    double sgnf;
    int nrep, npiv, nres;
  } repc_;

  extern struct
  {
    int mx, mxmc;
  } refactorc_;

  extern struct
  {
    int mxm1;
  } mxm1c_;

  extern struct
  {
    double rgnorm, vstep;
    int iter, npv, nfn, ngr;
  } infoc_;

  extern struct
  {
    int kk, ll, kkk, lll, mxws, mxlws;
  } wsc_;
} // end extern "C"

// external Fortran call
extern "C"
{
  void qlcpd_(int* n, int* m, int* k, int* kmax, int* maxg,
             double* a, int* la, double* x, double* bl, double* bu,
             double *f, double* fmin, double *g,
             double* r, double* w, double* e, int* ls, double* alp, int* lp,
             int* mlp, int* peq, double* ws, int* lws, double* v,
             int* nv, int* lin, double* rgtol, int* mode, int* ifail, int* mxgr,
             int* iprint, int* nout);

  // catches fortran exit
  void _gfortran_stop_string(const char*, int len)
  {
    throw std::runtime_error("Fortran exit");
  }
}
  
// Sole constructor
QLCPD::QLCPD()
{
}

/// Calls qlcpd
void QLCPD::qlcpd(int* n, int* m, int* k, int* kmax, int* maxg,
             double* a, int* la, double* x, double* bl, double* bu,
             double *f, double* fmin, double *g,
             double* r, double* w, double* e, int* ls, double* alp, int* lp,
             int* mlp, int* peq, double* ws, int* lws, double* v,
             int* nv, int* lin, double* rgtol, int* mode, int* ifail, int* mxgr,
             int* iprint, int* nout)
{
  qlcpd_(n, m, k, kmax, maxg,
             a, la,  x,  bl,  bu,
             f,  fmin, g,
              r,  w,  e, ls,  alp, lp,
             mlp, peq,  ws, lws,  v,
             nv, lin,  rgtol, mode, ifail, mxgr,
             iprint, nout);
}

/// Active-set QP algorithm for solving QPs 
/**
 * \param H the quadratic term
 * \param c the linear term
 * \param lb the lower bound
 * \param ub the lower bound
 * \param M the linear inequality constraints (M*z >= q)
 * \param q the linear inequality bound (M*z >= q)
 * \param A the linear equality constraint (A*z = b)
 * \param b the linear equality constraint (A*z = b)
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
void QLCPD::set_values(int n, int m)
{
  // set k to dummy value (should not be used since mode < 2)
  int k = 0;
  int kmax = n;
  int maxg = std::min(6,kmax+1);

  // setup 'common' variables (generally tolerance parameters)
  epsc_.emin = 0.0;
  epsc_.eps = std::numeric_limits<double>::epsilon();
  epsc_.t0l = epsc_.eps;
  repc_.sgnf = 1e-8;
  repc_.nrep = 100;  // number of repeats
  repc_.npiv = 100;  // maximum number of pivoting operations 
  repc_.nres = 100;  // number of restarts allowed
  refactorc_.mxmc = 100;  // unused b/c dense?
  mxm1c_.mxm1 = std::min(m+1,n);
  wsc_.kk = (n+1)*n;
  wsc_.ll = 0;
  int nmi = n + m;
  int kkk = nmi + n + maxg*(maxg+1)/2+maxg*(kmax+5);
  _mxws = wsc_.kk + kkk + mxm1c_.mxm1*(mxm1c_.mxm1+1)/2+3*n+mxm1c_.mxm1;
  _mxlws = wsc_.ll + wsc_.lll + n + mxm1c_.mxm1 + nmi;
  wsc_.mxws = _mxws;
  wsc_.mxlws = _mxlws;
}

