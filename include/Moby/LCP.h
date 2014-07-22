/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_LCP_H
#define _MOBY_LCP_H

#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>

namespace Moby {

class LCP
{
  public:
    LCP();
    bool lcp_lemke_regularized(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, int min_exp = -20, unsigned step_exp = 1, int max_exp = 1, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_lemke_regularized(const Ravelin::SparseMatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, int min_exp = -20, unsigned step_exp = 4, int max_exp = 20, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_lemke(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_lemke(const Ravelin::SparseMatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_fast(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double zero_tol = -1.0);
    bool fast_pivoting(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double eps = std::sqrt(std::numeric_limits<double>::epsilon()));

  private:
    static void log_failure(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q);
    static void set_basis(unsigned n, unsigned count, std::vector<unsigned>& bas, std::vector<unsigned>& nbas);
    static unsigned rand_min(const Ravelin::VectorNd& v, double zero_tol);

    // temporaries for regularized solver
    Ravelin::MatrixNd _MM;
    Ravelin::VectorNd _wx;

    // temporaries for fast pivoting solver
    Ravelin::VectorNd _z, _w, _qbas, _qprime;
    Ravelin::MatrixNd _Msub, _Mmix, _M;

    // temporaries for Lemke solver
    Ravelin::VectorNd _d, _Be, _u, _z0, _x, _dl, _xj, _dj, _wl, _result;
    Ravelin::VectorNd _restart_z0;
    Ravelin::MatrixNd _Bl, _Al, _t1, _t2;
    std::vector<unsigned> _all, _tlist, _bas, _nonbas, _j, _max_elm;

    // temporary for sparse Lemke solver
    Ravelin::SparseMatrixNd _sBl;
    Ravelin::SparseMatrixNd _MMs, _MMx, _eye, _zero, _diag_lambda;

    // linear algebra
    Ravelin::LinAlgd _LA;
}; // end class

} // end namespace 

#endif

