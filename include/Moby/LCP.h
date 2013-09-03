/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
    bool lcp_lemke_regularized(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, int min_exp = -20, unsigned step_exp = 4, int max_exp = 20, double piv_tol = -1.0, double zero_tol = -1.0);
bool lcp_lemke(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_lemke(const Ravelin::SparseMatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);

  private:
    // temporaries for regularized solver
    Ravelin::MatrixNd _MMx, _MM;
    Ravelin::VectorNd _qqx, _wx, _qq;

    // temporaries for Lemke solver
    Ravelin::VectorNd _d, _Be, _u, _z0, _x, _dl, _xj, _dj, _wl, _result;
    Ravelin::MatrixNd _Bl, _Al, _t1, _t2;
    std::vector<unsigned> _all, _tlist, _bas, _nonbas, _j, _max_elm;

    // temporary for sparse Lemke solver
    Ravelin::SparseMatrixNd _sBl;

    // linear algebra
    Ravelin::LinAlgd _LA;
}; // end class

} // end namespace 

#endif

