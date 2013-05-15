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
    bool lcp_lemke_regularized(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, int min_exp = -20, unsigned step_exp = 4, int max_exp = 20, double piv_tol = -1.0, double zero_tol = -1.0);
bool lcp_lemke(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);
    bool lcp_lemke(const Ravelin::SparseMatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z, double piv_tol = -1.0, double zero_tol = -1.0);
}; // end class

} // end namespace 

#endif

