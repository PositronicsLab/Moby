/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_QP_H
#define _MOBY_QP_H

#include <Ravelin/NonsquareMatrixException.h>

namespace Moby {

// class for solving positive definite QPs with box constraints
class QP
{
  private:
    Ravelin::LinAlgd _LA;

    // temporary variables for qp_gradproj and find_cauchy_point()
    Ravelin::VectorNd _p, _Gp, _Gx_c, _Gx, _H, _Wzz;
    Ravelin::MatrixNd _ZtGZ;
    Ravelin::VectorNd _xz, _cz, _lz, _uz, _dxz;
    Ravelin::VectorNd _r1, _r2, _w, _w1, _w2, _v, _xstar, _workv, _y;
    std::vector<unsigned> _inactive_set;
    std::vector<double> _alphas;

  public:
    #include "QP.inl"
}; // end class

} // end namespace

#endif

