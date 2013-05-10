/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_SPATIAL_H
#define _MOBY_SPATIAL_H

#include <vector>
#include <Ravelin/MissizeException.h>
#include <Ravelin/Twistd.h>
#include <Ravelin/Wrenchd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>

namespace Moby {

/// Converts an STL vector of wrenches to a matrix (type X)
template <class X>
X& to_matrix(const std::vector<Ravelin::Wrenchd>& w, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, w.size());
  double* data = m.data();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    Ravelin::Vector3d f = w[i].get_force();  
    Ravelin::Vector3d t = w[i].get_torque();
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
  }

  return m;
}

/// Computes the "spatial dot product" between a vector of twists and a wrench and returns the result in the matrix container (X)
template <class X>
X& transpose_mult(const std::vector<Ravelin::Twistd>& t, const Ravelin::Wrenchd& w, X& result)
{
  result.resize(t.size());
  double* data = result.data();
  for (unsigned i=0; i< t.size(); i++)
    data[i] = t[i].dot(w);

  return result;
}

Ravelin::Twistd spatial_cross(const Ravelin::Twistd& v1, const Ravelin::Twistd& v2);
Ravelin::VectorNd& concat(const Ravelin::VectorNd& v, const Ravelin::Wrenchd& w, Ravelin::VectorNd& result);
Ravelin::Twistd mult(const std::vector<Ravelin::Twistd>& t, const Ravelin::VectorNd& v);

} // end namespace Moby

#endif

