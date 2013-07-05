/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_SPATIAL_H
#define _MOBY_SPATIAL_H

#include <vector>
#include <Ravelin/MissizeException.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/SVelocityd.h>
#include <Ravelin/SForced.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/SpatialABInertiad.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>

namespace Moby {

/// Converts an STL vector of forces to a matrix (type X)
template <class X>
X& to_matrix(const std::vector<Ravelin::SForced>& w, X& m)
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

/// Converts an STL vector of momenta to a matrix (type X)
template <class X>
X& to_matrix(const std::vector<Ravelin::SMomentumd>& w, X& m)
{
  // TODO: check whether SMomentum should be defined with linear components on
  // top and angular on bottom
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, w.size());
  double* data = m.data();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    Ravelin::Vector3d f = w[i].get_linear();  
    Ravelin::Vector3d t = w[i].get_angular();
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
  }

  return m;
}

/// Converts an STL vector of spatial velocities to a force matrix (type X)
template <class X>
X& transpose_to_matrix(const std::vector<Ravelin::SVelocityd>& t, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, t.size());
  double* data = m.data();
  for (unsigned k=0, i=0; i< t.size(); i++)
  {
    Ravelin::Vector3d lin = t[i].get_linear();  
    Ravelin::Vector3d ang = t[i].get_angular();
    data[k++] = lin[0];  data[k++] = lin[1];  data[k++] = lin[2];
    data[k++] = ang[0];  data[k++] = ang[1];  data[k++] = ang[2];
  }

  return m;
}

/// Converts an STL vector of spatial accelerations to a force matrix (type X)
template <class X>
X& transpose_to_matrix(const std::vector<Ravelin::SAcceld>& t, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, t.size());
  double* data = m.data();
  for (unsigned k=0, i=0; i< t.size(); i++)
  {
    Ravelin::Vector3d lin = t[i].get_linear();  
    Ravelin::Vector3d ang = t[i].get_angular();
    data[k++] = lin[0];  data[k++] = lin[1];  data[k++] = lin[2];
    data[k++] = ang[0];  data[k++] = ang[1];  data[k++] = ang[2];
  }

  return m;
}

/// Computes the "spatial dot product" between a vector of velocities and a vector of forces and returns the result in the matrix container (X)
template <class X>
X& transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const std::vector<Ravelin::SVector6d>& v, X& result)
{
  result.resize(t.size(), v.size());
  double* data = result.data();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< v.size(); j++)
      data[k++] = t[i].dot(v[j]);

  return result;
}

/// Computes the "spatial dot product" between a vector of velocities and a vector of forces and returns the result in the matrix container (X)
template <class X>
X& transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const std::vector<Ravelin::SForced>& w, X& result)
{
  result.resize(t.size(), w.size());
  double* data = result.data();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< w.size(); j++)
      data[k++] = t[i].dot(w[j]);

  return result;
}

/// Computes the "spatial dot product" between a vector of velocities and a matrix or vector and returns the result in the matrix container (X)
template <class Y, class X>
X& transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const Y& y, X& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(t.size(), y.columns(), false);
  double* data = result.data();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< y.columns(); j++)
      data[k++] = t[i].dot(y.column(j));

  return result;
}

/// Computes the "spatial dot product" between a vector of velocities and a force and returns the result in the matrix container (X)
template <class X>
X& transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const Ravelin::SForced& w, X& result)
{
  result.resize(t.size(), 1, false);
  double* data = result.data();
  for (unsigned i=0; i< t.size(); i++)
    data[i] = t[i].dot(w);

  return result;
}

/// Computes the "spatial dot product" between a vector of forces and a twist and returns the result in the matrix container (X)
template <class X>
X& transpose_mult(const std::vector<Ravelin::SForced>& w, const Ravelin::SVelocityd& t, X& result)
{
  result.resize(w.size());
  double* data = result.data();
  for (unsigned i=0; i< w.size(); i++)
    data[i] = w[i].dot(t);

  return result;
}

Ravelin::MatrixNd& mult(const Ravelin::SpatialABInertiad& I, const std::vector<Ravelin::SVector6d>& s, Ravelin::MatrixNd& result);
Ravelin::MatrixNd& mult(const Ravelin::SpatialABInertiad& I, const std::vector<Ravelin::SVelocityd>& s, Ravelin::MatrixNd& result);
Ravelin::VectorNd& concat(const Ravelin::VectorNd& v, const Ravelin::SForced& w, Ravelin::VectorNd& result);
Ravelin::SVector6d mult(const std::vector<Ravelin::SVector6d>& t, const Ravelin::VectorNd& v);
Ravelin::MatrixNd& mult(const std::vector<Ravelin::SVector6d>& v, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result);
Ravelin::MatrixNd& mult(const std::vector<Ravelin::SVelocityd>& t, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result);
Ravelin::MatrixNd& mult(const std::vector<Ravelin::SAcceld>& t, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result);
Ravelin::MatrixNd& mult(const std::vector<Ravelin::SForced>& w, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result);
Ravelin::SVector6d mult(const std::vector<Ravelin::SVelocityd>& t, const Ravelin::VectorNd& v);
Ravelin::SVector6d mult(const std::vector<Ravelin::SAcceld>& t, const Ravelin::VectorNd& v);
Ravelin::SVector6d mult(const std::vector<Ravelin::SForced>& w, const Ravelin::VectorNd& v);

} // end namespace Moby

#endif

