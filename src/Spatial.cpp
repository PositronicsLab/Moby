/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/cblas.h>
#include <Moby/Spatial.h>

using namespace Ravelin;
using namespace Moby;
using std::vector;

/// Concates a vector with a wrench to make a new vector
VectorNd& concat(const VectorNd& v, const Wrenchd& w, VectorNd& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(v.size()+SPATIAL_DIM);
  result.set_sub_vec(0, result);
  result.set_sub_vec(v.size()+0, w.get_force());
  result.set_sub_vec(v.size()+3, w.get_torque());
  return result;
}

/// Multiplies a vector of spatial vectors by a matrix 
MatrixNd& mult(const vector<SVector6d>& v, const MatrixNd& m, MatrixNd& result)
{
  const unsigned SPATIAL_DIM = 6;

  if (v.size() != m.rows())
    throw MissizeException();

  // setup the result
  result.set_zero(SPATIAL_DIM, m.columns());

  // if the vector is empty, return now
  if (v.empty())
    return result;

  // verify that all twists are in the same pose
  for (unsigned i=1; i< v.size(); i++)
    if (v[i].pose != v[i-1].pose)
      throw FrameException(); 

  // finally, do the computation
  for (unsigned k=0; k< m.columns(); k++)
  {
    const double* mdata = m.column(k).data();
    for (unsigned j=0; j< SPATIAL_DIM; j++)
      for (unsigned i=0; i< v.size(); i++)
        result(k,j) += v[i][j]*mdata[i];
  }

  return result;
}

/// Multiplies a vector of spatial vectors by a vector
SVector6d mult(const vector<SVector6d>& t, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (t.size() != v.size())
    throw MissizeException();

  // setup the result
  SVector6d result = SVector6d::zero();

  // if the vector is empty, return now
  if (t.empty())
    return result;

  // verify that all twists are in the same pose
  result.pose = t.front().pose;
  for (unsigned i=1; i< t.size(); i++)
    if (t[i].pose != result.pose)
      throw FrameException(); 

  // finally, do the computation
  const double* vdata = v.data();
  for (unsigned j=0; j< SPATIAL_DIM; j++)
    for (unsigned i=0; i< t.size(); i++)
      result[j] += t[i][j]*vdata[i];

  return result;
}
/// Multiplies a vector of wrenches by a vector
SVector6d mult(const vector<Wrenchd>& w, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (w.size() != v.size())
    throw MissizeException();

  // setup the result
  SVector6d result = SVector6d::zero();

  // if the vector is empty, return now
  if (w.empty())
    return result;

  // verify that all twists are in the same pose
  result.pose = w.front().pose;
  for (unsigned i=1; i< w.size(); i++)
    if (w[i].pose != result.pose)
      throw FrameException(); 

  // finally, do the computation
  const double* vdata = v.data();
  for (unsigned j=0; j< SPATIAL_DIM; j++)
    for (unsigned i=0; i< w.size(); i++)
      result[j] += w[i][j]*vdata[i];

  return result;
}

/// Multiplies a vector of twists by a vector
SVector6d mult(const vector<Twistd>& t, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (t.size() != v.size())
    throw MissizeException();

  // setup the result
  SVector6d result = SVector6d::zero();

  // if the vector is empty, return now
  if (t.empty())
    return result;

  // verify that all twists are in the same pose
  result.pose = t.front().pose;
  for (unsigned i=1; i< t.size(); i++)
    if (t[i].pose != result.pose)
      throw FrameException(); 

  // finally, do the computation
  const double* vdata = v.data();
  for (unsigned j=0; j< SPATIAL_DIM; j++)
    for (unsigned i=0; i< t.size(); i++)
      result[j] += t[i][j]*vdata[i];

  return result;
}

Twistd spatial_cross(const Twistd& v1, const Twistd& v2)
{
  return SVector6d::spatial_cross(v1, v2);
}

