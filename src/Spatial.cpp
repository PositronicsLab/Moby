/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

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

/// Multiplies a vector of twists by a vector
Twistd mult(const vector<Twistd>& t, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (t.size() != v.size())
    throw MissizeException();

  // setup the result
  Twistd result = Twistd::zero();

  // if the vector is empty, return now
  if (t.empty())
    return result;

  // verify that all twists are in the same pose
  result.pose = t.front().pose;
  for (unsigned i=1; i< t.size(); i++)
    if (t[i].pose != result.pose)
      throw FrameException(); 

  // finally, do the computation
  for (unsigned j=0; j< SPATIAL_DIM; j++)
    for (unsigned i=0; i< t.size(); i++)
      result[j] += t[i][j]*v[i];

  return result;
}

Twistd spatial_cross(const Twistd& v1, const Twistd& v2)
{
  return SVector6d::spatial_cross(v1, v2);
}

