/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/cblas.h>
#include <Moby/Spatial.h>

using namespace Ravelin;
using std::vector;

namespace Moby {

/// Concates a vector with a force to make a new vector
VectorNd& concat(const VectorNd& v, const SForced& w, VectorNd& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(v.size()+SPATIAL_DIM);
  result.set_sub_vec(0, result);
  result.set_sub_vec(v.size()+0, w.get_force());
  result.set_sub_vec(v.size()+3, w.get_torque());
  return result;
}

/// Concates a vector with a momentum to make a new vector
VectorNd& concat(const VectorNd& v, const SMomentumd& w, VectorNd& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(v.size()+SPATIAL_DIM);
  result.set_sub_vec(0, result);
  result.set_sub_vec(v.size()+0, w.get_linear());
  result.set_sub_vec(v.size()+3, w.get_angular());
  return result;
}

/// Multiplies a vector of spatial momenta by a vector 
VectorNd& mult(const vector<SMomentumd>& Is, const VectorNd& v, VectorNd& result)
{
  const unsigned SPATIAL_DIM = 6;

  if (Is.size() != v.rows())
    throw MissizeException();

  // setup the result
  result.set_zero(SPATIAL_DIM);

  // if the vector is empty, return now
  if (Is.empty())
    return result;

  // verify that all twists are in the same pose
  for (unsigned i=1; i< Is.size(); i++)
    if (Is[i].pose != Is[i-1].pose)
      throw FrameException(); 

  // finally, do the computation
  const double* vdata = v.data();
  double* rdata = result.data();
  for (unsigned j=0; j< SPATIAL_DIM; j++)
    for (unsigned i=0; i< Is.size(); i++)
    {
      const double* Isdata = Is[i].data();
      rdata[j] += Isdata[j]*vdata[i];
    }

  return result;
}

/// Multiplies a vector of spatial momenta by a matrix 
MatrixNd& mult(const vector<SMomentumd>& Is, const MatrixNd& m, MatrixNd& result)
{
  const unsigned SPATIAL_DIM = 6;

  if (Is.size() != m.rows())
    throw MissizeException();

  // setup the result
  result.set_zero(SPATIAL_DIM, m.columns());

  // if the vector is empty, return now
  if (Is.empty())
    return result;

  // verify that all twists are in the same pose
  for (unsigned i=1; i< Is.size(); i++)
    if (Is[i].pose != Is[i-1].pose)
      throw FrameException(); 

  // finally, do the computation
  double* rdata = result.data();
  for (unsigned k=0; k< m.columns(); k++)
  {
    const double* mdata = m.column(k).data();
    for (unsigned j=0; j< SPATIAL_DIM; j++)
      for (unsigned i=0; i< Is.size(); i++)
      {
        const double* Isdata = Is[i].data();
        rdata[k+j*result.rows()] += Isdata[j]*mdata[i];
      }
  }

  return result;
}

/// Multiplies a spatial inertia by a vector of spatial axes
MatrixNd& mult(const SpatialABInertiad& I, const std::vector<SAxisd>& s, MatrixNd& result)
{
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(SPATIAL_DIM, s.size());

  // compute the individual momenta
  for (unsigned i=0; i< s.size(); i++)
  {
    SMomentumd m = I.mult(s[i]);
    SharedVectorNd col = result.column(i);
    m.to_vector(col);
  } 

  return result;
}

/// Multiplies a spatial inertia by a vector of spatial axes
vector<SMomentumd>& mult(const SpatialABInertiad& I, const std::vector<SAxisd>& s, vector<SMomentumd>& result)
{
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(s.size());

  // compute the individual momenta
  for (unsigned i=0; i< s.size(); i++)
    result[i] = I.mult(s[i]);

  return result;
}

/// Multiplies a vector of spatial axes by a vector
SVelocityd mult(const vector<SAxisd>& t, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (t.size() != v.size())
    throw MissizeException();

  // setup the result
  SVelocityd result = SVelocityd::zero();

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

} // end namespace

