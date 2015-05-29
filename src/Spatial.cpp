/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/cblas.h>
#include <Ravelin/Transform3d.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Spatial.h>

using namespace Ravelin;
using std::vector;

namespace Moby {

/// Duplication of Pose3d::get_r_E(.)
/// Gets r and E from a transform 
/**
 * \param T a transformation from frame A to frame B
 * \param r on return, the vector from A's origin to B's origin in A frame 
 * \param E rotates vectors in A's orientation to vectors in B's orientation 
 */ 
static void get_r_E(const Transform3d& T, Vector3d& r, Matrix3d& E)
{
  // x is the translation from frame A to frame B

  // p_a
  // bTa * p_a = bQa * p_a + bXa

  // note that x is translation from relative pose to this pose
  // q is rotation from vectors in this pose to relative pose 
  E = T.q;
  r = Vector3d(E.transpose_mult(-T.x), T.source);
}

/// transforms a spatial acceleration using precomputation *without accounting for moving frames*
static void transform_accel(boost::shared_ptr<const Pose3d> target, const SAcceld& w, const Vector3d& r, const Matrix3d& E, SAcceld& result)
{
  // get the components of w[i] 
  Vector3d top = w.get_upper();
  Vector3d bottom = w.get_lower();

  // do the calculations
  Vector3d Etop(E * Origin3d(top), target);
  Vector3d cross = Vector3d::cross(r, top);
  result.set_upper(Etop);
  result.set_lower(Vector3d(E * Origin3d(bottom - cross), target));
  result.pose = target;
} 

/// Special transformation of acceleration when moving aspect of pose accounted for elsewhere
SAcceld transform_accel(boost::shared_ptr<const Pose3d> target, const SAcceld& a)
{
  SAcceld s;

  // NOTE: this is a duplication of the transform_spatial(.) function in
  // Ravelin::Pose3x
  // setup the source pose 
  boost::shared_ptr<const Pose3d> source = a.pose;

  // quick check
  if (source == target)
  {
    s=a;
    return s;
  }

  // compute the relative transform
  Transform3d Tx = Pose3d::calc_relative_pose(source, target);

  // setup r and E
  Vector3d r;
  Matrix3d E;
  get_r_E(Tx, r, E);

  // get the components of a
  Vector3d top = a.get_upper();
  Vector3d bottom = a.get_lower();

  // do the calculations
  Vector3d Etop(E * Origin3d(top), target);
  Vector3d cross = Vector3d::cross(r, top);
  s.set_upper(Etop);
  s.set_lower(Vector3d(E * Origin3d(bottom - cross), target));
  s.pose = target;
  return s;
}

/// Special transformation of acceleration when moving aspect of pose accounted for elsewhere
std::vector<SAcceld>& transform_accel(boost::shared_ptr<const Pose3d> target, const std::vector<SAcceld>& t, std::vector<SAcceld>& result)
{
  // NOTE: this is a duplication of the transform_spatial(.) function in
  // Ravelin::Pose3x

  // look for empty vector (easy case)
  if (t.empty())
  {
    result.clear();
    return result;
  }

  // setup the source pose
  boost::shared_ptr<const Pose3d> source = t[0].pose; 

  #ifndef NEXCEPT
  for (unsigned i=1; i< t.size(); i++)
    if (source != t[i].pose)
      throw FrameException();
  #endif

  // quick check
  if (source == target)
    return (result = t);

  // compute the relative transform
  Transform3d Tx = Pose3d::calc_relative_pose(source, target);

  // setup r and E
  Vector3d r;
  Matrix3d E;
  get_r_E(Tx, r, E);

  // resize the result vector
  result.resize(t.size());

  // transform the individual vectors 
  for (unsigned i=0; i< t.size(); i++)
    transform_accel(target, t[i], r, E, result[i]);

  return result;
}

/// Concates a vector with a force to make a new vector
VectorNd& concat(const VectorNd& v, const SForced& w, VectorNd& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(v.size()+SPATIAL_DIM);
  result.set_sub_vec(0, v);
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
    for (unsigned j=0; j< SPATIAL_DIM; j++)  // j is row index for Is
      for (unsigned i=0; i< Is.size(); i++)  // i is column index for Is
      {
        const double* Isdata = Is[i].data();
        rdata[k*result.rows()+j] += Isdata[j]*mdata[i];
      }
  }

  return result;
}

/// Multiplies a spatial inertia by a vector of spatial axes
MatrixNd& mult(const SpatialABInertiad& I, const std::vector<SVelocityd>& s, MatrixNd& result)
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
vector<SMomentumd>& mult(const SpatialABInertiad& I, const std::vector<SVelocityd>& s, vector<SMomentumd>& result)
{
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(s.size());

  // compute the individual momenta
  for (unsigned i=0; i< s.size(); i++)
    result[i] = I.mult(s[i]);

  return result;
}

/// Multiplies a spatial inertia by a vector of spatial axes
MatrixNd& mult(const SpatialRBInertiad& I, const std::vector<SVelocityd>& s, MatrixNd& result)
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
vector<SMomentumd>& mult(const SpatialRBInertiad& I, const std::vector<SVelocityd>& s, vector<SMomentumd>& result)
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
SVelocityd mult(const vector<SVelocityd>& t, const VectorNd& v)
{
  const unsigned SPATIAL_DIM = 6;

  if (t.size() != v.size())
    throw MissizeException();

  // verify that the vector is not empty - we lose frame info!
  if (t.empty())
    throw std::runtime_error("loss of frame information");

  // setup the result
  SVelocityd result = SVelocityd::zero();

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

