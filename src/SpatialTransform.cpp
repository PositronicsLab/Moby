/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Constants.h>
#include <Moby/SpatialTransform.h>

using namespace Moby;

/// Default constructor -- constructs an identity transformation 
SpatialTransform::SpatialTransform()
{
  E = IDENTITY_3x3;
  r = ZEROS_3;
}

/// Computes the spatial transform from a global frame
SpatialTransform SpatialTransform::from_global(const Matrix4& T)
{
  Matrix3 R;

  // get the two rotation matrices that we will need
  T.get_rotation(&R);
  R.transpose();

  // setup the spatial transform
  SpatialTransform X;
  X.r = R.mult(T.get_translation());
  X.E = R;
  return X;
}

/// Computes the spatial transform to a global frame
SpatialTransform SpatialTransform::to_global(const Matrix4& T)
{
  Matrix3 R;

  // get the two rotation matrices that we will need
  T.get_rotation(&R);

  // setup the spatial transform
  SpatialTransform X;
  X.r = -T.get_translation();
  X.E = R;
  return X;
}

/// Calculates the spatial transform given a 3x3 orientation matrix, a 3-dimensional translation vector, and one 4x4 transformation matrix
SpatialTransform::SpatialTransform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix4& target)
{
  Matrix3 Rtarget;

  // get the two rotation matrices that we will need
  target.get_rotation(&Rtarget);

  // setup the spatial transform
  r = Rtarget.transpose_mult(target.get_translation() - sourceX);
  E = Rtarget.transpose_mult(sourceR);
}

/// Calculates the spatial transform given a 4x4 transformation matrix, a 3x3 orientation matrix, and a 3-dimensional translation vector
SpatialTransform::SpatialTransform(const Matrix4& source, const Matrix3& targetR, const Vector3& targetX)
{
  Matrix3 Rsource;

  // get the rotation matrix that we will need
  source.get_rotation(&Rsource);

  // setup the spatial transform
  r = targetR.transpose_mult(targetX - source.get_translation());
  E = targetR.transpose_mult(Rsource);
}

/// Calculates the spatial transform given two 4x4 transformation matrices
SpatialTransform::SpatialTransform(const Matrix4& source, const Matrix4& target)
{
  Matrix3 Rsource, Rtarget;

  // get the two rotation matrices that we will need
  source.get_rotation(&Rsource);
  target.get_rotation(&Rtarget);

  // setup the spatial transform
  r = Rtarget.transpose_mult(target.get_translation() - source.get_translation());
  E = Rtarget.transpose_mult(Rsource);
}

/// Calculates the spatial transform given two 3x3 rotation matrices and two translation vectors
SpatialTransform::SpatialTransform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix3& targetR, const Vector3& targetX)
{
  // setup the matrix
  r = targetR.transpose_mult(targetX - sourceX);
  E = targetR.transpose_mult(sourceR);
}

/// Copies a spatial matrix to this one 
SpatialTransform& SpatialTransform::operator=(const SpatialTransform& m)
{
  this->r = m.r;
  this->E = m.E;
  return *this;
}

/// Transforms a rigid body inertia in its frame
SpatialRBInertia SpatialTransform::transform(Real mass, const Matrix3& J) const
{
  // precompute some things we'll need
  SpatialRBInertia result;
  Vector3 mr = r * mass;
  Matrix3 rx = Matrix3::skew_symmetric(r);
  Matrix3 mrxrx = rx*Matrix3::skew_symmetric(mr);
  Matrix3 ET = Matrix3::transpose(E);

  // compute the result
  result.m = mass;
  result.h = -mr;
  result.J = (E*J*ET) - mrxrx;

  return result;
}

/// Transforms the spatial matrix
SMatrix6N SpatialTransform::transform(const SMatrix6N& m) const
{
  SMatrix6N result;
  transform(m, result);
  return result;
}

/// Transforms a spatial rigid body inertia
SpatialRBInertia SpatialTransform::transform(const SpatialRBInertia& m) const
{
  // precompute some things we'll need
  SpatialRBInertia result;
  Vector3 mr = r * m.m;
  Matrix3 rx = Matrix3::skew_symmetric(r);
  Matrix3 hx = Matrix3::skew_symmetric(m.h);
  Matrix3 mrxrx = rx*Matrix3::skew_symmetric(mr);
  Matrix3 ET = Matrix3::transpose(E);
  Matrix3 EhxETrx = E * hx * ET * rx;

  // compute the result
  result.m = m.m;
  result.h = E * m.h - mr;
  result.J = EhxETrx + Matrix3::transpose(EhxETrx) + (E*m.J*ET) - mrxrx;

  return result;
}

/// Transforms a spatial articulated body inertia 
SpatialABInertia SpatialTransform::transform(const SpatialABInertia& m) const
{
  // precompute some things we'll need
  Matrix3 ET = Matrix3::transpose(E);
  Matrix3 rx = Matrix3::skew_symmetric(r);
  Matrix3 HT = Matrix3::transpose(m.H);
  Matrix3 EJET = E * m.J * ET;
  Matrix3 rx_E_HT_ET = rx*E*HT*ET;
  Matrix3 EHET = E * m.H * ET;
  Matrix3 EMET = E * m.M * ET;
  Matrix3 rxEMET = rx * EMET;

  SpatialABInertia result;
  result.M = EMET;
  result.H = EHET - rxEMET;
  result.J = EJET - rx_E_HT_ET + ((EHET - rxEMET) * rx); 
  return result;
}

/// Transforms the spatial vector 
SVector6 SpatialTransform::transform(const SVector6& v) const
{
  // get the components of v
  Vector3 top = v.get_upper();
  Vector3 bottom = v.get_lower();

  // do the calculations
  Vector3 Etop = E * top;
  Vector3 cross = Vector3::cross(r, Etop);
  return SVector6(Etop, (E * bottom) - cross);
}

/// Combines (concatenates) two spatial transforms 
SpatialTransform SpatialTransform::concatenate(const SpatialTransform& m) const
{
  SpatialTransform result;
  result.E = E * m.E;
  result.r = r + (E * m.r);
  return result;
}

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
SMatrix6N& SpatialTransform::transform(const SMatrix6N& m, SMatrix6N& result) const
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of columns of m
  const unsigned NCOLS = m.columns();

  // look for empty result
  if (NCOLS == 0)
  {
    result.set_zero(SPATIAL_DIM, NCOLS);
    return result;
  }

  // resize the new matrix
  result.resize(SPATIAL_DIM, NCOLS);
  
  // do multiplication on each column of m
  for (unsigned i=0; i< NCOLS; i++)
  {
    SVector6 v = m.get_column(i);
    v = transform(v);
    result.set_column(i, v);
  }

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const SpatialTransform& m) 
{
  out << "spatial transform R: " << std::endl << m.E;
  out << "spatial transform r: " << " r: " << m.r << std::endl;
  out << "spatial transform as 6x6: " << std::endl;
  // top 3x6 matrix
  for (unsigned i=0; i< 3; i++)
  {
    for (unsigned j=0; j< 3; j++)
      out << m.E(i,j) << " ";
    out << "0 0 0" << std::endl;
  }

  // bottom 3x6 matrix
  Matrix3 X = Matrix3::skew_symmetric(m.r) * m.E;
  for (unsigned i=0; i< 3; i++)
  {
    for (unsigned j=0; j< 3; j++)
      out << X(i,j) << " ";
    for (unsigned j=0; j< 3; j++)
      out << " " << m.E(i,j);
    out << std::endl;
  }

  return out;
}

