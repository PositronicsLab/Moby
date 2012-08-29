/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/cblas.h>
#include <cstring>
#include <functional>
#include <algorithm>
#include <Moby/Matrix3.h>
#include <Moby/SVector6.h>

using namespace Moby;

/// Constructs this vector with the given values
SVector6::SVector6(Real x, Real y, Real z, Real a, Real b, Real c)
{
  _data[0] = x;
  _data[1] = y;
  _data[2] = z;
  _data[3] = a;
  _data[4] = b;
  _data[5] = c;
}

/// Constructs this vector from the given array
/**
 * \param array a 6-dimensional (or larger) array
 */
SVector6::SVector6(const Real* array)
{
  for (unsigned i=0; i< 6; i++)
    _data[i] = array[i];
}

/// Constructs the given spatial vector with given upper and lower components
SVector6::SVector6(const Vector3& upper, const Vector3& lower)
{
  set_upper(upper);
  set_lower(lower);
}

/// Constructs the given spatial vector using a VectorN
SVector6::SVector6(const VectorN& v)
{
  if (v.size() != 6)
    throw MissizeException();
  for (unsigned i=0; i< 6; i++)
    _data[i] = v.data()[i];
}

/// Computes the spatial cross product between two vectors
SVector6 SVector6::spatial_cross(const SVector6& v1, const SVector6& v2)
{
  // get the two skew symmetric matrices we need
  Matrix3 ax = Matrix3::skew_symmetric(v1.get_upper());
  Matrix3 bx = Matrix3::skew_symmetric(v1.get_lower());

  // multiply
  Vector3 v2top = v2.get_upper();
  Vector3 v2bot = v2.get_lower();
  Vector3 top = ax * v2top;
  Vector3 bot = (bx * v2top) + (ax * v2bot);
  return SVector6(top, bot);
}

/// Gets the lower 3-dimensional vector
Vector3 SVector6::get_lower() const
{
  return Vector3(_data[3], _data[4], _data[5]);
}

/// Gets the upper 3-dimensional vector
Vector3 SVector6::get_upper() const
{
  return Vector3(_data[0], _data[1], _data[2]);
}
 
/// Sets the lower 3-dimensional vector
void SVector6::set_lower(const Vector3& lower)
{
  _data[3] = lower[0];
  _data[4] = lower[1];
  _data[5] = lower[2];
}

/// Sets the upper 3-dimensional vector
void SVector6::set_upper(const Vector3& upper)
{
  _data[0] = upper[0];
  _data[1] = upper[1];
  _data[2] = upper[2];
}

/// Performs the dot product
Real SVector6::dot(const SVector6& v1, const SVector6& v2)
{
  return v1[3]*v2[0] + v1[4]*v2[1] + v1[5]*v2[2] + v1[0]*v2[3] + v1[1]*v2[4] + v1[2]*v2[5];
}

/// Sets this vector to its spatial vector transpose
void SVector6::transpose()
{
  std::swap(_data[0], _data[3]);
  std::swap(_data[1], _data[4]);
  std::swap(_data[2], _data[5]);
}

/// Returns the spatial transpose of the given vector
SVector6 SVector6::transpose(const SVector6& v)
{
  return SVector6(v.get_lower(), v.get_upper());
}

/// Copies this vector from another SVector6
SVector6& SVector6::operator=(const SVector6& v)
{
  const unsigned SZ = 6;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = v.data()[i];
  return *this;
}

/// Returns the negation of this vector
SVector6 SVector6::operator-() const
{
  SVector6 v;
  std::transform(begin(), end(), v.begin(), std::negate<Real>());

  return v;
}

/// Multiplies this vector by a scalar in place
SVector6& SVector6::operator*=(Real scalar)
{
  for (unsigned i=0; i< 6; i++)
    _data[i] *= scalar;

  return *this;
}

/// Adds another vector to this one in place
SVector6& SVector6::operator+=(const SVector6& v)
{
  std::transform(begin(), end(), v.begin(), begin(), std::plus<Real>());

  return *this;
}

/// Subtracts another vector from this one in place
SVector6& SVector6::operator-=(const SVector6& v)
{
  std::transform(begin(), end(), v.begin(), begin(), std::minus<Real>());

  return *this;
}

/// Adds this vector to another and returns the result in a new vector
SVector6 SVector6::operator+(const SVector6& v) const
{
  SVector6 result;
  std::transform(begin(), end(), v.begin(), result.begin(), std::plus<Real>());
  return result;
}

/// Subtracts another vector from this vector and returns the result in a new vector
SVector6 SVector6::operator-(const SVector6& v) const
{
  SVector6 result;
  std::transform(begin(), end(), v.begin(), result.begin(), std::minus<Real>());
  return result;
}

