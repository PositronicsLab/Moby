/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <Moby/Constants.h>
#include <Moby/MatrixNN.h>
#include <Moby/Vector3.h>

using namespace Moby;
using namespace boost;

/// Constructs this vector with the given values
Vector3::Vector3(Real x, Real y, Real z)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
}

/// Constructs this vector from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
Vector3::Vector3(const Real* array)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
}

/// Determines whether all components of this vector are finite
bool Vector3::is_finite() const
{
  return !std::isinf(_data[0]) && !std::isinf(_data[1]) && !std::isinf(_data[2]);
}

/// Swaps the contents of this vector with another
void Vector3::swap(Vector3& source)
{
  const unsigned X = 0, Y = 1, Z = 2;
  std::swap(_data[X], source[X]);
  std::swap(_data[Y], source[Y]);
  std::swap(_data[Z], source[Z]);
}

/// Computes the outer product of two Vector3 objects
Matrix3* Vector3::outer_prod(const Vector3& v1, const Vector3& v2, Matrix3* M)
{
  const unsigned X = 0, Y = 1, Z = 2;
  (*M)(X,X) = v1[X]*v2[X];
  (*M)(X,Y) = v1[X]*v2[Y];
  (*M)(X,Z) = v1[X]*v2[Z];
  (*M)(Y,X) = v1[Y]*v2[X];
  (*M)(Y,Y) = v1[Y]*v2[Y];
  (*M)(Y,Z) = v1[Y]*v2[Z];
  (*M)(Z,X) = v1[Z]*v2[X];
  (*M)(Z,Y) = v1[Z]*v2[Y];
  (*M)(Z,Z) = v1[Z]*v2[Z];
  return M;
}

/// Computes the cross-product of two vectors
Vector3 Vector3::cross(const Vector3& v1, const Vector3& v2)
{
  Vector3 w;
  w[0] = v1[1] * v2[2] - v1[2] * v2[1];
  w[1] = v1[2] * v2[0] - v1[0] * v2[2];
  w[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return w;
}

/// Determines a vector orthogonal to v
/**
 * \note the returned vector is not normalized
 */ 
Vector3 Vector3::determine_orthogonal_vec(const Vector3& v)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // make a vector of all ones
  Vector3 ones(1,1,1);
      
  // get the absolute values of the three components
  Real x = std::fabs(v[X]);
  Real y = std::fabs(v[Y]);
  Real z = std::fabs(v[Z]);
    
  // make the component zero that is equal to the largest component of v1 
  if (x > y)
  {
    if (x > z)
      ones[X] = 0;
    else
      ones[Z] = 0;
  }
  else
  {
    if (y > z)
      ones[Y] = 0;
    else
      ones[Z] = 0;
  }
      
  // compute the cross product of v and the ones vector
  return Vector3::cross(v, ones);
}

/// Computes an orthonormal basis, given a single vector
/**
 * \return an orthonormal basis constructed from a single 3D vector, v1; v1 x v2 = v3
 */
void Vector3::determine_orthonormal_basis(const Vector3& v1, Vector3& v2, Vector3& v3)
{
  // get a vector orthogonal to v1
  v2 = Vector3::normalize(determine_orthogonal_vec(v1));
  
  // compute the second vector
  v3 = Vector3::normalize(Vector3::cross(v1, v2));
}

/// Determines whether two numbers are relatively equal
bool Vector3::rel_equal(Real x, Real y)
{
  return (std::fabs(x-y) <= NEAR_ZERO * std::max(std::fabs(x), std::max(std::fabs(y), (Real) 1.0)));
}

/// Compares the two vectors lexographically
bool Vector3::operator<(const Vector3& v) const
{
  const unsigned LEN = 3;
  for (unsigned i=0; i< LEN; i++)
  {
    if (rel_equal(_data[i], v[i]))
      continue;
    return _data[i] < v[i];
  }

  // still here?  comparison was identical
  return false;
}

