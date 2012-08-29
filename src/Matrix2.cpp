/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
#include <cstring>
#include <iostream>
#include <cmath>
#include <Moby/Constants.h>
#include <Moby/Matrix2.h>

using namespace Moby;

/// Constructs a matrix from an array
/**
 * \param array an array of 4 Real values in row-major format
 */
Matrix2::Matrix2(const Real* array)
{
  for (unsigned i=0; i< 4; i++)
    _data[i] = array[i];
}

/// Constructs a matrix from 4 values
/**
 * The resulting matrix will appear as follows: <br />
 * m00 m01 <br />
 * m10 m11 
 */
Matrix2::Matrix2(Real m00, Real m01, Real m10, Real m11)
{
  const unsigned X = 0, Y = 1;
  operator()(X,X) = m00;
  operator()(X,Y) = m01;
  operator()(Y,X) = m10;
  operator()(Y,Y) = m11;
}

/// Inverts this matrix 
void Matrix2::inverse()
{
  *this = inverse(*this);
}

/// Calculates the determinant for a 2x2 matrix
Real Matrix2::det() const
{
  return _data[0]*_data[3] - _data[1]*_data[2];
}

/// Determines whether this is an orthonormal matrix
bool Matrix2::is_orthonormal() const
{
  Real determinant = det();
  return (std::fabs(determinant - 1.0) < NEAR_ZERO || std::fabs(determinant + 1.0) < NEAR_ZERO);
}

/// Determines the inverse of the given matrix
Matrix2 Matrix2::inverse(const Matrix2& m)
{
  // compute the determinant
  Real dt = 1.0/m.det();
  return Matrix2(m._data[3]*dt, -m._data[2]*dt, -m._data[1]*dt, m._data[0]*dt);
}

/// Sets this matrix to the rotation matrix of the specified angle around the Z axis
void Matrix2::set_rot_Z(Real angle)
{
  *this = rot_Z(angle);
}

/// Returns the rotation matrix of the specified angle around the Z axis
Matrix2 Matrix2::rot_Z(Real angle)
{
  Real sina = std::sin(angle);
  Real cosa = std::cos(angle);
  Matrix2 r;
  r(0,0) = cosa;
  r(0,1) = -sina;
  r(1,0) = sina;
  r(1,1) = cosa;
  return r;
}

/// Graham-Schmidt orthonormalization
bool Matrix2::orthonormalize(Vector2& a, Vector2& b)
{
  // Gram-Schmidt orthonormalization produces vectors u0, u1, and u2 as
  // follows:
  //
  //   u0 = v0/|v0|
  //   u1 = (v1 - Dot(u0,v1)*u0)/|v1 - Dot(u0,v1)*u0|

  // Compute u0.
  Real flen = a.norm();
  if (flen == 0.0)
    return false;

  // Compute u1.
  Real fd0 = Vector2::dot(a, b);
  b -= fd0*a;
  flen = b.norm();
  return (flen != 0.0);
}

/// Makes the matrix orthonormal using Gram-Schmidt orthogonalization
bool Matrix2::orthonormalize()
{
  const unsigned X = 0, Y = 1;

  // The algorithm uses Gram-Schmidt orthogonalization applied to the
  // columns of M.
  Vector2 a(operator()(X,X),operator()(X,Y));
  Vector2 b(operator()(Y,X),operator()(Y,Y));
  if (orthonormalize(a,b))
  {
    operator()(X,X) = a[X];
    operator()(X,Y) = a[Y];
    operator()(Y,X) = b[X];
    operator()(Y,Y) = b[Y];
    return true;
  }
  
  return false;
}

/// Sets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param v a 3D vector
 * \note the number of columns of this must be three!
 */
void Matrix2::set_row(unsigned i, const Vector2& v)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  _data[st] = v[0];
  st += SZ;
  _data[st] = v[1]; 
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
Vector2 Matrix2::get_row(unsigned i) const
{
  Vector2 v;
  get_row(i, v);
  return v;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param target an array of size columns() or greater
 */
void Matrix2::get_row(unsigned i, Real* target) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i;
  target[0] = _data[st];
  st += SZ;
  target[1] = _data[st];
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
void Matrix2::get_row(unsigned i, Vector2& result) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i;
  result[0] = _data[st];
  st += SZ;
  result[1] = _data[st];
}

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param v a 3D vector 
 */
void Matrix2::set_column(unsigned i, const Vector2& v)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  _data[st++] = v[0];
  _data[st] = v[1];
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 */
Vector2 Matrix2::get_column(unsigned i) const
{
  Vector2 v;
  get_column(i, v);
  return v;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param target an array of size rows() or greater
 */
void Matrix2::get_column(unsigned i, Real* target) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i*SZ;
  target[0] = _data[st++];
  target[1] = _data[st];
}

/// Gets the specified column of the matrix
/*
 * \param i the 0-index of the column
 */
void Matrix2::get_column(unsigned i, Vector2& result) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i*SZ;
  result[0] = _data[st++];
  result[1] = _data[st];
}

/// Sets this matrix to its transpose
void Matrix2::transpose()
{
  const unsigned YX = 1, XY = 3;
  std::swap(_data[XY], _data[YX]);
}

/// Determines the transpose of this matrix
Matrix2 Matrix2::transpose(const Matrix2& m)
{
  Matrix2 n = m;
  n.transpose();
  return n;
}

/// Copies a matrix to this one
Matrix2& Matrix2::operator=(const Matrix2& m)
{
  const unsigned SZ = 4;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m.begin()[i];
  return *this;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
Vector2 Matrix2::mult(const Vector2& v) const
{
  Vector2 result;
  result[0] = _data[0]*v[0] + _data[2]*v[1];
  result[1] = _data[1]*v[0] + _data[3]*v[1];
  return result;
}

/// Multiplies the transpose of this matrix by a vector and returns the result in a new vector
Vector2 Matrix2::transpose_mult(const Vector2& v) const
{
  Vector2 result;
  result[0] = _data[0]*v[0] + _data[1]*v[1];
  result[1] = _data[2]*v[0] + _data[3]*v[1];
  return result;
}

/// Multiplies the transpose of this matrix by a matrix and returns the result in a new matrix 
Matrix2 Matrix2::transpose_mult(const Matrix2& m) const
{
  Matrix2 result;
  result._data[0] = _data[0]*m._data[0] + _data[1]*m._data[1];
  result._data[1] = _data[2]*m._data[0] + _data[3]*m._data[1];
  result._data[2] = _data[0]*m._data[2] + _data[1]*m._data[3];
  result._data[3] = _data[2]*m._data[2] + _data[3]*m._data[3];
  return result;
}

/// Multiplies the transpose of this matrix by the transpose of a matrix and returns the result in a new matrix 
Matrix2 Matrix2::transpose_mult_transpose(const Matrix2& m) const
{
  Matrix2 result;
  result._data[0] = _data[0]*m._data[0] + _data[1]*m._data[2];
  result._data[1] = _data[2]*m._data[0] + _data[3]*m._data[2];
  result._data[2] = _data[0]*m._data[1] + _data[1]*m._data[3];
  result._data[3] = _data[2]*m._data[1] + _data[3]*m._data[3];
  return result;
}

/// Multiplies this matrix by the transpose of a matrix and returns the result in a new matrix 
Matrix2 Matrix2::mult_transpose(const Matrix2& m) const
{
  Matrix2 result;
  result._data[0] = _data[0]*m._data[0] + _data[2]*m._data[2];
  result._data[1] = _data[1]*m._data[0] + _data[3]*m._data[2];
  result._data[2] = _data[0]*m._data[1] + _data[2]*m._data[3];
  result._data[3] = _data[1]*m._data[1] + _data[3]*m._data[3];
  return result;
}

/// Multiplies this matrix by another 3x3 matrix
Matrix2 Matrix2::mult(const Matrix2& m) const
{
  Matrix2 result;
  result._data[0] = _data[0]*m._data[0] + _data[2]*m._data[1];
  result._data[1] = _data[1]*m._data[0] + _data[3]*m._data[1];
  result._data[2] = _data[0]*m._data[2] + _data[2]*m._data[3];
  result._data[3] = _data[1]*m._data[2] + _data[3]*m._data[3];
  return result;
}

/// Multiplies this matrix by a scalar in place
Matrix2& Matrix2::operator*=(Real scalar)
{
  _data[0] *= scalar;
  _data[1] *= scalar;
  _data[2] *= scalar;
  _data[3] *= scalar;
  return *this;
}

/// Returns the negation of this matrix
Matrix2 Matrix2::operator-() const
{
  Matrix2 m;
  std::transform(begin(), end(), m.begin(), std::negate<Real>());
  return m;
}

/// Adds m to this in place
Matrix2& Matrix2::operator+=(const Matrix2& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::plus<Real>());
  return *this;
}

/// Subtracts m from this in place
Matrix2& Matrix2::operator-=(const Matrix2& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::minus<Real>());
  return *this;
}

/// Sets this matrix to identity
void Matrix2::set_identity()
{
  const unsigned XX = 0, YY = 3;
  set_zero();
  _data[XX] = _data[YY] = 1.0;
}

/// Sets this matrix to zero
void Matrix2::set_zero()
{
  const unsigned N = 2*2;
  for (unsigned i=0; i< N; i++)
    _data[i] = (Real) 0.0;
}

/// Checks whether this matrix is symmetric
bool Matrix2::is_symmetric(Real tolerance) const
{
  const unsigned X = 0, Y = 1;

  // make sure tolerance is non-negative
  tolerance = std::fabs(tolerance);

  // check symmetric
  if (!rel_equal(operator()(X,Y), operator()(Y,X), tolerance))
    return false;

  return true;
}

/// Determines whether two numbers are relatively equal
bool Matrix2::rel_equal(Real x, Real y, Real tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), (Real) 1.0)));
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const Matrix2& m)
{
  const unsigned SZ = 2; 
  for (unsigned i=0; i< SZ; i++)
  {
    for (unsigned j=0; j< SZ-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,SZ-1) << std::endl;
  }
   
  return out;
}

