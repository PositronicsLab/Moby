/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/cblas.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstring>
#include <Moby/Constants.h>
#include <Moby/SMatrix6N.h>
#include <Moby/SMatrix6.h>

using namespace Moby;

/// Default constructor -- constructs an identity matrix
SMatrix6::SMatrix6()
{
  set_identity();
}

/// Constructs the 6x6 spatial matrix from an array of 36 values
SMatrix6::SMatrix6(const Real* source)
{
  const unsigned SZ = 6*6;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = source[i];
}

/// Constructs the 6x6 spatial matrix from a 6x6 dimensional MatrixN
/**
 * \throws MissizeException if MatrixN is not 6x6
 */
SMatrix6::SMatrix6(const MatrixN& m)
{
  if (m.rows() != 6 && m.columns() != 6)
    throw MissizeException();

  const unsigned SZ = 6*6;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m.data()[i];
}

/// Sets this matrix to identity
void SMatrix6::set_identity()
{
  const Matrix3 IDENTITY_3x3 = Matrix3::identity();
  set_upper_left(IDENTITY_3x3);
  set_upper_right(ZEROS_3x3);
  set_lower_left(ZEROS_3x3);
  set_lower_right(IDENTITY_3x3);
}

/// Gets the upper left 3x3 matrix
Matrix3 SMatrix6::get_upper_left() const
{
  Matrix3 m;
  
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, _data+i*6, 1, m.begin()+i*3, 1);

  return m; 
}

/// Gets the upper right 3x3 matrix
Matrix3 SMatrix6::get_upper_right() const
{
  Matrix3 m;
  
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, _data+(i+3)*6, 1, m.begin()+i*3, 1);

  return m; 
}

/// Gets the lower left 3x3 matrix
Matrix3 SMatrix6::get_lower_left() const
{
  Matrix3 m;
  
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, _data+i*6+3, 1, m.begin()+i*3, 1);

  return m; 
}

/// Gets the lower right 3x3 matrix
Matrix3 SMatrix6::get_lower_right() const
{
  Matrix3 m;
  
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, _data+(i+3)*6+3, 1, m.begin()+i*3, 1);

  return m; 
}

/// Sets the 3x3 upper left matrix
void SMatrix6::set_upper_left(const Matrix3& m)
{
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, m.begin()+i*3, 1, _data+i*6, 1);
}

/// Sets the 3x3 upper right matrix
void SMatrix6::set_upper_right(const Matrix3& m)
{
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, m.begin()+i*3, 1, _data+(i+3)*6, 1);
}

/// Sets the 3x3 lower left matrix
void SMatrix6::set_lower_left(const Matrix3& m)
{
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, m.begin()+i*3, 1, _data+i*6+3, 1);
}

/// Sets the 3x3 lower right matrix
void SMatrix6::set_lower_right(const Matrix3& m)
{
  for (unsigned i=0; i< 3; i++)
    CBLAS::copy(3, m.begin()+i*3, 1, _data+(i+3)*6+3, 1);
}

/// Calculates the spatial transform given two 4x4 transformation matrices
SMatrix6 SMatrix6::calc_spatial_transform(const Matrix4& source, const Matrix4& target)
{
  SMatrix6 X;
  Matrix3 Rsource, Rtarget;

  // get the two rotation matrices that we will need
  source.get_rotation(&Rsource);
  target.get_rotation(&Rtarget);

  // setup the matrix
  Vector3 r = Rtarget.transpose_mult(source.get_translation() - target.get_translation());
  Matrix3 R = Rtarget.transpose_mult(Rsource);
  X.set_upper_left(R);
  X.set_upper_right(ZEROS_3x3);
  X.set_lower_left(Matrix3::skew_symmetric(r) * R);
  X.set_lower_right(R);

  return X;
}

/// Calculates the spatial transform given two 3x3 rotation matrices and two translation vectors
SMatrix6 SMatrix6::calc_spatial_transform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix3& targetR, const Vector3& targetX)
{
  SMatrix6 X;
  
  // setup the matrix
  Vector3 r = targetR.transpose_mult(sourceX - targetX);
  Matrix3 R = targetR.transpose_mult(sourceR);
  X.set_upper_right(ZEROS_3x3);
  X.set_upper_left(R);
  X.set_lower_right(R);
  X.set_lower_left(Matrix3::skew_symmetric(r) * R);
  
  return X;
}

/// Special routine for calculating the inverse of a spatial inertia matrix
/**
 * Note that this does not calculate the inverse of a general 6x6 matrix.   
 */
SMatrix6 SMatrix6::inverse_inertia(const SMatrix6& I)
{
  const unsigned X = 0, Y = 1, Z = 2;
  SMatrix6 invI;

  Matrix3 Binv = -Matrix3::inverse(I.get_upper_right());
  Matrix3 tmp = I.get_lower_right() * Binv;
  invI.set_upper_right(Matrix3::inverse(tmp * I.get_upper_left() + I.get_lower_left()));
  tmp = invI.get_upper_right() * I.get_lower_right();
  invI.set_upper_left(tmp * Binv);
  invI.set_lower_right(Matrix3::transpose(invI.get_upper_left()));
  tmp = I.get_upper_left() * invI.get_upper_left();
  tmp(X,X)-= 1.0;
  tmp(Y,Y)-= 1.0;
  tmp(Z,Z)-= 1.0;
  invI.set_lower_left(Binv * tmp);
  
  return invI;
}

/// Computes the spatial cross operator, as defined in [Featherstone, 87], p. 29
SMatrix6 SMatrix6::spatial_cross(const SVector6& v)
{
  SMatrix6 X;

  Matrix3 ax = Matrix3::skew_symmetric(v.get_upper());
  Matrix3 bx = Matrix3::skew_symmetric(v.get_lower());

  X.set_upper_left(ax);
  X.set_upper_right(ZEROS_3x3);
  X.set_lower_left(bx);
  X.set_lower_right(ax);

  return X;
}

/// Computes the spatial transpose of this matrix, as defined in [Featherstone, 87], p. 52
void SMatrix6::transpose()
{
  Matrix3 A = Matrix3::transpose(get_upper_left());
  set_upper_left(Matrix3::transpose(get_lower_right()));
  set_upper_right(Matrix3::transpose(get_upper_right()));
  set_lower_left(Matrix3::transpose(get_lower_left()));
  set_lower_right(A);
}

/// Computes the spatial transpose of a 6x6 matrix, as defined in [Featherstone, 87], p. 52
SMatrix6 SMatrix6::transpose(const SMatrix6& m)
{
  SMatrix6 n;
  n.set_upper_left(Matrix3::transpose(m.get_lower_right()));
  n.set_upper_right(Matrix3::transpose(m.get_upper_right()));
  n.set_lower_left(Matrix3::transpose(m.get_lower_left()));
  n.set_lower_right(Matrix3::transpose(m.get_upper_left()));
  return n;
}

/// Copies a spatial matrix to this one 
SMatrix6& SMatrix6::operator=(const SMatrix6& m)
{
  const unsigned SZ = 6*6;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m.begin()[i];
  return *this;
}

/// Creates a zero matrix
void SMatrix6::set_zero()
{
  const unsigned N = 6*6;
  for (unsigned i=0; i< N; i++)
    _data[i] = (Real) 0.0;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
SVector6 SMatrix6::operator*(const SVector6& v) const
{
  SVector6 result; 
  CBLAS::gemv(CblasNoTrans, *this, v, (Real) 1.0, (Real) 0.0, result); 
  return result;
}

/// Multiplies this matrix by another 6x6 matrix
SMatrix6 SMatrix6::operator*(const SMatrix6& m) const
{
  SMatrix6 result;
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, *this, m, (Real) 1.0, (Real) 0.0, result); 
  return result;
}

/// Multiplies this matrix by a scalar in place
SMatrix6& SMatrix6::operator*=(Real scalar)
{
  const unsigned SZ = 6*6;
  CBLAS::scal(SZ, scalar, begin(), 1);
  return *this;
}

/// Returns the negation of this matrix
SMatrix6 SMatrix6::operator-() const
{
  SMatrix6 m;
  std::transform(begin(), end(), m.begin(), std::negate<Real>());
  return m;
}

/// Adds two spatial matrices
SMatrix6 SMatrix6::operator+(const SMatrix6& m) const
{
  SMatrix6 result;
  std::transform(begin(), end(), m.begin(), result.begin(), std::plus<Real>());
  return result;
}

/// Subtracts two spatial matrices
SMatrix6 SMatrix6::operator-(const SMatrix6& m) const
{
  SMatrix6 result;
  std::transform(begin(), end(), m.begin(), result.begin(), std::minus<Real>());
  return result;
}

/// Adds m to this in place
SMatrix6& SMatrix6::operator+=(const SMatrix6& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::plus<Real>());
  return *this;
}

/// Subtracts m from this in place
SMatrix6& SMatrix6::operator-=(const SMatrix6& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::minus<Real>());
  return *this;
}

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
SMatrix6N* SMatrix6::mult(const SMatrix6N* m, SMatrix6N* result) const
{
  const unsigned COLUMNS = 6, ROWS = 6;

  if (COLUMNS != m->rows())
    throw MissizeException();

  // look for empty result
  if (ROWS == 0 || COLUMNS == 0 || m->columns() == 0)
  {
    result->resize(ROWS, m->columns());
    result->set_zero();
    return result;
  }

  // resize the new matrix
  result->resize(ROWS, m->columns());

  // carry out multiplication using BLAS
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, *this, *m, (Real) 1.0, (Real) 0.0, *result); 

  return result;
}

/// Multiplies a spatial matrix by a standard matrix and returns the result in a standard matrix
MatrixN& SMatrix6::mult(const SMatrix6N* m, MatrixN& result) const
{
  const unsigned COLUMNS = 6, ROWS = 6;

  if (COLUMNS != m->rows())
    throw MissizeException();

  // look for empty result
  if (ROWS == 0 || COLUMNS == 0 || m->columns() == 0)
  {
    result.resize(ROWS, m->columns());
    result.set_zero();
    return result;
  }

  // resize the new matrix
  result.resize(ROWS, m->columns());

  // carry out multiplication using BLAS
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, *this, *m, (Real) 1.0, (Real) 0.0, result); 

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const SMatrix6& m) 
{
  const unsigned SZ = 6; 
  for (unsigned i=0; i< SZ; i++)
  {
    for (unsigned j=0; j< SZ-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,SZ-1) << std::endl;
  }
   
  return out;
}

