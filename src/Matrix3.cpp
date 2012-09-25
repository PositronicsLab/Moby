/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/cblas.h>
#include <iomanip>
#include <cstring>
#include <iostream>
#include <cmath>
#include <Moby/AAngle.h>
#include <Moby/Quat.h>
#include <Moby/MissizeException.h>
#include <Moby/MatrixN.h>
#include <Moby/Matrix3.h>

using namespace Moby;

/// Sets the orientation of this matrix to that represented by the given quaternion 
Matrix3::Matrix3(const Quat* q)
{
  Quat qn = Quat::normalize(*q);
  set(&qn);
}

/// Sets the orientation fo this matrix to that represented by the given axis-angle parameter
Matrix3::Matrix3(const AAngle* a)
{
  set(a);
}

/// Constructs a matrix from an array
/**
 * \param array an array of 9 Real values in row-major format
 */
Matrix3::Matrix3(const Real* array)
{
  for (unsigned i=0; i< 9; i++)
    _data[i] = array[i];
}

/// Constructs a matrix from 9 values
/**
 * The resulting matrix will appear as follows: <br />
 * m00 m01 m02 <br />
 * m10 m11 m12 <br />
 * m20 m21 m22
 */
Matrix3::Matrix3(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22)
{
  const unsigned X = 0, Y = 1, Z = 2;
  operator()(X,X) = m00;
  operator()(X,Y) = m01;
  operator()(X,Z) = m02;
  operator()(Y,X) = m10;
  operator()(Y,Y) = m11;
  operator()(Y,Z) = m12;
  operator()(Z,X) = m20;
  operator()(Z,Y) = m21;
  operator()(Z,Z) = m22;
}

/// Computes the l-infinity norm of this matrix
Real Matrix3::norm_inf() const
{
  const unsigned NELMS = 9;

  Real nrm = (Real) 0.0;
  for (unsigned i=0; i< NELMS; i++)
    nrm = std::max(nrm, std::fabs(_data[i]));

  return nrm;
}

/// Constructs the vector from a skew-symmetric matrix
Vector3 Matrix3::inverse_skew_symmetric(const Matrix3& R)
{
  const unsigned X = 0, Y = 1, Z = 2;

  Real vx = (R(Z,Y) - R(Y,Z)) * (Real) 0.5;
  Real vy = (R(X,Z) - R(Z,X)) * (Real) 0.5;
  Real vz = (R(Y,X) - R(X,Y)) * (Real) 0.5;

  return Vector3(vx, vy, vz);
}

/// Constructs a skew-symmetric matrix from the given values
/**
 * The skew symmetric matrix generated will be:
 * | 0   -z   y |
 * | z    0  -x |
 * | -y   x   0 |
 */
Matrix3 Matrix3::skew_symmetric(Real x, Real y, Real z)
{
  return Matrix3(0,-z,y,z,0,-x,-y,x,0);
}

/// Constructs a skew-symmetric matrix from the given values
/**
 * The skew symmetric matrix generated will be:
 * |   0     -v[2]   v[1] |
 * |  v[2]     0    -v[0] |
 * | -v[1]    v[0]    0   |
 */
Matrix3 Matrix3::skew_symmetric(const Vector3& v)
{
  return Matrix3(0,-v[2],v[1],v[2],0,-v[0],-v[1],v[0],0);
}

/// Inverts this matrix 
void Matrix3::inverse()
{
  *this = inverse(*this);
}

/// Calculates the determinant for a 3x3 matrix
Real Matrix3::det() const
{
  const unsigned X = 0, Y = 1, Z = 2;
   Real v1 = operator()(X,X) * (operator()(Y,Y) * operator()(Z,Z) - operator()(Y,Z) * operator()(Z,Y)); 
   Real v2 = operator()(Y,X) * (operator()(X,Z) * operator()(Z,Y) - operator()(X,Y) * operator()(Z,Z)); 
   Real v3 = operator()(Z,X) * (operator()(X,Y) * operator()(Y,Z) - operator()(X,Z) * operator()(Y,Y)); 
  return v1 + v2 + v3;
}

/// Determines whether this is an orthonormal matrix
bool Matrix3::is_orthonormal() const
{
  Real determinant = det();
  return (std::fabs(determinant - (Real) 1.0) < NEAR_ZERO || std::fabs(determinant + (Real) 1.0) < NEAR_ZERO);
}

/// Determines the inverse of the given matrix
Matrix3 Matrix3::inverse(const Matrix3& m)
{
  // compute the determinant
  Real determ = 1.0/m.det();
  Real m00 = determ * (m(1,1)*m(2,2) - m(2,1)*m(1,2));
  Real m01 = determ * (m(2,1)*m(0,2) - m(0,1)*m(2,2));
  Real m02 = determ * (m(0,1)*m(1,2) - m(1,1)*m(0,2));
  Real m10 = determ * (m(1,2)*m(2,0) - m(1,0)*m(2,2));
  Real m11 = determ * (m(0,0)*m(2,2) - m(2,0)*m(0,2));
  Real m12 = determ * (m(1,0)*m(0,2) - m(0,0)*m(1,2));
  Real m20 = determ * (m(1,0)*m(2,1) - m(2,0)*m(1,1));
  Real m21 = determ * (m(2,0)*m(0,1) - m(0,0)*m(2,1));
  Real m22 = determ * (m(0,0)*m(1,1) - m(0,1)*m(1,0));
  return Matrix3(m00, m01, m02, m10, m11, m12, m20, m21, m22);
}

/// Sets this matrix to the rotation matrix of the specified angle around the X axis
void Matrix3::set_rot_X(Real angle)
{
  *this = rot_X(angle);
}

/// Returns the rotation matrix of the specified angle around the X axis
Matrix3 Matrix3::rot_X(Real angle)
{
  Real sina = std::sin(angle);
  Real cosa = std::cos(angle);
  Matrix3 r;
  r(0,0) = 1;
  r(0,1) = 0;
  r(0,2) = 0;
  r(1,0) = 0;
  r(1,1) = cosa;
  r(1,2) = -sina;
  r(2,0) = 0;
  r(2,1) = sina;
  r(2,2) = cosa;
  return r;
}

/// Sets this matrix to the rotation matrix of the specified angle around the Y axis
void Matrix3::set_rot_Y(Real angle)
{
  *this = rot_Y(angle);
}

/// Returns the rotation matrix of the specified angle around the Y axis
Matrix3 Matrix3::rot_Y(Real angle)
{
  Real sina = std::sin(angle);
  Real cosa = std::cos(angle);
  Matrix3 r;
  r(0,0) = cosa;
  r(0,1) = 0;
  r(0,2) = sina;
  r(1,0) = 0;
  r(1,1) = 1;
  r(1,2) = 0;
  r(2,0) = -sina;
  r(2,1) = 0;
  r(2,2) = cosa;
  return r;
}

/// Sets this matrix to the rotation matrix of the specified angle around the Z axis
void Matrix3::set_rot_Z(Real angle)
{
  *this = rot_Z(angle);
}

/// Returns the rotation matrix of the specified angle around the Z axis
Matrix3 Matrix3::rot_Z(Real angle)
{
  Real sina = std::sin(angle);
  Real cosa = std::cos(angle);
  Matrix3 r;
  r(0,0) = cosa;
  r(0,1) = -sina;
  r(0,2) = 0;
  r(1,0) = sina;
  r(1,1) = cosa;
  r(1,2) = 0;
  r(2,0) = 0;
  r(2,1) = 0;
  r(2,2) = 1;
  return r;
}

/// Sets this rotation matrix from the given axis-angle representation
/**
 * Note: this method could be faster...
 */
void Matrix3::set(const AAngle* a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Real x = a->x;
  Real y = a->y;
  Real z = a->z;
  Real ca = cos(a->angle);
  Real sa = sin(a->angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const Real SOMEWHAT_NEAR_ZERO = 1e-2;
  Real va = (std::fabs(a->angle) > SOMEWHAT_NEAR_ZERO) ? 1 - ca : (sa*sa)/(1+ca);

  // setup the matrix
  operator()(X,X) = x*x*va + ca;
  operator()(X,Y) = x*y*va - z*sa;
  operator()(X,Z) = x*z*va + y*sa;
  operator()(Y,X) = x*y*va + z*sa;
  operator()(Y,Y) = y*y*va + ca;
  operator()(Y,Z) = y*z*va - x*sa;
  operator()(Z,X) = x*z*va - y*sa;
  operator()(Z,Y) = y*z*va + x*sa;
  operator()(Z,Z) = z*z*va + ca;
}

/// Sets this rotation matrix from the given quaternion
void Matrix3::set(const Quat* q)
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // verify that the quaternion is normalized
  assert(std::fabs(q->magnitude()) - 1 < NEAR_ZERO);

  // setup repeated products
  const Real xx = q->x*q->x;
  const Real xy = q->x*q->y;
  const Real xz = q->x*q->z;
  const Real xw = q->x*q->w;
  const Real yy = q->y*q->y;
  const Real yz = q->y*q->z;
  const Real yw = q->y*q->w;
  const Real zz = q->z*q->z;
  const Real zw = q->z*q->w; 
  const Real ww = q->w*q->w;

  operator()(X,X) = 2*(xx + ww) - 1;
  operator()(X,Y) = 2*(xy - zw);
  operator()(X,Z) = 2*(xz + yw);
  operator()(Y,X) = 2*(xy + zw);
  operator()(Y,Y) = 2*(yy + ww) - 1;
  operator()(Y,Z) = 2*(yz - xw);
  operator()(Z,X) = 2*(xz - yw);
  operator()(Z,Y) = 2*(yz + xw);
  operator()(Z,Z) = 2*(zz + ww) - 1;
}

/// Determines whether the rotation matrix represented by this matrix is equal to another (to within numerical tolerance epsilon)
bool Matrix3::epsilon_equals(const Matrix3& m, Real epsilon) const
{
  Quat q1(this);
  Quat q2(&m);
  return Quat::epsilon_equals(q1, q2, epsilon);
}

/// Determines whether the rotation matrix represented by two matrices are equal (to within numerical tolerance epsilon)
bool Matrix3::epsilon_equals(const Matrix3& m1, const Matrix3& m2, Real epsilon)
{
  Quat q1(&m1);
  Quat q2(&m2);
  return Quat::epsilon_equals(q1, q2, epsilon);
}

/// Checks whether this matrix represents a valid rotation and scale for a right-handed coordinate system
bool Matrix3::valid_rotation_scale(const Matrix3& R)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real TOLERANCE = std::sqrt(NEAR_ZERO);

  // create vectors for each column
  Vector3 cx, cy, cz;
  R.get_column(X, cx.begin());
  R.get_column(Y, cy.begin());
  R.get_column(Z, cz.begin());

  // get the norms for each column
  Real cx_norm = cx.norm();
  Real cy_norm = cy.norm();
  Real cz_norm = cz.norm();

  // check cross product for right-handedness
  Vector3 cp = Vector3::cross(cx/cx_norm, cy/cy_norm);
  Real cp_norm = cp.norm();

  return (std::fabs(Vector3::dot(cp/cp_norm, cz/cz_norm) - (Real) 1.0) < TOLERANCE);
}

/// Checks whether this matrix represents a valid rotation for a right-handed coordinate system
bool Matrix3::valid_rotation(const Matrix3& R)
{
  Real TOLERANCE = std::sqrt(NEAR_ZERO);

  // check the determinant of the matrix 
  return (std::fabs(R.det() - 1.0) < TOLERANCE);
}

/// Graham-Schmidt orthonormalization
bool Matrix3::orthonormalize(Vector3& a, Vector3& b, Vector3& c)
{
  // Gram-Schmidt orthonormalization produces vectors u0, u1, and u2 as
  // follows:
  //
  //   u0 = v0/|v0|
  //   u1 = (v1 - Dot(u0,v1)*u0)/|v1 - Dot(u0,v1)*u0|
  //   u2 = (v2 - Dot(u0,v2)*u0 - Dot(u1,v2)*u1) /
  //          |v2 - Dot(u0,v2)*u0 - Dot(u1,v2)*u1|

  // Compute u0.
  Real flen = a.norm();
  if (flen == 0.0)
    return false;
  a /= flen;

  // Compute u1.
  Real fd0 = Vector3::dot(a, b);
  b -= fd0*a;
  flen = b.norm();
  if (flen == 0.0)
    return false;
  b /= flen;

  // Compute u2.
  Real fd1 = Vector3::dot(b, c);
  fd0 = Vector3::dot(a,c);
  c -= fd0*a + fd1*b;
  flen = c.norm();
  if (flen == 0.0)
    return false;

  c /= flen;
  return true;
}

/// Makes the matrix orthonormal using Gram-Schmidt orthogonalization
bool Matrix3::orthonormalize()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // The algorithm uses Gram-Schmidt orthogonalization applied to the
  // columns of M.
  Vector3 a(operator()(X,X),operator()(X,Y),operator()(X,Z));
  Vector3 b(operator()(Y,X),operator()(Y,Y),operator()(Y,Z));
  Vector3 c(operator()(Z,X),operator()(Z,Y),operator()(Z,Z));
  if (orthonormalize(a,b,c))
  {
    operator()(X,X) = a[X];
    operator()(X,Y) = a[Y];
    operator()(X,Z) = a[Z];
    operator()(Y,X) = b[X];
    operator()(Y,Y) = b[Y];
    operator()(Y,Z) = b[Z];
    operator()(Z,X) = c[X];
    operator()(Z,Y) = c[Y];
    operator()(Z,Z) = c[Z];
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
void Matrix3::set_row(unsigned i, const Vector3& v)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  _data[st] = v[0];
  st += SZ;
  _data[st] = v[1]; 
  st += SZ;
  _data[st] = v[2]; 
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
Vector3 Matrix3::get_row(unsigned i) const
{
  Vector3 v;
  get_row(i, v);
  return v;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param target an array of size columns() or greater
 */
void Matrix3::get_row(unsigned i, Real* target) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  target[0] = _data[st];
  st += SZ;
  target[1] = _data[st];
  st += SZ;
  target[2] = _data[st];
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
void Matrix3::get_row(unsigned i, Vector3& result) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  result[0] = _data[st];
  st += SZ;
  result[1] = _data[st];
  st += SZ;
  result[2] = _data[st];
}

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param v a 3D vector 
 */
void Matrix3::set_column(unsigned i, const Vector3& v)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  _data[st++] = v[0];
  _data[st++] = v[1];
  _data[st] = v[2];
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 */
Vector3 Matrix3::get_column(unsigned i) const
{
  Vector3 v;
  get_column(i, v);
  return v;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param target an array of size rows() or greater
 */
void Matrix3::get_column(unsigned i, Real* target) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  target[0] = _data[st++];
  target[1] = _data[st++];
  target[2] = _data[st];
}

/// Gets the specified column of the matrix
/*
 * \param i the 0-index of the column
 */
void Matrix3::get_column(unsigned i, Vector3& result) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  result[0] = _data[st++];
  result[1] = _data[st++];
  result[2] = _data[st];
}

/// Sets this matrix to its transpose
void Matrix3::transpose()
{
  const unsigned YX = 1, ZX = 2, XY = 3, ZY = 5, XZ = 6, YZ = 7;
  std::swap(_data[XY], _data[YX]);
  std::swap(_data[XZ], _data[ZX]);
  std::swap(_data[YZ], _data[ZY]);
}

/// Determines the transpose of this matrix
Matrix3 Matrix3::transpose(const Matrix3& m)
{
  Matrix3 n = m;
  n.transpose();
  return n;
}

/// Copies a matrix to this one
Matrix3& Matrix3::operator=(const Matrix3& m)
{
  const unsigned SZ = 9;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m.begin()[i];
  return *this;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
Vector3 Matrix3::mult(const Vector3& v) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Vector3 result;
  CBLAS::gemv(CblasNoTrans, ROWS, COLUMNS, *this, ROWS, v, 1, (Real) 1.0, (Real) 0.0, result, 1);
  return result;
}

/// Multiplies the transpose of this matrix by a vector and returns the result in a new vector
Vector3 Matrix3::transpose_mult(const Vector3& v) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Vector3 result;
  CBLAS::gemv(CblasTrans, COLUMNS, ROWS, *this, ROWS, v, 1, (Real) 1.0, (Real) 0.0, result, 1);
  return result;
}

/// Multiplies the transpose of this matrix by a matrix and returns the result in a given matrix
MatrixN* Matrix3::transpose_mult(const MatrixN* m, MatrixN* result) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  // verify that m is the proper size
  if (m->rows() != size())
    throw MissizeException();

  // resize the result matrix
  result->resize(size(), m->columns());

  // do the multiplication
  CBLAS::gemm(CblasTrans, CblasNoTrans, COLUMNS, ROWS, m->rows(), *this, ROWS, *m, m->rows(), (Real) 1.0, (Real) 0.0, *result, COLUMNS);

  return result;
}

/// Multiplies the transpose of this matrix by a matrix and returns the result in a new matrix 
Matrix3 Matrix3::transpose_mult(const Matrix3& m) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Matrix3 result;
  CBLAS::gemm(CblasTrans, CblasNoTrans, COLUMNS, ROWS, m.columns(), *this, ROWS, m, m.rows(), (Real) 1.0, (Real) 0.0, result, COLUMNS);
  return result;
}

/// Multiplies the transpose of this matrix by the transpose of a matrix and returns the result in a given matrix
MatrixN* Matrix3::transpose_mult_transpose(const MatrixN* m, MatrixN* result) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  // verify that m is the proper size
  if (m->columns() != size())
    throw MissizeException();

  // resize the result matrix
  result->resize(size(), m->rows());

  // do the multiplication
  CBLAS::gemm(CblasTrans, CblasTrans, COLUMNS, ROWS, m->rows(), *this, ROWS, *m, m->rows(), (Real) 1.0, (Real) 0.0, *result, COLUMNS);

  return result;
}

/// Multiplies the transpose of this matrix by the transpose of a matrix and returns the result in a new matrix 
Matrix3 Matrix3::transpose_mult_transpose(const Matrix3& m) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Matrix3 result;
  CBLAS::gemm(CblasTrans, CblasTrans, COLUMNS, ROWS, m.rows(), *this, ROWS, m, m.rows(), (Real) 1.0, (Real) 0.0, result, COLUMNS);
  return result;
}

/// Multiplies this matrix by the transpose of a matrix and returns the result in a given matrix
MatrixN* Matrix3::mult_transpose(const MatrixN* m, MatrixN* result) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  // verify that m is the proper size
  if (m->columns() != size())
    throw MissizeException();

  // resize the result matrix
  result->resize(size(), m->rows());

  // do the multiplication
  CBLAS::gemm(CblasNoTrans, CblasTrans, ROWS, COLUMNS, m->rows(), *this, ROWS, *m, m->rows(), (Real) 1.0, (Real) 0.0, *result, ROWS);

  return result;
}

/// Multiplies this matrix by the transpose of a matrix and returns the result in a new matrix 
Matrix3 Matrix3::mult_transpose(const Matrix3& m) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Matrix3 result;
  CBLAS::gemm(CblasNoTrans, CblasTrans, ROWS, COLUMNS, m.rows(), *this, ROWS, m, m.rows(), (Real) 1.0, (Real) 0.0, result, ROWS);
  return result;
}

/// Multiplies this matrix by a matrix and returns the result in a given matrix
MatrixN* Matrix3::mult(const MatrixN* m, MatrixN* result) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  // verify that m is the proper size
  if (m->rows() != size())
    throw MissizeException();

  // resize the result matrix
  result->resize(size(), m->columns());

  // do the multiplication
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, ROWS, COLUMNS, m->rows(), *this, ROWS, *m, m->rows(), (Real) 1.0, (Real) 0.0, *result, ROWS);

  return result;
}

/// Multiplies this matrix by another 3x3 matrix
Matrix3 Matrix3::mult(const Matrix3& m) const
{
  const unsigned ROWS = 3, COLUMNS = 3;
  Matrix3 result;
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, ROWS, COLUMNS, m.rows(), *this, ROWS, m, m.rows(), (Real) 1.0, (Real) 0.0, result, ROWS);
  return result;
}

/// Multiplies this matrix by a scalar in place
Matrix3& Matrix3::operator*=(Real scalar)
{
  const unsigned SZ = 3*3;
  CBLAS::scal(SZ, scalar, begin(), 1);
  return *this;
}

/// Returns the negation of this matrix
Matrix3 Matrix3::operator-() const
{
  Matrix3 m;
  std::transform(begin(), end(), m.begin(), std::negate<Real>());
  return m;
}

/// Adds m to this in place
Matrix3& Matrix3::operator+=(const Matrix3& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::plus<Real>());
  return *this;
}

/// Subtracts m from this in place
Matrix3& Matrix3::operator-=(const Matrix3& m)
{
  std::transform(begin(), end(), m.begin(), begin(), std::minus<Real>());
  return *this;
}

/// Sets this matrix to identity
void Matrix3::set_identity()
{
  const unsigned XX = 0, YY = 4, ZZ = 8;
  set_zero();
  _data[XX] = _data[YY] = _data[ZZ] = 1.0;
}

/// Sets this matrix to zero
void Matrix3::set_zero()
{
  const unsigned N = 3*3;
  for (unsigned i=0; i< N; i++)
    _data[i] = (Real) 0.0;
}

/// Checks whether this matrix is symmetric
bool Matrix3::is_symmetric(Real tolerance) const
{
  // make sure tolerance is non-negative
  tolerance = std::fabs(tolerance);

  // check symmetric
  for (unsigned i=0; i< size()-1; i++)
    for (unsigned j=i+1; j< size(); j++)
      if (!rel_equal(operator()(i,j), operator()(j,i), tolerance))
        return false;

  return true;
}

/// Calculates the differential between two rotations
Vector3 Matrix3::calc_differential(const Matrix3& R1, const Matrix3& R2)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Vector3 R1x = R1.get_column(X);
  Vector3 R1y = R1.get_column(Y);
  Vector3 R1z = R1.get_column(Z);
  Vector3 R2x = R2.get_column(X);
  Vector3 R2y = R2.get_column(Y);
  Vector3 R2z = R2.get_column(Z);
  return (Vector3::cross(R1x,R2x) + Vector3::cross(R1y,R2y) + 
          Vector3::cross(R1z,R2z))*0.5;
}

/// Determines whether two numbers are relatively equal
bool Matrix3::rel_equal(Real x, Real y, Real tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), (Real) 1.0)));
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const Matrix3& m)
{
  const unsigned SZ = 3; 
  for (unsigned i=0; i< SZ; i++)
  {
    for (unsigned j=0; j< SZ-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,SZ-1) << std::endl;
  }
   
  return out;
}

