/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/cblas.h>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <Moby/Constants.h>
#include <Moby/Quat.h>
#include <Moby/AAngle.h>
#include <Moby/Matrix4.h>

using namespace Moby;

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
Matrix4::Matrix4()
{
  set_identity();
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion and translation vector
Matrix4::Matrix4(const Quat* q, const Vector3* v)
{
  Quat qn = Quat::normalize(*q);
  set(&qn, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion (for rotation) and zero translation
Matrix4::Matrix4(const Quat* q)
{
  Quat qn = Quat::normalize(*q);
  set(&qn, &ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and translation vector
Matrix4::Matrix4(const Matrix3* r, const Vector3* v) 
{
  set(r, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and zero translation
Matrix4::Matrix4(const Matrix3* r)
{
  set(r, &ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation and a translation vector
Matrix4::Matrix4(const AAngle* a, const Vector3* v)
{
  set(a, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation (for rotation) and zero translation
Matrix4::Matrix4(const AAngle* a)
{
  set(a, &ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix using identity orientation and a translation vector
Matrix4::Matrix4(const Vector3* v)
{
  const Matrix3 IDENTITY_3x3 = Matrix3::identity();
  set_translation(*v);
  set_rotation(&IDENTITY_3x3);
}

/// Constructs a matrix from an array
/**
 * \param array an array of 12 Real values in row-major format (all entries
 *        should be given except the last row
 */
Matrix4::Matrix4(const Real* array)
{
  const unsigned SZ = 3;
  CBLAS::copy(SZ, array, 1, _data, 1);
  CBLAS::copy(SZ, array+SZ, 1, _data+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ, 1, _data+SZ+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ+SZ, 1, _data+SZ+SZ+SZ, 1);
}

/// Constructs a matrix from 12 values
/**
 * The resulting matrix will appear as follows: <br />
 * m00 m01 m02 m03<br />
 * m10 m11 m12 m13<br />
 * m20 m21 m22 m23<br />
 * 0   0   0   1<br />
 */
Matrix4::Matrix4(Real m00, Real m01, Real m02, Real m03, Real m10, Real m11, Real m12, Real m13, Real m20, Real m21, Real m22, Real m23)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // set rotation component
  operator()(X,X) = m00;  operator()(X,Y) = m01;  operator()(X,Z) = m02;
  operator()(Y,X) = m10;  operator()(Y,Y) = m11;  operator()(Y,Z) = m12;
  operator()(Z,X) = m20;  operator()(Z,Y) = m21;  operator()(Z,Z) = m22;

  // set translation component
  operator()(X,W) = m03;  operator()(Y,W) = m13;  operator()(Z,W) = m23;
}

/// Interpolates between two 4x4 transforms using spherical linear interpolation
/**
 * \param T1 the matrix to use when t=0
 * \param T2 the matrix to use when t=1
 * \param t a real value in the interval [0,1]
 * \return the interpolated transform
 */
Matrix4 Matrix4::interpolate(const Matrix4& T1, const Matrix4& T2, Real t)
{
  // interpolate the positions
  Vector3 x = T1.get_translation()*(1-t) + T2.get_translation()*t;

  // interpolate the rotations
  Quat q1(&T1);
  Quat q2(&T2);
  Quat q = Quat::slerp(q1, q2, t);

  // return the new matrix
  return Matrix4(&q, &x);
}

/// Gets rotation 3x3 submatrix from a 4x4 homogeneous transform 
/**
 * \note identical to get_rotation()
 */
Matrix3 Matrix4::get_rotation() const 
{
  Matrix3 m;
  get_rotation(&m);
  return m;
}

/// Gets rotation 3x3 submatrix from a 4x4 homogeneous transform
void Matrix4::get_rotation(Matrix3* m) const 
{
  // copy each of three columns using BLAS
  const unsigned SZ = 3;
  for (unsigned i=0, inc = 0; i< 3; i++, inc +=3)
    CBLAS::copy(SZ, _data+inc, 1, m->begin()+inc, 1); 
}

/// Gets the translation component of a 4x4 homogeneous transform 
Vector3 Matrix4::get_translation() const
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  return Vector3(_data[TX], _data[TY], _data[TZ]);
}

/// Gets the translation component of a 4x4 homogeneous transform
/**
 * \param array a 3-dimensional (or larger) array
 */
void Matrix4::get_translation(Real* array) const
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  array[0] = _data[TX]; 
  array[1] = _data[TY]; 
  array[2] = _data[TZ]; 
}

/// Gets the translation component of a 4x4 homogeneous transform
void Matrix4::get_translation(Vector3& v) const
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  v[0] = _data[TX]; 
  v[1] = _data[TY]; 
  v[2] = _data[TZ]; 
}

/// Gets the translation component of a 4x4 homogeneous transform
/**
 * The values for x, y, and z are overwritten with the output.
 */
void Matrix4::get_translation(Real& x, Real& y, Real& z) const
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  x = _data[TX]; 
  y = _data[TY]; 
  z = _data[TZ]; 
}

/// Sets the matrix using an array of Real values
/**
 * \param array an array of 12 real values in row-major format (fourth row is 
 *        not modified)
 */
void Matrix4::set(const Real* array)
{
  const unsigned SZ = 3;
  CBLAS::copy(SZ, array, 1, _data, 1);
  CBLAS::copy(SZ, array+SZ, 1, _data+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ, 1, _data+SZ+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ+SZ, 1, _data+SZ+SZ+SZ, 1);
}

// Sets the matrix to be a 4x4 homogeneous transform from an axis-angle representation and a translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Matrix4::set(const AAngle* a, const Vector3*  v)
{
  set_rotation(a);
  set_translation(*v);
}

/// Sets the matrix to be a 4x4 homogeneous transform from a rotation matrix and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Matrix4::set(const Matrix3* m, const Vector3*  v)
{
  set_rotation(m);
  set_translation(*v);
}

/// Sets the matrix to be a 4x4 homogeneous transform from a unit quaternion and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Matrix4::set(const Quat* q, const Vector3*  v)
{
  set_rotation(q);
  set_translation(*v);
}

/// Sets the rotation component of a 4x4 homogeneous transform from the given axis-angle representation
void Matrix4::set_rotation(const AAngle* a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Real x = a->x;
  Real y = a->y;
  Real z = a->z;
  Real ca = std::cos(a->angle);
  Real sa = std::sin(a->angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const Real SOMEWHAT_NEAR_ZERO = 1e-2;
  Real va = (std::fabs(a->angle) > SOMEWHAT_NEAR_ZERO) ? 1 - ca : (sa*sa)/(1+ca);

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

/// Sets the rotation component of a 4x4 homogeneous transform from the given quaternion 
void Matrix4::set_rotation(const Quat* q)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Quat qn = Quat::normalize(*q);

  operator()(X,X) = 2*qn.x*qn.x + 2*qn.w*qn.w - 1;
  operator()(X,Y) = 2*qn.x*qn.y - 2*qn.z*qn.w;
  operator()(X,Z) = 2*qn.x*qn.z + 2*qn.y*qn.w;
  operator()(Y,X) = 2*qn.x*qn.y + 2*qn.z*qn.w;
  operator()(Y,Y) = 2*qn.y*qn.y + 2*qn.w*qn.w - 1;
  operator()(Y,Z) = 2*qn.y*qn.z - 2*qn.x*qn.w;
  operator()(Z,X) = 2*qn.x*qn.z - 2*qn.y*qn.w;
  operator()(Z,Y) = 2*qn.y*qn.z + 2*qn.x*qn.w;
  operator()(Z,Z) = 2*qn.z*qn.z + 2*qn.w*qn.w - 1;
}

/// Sets the rotation component of a 4x4 homogeneous transform from the given rotation matrix
void Matrix4::set_rotation(const Matrix3* m)
{
  const unsigned SZ = 9;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m->begin()[i];
}

/// Sets the translation component of a 4x4 homogeneous transform
void Matrix4::set_translation(Real x, Real y, Real z)
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  _data[TX] = x;
  _data[TY] = y;
  _data[TZ] = z;
}

/// Sets the translation component of a 4x4 homogeneous transform from the given 3-D array
void Matrix4::set_translation(const Real* array)
{
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;
  _data[TX] = array[X];
  _data[TY] = array[Y];
  _data[TZ] = array[Z];
}

/// Sets the translation componet of a 4x4 homogeneous transform from the given 3-D vector
void Matrix4::set_translation(const Vector3& v)
{
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;
  _data[TX] = v[X];
  _data[TY] = v[Y];
  _data[TZ] = v[Z];
}

/// Sets this matrix to identity
void Matrix4::set_identity()
{
  const unsigned XX = 0, YY = 4, ZZ = 8, N = 12;
  for (unsigned i=0; i< N; i++)
    _data[i] = (Real) 0.0;
  _data[XX] = _data[YY] = _data[ZZ] = 1.0;
}

/// Multiplies the 3x3 rotational component of this matrix by a vector (i.e., translation is not added)
Vector3 Matrix4::mult_vector(const Vector3& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3;

  // this is just regular 3x3 matrix multiplication
  Vector3 result;
  CBLAS::gemv(CblasColMajor, CblasNoTrans, ROWS, ROWS, (Real) 1.0, begin(), ROWS, v.begin(), 1, (Real) 0.0, result.begin(), 1);

  return result;
}

/// Multiplies the transpose of the 3x3 rotational component of this matrix by a vector (i.e., translation is not added)
Vector3 Matrix4::transpose_mult_vector(const Vector3& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3, X = 0, Y = 1, Z = 2;

  // setup starting rows indices
  const unsigned RX = 0, RY = 3, RZ = 6;

  // this is 3x3 matrix multiplication
  Vector3 result;
  result[X] = CBLAS::dot(ROWS, _data+RX, 1, v.data(), 1);
  result[Y] = CBLAS::dot(ROWS, _data+RY, 1, v.data(), 1);
  result[Z] = CBLAS::dot(ROWS, _data+RZ, 1, v.data(), 1); 

  return result;
}

/// Multiplies a 4x4 matrix by a point 
Vector3 Matrix4::mult_point(const Vector3& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3;
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;

  // matrix/vec multiplication is just regular 3x3 matrix multiplication
  Vector3 result;
  CBLAS::gemv(CblasColMajor, CblasNoTrans, ROWS, ROWS, (Real) 1.0, begin(), ROWS, v.begin(), 1, (Real) 0.0, result.begin(), 1);

  // add in the translational components
  result[X] += _data[TX];
  result[Y] += _data[TY];
  result[Z] += _data[TZ];

  return result;
}

/// Multiplies the inverse of the 4x4 matrix by a point
/**
 * \note if multiple inverse_mult_point() operations are performed, it may be
 *       faster to use the inverse_transform() function, followed by 
 *       mult_point().
 */ 
Vector3 Matrix4::inverse_mult_point(const Vector3& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3, X = 0, Y = 1, Z = 2;

  // setup starting rows indices
  const unsigned RX = 0, RY = 3, RZ = 6, RW = 9;

  // determine R' * v 
  Vector3 result1;
  result1[X] = CBLAS::dot(ROWS, _data+RX, 1, v.data(), 1);
  result1[Y] = CBLAS::dot(ROWS, _data+RY, 1, v.data(), 1);
  result1[Z] = CBLAS::dot(ROWS, _data+RZ, 1, v.data(), 1); 

  // determine R' * -t
  Vector3 result2;
  result2[X] = -CBLAS::dot(ROWS, _data+RX, 1, _data + RW, 1);
  result2[Y] = -CBLAS::dot(ROWS, _data+RY, 1, _data + RW, 1);
  result2[Z] = -CBLAS::dot(ROWS, _data+RZ, 1, _data + RW, 1); 

  // add the two together
  return result1 + result2;
}

/// Determines whether the transformation matrix represented by this matrix is equal to another (to within numerical tolerance epsilon)
bool Matrix4::epsilon_equals(const Matrix4& m, Real epsilon) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // first check position
  Vector3 p1 = get_translation();
  Vector3 p2 = m.get_translation();
  if (std::fabs(p1[X] - p2[X]) > epsilon ||
      std::fabs(p1[Y] - p2[Y]) > epsilon ||
      std::fabs(p1[Z] - p2[Z]))
    return false;
      
  // check orientation  
  Quat q1(this);
  Quat q2(&m);
  return Quat::epsilon_equals(q1, q2, epsilon);
}

/// Determines whether the transformation matrix represented by two matrices are equal (to within numerical tolerance epsilon)
bool Matrix4::epsilon_equals(const Matrix4& m1, const Matrix4& m2, Real epsilon)
{
  // first check position
  Vector3 p1 = m1.get_translation();
  Vector3 p2 = m2.get_translation();
  if (Vector3::norm(p1 - p2) > epsilon)
    return false;

  // check orientation
  Quat q1(&m1);
  Quat q2(&m2);
  return Quat::epsilon_equals(q1, q2, epsilon);
}

/// Determines whether the specified transform is valid
bool Matrix4::valid_transform(const Matrix4& T)
{
  // check that the rotation component is correct
  Matrix3 R;
  T.get_rotation(&R);
  return Matrix3::valid_rotation(R);
}

/// Special method for inverting a 4x4 transformation matrix in place
void Matrix4::invert_transform()
{
  // verify that it is a valid transform
  assert(valid_transform(*this));

  // get the new rotation matrix
  Matrix3 RT;
  get_rotation(&RT);
  RT.transpose();

  // determine the new translation
  Vector3 new_trans = RT * -get_translation();

  // set the new rotation and translation
  set_rotation(&RT);
  set_translation(new_trans);
}

/// Special method for inverseing a 4x4 transformation matrix
Matrix4 Matrix4::inverse_transform(const Matrix4& T)
{
  Matrix4 m(T);
  m.invert_transform();
  return m;
}

/// Copies a matrix
Matrix4& Matrix4::operator=(const Matrix4& m)
{
  const unsigned SZ = 12;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m.data()[i];
  return *this;
}

/// Multiplies this matrix by another
Matrix4 Matrix4::operator*(const Matrix4& m) const
{
  // determine the rotation part of the matrix
  const unsigned SZ = 3, XLAT_BEGIN = 9, TX = 9, TY = 10, TZ = 11;
  Matrix4 result;
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, SZ, SZ, SZ, (Real) 1.0, _data, SZ, m._data, SZ, (Real) 0.0, result._data, SZ);

  // determine the translation part of the matrix
  CBLAS::gemv(CblasColMajor, CblasNoTrans, SZ, SZ, (Real) 1.0, _data, SZ, m._data + XLAT_BEGIN, 1, (Real) 0.0, result._data + XLAT_BEGIN, 1);
  result._data[TX] += _data[TX];
  result._data[TY] += _data[TY];
  result._data[TZ] += _data[TZ];
  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const Matrix4& m)
{
  const unsigned ROWS = 3, COLUMNS = 4; 
  for (unsigned i=0; i< ROWS; i++)
  {
    for (unsigned j=0; j< COLUMNS-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,COLUMNS-1) << std::endl;
  }
  out << "0 0 0 1" << std::endl;
   
  return out;
}

