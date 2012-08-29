/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/AAngle.h>
#include <Moby/InvalidIndexException.h>
#include <Moby/MissizeException.h>
#include <Moby/Quat.h>

using namespace Moby;

/// Default constructor
/**
 * Quaternion constructed as [0,0,0] 1
 */
Quat::Quat()
{
  x = y = z = 0;
  w = 1;
}

/// Constructs a quaternion from an axis-angle object
Quat::Quat(const AAngle* a)
{
  set(a);
}

/// Constructs a quaternion from four Real values
/**
 * \param x the first component of the vector
 * \param y the second component of the vector
 * \param z the third component of the vector
 * \param w the scalar component of the quaternion
 */
Quat::Quat(Real x, Real y, Real z, Real w)
{
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

/// Constructs a quaternion from a vector
/**
 * \param v a four-dimensional vector (components 0..2 are the vector component, component 3 is the angular component)
 */
Quat::Quat(const VectorN& v)
{
  set(v);
}

/// Constructs a quaternion from a rotation matrix
Quat::Quat(const Matrix3* m)
{
  set(m);
}

/// Constructs a quaternion from a transformation matrix
Quat::Quat(const Matrix4* m)
{
  set(m);
}

/// Constructs a quaternion from a Real array of length 4
Quat::Quat(boost::shared_array<const Real>  q)
{
  x = q[0];
  y = q[1];
  z = q[2];
  w = q[3];
}

/// Copy constructor
Quat::Quat(const Quat& q)
{
  *this = q;
}

/// Constructs a zero quaternion
Quat Quat::zero()
{
  Quat q;
  q.x = q.y = q.z = q.w = 0;
  return q;
}

/// Gets the i'th component of the quaternion
/**
 * \param i for i=0, returns the first component of the vector (x), for
 *        i=1, returns the second component of the vector (y), for i=2,
 *        returns the third component of the vector (z), for i=3, returns
 *        the scalar (w).
 */
Real& Quat::operator[](unsigned i)
{
  switch (i)
  {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default:
      throw InvalidIndexException();
  }

  // make compiler happy
  return x;
}

/// Gets the i'th component of the quaternion
/**
 * \param i for i=0, returns the first component of the vector (x), for
 *        i=1, returns the second component of the vector (y), for i=2,
 *        returns the third component of the vector (z), for i=3, returns
 *        the scalar (w).
 */
Real Quat::operator[](unsigned i) const
{
  switch (i)
  {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default:
      throw InvalidIndexException();
  }

  // make compiler happy
  return x;
}


/// Copy operator
Quat& Quat::operator=(const Quat& q)
{
  x = q.x;
  y = q.y;
  z = q.z;
  w = q.w;
  return *this;
}

/// Negates the value of each of the x, y, and z coordinates in place
void Quat::conjugate()
{
  *this = conjugate(*this);
}

/// Negates the value of each of the x, y, and z coordinates of the given quaternion
/**
 * \param q the quaternion from which the conjugate is to be computed
 */
Quat Quat::conjugate(const Quat& q)
{
  return Quat(-q.x, -q.y, -q.z, q.w);
}

/// Computes the 3x4 matrix 'L' times a four dimensional array
Vector3 Quat::L_mult(const Quat& q) const
{
  Vector3 v;
  v[0] = -x*q[0] + w*q[1] + z*q[2] - y*q[3];
  v[1] = -y*q[0] - z*q[1] + w*q[2] + x*q[3];
  v[2] = -z*q[0] + y*q[1] - x*q[2] + w*q[3];
  return v;
}

/// Computes the matrix 'L' used for generalized coordinate calculations
MatrixN& Quat::determine_L(MatrixN& L) const
{
  L.resize(3,4);
  L(0,0) = -x;  L(0,1) = +w;  L(0,2) = +z;  L(0,3) = -y;
  L(1,0) = -y;  L(1,1) = -z;  L(1,2) = +w;  L(1,3) = +x;
  L(2,0) = -z;  L(2,1) = +y;  L(2,2) = -x;  L(2,3) = +w;
  return L;
}

/// Computes the matrix 'G' used for generalized coordinate calculations
MatrixN& Quat::determine_G(MatrixN& G) const
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  G.resize(3,4);
  G(X,X) = -x;  G(X,Y) = +w;  G(X,Z) = -z;  G(X,W) = +y;
  G(Y,X) = -y;  G(Y,Y) = +z;  G(Y,Z) = +w;  G(Y,W) = -x;
  G(Z,X) = -z;  G(Z,Y) = -y;  G(Z,Z) = +x;  G(Z,W) = +w;
  return G;
}

/// Multiplies the matrix 'G' by a quaternion vector
Vector3 Quat::G_mult(Real qw, Real qx, Real qy, Real qz) const
{
  Vector3 r;
  r[0] = -x*qw + w*qx - z*qy + y*qz;
  r[1] = -y*qw + z*qx + w*qy - x*qz;
  r[2] = -z*qw - y*qx + x*qy + w*qz;
  return r;
}

/// Multiplies the matrix 'G' (transpose) by a vector
Quat Quat::G_transpose_mult(const Vector3& v) const
{
  Quat output;
  output.w = -x*v[0] - y*v[1] - z*v[2];
  output.x = w*v[0] + z*v[1] - y*v[2];
  output.y = -z*v[0] + w*v[1] + x*v[2];
  output.z = y*v[0] - x*v[1] + w*v[2];
  return output;
}

/// Multiplies the matrix 'L' (transpose) by a vector
Quat Quat::L_transpose_mult(const Vector3& v) const
{
  Quat output;
  output.w = -x*v[0] - y*v[1] - z*v[2];
  output.x = +w*v[0] - z*v[1] + y*v[2];
  output.y = +z*v[0] + w*v[1] - x*v[2];
  output.z = -y*v[0] + x*v[1] + w*v[2];
  return output;
}

/// Computes the derivative of a quaternion given current orientation and angular velocity
/**
 * Note that the derivative is not generally a unit quaternion.
 * \param q the current orientation
 * \param w the angular velocity (in the global frame)
 * Uses the matrix: 
 *      |  -q.x  +q.w  -q.z  +q.y  |
 * G =  |  -q.y  +q.z  +q.w  -q.x  |
 *      |  -q.z  -q.y  +q.x  +q.w  |
 */
Quat Quat::deriv(const Quat& q, const Vector3& w)
{
  Quat qd;

  qd.w = .5 * (-q.x * w[0] - q.y * w[1] - q.z * w[2]); 
  qd.x = .5 * (+q.w * w[0] + q.z * w[1] - q.y * w[2]);
  qd.y = .5 * (-q.z * w[0] + q.w * w[1] + q.x * w[2]);
  qd.z = .5 * (+q.y * w[0] - q.x * w[1] + q.w * w[2]);

  return qd;
}

/// Computes the angular velocity of a body given the current quaternion orientation and the quaternion velocity
Vector3 Quat::to_omega(const Quat& q, const Quat& qd)
{
  Vector3 omega;
  omega[0] = 2 * (-q.x * qd.w + q.w * qd.x - q.z * qd.y + q.y * qd.z);
  omega[1] = 2 * (-q.y * qd.w + q.z * qd.x + q.w * qd.y - q.x * qd.z);
  omega[2] = 2 * (-q.z * qd.w - q.y * qd.x + q.x * qd.y + q.w * qd.z);
  return omega;
}

/// Calculates the second derivative of a quaternion
/**
 * \note alpha and omega are acceleration and velocity vectors in the global 
 *       frame
 */
Quat Quat::dderiv(const Quat& q, const Vector3& omega, const Vector3& alpha)
{
  Quat qdd = Quat::deriv(q, alpha) - (Real) 0.25 * omega.norm_sq() * q;
 
  return qdd;
}

/// Determines whether this quaternion is equal to another within the given tolerance
bool Quat::epsilon_equals(const Quat& q, Real epsilon)
{
  return epsilon_equals(*this, q, epsilon);
}

/// Determines whether to quaternions are equal to within the given tolerance
bool Quat::epsilon_equals(const Quat& q1, const Quat& q2, Real epsilon)
{
  if (std::fabs(q1.x - q2.x) < epsilon && std::fabs(q1.y - q2.y) < epsilon && std::fabs(q1.z - q2.z) < epsilon && std::fabs(q1.w - q2.w) < epsilon)
    return true;
  else
    return (std::fabs(-q1.x - q2.x) < epsilon && std::fabs(-q1.y - q2.y) < epsilon && std::fabs(-q1.z - q2.z) < epsilon && std::fabs(-q1.w - q2.w) < epsilon);
}

/// Performs spherical linear interpolation between this and q
/**
 * Sets the orientation represented by this quaternion to a linear interpolated
 * orientation between this and q.  If alpha is 0, then *this is unaltered.
 * If alpha is 1, then *this is set to q.
 * \param q the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 */
void Quat::slerp(const Quat& q, Real alpha)
{
  *this = slerp(*this, q, alpha);
}

/// Performs linear interpolation between this and q
/**
 * Sets the orientation represented by this quaternion to a linear interpolated
 * orientation between this and q.  If alpha is 0, then *this is unaltered.
 * If alpha is 1, then *this is set to q.
 * \param q the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 */
void Quat::lerp(const Quat& q, Real alpha)
{
  *this = lerp(*this, q, alpha);
}

/// Performs non-spherical linear interpolation between this and q
/**
 * Calculates the orientation of determined by linear interpolation between
 * q1 and q2.  If alpha is 0, then q1 is returned.  If alpha is 1, then 
 * q2 is returned. 
 * \param q1 the "initial" orientation
 * \param q2 the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 * \return the orientation linearly interpolated between q1 and q2 
 */
Quat Quat::lerp(const Quat& q1, const Quat& q2, Real alpha)
{
  if (alpha < 0.0 || alpha > 1.0)
    throw std::runtime_error("Attempting to interpolate using Quat::lerp() with t not in interval [0,1]");

  // compute q1'q2
  Real dot = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;

  // see whether we need to use the conjugate of q2
  if (dot < 0.0)  
  {
    Quat q = q1*(1.0 - alpha) - q2*alpha;
    q.normalize();
    return q;
  }
  else
  {
    Quat q = q1*(1.0 - alpha) + q2*alpha;
    q.normalize();
    return q;
  }
}

/// Performs spherical linear interpolation between this and q
/**
 * Calculates the orientation of determined by linear interpolation between
 * q1 and q2.  If alpha is 0, then q1 is returned.  If alpha is 1, then 
 * q2 is returned. 
 * \param q1 the "initial" orientation
 * \param q2 the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 * \return the orientation linearly interpolated between q1 and q2
 * \todo rewrite this function to avoid cancellation errors 
 */
Quat Quat::slerp(const Quat& q1, const Quat& q2, Real alpha)
{
  if (alpha < 0.0 || alpha > 1.0)
    throw std::runtime_error("Attempting to interpolate using Quat::slerp() with t not in interval [0,1]");

  // compute q1'q2
  Real dot = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;

  // see whether we need to use the conjugate of q2
  bool use_conj = (dot < 0.0);  

  // find the angle between the two
  Real theta = std::acos(std::fabs(dot));

  Real sin1at = std::sin((1-alpha)*theta);
  Real sinat = std::sin(alpha*theta);
  Real sint_i = 1.0/std::sin(theta);
  Quat qa(q1.x*sin1at, q1.y*sin1at, q1.z*sin1at, q1.w*sin1at);
  Quat qb(q2.x*sinat, q2.y*sinat, q2.z*sinat, q2.w*sinat);
  if (use_conj)
    qb = -qb;
  Quat qc = qa + qb;
  qc.x *= sint_i;
  qc.y *= sint_i;
  qc.z *= sint_i;
  qc.w *= sint_i;
  return qc; 
}

/// Computes the inverse orientation represented by this quaternion in place
void Quat::inverse()
{
  conjugate();
}

/// Computes the inverse orientation of a quaternion
Quat Quat::inverse(const Quat& q)
{
  return conjugate(q);
}

/// Normalizes this quaternion in place
void Quat::normalize()
{
  Real mag = magnitude();
  if (mag != 0.0)
  {
    Real magi = 1.0/mag;
    x *= magi;
    y *= magi;
    z *= magi;
    w *= magi;
  }
  else
  {
    x = y = z = 0;
    w = 1;
  }
}

/// Computes the normalized quaternion of q
Quat Quat::normalize(const Quat& q)
{
  Quat qn(q);
  qn.normalize();
  return qn;
}

/// Sets a quaternion from a vector
/**
 * \param v a four-dimensional vector (components 0..2 are the vector component, component 3 is the angular component)
 */
void Quat::set(const VectorN& v)
{
  // verify that the vector is four-dimensional
  if (v.size() != 4)
    throw MissizeException();
  
  // set the quaternion
  x = v[0];
  y = v[1];
  z = v[2];
  w = v[3];
}

/// Outputs a VectorN from a quaternion
VectorN Quat::to_vect(const Quat& q)
{
  VectorN v(4);
  v[0] = q.x;
  v[1] = q.y;
  v[2] = q.z;
  v[3] = q.w;
  return v;
}

/// Sets this to that represented by an axis-angle representation
void Quat::set(const AAngle* a)
{
  const Real half = a->angle*0.5;
  Real sina = std::sin(half);
  x = a->x * sina;
  y = a->y * sina;
  z = a->z * sina;
  w = std::cos(half);

  #ifndef NDEBUG
  if (!unit())
    std::cerr << "Quat::set() warning!  not a unit quaternion!" << std::endl;
  #endif
}

/// Determines whether this is a unit quaternion
bool Quat::unit() const
{
  Real mag = magnitude();
  return (std::fabs(mag - 1.0) < std::sqrt(NEAR_ZERO));
}

/// Sets this quaternion to that represented by a rotation matrix
/**
 * Code based on algorithm presented by Martin Baker
 */
void Quat::set(const Matrix3* m)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // core computation
  w = std::sqrt(std::max((Real) 0.0, (Real) 1.0 + (*m)(X,X) + (*m)(Y,Y) + (*m)(Z,Z))) * (Real) 0.5;
  x = std::sqrt(std::max((Real) 0.0, (Real) 1.0 + (*m)(X,X) - (*m)(Y,Y) - (*m)(Z,Z))) * (Real) 0.5;
  y = std::sqrt(std::max((Real) 0.0, (Real) 1.0 - (*m)(X,X) + (*m)(Y,Y) - (*m)(Z,Z))) * (Real) 0.5;
  z = std::sqrt(std::max((Real) 0.0, (Real) 1.0 - (*m)(X,X) - (*m)(Y,Y) + (*m)(Z,Z))) * (Real) 0.5;

  // sign computation
  if ((*m)(Z,Y) - (*m)(Y,Z) < (Real) 0.0)
    x = -x;
  if ((*m)(X,Z) - (*m)(Z,X) < (Real) 0.0)
    y = -y;
  if ((*m)(Y,X) - (*m)(X,Y) < (Real) 0.0)
    z = -z;

  #ifndef NDEBUG
  if (!unit())
    std::cerr << "Quat::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return;
}

/// Sets this quaternion to that represented by the rotation component of a transformation matrix
void Quat::set(const Matrix4* m)
{
  Matrix3 mr = m->get_rotation();
  set(&mr);
}

/// Constructs the conjugate of <b>this</b>
Quat Quat::operator-() const
{
  return Quat::conjugate(*this);
}

/// Subtracts a quaternion from <b>this</b>
Quat Quat::operator-(const Quat& q) const
{
  return Quat(x+q.x, y+q.y, z+q.z, w+q.w);
}

/// Subtracts one quaternion from another 
Quat& Quat::operator-=(const Quat& q)
{
  *this = *this - q;
  return *this;
}

/// Adds this quaternion to another and returns the result in a new quaternion
Quat Quat::operator+(const Quat& q) const
{
  return Quat(x+q.x, y+q.y, z+q.z, w+q.w);
}

/// Adds this quaon to another and stores the result in <b>this</b>
Quat& Quat::operator+=(const Quat& q)
{
  *this = *this + q;
  return *this;
}

/// Multiplies inv(q) by <b>this</b> and returns the result
Quat Quat::operator/(const Quat& q) const
{
  return conjugate(q) * (*this);
}

/// Multiplies <b>this</b> by q and returns the result 
Quat Quat::operator*(const Quat& q) const
{
  Quat qm;

  qm.w = w * q.w - x * q.x - y * q.y - z * q.z; 
  qm.x = w * q.x + x * q.w + y * q.z - z * q.y;
  qm.y = w * q.y + y * q.w + z * q.x - x * q.z;
  qm.z = w * q.z + z * q.w + x * q.y - y * q.x;

  return qm;
}

/// Multiplies a quaternion by a scalar
Quat Quat::operator*(Real scalar) const
{
  Quat qn;
  qn.x = x * scalar;
  qn.y = y * scalar;
  qn.z = z * scalar;
  qn.w = w * scalar;
  return qn;
}

/// Multiplies a quaternion by a scalar and stores the result in <b>this</b>
Quat& Quat::operator*=(Real scalar)
{
  x *= scalar;
  y *= scalar;
  z *= scalar;
  w *= scalar;
  return *this;
}

/// Multiplies the 3x3 matrix corresponding to this quaternion by a 3D vector and returns the result in a new 3D vector
Vector3 Quat::operator*(const Vector3& v) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real w2 = w*w;
  const Real x2 = x*x;
  const Real y2 = y*y;
  const Real z2 = z*z;
  const Real xy = x*y;
  const Real xz = x*z;
  const Real yz = y*z;
  const Real xw = x*w;
  const Real yw = y*w;
  const Real zw = z*w;
  return Vector3((-1.0+2.0*(w2+x2))*v[X] + 2.0*((xy-zw)*v[Y] + (yw+xz)*v[Z]), 2.0*((xy+zw)*v[X] + (-xw+yz)*v[Z]) + (-1.0+2.0*(w2+y2))*v[Y], 2.0*((-yw+xz)*v[X] + (xw+yz)*v[Y]) + (-1.0+2.0*(w2+z2))*v[Z]);
}

/// Multiplies <b>this</b> by q and stores the result in <b>this</b>
Quat& Quat::operator*=(const Quat& q)
{
  *this = *this * q;
  return *this;
}

/// Returns the l-infinity norm of the quaternion components
Real Quat::norm_inf() const
{
  return std::max(std::fabs(x), std::max(std::fabs(y), std::max(std::fabs(z),
                  std::fabs(w))));
}

/// Calculates the squared magnitude of a quaternion
Real Quat::norm_sq() const
{
  return w*w + x*x + y*y + z*z;
}

/// Calculates the magnitude of a quaternion
Real Quat::magnitude() const
{
  return safe_sqrt(w*w + x*x + y*y + z*z);
}

/// Computes the product between a quaternion and a scalar
Quat Moby::operator*(Real scalar, const Quat& q)
{
  return q * scalar;
}

/// Sends the quaternion to the output stream
std::ostream& Moby::operator<<(std::ostream& out, const Quat& q) 
{
  out << "<" << q.x << ", " << q.y << ", " << q.z << "> " << q.w;

  return out;
}
