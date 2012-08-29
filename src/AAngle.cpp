/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <Moby/Constants.h>
#include <Moby/Quat.h>
#include <Moby/AAngle.h>
#include <Moby/MissizeException.h>

using namespace Moby;

/// Default constructor
AAngle::AAngle()
{
  x = y = z = angle = 0;
}

/// Constructs an axis-angle object from an array of four Real values
/**
 * \param array a 4-dimensional (or larger) array; the first three values are the axis and the fourth is the angle
 * \note automatically normalizes the axis
 */
AAngle::AAngle(boost::shared_array<const Real> array)
{
  set(array);
}

/// Constructs an axis-angle object from a quaternion
AAngle::AAngle(const Quat* q)
{
  set(q);
}

/// Constructs an axis-angle object from four Real values
/*
 * \note automatically normalizes the axis
 */
AAngle::AAngle(Real x, Real y, Real z, Real angle)
{
  set(x, y, z, angle);
}

/// Constructs an axis-angle object from a vector and a Real value
/*
 * \note automatically normalizes the axis
 */
AAngle::AAngle(const Vector3*  v, Real angle)
{
  set(v, angle);
}

/// Constructs an axis-angle object from a 4D vector (components 0-2 are axis, component 3 is angle)
/*
 * \note automatically normalizes the axis
 */
AAngle::AAngle(const VectorN* v)
{
  set(v);
}

/// Constructs an axis-angle object from a 3x3 rotation matrix
AAngle::AAngle(const Matrix3* m)
{
  set(m);
}

/// Constructs an axis-angle object from a 3x3 rotation matrix and a desired axis
/*
 * \note automatically normalizes the axis
 */
AAngle::AAngle(const Matrix3* m, const Vector3* v)
{
  set(m, v);
}

/// Constructs an axis-angle object from a 4x4 transformation matrix
AAngle::AAngle(const Matrix4* m)
{
  set(m);
}

/// Determines whether this axis-angle object is equal to another within some tolerance epsilon
bool AAngle::epsilon_equals(const AAngle& a, Real epsilon) const
{
  return epsilon_equals(*this, a, epsilon);
}

/// Determines whether two axis-angle objects are equal to within some tolerance epsilon
bool AAngle::epsilon_equals(const AAngle& a1, const AAngle& a2, Real epsilon)
{
  // convert to quaternion and test there
  boost::shared_ptr<Quat> q1(new Quat(&a1));
  boost::shared_ptr<Quat> q2(new Quat(&a2));
  bool flag = q1->epsilon_equals(*q2, epsilon);
  return flag;
}

/// Sets this object from four Real values
void AAngle::set(Real x, Real y, Real z, Real angle)
{
  Real ilen = 1.0/safe_sqrt(x*x + y*y + z*z);
  this->x = x*ilen;
  this->y = y*ilen;
  this->z = z*ilen;
  this->angle = angle;
}

/// Constructs an axis-angle object from a 4D vector (components 0-2 are axis, component 3 is angle)
/**
 * \note automatically normalizes the axis
 * Prints an error message to stderr if v is not a 4D vector
 */
void AAngle::set(const VectorN* v)
{
  if (v->size() != 4)
    throw MissizeException();
  x = (*v)[0];
  y = (*v)[1];
  z = (*v)[2];
  Real ilen = 1.0/safe_sqrt(x*x + y*y + z*z);
  x *= ilen;
  y *= ilen;
  z *= ilen;
  angle = (*v)[3];
}

/// Sets this object from a three-dimensional vector and a Real value
/**
 * \note automatically normalizes the axis
 */
void AAngle::set(const Vector3* v, Real angle)
{
  Vector3 w = Vector3::normalize(*v);
  x = w[0];
  y = w[1];
  z = w[2];
  this->angle = angle;
}

/// Sets the object from a rotation matrix and rotation axis
/**
 * \note automatically normalizes the axis
 */
void AAngle::set(const Matrix3* m, const Vector3* axis)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // normalize the axis
  Vector3 axisn = Vector3::normalize(*axis);

  // get the arc-cosine of the angle (clip it as necessary)
  Real acosangle = ((*m)(X,X) + (*m)(Y,Y) + (*m)(Z,Z) - 1)/2; 
  if (acosangle > 1.0)
    acosangle = 1.0;
  else if (acosangle < -1.0)
    acosangle = -1.0;

  // compute angle of rotation
  angle = std::acos(acosangle);

  // set axis of rotation
  x = axisn[X];
  y = axisn[Y];
  z = axisn[Z];

  // if the angle of rotation is zero, can return now
  if (angle == 0.0)
    return;

  // if the angle of rotation is PI, solve for axis
  Real x,y,z;
  if (angle == M_PI)
  {
    x = safe_sqrt(((*m)(X,X)+1)/2);
    y = safe_sqrt(((*m)(Y,Y)+1)/2);
    z = safe_sqrt(((*m)(Z,Z)+1)/2);
  }
  else
  {
    Real constant = 1.0/(2.0*std::sin(angle));
    x = constant * ((*m)(Z,Y) - (*m)(Y,Z));
    y = constant * ((*m)(X,Z) - (*m)(Z,X));
    z = constant * ((*m)(Y,X) - (*m)(X,Y));
  }

  // make sure that x, y, and z are not NaN
  assert(!std::isnan(x));
  assert(!std::isnan(y));
  assert(!std::isnan(z));

  // get the length of determined axis [x y z]; if length is zero, angle is
  // zero...
  Real len = safe_sqrt(x*x + y*y + z*z);
  if (len == 0.0)
  {
    angle = 0.0;
    return;
  }

  // normalize vector [x y z] (generally not necessary, but safe...)
  if (std::fabs(len - 1.0) > NEAR_ZERO)
  {
    Real ilen = 1.0/len;
    x *= ilen;
    y *= ilen;
    z *= ilen;
  }

  // form the determined axis and verify it points in same or exact opposite
  // direction as axis
  Vector3 determine_axis(x,y,z);
  
  // reverse the angle, if necessary
  if (Vector3::dot(determine_axis, axisn) < 0)
    angle = -angle;
}

/// Sets the object from a rotation matrix
void AAngle::set(const Matrix3* m)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // get the arc-cosine of the angle (clip it as necessary)
  Real acosangle = ((*m)(X,X) + (*m)(Y,Y) + (*m)(Z,Z) - 1)/2; 
  if (acosangle > 1.0)
    acosangle = 1.0;
  else if (acosangle < -1.0)
    acosangle = -1.0;

  // compute angle of rotation
  angle = std::acos(acosangle);
  
  // if angle is 0, then axis is arbitrary
  if (angle < std::numeric_limits<Real>::epsilon())
  {
    x = 1;
    y = z = 0;
    return;
  }
  
  // if angle is pi then must solve for rx, ry, rz
  if (std::fabs(angle-M_PI) < std::numeric_limits<Real>::epsilon())
  {
    x = safe_sqrt(((*m)(X,X)+1)/2);
    y = safe_sqrt(((*m)(Y,Y)+1)/2);
    z = safe_sqrt(((*m)(Z,Z)+1)/2);
    assert(!std::isnan(x));
    assert(!std::isnan(y));
    assert(!std::isnan(z));
    return;
  }
  
  // standard case
  Real constant = 1.0/(2.0*std::sin(angle));
  x = constant * ((*m)(Z,Y) - (*m)(Y,Z));
  y = constant * ((*m)(X,Z) - (*m)(Z,X));
  z = constant * ((*m)(Y,X) - (*m)(X,Y));
  
  // normalize the axis (generally not necessary, but safe...)
  Real len = safe_sqrt(x * x + y * y + z * z);
  assert(len != (Real) 0.0);
  if (std::fabs(len - (Real) 1.0) > NEAR_ZERO)
  {
    Real ilen = (Real) 1.0/len;
    x *= ilen;
    y *= ilen;
    z *= ilen;
  }

  assert(!std::isnan(angle));
  assert(!std::isnan(x));
  assert(!std::isnan(y));
  assert(!std::isnan(z));
}

/// Sets the object from a transformation matrix
void AAngle::set(const Matrix4* m)
{
  // get the rotation component
  Matrix3 r = m->get_rotation();
  set(&r);
}

/// Sets this object from an array of four Reals
/**
 * \param array a 4-dimensional (or larger) array; the first three values are the axis and the fourth is the angle
 * \note automatically normalizes the axis
 */
void AAngle::set(boost::shared_array<const Real> array)
{
  Real ilen = 1.0/safe_sqrt(array[0]*array[0]+array[1]*array[1]+array[2]*array[2]);
  x = array[0]*ilen;
  y = array[1]*ilen;
  z = array[2]*ilen;
  angle = array[3];
}

/// Sets this object from a quaternion object
/**
 * \todo test this method
 */
void AAngle::set(const Quat* q)
{
  Quat qn = Quat::normalize(*q);
  x = qn.x;
  y = qn.y;
  z = qn.z;
  angle = std::acos(qn.w) * 2.0;

  // verify that [x y z] normalized
  Real nrm = safe_sqrt(x*x + y*y + z*z);
  if (std::fabs(nrm) < NEAR_ZERO)
  {
    // axis is zero; set it to [1 0 0] arbitrarily
    x = 1.0;
    y = z = 0.0;
  }
  else if (std::fabs(nrm - 1.0) > NEAR_ZERO)
  {
    x /= nrm;
    y /= nrm;
    z /= nrm;
  }
}

/// Copies a axis-angle object
void AAngle::operator=(const AAngle& a)
{
  x = a.x;
  y = a.y;
  z = a.z;
  angle = a.angle;
}

/// Multiplies two axis-angle representations
/**
 * This has the effect of performing a rotation by a, then performing a 
 * rotation by <b>this</b>.
 */
AAngle AAngle::operator*(const AAngle& a) const
{
  boost::shared_ptr<Quat> q1(new Quat(this));
  boost::shared_ptr<Quat> q2(new Quat(&a));
  *q1 *= *q2;
  return AAngle(q1.get());
}

/// Multiplies this axis-angle representation by another and stores the result in <b>this</b>
/**
 * This has the effect of performing a rotation by a, then performing a
 * rotation by <b>this</b>.
 */
void AAngle::operator*=(const AAngle& a)
{
  boost::shared_ptr<Quat> q1(new Quat(this));
  boost::shared_ptr<Quat> q2(new Quat(&a));
  *q1 *= *q2;
  set(q1.get());
}

/// Sends the representation to the specified stream
std::ostream& Moby::operator<<(std::ostream& out, const AAngle& a)
{
  out << "[ " << a.x << ' ' << a.y << ' ' << a.z << " ] "  << a.angle << "  ";
  return out;
}

