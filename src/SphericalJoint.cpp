/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/NullPointerException.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/SphericalJoint.h>

using namespace Moby;
using boost::dynamic_pointer_cast;

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
SphericalJoint::SphericalJoint() : Joint()
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1] = ZEROS_3;
  _u[eAxis2] = ZEROS_3;
  _u[eAxis3] = ZEROS_3;

  // set rotation to z axis matrix to zero (for debugging purposes)
  _R = ZEROS_3x3;

  // set the translation in the transformation to zero
  _T.set_translation(ZEROS_3);

  // setup the spatial axis derivative to zero
  _si_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
  _s0_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());

  // assign the spherical joint tolerance
  SINGULAR_TOL = (double) 1e-2;
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
SphericalJoint::SphericalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1] = ZEROS_3;
  _u[eAxis2] = ZEROS_3;
  _u[eAxis3] = ZEROS_3;

  // set the translation induced by the joint to zero
  _T.set_translation(ZEROS_3);

  // set rotation to z axis matrix to zero (for debugging purposes)
  _R = ZEROS_3x3;

  // setup the spatial axis derivative to zero
  _si_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
  _s0_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());

  // assign the spherical joint tolerance
  SINGULAR_TOL = (double) 1e-2;
}  

/// Determines whether two values are relatively equal
bool SphericalJoint::rel_equal(double x, double y)
{
  return (std::fabs(x - y) <= NEAR_ZERO * std::max(std::fabs(x), std::max(std::fabs(y), (double) 1.0)));
}

/// Gets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
Vector3 SphericalJoint::get_axis_global(Axis a) const
{
  const unsigned X = 0;

  // get the inboard link
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "SphericalJoint::get_axis_global() - attempt to get axis w/o inboard link!" << std::endl;
    return ZEROS_3;
  }

  // get the transform for the inboard link
  Matrix3 R;
  inboard_link->get_transform().get_rotation(&R);
  
  // axis one is easy 
  if (a == eAxis1)
  {
    Vector3 axis;
    _R.get_column(X, axis.begin());
    return R * axis;
  }

  // for both axes 2 and 3 we need cos and sin of q(1)
  const double c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);

  // axis two is obtained by multiplying rotation matrix around x by axis [0,1,0]
  if (a == eAxis2)
    return R * _R * Vector3(0,c1,s1);
    
  // axis 3, requires the rotation matrix induced by the axis-angle
  // representation for axis 2..  much simpler to just use the rotation matrix from the
  // universal joint induced transform
  const double c2 = std::cos(q[DOF_2]+_q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+_q_tare[DOF_2]);
  assert (a == eAxis3);
  return R * _R * Vector3(s2, -c2*s1, c1*c2);    
}

/// Sets the local axis of translation for this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void SphericalJoint::set_axis_local(const Vector3& axis, Axis a) 
{ 
  // normalize the axis in case the caller did not
  Vector3 naxis = Vector3::normalize(axis);

  _u[a] = naxis; 
  update_spatial_axes(); 
}        

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \note set_axis_global() should be called on all axes before setting q! 
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void SphericalJoint::set_axis_global(const Vector3& axis, Axis a)
{
  // normalize the axis in case the caller did not
  Vector3 naxis = Vector3::normalize(axis);

  // get the inboard link
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "SphericalJoint::set_axis_global() - attempt to set axis w/o inboard link!" << std::endl;
    return;
  }

  // get the orientation for the inboard link
  Matrix3 R;
  inboard_link->get_transform().get_rotation(&R);

  // transform the axis into local coordinates
  _u[a] = Vector3::normalize(R * naxis);

  // update the spatial axis
  update_spatial_axes();
}

/// Updates the spatial axis for this joint
void SphericalJoint::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // call the parent method 
  Joint::update_spatial_axes();

  // assign axes, if possible
  if (!assign_axes())
    return;

  // verify that all are unit-length and they are orthogonal
  assert(std::fabs(_u[eAxis1].norm() - 1.0) < NEAR_ZERO);
  assert(std::fabs(_u[eAxis2].norm() - 1.0) < NEAR_ZERO);
  assert(std::fabs(_u[eAxis3].norm() - 1.0) < NEAR_ZERO);
  assert(std::fabs(Vector3::dot(_u[eAxis1], _u[eAxis2])) < NEAR_ZERO);
  assert(std::fabs(Vector3::dot(_u[eAxis1], _u[eAxis3])) < NEAR_ZERO);
  assert(std::fabs(Vector3::dot(_u[eAxis2], _u[eAxis3])) < NEAR_ZERO);

  // set the axes
  _R.set_column(X, _u[eAxis1]);
  _R.set_column(Y, _u[eAxis2]);
  _R.set_column(Z, _u[eAxis3]);
}

/// Attempts to compute unassigned axes
bool SphericalJoint::assign_axes()
{
  // compute the rotation matrix necessary to transform the first axis to the
  // z-axis; can only compute this matrix if all three axes are set
  if (_u[eAxis1].norm() < NEAR_ZERO)
  {
    // check whether axis2 and/or 3 are set
    if (_u[eAxis2].norm() < NEAR_ZERO)
    {
      // axis2 is not set; if axis3 is not set, we can do nothing...
      if (_u[eAxis3].norm() < NEAR_ZERO)
        return false;

      // otherwise, set axes1 and 2 from 3
      _u[eAxis3].normalize();
      Vector3::determine_orthonormal_basis(_u[eAxis3], _u[eAxis1], _u[eAxis2]);
    }
    else
    {
      // check for axis2 set and axes 1 and 3 not set
      if (_u[eAxis3].norm() < NEAR_ZERO)
      {
        // axis3 is not set; set axes1/3 from 2
        _u[eAxis2].normalize();
        Vector3::determine_orthonormal_basis(_u[eAxis2], _u[eAxis3], _u[eAxis1]);
      }
      else
      {
        // axis2 and axis3 are set; set axis1
        _u[eAxis2].normalize();
        _u[eAxis3].normalize();
        _u[eAxis1] = Vector3::cross(_u[eAxis2], _u[eAxis3]);      
      }
    }
  }
  else
  {
    // check whether axis2 and/or 3 are set
    if (_u[eAxis2].norm() < NEAR_ZERO)
    {
      // axis2 is not set; if axis3 is not set, set them arbitrarily
      if (_u[eAxis3].norm() < NEAR_ZERO)
      {
        _u[eAxis1].normalize();
        Vector3::determine_orthonormal_basis(_u[eAxis1], _u[eAxis2], _u[eAxis3]);
      }
      else
      {
        // axis3 is set; just need to set axis2
        _u[eAxis1].normalize();
        _u[eAxis3].normalize();
        _u[eAxis2] = Vector3::cross(_u[eAxis3], _u[eAxis1]);
      }
    }
    else
    {
      // axis1 and 2 are set; see whether we need to set axis3
      if (_u[eAxis3].norm() < NEAR_ZERO)
      {
        _u[eAxis1].normalize();        
        _u[eAxis2].normalize();        
        _u[eAxis3] = Vector3::cross(_u[eAxis1], _u[eAxis2]);
      }
    }
  }

  return true;
} 

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const SMatrix6N& SphericalJoint::get_spatial_axes(ReferenceFrameType rftype)
{
  const unsigned X = 0, Y = 1, Z = 2;

  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes() called with NULL outboard link");

  // get current values of q
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;

  // get the outboard link's joint to com vector in outer link coordinates
  const Vector3& p = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

  // get the set of spatial axes
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);

  // form untransformed spatial axes -- this are the vectors describing each axis, after
  // rotation by preceding axis/axes; note that first axis always points toward 1,0,0
  Vector3 uu2(0, c1, s1);
  Vector3 uu3(s2, -c2*s1, c1*c2);

  // transform the spatial axes into the outer link frame
  Vector3 u1;
  _R.get_column(X, u1.begin());
  Vector3 u2 = _R * uu2;
  Vector3 u3 = _R * uu3;

  // update the spatial axis in link coordinates
  SVector6 si1, si2, si3;
  si1.set_upper(u1);
  si2.set_upper(u2);
  si3.set_upper(u3);
  si1.set_lower(Vector3::cross(u1, p));
  si2.set_lower(Vector3::cross(u2, p));
  si3.set_lower(Vector3::cross(u3, p));
  _si.set_column(eAxis1, si1);
  _si.set_column(eAxis2, si2);
  _si.set_column(eAxis3, si3);

  // setup the complement of the spatial axes in link coordinates
  calc_s_bar_from_si();

  // use the Joint function to do the rest
  return Joint::get_spatial_axes(rftype);
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const SMatrix6N& SphericalJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get q, _q_tare, and qd
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;
  const VectorN& qd = this->qd;

  // get the two transformed axes
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  double qd1 = qd[DOF_1];
  double qd2 = qd[DOF_2];

  // form the time derivatives of the non-constant spatial axes (untransformed) 
  Vector3 uu2(0, -s1*qd1, c1*qd1);
  Vector3 uu3(c2*qd2, -c2*c1*qd1 + s2*s1*qd2, -c1*s2*qd2 - s1*c2*qd1);

  // transform the axes into outer link coordinates
  Vector3 u2 = _R * uu2; 
  Vector3 u3 = _R * uu3; 

  // get the outboard link's joint to com vector in outer link coordinates
  const Vector3& p = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

  // update the spatial axis in link coordinates; note that third column of spatial axis
  // derivative set to zero in constructor and is never modified
  SVector6 si2, si3;
  si2.set_upper(u2);
  si3.set_upper(u3);
  si2.set_lower(Vector3::cross(u2, p));
  si3.set_lower(Vector3::cross(u3, p));
  _si_dot.set_column(eAxis2, si2);
  _si_dot.set_column(eAxis3, si3);

  // transform to global coordinates
  SpatialTransform X_0_i = outboard->get_spatial_transform_link_to_global();
  X_0_i.transform(_si_dot, _s0_dot);

  return (rftype == eLink) ? _si_dot : _s0_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void SphericalJoint::determine_q(VectorN& q)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the inboard and outboard links
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();

  // verify that the inboard and outboard links are set
  if (!inboard || !outboard)
    throw std::runtime_error("determine_q() called on NULL inboard and/or outboard links!");

  // if any of the axes are not defined, can't use this method
  if (std::fabs(_u[0].norm_sq() - 1.0) > NEAR_ZERO ||
      std::fabs(_u[1].norm_sq() - 1.0) > NEAR_ZERO ||
      std::fabs(_u[2].norm_sq() - 1.0) > NEAR_ZERO)
    throw UndefinedAxisException();

  // set proper size for q
  q.resize(num_dof());

  // get the link transforms
  Matrix3 R_inboard, R_outboard;
  inboard->get_transform().get_rotation(&R_inboard);
  outboard->get_transform().get_rotation(&R_outboard);

  // determine the joint transformation
  Matrix3 R_local = R_inboard.transpose_mult(R_outboard);

  // back out the transformation to z-axis
  Matrix3 RU = _R.transpose_mult(R_local * _R);

  // determine cos and sin values for q1, q2,  and q3
  double s2 = RU(X,Z);
  double c2 = std::cos(std::asin(s2));
  double s1, c1, s3, c3;
  if (std::fabs(c2) > NEAR_ZERO)
  {
    s1 = -RU(Y,Z)/c2;
    c1 = RU(Z,Z)/c2;
    s3 = -RU(X,Y)/c2;
    c3 = RU(X,X)/c2;
    assert(!std::isnan(s1));
    assert(!std::isnan(c1));
    assert(!std::isnan(s3));
    assert(!std::isnan(c3));
  }
  else
  {
    // singular, we can pick any value for s1, c1, s3, c3 as long as the
    // following conditions are satisfied
    // c1*s3 + s1*c3*s2 = RU(Y,X)
    // c1*c3 - s1*s3*s2 = RU(Y,Y)
    // s1*s3 - c1*c3*s2 = RU(Z,X)
    // s1*c3 + c1*s3*s2 = RU(Z,Y)
    // so, we'll set q1 to zero (arbitrarily) and obtain
    s1 = 0;
    c1 = 1;
    s3 = RU(Y,X);
    c3 = RU(Y,Y);
  }

  // now determine q; only q2 can be determined without ambiguity
  if (std::fabs(s1) < NEAR_ZERO)
    q[DOF_2] = std::atan2(RU(X,Z), RU(Z,Z)/c1);
  else
    q[DOF_2] = std::atan2(RU(X,Z), -RU(Y,Z)/s1);
  assert(!std::isnan(q[DOF_2]));

  // if cos(q2) is not singular, proceed easily from here..
  if (std::fabs(c2) > NEAR_ZERO)
  {
    q[DOF_1] = std::atan2(-RU(Y,Z)/c2, RU(Z,Z)/c2);
    q[DOF_3] = std::atan2(-RU(X,Y)/c2, RU(X,X)/c2);
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
  else
  {
    if (std::fabs(c1) > NEAR_ZERO)
      q[DOF_3] = std::atan2((RU(Y,X) - s1*s2*c3)/c1, (RU(Y,Y) + s1*s2*s3)/c1);
    else
      q[DOF_3] = std::atan2((RU(Z,X) + c1*s2*c3)/s1, (RU(Z,Y) - c1*s2*s3)/s1);
    if (std::fabs(c3) > NEAR_ZERO)
      q[DOF_1] = std::atan2((RU(Y,X) - c1*s3)/(s2*c3), (-RU(Y,X) + s1*s3)/(s2*c3));
    else
      q[DOF_1] = std::atan2((-RU(Y,Y) + c1*c3)/(s2*s3), (RU(Z,Y) - s1*c3)/(s2*s3));
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
}

/// Gets the (local) rotation induced by this joint
Matrix3 SphericalJoint::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q, _q_tare
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;

  // compute some needed quantities
  const double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  const double c3 = std::cos(q[DOF_3]+q_tare[DOF_3]);
  const double s3 = std::sin(q[DOF_3]+q_tare[DOF_3]);

  // determine untransformed rotation
  // this is just the rotation matrix induced by using Tait-Bryan angles
  Matrix3 RU;
  RU(X,X) = c2*c3;              RU(X,Y) = -c2*s3;              RU(X,Z) = s2;
  RU(Y,X) = s1*s2*c3 + c1*s3;   RU(Y,Y) = -s1*s2*s3 + c1*c3;   RU(Y,Z) = -c2*s1;
  RU(Z,X) = -c1*s2*c3 + s1*s3;  RU(Z,Y) = c1*s2*s3 + s1*c3;    RU(Z,Z) = c2*c1;

  // transform the rotation
  Matrix3 R_local = _R * RU.mult_transpose(_R);

  return R_local;
}

/// Gets the (local) transform for this joint
const Matrix4& SphericalJoint::get_transform()
{
  // note that translation is set to zero in the constructors
  Matrix3 R = get_rotation();
  _T.set_rotation(&R);

  return _T;
}

/// Evaluates the constraint equations
void SphericalJoint::evaluate_constraints(double C[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], pp. 430-431; variable
  // names have been altered, however

  // determine the global positions of the attachment points and subtract them
  Vector3 r1 = get_position_global(false);
  Vector3 r2 = get_position_global(true);
  Vector3 r12 = r1 - r2; 

  // copy values
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
}

/// Computes the constraint jacobian with respect to a body
void SphericalJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup the constraint equations (from Shabana, p. 432)
  if (body == inner)
  {
    // get the information necessary to compute the constraint equations
    const Quat& q = inner->get_orientation();
    const Vector3& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
    const double qx = q.x;
    const double qy = q.y;
    const double qz = q.z;
    const double qw = q.w;
    const double px = p[X];
    const double py = p[Y];
    const double pz = p[Z];

    switch (index)
    {
      case 0:
        Cq[0] = 1.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*px*qw + 2*pz*qy - 2*py*qz;
        Cq[4] = 4*qx*px + 2*qy*py + 2*qz*pz;
        Cq[5] = 2*pz*qw + 2*py*qx;
        Cq[6] = 2*pz*qx - 2*py*qw;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 1.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*py*qw - 2*pz*qx + 2*px*qz;
        Cq[4] = 2*qy*px - 2*qw*pz;
        Cq[5] = 2*px*qx + 4*py*qy + 2*pz*qz;
        Cq[6] = 2*px*qw + 2*pz*qy;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 1.0;
        Cq[3] = 4*pz*qw + 2*py*qx - 2*px*qy;
        Cq[4] = 2*qz*px + 2*qw*py;
        Cq[5] = 2*py*qz - 2*px*qw;
        Cq[6] = 4*pz*qz + 2*py*qy + 2*px*qx;
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // get the information necessary to compute the constraint equations
    const Quat& q = outer->get_orientation();
    const Vector3& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
    const double qx = q.x;
    const double qy = q.y;
    const double qz = q.z;
    const double qw = q.w;
    const double px = -p[X];
    const double py = -p[Y];
    const double pz = -p[Z];

    switch (index)
    {
      case 0:
        Cq[0] = -1.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -(4*px*qw + 2*pz*qy - 2*py*qz);
        Cq[4] = -(4*qx*px + 2*qy*py + 2*qz*pz);
        Cq[5] = -(2*pz*qw + 2*py*qx);
        Cq[6] = -(2*pz*qx - 2*py*qw);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = -1.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*py*qw - 2*pz*qx + 2*px*qz);
        Cq[4] = -(2*qy*px - 2*qw*pz);
        Cq[5] = -(2*px*qx + 4*py*qy + 2*pz*qz);
        Cq[6] = -(2*px*qw + 2*pz*qy);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = -1.0;
        Cq[3] = -(4*pz*qw + 2*py*qx - 2*px*qy);
        Cq[4] = -(2*qz*px + 2*qw*py);
        Cq[5] = -(2*py*qz - 2*px*qw);
        Cq[6] = -(4*pz*qz + 2*py*qy + 2*px*qx);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Computes the time derivative of the constraint jacobian with respect to a body
void SphericalJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup the constraint equations (from Shabana, p. 432)
  if (body == inner)
  {
    // get the information necessary to compute the constraint equations
    const Quat& q = inner->get_orientation();
    const Vector3& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
    const double px = p[X];
    const double py = p[Y];
    const double pz = p[Z];
    Quat qd = Quat::deriv(q, inner->get_avel());
    const double dqw = qd.w;
    const double dqx = qd.x;
    const double dqy = qd.y;
    const double dqz = qd.z;

    switch (index)
    {
      case 0:
        Cq[0] = 0.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*px*dqw + 2*pz*dqy - 2*py*dqz;
        Cq[4] = 4*dqx*px + 2*dqy*py + 2*dqz*pz;
        Cq[5] = 2*pz*dqw + 2*py*dqx;
        Cq[6] = 2*pz*dqx - 2*py*dqw;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*py*dqw - 2*pz*dqx + 2*px*dqz;
        Cq[4] = 2*dqy*px - 2*dqw*pz;
        Cq[5] = 2*px*dqx + 4*py*dqy + 2*pz*dqz;
        Cq[6] = 2*px*dqw + 2*pz*dqy;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = 4*pz*dqw + 2*py*dqx - 2*px*dqy;
        Cq[4] = 2*dqz*px + 2*dqw*py;
        Cq[5] = 2*py*dqz - 2*px*dqw;
        Cq[6] = 4*pz*dqz + 2*py*dqy + 2*px*dqx;
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // get the information necessary to compute the constraint equations
    const Quat& q = outer->get_orientation();
    const Vector3& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
    const double px = -p[X];
    const double py = -p[Y];
    const double pz = -p[Z];
    Quat qd = Quat::deriv(q, outer->get_avel());
    const double dqw = qd.w;
    const double dqx = qd.x;
    const double dqy = qd.y;
    const double dqz = qd.z;

    switch (index)
    {
      case 0:
        Cq[0] = 0.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -(4*px*dqw + 2*pz*dqy - 2*py*dqz);
        Cq[4] = -(4*dqx*px + 2*dqy*py + 2*dqz*pz);
        Cq[5] = -(2*pz*dqw + 2*py*dqx);
        Cq[6] = -(2*pz*dqx - 2*py*dqw);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = 0.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*py*dqw - 2*pz*dqx + 2*px*dqz);
        Cq[4] = -(2*dqy*px - 2*dqw*pz);
        Cq[5] = -(2*px*dqx + 4*py*dqy + 2*pz*dqz);
        Cq[6] = -(2*px*dqw + 2*pz*dqy);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = -(4*pz*dqw + 2*py*dqx - 2*px*dqy);
        Cq[4] = -(2*dqz*px + 2*dqw*py);
        Cq[5] = -(2*py*dqz - 2*px*dqw);
        Cq[6] = -(4*pz*dqz + 2*py*dqy + 2*px*dqx);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Implements Base::load_from_xml()
void SphericalJoint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "SphericalJoint") == 0);

  // read local joint axes
  const XMLAttrib* laxis_attrib = node->get_attrib("local-axis");
  if (laxis_attrib)
  {
    Vector3 laxis;
    laxis_attrib->get_vector_value(laxis);
    set_axis_local(laxis, eAxis1);
  }
  const XMLAttrib* laxis1_attrib = node->get_attrib("local-axis1");
  if (laxis1_attrib)
  {
    Vector3 laxis1;
    laxis1_attrib->get_vector_value(laxis1);
    set_axis_local(laxis1, eAxis1);
  }
  const XMLAttrib* laxis2_attrib = node->get_attrib("local-axis2");
  if (laxis2_attrib)
  {
    Vector3 laxis2;
    laxis2_attrib->get_vector_value(laxis2);
    set_axis_local(laxis2, eAxis2);
  }

  // read the global joint axes, if given
  const XMLAttrib* gaxis_attrib = node->get_attrib("global-axis");
  if (gaxis_attrib)
  {
    Vector3 gaxis;
    gaxis_attrib->get_vector_value(gaxis);
    set_axis_global(gaxis, eAxis1);  
  }
  const XMLAttrib* gaxis1_attrib = node->get_attrib("global-axis1");
  if (gaxis1_attrib)
  {
    Vector3 gaxis1;
    gaxis1_attrib->get_vector_value(gaxis1);
    set_axis_global(gaxis1, eAxis1);  
  }
  const XMLAttrib* gaxis2_attrib = node->get_attrib("global-axis2");
  if (gaxis2_attrib)
  {
    Vector3 gaxis2;
    gaxis2_attrib->get_vector_value(gaxis2);
    set_axis_global(gaxis2, eAxis2);  
  }

  // compute _q_tare if necessary
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void SphericalJoint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "SphericalJoint";

  // save the local joint axes
  node->attribs.insert(XMLAttrib("local-axis1", _u[eAxis1]));
  node->attribs.insert(XMLAttrib("local-axis2", _u[eAxis1]));
}

