/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/NullPointerException.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/SphericalJoint.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using boost::dynamic_pointer_cast;
using namespace Moby;

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
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _u[eAxis3].pose = _F;

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;

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
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _u[eAxis3].pose = _F;

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;

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
Vector3d SphericalJoint::get_axis(Axis a) const
{
  const unsigned X = 0;

  // axis one is easy 
  if (a == eAxis1)
  {
    return _u[0];
  }

  // for both axes 2 and 3 we need cos and sin of q(1)
  const double c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);

  // axis two is obtained by rotating around x axis
  if (a == eAxis2)
  {
    Vector3d u(0,c1,s1);
    u.pose = _u[0].pose;
    return u;
  }
    
  // axis 3, requires the rotation matrix induced by the axis-angle
  // representation for axis 2..  much simpler to just use the rotation matrix from the
  // universal joint induced transform
  const double c2 = std::cos(q[DOF_2]+_q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+_q_tare[DOF_2]);
  assert (a == eAxis3);
  Vector3d u(s2, -c2*s1, c1*c2);
  u.pose = _u[0].pose;
  return u;
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
void SphericalJoint::set_axis(const Vector3d& axis, Axis a) 
{ 
  // normalize the axis in case the caller did not
  Vector3d naxis = Vector3d::normalize(axis);
  _u[a] = Pose3d::transform(naxis.pose, get_pose(), naxis); 
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
  assert(std::fabs(Vector3d::dot(_u[eAxis1], _u[eAxis2])) < NEAR_ZERO);
  assert(std::fabs(Vector3d::dot(_u[eAxis1], _u[eAxis3])) < NEAR_ZERO);
  assert(std::fabs(Vector3d::dot(_u[eAxis2], _u[eAxis3])) < NEAR_ZERO);
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
      Vector3d::determine_orthonormal_basis(_u[eAxis3], _u[eAxis1], _u[eAxis2]);
    }
    else
    {
      // check for axis2 set and axes 1 and 3 not set
      if (_u[eAxis3].norm() < NEAR_ZERO)
      {
        // axis3 is not set; set axes1/3 from 2
        _u[eAxis2].normalize();
        Vector3d::determine_orthonormal_basis(_u[eAxis2], _u[eAxis3], _u[eAxis1]);
      }
      else
      {
        // axis2 and axis3 are set; set axis1
        _u[eAxis2].normalize();
        _u[eAxis3].normalize();
        _u[eAxis1] = Vector3d::cross(_u[eAxis2], _u[eAxis3]);      
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
        Vector3d::determine_orthonormal_basis(_u[eAxis1], _u[eAxis2], _u[eAxis3]);
      }
      else
      {
        // axis3 is set; just need to set axis2
        _u[eAxis1].normalize();
        _u[eAxis3].normalize();
        _u[eAxis2] = Vector3d::cross(_u[eAxis3], _u[eAxis1]);
      }
    }
    else
    {
      // axis1 and 2 are set; see whether we need to set axis3
      if (_u[eAxis3].norm() < NEAR_ZERO)
      {
        _u[eAxis1].normalize();        
        _u[eAxis2].normalize();        
        _u[eAxis3] = Vector3d::cross(_u[eAxis1], _u[eAxis2]);
      }
    }
  }

  return true;
} 

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<Twistd>& SphericalJoint::get_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;

  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes() called with NULL outboard link");

  // get current values of q
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;

  // get the outboard link's joint to com vector in outer link coordinates
  const Vector3d& p = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

  // get the set of spatial axes
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);

  // form untransformed spatial axes -- this are the vectors describing each axis, after
  // rotation by preceding axis/axes; note that first axis always points toward 1,0,0
  Vector3d uu2(0, c1, s1);
  Vector3d uu3(s2, -c2*s1, c1*c2);

  // transform the spatial axes into the joint frame 
  Vector3d u1 = _u[0];
  Vector3d u2 =  uu2;
  Vector3d u3 =  uu3;

  // setup relative poses for all three
  u1.pose = _F;
  u2.pose = _F;
  u3.pose = _F;

  // update the spatial axis in link coordinates
  _s[0].set_angular(u1);
  _s[0].set_linear(ZEROS_3);
  _s[1].set_angular(u2);
  _s[1].set_linear(ZEROS_3);
  _s[2].set_angular(u3);
  _s[2].set_linear(ZEROS_3);

  // setup the complement of the spatial axes in link coordinates
  calc_s_bar_from_s();

  // use the Joint function to do the rest
  return Joint::get_spatial_axes();
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<Twistd>& SphericalJoint::get_spatial_axes_dot()
{
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("SphericalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get q, _q_tare, and qd
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;
  const VectorNd& qd = this->qd;

  // get the two transformed axes
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  double qd1 = qd[DOF_1];
  double qd2 = qd[DOF_2];

  // form the time derivatives of the non-constant spatial axes (untransformed) 
  Vector3d uu2(0, -s1*qd1, c1*qd1);
  Vector3d uu3(c2*qd2, -c2*c1*qd1 + s2*s1*qd2, -c1*s2*qd2 - s1*c2*qd1);

  // transform the axes into outer link coordinates
  Vector3d u2 = uu2; 
  Vector3d u3 = uu3; 

  // update the spatial axis in joint coordinates; note that third column of spatial axis
  // derivative set to zero in constructor and is never modified
  _s_dot[1].set_angular(u2);
  _s_dot[1].set_lower(ZEROS_3);
  _s_dot[2].set_angular(u3);
  _s_dot[2].set_linear(ZEROS_3);

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void SphericalJoint::determine_q(VectorNd& q)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the outboard link
  RigidBodyPtr outboard = get_outboard_link();

  // verify that the outboard link is set
  if (!outboard)
    throw std::runtime_error("determine_q() called on NULL outboard link!");

  // if any of the axes are not defined, can't use this method
  if (std::fabs(_u[0].norm_sq() - 1.0) > NEAR_ZERO ||
      std::fabs(_u[1].norm_sq() - 1.0) > NEAR_ZERO ||
      std::fabs(_u[2].norm_sq() - 1.0) > NEAR_ZERO)
    throw UndefinedAxisException();

  // set proper size for q
  q.resize(num_dof());

  // get the poses of the joint and outboard link
  shared_ptr<const Pose3d> Fj = get_pose();
  shared_ptr<const Pose3d> Fo = outboard->get_pose();
  shared_ptr<Pose3d> Fjprime;

  // Fo will be relative to Fj' which will be relative to Fj

  // determine the joint transformation
  Matrix3d R = Fjprime->q;

  // determine cos and sin values for q1, q2,  and q3
  double s2 = R(X,Z);
  double c2 = std::cos(std::asin(s2));
  double s1, c1, s3, c3;
  if (std::fabs(c2) > NEAR_ZERO)
  {
    s1 = -R(Y,Z)/c2;
    c1 = R(Z,Z)/c2;
    s3 = -R(X,Y)/c2;
    c3 = R(X,X)/c2;
    assert(!std::isnan(s1));
    assert(!std::isnan(c1));
    assert(!std::isnan(s3));
    assert(!std::isnan(c3));
  }
  else
  {
    // singular, we can pick any value for s1, c1, s3, c3 as long as the
    // following conditions are satisfied
    // c1*s3 + s1*c3*s2 = R(Y,X)
    // c1*c3 - s1*s3*s2 = R(Y,Y)
    // s1*s3 - c1*c3*s2 = R(Z,X)
    // s1*c3 + c1*s3*s2 = R(Z,Y)
    // so, we'll set q1 to zero (arbitrarily) and obtain
    s1 = 0;
    c1 = 1;
    s3 = R(Y,X);
    c3 = R(Y,Y);
  }

  // now determine q; only q2 can be determined without ambiguity
  if (std::fabs(s1) < NEAR_ZERO)
    q[DOF_2] = std::atan2(R(X,Z), R(Z,Z)/c1);
  else
    q[DOF_2] = std::atan2(R(X,Z), -R(Y,Z)/s1);
  assert(!std::isnan(q[DOF_2]));

  // if cos(q2) is not singular, proceed easily from here..
  if (std::fabs(c2) > NEAR_ZERO)
  {
    q[DOF_1] = std::atan2(-R(Y,Z)/c2, R(Z,Z)/c2);
    q[DOF_3] = std::atan2(-R(X,Y)/c2, R(X,X)/c2);
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
  else
  {
    if (std::fabs(c1) > NEAR_ZERO)
      q[DOF_3] = std::atan2((R(Y,X) - s1*s2*c3)/c1, (R(Y,Y) + s1*s2*s3)/c1);
    else
      q[DOF_3] = std::atan2((R(Z,X) + c1*s2*c3)/s1, (R(Z,Y) - c1*s2*s3)/s1);
    if (std::fabs(c3) > NEAR_ZERO)
      q[DOF_1] = std::atan2((R(Y,X) - c1*s3)/(s2*c3), (-R(Y,X) + s1*s3)/(s2*c3));
    else
      q[DOF_1] = std::atan2((-R(Y,Y) + c1*c3)/(s2*s3), (R(Z,Y) - s1*c3)/(s2*s3));
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
}

/// Gets the (local) rotation induced by this joint
Matrix3d SphericalJoint::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q, _q_tare
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;

  // compute some needed quantities
  const double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  const double c3 = std::cos(q[DOF_3]+q_tare[DOF_3]);
  const double s3 = std::sin(q[DOF_3]+q_tare[DOF_3]);

  // determine rotation
  // this is just the rotation matrix induced by using Tait-Bryan angles
  Matrix3d R;
  R(X,X) = c2*c3;              R(X,Y) = -c2*s3;              R(X,Z) = s2;
  R(Y,X) = s1*s2*c3 + c1*s3;   R(Y,Y) = -s1*s2*s3 + c1*c3;   R(Y,Z) = -c2*s1;
  R(Z,X) = -c1*s2*c3 + s1*s3;  R(Z,Y) = c1*s2*s3 + s1*c3;    R(Z,Z) = c2*c1;

  return R;
}

/// Gets the (local) transform for this joint
shared_ptr<const Pose3d> SphericalJoint::get_induced_pose()
{
  // note that translation is zero by default 
  _Fprime->q = get_rotation();
  return _Fprime;
}

/*
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
  Vector3d r1 = get_position_global(false);
  Vector3d r2 = get_position_global(true);
  Vector3d r12 = r1 - r2; 

  // copy values
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
}

/// Computes the constraint jacobian with respect to a body
void SphericalJoint::calc_constraint_jacobian_euler(RigidBodyPtr body, unsigned index, double Cq[7])
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
    const Vector3d& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
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
    const Vector3d& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
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
void SphericalJoint::calc_constraint_jacobian_dot_euler(RigidBodyPtr body, unsigned index, double Cq[7])
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
    const Vector3d& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
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
    const Vector3d& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
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
*/

/// Implements Base::load_from_xml()
void SphericalJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "SphericalJoint") == 0);

  // read the global joint axes, if given
  const XMLAttrib* axis_attrib = node->get_attrib("axis");
  if (axis_attrib)
  {
    Vector3d axis;
    axis_attrib->get_vector_value(axis);
    set_axis(axis, eAxis1);  
  }
  const XMLAttrib* axis1_attrib = node->get_attrib("axis1");
  if (axis1_attrib)
  {
    Vector3d axis1;
    axis1_attrib->get_vector_value(axis1);
    set_axis(axis1, eAxis1);  
  }
  const XMLAttrib* axis2_attrib = node->get_attrib("axis2");
  if (axis2_attrib)
  {
    Vector3d axis2;
    axis2_attrib->get_vector_value(axis2);
    set_axis(axis2, eAxis2);  
  }

  // compute _q_tare if necessary
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void SphericalJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  Vector3d u0;

  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "SphericalJoint";

  // convert local axes to global axes and save
  u0 = Pose3d::transform(get_pose(), shared_ptr<const Pose3d>(), _u[eAxis1]);
  node->attribs.insert(XMLAttrib("axis1", u0));
  u0 = Pose3d::transform(get_pose(), shared_ptr<const Pose3d>(), _u[eAxis2]);
  node->attribs.insert(XMLAttrib("axis2", u0));
}

