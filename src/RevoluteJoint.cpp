/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/RevoluteJoint.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The axis of rotation is set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
RevoluteJoint::RevoluteJoint() : Joint()
{
  // init joint data
  init_data();

  // init the joint axes
  _u.set_zero(_F);
  _ui.set_zero();
  _uj.set_zero();
  _v2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].  
 */
RevoluteJoint::RevoluteJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  // init joint data
  init_data();

  // init the joint axes
  _u.set_zero();
  _ui.set_zero();
  _uj.set_zero();
  _v2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}  

/// Sets the axis of rotation for this joint
void RevoluteJoint::set_axis(const Vector3d& axis) 
{ 
  // normalize the joint axis, in case the caller didn't 
  Vector3d naxis = Vector3d::normalize(axis);

  // transform the axis as necessary
  _u = Pose3d::transform_vector(_F, naxis);

  // update the spatial axes
  update_spatial_axes(); 

  // setup ui and uj
  Vector3d::determine_orthonormal_basis(_u, _ui, _uj);
 
  // get the two links
  RigidBodyPtr b1 = get_inboard_link();
  RigidBodyPtr b2 = get_outboard_link();
  if (!b1 || !b2)
    throw std::runtime_error("Attempt to set joint axis without setting inboard and outboard links first!");
 
  // compute joint axis in outer link frame 
  _v2 = Pose3d::transform_vector(b2->get_pose(), naxis); 
}        

/// Updates the spatial axis for this joint
void RevoluteJoint::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Vector3d ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // call parent method
  Joint::update_spatial_axes();

  // update the spatial axis in joint coordinates
  _s[0].set_angular(_u);
  _s[0].set_linear(ZEROS_3);

  // setup s_bar
  calc_s_bar_from_s();
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void RevoluteJoint::determine_q(VectorNd& q)
{
  // get the outboard link pointer
  RigidBodyPtr outboard = get_outboard_link();
  
  // verify that the outboard link is set
  if (!outboard)
    throw std::runtime_error("determine_q() called on NULL outboard link!");

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - 1.0) > NEAR_ZERO)
    throw UndefinedAxisException();

  // get the poses of the joint and outboard link
  shared_ptr<const Pose3d> Fj = get_pose();
  shared_ptr<const Pose3d> Fo = outboard->get_pose();

  // compute transforms
  Transform3d wTo = Pose3d::calc_relative_pose(Fo, GLOBAL); 
  Transform3d jTw = Pose3d::calc_relative_pose(GLOBAL, Fj);
  Transform3d jTo = jTw * wTo;

  // determine the joint transformation
  Matrix3d R = jTo.q;
  AAngled aa(R, _u);

  // set q 
  q.resize(num_dof());
  q[DOF_1] = aa.angle;
}

/// Gets the pose for this joint
shared_ptr<const Pose3d> RevoluteJoint::get_induced_pose()
{
  assert(std::fabs(_u.norm() - 1.0) < NEAR_ZERO);

  // invalidate pose quantities for the outer link
  invalidate_pose_vectors();

  // note that translation is set to zero in the constructors
  _Fprime->q = AAngled(_u, this->q[DOF_1]+this->_q_tare[DOF_1]);

  return _Fprime;
}

/// Gets the derivative for the spatial axes for this joint
const std::vector<SVelocityd>& RevoluteJoint::get_spatial_axes_dot()
{
  return _s_dot;
}

/// Computes the constraint Jacobian with respect to a body
void RevoluteJoint::calc_constraint_jacobian(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr b1 = get_inboard_link();
  RigidBodyPtr b2 = get_outboard_link();

  // make sure that _u (and by extension _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Revolute joint axis has not been set; set before calling dynamics functions.");

  // mke sure that body is one of the links
  if (b1 != body && b2 != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup constants for calculations
  const Quatd& q1 = b1->get_orientation();
  const Quatd& q2 = b2->get_orientation();
  const Vector3d& p1 = b1->get_outer_joint_data(b2).com_to_joint_vec;
  const Vector3d& p2 = b2->get_inner_joint_data(b1).joint_to_com_vec_of;
  const double p1x = p1[X];
  const double p1y = p1[Y];
  const double p1z = p1[Z];
  const double p2x = -p2[X];
  const double p2y = -p2[Y];
  const double p2z = -p2[Z];
  const double qw1 = q1.w;
  const double qx1 = q1.x;
  const double qy1 = q1.y;
  const double qz1 = q1.z;
  const double qw2 = q2.w;
  const double qx2 = q2.x;
  const double qy2 = q2.y;
  const double qz2 = q2.z;
  const double uix = _ui[X];
  const double uiy = _ui[Y];
  const double uiz = _ui[Z];
  const double ujx = _uj[X];
  const double ujy = _uj[Y];
  const double ujz = _uj[Z];
  const double v2x = _v2[X];
  const double v2y = _v2[Y];
  const double v2z = _v2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == b1)
  {
    switch (index)
    {
      case 0:
        Cq[0] = 1.0;     
        Cq[1] = 0.0;     
        Cq[2] = 0.0;     
        Cq[3] = 4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1; 
        Cq[4] = 4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1; 
        Cq[5] = 2*p1z*qw1 + 2*p1y*qx1; 
        Cq[6] = -2*p1y*qw1 + 2*p1z*qx1; 
        break;

      case 1:
        Cq[0] = 0.0;     
        Cq[1] = 1.0;     
        Cq[2] = 0.0;     
        Cq[3] = 4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1; 
        Cq[4] = -2*p1z*qw1 + 2*p1x*qy1; 
        Cq[5] = 2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1; 
        Cq[6] = 2*p1x*qw1 + 2*p1z*qy1; 
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 1.0;
        Cq[3] = 4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1;
        Cq[4] = 2*p1y*qw1 + 2*p1x*qz1;
        Cq[5] = -2*p1x*qw1 + 2*p1y*qz1;
        Cq[6] = 2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1;
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (4*qw1*uix - 2*qz1*uiy + 2*qy1*uiz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qz1*uix + 4*qw1*uiy - 2*qx1*uiz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (-2*qy1*uix + 2*qx1*uiy + 4*qw1*uiz)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[4] = (4*qx1*uix + 2*qy1*uiy + 2*qz1*uiz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qy1*uix - 2*qw1*uiz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (2*qz1*uix + 2*qw1*uiy)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[5] = (2*qx1*uiy + 2*qw1*uiz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qx1*uix + 4*qy1*uiy + 2*qz1*uiz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (-2*qw1*uix + 2*qz1*uiy)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[6] = (-2*qw1*uiy + 2*qx1*uiz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qw1*uix + 2*qy1*uiz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (2*qx1*uix + 2*qy1*uiy + 4*qz1*uiz)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        break;
 
      case 4:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (4*qw1*ujx - 2*qz1*ujy + 2*qy1*ujz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qz1*ujx + 4*qw1*ujy - 2*qx1*ujz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (-2*qy1*ujx + 2*qx1*ujy + 4*qw1*ujz)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[4] = (4*qx1*ujx + 2*qy1*ujy + 2*qz1*ujz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qy1*ujx - 2*qw1*ujz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (2*qz1*ujx + 2*qw1*ujy)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[5] = (2*qx1*ujy + 2*qw1*ujz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qx1*ujx + 4*qy1*ujy + 2*qz1*ujz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (-2*qw1*ujx + 2*qz1*ujy)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[6] = (-2*qw1*ujy + 2*qx1*ujz)*
     ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
       2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
    (2*qw1*ujx + 2*qy1*ujz)*
     (2*(qx2*qy2 + qw2*qz2)*v2x + 
       (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
       2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
    (2*qx1*ujx + 2*qy1*ujy + 4*qz1*ujz)*
     (2*(-(qw2*qy2) + qx2*qz2)*v2x + 
       2*(qw2*qx2 + qy2*qz2)*v2y + 
       (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    switch (index)
    {
      case 0:
        Cq[0] = -1.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2; 
        Cq[4] = -4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2; 
        Cq[5] = -2*p2z*qw2 - 2*p2y*qx2; 
        Cq[6] = 2*p2y*qw2 - 2*p2z*qx2; 
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = -1.0;     
        Cq[2] = 0.0;      
        Cq[3] = -4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2; 
        Cq[4] = 2*p2z*qw2 - 2*p2x*qy2; 
        Cq[5] = -2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2; 
        Cq[6] = -2*p2x*qw2 - 2*p2z*qy2; 
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = -1.0;
        Cq[3] = -4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2;
        Cq[4] = -2*p2y*qw2 - 2*p2x*qz2;
        Cq[5] = 2*p2x*qw2 - 2*p2y*qz2;
        Cq[6] = -2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2;
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
     (-2*qy2*v2x + 2*qx2*v2y + 4*qw2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz)*
     (2*qz2*v2x + 4*qw2*v2y - 2*qx2*v2z) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
     (4*qw2*v2x - 2*qz2*v2y + 2*qy2*v2z);
        Cq[4] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
     (2*qz2*v2x + 2*qw2*v2y) + 
    (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*qy2*v2x - 2*qw2*v2z) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
     (4*qx2*v2x + 2*qy2*v2y + 2*qz2*v2z);
        Cq[5] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
     (-2*qw2*v2x + 2*qz2*v2y) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
     (2*qx2*v2y + 2*qw2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz)*
     (2*qx2*v2x + 4*qy2*v2y + 2*qz2*v2z);
        Cq[6] = ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
     (-2*qw2*v2y + 2*qx2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*qw2*v2x + 2*qy2*v2z) + 
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
     (2*qx2*v2x + 2*qy2*v2y + 4*qz2*v2z);
        break;
 
      case 4:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
     (-2*qy2*v2x + 2*qx2*v2y + 4*qw2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz)*
     (2*qz2*v2x + 4*qw2*v2y - 2*qx2*v2z) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
     (4*qw2*v2x - 2*qz2*v2y + 2*qy2*v2z);
        Cq[4] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
     (2*qz2*v2x + 2*qw2*v2y) + 
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz)*(2*qy2*v2x - 2*qw2*v2z) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
     (4*qx2*v2x + 2*qy2*v2y + 2*qz2*v2z);
        Cq[5] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
     (-2*qw2*v2x + 2*qz2*v2y) + 
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
     (2*qx2*v2y + 2*qw2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz)*
     (2*qx2*v2x + 4*qy2*v2y + 2*qz2*v2z);
        Cq[6] = ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
     (-2*qw2*v2y + 2*qx2*v2z) + 
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz)*(2*qw2*v2x + 2*qy2*v2z) + 
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
     (2*qx2*v2x + 2*qy2*v2y + 4*qz2*v2z);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
*/
}

/// Computes the time derivative of the constraint Jacobian with respect to a body
void RevoluteJoint::calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr b1 = get_inboard_link();
  RigidBodyPtr b2 = get_outboard_link();

  // setup frames
  shared_ptr<const Pose3d> F1(new Pose3d(Quatd::identity()), b1->get_pose()->x);
  shared_ptr<const Pose3d> F2(new Pose3d(Quatd::identity()), b2->get_pose()->x);

  // need velocities in the proper frame
  SVelocityd v1 = Pose3d::transform(b1->velocity().pose, F1, b1->velocity());
  SVelocityd v2 = Pose3d::transform(b2->velocity().pose, F2, b2->velocity());

  // make sure that _u (and by extension _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Revolute joint axis has not been set; set before calling dynamics functions.");

  // mke sure that body is one of the links
  if (b1 != body && b2 != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup constants for calculations
  const Quatd& q1 = b1->get_pose()->q;
  const Quatd& q2 = b2->get_pose()->q;
  const Quatd qd1 = Quatd::deriv(q1, v1.get_angular());
  const Quatd qd2 = Quatd::deriv(q2, v2.get_angular());
  const Vector3d& p1 = b1->get_outer_joint_data(b2).com_to_joint_vec;
  const Vector3d& p2 = b2->get_inner_joint_data(b1).joint_to_com_vec_of;
  const double dqw1 = qd1.w;
  const double dqx1 = qd1.x;
  const double dqy1 = qd1.y;
  const double dqz1 = qd1.z;
  const double dqw2 = qd2.w;
  const double dqx2 = qd2.x;
  const double dqy2 = qd2.y;
  const double dqz2 = qd2.z;
  const double p1x = p1[X];
  const double p1y = p1[Y];
  const double p1z = p1[Z];
  const double p2x = -p2[X];
  const double p2y = -p2[Y];
  const double p2z = -p2[Z];
  const double qw1 = q1.w;
  const double qx1 = q1.x;
  const double qy1 = q1.y;
  const double qz1 = q1.z;
  const double qw2 = q2.w;
  const double qx2 = q2.x;
  const double qy2 = q2.y;
  const double qz2 = q2.z;
  const double uix = _ui[X];
  const double uiy = _ui[Y];
  const double uiz = _ui[Z];
  const double ujx = _uj[X];
  const double ujy = _uj[Y];
  const double ujz = _uj[Z];
  const double v2x = _v2[X];
  const double v2y = _v2[Y];
  const double v2z = _v2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == b1)
  {
    switch (index)
    {
      case 0:
        Cq[0] = (double) 0.0;     
        Cq[1] = (double) 0.0;     
        Cq[2] = (double) 0.0;     
        Cq[3] = 4*p1x*dqw1 + 2*p1z*dqy1 - 2*p1y*dqz1; 
        Cq[4] = 4*p1x*dqx1 + 2*p1y*dqy1 + 2*p1z*dqz1; 
        Cq[5] = 2*p1z*dqw1 + 2*p1y*dqx1; 
        Cq[6] = -2*p1y*dqw1 + 2*p1z*dqx1; 
        break;

      case 1:
        Cq[0] = (double) 0.0;     
        Cq[1] = (double) 0.0;     
        Cq[2] = (double) 0.0;     
        Cq[3] = 4*p1y*dqw1 - 2*p1z*dqx1 + 2*p1x*dqz1; 
        Cq[4] = -2*p1z*dqw1 + 2*p1x*dqy1; 
        Cq[5] = 2*p1x*dqx1 + 4*p1y*dqy1 + 2*p1z*dqz1; 
        Cq[6] = 2*p1x*dqw1 + 2*p1z*dqy1; 
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*p1z*dqw1 + 2*p1y*dqx1 - 2*p1x*dqy1;
        Cq[4] = 2*p1y*dqw1 + 2*p1x*dqz1;
        Cq[5] = -2*p1x*dqw1 + 2*p1y*dqz1;
        Cq[6] = 2*p1x*dqx1 + 2*p1y*dqy1 + 4*p1z*dqz1;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (4*qw1*uix - 2*qz1*uiy + 2*qy1*uiz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qz1*uix + 4*qw1*uiy - 2*qx1*uiz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (-2*qy1*uix + 2*qx1*uiy + 4*qw1*uiz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (4*dqw1*uix - 2*dqz1*uiy + 2*dqy1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqz1*uix + 4*dqw1*uiy - 2*dqx1*uiz)*
    (2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (-2*dqy1*uix + 2*dqx1*uiy + 4*dqw1*uiz)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[4] = (4*qx1*uix + 2*qy1*uiy + 2*qz1*uiz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qy1*uix - 2*qw1*uiz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (2*qz1*uix + 2*qw1*uiy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (4*dqx1*uix + 2*dqy1*uiy + 2*dqz1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqy1*uix - 2*dqw1*uiz)*(2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (2*dqz1*uix + 2*dqw1*uiy)*(2*(-(qw2*qy2) + qx2*qz2)*v2x + 
      2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[5] = (2*qx1*uiy + 2*qw1*uiz)*((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qx1*uix + 4*qy1*uiy + 2*qz1*uiz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (-2*qw1*uix + 2*qz1*uiy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (2*dqx1*uiy + 2*dqw1*uiz)*((-1 + 2*qw2*qw2 + 2*qx2*qx2)*
       v2x + 2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqx1*uix + 4*dqy1*uiy + 2*dqz1*uiz)*
    (2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (-2*dqw1*uix + 2*dqz1*uiy)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[6] = (-2*qw1*uiy + 2*qx1*uiz)*((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qw1*uix + 2*qy1*uiz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (2*qx1*uix + 2*qy1*uiy + 4*qz1*uiz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (-2*dqw1*uiy + 2*dqx1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqw1*uix + 2*dqy1*uiz)*(2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (2*dqx1*uix + 2*dqy1*uiy + 4*dqz1*uiz)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        break;
 
      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (4*qw1*ujx - 2*qz1*ujy + 2*qy1*ujz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qz1*ujx + 4*qw1*ujy - 2*qx1*ujz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (-2*qy1*ujx + 2*qx1*ujy + 4*qw1*ujz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (4*dqw1*ujx - 2*dqz1*ujy + 2*dqy1*ujz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqz1*ujx + 4*dqw1*ujy - 2*dqx1*ujz)*
    (2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (-2*dqy1*ujx + 2*dqx1*ujy + 4*dqw1*ujz)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[4] = (4*qx1*ujx + 2*qy1*ujy + 2*qz1*ujz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qy1*ujx - 2*qw1*ujz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (2*qz1*ujx + 2*qw1*ujy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (4*dqx1*ujx + 2*dqy1*ujy + 2*dqz1*ujz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqy1*ujx - 2*dqw1*ujz)*(2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (2*dqz1*ujx + 2*dqw1*ujy)*(2*(-(qw2*qy2) + qx2*qz2)*v2x + 
      2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[5] = (2*qx1*ujy + 2*qw1*ujz)*((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qx1*ujx + 4*qy1*ujy + 2*qz1*ujz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (-2*qw1*ujx + 2*qz1*ujy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (2*dqx1*ujy + 2*dqw1*ujz)*((-1 + 2*qw2*qw2 + 2*qx2*qx2)*
       v2x + 2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqx1*ujx + 4*dqy1*ujy + 2*dqz1*ujz)*
    (2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (-2*dqw1*ujx + 2*dqz1*ujy)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        Cq[6] = (-2*qw1*ujy + 2*qx1*ujz)*((4*dqw2*qw2 + 4*dqx2*qx2)*v2x + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*v2y + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*v2z) + 
   (2*qw1*ujx + 2*qy1*ujz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*v2x + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*v2y + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2z) + 
   (2*qx1*ujx + 2*qy1*ujy + 4*qz1*ujz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*v2x + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*v2y + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*v2z) + 
   (-2*dqw1*ujy + 2*dqx1*ujz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*v2x + 
      2*(qx2*qy2 - qw2*qz2)*v2y + 2*(qw2*qy2 + qx2*qz2)*v2z) + 
   (2*dqw1*ujx + 2*dqy1*ujz)*(2*(qx2*qy2 + qw2*qz2)*v2x + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*v2y + 
      2*(-(qw2*qx2) + qy2*qz2)*v2z) + 
   (2*dqx1*ujx + 2*dqy1*ujy + 4*dqz1*ujz)*
    (2*(-(qw2*qy2) + qx2*qz2)*v2x + 2*(qw2*qx2 + qy2*qz2)*v2y + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*v2z);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    switch (index)
    {
      case 0:
        Cq[0] = (double) 0.0;     
        Cq[1] = (double) 0.0;      
        Cq[2] = (double) 0.0;      
        Cq[3] = -4*p2x*dqw2 - 2*p2z*dqy2 + 2*p2y*dqz2; 
        Cq[4] = -4*p2x*dqx2 - 2*p2y*dqy2 - 2*p2z*dqz2; 
        Cq[5] = -2*p2z*dqw2 - 2*p2y*dqx2; 
        Cq[6] = 2*p2y*dqw2 - 2*p2z*dqx2; 
        break;

      case 1:
        Cq[0] = (double) 0.0;      
        Cq[1] = (double) 0.0;     
        Cq[2] = (double) 0.0;      
        Cq[3] = -4*p2y*dqw2 + 2*p2z*dqx2 - 2*p2x*dqz2; 
        Cq[4] = 2*p2z*dqw2 - 2*p2x*dqy2; 
        Cq[5] = -2*p2x*dqx2 - 4*p2y*dqy2 - 2*p2z*dqz2; 
        Cq[6] = -2*p2x*dqw2 - 2*p2z*dqy2; 
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -4*p2z*dqw2 - 2*p2y*dqx2 + 2*p2x*dqy2;
        Cq[4] = -2*p2y*dqw2 - 2*p2x*dqz2;
        Cq[5] = 2*p2x*dqw2 - 2*p2y*dqz2;
        Cq[6] = -2*p2x*dqx2 - 2*p2y*dqy2 - 4*p2z*dqz2;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (-2*dqy2*v2x + 2*dqx2*v2y + 4*dqw2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*
    (2*dqz2*v2x + 4*dqw2*v2y - 2*dqx2*v2z) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (4*dqw2*v2x - 2*dqz2*v2y + 2*dqy2*v2z) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*
    (-2*qy2*v2x + 2*qx2*v2y + 4*qw2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qz2*v2x + 4*qw2*v2y - 2*qx2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (4*qw2*v2x - 2*qz2*v2y + 2*qy2*v2z);
        Cq[4] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (2*dqz2*v2x + 2*dqw2*v2y) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*(2*qz2*v2x + 2*qw2*v2y) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*dqy2*v2x - 2*dqw2*v2z) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (4*dqx2*v2x + 2*dqy2*v2y + 2*dqz2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qy2*v2x - 2*qw2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (4*qx2*v2x + 2*qy2*v2y + 2*qz2*v2z);
        Cq[5] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (-2*dqw2*v2x + 2*dqz2*v2y) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*(-2*qw2*v2x + 2*qz2*v2y) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (2*dqx2*v2y + 2*dqw2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*
    (2*dqx2*v2x + 4*dqy2*v2y + 2*dqz2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (2*qx2*v2y + 2*qw2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qx2*v2x + 4*qy2*v2y + 2*qz2*v2z);
        Cq[6] = ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (-2*dqw2*v2y + 2*dqx2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*dqw2*v2x + 2*dqy2*v2z) + 
   (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (2*dqx2*v2x + 2*dqy2*v2y + 4*dqz2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (-2*qw2*v2y + 2*qx2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qw2*v2x + 2*qy2*v2z) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*
    (2*qx2*v2x + 2*qy2*v2y + 4*qz2*v2z);
        break;
 
      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
    (-2*dqy2*v2x + 2*dqx2*v2y + 4*dqw2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz)*
    (2*dqz2*v2x + 4*dqw2*v2y - 2*dqx2*v2z) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
    (4*dqw2*v2x - 2*dqz2*v2y + 2*dqy2*v2z) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz)*
    (-2*qy2*v2x + 2*qx2*v2y + 4*qw2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz)*
    (2*qz2*v2x + 4*qw2*v2y - 2*qx2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz)*
    (4*qw2*v2x - 2*qz2*v2y + 2*qy2*v2z);
        Cq[4] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
    (2*dqz2*v2x + 2*dqw2*v2y) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz)*(2*qz2*v2x + 2*qw2*v2y) + 
   (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz)*(2*dqy2*v2x - 2*dqw2*v2z) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
    (4*dqx2*v2x + 2*dqy2*v2y + 2*dqz2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz)*
    (2*qy2*v2x - 2*qw2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz)*
    (4*qx2*v2x + 2*qy2*v2y + 2*qz2*v2z);
        Cq[5] = (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
    (-2*dqw2*v2x + 2*dqz2*v2y) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz)*(-2*qw2*v2x + 2*qz2*v2y) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
    (2*dqx2*v2y + 2*dqw2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz)*
    (2*dqx2*v2x + 4*dqy2*v2y + 2*dqz2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz)*
    (2*qx2*v2y + 2*qw2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz)*
    (2*qx2*v2x + 4*qy2*v2y + 2*qz2*v2z);
        Cq[6] = ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz)*
    (-2*dqw2*v2y + 2*dqx2*v2z) + 
   (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz)*(2*dqw2*v2x + 2*dqy2*v2z) + 
   (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz)*
    (2*dqx2*v2x + 2*dqy2*v2y + 4*dqz2*v2z) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz)*
    (-2*qw2*v2y + 2*qx2*v2z) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz)*
    (2*qw2*v2x + 2*qy2*v2z) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz)*
    (2*qx2*v2x + 2*qy2*v2y + 4*qz2*v2z);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
*/
}

/// Evaluates the constraint equations
void RevoluteJoint::evaluate_constraints(double C[])
{
/*
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr b1 = get_inboard_link();
  RigidBodyPtr b2 = get_outboard_link();

  // This code was developed using [Shabana, 2003], p. 435-436; variable names
  // have been altered however

  // determine v1, v1i, v1j, and v2 (all in global coordinates)
  Vector3d v1i, v1j;
  Vector3d v1 = Pose3d::transform(get_pose(), GLOBAL, _u);
  Vector3d v2 = Pose3d::transform(b2->get_pose(), GLOBAL, _v2);
  Vector3d::determine_orthonormal_basis(v1, v1i, v1j);

  // determine the global positions of the attachment points and subtract them
  Vector3d r1 = get_location(false);
  Vector3d r2 = get_location(true);
  Vector3d r12 = r1 - r2; 

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = v1i.dot(v2);
  C[4] = v1j.dot(v2); 
*/
}

/// Implements Base::load_from_xml()
void RevoluteJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "RevoluteJoint") == 0);

  // read the global joint axis, if given
  XMLAttrib* axis_attrib = node->get_attrib("axis");
  if (axis_attrib)
  {
    Vector3d axis;
    axis_attrib->get_vector_value(axis);
    set_axis(axis);  
  }

  // compute _q_tare if necessary 
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void RevoluteJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "RevoluteJoint";

  // save the joint axis (global coordinates)
  Vector3d u0 = Pose3d::transform_vector(shared_ptr<const Pose3d>(), _u);
  node->attribs.insert(XMLAttrib("axis", u0));
}

