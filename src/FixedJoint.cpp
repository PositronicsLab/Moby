/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/FixedJoint.h>

using namespace Moby;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
FixedJoint::FixedJoint() : Joint()
{
  // init joint data
  init_data();

  // setup the spatial axis derivative to zero
  _si_deriv = SMatrix6N::zero(6,0);
}

/// Initializes the joint with the specified inboard and outboard links
FixedJoint::FixedJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  // init joint data
  init_data();

  // setup the spatial axis and its derivative to zero
  _si_deriv = SMatrix6N::zero(6,0);
}  

/// Sets spatial axes to zero 
void FixedJoint::update_spatial_axes()
{
  // call parent method
  Joint::update_spatial_axes();

  // setup si
  _si.set_zero();

  // setup complement of si
  _s_bar.set_zero();
  for (unsigned i=0; i< 6; i++)
    _s_bar(i,i) = (Real) 1.0;
}

/// Setup the transform from one joint to another
void FixedJoint::setup_joint()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  assert(inboard);
  assert(outboard);

  // get the transforms
  const Matrix4& Ti = inboard->get_transform();
  const Matrix4& To = outboard->get_transform();

  // get the rotation matrices
  Matrix3 Ri = Ti.get_rotation();
  Matrix3 Ro = To.get_rotation();

  // compute the relative transform
  Matrix3 Rrel = Matrix3::transpose(Ro) * Ri;
  _T.set_rotation(&Rrel);
  _T.set_translation(ZEROS_3);

  // compute the vector from the inner link to the outer link in inner link
  // frame
  _ui = Ti.inverse_mult_point(To.get_translation());

  // compute the constant orientation term
  _rconst[X] = Vector3::dot(Ri.get_row(X), Ro.get_row(X));
  _rconst[Y] = Vector3::dot(Ri.get_row(Y), Ro.get_row(Y));
  _rconst[Z] = Vector3::dot(Ri.get_row(Z), Ro.get_row(Z));
}

/// Sets the inboard link
void FixedJoint::set_inboard_link(RigidBodyPtr link)
{
  // call parent method since it does all of the work
  Joint::set_inboard_link(link);

  // see whether there is already an outboard link
  RigidBodyPtr outboard = get_outboard_link();

  // if there is an outboard link set, compute relative, fixed generalized
  // coordinates
  if (outboard)
    setup_joint();
}

/// Sets the outboard link
void FixedJoint::set_outboard_link(RigidBodyPtr link)
{
  // call parent method since it does all of the work
  Joint::set_outboard_link(link);

  // see whether there is already an inboard link
  RigidBodyPtr inboard = get_inboard_link();

  // if there is an inboard link set, compute relative, fixed generalized
  // coordinates
  if (inboard)
    setup_joint();
}

/// Gets the (local) transform for this joint (constant)
const Matrix4& FixedJoint::get_transform()
{
  // get the link transforms
  return _T;
}

/// Gets the derivative for the spatial axes for this joint
const SMatrix6N& FixedJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  return _si_deriv;
}

/// Computes the constraint Jacobian with respect to a body
void FixedJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned idx, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (idx >= num_constraint_eqns())
    throw std::runtime_error("Invalid constraint index specified!");

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
    return;

  // setup necessary variables
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Real qw1 = q1.w;
  const Real qx1 = q1.x;
  const Real qy1 = q1.y;
  const Real qz1 = q1.z;
  const Real qw2 = q2.w;
  const Real qx2 = q2.x;
  const Real qy2 = q2.y;
  const Real qz2 = q2.z;
  const Real p1x = _ui[X];
  const Real p1y = _ui[Y];
  const Real p1z = _ui[Z];
  const Real p2x = (Real) 0.0;
  const Real p2y = (Real) 0.0;
  const Real p2z = (Real) 0.0;

  // setup the constraint equations (derived from Mathematica)
  if (body == inner)
  {
    switch (idx)
    {
      case 0:
        Cq[0] = (Real) 1.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1;
        Cq[4] = 4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1;
        Cq[5] = 2*p1z*qw1 + 2*p1y*qx1;
        Cq[6] = -2*p1y*qw1 + 2*p1z*qx1;
        break;

      case 1:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 1.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1;
        Cq[4] = -2*p1z*qw1 + 2*p1x*qy1;
        Cq[5] = 2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1;
        Cq[6] = 2*p1x*qw1 + 2*p1z*qy1;
        break;

      case 2:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 1.0;
        Cq[3] = 4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1;
        Cq[4] = 2*p1y*qw1 + 2*p1x*qz1;
        Cq[5] = -2*p1x*qw1 + 2*p1y*qz1;
        Cq[6] = 2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1;
        break;

      case 3:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) - 
    4*qz1*(qx2*qy2 - qw2*qz2) + 4*qy1*(qw2*qy2 + qx2*qz2);
        Cq[4] = 4*qx1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
    4*qy1*(qx2*qy2 - qw2*qz2) + 4*qz1*(qw2*qy2 + qx2*qz2);
        Cq[5] = 4*qx1*(qx2*qy2 - qw2*qz2) + 4*qw1*(qw2*qy2 + qx2*qz2);
        Cq[6] = -4*qw1*(qx2*qy2 - qw2*qz2) + 4*qx1*(qw2*qy2 + qx2*qz2);
        break;

      case 4:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
    4*qz1*(qx2*qy2 + qw2*qz2) - 4*qx1*(-(qw2*qx2) + qy2*qz2);
        Cq[4] = 4*qy1*(qx2*qy2 + qw2*qz2) - 4*qw1*(-(qw2*qx2) + qy2*qz2);
        Cq[5] = 4*qy1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
    4*qx1*(qx2*qy2 + qw2*qz2) + 4*qz1*(-(qw2*qx2) + qy2*qz2);
        Cq[6] = 4*qw1*(qx2*qy2 + qw2*qz2) + 4*qy1*(-(qw2*qx2) + qy2*qz2);
        break;

      case 5:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*qy1*(-(qw2*qy2) + qx2*qz2) + 
    4*qx1*(qw2*qx2 + qy2*qz2) + 
    4*qw1*(-1 + 2*(qw2*qw2 + qz2*qz2));
        Cq[4] = 4*qz1*(-(qw2*qy2) + qx2*qz2) + 4*qw1*(qw2*qx2 + qy2*qz2);
        Cq[5] = -4*qw1*(-(qw2*qy2) + qx2*qz2) + 4*qz1*(qw2*qx2 + qy2*qz2);
        Cq[6] = 4*qx1*(-(qw2*qy2) + qx2*qz2) + 4*qy1*(qw2*qx2 + qy2*qz2) + 
    4*qz1*(-1 + 2*(qw2*qw2 + qz2*qz2));
        break;
    }
  }
  else
  {
    switch (idx)
    {
      case 0:
        Cq[0] = (Real) -1.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2;
        Cq[4] = -4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2;
        Cq[5] = -2*p2z*qw2 - 2*p2y*qx2;
        Cq[6] = 2*p2y*qw2 - 2*p2z*qx2;
        break;

      case 1:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) -1.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2;
        Cq[4] = 2*p2z*qw2 - 2*p2x*qy2;
        Cq[5] = -2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2;
        Cq[6] = -2*p2x*qw2 - 2*p2z*qy2;
        break;

      case 2:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) -1.0;
        Cq[3] = -4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2;
        Cq[4] = -2*p2y*qw2 - 2*p2x*qz2;
        Cq[5] = 2*p2x*qw2 - 2*p2y*qz2;
        Cq[6] = -2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2;
        break;

      case 3:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw2*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
    4*qy2*(qw1*qy1 + qx1*qz1) - 4*(qx1*qy1 - qw1*qz1)*qz2;
        Cq[4] = 4*(-1 + 2*qw1*qw1 + 2*qx1*qx1)*qx2 + 
    4*qy2*(qx1*qy1 - qw1*qz1) + 4*(qw1*qy1 + qx1*qz1)*qz2;
        Cq[5] = 4*qx2*(qx1*qy1 - qw1*qz1) + 4*qw2*(qw1*qy1 + qx1*qz1);
        Cq[6] = -4*qw2*(qx1*qy1 - qw1*qz1) + 4*qx2*(qw1*qy1 + qx1*qz1);
        break;

      case 4:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw2*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
    4*qx2*(-(qw1*qx1) + qy1*qz1) + 4*(qx1*qy1 + qw1*qz1)*qz2;
        Cq[4] = 4*qy2*(qx1*qy1 + qw1*qz1) - 4*qw2*(-(qw1*qx1) + qy1*qz1);
        Cq[5] = 4*(-1 + 2*(qw1*qw1 + qy1*qy1))*qy2 + 
    4*qx2*(qx1*qy1 + qw1*qz1) + 4*(-(qw1*qx1) + qy1*qz1)*qz2;
        Cq[6] = 4*qw2*(qx1*qy1 + qw1*qz1) + 4*qy2*(-(qw1*qx1) + qy1*qz1);
        break;

      case 5:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*qy2*(-(qw1*qy1) + qx1*qz1) + 4*qx2*(qw1*qx1 + qy1*qz1) + 
    4*qw2*(-1 + 2*(qw1*qw1 + qz1*qz1));
        Cq[4] = 4*qw2*(qw1*qx1 + qy1*qz1) + 4*(-(qw1*qy1) + qx1*qz1)*qz2;
        Cq[5] = -4*qw2*(-(qw1*qy1) + qx1*qz1) + 4*(qw1*qx1 + qy1*qz1)*qz2;
        Cq[6] = 4*qx2*(-(qw1*qy1) + qx1*qz1) + 4*qy2*(qw1*qx1 + qy1*qz1) + 
    4*(-1 + 2*(qw1*qw1 + qz1*qz1))*qz2;
        break;
    }
  }
}

/// Computes the time derivative of the constraint Jacobian with respect to a body
void FixedJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned idx, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (idx >= num_constraint_eqns())
    throw std::runtime_error("Invalid constraint index specified!");

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
    return;

  // setup necessary variables
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Real qw1 = q1.w;
  const Real qx1 = q1.x;
  const Real qy1 = q1.y;
  const Real qz1 = q1.z;
  const Real qw2 = q2.w;
  const Real qx2 = q2.x;
  const Real qy2 = q2.y;
  const Real qz2 = q2.z;
  const Real p1x = _ui[X];
  const Real p1y = _ui[Y];
  const Real p1z = _ui[Z];
  const Real p2x = (Real) 0.0;
  const Real p2y = (Real) 0.0;
  const Real p2z = (Real) 0.0;
  Quat q1d = Quat::deriv(q1, inner->get_avel());
  Quat q2d = Quat::deriv(q2, outer->get_avel());
  const Real dqw1 = q1d.w;
  const Real dqx1 = q1d.x;
  const Real dqy1 = q1d.y;
  const Real dqz1 = q1d.z;
  const Real dqw2 = q2d.w;
  const Real dqx2 = q2d.x;
  const Real dqy2 = q2d.y;
  const Real dqz2 = q2d.z;

  // setup the constraint equations (derived from Mathematica)
  if (body == inner)
  {
    switch (idx)
    {
      case 0:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*p1x*dqw1 + 2*p1z*dqy1 - 2*p1y*dqz1;
        Cq[4] = 4*p1x*dqx1 + 2*p1y*dqy1 + 2*p1z*dqz1;
        Cq[5] = 2*p1z*dqw1 + 2*p1y*dqx1;
        Cq[6] = -2*p1y*dqw1 + 2*p1z*dqx1;
        break;

      case 1:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*p1y*dqw1 - 2*p1z*dqx1 + 2*p1x*dqz1;
        Cq[4] = -2*p1z*dqw1 + 2*p1x*dqy1;
        Cq[5] = 2*p1x*dqx1 + 4*p1y*dqy1 + 2*p1z*dqz1;
        Cq[6] = 2*p1x*dqw1 + 2*p1z*dqy1;
        break;

      case 2:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*p1z*dqw1 + 2*p1y*dqx1 - 2*p1x*dqy1;
        Cq[4] = 2*p1y*dqw1 + 2*p1x*dqz1;
        Cq[5] = -2*p1x*dqw1 + 2*p1y*dqz1;
        Cq[6] = 2*p1x*dqx1 + 2*p1y*dqy1 + 4*p1z*dqz1;
        break;

      case 3:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw1*(4*dqw2*qw2 + 4*dqx2*qx2) + 
   4*dqw1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) - 
   4*qz1*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
   4*qy1*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2) - 
   4*dqz1*(qx2*qy2 - qw2*qz2) + 4*dqy1*(qw2*qy2 + qx2*qz2);
        Cq[4] = 4*qx1*(4*dqw2*qw2 + 4*dqx2*qx2) + 
   4*dqx1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
   4*qy1*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
   4*qz1*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2) + 
   4*dqy1*(qx2*qy2 - qw2*qz2) + 4*dqz1*(qw2*qy2 + qx2*qz2);
        Cq[5] = 4*qx1*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
   4*qw1*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2) + 
   4*dqx1*(qx2*qy2 - qw2*qz2) + 4*dqw1*(qw2*qy2 + qx2*qz2);
        Cq[6] = -4*qw1*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
   4*qx1*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2) - 
   4*dqw1*(qx2*qy2 - qw2*qz2) + 4*dqx1*(qw2*qy2 + qx2*qz2);
        break;

      case 4:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 8*qw1*(2*dqw2*qw2 + 2*dqy2*qy2) + 
   4*dqw1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
   4*qz1*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
   4*qx1*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   4*dqz1*(qx2*qy2 + qw2*qz2) - 4*dqx1*(-(qw2*qx2) + qy2*qz2);
        Cq[4] = 4*qy1*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
   4*qw1*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   4*dqy1*(qx2*qy2 + qw2*qz2) - 4*dqw1*(-(qw2*qx2) + qy2*qz2);
        Cq[5] = 8*qy1*(2*dqw2*qw2 + 2*dqy2*qy2) + 
   4*dqy1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
   4*qx1*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
   4*qz1*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   4*dqx1*(qx2*qy2 + qw2*qz2) + 4*dqz1*(-(qw2*qx2) + qy2*qz2);
        Cq[6] = 4*qw1*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
   4*qy1*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   4*dqw1*(qx2*qy2 + qw2*qz2) + 4*dqy1*(-(qw2*qx2) + qy2*qz2);
        break;

      case 5:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*qy1*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
   4*qx1*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   8*qw1*(2*dqw2*qw2 + 2*dqz2*qz2) - 4*dqy1*(-(qw2*qy2) + qx2*qz2) + 
   4*dqx1*(qw2*qx2 + qy2*qz2) + 
   4*dqw1*(-1 + 2*(qw2*qw2 + qz2*qz2));
        Cq[4] = 4*qz1*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
   4*qw1*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   4*dqz1*(-(qw2*qy2) + qx2*qz2) + 4*dqw1*(qw2*qx2 + qy2*qz2);
        Cq[5] = -4*qw1*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
   4*qz1*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
   4*dqw1*(-(qw2*qy2) + qx2*qz2) + 4*dqz1*(qw2*qx2 + qy2*qz2);
        Cq[6] = 4*qx1*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
   4*qy1*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
   8*qz1*(2*dqw2*qw2 + 2*dqz2*qz2) + 4*dqx1*(-(qw2*qy2) + qx2*qz2) + 
   4*dqy1*(qw2*qx2 + qy2*qz2) + 
   4*dqz1*(-1 + 2*(qw2*qw2 + qz2*qz2)); 
        break;
    }
  }
  else
  {
    switch (idx)
    {
      case 0:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*dqw2*p2x + 2*dqz2*p2y - 2*dqy2*p2z;
        Cq[4] = -4*dqx2*p2x - 2*dqy2*p2y - 2*dqz2*p2z;
        Cq[5] = -2*dqx2*p2y - 2*dqw2*p2z;
        Cq[6] = 2*dqw2*p2y - 2*dqx2*p2z;
        break;

      case 1:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -2*dqz2*p2x - 4*dqw2*p2y + 2*dqx2*p2z;
        Cq[4] = -2*dqy2*p2x + 2*dqw2*p2z;
        Cq[5] = -2*dqx2*p2x - 4*dqy2*p2y - 2*dqz2*p2z;
        Cq[6] = -2*dqw2*p2x - 2*dqy2*p2z;
        break;

      case 2:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 2*dqy2*p2x - 2*dqx2*p2y - 4*dqw2*p2z;
        Cq[4] = -2*dqz2*p2x - 2*dqw2*p2y;
        Cq[5] = 2*dqw2*p2x - 2*dqz2*p2y;
        Cq[6] = -2*dqx2*p2x - 2*dqy2*p2y - 4*dqz2*p2z;
        break;

      case 3:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 4*qw2*(4*dqw1*qw1 + 4*dqx1*qx1) + 
   4*dqw2*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
   4*qy2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
   4*dqz2*(qx1*qy1 - qw1*qz1) + 4*dqy2*(qw1*qy1 + qx1*qz1) - 
   4*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*qz2;
        Cq[4] = 4*dqx2*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
   4*(4*dqw1*qw1 + 4*dqx1*qx1)*qx2 + 
   4*qy2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
   4*dqy2*(qx1*qy1 - qw1*qz1) + 4*dqz2*(qw1*qy1 + qx1*qz1) + 
   4*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*qz2;
        Cq[5] = 4*qx2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
   4*qw2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) + 
   4*dqx2*(qx1*qy1 - qw1*qz1) + 4*dqw2*(qw1*qy1 + qx1*qz1);
        Cq[6] = -4*qw2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
   4*qx2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
   4*dqw2*(qx1*qy1 - qw1*qz1) + 4*dqx2*(qw1*qy1 + qx1*qz1);
        break;

      case 4:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = 8*qw2*(2*dqw1*qw1 + 2*dqy1*qy1) + 
   4*dqw2*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
   4*qx2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   4*dqz2*(qx1*qy1 + qw1*qz1) - 4*dqx2*(-(qw1*qx1) + qy1*qz1) + 
   4*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*qz2;
        Cq[4] = 4*qy2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) - 
   4*qw2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   4*dqy2*(qx1*qy1 + qw1*qz1) - 4*dqw2*(-(qw1*qx1) + qy1*qz1);
        Cq[5] = 4*dqy2*(-1 + 2*(qw1*qw1 + qy1*qy1)) + 
   8*(2*dqw1*qw1 + 2*dqy1*qy1)*qy2 + 
   4*qx2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
   4*dqx2*(qx1*qy1 + qw1*qz1) + 4*dqz2*(-(qw1*qx1) + qy1*qz1) + 
   4*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*qz2;
        Cq[6] = 4*qw2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
   4*qy2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   4*dqw2*(qx1*qy1 + qw1*qz1) + 4*dqy2*(-(qw1*qx1) + qy1*qz1);
        break;

      case 5:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*qy2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
   4*qx2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   8*qw2*(2*dqw1*qw1 + 2*dqz1*qz1) - 4*dqy2*(-(qw1*qy1) + qx1*qz1) + 
   4*dqx2*(qw1*qx1 + qy1*qz1) + 
   4*dqw2*(-1 + 2*(qw1*qw1 + qz1*qz1));
        Cq[4] = 4*qw2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   4*dqz2*(-(qw1*qy1) + qx1*qz1) + 4*dqw2*(qw1*qx1 + qy1*qz1) + 
   4*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*qz2;
        Cq[5] = -4*qw2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) - 
   4*dqw2*(-(qw1*qy1) + qx1*qz1) + 4*dqz2*(qw1*qx1 + qy1*qz1) + 
   4*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*qz2;
        Cq[6] = 4*qx2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
   4*qy2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
   4*dqx2*(-(qw1*qy1) + qx1*qz1) + 4*dqy2*(qw1*qx1 + qy1*qz1) + 
   4*dqz2*(-1 + 2*(qw1*qw1 + qz1*qz1)) + 
   8*(2*dqw1*qw1 + 2*dqz1*qz1)*qz2;
        break;
    }
  }
}

/// Evaluates the constraint equations
void FixedJoint::evaluate_constraints(Real C[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // get the transforms and orientations for the two links
  const Matrix4& Ti = inner->get_transform();
  const Matrix4& To = outer->get_transform();
  Matrix3 Ri = Ti.get_rotation();
  Matrix3 Ro = To.get_rotation();

  // evaluate the relative position
  Vector3 rpos = Ti.mult_point(_ui) - To.get_translation(); 

  // setup C
  C[0] = rpos[0];
  C[1] = rpos[1];
  C[2] = rpos[2];
  C[3] = Vector3::dot(Ri.get_row(X), Ro.get_row(X)) - _rconst[X];
  C[4] = Vector3::dot(Ri.get_row(Y), Ro.get_row(Y)) - _rconst[Y];
  C[5] = Vector3::dot(Ri.get_row(Z), Ro.get_row(Z)) - _rconst[Z];
}

/// Implements Base::load_from_xml()
void FixedJoint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "FixedJoint") == 0);

  // set the relative joint transformation
  determine_q_tare();
}

/// Implements Base::save_to_xml()
void FixedJoint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "FixedJoint";
}

