/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/FixedJoint.h>

using namespace Ravelin;
using namespace Moby;
using std::vector;
using boost::shared_ptr;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
FixedJoint::FixedJoint() : Joint()
{
  // init joint data
  init_data();

  // setup the spatial axis derivative to zero
  _si_deriv.resize(6);
  for (unsigned i=0; i< _si_deriv.size(); i++)
    _si_deriv[i].set_zero();
}

/// Initializes the joint with the specified inboard and outboard links
FixedJoint::FixedJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  // init joint data
  init_data();

  // setup the spatial axis and its derivative to zero
  _si_deriv.resize(6);
  for (unsigned i=0; i< _si_deriv.size(); i++)
    _si_deriv[i].set_zero();
}  

/// Sets spatial axes to zero 
void FixedJoint::update_spatial_axes()
{
  // call parent method
  Joint::update_spatial_axes();

  // setup si
  for (unsigned i=0; i< _si.size(); i++)
    _si[i].set_zero();

  // setup complement of si
  _s_bar.resize(6);
  for (unsigned i=0; i< _s_bar.size(); i++)
  {
    _s_bar[i].set_zero();
    _s_bar[i][i] = (double) 1.0;
  }
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
  shared_ptr<const Pose3d> Ti = inboard->get_transform();
  shared_ptr<const Pose3d> To = outboard->get_transform();

  // get the rotation matrices
  Matrix3d Ri = Ti->q;
  Matrix3d Ro = To->q;

  // compute the relative transform
  Matrix3d Rrel = Matrix3d::transpose(Ro) * Ri;
  _T->q = Rrel;
  _T->x.set_zero();

  // compute the vector from the inner link to the outer link in inner link
  // frame
  _ui = Ti->inverse_transform(Point3d(To->x));

  // compute the constant orientation term
  _rconst[X] = Vector3d::dot(Ri.get_row(X), Ro.get_row(X));
  _rconst[Y] = Vector3d::dot(Ri.get_row(Y), Ro.get_row(Y));
  _rconst[Z] = Vector3d::dot(Ri.get_row(Z), Ro.get_row(Z));
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
shared_ptr<const Pose3d> FixedJoint::get_transform()
{
  // get the link transforms
  return _T;
}

/// Gets the derivative for the spatial axes for this joint
const vector<Twistd>& FixedJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  return _si_deriv;
}

/// Computes the constraint Jacobian with respect to a body
void FixedJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned idx, double Cq[7])
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
  const Quatd& q1 = inner->get_transform()->q;
  const Quatd& q2 = outer->get_transform()->q;
  const double qw1 = q1.w;
  const double qx1 = q1.x;
  const double qy1 = q1.y;
  const double qz1 = q1.z;
  const double qw2 = q2.w;
  const double qx2 = q2.x;
  const double qy2 = q2.y;
  const double qz2 = q2.z;
  const double p1x = _ui[X];
  const double p1y = _ui[Y];
  const double p1z = _ui[Z];
  const double p2x = (double) 0.0;
  const double p2y = (double) 0.0;
  const double p2z = (double) 0.0;

  // setup the constraint equations (derived from Mathematica)
  if (body == inner)
  {
    switch (idx)
    {
      case 0:
        Cq[0] = (double) 1.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1;
        Cq[4] = 4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1;
        Cq[5] = 2*p1z*qw1 + 2*p1y*qx1;
        Cq[6] = -2*p1y*qw1 + 2*p1z*qx1;
        break;

      case 1:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 1.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1;
        Cq[4] = -2*p1z*qw1 + 2*p1x*qy1;
        Cq[5] = 2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1;
        Cq[6] = 2*p1x*qw1 + 2*p1z*qy1;
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 1.0;
        Cq[3] = 4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1;
        Cq[4] = 2*p1y*qw1 + 2*p1x*qz1;
        Cq[5] = -2*p1x*qw1 + 2*p1y*qz1;
        Cq[6] = 2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*qw1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) - 
    4*qz1*(qx2*qy2 - qw2*qz2) + 4*qy1*(qw2*qy2 + qx2*qz2);
        Cq[4] = 4*qx1*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
    4*qy1*(qx2*qy2 - qw2*qz2) + 4*qz1*(qw2*qy2 + qx2*qz2);
        Cq[5] = 4*qx1*(qx2*qy2 - qw2*qz2) + 4*qw1*(qw2*qy2 + qx2*qz2);
        Cq[6] = -4*qw1*(qx2*qy2 - qw2*qz2) + 4*qx1*(qw2*qy2 + qx2*qz2);
        break;

      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*qw1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
    4*qz1*(qx2*qy2 + qw2*qz2) - 4*qx1*(-(qw2*qx2) + qy2*qz2);
        Cq[4] = 4*qy1*(qx2*qy2 + qw2*qz2) - 4*qw1*(-(qw2*qx2) + qy2*qz2);
        Cq[5] = 4*qy1*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
    4*qx1*(qx2*qy2 + qw2*qz2) + 4*qz1*(-(qw2*qx2) + qy2*qz2);
        Cq[6] = 4*qw1*(qx2*qy2 + qw2*qz2) + 4*qy1*(-(qw2*qx2) + qy2*qz2);
        break;

      case 5:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
        Cq[0] = (double) -1.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2;
        Cq[4] = -4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2;
        Cq[5] = -2*p2z*qw2 - 2*p2y*qx2;
        Cq[6] = 2*p2y*qw2 - 2*p2z*qx2;
        break;

      case 1:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) -1.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2;
        Cq[4] = 2*p2z*qw2 - 2*p2x*qy2;
        Cq[5] = -2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2;
        Cq[6] = -2*p2x*qw2 - 2*p2z*qy2;
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) -1.0;
        Cq[3] = -4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2;
        Cq[4] = -2*p2y*qw2 - 2*p2x*qz2;
        Cq[5] = 2*p2x*qw2 - 2*p2y*qz2;
        Cq[6] = -2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*qw2*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
    4*qy2*(qw1*qy1 + qx1*qz1) - 4*(qx1*qy1 - qw1*qz1)*qz2;
        Cq[4] = 4*(-1 + 2*qw1*qw1 + 2*qx1*qx1)*qx2 + 
    4*qy2*(qx1*qy1 - qw1*qz1) + 4*(qw1*qy1 + qx1*qz1)*qz2;
        Cq[5] = 4*qx2*(qx1*qy1 - qw1*qz1) + 4*qw2*(qw1*qy1 + qx1*qz1);
        Cq[6] = -4*qw2*(qx1*qy1 - qw1*qz1) + 4*qx2*(qw1*qy1 + qx1*qz1);
        break;

      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*qw2*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
    4*qx2*(-(qw1*qx1) + qy1*qz1) + 4*(qx1*qy1 + qw1*qz1)*qz2;
        Cq[4] = 4*qy2*(qx1*qy1 + qw1*qz1) - 4*qw2*(-(qw1*qx1) + qy1*qz1);
        Cq[5] = 4*(-1 + 2*(qw1*qw1 + qy1*qy1))*qy2 + 
    4*qx2*(qx1*qy1 + qw1*qz1) + 4*(-(qw1*qx1) + qy1*qz1)*qz2;
        Cq[6] = 4*qw2*(qx1*qy1 + qw1*qz1) + 4*qy2*(-(qw1*qx1) + qy1*qz1);
        break;

      case 5:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
void FixedJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned idx, double Cq[7])
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
  const Quatd& q1 = inner->get_transform()->q;
  const Quatd& q2 = outer->get_transform()->q;
  const double qw1 = q1.w;
  const double qx1 = q1.x;
  const double qy1 = q1.y;
  const double qz1 = q1.z;
  const double qw2 = q2.w;
  const double qx2 = q2.x;
  const double qy2 = q2.y;
  const double qz2 = q2.z;
  const double p1x = _ui[X];
  const double p1y = _ui[Y];
  const double p1z = _ui[Z];
  const double p2x = (double) 0.0;
  const double p2y = (double) 0.0;
  const double p2z = (double) 0.0;
  Quatd q1d = Quatd::deriv(q1, inner->get_velocity().get_linear());
  Quatd q2d = Quatd::deriv(q2, outer->get_velocity().get_angular());
  const double dqw1 = q1d.w;
  const double dqx1 = q1d.x;
  const double dqy1 = q1d.y;
  const double dqz1 = q1d.z;
  const double dqw2 = q2d.w;
  const double dqx2 = q2d.x;
  const double dqy2 = q2d.y;
  const double dqz2 = q2d.z;

  // setup the constraint equations (derived from Mathematica)
  if (body == inner)
  {
    switch (idx)
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
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -4*dqw2*p2x + 2*dqz2*p2y - 2*dqy2*p2z;
        Cq[4] = -4*dqx2*p2x - 2*dqy2*p2y - 2*dqz2*p2z;
        Cq[5] = -2*dqx2*p2y - 2*dqw2*p2z;
        Cq[6] = 2*dqw2*p2y - 2*dqx2*p2z;
        break;

      case 1:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -2*dqz2*p2x - 4*dqw2*p2y + 2*dqx2*p2z;
        Cq[4] = -2*dqy2*p2x + 2*dqw2*p2z;
        Cq[5] = -2*dqx2*p2x - 4*dqy2*p2y - 2*dqz2*p2z;
        Cq[6] = -2*dqw2*p2x - 2*dqy2*p2z;
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 2*dqy2*p2x - 2*dqx2*p2y - 4*dqw2*p2z;
        Cq[4] = -2*dqz2*p2x - 2*dqw2*p2y;
        Cq[5] = 2*dqw2*p2x - 2*dqz2*p2y;
        Cq[6] = -2*dqx2*p2x - 2*dqy2*p2y - 4*dqz2*p2z;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
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
void FixedJoint::evaluate_constraints(double C[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // get the transforms and orientations for the two links
  shared_ptr<const Pose3d> Ti = inner->get_transform();
  shared_ptr<const Pose3d> To = outer->get_transform();
  Matrix3d Ri = Ti->q;
  Matrix3d Ro = To->q;

  // evaluate the relative position
  Point3d rpos = Ti->transform(_ui) + Ti->x - To->x; 

  // setup C
  C[0] = rpos[0];
  C[1] = rpos[1];
  C[2] = rpos[2];
  C[3] = Vector3d::dot(Ri.get_row(X), Ro.get_row(X)) - _rconst[X];
  C[4] = Vector3d::dot(Ri.get_row(Y), Ro.get_row(Y)) - _rconst[Y];
  C[5] = Vector3d::dot(Ri.get_row(Z), Ro.get_row(Z)) - _rconst[Z];
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

