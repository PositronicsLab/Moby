/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/XMLTree.h>
#include <Moby/RigidBody.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/PrismaticJoint.h>

using namespace Moby;

/// Initializes the joint
/**
 * The axis of rotation is set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
PrismaticJoint::PrismaticJoint() : Joint()
{
  // init the joint data
  init_data();

  // set the axes and associated vectors to zeros initially
  _u = ZEROS_3;
  _v2 = ZEROS_3;
  _ui = ZEROS_3;
  _uj = ZEROS_3;

  // set the rotation for the transform to identity
  const Matrix3 IDENTITY_3x3 = Matrix3::identity();
  _T.set_rotation(&IDENTITY_3x3);

  // setup the spatial axis derivative to zero
  _si_deriv.clear();
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
PrismaticJoint::PrismaticJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  // init the joint data
  init_data();

  // set the axis to zeros initially
  _u = ZEROS_3;
  _v2 = ZEROS_3;
  _ui = ZEROS_3;
  _uj = ZEROS_3;

  // set the rotation for the transform to identity
  const Matrix3 IDENTITY_3x3 = Matrix3::identity();
  _T.set_rotation(&IDENTITY_3x3);

  // setup the spatial axis derivative to zero
  _si_deriv.clear();
}  

/// Gets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
Vector3 PrismaticJoint::get_axis_global() const
{
  // make sure that the inboard link has been set
  RigidBodyPtr inboard_link(get_inboard_link());
  if (!inboard_link)
  {
    std::cerr << "PrismaticJoint::get_axis_global() - attempt to get axis w/o inboard link!" << std::endl;
    return ZEROS_3;
  }

  // transform into global coordinates and return
  return inboard_link->get_transform().mult_vector(_u);
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
void PrismaticJoint::set_axis_local(const Vector3& axis) 
{
  // check that axis is ok 
  if (std::fabs(axis.norm() - (double) 1.0) < NEAR_ZERO)
    throw UndefinedAxisException(); 
 
  // normalize the axis, in case caller did not 
  Vector3 naxis = Vector3::normalize(axis); 

  // get the inboard and outboard links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that the links have been set
  if (!inner || !outer)
    throw std::runtime_error("Unable to set joint axis without links being set!");

  // set the joint axis in the inner link frame
  _u = naxis; 
  update_spatial_axes(); 

  // set the joint axis in the outer link frame and setup associated
  // vectors needed for maximal coordinate articulated bodies
  _v2 = inner->get_transform().mult_vector(naxis);
  _v2 = outer->get_transform().transpose_mult_vector(_v2);
  Vector3::determine_orthonormal_basis(_u, _ui, _uj);
}        

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void PrismaticJoint::set_axis_global(const Vector3& axis)
{
  // check that axis is unit vector 
  if (std::fabs(axis.norm() - (double) 1.0) < NEAR_ZERO)
    throw UndefinedAxisException(); 

   // normalize the axis, in case caller did not 
  Vector3 naxis = Vector3::normalize(axis); 

  // get the inboard and outboard links
  RigidBodyPtr inboard_link = get_inboard_link();
  RigidBodyPtr outboard_link = get_outboard_link();

  // transform the axis into local coordinates and update spatial axis
  _u = inboard_link->get_transform().transpose_mult_vector(naxis);
  update_spatial_axes();

  // set the joint axis in the outer link frame and setup associated
  // vectors needed for maximal coordinate articulated bodies
  _v2 = outboard_link->get_transform().transpose_mult_vector(naxis);
  Vector3::determine_orthonormal_basis(_u, _ui, _uj);
}

/// Updates the spatial axis for this joint
void PrismaticJoint::update_spatial_axes()
{
  // call parent method
  Joint::update_spatial_axes();

  // get the inboard and outboard links
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard || !outboard)
    return;

  // if the axis is not normal, return
  if (std::fabs(_u.norm_sq() - (double) 1.0) > NEAR_ZERO)
    return;

  // get the axis in outer link coordinates
  Vector3 u0 = inboard->get_transform().mult_vector(_u);
  Vector3 ui = outboard->get_transform().transpose_mult_vector(u0);

  // update the spatial axis in link coordinates
  SVector6 si_vec;
  si_vec.set_upper(ZEROS_3);
  si_vec.set_lower(ui);
  _si.set_column(0, si_vec);

  // update the complement of the spatial axis in link coordinates
  calc_s_bar_from_si();
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void PrismaticJoint::determine_q(VectorN& q)
{
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();

  // verify that the inboard and outboard links are set
  if (!inboard || !outboard)
    throw std::runtime_error("determine_q() called on NULL inboard and/or outboard links!");

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - 1.0) > NEAR_ZERO)
    throw UndefinedAxisException();

  // get the attachment points on the link (global coords)
  Vector3 p1 = get_position_global(false);
  Vector3 p2 = get_position_global(true);

  // get the joint axis in the global frame
  Vector3 ug = inboard->get_transform().mult_vector(_u);

  // now, we'll project p2 onto the axis ug; points will be setup so that
  // ug passes through origin on inboard
  q.resize(num_dof());
  q[DOF_1] = ug.dot(p2-p1);
}

/// Gets the (local) transform for this joint
const Matrix4& PrismaticJoint::get_transform()
{
  // note that rotation is set to identity in the constructors
  _T.set_translation(_u * (this->q[DOF_1] + this->_q_tare[DOF_1]));

  return _T;
}

/// Gets the derivative fo the spatial axes for this joint
vector<Twistd>& PrismaticJoint::get_spatial_axes_dot()
{
  return _si_deriv;
}

/// Calculates the constraint Jacobian
void PrismaticJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup constants for calculations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
  const double x1 = inner->get_position()[X];
  const double y1 = inner->get_position()[Y];
  const double z1 = inner->get_position()[Z];
  const double x2 = outer->get_position()[X];
  const double y2 = outer->get_position()[Y];
  const double z2 = outer->get_position()[Z];
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

  // setup the constraint equations (from Shabana, p. 437), eq. 7.179
  if (body == inner)
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
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

      case 1:
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

      case 2:
        Cq[0] = (-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
    2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz;
        Cq[1] = 2*(qx1*qy1 + qw1*qz1)*uix + 
    (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
    2*(-(qw1*qx1) + qy1*qz1)*uiz;
        Cq[2] = 2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
    (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz;
        Cq[3] = (4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
    (4*qw1*uix - 2*qz1*uiy + 2*qy1*uiz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qz1*uix + 4*qw1*uiy - 2*qx1*uiz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (-2*qy1*uix + 2*qx1*uiy + 4*qw1*uiz)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[4] = (4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (-2*p1z*qw1 + 2*p1x*qy1)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (2*p1y*qw1 + 2*p1x*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
    (4*qx1*uix + 2*qy1*uiy + 2*qz1*uiz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qy1*uix - 2*qw1*uiz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (2*qz1*uix + 2*qw1*uiy)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[5] = (2*p1z*qw1 + 2*p1y*qx1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (-2*p1x*qw1 + 2*p1y*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
    (2*qx1*uiy + 2*qw1*uiz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qx1*uix + 4*qy1*uiy + 2*qz1*uiz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (-2*qw1*uix + 2*qz1*uiy)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[6] = (-2*p1y*qw1 + 2*p1z*qx1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (2*p1x*qw1 + 2*p1z*qy1)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
    (-2*qw1*uiy + 2*qx1*uiz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qw1*uix + 2*qy1*uiz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (2*qx1*uix + 2*qy1*uiy + 4*qz1*uiz)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        break;

      case 3:
        Cq[0] = (-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
    2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz;
        Cq[1] = 2*(qx1*qy1 + qw1*qz1)*ujx + 
    (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
    2*(-(qw1*qx1) + qy1*qz1)*ujz;
        Cq[2] = 2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
    (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz;
        Cq[3] = (4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
    (4*qw1*ujx - 2*qz1*ujy + 2*qy1*ujz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qz1*ujx + 4*qw1*ujy - 2*qx1*ujz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (-2*qy1*ujx + 2*qx1*ujy + 4*qw1*ujz)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2); 
        Cq[4] = (4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (-2*p1z*qw1 + 2*p1x*qy1)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (2*p1y*qw1 + 2*p1x*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
    (4*qx1*ujx + 2*qy1*ujy + 2*qz1*ujz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qy1*ujx - 2*qw1*ujz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (2*qz1*ujx + 2*qw1*ujy)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2); 
        Cq[5] = (2*p1z*qw1 + 2*p1y*qx1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (-2*p1x*qw1 + 2*p1y*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
    (2*qx1*ujy + 2*qw1*ujz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qx1*ujx + 4*qy1*ujy + 2*qz1*ujz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (-2*qw1*ujx + 2*qz1*ujy)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2); 
        Cq[6] = (-2*p1y*qw1 + 2*p1z*qx1)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (2*p1x*qw1 + 2*p1z*qy1)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
    (-2*qw1*ujy + 2*qx1*ujz)*
     (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
       p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
       2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + 
       x1 - x2) + (2*qw1*ujx + 2*qy1*ujz)*
     (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
       p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*p1x*(qx1*qy1 + qw1*qz1) + 
       2*p1z*(-(qw1*qx1) + qy1*qz1) - 
       2*p2x*(qx2*qy2 + qw2*qz2) - 
       2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2) + 
    (2*qx1*ujx + 2*qy1*ujy + 4*qz1*ujz)*
     (2*p1x*(-(qw1*qy1) + qx1*qz1) + 
       2*p1y*(qw1*qx1 + qy1*qz1) + 
       p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
       2*p2x*(-(qw2*qy2) + qx2*qz2) - 
       2*p2y*(qw2*qx2 + qy2*qz2) - 
       p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2); 
        break;

      case 4:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (4*uix*qw1 + 2*uiz*qy1 - 2*uiy*qz1)*
     (ujx*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*ujy*(qx2*qy2 - qw2*qz2) + 2*ujz*(qw2*qy2 + qx2*qz2)) + 
    (4*uiy*qw1 - 2*uiz*qx1 + 2*uix*qz1)*
     (ujy*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*ujx*(qx2*qy2 + qw2*qz2) + 2*ujz*(-(qw2*qx2) + qy2*qz2))
      + (4*uiz*qw1 + 2*uiy*qx1 - 2*uix*qy1)*
     (2*ujx*(-(qw2*qy2) + qx2*qz2) + 
       2*ujy*(qw2*qx2 + qy2*qz2) + 
       ujz*(-1 + 2*(qw2*qw2 + qz2*qz2))); 
        Cq[4] = (4*uix*qx1 + 2*uiy*qy1 + 2*uiz*qz1)*
     (ujx*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*ujy*(qx2*qy2 - qw2*qz2) + 2*ujz*(qw2*qy2 + qx2*qz2)) + 
    (-2*uiz*qw1 + 2*uix*qy1)*
     (ujy*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*ujx*(qx2*qy2 + qw2*qz2) + 2*ujz*(-(qw2*qx2) + qy2*qz2))
      + (2*uiy*qw1 + 2*uix*qz1)*
     (2*ujx*(-(qw2*qy2) + qx2*qz2) + 
       2*ujy*(qw2*qx2 + qy2*qz2) + 
       ujz*(-1 + 2*(qw2*qw2 + qz2*qz2)));
        Cq[5] = (2*uiz*qw1 + 2*uiy*qx1)*
     (ujx*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*ujy*(qx2*qy2 - qw2*qz2) + 2*ujz*(qw2*qy2 + qx2*qz2)) + 
    (2*uix*qx1 + 4*uiy*qy1 + 2*uiz*qz1)*
     (ujy*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*ujx*(qx2*qy2 + qw2*qz2) + 2*ujz*(-(qw2*qx2) + qy2*qz2))
      + (-2*uix*qw1 + 2*uiy*qz1)*
     (2*ujx*(-(qw2*qy2) + qx2*qz2) + 
       2*ujy*(qw2*qx2 + qy2*qz2) + 
       ujz*(-1 + 2*(qw2*qw2 + qz2*qz2)));
        Cq[6] = (-2*uiy*qw1 + 2*uiz*qx1)*
     (ujx*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
       2*ujy*(qx2*qy2 - qw2*qz2) + 2*ujz*(qw2*qy2 + qx2*qz2)) + 
    (2*uix*qw1 + 2*uiz*qy1)*
     (ujy*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
       2*ujx*(qx2*qy2 + qw2*qz2) + 2*ujz*(-(qw2*qx2) + qy2*qz2))
      + (2*uix*qx1 + 2*uiy*qy1 + 4*uiz*qz1)*
     (2*ujx*(-(qw2*qy2) + qx2*qz2) + 
       2*ujy*(qw2*qx2 + qy2*qz2) + 
       ujz*(-1 + 2*(qw2*qw2 + qz2*qz2)));
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
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

      case 1:
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

      case 2:
        Cq[0] = -((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix) - 
    2*(qx1*qy1 - qw1*qz1)*uiy - 2*(qw1*qy1 + qx1*qz1)*uiz;
        Cq[1] = -2*(qx1*qy1 + qw1*qz1)*uix - 
    (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy - 
    2*(-(qw1*qx1) + qy1*qz1)*uiz;
        Cq[2] = -2*(-(qw1*qy1) + qx1*qz1)*uix - 
    2*(qw1*qx1 + qy1*qz1)*uiy - 
    (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz;
        Cq[3] = (-4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (-4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (-4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[4] = (-4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (2*p2z*qw2 - 2*p2x*qy2)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (-2*p2y*qw2 - 2*p2x*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[5] = (-2*p2z*qw2 - 2*p2y*qx2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (-2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (2*p2x*qw2 - 2*p2y*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[6] = (2*p2y*qw2 - 2*p2z*qx2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
       2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
    (-2*p2x*qw2 - 2*p2z*qy2)*
     (2*(qx1*qy1 + qw1*qz1)*uix + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
       2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
    (-2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*uix + 
       2*(qw1*qx1 + qy1*qz1)*uiy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        break;

      case 3:
        Cq[0] = -((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx) - 
    2*(qx1*qy1 - qw1*qz1)*ujy - 2*(qw1*qy1 + qx1*qz1)*ujz;
        Cq[1] = -2*(qx1*qy1 + qw1*qz1)*ujx - 
    (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy - 
    2*(-(qw1*qx1) + qy1*qz1)*ujz;
        Cq[2] = -2*(-(qw1*qy1) + qx1*qz1)*ujx - 
    2*(qw1*qx1 + qy1*qz1)*ujy - 
    (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz;
        Cq[3] = (-4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (-4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (-4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz); 
        Cq[4] = (-4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (2*p2z*qw2 - 2*p2x*qy2)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (-2*p2y*qw2 - 2*p2x*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz); 
        Cq[5] = (-2*p2z*qw2 - 2*p2y*qx2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (-2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (2*p2x*qw2 - 2*p2y*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz); 
        Cq[6] = (2*p2y*qw2 - 2*p2z*qx2)*
     ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
       2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
    (-2*p2x*qw2 - 2*p2z*qy2)*
     (2*(qx1*qy1 + qw1*qz1)*ujx + 
       (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
       2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
    (-2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2)*
     (2*(-(qw1*qy1) + qx1*qz1)*ujx + 
       2*(qw1*qx1 + qy1*qz1)*ujy + 
       (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz); 
        break;

      case 4:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = (4*ujz*qw2 + 2*ujy*qx2 - 2*ujx*qy2)*
     (2*uix*(-(qw1*qy1) + qx1*qz1) + 
       2*uiy*(qw1*qx1 + qy1*qz1) + 
       uiz*(-1 + 2*(qw1*qw1 + qz1*qz1))) + 
    (uiy*(-1 + 2*(qw1*qw1 + qy1*qy1)) + 
       2*uix*(qx1*qy1 + qw1*qz1) + 2*uiz*(-(qw1*qx1) + qy1*qz1))
      *(4*ujy*qw2 - 2*ujz*qx2 + 2*ujx*qz2) + 
    (uix*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
       2*uiy*(qx1*qy1 - qw1*qz1) + 2*uiz*(qw1*qy1 + qx1*qz1))*
     (4*ujx*qw2 + 2*ujz*qy2 - 2*ujy*qz2); 
        Cq[4] = (-2*ujz*qw2 + 2*ujx*qy2)*
     (uiy*(-1 + 2*(qw1*qw1 + qy1*qy1)) + 
       2*uix*(qx1*qy1 + qw1*qz1) + 2*uiz*(-(qw1*qx1) + qy1*qz1))
      + (2*uix*(-(qw1*qy1) + qx1*qz1) + 
       2*uiy*(qw1*qx1 + qy1*qz1) + 
       uiz*(-1 + 2*(qw1*qw1 + qz1*qz1)))*
     (2*ujy*qw2 + 2*ujx*qz2) + 
    (uix*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
       2*uiy*(qx1*qy1 - qw1*qz1) + 2*uiz*(qw1*qy1 + qx1*qz1))*
     (4*ujx*qx2 + 2*ujy*qy2 + 2*ujz*qz2);
        Cq[5] = (2*ujz*qw2 + 2*ujy*qx2)*
     (uix*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
       2*uiy*(qx1*qy1 - qw1*qz1) + 2*uiz*(qw1*qy1 + qx1*qz1)) + 
    (2*uix*(-(qw1*qy1) + qx1*qz1) + 2*uiy*(qw1*qx1 + qy1*qz1) + 
       uiz*(-1 + 2*(qw1*qw1 + qz1*qz1)))*
     (-2*ujx*qw2 + 2*ujy*qz2) + 
    (uiy*(-1 + 2*(qw1*qw1 + qy1*qy1)) + 
       2*uix*(qx1*qy1 + qw1*qz1) + 2*uiz*(-(qw1*qx1) + qy1*qz1))
      *(2*ujx*qx2 + 4*ujy*qy2 + 2*ujz*qz2);
        Cq[6] = (-2*ujy*qw2 + 2*ujz*qx2)*
     (uix*(-1 + 2*qw1*qw1 + 2*qx1*qx1) + 
       2*uiy*(qx1*qy1 - qw1*qz1) + 2*uiz*(qw1*qy1 + qx1*qz1)) + 
    (2*ujx*qw2 + 2*ujz*qy2)*
     (uiy*(-1 + 2*(qw1*qw1 + qy1*qy1)) + 
       2*uix*(qx1*qy1 + qw1*qz1) + 2*uiz*(-(qw1*qx1) + qy1*qz1))
      + (2*uix*(-(qw1*qy1) + qx1*qz1) + 
       2*uiy*(qw1*qx1 + qy1*qz1) + 
       uiz*(-1 + 2*(qw1*qw1 + qz1*qz1)))*
     (2*ujx*qx2 + 2*ujy*qy2 + 4*ujz*qz2);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Calculates the time derivative of the constraint Jacobian
void PrismaticJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // setup constants for calculations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Quat qd1 = Quat::deriv(q1, inner->get_avel());
  const Quat qd2 = Quat::deriv(q2, outer->get_avel());
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
  const double x1 = inner->get_position()[X];
  const double y1 = inner->get_position()[Y];
  const double z1 = inner->get_position()[Z];
  const double x2 = outer->get_position()[X];
  const double y2 = outer->get_position()[Y];
  const double z2 = outer->get_position()[Z];
  const double dx1 = inner->get_lvel()[X];
  const double dy1 = inner->get_lvel()[Y];
  const double dz1 = inner->get_lvel()[Z];
  const double dx2 = outer->get_lvel()[X];
  const double dy2 = outer->get_lvel()[Y];
  const double dz2 = outer->get_lvel()[Z];
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
  const double dqw1 = qd1.w;
  const double dqx1 = qd1.x;
  const double dqy1 = qd1.y;
  const double dqz1 = qd1.z;
  const double dqw2 = qd2.w;
  const double dqx2 = qd2.x;
  const double dqy2 = qd2.y;
  const double dqz2 = qd2.z;
  const double uix = _ui[X];
  const double uiy = _ui[Y];
  const double uiz = _ui[Z];
  const double ujx = _uj[X];
  const double ujy = _uj[Y];
  const double ujz = _uj[Z];
  const double v2x = _v2[X];
  const double v2y = _v2[Y];
  const double v2z = _v2[Z];

  // setup the constraint equations (from Shabana, p. 437), eq. 7.179
  if (body == inner)
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
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

      case 1:
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

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*
    (-2*qy1*uix + 2*qx1*uiy + 4*qw1*uiz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qz1*uix + 4*qw1*uiy - 2*qx1*uiz) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qw1*uix - 2*qz1*uiy + 2*qy1*uiz) + 
   (4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (4*dqw1*p1x - 2*dqz1*p1y + 2*dqy1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (2*dqz1*p1x + 4*dqw1*p1y - 2*dqx1*p1z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (-2*dqy1*p1x + 2*dqx1*p1y + 4*dqw1*p1z)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
   (4*dqw1*uix - 2*dqz1*uiy + 2*dqy1*uiz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqz1*uix + 4*dqw1*uiy - 2*dqx1*uiz)*
    (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (-2*dqy1*uix + 2*dqx1*uiy + 4*dqw1*uiz)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[4] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qz1*uix + 2*qw1*uiy) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qy1*uix - 2*qw1*uiz) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qx1*uix + 2*qy1*uiy + 2*qz1*uiz) + 
   (4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (-2*p1z*qw1 + 2*p1x*qy1)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (2*p1y*qw1 + 2*p1x*qz1)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (4*dqx1*p1x + 2*dqy1*p1y + 2*dqz1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (2*dqy1*p1x - 2*dqw1*p1z)*(2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (2*dqz1*p1x + 2*dqw1*p1y)*(2*(-(qw1*qy1) + qx1*qz1)*uix + 
      2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
   (4*dqx1*uix + 2*dqy1*uiy + 2*dqz1*uiz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqy1*uix - 2*dqw1*uiz)*(p1y*
       (-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (2*dqz1*uix + 2*dqw1*uiy)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[5] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qw1*uix + 2*qz1*uiy) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (2*qx1*uiy + 2*qw1*uiz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qx1*uix + 4*qy1*uiy + 2*qz1*uiz) + 
   (2*p1z*qw1 + 2*p1y*qx1)*((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (-2*p1x*qw1 + 2*p1y*qz1)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (2*dqx1*p1y + 2*dqw1*p1z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       uix + 2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (2*dqx1*p1x + 4*dqy1*p1y + 2*dqz1*p1z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (-2*dqw1*p1x + 2*dqz1*p1y)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
   (2*dqx1*uiy + 2*dqw1*uiz)*(p1x*
       (-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqx1*uix + 4*dqy1*uiy + 2*dqz1*uiz)*
    (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (-2*dqw1*uix + 2*dqz1*uiy)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[6] = (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (-2*qw1*uiy + 2*qx1*uiz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qw1*uix + 2*qy1*uiz) + 
   (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*
    (2*qx1*uix + 2*qy1*uiy + 4*qz1*uiz) + 
   (-2*p1y*qw1 + 2*p1z*qx1)*((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (2*p1x*qw1 + 2*p1z*qy1)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (-2*dqw1*p1y + 2*dqx1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (2*dqw1*p1x + 2*dqy1*p1z)*(2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (2*dqx1*p1x + 2*dqy1*p1y + 4*dqz1*p1z)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz) + 
   (-2*dqw1*uiy + 2*dqx1*uiz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqw1*uix + 2*dqy1*uiz)*(p1y*
       (-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (2*dqx1*uix + 2*dqy1*uiy + 4*dqz1*uiz)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*
    (-2*qy1*ujx + 2*qx1*ujy + 4*qw1*ujz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qz1*ujx + 4*qw1*ujy - 2*qx1*ujz) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qw1*ujx - 2*qz1*ujy + 2*qy1*ujz) + 
   (4*p1x*qw1 + 2*p1z*qy1 - 2*p1y*qz1)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (4*p1y*qw1 - 2*p1z*qx1 + 2*p1x*qz1)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (4*p1z*qw1 + 2*p1y*qx1 - 2*p1x*qy1)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (4*dqw1*p1x - 2*dqz1*p1y + 2*dqy1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (2*dqz1*p1x + 4*dqw1*p1y - 2*dqx1*p1z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (-2*dqy1*p1x + 2*dqx1*p1y + 4*dqw1*p1z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
   (4*dqw1*ujx - 2*dqz1*ujy + 2*dqy1*ujz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqz1*ujx + 4*dqw1*ujy - 2*dqx1*ujz)*
    (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (-2*dqy1*ujx + 2*dqx1*ujy + 4*dqw1*ujz)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[4] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qz1*ujx + 2*qw1*ujy) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qy1*ujx - 2*qw1*ujz) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qx1*ujx + 2*qy1*ujy + 2*qz1*ujz) + 
   (4*p1x*qx1 + 2*p1y*qy1 + 2*p1z*qz1)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (-2*p1z*qw1 + 2*p1x*qy1)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (2*p1y*qw1 + 2*p1x*qz1)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (4*dqx1*p1x + 2*dqy1*p1y + 2*dqz1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (2*dqy1*p1x - 2*dqw1*p1z)*(2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (2*dqz1*p1x + 2*dqw1*p1y)*(2*(-(qw1*qy1) + qx1*qz1)*ujx + 
      2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
   (4*dqx1*ujx + 2*dqy1*ujy + 2*dqz1*ujz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqy1*ujx - 2*dqw1*ujz)*(p1y*
       (-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (2*dqz1*ujx + 2*dqw1*ujy)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[5] = (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qw1*ujx + 2*qz1*ujy) + 
   (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (2*qx1*ujy + 2*qw1*ujz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qx1*ujx + 4*qy1*ujy + 2*qz1*ujz) + 
   (2*p1z*qw1 + 2*p1y*qx1)*((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (2*p1x*qx1 + 4*p1y*qy1 + 2*p1z*qz1)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (-2*p1x*qw1 + 2*p1y*qz1)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (2*dqx1*p1y + 2*dqw1*p1z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       ujx + 2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (2*dqx1*p1x + 4*dqy1*p1y + 2*dqz1*p1z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (-2*dqw1*p1x + 2*dqz1*p1y)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
   (2*dqx1*ujy + 2*dqw1*ujz)*(p1x*
       (-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqx1*ujx + 4*dqy1*ujy + 2*dqz1*ujz)*
    (p1y*(-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (-2*dqw1*ujx + 2*dqz1*ujy)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        Cq[6] = (dx1 - dx2 + p1x*(4*dqw1*qw1 + 4*dqx1*qx1) - 
      p2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*p1y*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1) + 
      2*p1z*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1) - 
      2*p2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) - 
      2*p2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (-2*qw1*ujy + 2*qx1*ujz) + 
   (dy1 - dy2 + 2*p1y*(2*dqw1*qw1 + 2*dqy1*qy1) - 
      2*p2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*p1x*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1) + 
      2*p1z*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1) - 
      2*p2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) - 
      2*p2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qw1*ujx + 2*qy1*ujz) + 
   (dz1 - dz2 + 2*p1x*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1) + 
      2*p1y*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1) + 
      2*p1z*(2*dqw1*qw1 + 2*dqz1*qz1) - 
      2*p2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) - 
      2*p2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) - 
      2*p2z*(2*dqw2*qw2 + 2*dqz2*qz2))*
    (2*qx1*ujx + 2*qy1*ujy + 4*qz1*ujz) + 
   (-2*p1y*qw1 + 2*p1z*qx1)*((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (2*p1x*qw1 + 2*p1z*qy1)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (2*p1x*qx1 + 2*p1y*qy1 + 4*p1z*qz1)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (-2*dqw1*p1y + 2*dqx1*p1z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (2*dqw1*p1x + 2*dqy1*p1z)*(2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (2*dqx1*p1x + 2*dqy1*p1y + 4*dqz1*p1z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz) + 
   (-2*dqw1*ujy + 2*dqx1*ujz)*
    (p1x*(-1 + 2*qw1*qw1 + 2*qx1*qx1) - 
      p2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*p1y*(qx1*qy1 - qw1*qz1) + 2*p1z*(qw1*qy1 + qx1*qz1) - 
      2*p2y*(qx2*qy2 - qw2*qz2) - 2*p2z*(qw2*qy2 + qx2*qz2) + x1 - x2) + 
   (2*dqw1*ujx + 2*dqy1*ujz)*(p1y*
       (-1 + 2*(qw1*qw1 + qy1*qy1)) - 
      p2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*p1x*(qx1*qy1 + qw1*qz1) + 2*p1z*(-(qw1*qx1) + qy1*qz1) - 
      2*p2x*(qx2*qy2 + qw2*qz2) - 2*p2z*(-(qw2*qx2) + qy2*qz2) + y1 - y2)
     + (2*dqx1*ujx + 2*dqy1*ujy + 4*dqz1*ujz)*
    (2*p1x*(-(qw1*qy1) + qx1*qz1) + 2*p1y*(qw1*qx1 + qy1*qz1) + 
      p1z*(-1 + 2*(qw1*qw1 + qz1*qz1)) - 
      2*p2x*(-(qw2*qy2) + qx2*qz2) - 2*p2y*(qw2*qx2 + qy2*qz2) - 
      p2z*(-1 + 2*(qw2*qw2 + qz2*qz2)) + z1 - z2);
        break;

      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (4*qw1*uix - 2*qz1*uiy + 2*qy1*uiz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*ujx + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*ujy + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*ujz) + 
   (2*qz1*uix + 4*qw1*uiy - 2*qx1*uiz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*ujx + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*ujy + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujz) + 
   (-2*qy1*uix + 2*qx1*uiy + 4*qw1*uiz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*ujx + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujy + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*ujz) + 
   (4*dqw1*uix - 2*dqz1*uiy + 2*dqy1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*ujx + 
      2*(qx2*qy2 - qw2*qz2)*ujy + 2*(qw2*qy2 + qx2*qz2)*ujz) + 
   (2*dqz1*uix + 4*dqw1*uiy - 2*dqx1*uiz)*
    (2*(qx2*qy2 + qw2*qz2)*ujx + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*ujy + 
      2*(-(qw2*qx2) + qy2*qz2)*ujz) + 
   (-2*dqy1*uix + 2*dqx1*uiy + 4*dqw1*uiz)*
    (2*(-(qw2*qy2) + qx2*qz2)*ujx + 2*(qw2*qx2 + qy2*qz2)*ujy + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*ujz);
        Cq[4] = (4*qx1*uix + 2*qy1*uiy + 2*qz1*uiz)*
    ((4*dqw2*qw2 + 4*dqx2*qx2)*ujx + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*ujy + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*ujz) + 
   (2*qy1*uix - 2*qw1*uiz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*ujx + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*ujy + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujz) + 
   (2*qz1*uix + 2*qw1*uiy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*ujx + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujy + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*ujz) + 
   (4*dqx1*uix + 2*dqy1*uiy + 2*dqz1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*ujx + 
      2*(qx2*qy2 - qw2*qz2)*ujy + 2*(qw2*qy2 + qx2*qz2)*ujz) + 
   (2*dqy1*uix - 2*dqw1*uiz)*(2*(qx2*qy2 + qw2*qz2)*ujx + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*ujy + 
      2*(-(qw2*qx2) + qy2*qz2)*ujz) + 
   (2*dqz1*uix + 2*dqw1*uiy)*(2*(-(qw2*qy2) + qx2*qz2)*ujx + 
      2*(qw2*qx2 + qy2*qz2)*ujy + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*ujz);
        Cq[5] = (2*qx1*uiy + 2*qw1*uiz)*((4*dqw2*qw2 + 4*dqx2*qx2)*ujx + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*ujy + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*ujz) + 
   (2*qx1*uix + 4*qy1*uiy + 2*qz1*uiz)*
    (2*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*ujx + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*ujy + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujz) + 
   (-2*qw1*uix + 2*qz1*uiy)*(2*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*ujx + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujy + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*ujz) + 
   (2*dqx1*uiy + 2*dqw1*uiz)*((-1 + 2*qw2*qw2 + 2*qx2*qx2)*
       ujx + 2*(qx2*qy2 - qw2*qz2)*ujy + 2*(qw2*qy2 + qx2*qz2)*ujz) + 
   (2*dqx1*uix + 4*dqy1*uiy + 2*dqz1*uiz)*
    (2*(qx2*qy2 + qw2*qz2)*ujx + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*ujy + 
      2*(-(qw2*qx2) + qy2*qz2)*ujz) + 
   (-2*dqw1*uix + 2*dqz1*uiy)*
    (2*(-(qw2*qy2) + qx2*qz2)*ujx + 2*(qw2*qx2 + qy2*qz2)*ujy + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*ujz);
        Cq[6] = (-2*qw1*uiy + 2*qx1*uiz)*((4*dqw2*qw2 + 4*dqx2*qx2)*ujx + 
      2*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2)*ujy + 
      2*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2)*ujz) + 
   (2*qw1*uix + 2*qy1*uiz)*(2*
       (dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2)*ujx + 
      2*(2*dqw2*qw2 + 2*dqy2*qy2)*ujy + 
      2*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujz) + 
   (2*qx1*uix + 2*qy1*uiy + 4*qz1*uiz)*
    (2*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2)*ujx + 
      2*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2)*ujy + 
      2*(2*dqw2*qw2 + 2*dqz2*qz2)*ujz) + 
   (-2*dqw1*uiy + 2*dqx1*uiz)*
    ((-1 + 2*qw2*qw2 + 2*qx2*qx2)*ujx + 
      2*(qx2*qy2 - qw2*qz2)*ujy + 2*(qw2*qy2 + qx2*qz2)*ujz) + 
   (2*dqw1*uix + 2*dqy1*uiz)*(2*(qx2*qy2 + qw2*qz2)*ujx + 
      (-1 + 2*(qw2*qw2 + qy2*qy2))*ujy + 
      2*(-(qw2*qx2) + qy2*qz2)*ujz) + 
   (2*dqx1*uix + 2*dqy1*uiy + 4*dqz1*uiz)*
    (2*(-(qw2*qy2) + qx2*qz2)*ujx + 2*(qw2*qx2 + qy2*qz2)*ujy + 
      (-1 + 2*(qw2*qw2 + qz2*qz2))*ujz);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
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

      case 1:
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

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (-4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (-4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (-4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (-4*dqw2*p2x + 2*dqz2*p2y - 2*dqy2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (-2*dqz2*p2x - 4*dqw2*p2y + 2*dqx2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (2*dqy2*p2x - 2*dqx2*p2y - 4*dqw2*p2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[4] = (-4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (2*p2z*qw2 - 2*p2x*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (-2*p2y*qw2 - 2*p2x*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (-4*dqx2*p2x - 2*dqy2*p2y - 2*dqz2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (-2*dqy2*p2x + 2*dqw2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (-2*dqz2*p2x - 2*dqw2*p2y)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[5] = (-2*p2z*qw2 - 2*p2y*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (-2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (2*p2x*qw2 - 2*p2y*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (-2*dqx2*p2y - 2*dqw2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (-2*dqx2*p2x - 4*dqy2*p2y - 2*dqz2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (2*dqw2*p2x - 2*dqz2*p2y)*(2*(-(qw1*qy1) + qx1*qz1)*uix + 
      2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        Cq[6] = (2*p2y*qw2 - 2*p2z*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz) + 
   (-2*p2x*qw2 - 2*p2z*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz) + 
   (-2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz) + 
   (2*dqw2*p2y - 2*dqx2*p2z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       uix + 2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz) + 
   (-2*dqw2*p2x - 2*dqy2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz) + 
   (-2*dqx2*p2x - 2*dqy2*p2y - 4*dqz2*p2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz);
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (-4*p2x*qw2 - 2*p2z*qy2 + 2*p2y*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (-4*p2y*qw2 + 2*p2z*qx2 - 2*p2x*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (-4*p2z*qw2 - 2*p2y*qx2 + 2*p2x*qy2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (-4*dqw2*p2x + 2*dqz2*p2y - 2*dqy2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (-2*dqz2*p2x - 4*dqw2*p2y + 2*dqx2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (2*dqy2*p2x - 2*dqx2*p2y - 4*dqw2*p2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz);
        Cq[4] = (-4*p2x*qx2 - 2*p2y*qy2 - 2*p2z*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (2*p2z*qw2 - 2*p2x*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (-2*p2y*qw2 - 2*p2x*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (-4*dqx2*p2x - 2*dqy2*p2y - 2*dqz2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (-2*dqy2*p2x + 2*dqw2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (-2*dqz2*p2x - 2*dqw2*p2y)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz);
        Cq[5] = (-2*p2z*qw2 - 2*p2y*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (-2*p2x*qx2 - 4*p2y*qy2 - 2*p2z*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (2*p2x*qw2 - 2*p2y*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (-2*dqx2*p2y - 2*dqw2*p2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ujx + 
      2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (-2*dqx2*p2x - 4*dqy2*p2y - 2*dqz2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (2*dqw2*p2x - 2*dqz2*p2y)*(2*(-(qw1*qy1) + qx1*qz1)*ujx + 
      2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz);
        Cq[6] = (2*p2y*qw2 - 2*p2z*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ujx + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*ujy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*ujz) + 
   (-2*p2x*qw2 - 2*p2z*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ujx + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*ujy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujz) + 
   (-2*p2x*qx2 - 2*p2y*qy2 - 4*p2z*qz2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ujx + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*ujy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*ujz) + 
   (2*dqw2*p2y - 2*dqx2*p2z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       ujx + 2*(qx1*qy1 - qw1*qz1)*ujy + 2*(qw1*qy1 + qx1*qz1)*ujz) + 
   (-2*dqw2*p2x - 2*dqy2*p2z)*
    (2*(qx1*qy1 + qw1*qz1)*ujx + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*ujy + 
      2*(-(qw1*qx1) + qy1*qz1)*ujz) + 
   (-2*dqx2*p2x - 2*dqy2*p2y - 4*dqz2*p2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ujx + 2*(qw1*qx1 + qy1*qz1)*ujy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*ujz);
        break;

      case 4:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (-2*dqy2*ujx + 2*dqx2*ujy + 4*dqw2*ujz) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*
    (2*dqz2*ujx + 4*dqw2*ujy - 2*dqx2*ujz) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (4*dqw2*ujx - 2*dqz2*ujy + 2*dqy2*ujz) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*
    (-2*qy2*ujx + 2*qx2*ujy + 4*qw2*ujz) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qz2*ujx + 4*qw2*ujy - 2*qx2*ujz) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (4*qw2*ujx - 2*qz2*ujy + 2*qy2*ujz);
        Cq[4] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (2*dqz2*ujx + 2*dqw2*ujy) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*(2*qz2*ujx + 2*qw2*ujy) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*dqy2*ujx - 2*dqw2*ujz) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (4*dqx2*ujx + 2*dqy2*ujy + 2*dqz2*ujz) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qy2*ujx - 2*qw2*ujz) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (4*qx2*ujx + 2*qy2*ujy + 2*qz2*ujz);
        Cq[5] = (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (-2*dqw2*ujx + 2*dqz2*ujy) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*(-2*qw2*ujx + 2*qz2*ujy) + 
   ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (2*dqx2*ujy + 2*dqw2*ujz) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*
    (2*dqx2*ujx + 4*dqy2*ujy + 2*dqz2*ujz) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (2*qx2*ujy + 2*qw2*ujz) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qx2*ujx + 4*qy2*ujy + 2*qz2*ujz);
        Cq[6] = ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*uix + 
      2*(qx1*qy1 - qw1*qz1)*uiy + 2*(qw1*qy1 + qx1*qz1)*uiz)*
    (-2*dqw2*ujy + 2*dqx2*ujz) + 
   (2*(qx1*qy1 + qw1*qz1)*uix + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uiy + 
      2*(-(qw1*qx1) + qy1*qz1)*uiz)*(2*dqw2*ujx + 2*dqy2*ujz) + 
   (2*(-(qw1*qy1) + qx1*qz1)*uix + 2*(qw1*qx1 + qy1*qz1)*uiy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uiz)*
    (2*dqx2*ujx + 2*dqy2*ujy + 4*dqz2*ujz) + 
   ((4*dqw1*qw1 + 4*dqx1*qx1)*uix + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uiy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uiz)*
    (-2*qw2*ujy + 2*qx2*ujz) + 
   (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*uix + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uiy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiz)*
    (2*qw2*ujx + 2*qy2*ujz) + 
   (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*uix + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uiy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uiz)*
    (2*qx2*ujx + 2*qy2*ujy + 4*qz2*ujz);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }    
  }
}

/// Evaluates the constraint equations
void PrismaticJoint::evaluate_constraints(double C[])
{
  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], p. 437; some variable names
  // have been altered

  // get v1 in global coordinates
  Vector3 v1 = get_axis_global();

  // determine axis in global coordinates
  Vector3 v2 = outer->get_transform().mult_vector(_v2);

  // determine v1i, v1j
  Vector3 v1i, v1j;
  Vector3::determine_orthonormal_basis(v1, v1i, v1j);

  // determine h1 and h2
  Vector3 h1 = inner->get_transform().mult_vector(_ui);
  Vector3 h2 = outer->get_transform().mult_vector(_uj);

  // determine the global positions of the attachment points and subtract them
  const Vector3& p1 = get_position_global(false); 
  const Vector3& p2 = get_position_global(true); 
  Vector3 r12 = p1 - p2; 

  // evaluate the constraint equations
  C[0] = v1i.dot(v2);
  C[1] = v1j.dot(v2);
  C[2] = v1i.dot(r12);
  C[3] = v1j.dot(r12);
  C[4] = h1.dot(h2);
}

/// Implements Base::load_from_xml()
void PrismaticJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "PrismaticJoint") == 0);

  // read the local joint axis
  const XMLAttrib* laxis_attrib = node->get_attrib("local-axis");
  if (laxis_attrib)
  {
    Vector3 laxis;
    laxis_attrib->get_vector_value(laxis);
    set_axis_local(laxis);
  }

  // read the global joint axis, if given
  const XMLAttrib* gaxis_attrib = node->get_attrib("global-axis");
  if (gaxis_attrib)
  {
    Vector3 gaxis;
    gaxis_attrib->get_vector_value(gaxis);
    set_axis_global(gaxis);  
  }

  // set the joint tare
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void PrismaticJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get the majority of the info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "PrismaticJoint";

  // save the local joint axis
  node->attribs.insert(XMLAttrib("local-axis", _u));
}

