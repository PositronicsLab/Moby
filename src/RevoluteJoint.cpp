/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/RevoluteJoint.h>

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
  _u = ZEROS_3;
  _ui = ZEROS_3;
  _uj = ZEROS_3;
  _v2 = ZEROS_3;

  // set the translation to zero
  _T.set_translation(ZEROS_3);

  // setup the spatial axis derivative to zero
  _si_deriv = SMatrix6N(1);
  _si_deriv.set_column(0, SVector6(0,0,0,0,0,0));
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
  _u = ZEROS_3;
  _ui = ZEROS_3;
  _uj = ZEROS_3;
  _v2 = ZEROS_3;

  // set the translation to zero
  _T.set_translation(ZEROS_3);

  // setup the spatial axis derivative to zero
  _si_deriv = SMatrix6N(1);
  _si_deriv.set_column(0, SVector6(0,0,0,0,0,0));
}  

/// Gets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
Vector3 RevoluteJoint::get_axis_global() const
{
  // get the inboard link 
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "RevoluteJoint::get_axis_global() - attempt to get axis w/o inboard link!" << std::endl;
    return ZEROS_3;
  }

  // get the transform for the inboard link
  const Matrix4& T = inboard_link->get_transform();
  
  // transform into global coordinates and return
  return T.transpose_mult_vector(_u);
}

/// Sets the local axis of rotation for this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void RevoluteJoint::set_axis_local(const Vector3& axis) 
{ 
  // normalize the joint axis, in case the caller didn't 
  Vector3 naxis = Vector3::normalize(axis);

  // get the inboard and outboard links 
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();
  if (!inner || !outer)
    throw std::runtime_error("Unable to set joint axis without both links set!");

  // set joint axis in inner link frame
  _u = naxis; 
  update_spatial_axes(); 

  // setup ui and uj
  Vector3::determine_orthonormal_basis(_u, _ui, _uj);
  
  // set joint axis in outer link frame
  _v2 = inner->get_transform().mult_vector(naxis);
  _v2 = outer->get_transform().transpose_mult_vector(_v2);
}        

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void RevoluteJoint::set_axis_global(const Vector3& axis)
{
  // normalize the joint axis, in case the caller didn't 
  Vector3 naxis = Vector3::normalize(axis);

  // get the inboard and outboard links
  RigidBodyPtr inboard_link = get_inboard_link();
  RigidBodyPtr outboard_link = get_outboard_link();

  // make sure that the links have been set
  if (!inboard_link || !outboard_link)
    throw std::runtime_error("Unable to set joint axis without links being set!");

  // transform the axis into local coordinates
  _u = inboard_link->get_transform().transpose_mult_vector(naxis);
  _v2 = outboard_link->get_transform().transpose_mult_vector(naxis);

  // setup ui and uj
  Vector3::determine_orthonormal_basis(_u, _ui, _uj);
  
  // update the spatial axis
  update_spatial_axes();
}

/// Updates the spatial axis for this joint
void RevoluteJoint::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // call parent method
  Joint::update_spatial_axes();

  // make sure the outboard link exists
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard || !outboard)
    return;

  try
  {
    // get the joint to com vector in outer link coordinates
    const Vector3& di = outboard->get_inner_joint_data(inboard).joint_to_com_vec;

    // get the joint axis in outer link coordinates
    Vector3 u0 = inboard->get_transform().mult_vector(_u);
    Vector3 ui = outboard->get_transform().transpose_mult_vector(u0);

    // update the spatial axis in link coordinates
    Vector3 x = Vector3::cross(ui, di);
    SVector6 si_vec;
    si_vec.set_upper(ui);
    si_vec.set_lower(x);
    _si.set_column(0, si_vec);

    // setup s_bar
    calc_s_bar_from_si();
  }
  catch (std::runtime_error e)
  {
    // do nothing -- joint data has not yet been set in the link
  }
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void RevoluteJoint::determine_Q()
{
  // get the inboard and outboard link pointers
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  
  // verify that the inboard and outboard links are set
  if (!inboard || !outboard)
  {
    std::cerr << "RevoluteJoint::determine_Q() called on NULL inboard and/or outboard links!" << std::endl;
    assert(false);
    return;
  }

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - 1.0) > NEAR_ZERO)
  {
    std::cerr << "RevoluteJoint::determine_Q() warning: some axes undefined; aborting..." << std::endl;
    return;
  }

  // get the link transforms
  Matrix3 R_inboard, R_outboard;
  inboard->get_transform().get_rotation(&R_inboard);
  outboard->get_transform().get_rotation(&R_outboard);

  // determine the joint transformation
  Matrix3 R_local = R_inboard.transpose_mult(R_outboard);
  AAngle aa(&R_local, &_u);

  // set Q
  this->_q_tare.resize(num_dof());
  this->_q_tare[DOF_1] = aa.angle;
}

/// Gets the (local) transform for this joint
const Matrix4& RevoluteJoint::get_transform()
{
  // note that translation is set to zero in the constructors
  AAngle a(&_u, this->q[DOF_1]+this->_q_tare[DOF_1]);
  _T.set_rotation(&a);
  assert(_T.get_translation().norm() < NEAR_ZERO);

  return _T;
}

/// Gets the derivative for the spatial axes for this joint
const SMatrix6N& RevoluteJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  return _si_deriv;
}

/// Computes the constraint Jacobian with respect to a body
void RevoluteJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _v2) is set
  if (_u.norm_sq() < std::numeric_limits<Real>::epsilon())
    throw std::runtime_error("Revolute joint axis has not been set; set before calling dynamics functions.");

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (Real) 0.0;
    return;
  }

  // setup constants for calculations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec;
  const Real p1x = p1[X];
  const Real p1y = p1[Y];
  const Real p1z = p1[Z];
  const Real p2x = -p2[X];
  const Real p2y = -p2[Y];
  const Real p2z = -p2[Z];
  const Real qw1 = q1.w;
  const Real qx1 = q1.x;
  const Real qy1 = q1.y;
  const Real qz1 = q1.z;
  const Real qw2 = q2.w;
  const Real qx2 = q2.x;
  const Real qy2 = q2.y;
  const Real qz2 = q2.z;
  const Real uix = _ui[X];
  const Real uiy = _ui[Y];
  const Real uiz = _ui[Z];
  const Real ujx = _uj[X];
  const Real ujy = _uj[Y];
  const Real ujz = _uj[Z];
  const Real v2x = _v2[X];
  const Real v2y = _v2[Y];
  const Real v2z = _v2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
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
}

/// Computes the time derivative of the constraint Jacobian with respect to a body
void RevoluteJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _v2) is set
  if (_u.norm_sq() < std::numeric_limits<Real>::epsilon())
    throw std::runtime_error("Revolute joint axis has not been set; set before calling dynamics functions.");

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (Real) 0.0;
    return;
  }

  // setup constants for calculations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Quat qd1 = Quat::deriv(q1, inner->get_avel());
  const Quat qd2 = Quat::deriv(q2, outer->get_avel());
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec;
  const Real dqw1 = qd1.w;
  const Real dqx1 = qd1.x;
  const Real dqy1 = qd1.y;
  const Real dqz1 = qd1.z;
  const Real dqw2 = qd2.w;
  const Real dqx2 = qd2.x;
  const Real dqy2 = qd2.y;
  const Real dqz2 = qd2.z;
  const Real p1x = p1[X];
  const Real p1y = p1[Y];
  const Real p1z = p1[Z];
  const Real p2x = -p2[X];
  const Real p2y = -p2[Y];
  const Real p2z = -p2[Z];
  const Real qw1 = q1.w;
  const Real qx1 = q1.x;
  const Real qy1 = q1.y;
  const Real qz1 = q1.z;
  const Real qw2 = q2.w;
  const Real qx2 = q2.x;
  const Real qy2 = q2.y;
  const Real qz2 = q2.z;
  const Real uix = _ui[X];
  const Real uiy = _ui[Y];
  const Real uiz = _ui[Z];
  const Real ujx = _uj[X];
  const Real ujy = _uj[Y];
  const Real ujz = _uj[Z];
  const Real v2x = _v2[X];
  const Real v2y = _v2[Y];
  const Real v2z = _v2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
  {
    switch (index)
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
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
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
        Cq[0] = (Real) 0.0;     
        Cq[1] = (Real) 0.0;      
        Cq[2] = (Real) 0.0;      
        Cq[3] = -4*p2x*dqw2 - 2*p2z*dqy2 + 2*p2y*dqz2; 
        Cq[4] = -4*p2x*dqx2 - 2*p2y*dqy2 - 2*p2z*dqz2; 
        Cq[5] = -2*p2z*dqw2 - 2*p2y*dqx2; 
        Cq[6] = 2*p2y*dqw2 - 2*p2z*dqx2; 
        break;

      case 1:
        Cq[0] = (Real) 0.0;      
        Cq[1] = (Real) 0.0;     
        Cq[2] = (Real) 0.0;      
        Cq[3] = -4*p2y*dqw2 + 2*p2z*dqx2 - 2*p2x*dqz2; 
        Cq[4] = 2*p2z*dqw2 - 2*p2x*dqy2; 
        Cq[5] = -2*p2x*dqx2 - 4*p2y*dqy2 - 2*p2z*dqz2; 
        Cq[6] = -2*p2x*dqw2 - 2*p2z*dqy2; 
        break;

      case 2:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
        Cq[3] = -4*p2z*dqw2 - 2*p2y*dqx2 + 2*p2x*dqy2;
        Cq[4] = -2*p2y*dqw2 - 2*p2x*dqz2;
        Cq[5] = 2*p2x*dqw2 - 2*p2y*dqz2;
        Cq[6] = -2*p2x*dqx2 - 2*p2y*dqy2 - 4*p2z*dqz2;
        break;

      case 3:
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
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
        Cq[0] = (Real) 0.0;
        Cq[1] = (Real) 0.0;
        Cq[2] = (Real) 0.0;
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
}

/// Evaluates the constraint equations
void RevoluteJoint::evaluate_constraints(Real C[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], p. 435-436; variable names
  // have been altered however

  // determine v1, v1i, v1j, and v2 (all in global coordinates)
  Vector3 v1i, v1j;
  Vector3 v1 = inner->get_transform().mult_vector(_u);
  Vector3 v2 = outer->get_transform().mult_vector(_v2);
  Vector3::determine_orthonormal_basis(v1, v1i, v1j);

  // determine the global positions of the attachment points and subtract them
  Vector3 r1 = get_position_global(false);
  Vector3 r2 = get_position_global(true);
  Vector3 r12 = r1 - r2; 

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = v1i.dot(v2);
  C[4] = v1j.dot(v2); 
}

/// Implements Base::load_from_xml()
void RevoluteJoint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "RevoluteJoint") == 0);

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

  // reset the joint position -- this will override any value of Q specified
  // in Joint::load_from_xml()
  determine_Q();
}

/// Implements Base::save_to_xml()
void RevoluteJoint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "RevoluteJoint";

  // save the local joint axis
  node->attribs.insert(XMLAttrib("local-axis", _u));
}

