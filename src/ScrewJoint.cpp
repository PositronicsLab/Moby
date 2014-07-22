/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>

using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The axis is set to [0 0 0]. The pitch is set to 1.0. The inboard and 
 * outboard links are set to NULL.
 */
ScrewJoint::ScrewJoint() : Joint()
{
  // init joint data
  init_data();

  // set the pitch
  _pitch = (double) 1.0;

  // init the joint axes
  _u.set_zero();
  _v2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis is set to [0 0 0] and the pitch is set to 1.0.  
 */
ScrewJoint::ScrewJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  // init joint data
  init_data();

  // set the pitch
  _pitch = (double) 1.0;

  // init the joint axes
  _u.set_zero();
  _v2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}  

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link 
 * into account; thus, if the orientation of the inboard link changes, then 
 * the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void ScrewJoint::set_axis(const Vector3d& axis)
{
  // normalize the joint axis, in case the caller didn't 
  Vector3d naxis = Vector3d::normalize(axis);

  // transform the axis as necessary
  _u = Pose3d::transform(axis.pose, _F, naxis);

  // get the inboard and outboard links
  RigidBodyPtr inboard_link = get_inboard_link();
  RigidBodyPtr outboard_link = get_outboard_link();

  // make sure that the links have been set
  if (!inboard_link || !outboard_link)
    throw std::runtime_error("Unable to set joint axis without links being set!");

  // transform the axis into local coordinates
  _u = inboard_link->get_transform().transpose_mult_vector(naxis);
  _v2 = outboard_link->get_transform().transpose_mult_vector(naxis);

  // update the spatial axis
  update_spatial_axes();
}

/// Updates the spatial axis for this joint
void ScrewJoint::update_spatial_axes()
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
    const Vector3d& di = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

    // get the joint axis in outer link coordinates
    Vector3d u0 = inboard->get_transform().mult_vector(_u);
    Vector3d ui = outboard->get_transform().transpose_mult_vector(u0);

    // update the spatial axis in link coordinates
    Vector3d x = Vector3d::cross(ui, di);
    _s[0].set_angular(ui);
    _s[0].set_lower(ui*_pitch + x);

    // setup s_bar
    calc_s_bar_from_s();
  }
  catch (std::runtime_error e)
  {
    // do nothing -- joint data has not yet been set in the link
  }
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void ScrewJoint::determine_q(VectorN& q)
{
  // get the inboard and outboard link pointers
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  
  // verify that the inboard and outboard links are set
  if (!inboard || !outboard)
  {
    std::cerr << "ScrewJoint::determine_Q() called on NULL inboard and/or outboard links!" << std::endl;
    assert(false);
    return;
  }

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - 1.0) > NEAR_ZERO)
  {
    std::cerr << "ScrewJoint::determine_Q() warning: some axes undefined; aborting..." << std::endl;
    return;
  }

  // get the attachment points on the link (global coords)
  Vector3d p1 = get_position_global(false);
  Vector3d p2 = get_position_global(true);

  // get the joint axis in the global frame
  Vector3d ug = inboard->get_transform().mult_vector(_u);

  // now, we'll project p2 onto the axis ug; points will be setup so that
  // ug passes through origin on inboard
  q.resize(num_dof());
  q[DOF_1] = ug.dot(p2-p1)/_pitch;
}

/// Gets the (local) transform for this joint
shared_ptr<const Pose3d> ScrewJoint::get_induced_pose()
{
/*
  // setup rotation 
  AAngle a(&_u, this->q[DOF_1]+this->_q_tare[DOF_1]);
  _T.set_rotation(&a);

  // setup translation
  _T.set_translation(_u*(this->q[DOF_1]+this->_q_tare[DOF_1]));
*/
  return _Fprime;
}

/// Gets the derivative for the spatial axes for this joint
const vector<SVelocityd>& ScrewJoint::get_spatial_axes_dot()
{
  return _s_dot;
}

/// Implements Base::load_from_xml()
void ScrewJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "ScrewJoint") == 0);

  // read the pitch
  XMLAttrib* pitch_attrib = node->get_attrib("pitch");
  if (pitch_attrib)
    _pitch = pitch_attrib->get_real_value(); 

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
void ScrewJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "ScrewJoint";

  // save the joint axis (global coordinates)
  Vector3dd u0 = Pose3d::transform(_u.pose, shared_ptr<const Pose3d>(), _u);
  node->attribs.insert(XMLAttrib("axis", u0));

  // save the pitch
  node->attribs.insert(XMLAttrib("pitch", _pitch));
}

/// Calculates the constraint Jacobian
void ScrewJoint::calc_constraint_jacobian(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Screw joint axis has not been set; set before calling dynamics functions.");

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
  const Vector3d& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3d& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
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

  // setup the constraint equations (from Shabana, p. 101), eq. 3.8 
  if (body == inner)
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
        break;

      case 1:
        break;

      case 2:
        break;

      case 3:
        break;

      case 4:
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
        break;

      case 1:
        break;

      case 2:
        break;

      case 3:
        break;

      case 4:
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
*/
}

/// Calculates the time derivative of the constraint Jacobian
// TODO: implement this properly
void ScrewJoint::calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Screw joint axis has not been set; set before calling dynamics functions.");

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
  const Vector3d& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3d& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
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

  // setup the constraint equations (from Shabana, p. 101), eq. 3.8 
  if (body == inner)
  {
    // now setup the constraint equations
    switch (index)
    {
      case 0:
        break;

      case 1:
        break;

      case 2:
        break;

      case 3:
        break;

      case 4:
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
        break;

      case 1:
        break;

      case 2:
        break;

      case 3:
        break;

      case 4:
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }    
  }
*/
}

/// Evaluates the constraint equations
// TODO: implement this properly
void ScrewJoint::evaluate_constraints(double C[])
{
/*
  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2001], p. 101
  // some variable names have been altered

  // get v1 in global coordinates
  Vector3d v1 = get_axis_global();

  // determine axis in global coordinates
  Vector3d v2 = outer->get_transform().mult_vector(_v2);

  // determine v1i, v1j
  Vector3d v1i, v1j;
  Vector3d::determine_orthonormal_basis(v1, v1i, v1j);

  // determine the global positions of the attachment points and subtract them
  const Vector3d& p1 = get_position_global(false); 
  const Vector3d& p2 = get_position_global(true); 
  Vector3d r12 = p1 - p2; 

  // evaluate the constraint equations
  // 1: h1 x h2 = 0   (joint axis according to body 1 and body2)
  // 2: h2 x s12 = 0  (s12 is a vector defined along the joint axis)
  // 3: tij - pitch*thetaij = 0 (tij/thetaij is relative translation/rotation)
*/
}


