/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/ScrewJoint.h>

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
  _pitch = (Real) 1.0;

  // init the joint axes
  _u = ZEROS_3;
  _v2 = ZEROS_3;

  // set the transformation to identity 
  _T = IDENTITY_4x4;;

  // setup the spatial axis derivative to zero
  _si_deriv = SMatrix6N(1);
  _si_deriv.set_column(0, SVector6(0,0,0,0,0,0));
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
  _pitch = (Real) 1.0;

  // init the joint axes
  _u = ZEROS_3;
  _v2 = ZEROS_3;

  // set the transformation to identity 
  _T = IDENTITY_4x4;

  // setup the spatial axis derivative to zero
  _si_deriv = SMatrix6N(1);
  _si_deriv.set_column(0, SVector6(0,0,0,0,0,0));
}  

/// Gets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link 
 * into account; thus, if the orientation of the inboard link changes, then 
 * the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
Vector3 ScrewJoint::get_axis_global() const
{
  // get the inboard link 
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "ScrewJoint::get_axis_global() - attempt to get axis w/o inboard link!" << std::endl;
    return ZEROS_3;
  }

  // get the transform for the inboard link
  const Matrix4& T = inboard_link->get_transform();
  
  // transform into global coordinates and return
  return T.transpose_mult_vector(_u);
}

/// Sets the local axis for this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void ScrewJoint::set_axis_local(const Vector3& axis) 
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

  // set joint axis in outer link frame
  _v2 = inner->get_transform().mult_vector(naxis);
  _v2 = outer->get_transform().transpose_mult_vector(_v2);
}        

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link 
 * into account; thus, if the orientation of the inboard link changes, then 
 * the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void ScrewJoint::set_axis_global(const Vector3& axis)
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
    const Vector3& di = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

    // get the joint axis in outer link coordinates
    Vector3 u0 = inboard->get_transform().mult_vector(_u);
    Vector3 ui = outboard->get_transform().transpose_mult_vector(u0);

    // update the spatial axis in link coordinates
    Vector3 x = Vector3::cross(ui, di);
    SVector6 si_vec;
    si_vec.set_upper(ui);
    si_vec.set_lower(ui*_pitch + x);
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
  Vector3 p1 = get_position_global(false);
  Vector3 p2 = get_position_global(true);

  // get the joint axis in the global frame
  Vector3 ug = inboard->get_transform().mult_vector(_u);

  // now, we'll project p2 onto the axis ug; points will be setup so that
  // ug passes through origin on inboard
  q.resize(num_dof());
  q[DOF_1] = ug.dot(p2-p1)/_pitch;
}

/// Gets the (local) transform for this joint
const Matrix4& ScrewJoint::get_transform()
{
  // setup rotation 
  AAngle a(&_u, this->q[DOF_1]+this->_q_tare[DOF_1]);
  _T.set_rotation(&a);

  // setup translation
  _T.set_translation(_u*(this->q[DOF_1]+this->_q_tare[DOF_1]));

  return _T;
}

/// Gets the derivative for the spatial axes for this joint
const SMatrix6N& ScrewJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  return _si_deriv;
}

/// Implements Base::load_from_xml()
void ScrewJoint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "ScrewJoint") == 0);

  // read the local joint axis
  const XMLAttrib* laxis_attrib = node->get_attrib("local-axis");
  if (laxis_attrib)
  {
    Vector3 laxis;
    laxis_attrib->get_vector_value(laxis);
    set_axis_local(laxis);
  }

  // read the pitch
  const XMLAttrib* pitch_attrib = node->get_attrib("pitch");
  if (pitch_attrib)
    _pitch = pitch_attrib->get_real_value(); 

  // read the global joint axis, if given
  const XMLAttrib* gaxis_attrib = node->get_attrib("global-axis");
  if (gaxis_attrib)
  {
    Vector3 gaxis;
    gaxis_attrib->get_vector_value(gaxis);
    set_axis_global(gaxis);  
  }

  // compute _q_tare if necessary 
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void ScrewJoint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "ScrewJoint";

  // save the local joint axis
  node->attribs.insert(XMLAttrib("local-axis", _u));

  // save the pitch
  node->attribs.insert(XMLAttrib("pitch", _pitch));
}

/// Calculates the constraint Jacobian
// TODO: implement this properly
void ScrewJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<Real>::epsilon())
    throw std::runtime_error("Screw joint axis has not been set; set before calling dynamics functions.");

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
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
  const Real x1 = inner->get_position()[X];
  const Real y1 = inner->get_position()[Y];
  const Real z1 = inner->get_position()[Z];
  const Real x2 = outer->get_position()[X];
  const Real y2 = outer->get_position()[Y];
  const Real z2 = outer->get_position()[Z];
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
}

/// Calculates the time derivative of the constraint Jacobian
// TODO: implement this properly
void ScrewJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _ui, _uj, _v2) is set
  if (_u.norm_sq() < std::numeric_limits<Real>::epsilon())
    throw std::runtime_error("Screw joint axis has not been set; set before calling dynamics functions.");

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
  const Vector3& p2 = outer->get_inner_joint_data(inner).joint_to_com_vec_of;
  const Real x1 = inner->get_position()[X];
  const Real y1 = inner->get_position()[Y];
  const Real z1 = inner->get_position()[Z];
  const Real x2 = outer->get_position()[X];
  const Real y2 = outer->get_position()[Y];
  const Real z2 = outer->get_position()[Z];
  const Real dx1 = inner->get_lvel()[X];
  const Real dy1 = inner->get_lvel()[Y];
  const Real dz1 = inner->get_lvel()[Z];
  const Real dx2 = outer->get_lvel()[X];
  const Real dy2 = outer->get_lvel()[Y];
  const Real dz2 = outer->get_lvel()[Z];
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
  const Real dqw1 = qd1.w;
  const Real dqx1 = qd1.x;
  const Real dqy1 = qd1.y;
  const Real dqz1 = qd1.z;
  const Real dqw2 = qd2.w;
  const Real dqx2 = qd2.x;
  const Real dqy2 = qd2.y;
  const Real dqz2 = qd2.z;
  const Real uix = _ui[X];
  const Real uiy = _ui[Y];
  const Real uiz = _ui[Z];
  const Real ujx = _uj[X];
  const Real ujy = _uj[Y];
  const Real ujz = _uj[Z];
  const Real v2x = _v2[X];
  const Real v2y = _v2[Y];
  const Real v2z = _v2[Z];

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
}

/// Evaluates the constraint equations
// TODO: implement this properly
void ScrewJoint::evaluate_constraints(Real C[])
{
  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2001], p. 101
  // some variable names have been altered

  // get v1 in global coordinates
  Vector3 v1 = get_axis_global();

  // determine axis in global coordinates
  Vector3 v2 = outer->get_transform().mult_vector(_v2);

  // determine v1i, v1j
  Vector3 v1i, v1j;
  Vector3::determine_orthonormal_basis(v1, v1i, v1j);

  // determine the global positions of the attachment points and subtract them
  const Vector3& p1 = get_position_global(false); 
  const Vector3& p2 = get_position_global(true); 
  Vector3 r12 = p1 - p2; 

  // evaluate the constraint equations
  // 1: h1 x h2 = 0   (joint axis according to body 1 and body2)
  // 2: h2 x s12 = 0  (s12 is a vector defined along the joint axis)
  // 3: tij - pitch*thetaij = 0 (tij/thetaij is relative translation/rotation)

}


