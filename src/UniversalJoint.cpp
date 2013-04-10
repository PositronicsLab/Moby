/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/NullPointerException.h>
#include <Moby/Constants.h>
#include <Moby/AAngle.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/UniversalJoint.h>

using namespace Moby;

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
UniversalJoint::UniversalJoint() : Joint()
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1] = ZEROS_3;
  _u[eAxis2] = ZEROS_3;
  _h2 = ZEROS_3;

  // set rotation to z axis matrix to zero (for debugging purposes)
  _R = Matrix3::zero();

  // set the translation in the transformation to zero
  _T.set_translation(ZEROS_3);

  // setup the spatial axis derivative to zero
  _si_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
  _s0_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
UniversalJoint::UniversalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1] = ZEROS_3;
  _u[eAxis2] = ZEROS_3;
  _h2 = ZEROS_3;

  // set the translation induced by the joint to zero
  _T.set_translation(ZEROS_3);

  // set rotation to z axis matrix to zero (for debugging purposes)
  _R = Matrix3::zero();

  // setup the spatial axis derivative to zero
  _si_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
  _s0_dot = SMatrix6N::zero(SPATIAL_DIM,num_dof());
}  

/// Determines whether two values are relatively equal
bool UniversalJoint::rel_equal(double x, double y)
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
Vector3 UniversalJoint::get_axis_global(Axis a) const
{
  const unsigned X = 0;

  // get the inboard link
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "UniversalJoint::get_axis_global() - attempt to get axis w/o inboard link!" << std::endl;
    return ZEROS_3;
  }

  // get the transform for the inboard link
  Matrix3 R;
  inboard_link->get_transform().get_rotation(&R);

  // axis one is just R times the first column of _R
  if (a == eAxis1)
  {
    Vector3 axis;
    _R.get_column(X, axis.begin());
    return R * axis;
  }

  // axis two is obtained by multiplying rotation matrix around x by y-axis 
  assert(a == eAxis2);
  const double c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);
  return R * _R * Vector3(0,c1,s1);
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
void UniversalJoint::set_axis_local(const Vector3& axis, Axis a) 
{ 
  // normalize the axis in case the user did not
  Vector3 naxis = Vector3::normalize(axis); 

  // set the axis
  _u[a] = naxis; 

  // update the spatial axes
  update_spatial_axes(); 

  // set the second axis in outer link coordinates, if we just changed it
  if (a == eAxis2)
  {
    RigidBodyPtr inboard = get_inboard_link();
    RigidBodyPtr outboard = get_outboard_link();
    Vector3 h2_g = inboard->get_transform().mult_vector(naxis);
    _h2 = outboard->get_transform().transpose_mult_vector(h2_g);
  }
}        

/// Sets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \note set_axis_global() should be called on all axes before setting q! 
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
void UniversalJoint::set_axis_global(const Vector3& axis, Axis a)
{
  // normalize the axis in case the user did not
  Vector3 naxis = Vector3::normalize(axis); 

  // get the inboard link 
  RigidBodyPtr inboard_link = get_inboard_link();

  // make sure that the inboard link has been set
  if (!inboard_link)
  {
    std::cerr << "UniversalJoint::set_axis_global() - attempt to set axis w/o inboard link!" << std::endl;
    return;
  }

  // get the orientation for the inboard link
  Matrix3 R;
  inboard_link->get_transform().get_rotation(&R);

  // transform the axis into inboard link coordinates
  _u[a] = R * naxis;

  // update the spatial axes
  update_spatial_axes();

  // set the second axis in outboard link coordinates, if we just changed it
  if (a == eAxis2)
  {
    RigidBodyPtr outboard = get_outboard_link();
    _h2 = outboard->get_transform().transpose_mult_vector(naxis);
  }
}

/// Updates the spatial axis for this joint
void UniversalJoint::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // call the parent method
  Joint::update_spatial_axes();

  // compute the rotation matrix necessary to transform the first axis to the
  // z-axis
  if (_u[eAxis1].norm() < NEAR_ZERO || _u[eAxis2].norm() < NEAR_ZERO)
    return;

  // first, verify that all are unit-length, and they are orthogonal in seq.
  assert(std::fabs(_u[eAxis1].norm() - 1.0) < NEAR_ZERO);
  assert(std::fabs(_u[eAxis2].norm() - 1.0) < NEAR_ZERO);
  assert(std::fabs(Vector3::dot(_u[eAxis1], _u[eAxis2])) < NEAR_ZERO);
  Vector3 axis3 = Vector3::normalize(Vector3::cross(_u[eAxis1], _u[eAxis2]));
  _R.set_column(X, _u[eAxis1]);
  _R.set_column(Y, _u[eAxis2]);
  _R.set_column(Z, axis3);
}

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const SMatrix6N& UniversalJoint::get_spatial_axes(ReferenceFrameType rftype)
{
  const unsigned X = 0, Y = 1, Z = 2;

  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get current values of q
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;

  // get the outboard link's joint to com vector in link coordinates
  const Vector3& p = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

  // get the axes of the joint transformed into the inner link frame 
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  Vector3 u1;
  _R.get_column(X, u1.begin());
  Vector3 u2 = _R * Vector3(0, c1, s1);

  // update the spatial axes in link coordinates
  SVector6 si1, si2;
  si1.set_upper(u1);
  si2.set_upper(u2);
  si1.set_lower(Vector3::cross(u1, p));
  si2.set_lower(Vector3::cross(u2, p));
  _si.set_column(eAxis1, si1);
  _si.set_column(eAxis2, si2);

  // setup s_bar 
  calc_s_bar_from_si();

  // use the Joint function to do the rest
  return Joint::get_spatial_axes(rftype);
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const SMatrix6N& UniversalJoint::get_spatial_axes_dot(ReferenceFrameType rftype)
{
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get q and qd
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;
  const VectorN& qd = this->qd;

  // form the time derivative of the spatial axis for the second DOF; note that spatial
  // axis for first DOF is constant, so time-derivative is zero 
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double qd1 = qd[DOF_1];
  Vector3 axis(0,-s1*qd1,c1*qd1);

  // get the axis transformed into link coordinates
  Vector3 u = _R * axis;

  // get the outboard link's joint to com vector in link coordinates
  const Vector3& p = outboard->get_inner_joint_data(inboard).joint_to_com_vec_of;

  // update the spatial axis in link coordinates; note that axis 1 is always
  // set to zero (init'd in constructor)
  SVector6 si;
  si.set_upper(u);
  si.set_lower(Vector3::cross(u, p));
  _si_dot.set_column(eAxis2, si);

  // transform to global coordinates
  SpatialTransform X_0_i = outboard->get_spatial_transform_link_to_global();
  X_0_i.transform(_si_dot, _s0_dot);

  return (rftype == eLink) ? _si_dot : _s0_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void UniversalJoint::determine_q(VectorN& q)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // set proper size for q
  this->q.resize(num_dof());

  // get the inboard and outboard links
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();

  // verify that the inboard and outboard links are set
  if (!inboard || !outboard)
    throw std::runtime_error("determine_q() called on NULL inboard and/or outboard links!");

  // if any of the axes are not defined, can't use this method
  if (std::fabs(_u[0].norm() - 1.0) > NEAR_ZERO ||
      std::fabs(_u[1].norm() - 1.0) > NEAR_ZERO)
    throw UndefinedAxisException();

  // get the link transforms
  Matrix3 R_inboard, R_outboard;
  inboard->get_transform().get_rotation(&R_inboard);
  outboard->get_transform().get_rotation(&R_outboard);

  // determine the joint transformation
  Matrix3 R_local = R_inboard.transpose_mult(R_outboard);

  // back out the transformation to z-axis
  Matrix3 RU = _R.transpose_mult(R_local * _R);

  // determine q1 and q2 -- they are uniquely determined by examining the rotation matrix
  // (see get_rotation())
  q.resize(num_dof());
  q[DOF_1] = std::atan2(RU(Z,Y), RU(Y,Y));
  q[DOF_2] = std::atan2(RU(X,Z), RU(X,X));   
}

/// Gets the (local) transform for this joint
Matrix3 UniversalJoint::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q and _q_tare
  const VectorN& q = this->q;
  const VectorN& q_tare = this->_q_tare;

  // compute some needed quantities
  const double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);

  // determine untransformed rotation; this rotation matrix is obtained by
  // using Tait-Bryan angles without a final rotation
  Matrix3 RU;
  RU(X,X) = c2;      RU(X,Y) = 0;    RU(X,Z) = s2;
  RU(Y,X) = s1*s2;  RU(Y,Y) = c1;    RU(Y,Z) = -c2*s1;
  RU(Z,X) = -c1*s2;  RU(Z,Y) = s1;    RU(Z,Z) = c1*c2;

  // transform rotation
  Matrix3 R_local = _R * RU.mult_transpose(_R);

  return R_local;
}

/// Gets the transform induced by this joint
const Matrix4& UniversalJoint::get_transform()
{
  // get the rotation
  Matrix3 R = get_rotation();

  // note that translation is set to zero in the constructors
  _T.set_rotation(&R);

  return _T;
}

/// Computes the constraint jacobian
void UniversalJoint::calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _h2) is set
  if (_u[eAxis1].norm_sq() < std::numeric_limits<double>::epsilon() ||
      _u[eAxis2].norm_sq() < std::numeric_limits<double>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // get the information necessary to compute the constraint equations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
  const double q1x = q1.x;
  const double q1y = q1.y;
  const double q1z = q1.z;
  const double q1w = q1.w;
  const double p1x = p1[X];
  const double p1y = p1[Y];
  const double p1z = p1[Z];
  const double q2x = q2.x;
  const double q2y = q2.y;
  const double q2z = q2.z;
  const double q2w = q2.w;
  const double p2x = -p2[X];
  const double p2y = -p2[Y];
  const double p2z = -p2[Z];
  const double u0x = _u[0][X];
  const double u0y = _u[0][Y];
  const double u0z = _u[0][Z];
  const double h2x = _h2[X];
  const double h2y = _h2[Y];
  const double h2z = _h2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
  {
    switch (index)
    {
      case 0:
        Cq[0] = 1.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*p1x*q1w + 2*p1z*q1y - 2*p1y*q1z;
        Cq[4] = 4*q1x*p1x + 2*q1y*p1y + 2*q1z*p1z;
        Cq[5] = 2*p1z*q1w + 2*p1y*q1x;
        Cq[6] = 2*p1z*q1x - 2*p1y*q1w;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 1.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*p1y*q1w - 2*p1z*q1x + 2*p1x*q1z;
        Cq[4] = 2*q1y*p1x - 2*q1w*p1z;
        Cq[5] = 2*p1x*q1x + 4*p1y*q1y + 2*p1z*q1z;
        Cq[6] = 2*p1x*q1w + 2*p1z*q1y;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 1.0;
        Cq[3] = 4*p1z*q1w + 2*p1y*q1x - 2*p1x*q1y;
        Cq[4] = 2*q1z*p1x + 2*q1w*p1y;
        Cq[5] = 2*p1y*q1z - 2*p1x*q1w;
        Cq[6] = 4*p1z*q1z + 2*p1y*q1y + 2*p1x*q1x;
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                2*(q2x*q2y + q2w*q2z)*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                (-1 + 2*q2w*q2w + 2*q2x*q2x)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z)) + 
                h2y*(2*(q2w*q2x + q2y*q2z)*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                (-1 + 2*(q2w*q2w + q2y*q2y))*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                2*(q2x*q2y - q2w*q2z)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z)) + 
                h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                2*(-(q2w*q2x) + q2y*q2z)*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                2*(q2w*q2y + q2x*q2z)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z));
        Cq[4] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
                (2*q1z*u0x + 2*q1w*u0y) + 
                2*(q2x*q2y + q2w*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
                (-1 + 2*q2w*q2w + 2*q2x*q2x)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
                h2y*(2*(q2w*q2x + q2y*q2z)*(2*q1z*u0x + 2*q1w*u0y) + 
                (-1 + 2*(q2w*q2w + q2y*q2y))*
                (2*q1y*u0x - 2*q1w*u0z) + 
                2*(q2x*q2y - q2w*q2z)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
                h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
                (2*q1z*u0x + 2*q1w*u0y) + 
                2*(-(q2w*q2x) + q2y*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
                2*(q2w*q2y + q2x*q2z)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z));
        Cq[5] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
               (2*q1z*u0x + 2*q1w*u0y) + 
               2*(q2x*q2y + q2w*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
               (-1 + 2*q2w*q2w + 2*q2x*q2x)*
               (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
               h2y*(2*(q2w*q2x + q2y*q2z)*(2*q1z*u0x + 2*q1w*u0y) + 
               (-1 + 2*(q2w*q2w + q2y*q2y))*
               (2*q1y*u0x - 2*q1w*u0z) + 
               2*(q2x*q2y - q2w*q2z)*
               (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
               h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
              (2*q1z*u0x + 2*q1w*u0y) + 
              2*(-(q2w*q2x) + q2y*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
              2*(q2w*q2y + q2x*q2z)*
              (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z));
        Cq[6] = h2x*((-1 + 2*q2w*q2w + 2*q2x*q2x)*
              (-2*q1w*u0y + 2*q1x*u0z) + 
              2*(q2x*q2y + q2w*q2z)*(2*q1w*u0x + 2*q1y*u0z) + 
              2*(-(q2w*q2y) + q2x*q2z)*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z)) + 
              h2y*(2*(q2x*q2y - q2w*q2z)*(-2*q1w*u0y + 2*q1x*u0z) + 
              (-1 + 2*(q2w*q2w + q2y*q2y))*
              (2*q1w*u0x + 2*q1y*u0z) + 
              2*(q2w*q2x + q2y*q2z)*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z)) + 
              h2z*(2*(q2w*q2y + q2x*q2z)*(-2*q1w*u0y + 2*q1x*u0z) + 
              2*(-(q2w*q2x) + q2y*q2z)*(2*q1w*u0x + 2*q1y*u0z) + 
              (-1 + 2*(q2w*q2w + q2z*q2z))*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z));

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
        Cq[3] = -(4*p2x*q2w + 2*p2z*q2y - 2*p2y*q2z);
        Cq[4] = -(4*q2x*p2x + 2*q2y*p2y + 2*q2z*p2z);
        Cq[5] = -(2*p2z*q2w + 2*p2y*q2x);
        Cq[6] = -(2*p2z*q2x - 2*p2y*q2w);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = -1.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*p2y*q2w - 2*p2z*q2x + 2*p2x*q2z);
        Cq[4] = -(2*q2y*p2x - 2*q2w*p2z);
        Cq[5] = -(2*p2x*q2x + 4*p2y*q2y + 2*p2z*q2z);
        Cq[6] = -(2*p2x*q2w + 2*p2z*q2y);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = -1.0;
        Cq[3] = -(4*p2z*q2w + 2*p2y*q2x - 2*p2x*q2y);
        Cq[4] = -(2*q2z*p2x + 2*q2w*p2y);
        Cq[5] = -(2*p2y*q2z - 2*p2x*q2w);
        Cq[6] = -(4*p2z*q2z + 2*p2y*q2y + 2*p2x*q2x);
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = h2z*(2*q2y*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) - 
                2*q2x*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                4*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(-2*q2z*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                4*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2x*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2x*(4*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2z*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) - 
                2*q2y*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[4] = h2z*(2*q2z*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) - 
                2*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z)) + 
                h2y*(2*q2y*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2x*(4*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[5] = h2z*(2*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2z*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z)) + 
                h2x*(2*q2x*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) - 
                2*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(2*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                4*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[6] = h2x*(2*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2x*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(-2*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2z*(2*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                4*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Computes the constraint jacobian
void UniversalJoint::calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // make sure that _u (and by extension _h2) is set
  if (_u[eAxis1].norm_sq() < std::numeric_limits<double>::epsilon() ||
      _u[eAxis2].norm_sq() < std::numeric_limits<double>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (double) 0.0;
    return;
  }

  // get the information necessary to compute the constraint equations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Quat qd1 = Quat::deriv(q1, inner->get_avel());
  const Quat qd2 = Quat::deriv(q2, outer->get_avel());
  const Vector3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
  const double qx1 = q1.x;
  const double qy1 = q1.y;
  const double qz1 = q1.z;
  const double qw1 = q1.w;
  const double p1x = p1[X];
  const double p1y = p1[Y];
  const double p1z = p1[Z];
  const double qx2 = q2.x;
  const double qy2 = q2.y;
  const double qz2 = q2.z;
  const double qw2 = q2.w;
  const double p2x = -p2[X];
  const double p2y = -p2[Y];
  const double p2z = -p2[Z];
  const double ux = _u[0][X];
  const double uy = _u[0][Y];
  const double uz = _u[0][Z];
  const double h2x = _h2[X];
  const double h2y = _h2[Y];
  const double h2z = _h2[Z];
  const double dqw1 = qd1.w;
  const double dqx1 = qd1.x;
  const double dqy1 = qd1.y;
  const double dqz1 = qd1.z;
  const double dqw2 = qd2.w;
  const double dqx2 = qd2.x;
  const double dqy2 = qd2.y;
  const double dqz2 = qd2.z;

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
  {
    switch (index)
    {
      case 0:
        Cq[0] = (double) 0.0;    
        Cq[1] = (double) 0.0;    
        Cq[2] = (double) 0.0;    
        Cq[3] = 4*p1x*dqw1 + 2*p1z*dqy1 - 2*p1y*dqz1;
        Cq[4] = 4*dqx1*p1x + 2*dqy1*p1y + 2*dqz1*p1z;
        Cq[5] = 2*p1z*dqw1 + 2*p1y*dqx1;
        Cq[6] = 2*p1z*dqx1 - 2*p1y*dqw1;
        break;

      case 1:
        Cq[0] = (double) 0.0;    
        Cq[1] = (double) 0.0;    
        Cq[2] = (double) 0.0;    
        Cq[3] = 4*p1y*dqw1 - 2*p1z*dqx1 + 2*p1x*dqz1;
        Cq[4] = 2*dqy1*p1x - 2*dqw1*p1z;
        Cq[5] = 2*p1x*dqx1 + 4*p1y*dqy1 + 2*p1z*dqz1;
        Cq[6] = 2*p1x*dqw1 + 2*p1z*dqy1;
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = 4*p1z*dqw1 + 2*p1y*dqx1 - 2*p1x*dqy1;
        Cq[4] = 2*dqz1*p1x + 2*dqw1*p1y;
        Cq[5] = 2*p1y*dqz1 - 2*p1x*dqw1;
        Cq[6] = 4*p1z*dqz1 + 2*p1y*dqy1 + 2*p1x*dqx1;
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (-2*dqy1*ux + 2*dqx1*uy + 4*dqw1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqz1*ux + 4*dqw1*uy - 2*dqx1*uz) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (4*dqw1*ux - 2*dqz1*uy + 2*dqy1*uz) + 
   (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qy1*ux + 2*qx1*uy + 4*qw1*uz)\
    + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qz1*ux + 4*qw1*uy - 2*qx1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qw1*ux - 2*qz1*uy + 2*qy1*uz);
        Cq[4] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*(2*dqz1*ux + 2*dqw1*uy)
     + (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qz1*ux + 2*qw1*uy) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqy1*ux - 2*dqw1*uz) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (4*dqx1*ux + 2*dqy1*uy + 2*dqz1*uz) + 
   (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qy1*ux - 2*qw1*uz) + (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qx1*ux + 2*qy1*uy + 2*qz1*uz);
        Cq[5] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (-2*dqw1*ux + 2*dqz1*uy) + 
   (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qw1*ux + 2*qz1*uy) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (2*dqx1*uy + 2*dqw1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqx1*ux + 4*dqy1*uy + 2*dqz1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (2*qx1*uy + 2*qw1*uz) + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qx1*ux + 4*qy1*uy + 2*qz1*uz);
        Cq[6] = (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (-2*dqw1*uy + 2*dqx1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqw1*ux + 2*dqy1*uz) + 
   (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (2*dqx1*ux + 2*dqy1*uy + 4*dqz1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (-2*qw1*uy + 2*qx1*uz) + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qw1*ux + 2*qy1*uz) + (2*h2x*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qx1*ux + 2*qy1*uy + 4*qz1*uz);

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
        Cq[3] = -(4*p2x*dqw2 + 2*p2z*dqy2 - 2*p2y*dqz2);
        Cq[4] = -(4*dqx2*p2x + 2*dqy2*p2y + 2*dqz2*p2z);
        Cq[5] = -(2*p2z*dqw2 + 2*p2y*dqx2);
        Cq[6] = -(2*p2z*dqx2 - 2*p2y*dqw2);
        break;

      case 1:
        Cq[0] = (double) 0.0;      
        Cq[1] = (double) 0.0;     
        Cq[2] = (double) 0.0;      
        Cq[3] = -(4*p2y*dqw2 - 2*p2z*dqx2 + 2*p2x*dqz2);
        Cq[4] = -(2*dqy2*p2x - 2*dqw2*p2z);
        Cq[5] = -(2*p2x*dqx2 + 4*p2y*dqy2 + 2*p2z*dqz2);
        Cq[6] = -(2*p2x*dqw2 + 2*p2z*dqy2);
        break;

      case 2:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = -(4*p2z*dqw2 + 2*p2y*dqx2 - 2*p2x*dqy2);
        Cq[4] = -(2*dqz2*p2x + 2*dqw2*p2y);
        Cq[5] = -(2*p2y*dqz2 - 2*p2x*dqw2);
        Cq[6] = -(4*p2z*dqz2 + 2*p2y*dqy2 + 2*p2x*dqx2);
        break;

      case 3:
        Cq[0] = (double) 0.0;
        Cq[1] = (double) 0.0;
        Cq[2] = (double) 0.0;
        Cq[3] = (4*h2x*qw2 + 2*h2z*qy2 - 2*h2y*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (4*h2y*qw2 - 2*h2z*qx2 + 2*h2x*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (4*h2z*qw2 + 2*h2y*qx2 - 2*h2x*qy2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (4*dqw2*h2x - 2*dqz2*h2y + 2*dqy2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqz2*h2x + 4*dqw2*h2y - 2*dqx2*h2z)*
    (2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (-2*dqy2*h2x + 2*dqx2*h2y + 4*dqw2*h2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[4] = (4*h2x*qx2 + 2*h2y*qy2 + 2*h2z*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (-2*h2z*qw2 + 2*h2x*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (2*h2y*qw2 + 2*h2x*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (4*dqx2*h2x + 2*dqy2*h2y + 2*dqz2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqy2*h2x - 2*dqw2*h2z)*(2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (2*dqz2*h2x + 2*dqw2*h2y)*(2*(-(qw1*qy1) + qx1*qz1)*ux + 
      2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[5] = (2*h2z*qw2 + 2*h2y*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (2*h2x*qx2 + 4*h2y*qy2 + 2*h2z*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (-2*h2x*qw2 + 2*h2y*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (2*dqx2*h2y + 2*dqw2*h2z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       ux + 2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqx2*h2x + 4*dqy2*h2y + 2*dqz2*h2z)*
    (2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (-2*dqw2*h2x + 2*dqz2*h2y)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[6] = (-2*h2y*qw2 + 2*h2z*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (2*h2x*qw2 + 2*h2z*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (2*h2x*qx2 + 2*h2y*qy2 + 4*h2z*qz2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (-2*dqw2*h2y + 2*dqx2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqw2*h2x + 2*dqy2*h2z)*(2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (2*dqx2*h2x + 2*dqy2*h2y + 4*dqz2*h2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Evaluates the constraint equations
void UniversalJoint::evaluate_constraints(double C[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], p. 438; variable names
  // have been altered however

  // determine h1 and h2 in global coordinates
  Vector3 h1 = inner->get_transform().mult_vector(_u[0]);
  Vector3 h2 = outer->get_transform().mult_vector(_h2);

  // determine the global positions of the attachment points and subtract them
  Vector3 r1 = get_position_global(false);
  Vector3 r2 = get_position_global(true);
  Vector3 r12 = r1 - r2;

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = h1.dot(h2);
}

/// Implements Base::load_from_xml()
void UniversalJoint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "UniversalJoint") == 0);

  // read local joint axes
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

  // determine _q_tare if necessary 
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void UniversalJoint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "UniversalJoint";

  // save the local joint axes
  node->attribs.insert(XMLAttrib("local-axis1", _u[eAxis1]));
  node->attribs.insert(XMLAttrib("local-axis2", _u[eAxis1]));
}

