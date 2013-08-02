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
#include <Moby/UndefinedAxisException.h>
#include <Moby/UniversalJoint.h>

using boost::shared_ptr;
using std::vector;
using namespace Ravelin;
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
  _u[eAxis1].set_zero();
  _u[eAxis2].set_zero();
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _h2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;
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
  _u[eAxis1].set_zero();
  _u[eAxis2].set_zero();
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _h2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;
}  

/// Determines whether two values are relatively equal
bool UniversalJoint::rel_equal(double x, double y)
{
  return (std::fabs(x - y) <= NEAR_ZERO * std::max(std::fabs(x), std::max(std::fabs(y), (double) 1.0)));
}

/// Gets the axis for this joint
Vector3d UniversalJoint::get_axis(Axis a) const
{
  const unsigned X = 0;

  // axis one is already set 
  if (a == eAxis1)
  {
    return _u[0];
  }

  // axis two is obtained by multiplying rotation matrix around x by y-axis 
  assert(a == eAxis2);
  const double c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);
  Vector3d u(0,c1,s1);
  u.pose = get_pose();
  return u;
}

/// Sets an axis of this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void UniversalJoint::set_axis(const Vector3d& axis, Axis a) 
{ 
  // normalize the axis in case the user did not
  Vector3d naxis = Vector3d::normalize(axis); 

  // set the axis
  _u[a] = Pose3d::transform(get_pose(), naxis); 

  // update the spatial axes
  update_spatial_axes(); 

  // set the second axis in outer link coordinates, if we just changed it
  if (a == eAxis2)
  {
    RigidBodyPtr inboard = get_inboard_link();
    RigidBodyPtr outboard = get_outboard_link();
// TODO: re-enable this with fix to joint equations
/*
    Vector3d h2_g = inboard->get_transform().mult_vector(naxis);
    _h2 = outboard->get_transform().transpose_mult_vector(h2_g);
*/
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
  assert(std::fabs(Vector3d::dot(_u[eAxis1], _u[eAxis2])) < NEAR_ZERO);
  Vector3d axis3 = Vector3d::normalize(Vector3d::cross(_u[eAxis1], _u[eAxis2]));
}

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SAxisd>& UniversalJoint::get_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Vector3d ZEROS_3(0.0, 0.0, 0.0, get_pose());

  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get current values of q
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;

  // get the axes of the joint transformed into the inner link frame 
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  Vector3d u1 = _u[0];
  Vector3d u2(0, c1, s1);
  u2.pose = u1.pose;

  // update the spatial axes in link coordinates
  _s[0].set_angular(u1);
  _s[0].set_linear(ZEROS_3);
  _s[1].set_angular(u2);
  _s[1].set_linear(ZEROS_3);

  // setup s_bar 
  calc_s_bar_from_s();

  // use the Joint function to do the rest
  return Joint::get_spatial_axes();
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SAxisd>& UniversalJoint::get_spatial_axes_dot()
{
  const Vector3d ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get the inboard and outboard links
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  if (!inboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UniversalJoint::get_spatial_axes_dot() called with NULL outboard link");

  // get q and qd
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;
  const VectorNd& qd = this->qd;

  // form the time derivative of the spatial axis for the second DOF; note that spatial
  // axis for first DOF is constant, so time-derivative is zero 
  double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  double qd1 = qd[DOF_1];
  Vector3d u(0,-s1*qd1,c1*qd1);
  u.pose = get_pose();

  // update the spatial axis in link coordinates; note that axis 1 is always
  // set to zero (init'd in constructor)
  _s_dot[1].set_upper(u);
  _s_dot[1].set_lower(ZEROS_3);

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void UniversalJoint::determine_q(VectorNd& q)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the outboard link
  RigidBodyPtr outboard = get_outboard_link();

  // verify that the outboard link is set
  if (!outboard)
    throw std::runtime_error("determine_q() called on NULL outboard link!");

  // set proper size for q
  this->q.resize(num_dof());

  // get the poses of the joint and outboard link
  shared_ptr<const Pose3d> Fj = get_pose();
  shared_ptr<const Pose3d> Fo = outboard->get_pose();

  // compute transforms
  Transform3d oT0 = Pose3d::calc_relative_pose(GLOBAL, Fo); 
  Transform3d jT0 = Pose3d::calc_relative_pose(GLOBAL, Fj);
  Transform3d oTj = oT0 * jT0.inverse();

  // determine the joint transformation
  Matrix3d R = oTj.q;

  // determine q1 and q2 -- they are uniquely determined by examining the rotation matrix
  // (see get_rotation())
  q.resize(num_dof());
  q[DOF_1] = std::atan2(R(Z,Y), R(Y,Y));
  q[DOF_2] = std::atan2(R(X,Z), R(X,X));   
}

/// Gets the (local) transform for this joint
Matrix3d UniversalJoint::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q and _q_tare
  const VectorNd& q = this->q;
  const VectorNd& q_tare = this->_q_tare;

  // compute some needed quantities
  const double c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const double s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  const double s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);

  // determine untransformed rotation; this rotation matrix is obtained by
  // using Tait-Bryan angles without a final rotation
  Matrix3d R;
  R(X,X) = c2;      R(X,Y) = 0;    R(X,Z) = s2;
  R(Y,X) = s1*s2;  R(Y,Y) = c1;    R(Y,Z) = -c2*s1;
  R(Z,X) = -c1*s2;  R(Z,Y) = s1;   R(Z,Z) = c1*c2;

  return R;
}

/// Gets the transform induced by this joint
shared_ptr<const Pose3d> UniversalJoint::get_induced_pose()
{
  // invalidate pose quantities for the outer link
  invalidate_pose_vectors();

  // get the rotation
  _Fprime->q = get_rotation();
  return _Fprime;
}

/// Computes the constraint jacobian
void UniversalJoint::calc_constraint_jacobian(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
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
  const Vector3d& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3d& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
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
*/
}

/// Computes the constraint jacobian
void UniversalJoint::calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, double Cq[7])
{
/*
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
  const Vector3d& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const Vector3d& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
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
*/
}

/// Evaluates the constraint equations
void UniversalJoint::evaluate_constraints(double C[])
{
/*
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], p. 438; variable names
  // have been altered however

  // determine h1 and h2 in global coordinates
  Vector3d h1 = inner->get_transform().mult_vector(_u[0]);
  Vector3d h2 = outer->get_transform().mult_vector(_h2);

  // determine the global positions of the attachment points and subtract them
  Vector3d r1 = get_position_global(false);
  Vector3d r2 = get_position_global(true);
  Vector3d r12 = r1 - r2;

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = h1.dot(h2);
*/
}

/// Implements Base::load_from_xml()
void UniversalJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "UniversalJoint") == 0);

  // read the global joint axes, if given
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

  // determine _q_tare if necessary 
  if (_determine_q_tare)
    determine_q_tare();
}

/// Implements Base::save_to_xml()
void UniversalJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  Vector3d u0;

  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "UniversalJoint";

  // convert local axes to global axes and save
  u0 = Pose3d::transform(shared_ptr<const Pose3d>(), _u[eAxis1]);
  node->attribs.insert(XMLAttrib("axis1", u0));
  u0 = Pose3d::transform(shared_ptr<const Pose3d>(), _u[eAxis2]);
  node->attribs.insert(XMLAttrib("axis2", u0));
  node->attribs.insert(XMLAttrib("axis1", _u[eAxis1]));
  node->attribs.insert(XMLAttrib("axis2", _u[eAxis1]));
}

