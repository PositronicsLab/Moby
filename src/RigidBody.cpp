/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <queue>
#include <iostream>
#include <iomanip>
#include <limits>
#include <Moby/CollisionGeometry.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/XMLTree.h>
#include <Moby/Joint.h>
#include <Moby/Log.h>
#include <Moby/Spatial.h>
#include <Moby/RigidBody.h>
#ifdef USE_OSG
#include "Color.h"
#endif

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::cerr;
using std::endl;
using std::map;
using std::list;
using std::queue;

/// Default constructor
/**
 * Constructs a rigid body with zero mass, zero inertia tensor, and center
 * of mass at [0,0,0] with position at [0,0,0], identity orientation, and zero
 * linear and angular velocity.  Body is enabled by default.
 */
RigidBody::RigidBody()
{
  const unsigned SPATIAL_DIM = 6;

  // setup reference pose
  _F = shared_ptr<Pose3d>(new Pose3d(Pose3d::identity()));
  _Ji.pose = _xdi.pose = _xddi.pose = _forcei.pose = _F;

  // setup inertial pose
  _jF = shared_ptr<Pose3d>(new Pose3d(Pose3d::identity()));
  _jF->rpose = _F;
  _Jm.pose = _xdm.pose = _xddm.pose = _forcem.pose = _jF;

  // setup visualization pose
  _vF->rpose = _F;

  // setup c.o.m. frame link pose
  _F2 = shared_ptr<Pose3d>(new Pose3d);
  _Jcom.pose = _xdcom.pose = _xddcom.pose = _forcecom.pose = _F2;

  // setup poses for acceleration bounds
  _vel_limit_lo.set_zero(_F);
  _vel_limit_hi.set_zero(_F);

  // invalidate everything
  _forcei_valid = false;
  _forcej_valid = false;
  _forcem_valid = false;
  _force0_valid = false;
  _xdi_valid = false;
  _xdj_valid = false;
  _xdm_valid = false;
  _xd0_valid = false;
  _xddi_valid = false;
  _xddj_valid = false;
  _xddm_valid = false;
  _xdd0_valid = false;
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // use link c.o.m. frame by default
  _rftype = eLinkCOM;

  // set everything else
  _enabled = true;
  _link_idx = std::numeric_limits<unsigned>::max();
  viscous_coeff = VectorNd::zero(SPATIAL_DIM);

  // indicate velocity limit has been exceeded (safe initialization)
 _vel_limit_exceeded = true;

  // setup the default limit bound expansion
  limit_bound_expansion = 0.15;
  compliance = eRigid;
}

/// Resets the acceleration limit estimates
void RigidBody::reset_limit_estimates()
{
  const unsigned SPATIAL_DIM = 6;

  // mark velocity limits as not exceeded
 _vel_limit_exceeded = false;

  SVelocityd v = Pose3d::transform(_F, get_velocity());
  for (unsigned i=0; i< SPATIAL_DIM; i++)
  {
    _vel_limit_lo[i] = v[i];
    _vel_limit_hi[i] = v[i];
    if (v[i] < 0.0)
    {
      _vel_limit_lo[i] *= (1.0+limit_bound_expansion);
      _vel_limit_hi[i] *= (1.0-limit_bound_expansion);
    }
    else
    {
      _vel_limit_lo[i] *= (1.0-limit_bound_expansion);
      _vel_limit_hi[i] *= (1.0+limit_bound_expansion);
    }
  }
}

/// Updates this rigid body's velocity limit
void RigidBody::update_vel_limits()
{
  const unsigned SPATIAL_DIM = 6;

  // mark limit as not exceeded
  _vel_limit_exceeded = false;

  SVelocityd v = Pose3d::transform(_F, get_velocity());
  for (unsigned i=0; i< SPATIAL_DIM; i++)
  {
    if (v[i] < _vel_limit_lo[i])
    {
      _vel_limit_lo[i] = v[i];
      if (v[i] < 0.0)
        _vel_limit_lo[i] *= (1.0+limit_bound_expansion);
      else
        _vel_limit_lo[i] *= (1.0-limit_bound_expansion);
    }
    if (v[i] > _vel_limit_hi[i])
    {
      _vel_limit_hi[i] = v[i];
      if (v[i] < 0.0)
        _vel_limit_hi[i] *= (1.0-limit_bound_expansion);
      else
        _vel_limit_hi[i] *= (1.0+limit_bound_expansion);
    }

    assert(_vel_limit_lo[i] <= _vel_limit_hi[i]);
  }
}

/// Checks whether a rigid body has exceeded its velocity limit and updates the velocity limit
void RigidBody::check_vel_limit_exceeded_and_update()
{
  const unsigned SPATIAL_DIM = 6;

  SVelocityd v = Pose3d::transform(_F, get_velocity());
  for (unsigned i=0; i< SPATIAL_DIM; i++)
  {
    if (v[i] < _vel_limit_lo[i])
    {
      _vel_limit_lo[i] = v[i];
       if (v[i] < 0.0)
         _vel_limit_lo[i] *= (1.0+limit_bound_expansion);
       else
         _vel_limit_lo[i] *= (1.0-limit_bound_expansion);
      _vel_limit_exceeded = true;
    }
    if (v[i] > _vel_limit_hi[i])
    {
      _vel_limit_hi[i] = v[i];
       if (v[i] < 0.0)
         _vel_limit_hi[i] *= (1.0-limit_bound_expansion);
       else
         _vel_limit_hi[i] *= (1.0+limit_bound_expansion);
      _vel_limit_exceeded = true;
    }
    assert(_vel_limit_lo[i] <= _vel_limit_hi[i]);
  }
}

/// Gets the frame in which kinematics and dynamics computations occur
shared_ptr<const Pose3d> RigidBody::get_computation_frame() const
{
  switch (_rftype)
  {
    case eLink:
      return _F;

    case eLinkInertia:
      return _jF;

    case eLinkCOM:
      return _F2;

    case eGlobal:
      return shared_ptr<const Pose3d>();

    case eJoint:
      return (_abody.expired() || is_base()) ? _F : get_inner_joint_explicit()->get_pose();

    default:
      assert(false);
  }

  return shared_ptr<const Pose3d>();
}

/// Rotates the rigid body
void RigidBody::rotate(const Quatd& q)
{
  // save the current relative pose
  shared_ptr<const Pose3d> Frel = _F->rpose;

  // update the rotation
  _F->update_relative_pose(GLOBAL);
  _F->q *= q;
  _F->update_relative_pose(Frel);

  // invalidate vector quantities
  _forcei_valid = _forcem_valid = _force0_valid = false;
  _xdi_valid = _xdm_valid = _xd0_valid = false;
  _xddi_valid = _xddm_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<RigidBodyPtr> outer;
  BOOST_FOREACH(JointPtr j, _outer_joints)
  {
    if (j->get_constraint_type() == Joint::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<RigidBodyPtr>::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<RigidBodyPtr>::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// Computes the Jacobian
MatrixNd& RigidBody::calc_jacobian_dot(shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  if (body != shared_from_this())
    throw std::runtime_error("RigidBody::calc_jacobian_dot() called with wrong body!");

  if (!is_enabled())
  {
    J.set_zero(SPATIAL_DIM, 0);
    return J;
  }

  J.set_zero(SPATIAL_DIM, SPATIAL_DIM);

  return J;
}

/// Computes the Jacobian
MatrixNd& RigidBody::calc_jacobian(shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  if (body != shared_from_this())
    throw std::runtime_error("RigidBody::calc_jacobian() called with wrong body!");

  // if the body is disabled, do not compute a Jacobian
  if (!is_enabled())
  {
    J.set_zero(6, 0);
    return J;
  }

  // construct the spatial transform
  Pose3d::spatial_transform_to_matrix2(_F2, frame, J);

  FILE_LOG(LOG_DYNAMICS) << "RigidBody::calc_jacobian() entered" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  pose: " << ((frame) ? Pose3d(*frame).update_relative_pose(GLOBAL) : GLOBAL) << std::endl;

  return J;
}

/// Translates the rigid body
void RigidBody::translate(const Origin3d& x)
{
  // save the current relative pose
  shared_ptr<const Pose3d> Frel = _F->rpose;

  // update the translation
  _F->update_relative_pose(GLOBAL);
  _F->x += x;
  _F->update_relative_pose(Frel);

  // update the mixed pose
  _F2->set_identity();
  _F2->rpose = _F;
  _F2->update_relative_pose(GLOBAL);
  _F2->q.set_identity();

  // invalidate vector quantities
  _forcei_valid = _forcem_valid = _force0_valid = false;
  _xdi_valid = _xdm_valid = _xd0_valid = false;
  _xddi_valid = _xddm_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<RigidBodyPtr> outer;
  BOOST_FOREACH(JointPtr j, _outer_joints)
  {
    if (j->get_constraint_type() == Joint::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<RigidBodyPtr>::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<RigidBodyPtr>::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// (Re)sets the computation frame
void RigidBody::set_computation_frame_type(ReferenceFrameType rftype)
{
  // correct rftype if necessary
  if (_abody.expired() && rftype == eJoint)
    rftype = eLink;

  // store the new reference frame type
  _rftype = rftype;
}

/// Computes the forward dynamics for this body
void RigidBody::calc_fwd_dyn()
{
  // if the body is free, just compute linear and angular acceleration via
  // Newton's and Euler's laws
  if (_abody.expired())
  {
    // don't do anything if the body is enabled
    if (!is_enabled())
      return;

    // make sure that the inertia is reasonable
    const SpatialRBInertiad& J = get_inertia();
    #ifndef NDEBUG
    if (J.m <= 0.0 || J.J.norm_inf() <= 0.0)
      throw std::runtime_error("Tried to calculate forward dynamics on body with zero mass/inertia");
    #endif

    // otherwise, calculate forward dynamics
    SForced f = sum_forces() - calc_euler_torques();
    SAcceld xdd = J.inverse_mult(f);

    // set the acceleration
    switch (_rftype)
    {
      case eGlobal:
        _xdd0 = xdd;
        _xddcom = Pose3d::transform(_F2, xdd);
        _xddi_valid = _xddj_valid = _xddm_valid = false;
        break;

      case eLinkCOM:
        _xddcom = xdd;
        _xddi_valid = _xddj_valid = _xddm_valid = _xdd0_valid = false;
        break;

      case eLink:
        _xddi = xdd;
        _xddcom = Pose3d::transform(_F2, xdd);
        _xddi_valid = true;
        _xddj_valid = _xddm_valid = _xdd0_valid = false;
        break;

      case eLinkInertia:
        _xddm = xdd;
        _xddcom = Pose3d::transform(_F2, xdd);
        _xddm_valid = true;
        _xddi_valid = _xddj_valid = _xdd0_valid = false;
        break;

      case eJoint:
        _xddj = xdd;
        _xddcom = Pose3d::transform(_F2, xdd);
        _xddj_valid = true;
        _xddi_valid = _xddm_valid = _xdd0_valid = false;
        break;

      default:
        assert(false);
    }
  }
  else
  {
    // otherwise, need to call forward dynamics on the articulated body
    ArticulatedBodyPtr abody(_abody);

    // calculate forward dynamics on it
    abody->calc_fwd_dyn();
  }
}

/// Sets the body to enabled / disabled.
/**
 * If the body is disabled, the linear and angular velocity are set to zero,
 * and the body will not be updated if it is attempted to integrate its
 * equations of motion.
 */
void RigidBody::set_enabled(bool flag)
{
  // mark as enabled / disabled
  _enabled = flag;

  // if disabled, then zero the velocities and accelerations
  if (!_enabled)
  {
    _xdcom.set_zero();
    _xddcom.set_zero();
    _xd0.set_zero();
    _xdd0.set_zero();
    _xdi.set_zero();
    _xddi.set_zero();
    _xdj.set_zero();
    _xddj.set_zero();
    _xdm.set_zero();
    _xddm.set_zero();
    _xdi_valid = _xdj_valid = _xdm_valid = _xd0_valid = true;
    _xddi_valid = _xddj_valid = _xddm_valid = _xdd0_valid = true;
  }
}

/// Sets the velocity of this body
void RigidBody::set_velocity(const SVelocityd& xd)
{
  // set the velocity
  _xdcom = Pose3d::transform(_F2, xd);

  // invalidate the remaining velocities
  _xdi_valid = _xdj_valid = _xdm_valid = _xd0_valid = false;

  // see whether we can re-validate a velocity
  if (xd.pose == _F)
  {
    _xdi_valid = true;
    _xdi = xd;
  }
  else if (xd.pose == _jF)
  {
    _xdm_valid = true;
    _xdm = xd;
  }
  else if (!is_base() && xd.pose == get_inner_joint_explicit()->get_pose())
  {
    _xdj_valid = true;
    _xdj = xd;
  }
  else if (xd.pose == GLOBAL)
  {
    _xd0_valid = true;
    _xd0 = xd;
  }
}

/// Sets the acceleration of this body
void RigidBody::set_accel(const SAcceld& xdd)
{
  // set the acceleration
  _xddcom = Pose3d::transform(_F2, xdd);

  // invalidate the remaining accelerations
  _xddi_valid = _xddj_valid = _xddm_valid = _xdd0_valid = false;

  // see whether we can re-validate an acceleration
  if (xdd.pose == _F)
  {
    _xddi_valid = true;
    _xddi = xdd;
  }
  else if (xdd.pose == _jF)
  {
    _xddm_valid = true;
    _xddm = xdd;
  }
  else if (!is_base() && xdd.pose == get_inner_joint_explicit()->get_pose())
  {
    _xddj_valid = true;
    _xddj = xdd;
  }
  else if (xdd.pose == GLOBAL)
  {
    _xdd0_valid = true;
    _xdd0 = xdd;
  }
}

/// Sets the rigid body inertia for this body
void RigidBody::set_inertia(const SpatialRBInertiad& inertia)
{
  // set the inertia
  _Jm = Pose3d::transform(_jF, inertia);

  // invalidate the remaining inertias
  _Ji_valid = _Jj_valid = _J0_valid = _Jcom_valid = false;

  // see whether we can re-validate an acceleration
  if (inertia.pose == _F)
  {
    _Ji_valid = true;
    _Ji = inertia;
  }
  else if (inertia.pose == _F2)
  {
    _Jcom_valid = true;
    _Jcom = inertia;
  }
  else if (!is_base() && inertia.pose == get_inner_joint_explicit()->get_pose())
  {
    _Jj_valid = true;
    _Jj = inertia;
  }
  else if (inertia.pose == GLOBAL)
  {
    _J0_valid = true;
    _J0 = inertia;
  }
}

/// Sets the inertial pose for this rigid body
/**
 * Inertial pose should be defined relative to the rigid body pose
 */
void RigidBody::set_inertial_pose(const Pose3d& P)
{
  // verify that pose is set relative to body pose
  if (P.rpose != _F)
    throw std::runtime_error("RigidBody::set_inertial_pose() - inertial pose not defined relative to body pose");

  // set the inertial pose
  *_jF = P;

  // invalidate vectors using inertial frame
  _xdm_valid = _xddm_valid = _forcem_valid = false;
}

/// Gets the current sum of forces on this body
const SForced& RigidBody::sum_forces()
{
  switch (_rftype)
  {
    case eGlobal:
      if (!_force0_valid)
        _force0 = Pose3d::transform(GLOBAL, _forcecom);
      _force0_valid = true;
      return _force0;

    case eLink:
      if (!_forcei_valid)
        _forcei = Pose3d::transform(_F, _forcecom);
      _forcei_valid = true;
      return _forcei;

    case eLinkInertia:
      if (!_forcem_valid)
        _forcem = Pose3d::transform(_jF, _forcecom);
      _forcem_valid = true;
      return _forcem;

    case eLinkCOM:
      return _forcecom;

    case eJoint:
      if (!_forcej_valid)
        _forcej = Pose3d::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _forcecom);
      _forcej_valid = true;
      return _forcej;

    default:
      assert(false);
  }
}

/// Gets the current body velocity
const SVelocityd& RigidBody::get_velocity()
{
  switch (_rftype)
  {
    case eLinkCOM:
      return _xdcom;

    case eLink:
      if (!_xdi_valid)
        _xdi = Pose3d::transform(_F, _xdcom);
      _xdi_valid = true;
      return _xdi;

    case eLinkInertia:
      if (!_xdm_valid)
        _xdm = Pose3d::transform(_jF, _xdcom);
      _xdm_valid = true;
      return _xdm;

    case eJoint:
      if (!_xdj_valid)
        _xdj = Pose3d::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _xdcom);
      _xdj_valid = true;
      return _xdj;

    case eGlobal:
      if (!_xd0_valid)
        _xd0 = Pose3d::transform(GLOBAL, _xdcom);
      _xd0_valid = true;
      return _xd0;

    default:
      assert(false);
  }
}

/// Gets the body inertia
const SpatialRBInertiad& RigidBody::get_inertia()
{
  switch (_rftype)
  {
    case eGlobal:
      if (!_J0_valid)
        _J0 = Pose3d::transform(GLOBAL, _Jm);
      _J0_valid = true;
      return _J0;

    case eLink:
      if (!_Ji_valid)
        _Ji = Pose3d::transform(_F, _Jm);
      _Ji_valid = true;
      return _Ji;

    case eLinkCOM:
      if (!_Jcom_valid)
        _Jcom = Pose3d::transform(_F2, _Jm);
      _Jcom_valid = true;
      return _Jcom;

    case eLinkInertia:
      return _Jm;

    case eJoint:
      if (!_Jj_valid)
        _Jj = Pose3d::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _Jm);
      _Jj_valid = true;
      return _Jj;

    default:
      assert(false);
  }
}

/// Gets the current body acceleration
const SAcceld& RigidBody::get_accel()
{
  // do simplified case where body is disabled
  if (!is_enabled())
  {
    _xddi_valid = _xddm_valid = _xddj_valid = _xdd0_valid = true;
  }

  switch (_rftype)
  {
    case eLinkCOM:
      return _xddcom;

    case eLink:
      if (!_xddi_valid)
        _xddi = Pose3d::transform(_F, _xddcom);
      _xddi_valid = true;
      return _xddi;

    case eLinkInertia:
      if (!_xddm_valid)
        _xddm = Pose3d::transform(_jF, _xddcom);
      _xddm_valid = true;
      return _xddm;

    case eJoint:
      if (!_xddj_valid)
        _xddj = Pose3d::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _xddcom);
      _xddj_valid = true;
      return _xddj;

    case eGlobal:
      if (!_xdd0_valid)
        _xdd0 = Pose3d::transform(GLOBAL, _xddcom);
      _xdd0_valid = true;
      return _xdd0;

    default:
      assert(false);
  }
}

/// Resets the force accumulators on this body
void RigidBody::reset_accumulators()
{
  // clear forces
  _force0.set_zero();
  _forcei.set_zero();
  _forcem.set_zero();
  _forcej.set_zero();
  _forcecom.set_zero();

  // validate all forces
  _forcei_valid = true;
  _forcem_valid = true;
  _forcej_valid = true;
  _force0_valid = true;
}

/// Computes the torques (w x Jw) that come from the Euler component of the Newton-Euler equations
SForced RigidBody::calc_euler_torques()
{
  FILE_LOG(LOG_DYNAMICS) << "Calculating Euler torques for " << id << std::endl;
  const SVelocityd& xd = get_velocity();
  return xd.cross(get_inertia() * xd);
}

/// Sets the current 3D pose for this rigid body
/**
 * Also updates the transforms for associated visualization and collision data.
 */
void RigidBody::set_pose(const Pose3d& p)
{
  // verify that the two poses are relative to the same pose
  if (p.rpose != _F->rpose)
    throw std::runtime_error("RigidBody::set_pose() - relative pose is not correct");

  // update the pose
  *_F = p;

  // update the mixed pose
  _F2->set_identity();
  _F2->rpose = _F;
  _F2->update_relative_pose(GLOBAL);
  _F2->q.set_identity();

  // invalidate pose vectors
  invalidate_pose_vectors();
}

/// Invalidates pose quantities
void RigidBody::invalidate_pose_vectors()
{
  // invalidate vector quantities
  _forcei_valid = _forcem_valid = _force0_valid = false;
  _xdi_valid = _xdm_valid = _xd0_valid = false;
  _xddi_valid = _xddm_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<RigidBodyPtr> outer;
  BOOST_FOREACH(JointPtr j, _outer_joints)
  {
    if (j->get_constraint_type() == Joint::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<RigidBodyPtr>::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<RigidBodyPtr>::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// Gets the desired child link
RigidBodyPtr RigidBody::get_child_link(JointPtr j) const
{
  assert(_outer_joints.find(j) != _outer_joints.end());
  return j->get_outboard_link();
}

/// Sets the force on the body (this function is private b/c I can't imagine where it should be called by the user)
void RigidBody::set_force(const SForced& w)
{
  // do not add forces to disabled bodies
  if (!_enabled)
    return;

  // update the force
  _forcecom = Pose3d::transform(_F2, w);

  // see whether we update a force
  if (w.pose == _F)
  {
    _forcei_valid = true;
    _forcei = w;

    // invalidate the remaining forces
    _forcej_valid = _forcem_valid = _force0_valid = false;
  }
  else if (w.pose == _jF)
  {
    _forcem_valid = true;
    _forcem = w;

    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _force0_valid = false;
  }
  else if (!is_base() && w.pose == get_inner_joint_explicit()->get_pose())
  {
    _forcej_valid = true;
    _forcej = w;

    // invalidate the remaining forces
    _forcei_valid = _forcem_valid = _force0_valid = false;
  }
  else if (w.pose == GLOBAL)
  {
    _force0_valid = true;
    _force0 = w;

    // invalidate the remaining forces
    _forcei_valid = _forcem_valid = _forcej_valid = false;
  }
  else
    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _forcem_valid = _force0_valid = false;
}

/// Adds a force to the body
void RigidBody::add_force(const SForced& w)
{
  // do not add forces to disabled bodies
  if (!_enabled)
    return;

  // update the force
  _forcecom += Pose3d::transform(_F2, w);

  // see whether we update a force
  if (w.pose == _F)
  {
    if (_forcei_valid)
      _forcei += w;

    // invalidate the remaining forces
    _forcej_valid = _forcem_valid = _force0_valid = false;
  }
  else if (w.pose == _jF)
  {
    if (_forcem_valid)
      _forcem += w;

    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _force0_valid = false;
  }
  else if (!is_base() && w.pose == get_inner_joint_explicit()->get_pose())
  {
    if (_forcej_valid)
      _forcej += w;

    // invalidate the remaining forces
    _forcei_valid = _forcem_valid = _force0_valid  = false;
  }
  else if (w.pose == GLOBAL)
  {
    if (_force0_valid)
      _force0 += w;

    // invalidate the remaining forces
    _forcei_valid = _forcem_valid = _forcej_valid  = false;
  }
  else
    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _forcem_valid = _force0_valid  = false;
}

/// Calculates the velocity of a point on this rigid body in the body frame
Vector3d RigidBody::calc_point_vel(const Point3d& point) const
{
  // if the body is disabled, point velocity is zero
  if (!_enabled)
    return Vector3d::zero(_F);

  // convert point to a vector in the body frame
  Vector3d r = Pose3d::transform_point(_F, point);

  // get the velocity in the body frame
  SVelocityd xd = Pose3d::transform(_F, _xd0);

  // compute the point velocity - in the body frame
  Vector3d pv = xd.get_linear() + Vector3d::cross(xd.get_angular(), r);
  pv.pose = _F;
  return pv;
}

/// Implements Base::load_from_xml()
void RigidBody::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  const unsigned X = 0, Y = 1, Z = 2;
  map<std::string, BasePtr>::const_iterator id_iter;

  // load parent data
  SingleBody::load_from_xml(node, id_map);

  // ***********************************************************************
  // don't verify that the node is correct, b/c RigidBody can be subclassed
  // ***********************************************************************

#ifdef USE_OSG
  /// Color to add to the rigid body when rendered
  Ravelin::VectorNd color_rgba;

  // read the viscous dampening coefficient, if provided
  XMLAttrib* color_attr = node->get_attrib("color");
  if (color_attr){
    color_attr->get_vector_value(color_rgba);

    osg::Group* this_group = _vizdata->get_group();

    for(int i=0;i<this_group->getNumChildren();i++){
      osg::Node* n = this_group->getChild(i);
      CcolorVisitor  newColor;
      newColor.setColor( color_rgba[0], color_rgba[1], color_rgba[2], color_rgba[3] );
      n->accept( newColor );
    }
  }
#endif
  // read the viscous dampening coefficient, if provided
  XMLAttrib* viscous_coeff_attr = node->get_attrib("viscous-dampening-coeff");
  if (viscous_coeff_attr)
    viscous_coeff_attr->get_vector_value(viscous_coeff);

  // read whether the body is enabled, if provided
  XMLAttrib* enabled_attr = node->get_attrib("enabled");
  if (enabled_attr)
    _enabled = enabled_attr->get_bool_value();

  // read whether the body is compliant, if provided
  XMLAttrib* compliant_attr = node->get_attrib("compliant");
  if (compliant_attr)
    compliance = (compliant_attr->get_bool_value()) ? eCompliant : eRigid;

  // read the mass, if provided
  XMLAttrib* mass_attr = node->get_attrib("mass");
  if (mass_attr)
  {
    SpatialRBInertiad J = Pose3d::transform(_jF, get_inertia());
    J.m = mass_attr->get_real_value();
    set_inertia(J);
  }

  // read the limit bound expansion, if provided
  XMLAttrib* lbe_attr = node->get_attrib("limit-bound-expansion");
  if (lbe_attr)
    limit_bound_expansion = lbe_attr->get_real_value();

  // read the inertia matrix, if provided
  XMLAttrib* inertia_attr = node->get_attrib("inertia");
  if (inertia_attr)
  {
    SpatialRBInertiad J = Pose3d::transform(_jF, get_inertia());
    inertia_attr->get_matrix_value(J.J);
    set_inertia(J);
  }

  // read the position and orientation, if provided
  XMLAttrib* position_attr = node->get_attrib("position");
  XMLAttrib* rpy_attr = node->get_attrib("rpy");
  XMLAttrib* quat_attr = node->get_attrib("quat");
  if (position_attr || rpy_attr || quat_attr)
  {
    Pose3d T;
    if (position_attr)
      T.x = position_attr->get_origin_value();
    if (quat_attr)
      T.q = quat_attr->get_quat_value();
    else if (rpy_attr)
      T.q = rpy_attr->get_rpy_value();
    set_pose(T);
  }

  // read the inertial frame here...
  XMLAttrib* com_attr = node->get_attrib("inertial-relative-com");
  XMLAttrib* J_rpy_attr = node->get_attrib("inertial-relative-rpy");
  XMLAttrib* J_quat_attr = node->get_attrib("inertial-relative-quat");
  if (com_attr || J_rpy_attr || J_quat_attr)
  {
    // reset the inertial frame
    _jF->set_identity();

    // read the com
    if (com_attr)
      _jF->x = com_attr->get_origin_value();
    if (J_quat_attr)
      _jF->q = J_quat_attr->get_quat_value();
    else if (J_rpy_attr)
      _jF->q = J_rpy_attr->get_rpy_value();
  }

  // set the collision geometries, if provided
  list<shared_ptr<const XMLTree> > cg_nodes = node->find_child_nodes("CollisionGeometry");
  if (!cg_nodes.empty())
  {
    // ok to clear the set of geometries
    geometries.clear();

    // read in the collision geometries
    for (list<shared_ptr<const XMLTree> >::const_iterator i = cg_nodes.begin(); i != cg_nodes.end(); i++)
    {
      // create a new CollisionGeometry object
      CollisionGeometryPtr cg(new CollisionGeometry());

      // set the single body for the geometry
      cg->set_single_body(get_this());

      // populate the CollisionGeometry object
      cg->load_from_xml(*i, id_map);

      // add the collision geometry
      geometries.push_back(cg);
    }
  }

  // look for a inertia from primitives nodes
  // NOTE: we must do this step *before* setting velocities b/c setting
  // velocities updates momenta!
  list<shared_ptr<const XMLTree> > ifp_nodes = node->find_child_nodes("InertiaFromPrimitive");
  if (!ifp_nodes.empty())
  {
    // set inertia to zero initially
    SpatialRBInertiad J(_jF);

    // loop over all InertiaFromPrimitive nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = ifp_nodes.begin(); i != ifp_nodes.end(); i++)
    {
      // make sure the child node has the ID
      XMLAttrib* pid_attr = (*i)->get_attrib("primitive-id");
      if (!pid_attr)
      {
        cerr << "RigidBody::load_from_xml() - InertiaFromPrimitive node ";
        cerr << "has no" << endl << "  primitive-id attribute!";
        cerr << endl << "  offending node: " << endl << *node;
        continue;
      }

      // get the ID
      const std::string& ID = pid_attr->get_string_value();

      // attempt to find the ID
      if ((id_iter = id_map.find(ID)) == id_map.end())
      {
        cerr << "RigidBody::load_from_xml() - Primitive id: ";
        cerr << ID << " not found!" << endl << "  offending node: ";
        cerr << endl << *node;
        continue;
      }

      // get the primitive
      PrimitivePtr primitive = dynamic_pointer_cast<Primitive>(id_iter->second);

      // get the inertia from the primitive
      SpatialRBInertiad Jx = primitive->get_inertia();

      // convert the primitive's inertial frame
      // we want to treat the primitive's inertial frame as relative to the
      // rigid body's inertial frame
      shared_ptr<const Pose3d> Fx = primitive->get_inertial_pose();
      shared_ptr<Pose3d> Fxx(new Pose3d(*Fx));
      Fxx->update_relative_pose(GLOBAL);  // account for relative pose chain

      // now make the relative pose for Fxx be the inertial frame for this
      Fxx->rpose = _jF;

      // set the relative pose initially to identity for this primitive
      shared_ptr<Pose3d> rTR(new Pose3d);

      // read the relative transformation, if specified
      XMLAttrib* rel_origin_attr = (*i)->get_attrib("relative-origin");
      XMLAttrib* rel_rpy_attr = (*i)->get_attrib("relative-rpy");
      if (rel_origin_attr)
        rTR->x = rel_origin_attr->get_origin_value();
      if (rel_rpy_attr)
        rTR->q = rel_rpy_attr->get_rpy_value();
      rTR->rpose = Fxx;
      Jx.pose = rTR;

      // transform the inertia and update the inertia for this
      J += Pose3d::transform(_jF, Jx);
    }

    // set the mass and inertia of the RigidBody additively
    set_inertia(J);
  }

  // read the linear and/or velocity of the body, if provided
  XMLAttrib* lvel_attr = node->get_attrib("linear-velocity");
  XMLAttrib* avel_attr = node->get_attrib("angular-velocity");
  if (lvel_attr || avel_attr)
  {
    Vector3d lv = Vector3d::zero(), av = Vector3d::zero();
    shared_ptr<Pose3d> TARGET(new Pose3d);
    TARGET->rpose = _F;
    TARGET->q = Quatd::invert(_F->q);
    SVelocityd v;
    v.pose = TARGET;
    if (lvel_attr) lvel_attr->get_vector_value(lv);
    if (avel_attr) avel_attr->get_vector_value(av);
    v.set_linear(lv);
    v.set_angular(av);
    set_velocity(v);
  }

/*
  // read in the vector from the inner joint to the com in link coordinates
  XMLAttrib* d_attr = node->get_attrib("inner-joint-to-com-vector-link");
  if (d_attr)
  {
    Vector3 d;
    d_attr->get_vector_value(d);
    set_inner_joint_to_com_vector_link(d);
  }
*/
  // read the articulated body, if given
  XMLAttrib* ab_attr = node->get_attrib("articulated-body-id");
  if (ab_attr)
  {
    // get the ID
    const std::string& ID = ab_attr->get_string_value();

    // look for the ID -- only warn if it is not found
    if ((id_iter = id_map.find(ID)) == id_map.end())
    {
      FILE_LOG(LOG_DYNAMICS) << "RigidBody::load_from_xml() warning - ";
      FILE_LOG(LOG_DYNAMICS) << "articulated body" << endl << "  '" << ID << "' not ";
      FILE_LOG(LOG_DYNAMICS) << "found" << endl << "  ** This warning could result ";
      FILE_LOG(LOG_DYNAMICS) << "from links being constructed before articulated bodies ";
      FILE_LOG(LOG_DYNAMICS) << endl << "    and may not be serious..." << endl;
      FILE_LOG(LOG_DYNAMICS) << "  offending node: " << endl << *node;
    }
    else
      set_articulated_body(dynamic_pointer_cast<ArticulatedBody>(id_iter->second));
  }
}

/// Implements Base::save_to_xml()
void RigidBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // save parent data
  SingleBody::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "RigidBody";

  // save whether the body is enabled
  node->attribs.insert(XMLAttrib("enabled", _enabled));

  // save whether the body is compliant
  node->attribs.insert(XMLAttrib("compliant", compliance == eCompliant));

  // save the mass
  node->attribs.insert(XMLAttrib("mass", _Jm.m));

  // write the limit bound expansion
  node->attribs.insert(XMLAttrib("limit-bound-expansion", limit_bound_expansion));

  // save the inertia
  node->attribs.insert(XMLAttrib("inertia", _Jm.J));

  // convert the current pose to be with respect to global coordinates
  Pose3d F0 = *_F;
  F0.update_relative_pose(GLOBAL);
  node->attribs.insert(XMLAttrib("position", F0.x));
  node->attribs.insert(XMLAttrib("quat", F0.q));

  // save the inertial frame
  node->attribs.insert(XMLAttrib("inertial-relative-com", _jF->x));
  node->attribs.insert(XMLAttrib("inertial-relative-quat", _jF->q));

  // save the linear and angular velocities
  shared_ptr<Pose3d> TARGET(new Pose3d);
  TARGET->rpose = _F;
  TARGET->q = Quatd::invert(_F->q);
  SVelocityd v = Pose3d::transform(TARGET, _xd0);
  node->attribs.insert(XMLAttrib("linear-velocity", v.get_linear()));
  node->attribs.insert(XMLAttrib("angular-velocity", v.get_angular()));

  // save the dampening coefficients
  node->attribs.insert(XMLAttrib("viscous-coeff", viscous_coeff));

  // save all collision geometries
  BOOST_FOREACH(CollisionGeometryPtr g, geometries)
  {
    XMLTreePtr geom_node(new XMLTree("CollisionGeometry"));
    node->add_child(geom_node);
    g->save_to_xml(geom_node, shared_objects);
  }

  // save the ID articulated body (if any)
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody(_abody);
    node->attribs.insert(XMLAttrib("articulated-body-id", abody->id));
  }
/*
  // save the IDs of all child links and the vectors from the com to the outer
  // joints
  for (map<RigidBodyPtr, Vector3>::const_iterator i = _child_links.begin(); i != _child_links.end(); i++)
  {
    XMLTreePtr child_node(new XMLTree("ChildLink"));
    node->add_child(child_node);
    child_node->attribs.insert(XMLAttrib("link-id", i->first->id));
    child_node->attribs.insert(XMLAttrib("com-to-outer-vec-link", i->second));
  }
*/
}

/// Adds an inner joint for this link
/**
 * \param parent the outer link of the parent
 * \param j the joint connecting parent and this
 */
void RigidBody::add_inner_joint(JointPtr j)
{
  _inner_joints.insert(j);

  // update the spatial axes
  j->update_spatial_axes();

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1 = j->get_articulated_body();
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Adds an outer joint for this link
/**
 * \param j the joint connecting this and child
 * \note replaces the outer joint if it is already attached to this link
 */
void RigidBody::add_outer_joint(JointPtr j)
{
  // add the outer joint
  _outer_joints.insert(j);

  // update the spatial axes
  j->update_spatial_axes();

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1 = j->get_articulated_body();
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Determines whether the given link is a child link of this
bool RigidBody::is_child_link(shared_ptr<const RigidBody> query) const
{
  BOOST_FOREACH(JointPtr j, _outer_joints)
    if (RigidBodyPtr(j->get_outboard_link()) == query)
      return true;

  return false;
}

/// Determines whether the given link is a descendant of this
/**
 * \note returns <b>true</b> if query == this
 */
bool RigidBody::is_descendant_link(shared_ptr<const RigidBody> query) const
{
  queue<shared_ptr<const RigidBody> > q;

  // check for query == this
  if (query == shared_from_this())
    return true;

  // add all children to the queue
  BOOST_FOREACH(JointPtr j, _outer_joints)
    q.push(shared_ptr<const RigidBody>(j->get_outboard_link()));

  // continue processing children until no more children are able to be processed
  while (!q.empty())
  {
    shared_ptr<const RigidBody> link = q.front();
    q.pop();
    if (link == query)
      return true;
    BOOST_FOREACH(JointPtr j, link->_outer_joints)
      q.push(shared_ptr<const RigidBody>(j->get_outboard_link()));
  }

  return false;
}

/// Removes the specified outer joint from this link
void RigidBody::remove_outer_joint(JointPtr joint)
{
  _outer_joints.erase(joint);
}

/// Removes the specified outer joint from this link
/**
 * Returns true if the link was found.
 */
void RigidBody::remove_inner_joint(JointPtr joint)
{
  _inner_joints.erase(joint);
}

/// Applies a impulse to this link
/**
 * \param w the impulse as a force
 */
void RigidBody::apply_impulse(const SMomentumd& w)
{
  // if this is not an articulated body, just update linear and angular
  // momenta and velocites
  if (_abody.expired())
  {
    if (!_enabled)
      return;

    // get velocity update
    SMomentumd wx = Pose3d::transform(get_computation_frame(), w);
    SVelocityd dxd = get_inertia().inverse_mult(wx);

    // update linear and angular velocities
    _xdcom += Pose3d::transform(_F2, dxd);

    // see whether we can update any velocities
    if (dxd.pose == _F)
    {
      if (_xdi_valid)
        _xdi += dxd;
      _xdj_valid = _xdm_valid = _xd0_valid = false;
    }
    else if (dxd.pose == _jF)
    {
      if (_xdm_valid)
        _xdm += dxd;
      _xdj_valid = _xdi_valid = _xd0_valid = false;
    }
    else if (!is_base() && dxd.pose == get_inner_joint_explicit()->get_pose())
    {
      if (_xdj_valid)
        _xdj += dxd;
      _xdm_valid = _xdi_valid = _xd0_valid = false;
    }
    else if (dxd.pose == GLOBAL)
    {
      if (_xd0_valid)
        _xd0 += dxd;
      _xdm_valid = _xdi_valid = _xdj_valid = false;
    }
    else
      _xdm_valid = _xdi_valid = _xdj_valid = _xd0_valid = false;


    // reset the force and torque accumulators for this body
    reset_accumulators();
  }
  else
  {
    // get the articulated body
    ArticulatedBodyPtr abody(_abody);

    // apply the impulse to the articulated body
    abody->apply_impulse(w, get_this());
  }
}

/// Gets the generalized inertia of this rigid body
unsigned RigidBody::num_generalized_coordinates(GeneralizedCoordinateType gctype) const
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->num_generalized_coordinates(gctype);
  }
  else
    return num_generalized_coordinates_single(gctype);
}

/// Sets the generalized forces on the rigid body
void RigidBody::set_generalized_forces(const Ravelin::SharedVectorNd& gf)
{
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->add_generalized_force(gf);
    return;
  }

  // if we're still here, this is only an individual body
  assert(gf.size() == num_generalized_coordinates(DynamicBody::eSpatial));
  SForced w;

  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  // set the pose for w
  w.pose = _F2;

  // get the force and torque
  w.set_force(Vector3d(gf[0], gf[1], gf[2]));
  w.set_torque(Vector3d(gf[3], gf[4], gf[5]));

  // add the force to the sum of forces
  set_force(w);
}

/// Adds a generalized force to this rigid body
void RigidBody::add_generalized_force(const SharedVectorNd& gf)
{
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->add_generalized_force(gf);
    return;
  }

  // if we're still here, this is only an individual body
  assert(gf.size() == num_generalized_coordinates(DynamicBody::eSpatial));
  SForced w;

  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  // set the pose for w
  w.pose = _F2;

  // get the force and torque
  w.set_force(Vector3d(gf[0], gf[1], gf[2]));
  w.set_torque(Vector3d(gf[3], gf[4], gf[5]));

  // add the force to the sum of forces
  add_force(w);
}

/// Applies a generalized impulse to this rigid body
void RigidBody::apply_generalized_impulse(const SharedVectorNd& gj)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->apply_generalized_impulse(gj);
    return;
  }
  else
    apply_generalized_impulse_single(gj);
}

/// Applies a generalized impulse to this rigid body
void RigidBody::apply_generalized_impulse_single(const SharedVectorNd& gj)
{
  SMomentumd w;

  // don't do anything if this body is disabled
  if (!_enabled)
    return;

  // simple error check...
  assert(gj.size() == num_generalized_coordinates(DynamicBody::eSpatial));

  // clear the force accumulators (and validate them all)
  reset_accumulators();

  // get the impulses
  w.set_linear(Vector3d(gj[0], gj[1], gj[2]));
  w.set_angular(Vector3d(gj[3], gj[4], gj[5]));
  w.pose = _F2;

  // get the inertia in the link COM frame
  if (!_Jcom_valid)
  {
    _Jcom = Pose3d::transform(_F2, _Jm);
    _Jcom_valid = true;
  }

  // get the current velocity in the inertial frame
  SVelocityd v = _xdcom;

  // update the velocity
  v += _Jcom.inverse_mult(w);

  set_velocity(v);
}

/// Solves using the generalized inertia matrix
SharedMatrixNd& RigidBody::transpose_solve_generalized_inertia(const SharedMatrixNd& B, SharedMatrixNd& X)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->transpose_solve_generalized_inertia(B, X);
  }
  else
  {
    // get proper generalized inertia matrix
    const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
    MatrixNd M;
    M.resize(NGC, NGC);
    SharedMatrixNd Mshared = M.block(0, NGC, 0, NGC);
    get_generalized_inertia_inverse(Mshared);
    M.mult_transpose(B, X);
    return X;
  }
}

/// Solves using the generalized inertia matrix
SharedMatrixNd& RigidBody::solve_generalized_inertia(const SharedMatrixNd& B, SharedMatrixNd& X)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->solve_generalized_inertia(B, X);
  }
  else
    return solve_generalized_inertia_single(B, X);
}

/// Solves using the generalized inertia matrix (does not call articulated body version)
SharedMatrixNd& RigidBody::transpose_solve_generalized_inertia_single(const SharedMatrixNd& B, SharedMatrixNd& X)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  MatrixNd M;
  M.resize(NGC, NGC);
  SharedMatrixNd Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult_transpose(B, X);

  return X;
}

/// Solves using the generalized inertia matrix (does not call articulated body version)
SharedMatrixNd& RigidBody::solve_generalized_inertia_single(const SharedMatrixNd& B, SharedMatrixNd& X)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  MatrixNd M;
  M.resize(NGC, NGC);
  SharedMatrixNd Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult(B, X);

  return X;
}

/// Solves using the generalized inertia matrix
SharedVectorNd& RigidBody::solve_generalized_inertia(const SharedVectorNd& b, SharedVectorNd& x)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->solve_generalized_inertia(b, x);
  }
  else
    return solve_generalized_inertia_single(b, x);
}

/// Solves using the generalized inertia matrix
SharedVectorNd& RigidBody::solve_generalized_inertia_single(const SharedVectorNd& b, SharedVectorNd& x)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  MatrixNd M;
  M.resize(NGC, NGC);
  SharedMatrixNd Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult(b, x);

  return x;
}

/// Gets the generalized position of this rigid body
SharedVectorNd& RigidBody::get_generalized_coordinates(GeneralizedCoordinateType gctype, SharedVectorNd& gc)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->get_generalized_coordinates(gctype, gc);
  }
  else
    get_generalized_coordinates_generic(gctype, gc);

  return gc;
}

/// Sets the generalized coordinates of this rigid body
void RigidBody::set_generalized_coordinates(GeneralizedCoordinateType gctype, const SharedVectorNd& gc)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->set_generalized_coordinates(gctype, gc);
  }
  else
    set_generalized_coordinates_generic(gctype, gc);
}

/// Sets the generalized velocity of this rigid body
void RigidBody::set_generalized_velocity(GeneralizedCoordinateType gctype, const SharedVectorNd& gv)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->set_generalized_velocity(gctype, gv);
  }
  else
    set_generalized_velocity_generic(gctype, gv);
}

/// Gets the generalized velocity of this rigid body
SharedVectorNd& RigidBody::get_generalized_velocity(GeneralizedCoordinateType gctype, SharedVectorNd& gv)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->get_generalized_velocity(gctype, gv);
  }
  else
    get_generalized_velocity_generic(gctype, gv);

  return gv;
}

/// Gets the generalized acceleration of this body
SharedVectorNd& RigidBody::get_generalized_acceleration(SharedVectorNd& ga)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->get_generalized_acceleration(ga);
  }
  else
    get_generalized_acceleration_generic(ga);

  return ga;
}

/// Gets the generalized inertia of this rigid body
SharedMatrixNd& RigidBody::get_generalized_inertia(SharedMatrixNd& M)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->get_generalized_inertia(M);
  }
  else
    return get_generalized_inertia_single(M);
}

/// Gets the generalized inertia of this rigid body (does not call articulated body version)
SharedMatrixNd& RigidBody::get_generalized_inertia_single(SharedMatrixNd& M)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;

  // special case: disabled body
  if (!_enabled)
    return M.resize(0,0);

  // get the inertia
  SpatialRBInertiad J = Pose3d::transform(_F2, _Jm);

  // precompute some things
  Matrix3d hxm = Matrix3d::skew_symmetric(J.h * J.m);
  Matrix3d hxhxm = Matrix3d::skew_symmetric(J.h) * hxm;

  // arrange the matrix the way we want it: mass upper left, inertia lower right
  M.resize(SPATIAL_DIM, SPATIAL_DIM);
  M.set_sub_mat(0, 0, Matrix3d(J.m, 0, 0, 0, J.m, 0, 0, 0, J.m));
  M.set_sub_mat(3, 0, hxm);
  M.set_sub_mat(0, 3, hxm, Ravelin::eTranspose);
  M.set_sub_mat(3, 3, J.J - hxhxm);

  return M;
}

/// Gets the generalized inertia of this rigid body (does not call articulated body version)
SharedMatrixNd& RigidBody::get_generalized_inertia_inverse(SharedMatrixNd& M) const
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  static LinAlgd _LA;

  // don't invert inertia for disabled bodies
  if (!_enabled)
  {
    M.resize(0,0);
    return M;
  }

  // get the inertia
  SpatialRBInertiad J = Pose3d::transform(_F2, _Jm);

  // setup the matrix
  Matrix3d hx = Matrix3d::skew_symmetric(J.h);
  Matrix3d hxm = Matrix3d::skew_symmetric(J.h*J.m);
  M.resize(6,6);
  M.set_sub_mat(0,3, hxm, eTranspose);
  M.set_sub_mat(3,3, J.J - hx*hxm);
  M.set_sub_mat(0,0, Matrix3d(J.m, 0, 0, 0, J.m, 0, 0, 0, J.m));
  M.set_sub_mat(3,0, hxm);

  // invert the matrix
  _LA.invert(M);

  return M;
}

/// Gets the generalized inertia of this rigid body
SharedVectorNd& RigidBody::get_generalized_forces(SharedVectorNd& gf)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->get_generalized_forces(gf);
  }
  else
    return get_generalized_forces_single(gf);
}

/// Gets the generalized external forces (does not call articulated body version)
SharedVectorNd& RigidBody::get_generalized_forces_single(SharedVectorNd& gf)
{
  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // resize the generalized forces vector
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  gf.resize(NGC);

  // get force and torque
  Vector3d f = _forcecom.get_force();
  Vector3d t = _forcecom.get_torque();

  // setup the linear components of f
  gf[0] = f[0];
  gf[1] = f[1];
  gf[2] = f[2];
  gf[3] = t[0];
  gf[4] = t[1];
  gf[5] = t[2];

  return gf;
}

/// Converts a force to a generalized force
SharedVectorNd& RigidBody::convert_to_generalized_force(SingleBodyPtr body, const SForced& w, SharedVectorNd& gf)
{
  // if this belongs to an articulated body, call the articulated body method
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->convert_to_generalized_force(body, w, gf);
  }
  else
    return convert_to_generalized_force_single(body, w, gf);
}

/// Converts a force to a generalized force (does not call articulated body version)
SharedVectorNd& RigidBody::convert_to_generalized_force_single(SingleBodyPtr body, const SForced& w, SharedVectorNd& gf)
{
  // verify that body == this
  assert(body.get() == this);

  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // transform w to computation frame
  SForced wt = Pose3d::transform(_F2, w);

  // get linear and angular components of wt
  Vector3d f = wt.get_force();
  Vector3d t = wt.get_torque();

  // resize gf
  gf.resize(num_generalized_coordinates(DynamicBody::eSpatial));

  // setup the linear components
  gf[0] = f[0];
  gf[1] = f[1];
  gf[2] = f[2];
  gf[3] = t[0];
  gf[4] = t[1];
  gf[5] = t[2];

  return gf;
}

/// Calculates the kinetic energy of the body
double RigidBody::calc_kinetic_energy()
{
  if (!_enabled)
    return (double) 0.0;

/*
  const SVelocityd& xd = get_velocity();
  const SpatialRBInertiad& J = get_inertia();
*/
  SVelocityd xd = Pose3d::transform(get_gc_pose(), get_velocity());
  SpatialRBInertiad J = Pose3d::transform(get_gc_pose(), get_inertia());

  Vector3d v = xd.get_linear();
  Vector3d w = xd.get_angular();
  Vector3d wx = Vector3d(J.J*Origin3d(w), w.pose);
  return (v.norm_sq()*J.m + w.dot(wx))*0.5;

// return xd.dot(J * xd) * 0.5;
}

/// Gets the number of generalized coordinates
unsigned RigidBody::num_generalized_coordinates_single(DynamicBody::GeneralizedCoordinateType gctype) const
{
  const unsigned NGC_EULER = 7, NGC_SPATIAL = 6;

  // no generalized coordinates if this body is disabled
  if (!_enabled)
    return 0;

  // return the proper number of coordinates
  switch (gctype)
  {
    case DynamicBody::eEuler:
      return NGC_EULER;

    case DynamicBody::eSpatial:
      return NGC_SPATIAL;

    default:
      assert(false);
  }

  // make compiler happy
  assert(false);
  return 0;
}

/// Gets the first parent link of this link; returns NULL if there is no parent
RigidBodyPtr RigidBody::get_parent_link() const
{
  if (_inner_joints.size() > 1)
    throw std::runtime_error("Called RigidBody::get_parent_link() when multiple parent links present! It's not reasonable to call this method for links in maximal-coordinate articulated bodies.");

  // special case (no parent!)
  if (_inner_joints.empty())
    return RigidBodyPtr();

  JointPtr inner = *_inner_joints.begin();
  return RigidBodyPtr(inner->get_inboard_link());
}

/// Gets the explicit inner joint of this link; returns NULL if there is no explicit inner joint
/**
 * Throws an exception if this link has multiple explicit inner joints
 */
JointPtr RigidBody::get_inner_joint_explicit() const
{
  JointPtr ij;
  BOOST_FOREACH(JointPtr j, _inner_joints)
  {
    if (j->get_constraint_type() == Joint::eExplicit)
    {
      if (ij)
        throw std::runtime_error("Multiple explicit joints detected for a single link!");
      else
        ij = j;
    }
  }

  return ij;
}

/// Determines whether this link is a "ground" (fixed link)
bool RigidBody::is_ground() const
{
  // clear easy cases
  if (!_enabled)
    return true;

  // can't be a ground if not disabled and not part of an articulated body
  if (_abody.expired())
    return false;

  // now, case will differ depending on what type of articulated body this is
  ArticulatedBodyPtr ab(_abody);
  RCArticulatedBodyPtr rcab = dynamic_pointer_cast<RCArticulatedBody>(ab);
  if (rcab)
  {
    // check whether inner explicit joints are present (if none are present,
    // this is a base link)
    bool is_base = true;
    BOOST_FOREACH(JointPtr j, _inner_joints)
    {
      if (j->get_constraint_type() == Joint::eExplicit)
      {
        is_base = false;
        break;
      }
    }

    // if this link is a base and the base is fixed, it is a ground
    if (is_base && !rcab->is_floating_base())
      return true;
  }

  // still here? can't be a ground link
  return false;
}

/// Determines whether this link is the base
bool RigidBody::is_base() const
{
  // clear easy cases
  if (_abody.expired())
    return true;

  // check whether no explicit joints are present
  BOOST_FOREACH(JointPtr j, _inner_joints)
  {
    if (j->get_constraint_type() == Joint::eExplicit)
      return false;
  }

  // no explicit joints... it's the base
  return true;
}

/// Returns the ODE's for position and velocity (concatenated into x)
void RigidBody::ode_noexcept(SharedConstVectorNd& x, double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the number of generalized coordinates
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the shared pointer to this
  RigidBodyPtr shared_this = dynamic_pointer_cast<RigidBody>(shared_from_this());

  // get the generalized coordinates and velocity
  const SharedVectorNd gc = x.segment(0, NGC_EUL).get();
  const SharedVectorNd gv = x.segment(NGC_EUL, x.size()).get();

  // get the derivative of generalized coordinates and velocity
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, x.size());

  // set the state and velocity
  set_generalized_coordinates(DynamicBody::eEuler, gc);
  set_generalized_velocity(DynamicBody::eSpatial, gv);

  // check whether velocity limits have been exceeded
  if (!_vel_limit_exceeded)
    check_vel_limit_exceeded_and_update();

  // we need the generalized velocity as Rodrigues coordinates
  get_generalized_velocity(DynamicBody::eEuler, dgc);

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
  {
    FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << id << std::endl;
    (*controller)(shared_this, t, controller_arg);
  }

  // calculate forward dynamics at state x
  calc_fwd_dyn();
  get_generalized_acceleration(dgv);
}

/// Prepares to compute the ODE
void RigidBody::prepare_to_calc_ode(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the number of generalized coordinates
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the shared pointer to this
  RigidBodyPtr shared_this = dynamic_pointer_cast<RigidBody>(shared_from_this());

  // get the generalized coordinates and velocity
  const SharedVectorNd gc = x.segment(0, NGC_EUL).get();
  const SharedVectorNd gv = x.segment(NGC_EUL, x.size()).get();

  // set the state and velocity
  set_generalized_coordinates(DynamicBody::eEuler, gc);
  set_generalized_velocity(DynamicBody::eSpatial, gv);

  // check whether velocity limits have been exceeded
  check_vel_limit_exceeded_and_update();

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
  {
    FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << id << std::endl;
    (*controller)(shared_this, t, controller_arg);
  }
}

/// Computes the ODE
void RigidBody::ode(double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the number of generalized coordinates
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the shared pointer to this
  RigidBodyPtr shared_this = dynamic_pointer_cast<RigidBody>(shared_from_this());

  // get the derivative of generalized coordinates and velocity
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, dx.size());

  // we need the generalized velocity as Rodrigues coordinates
  get_generalized_velocity(DynamicBody::eEuler, dgc);

  // get the generalized acceleration
  get_generalized_acceleration(dgv);
}

/// Outputs the object state to the specified stream
/**
 * This method outputs all of the low-level details to the stream
 */
std::ostream& Moby::operator<<(std::ostream& out, Moby::RigidBody& rb)
{
  // write the ID of the object
  out << "rigid body ID: " << rb.id << endl;

  // indicate whether the body is enabled
  out << "  enabled? " << rb.is_enabled() << endl;

  // write the computation frame
  out << "  computation frame: ";
  switch (rb.get_computation_frame_type())
  {
    case eGlobal:        out << "global" << endl; break;
    case eLink:          out << "link inertia" << endl; break;
    case eLinkInertia:   out << "link" << endl; break;
    case eLinkCOM:       out << "link c.o.m." << endl; break;
    case eJoint:         out << "joint" << endl; break;
    default:
      assert(false);
  }

  out << "  Compliance type: ";
  switch (rb.compliance)
  {
    case RigidBody::eRigid:        out << "rigid" << endl; break;
    case RigidBody::eCompliant:    out << "compliant" << endl; break;
    default:
      assert(false);
  }

  // write inertial info
  shared_ptr<const Pose3d> jF = rb.get_inertial_pose();
  out << "  relative c.o.m.: " << jF->x << endl;
  out << "  relative inertial frame: " << AAngled(jF->q) << endl;
  out << "  mass: " << rb.get_inertia().m << endl;
  out << "  inertia: " << endl << rb.get_inertia().J;

  // write positions, velocities, and accelerations
  Pose3d F0 = *rb.get_pose();
  F0.update_relative_pose(GLOBAL);
  out << "  position: " << F0.x << endl;
  out << "  orientation: " << AAngled(F0.q) << endl;
  out << "  velocity (twist): " << rb.get_velocity() << endl;

  // write sum of forces
  out << "  forces: " << rb.sum_forces() << endl;

  // write the articulated body
  ArticulatedBodyPtr ab = rb.get_articulated_body();
  out << "  articulated body: " << ab << endl;

  // write all collision geometries
  out << "  collision geometries: ";
  BOOST_FOREACH(CollisionGeometryPtr g, rb.geometries)
    out << "    " << g << endl;

  return out;
}

