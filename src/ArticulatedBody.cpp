/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iostream>
#include <boost/foreach.hpp>
#include <queue>
#include <Moby/XMLTree.h>
#include <Moby/Joint.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/NumericalException.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/InvalidStateException.h>
#include <Moby/URDFReader.h>

using namespace Moby;
using namespace Ravelin;
using std::set;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::list;
using std::vector;
using std::map;
using std::string;
using std::queue;

ArticulatedBody::ArticulatedBody()
{
  // setup the default limit bound expansion
  limit_bound_expansion = 0.15;
}

/// Integrates a dynamic body
/*
void ArticulatedBody::integrate(double t, double h, shared_ptr<Integrator> integrator)
{
  FILE_LOG(LOG_DYNAMICS) << "ArticulatedBody::integrate() - integrating from " << t << " by " << h << std::endl;

  // reset the acceleration events exceeded
  _vel_limit_exceeded = false;

  if (_kinematic_update)
  {
    FILE_LOG(LOG_DYNAMICS) << " -- body set to kinematic update --" << std::endl;
    if (controller)
      (*controller)(dynamic_pointer_cast<ArticulatedBody>(shared_from_this()), t, controller_arg);
    return;
  }

  shared_ptr<ArticulatedBody> shared_this = dynamic_pointer_cast<ArticulatedBody>(shared_from_this());

  get_generalized_coordinates(eEuler, gc);
  get_generalized_velocity(eSpatial, gv); // gv depends on gc
  gcgv.resize(gc.size()+gv.size());
  gcgv.set_sub_vec(0, gc);
  gcgv.set_sub_vec(gc.size(), gv);
  integrator->integrate(gcgv, &ArticulatedBody::ode_both, t, h, (void*) &shared_this);
  gcgv.get_sub_vec(0, gc.size(), gc);
  gcgv.get_sub_vec(gc.size(), gcgv.size(), gv);
  // NOTE: velocity must be set first (it's computed w.r.t. old frame)
  set_generalized_coordinates(eEuler, gc);
  set_generalized_velocity(eSpatial, gv);
}
*/

void ArticulatedBody::reset_limit_estimates()
{
  // reset the acceleration events exceeded
  _vel_limit_exceeded = false;

  // clear the force limits on the individual rigid bodies
  BOOST_FOREACH(RigidBodyPtr link, _links)
    link->reset_limit_estimates();

  // clear the velocity limits
  const unsigned N = num_joint_dof();
  _vel_limits_lo.resize(N);
  _vel_limits_hi.resize(N);

  // update estimates
  for (unsigned i=0, j=0; i< _joints.size(); i++)
  {
    for (unsigned k=0; k< _joints[i]->num_dof(); k++, j++)
    {
      _vel_limits_lo[j] = _joints[i]->qd[k];
      if (_vel_limits_lo[k] < 0.0)
        _vel_limits_lo[k] *= (1.0 + limit_bound_expansion);
      else
        _vel_limits_lo[k] *= (1.0 - limit_bound_expansion);
      _vel_limits_hi[j] = _joints[i]->qd[k];
      if (_vel_limits_hi[k] < 0.0)
        _vel_limits_hi[k] *= (1.0 - limit_bound_expansion);
      else
        _vel_limits_hi[k] *= (1.0 + limit_bound_expansion);
    }
  }
}

/// Returns the ODE's for position and velocity (concatenated into x) without throwing an exception
void ArticulatedBody::ode_noexcept(SharedConstVectorNd& x, double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the shared pointer to this
  ArticulatedBodyPtr shared_this = dynamic_pointer_cast<ArticulatedBody>(shared_from_this());

  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the coordinates and velocity from x
  SharedConstVectorNd gc = x.segment(0, NGC_EUL);
  SharedConstVectorNd gv = x.segment(NGC_EUL, x.size());

  // set the state
  set_generalized_coordinates(DynamicBody::eEuler, gc);

  // set the velocity 
  set_generalized_velocity(DynamicBody::eSpatial, gv);

  // get the derivatives of coordinates and velocity from dx
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, x.size());

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

  // check whether velocity limits on the individual rigid bodies have been
  // violated
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->check_vel_limit_exceeded_and_update();

  // check whether joint velocities have been violated 
  check_joint_vel_limit_exceeded_and_update();
}

/// Prepares to compute the ODE  
void ArticulatedBody::prepare_to_calc_ode_sustained_constraints(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the shared pointer to this
  ArticulatedBodyPtr shared_this = dynamic_pointer_cast<ArticulatedBody>(shared_from_this());

  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the coordinates and velocity from x
  SharedConstVectorNd gc = x.segment(0, NGC_EUL);
  SharedConstVectorNd gv = x.segment(NGC_EUL, x.size());

  // set the state
  set_generalized_coordinates(DynamicBody::eEuler, gc);
  if (is_joint_constraint_violated())
    throw InvalidStateException();

  // set the velocity 
  set_generalized_velocity(DynamicBody::eSpatial, gv);

  // check whether velocity limits on the individual rigid bodies have been
  // violated
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->check_vel_limit_exceeded_and_update();

  // check whether joint velocities have been violated 
  check_joint_vel_limit_exceeded_and_update();

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

/// Prepares to compute the ODE  
void ArticulatedBody::prepare_to_calc_ode(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the shared pointer to this
  ArticulatedBodyPtr shared_this = dynamic_pointer_cast<ArticulatedBody>(shared_from_this());

  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the coordinates and velocity from x
  SharedConstVectorNd gc = x.segment(0, NGC_EUL);
  SharedConstVectorNd gv = x.segment(NGC_EUL, x.size());

  // set the state
  set_generalized_coordinates(DynamicBody::eEuler, gc);
  if (is_joint_constraint_violated())
    throw InvalidStateException();

  // set the velocity 
  set_generalized_velocity(DynamicBody::eSpatial, gv);

  // check whether velocity limits on the individual rigid bodies have been
  // violated
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    if (!rb->_vel_limit_exceeded)
      rb->check_vel_limit_exceeded_and_update();

  // check whether joint velocities have been violated 
  if (!_vel_limit_exceeded)
    check_joint_vel_limit_exceeded_and_update();

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
    (*controller)(shared_this, t, controller_arg);
}

/// Returns the ODE's for position and velocity (concatenated into x)
void ArticulatedBody::ode(double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the derivatives of coordinates and velocity from dx
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, dx.size());

  // we need the generalized velocity as Rodrigues coordinates
  get_generalized_velocity(DynamicBody::eEuler, dgc);

  // calculate forward dynamics 
  get_generalized_acceleration(dgv);
}

// TODO: remove this
/// Updates joint constraint violation (after integration)
void ArticulatedBody::update_joint_constraint_violations()
{
  // update the size of the constraint violation vector
  _cvio.resize(num_joint_dof());

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      _cvio[k] = 0.0;
      if (_joints[i]->q[j] < _joints[i]->lolimit[j])
        _cvio[k] = _joints[i]->lolimit[j] - _joints[i]->q[j];
      else if (_joints[i]->q[j] > _joints[i]->hilimit[j])
        _cvio[k] = _joints[i]->q[j] - _joints[i]->hilimit[j]; 
    }
  }
}

/// Checks for a joint constraint violation
bool ArticulatedBody::is_joint_constraint_violated() const
{
  // obvious check
  if (_cvio.size() != num_joint_dof())
    return false;

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      if (_joints[i]->q[j] + _cvio[k] < _joints[i]->lolimit[j] ||
          _joints[i]->q[j] - _cvio[k] > _joints[i]->hilimit[j])
        return true; 
    }
  }

  // no violation if here
  return false;
}

/// "Compiles" the articulated body
void ArticulatedBody::compile()
{
  _vel_limits_lo.resize(num_joint_dof(), 0.0);
  _vel_limits_hi.resize(num_joint_dof(), 0.0);
}

/// Finds the next event time on this articulated body for joint constraints
/**
 * \note skips current events
 */
double ArticulatedBody::find_next_joint_limit_time() const
{
  const double INF = std::numeric_limits<double>::max();
  const double VEL_TOL = NEAR_ZERO;

  // setup the maximum integration time
  double dt = std::numeric_limits<double>::max();

  // loop over all joints
  const vector<JointPtr>& joints = get_joints();
  for (unsigned i=0; i< joints.size(); i++)
  {
    // get the coordinate index for the joint
    unsigned k = joints[i]->get_coord_index();

    for (unsigned j=0; j< joints[i]->num_dof(); j++, k++)
    {
      // get the joint data
      const double q = joints[i]->q[j];
      const double qd = joints[i]->qd[j];
      const double l = joints[i]->lolimit[j];
      const double u = joints[i]->hilimit[j];

      // skip lower limit of DOF j of joint i if lower limit = -INF
      if (l > -INF)
      {
        // first check whether the limit is already exceeded
        if (q <= l)
          continue;

        // otherwise, determine when the joint limit will be met
        if (qd < -VEL_TOL)
          dt = std::min((l-q)/qd, dt);
      }

      // skip upper limit of DOF j of joint i if upper limit = INF
      if (u < INF)
      {
        // first check whether the limit is already exceeded
        if (q >= u)
          continue;

        // otherwise, determine when the joint limit will be met
        if (qd > VEL_TOL)
          dt = std::min((u-q)/qd, dt); 
      }
    }
  }

  return dt;
}

/// Computes the conservative advancement time on this articulated body for joint constraints
double ArticulatedBody::calc_CA_time_for_joints() const
{
  const double INF = std::numeric_limits<double>::max();

  // setup the maximum integration time
  double dt = std::numeric_limits<double>::max();

  // get the lower and upper velocity limits
  const vector<double>& vel_lo = _vel_limits_lo;
  const vector<double>& vel_hi = _vel_limits_hi; 
  assert(vel_lo.size() == num_joint_dof());
  assert(vel_hi.size() == num_joint_dof());

  // loop over all joints
  const vector<JointPtr>& joints = get_joints();
  for (unsigned i=0; i< joints.size(); i++)
  {
    // get the coordinate index for the joint
    unsigned k = joints[i]->get_coord_index();

    for (unsigned j=0; j< joints[i]->num_dof(); j++, k++)
    {
      // get the joint data
      const double q = joints[i]->q[j];
      const double qd_lo = vel_lo[k];
      const double qd_hi = vel_hi[k];
      const double l = joints[i]->lolimit[j];
      const double u = joints[i]->hilimit[j];

      // skip lower limit of DOF j of joint i if lower limit = -INF
      if (l > -INF)
      {
        // first check whether the limit is already exceeded
        if (q <= l)
          return 0.0;

        // find when lower limit would be exceeded: q + qd*t = l
        double t = (l-q)/qd_lo;
        if (t > 0.0)
          dt = std::min(t, dt);
      }

      // skip upper limit of DOF j of joint i if upper limit = INF
      if (u < INF)
      {
        // first check whether the limit is already exceeded
        if (q >= u)
          return 0.0;

        // find when upper limit would be exceeded: q + qd*t = u
        double t = (u - q)/qd_hi;
        if (t > 0.0)
          dt = std::min(t, dt);
      }
    }
  }

  return dt;
}

/// Validates the limit estimates
void ArticulatedBody::validate_limit_estimates()
{
  _vel_limit_exceeded = false;
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->validate_limit_estimates();
}

/// Determines whether the joint acceleration or link force estimates have been exceeded
bool ArticulatedBody::limit_estimates_exceeded() const
{
  // first look for the joint acceleration being exceeded
  if (_vel_limit_exceeded)
    return true;

  // now look for the link force estimates being exceeded
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    if (rb->limit_estimates_exceeded())
      return true;

  return false;
}

/// Updates the joint velocity limits
void ArticulatedBody::update_joint_vel_limits()
{
  _vel_limits_lo.resize(num_joint_dof());
  _vel_limits_hi.resize(num_joint_dof());

  // reset the velocity limits exceeded
  _vel_limit_exceeded = false;

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      if (_joints[i]->qd[j] < _vel_limits_lo[k])
      {
        _vel_limits_lo[k] = _joints[i]->qd[j];
        if (_vel_limits_lo[k] < 0.0)
          _vel_limits_lo[k] *= (1.0 + limit_bound_expansion);
        else
          _vel_limits_lo[k] *= (1.0 - limit_bound_expansion);
      }
      if (_joints[i]->qd[j] > _vel_limits_hi[k])
      {
        _vel_limits_hi[k] = _joints[i]->qd[j];
        if (_vel_limits_hi[k] < 0.0)
          _vel_limits_hi[k] *= (1.0 - limit_bound_expansion);
        else
          _vel_limits_hi[k] *= (1.0 + limit_bound_expansion);
      }
    }
  }
}

/// Checks whether a joint velocity exceeded the given limits, and updates the limits if necessary
void ArticulatedBody::check_joint_vel_limit_exceeded_and_update()
{
  // obvious check
  if (_vel_limits_lo.size() != num_joint_dof() ||
      _vel_limits_hi.size() != num_joint_dof())
    throw std::runtime_error("Joint velocities have not been setup!");

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      if (_joints[i]->qd[j] < _vel_limits_lo[k])
      {
        _vel_limits_lo[k] = _joints[i]->qd[j];
        if (_vel_limits_lo[k] < 0.0)
          _vel_limits_lo[k] *= (1.0 + limit_bound_expansion);
        else
          _vel_limits_lo[k] *= (1.0 - limit_bound_expansion);
        _vel_limit_exceeded = true;
      }
      if (_joints[i]->qd[j] > _vel_limits_hi[k])
      {
        _vel_limits_hi[k] = _joints[i]->qd[j];
        if (_vel_limits_hi[k] < 0.0)
          _vel_limits_hi[k] *= (1.0 - limit_bound_expansion);
        else
          _vel_limits_hi[k] *= (1.0 + limit_bound_expansion);
        _vel_limit_exceeded = true;
      } 
    }
  }
}

/// Gets the time-derivative of the Jacobian
/**
 * Columns correspond to joint coordinate indices.
 */
MatrixNd& ArticulatedBody::calc_jacobian_dot(boost::shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NEXP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NEXP_DOF + SPATIAL_DIM : NEXP_DOF;

  // setup the Jacobian
  J.set_zero(SPATIAL_DIM, NDOF); 

  // get the current link
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  
  // get the base link
  RigidBodyPtr base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    JointPtr joint = link->get_inner_joint_explicit();

    // get the parent link
    RigidBodyPtr parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVelocityd>& s = joint->get_spatial_axes_dot();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SharedVectorNd v = J.column(CIDX+i);
      Pose3d::transform(frame, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  // NOTE: we do not even check for a floating base, because the 
  // time-derivative of its Jacobian will always be zero

  return J;
}

/// Gets the Jacobian
/**
 * Columns correspond to joint coordinate indices.
 */
MatrixNd& ArticulatedBody::calc_jacobian(boost::shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NEXP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NEXP_DOF + SPATIAL_DIM : NEXP_DOF;

  // setup the Jacobian
  J.set_zero(SPATIAL_DIM, NDOF); 

  // get the current link
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  
  // get the base link
  RigidBodyPtr base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    JointPtr joint = link->get_inner_joint_explicit();

    // get the parent link
    RigidBodyPtr parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVelocityd>& s = joint->get_spatial_axes();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SharedVectorNd v = J.column(CIDX+i);
      Pose3d::transform(frame, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  // if base is floating, setup Jacobian columns at the end
  if (is_floating_base())
  {
    shared_ptr<const Pose3d> bpose = base->get_gc_pose();
    SharedMatrixNd Jbase = J.block(0, SPATIAL_DIM, NEXP_DOF, NEXP_DOF+SPATIAL_DIM);
    Pose3d::spatial_transform_to_matrix2(bpose, frame, Jbase);
  }

  return J;
}

/// Gets the maximum angular speed of the links of this articulated body
double ArticulatedBody::get_aspeed()
{
  double max_aspeed = (double) 0.0;
  for (unsigned i=0; i< _links.size(); i++)
  {
    double aspeed = _links[i]->get_aspeed();
    if (aspeed > max_aspeed)
      max_aspeed = aspeed;
  }

  return max_aspeed;
}

/// Determines the loop indices corresponding to each joint and the vector of links for each joint
void ArticulatedBody::find_loops(vector<unsigned>& loop_indices, vector<vector<unsigned> >& loop_links) const
{
  SAFESTATIC vector<JointPtr> loop_joints, implicit_joints;
  queue<RigidBodyPtr> q;

  // clear vectors
  loop_indices.resize(_joints.size());
  implicit_joints.clear();

  // get all implicit joints
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      implicit_joints.push_back(_joints[i]);

  // set all loop indices to INF (indicates no loop) initially
  for (unsigned i=0; i< _joints.size(); i++)
    loop_indices[i] = std::numeric_limits<unsigned>::max();

  // look for early exit
  if (implicit_joints.empty())
    return;

  // two cases: 1) body uses *only* implicit joints and 2) body uses 
  // explicit and implicit joints
  if (_joints.size() == implicit_joints.size())
  {
    // we're going to reset implicit_joints to hold only the joints that
    // complete loops
    implicit_joints.clear();
    for (unsigned i=0; i< _joints.size(); i++)
    {
      RigidBodyPtr inboard = _joints[i]->get_inboard_link();
      RigidBodyPtr outboard = _joints[i]->get_outboard_link();

      // check for obvious loop closure
      if (inboard->get_index() > outboard->get_index())
      {
        implicit_joints.push_back(_joints[i]);
        continue;
      }

      // check for not-so-obvious loop closure
      if (!outboard->is_enabled())
      {
        // outboard is fixed; look back to see whether one of the predecessor
        // links is fixed as well
        q.push(inboard);
        while (!q.empty())
        {
          RigidBodyPtr link = q.front();
          q.pop();
          if (!link->is_enabled())
          {
            implicit_joints.push_back(_joints[i]);
            break;
          }
          const set<JointPtr>& ij = link->get_inner_joints();
          BOOST_FOREACH(JointPtr j, ij)
            q.push(RigidBodyPtr(j->get_inboard_link()));
         }
       }
     }
  }

  // reset loop links
  loop_links.clear();
  loop_links.resize(implicit_joints.size());

  // for every kinematic loop
  for (unsigned k=0; k< implicit_joints.size(); k++)
  {
    // get the implicit joint
    JointPtr ejoint = implicit_joints[k];
    RigidBodyPtr outboard = ejoint->get_outboard_link();
    bool ground_outboard = outboard->is_ground();

    // determine all joints and links in the loop by iterating backward until
    // we get back to the first link in the loop
    loop_joints.clear();
    loop_joints.push_back(ejoint);
    loop_links[k].push_back(outboard->get_index());
    RigidBodyPtr inboard = ejoint->get_inboard_link();
    while (true)
    {
      JointPtr jx = inboard->get_inner_joint_explicit();
      loop_joints.push_back(jx);
      loop_links[k].push_back(inboard->get_index());
      inboard = jx->get_inboard_link();
 
      // check for loop termination
      if (inboard == outboard)
        break;
      if (ground_outboard && inboard->is_ground())
        break;
    }

    // reverse the vector of links so that it is (almost) sorted (all but
    // last link)
    std::reverse(loop_links[k].begin(), loop_links[k].end());
    #ifndef NDEBUG
    for (unsigned i=1; i< loop_links[k].size()-1; i++)
      assert(loop_links[k][i] > loop_links[k][i-1]);
    #endif

    // setup loop indices for each joint in the loop
    for (unsigned i=0; i< loop_joints.size(); i++)
      loop_indices[loop_joints[i]->get_index()] = k;
  }
}

/// Sets the vectors of links and joints
void ArticulatedBody::set_links_and_joints(const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  // copy the vector
  _links = links;

  // setup the link in the map 
  for (unsigned i=0; i< _links.size(); i++)
  {
    _links[i]->set_index(i);
    _links[i]->set_articulated_body(get_this());
  }

  // set vector of joints
  _joints = joints;

  // iterate over each joint
  for (unsigned i=0; i< _joints.size(); i++)
  {
    _joints[i]->set_index(i);
    _joints[i]->set_articulated_body(get_this());
  }

  compile();
}

/// Gets the number of explicit joint constraint equations
unsigned ArticulatedBody::num_constraint_eqns_explicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eExplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of implicit joint constraint equations
unsigned ArticulatedBody::num_constraint_eqns_implicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of joint degrees of freedom permitted by both implicit and explicit joint constraints
unsigned ArticulatedBody::num_joint_dof() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    ndof += _joints[i]->num_dof();
  return ndof;
}

/// Finds the joint with the given name
/**
 * \return NULL if the joint wasshared_ptr<void> not found
 */
JointPtr ArticulatedBody::find_joint(const string& jointname) const
{
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->id == jointname)
      return _joints[i];
      
  return JointPtr();
}

/// Gets the adjacent links
void ArticulatedBody::get_adjacent_links(list<sorted_pair<RigidBodyPtr> >& links) const
{
  for (unsigned i=0; i< _joints.size(); i++)
  {
    RigidBodyPtr ib = _joints[i]->get_inboard_link();
    RigidBodyPtr ob = _joints[i]->get_outboard_link();
    if (ib && ob)
      links.push_back(make_sorted_pair(ib, ob));
  }
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ArticulatedBody::translate(const Origin3d& x)
{
  // apply transform to all links
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->translate(x);
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ArticulatedBody::rotate(const Quatd& q)
{
  // apply transform to all links
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->rotate(q);
}

/// Calculates the combined kinetic energy of all links in this body
double ArticulatedBody::calc_kinetic_energy() 
{
  double KE = 0;
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    KE += rb->calc_kinetic_energy();

  return KE;
}

/// Resets force and torque accumulators on the body
/**
 * Force and torque accumulators on all links are reset.
 */
void ArticulatedBody::reset_accumulators()
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->reset_accumulators();
}

/// Finds the link with the given name
RigidBodyPtr ArticulatedBody::find_link(const string& linkid) const
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    if (rb->id == linkid)
      return rb;

  return RigidBodyPtr();
}

/// Updates visualization for the body
void ArticulatedBody::update_visualization()
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->update_visualization(); 
}

/// Loads a MCArticulatedBody object from an XML node
void ArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, std::map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // call parent method
  DynamicBody::load_from_xml(node, id_map);

  // don't verify the node name -- this class has derived classes
  // assert(strcasecmp(node->name().c_str(), "MCArticulatedBody") == 0);

  // see whether to load the model from a URDF file
  XMLAttrib* urdf_attr = node->get_attrib("urdf");
  if (urdf_attr)
  {
    // get the URDF filename
    std::string urdf_fname = urdf_attr->get_string_value();

    // load robots from the URDF
    std::string robot_name;
    std::vector<RigidBodyPtr> links;
    std::vector<JointPtr> joints; 
    if (URDFReader::read(urdf_fname, robot_name, links, joints))
      set_links_and_joints(links, joints);
    else
      std::cerr << "ArticulatedBody::load_from_xml()- unable to process URDF " << urdf_fname << std::endl;

    // do no more processing
    return;
  }

  // setup a list of joint nodes to find
  list<string> joint_node_names;
  joint_node_names.push_back("RevoluteJoint");
  joint_node_names.push_back("PrismaticJoint");
  joint_node_names.push_back("SphericalJoint");
  joint_node_names.push_back("UniversalJoint");
  joint_node_names.push_back("FixedJoint");
  joint_node_names.push_back("JointPlugin");

  // read the set of joint nodes and concatenate them into a single list
  list<shared_ptr<const XMLTree> > joint_nodes = node->find_child_nodes(joint_node_names);
  
  // read the set of link nodes
  list<shared_ptr<const XMLTree> > link_nodes = node->find_child_nodes("RigidBody");

  // if there were links read or joints read, add them 
  if (!joint_nodes.empty() || !link_nodes.empty())
  {
    // setup new lists for joints and links
    list<JointPtr> joints;
    list<RigidBodyPtr> links;

    // process all link nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = link_nodes.begin(); i != link_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body links are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the link was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the link
      links.push_back(dynamic_pointer_cast<RigidBody>(id_iter->second));
    }

    // process all joint nodes in the same manner
    for (list<shared_ptr<const XMLTree> >::const_iterator i = joint_nodes.begin(); i != joint_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body joints are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the joint was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the joints
      joints.push_back(dynamic_pointer_cast<Joint>(id_iter->second));
    }

    // set the joints and links
    set_links_and_joints(vector<RigidBodyPtr>(links.begin(), links.end()), vector<JointPtr>(joints.begin(), joints.end()));
  }

  // read the joint limit bound expansion, if provided
  XMLAttrib* jlbe_attr = node->get_attrib("joint-limit-bound-expansion");
  if (jlbe_attr)
    limit_bound_expansion = jlbe_attr->get_real_value();

  // read the link limit bound expansion, if provided
  XMLAttrib* llbe_attr = node->get_attrib("link-limit-bound-expansion");
  if (llbe_attr)
    BOOST_FOREACH(RigidBodyPtr rb, _links)
      rb->limit_bound_expansion = llbe_attr->get_real_value();
}

/// Saves this object to a XML tree
void ArticulatedBody::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call parent method
  DynamicBody::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "ArticulatedBody";

  // add all links
  for (unsigned i=0; i< _links.size(); i++)
  {
    // create the new node in the tree
    XMLTreePtr child_node(new XMLTree("RigidBody"));
    node->add_child(child_node);

    // write to this node
    _links[i]->save_to_xml(child_node, shared_objects);
  }

  // add all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // create the new node in the tree -- note that save_to_xml() should set
    // the node name correctly
    XMLTreePtr child_node(new XMLTree("Joint"));
    node->add_child(child_node);

    // write to this node
    _joints[i]->save_to_xml(child_node, shared_objects);
  }

  // write the limit bound expansion
  node->attribs.insert(XMLAttrib("limit-bound-expansion", limit_bound_expansion));
}

