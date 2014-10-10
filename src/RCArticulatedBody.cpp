/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <stack>
#include <queue>
#include <Moby/Log.h>
#include <Moby/Joint.h>
#include <Moby/RigidBody.h>
#include <Moby/CRBAlgorithm.h>
#include <Moby/XMLTree.h>
#include <Moby/FSABAlgorithm.h>
#include <Moby/Spatial.h>
#include <Moby/NumericalException.h>
#include <Moby/RCArticulatedBody.h>

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::queue;
using std::list;
using std::map;
using std::string;
using namespace Ravelin;
using namespace Moby;

/// Default constructor
/**
 * Constructs a reduced-coordinate articulated body with no joints and no links.
 */
RCArticulatedBody::RCArticulatedBody()
{
  _floating_base = false;
  _n_joint_DOF_explicit = 0;

  // create the linear algebra object
  _LA = shared_ptr<LinAlgd>(new LinAlgd);
  _fsab._LA = _LA;
  _crb._LA = _LA;

  // set default algorithm to FSAB and computation frame to link c.o.m. 
  algorithm_type = eFeatherstone;
  set_computation_frame_type(eLinkCOM);

  // setup baumgarte parameters
  b_alpha = (double) 0.0;
  b_beta = (double) 0.0;

  // invalidate position quanitites
  _position_invalidated = true;
}

/// Validates position variables
void RCArticulatedBody::validate_position_variables()
{
  _position_invalidated = false;
}

/// Gets the frame used for generalized coordinate calculations
shared_ptr<const Pose3d> RCArticulatedBody::get_gc_pose() const
{
  if (_links.empty())
    throw std::runtime_error("Cannot return gc pose when articulated body has no links");

  return _links.front()->get_gc_pose();
}

/// Sets the computation frame type
void RCArticulatedBody::set_computation_frame_type(ReferenceFrameType rftype)
{
  // set the reference frame
  _rftype = rftype;

  // invalidate
  _position_invalidated = true;

  // set the reference frame type for all links
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->set_computation_frame_type(rftype);
}

/// Determines whether all of the children of a link have been processed
bool RCArticulatedBody::all_children_processed(RigidBodyPtr link) const
{
  const std::set<JointPtr>& joints = link->get_outer_joints();
  BOOST_FOREACH(JointPtr j, joints)
  {
    RigidBodyPtr child = j->get_outboard_link();
    if (!_processed[child->get_index()])
      return false;
  }

  return true;
}

/// Gets the number of generalized coordinates for this body
unsigned RCArticulatedBody::num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const
{
  // look for trivial case
  if (_links.empty())
    return 0;

  if (!_floating_base)
    return _n_joint_DOF_explicit;
  else
    return _n_joint_DOF_explicit + _links.front()->num_generalized_coordinates_single(gctype);
}

/// Updates inverse generalized inertia matrix, as necessary
void RCArticulatedBody::update_factorized_generalized_inertia()
{
  // see whether we need to update
  if (!_position_invalidated)
    return;

  // get the body
  RCArticulatedBodyPtr body = dynamic_pointer_cast<RCArticulatedBody>(shared_from_this());

  // do precalculation on the body
  if (algorithm_type == eFeatherstone)
    _fsab.calc_spatial_inertias(body);
  else
    _crb.precalc(body);

  // indicate factorized inertia is valid
  validate_position_variables(); 
}

/// Solves using a generalized inertia matrix
SharedVectorNd& RCArticulatedBody::solve_generalized_inertia(const SharedVectorNd& v, SharedVectorNd& result)
{
  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // make x/b one vector
    result = v;

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // store the body's computation reference frame type
    ReferenceFrameType rftype = get_computation_frame_type();

    // set the reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(eLinkCOM);

    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // make x/b one vector
    result = v;

    // solve once
    _crb.M_solve_noprecalc(result);

    // revert the link reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(rftype); 
  }

  return result;
}

/// Solves the transpose using a generalized inertia matrix
SharedMatrixNd& RCArticulatedBody::transpose_solve_generalized_inertia(const SharedMatrixNd& m, SharedMatrixNd& result)
{
  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    MatrixNd::transpose(m, result);

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // store the body's computation reference frame type
    ReferenceFrameType rftype = get_computation_frame_type();

    // set the reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(eLinkCOM);

    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    MatrixNd::transpose(m, result);

    // solve
    _crb.M_solve_noprecalc(result);

    // revert the link reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(rftype); 
  }

  return result;
}

/// Solves using a generalized inertia matrix
SharedMatrixNd& RCArticulatedBody::solve_generalized_inertia(const SharedMatrixNd& m, SharedMatrixNd& result)
{
  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    result = m;

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // store the body's computation reference frame type
    ReferenceFrameType rftype = get_computation_frame_type();

    // set the reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(eLinkCOM);

    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    result = m;

    // solve
    _crb.M_solve_noprecalc(result);

    // revert the link reference frame type
    if (rftype != eLinkCOM)
      set_computation_frame_type(rftype); 
  }

  return result;
}

/// Applies a generalized impulse to the articulated body
void RCArticulatedBody::apply_generalized_impulse(const SharedVectorNd& gj)
{
  if (algorithm_type == eFeatherstone)
    _fsab.apply_generalized_impulse(gj);
  else
  {
    assert(algorithm_type == eCRB);

    // setup work variables
    SAFESTATIC VectorNd gv, gv_delta;

    // get the current generalized velocity
    get_generalized_velocity(DynamicBody::eSpatial, gv);

    // we'll solve for the change in generalized velocity
    DynamicBody::solve_generalized_inertia(gj, gv_delta);

    // apply the change in generalized velocity
    gv += gv_delta;
    set_generalized_velocity(DynamicBody::eSpatial, gv);
  }

  // reset the force and torque accumulators
  reset_accumulators();
}

/// Sets the generalized forces for the articulated body
void RCArticulatedBody::set_generalized_forces(const SharedVectorNd& gf)
{
  unsigned index = 0;
  SForced f0;

  if (_floating_base)
  {
    // get the base
    RigidBodyPtr base = _links.front();

    // first, get the force on the base link
    gf.get_sub_vec(num_joint_dof_explicit(), gf.size(), f0);

    // add the force to the base
    SForced fx = Pose3d::transform(base->get_gc_pose(), f0);
    base->set_force(fx);
  }

  // add to joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    SharedConstVectorNd f = gf.segment(idx, idx+_ejoints[i]->num_dof());
    _ejoints[i]->force = f;
  }
}

/// Adds a generalized force to the articulated body
void RCArticulatedBody::add_generalized_force(const SharedVectorNd& gf)
{
  unsigned index = 0;

  if (_floating_base)
  {
    // get the base
    RigidBodyPtr base = _links.front();

    // first, get the force on the base link
    SForced f0;
    gf.get_sub_vec(num_joint_dof_explicit(), gf.size(), f0);

    // add the force to the base
    f0.pose = base->get_gc_pose();
    base->add_force(f0);
  }

  // add to joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    SharedConstVectorNd f = gf.segment(idx, idx+_ejoints[i]->num_dof());
    _ejoints[i]->force += f;
  }
}

/// Determines whether the link is effectively a leaf link
bool RCArticulatedBody::treat_link_as_leaf(RigidBodyPtr link) const
{
  // if all children have lower link indices, link treatable as end-effector
  const std::set<JointPtr>& joints = link->get_outer_joints();
  BOOST_FOREACH(JointPtr j, joints)
  {
    RigidBodyPtr child = j->get_outboard_link();
    if (child->get_index() > link->get_index())
      return false;
  }

  return true;
}

/// Sets whether the base of this body is "floating" (or fixed)
void RCArticulatedBody::set_floating_base(bool flag)
{
  _floating_base = flag;
  compile();
}

/// Compiles this body (updates the link transforms and velocities)
void RCArticulatedBody::compile()
{
  // call parent method first
  ArticulatedBody::compile();

  // verify all links are enabled
  if (!is_floating_base())
  {
    for (unsigned i=1; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("Only first link can be disabled in a reduced coordinate body with fixed-base!");
  }
  else
  {
    for (unsigned i=0; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("No links can be disabled in a reduced coordinate body with floating-base!");
  }

  // update processed vector size
  _processed.resize(_links.size());

  // setup explicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx = 0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->set_coord_index(cidx);
    _ejoints[i]->set_constraint_index(ridx);
    cidx += _ejoints[i]->num_dof();
    ridx += _ejoints[i]->num_constraint_eqns();
  }

  // setup explicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->set_coord_index(cidx);
    _ijoints[i]->set_constraint_index(ridx);
    cidx += _ijoints[i]->num_dof();
    ridx += _ijoints[i]->num_constraint_eqns();
  }

  // point both algorithms to this body
  _crb.set_body(get_this());
  _fsab.set_body(get_this());

  // update link transforms and velocities
  update_link_poses();
  update_link_velocities();
}

/// Sets the vector of links and joints
void RCArticulatedBody::set_links_and_joints(const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  // setup the processed vector
  _processed.resize(links.size());

  // clear the vectors of joints
  _ejoints.clear();
  _ijoints.clear();

  // set the computation frame type on all links
  for (unsigned i=0; i< links.size(); i++)
    links[i]->set_computation_frame_type(_rftype);

  // find the base link and setup a processed map
  map<RigidBodyPtr, bool> processed;
  RigidBodyPtr base;
  for (unsigned i=0; i< joints.size(); i++)
  {
    RigidBodyPtr inboard = joints[i]->get_inboard_link();
    RigidBodyPtr outboard = joints[i]->get_outboard_link();
    RigidBodyPtr parent = inboard->get_parent_link();
    if (!parent)
    {
      if (base && inboard != base)
        throw std::runtime_error("Multiple base links detected!");
      base = inboard;
    }
    processed[inboard] = false;
    processed[outboard] = false;
  }

  // if there is no clearly defined base link, there are no links *and* there
  // are joints, we can't do anything
  if (joints.empty() && links.size() == 1)
    base = links.front();
  if (!base)
    throw std::runtime_error("Could not find base link!");

  // start processed at the base link
  queue<RigidBodyPtr> link_queue;
  link_queue.push(base);
  while (!link_queue.empty())
  {
    // get the link at the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // if the link has already been processed, no need to process it again
    if (processed[link])
      continue;

    // get all outer joints for this link
    BOOST_FOREACH(JointPtr joint, link->get_outer_joints())
    {
      // see whether the child has already been processed
      RigidBodyPtr child(joint->get_outboard_link());
      if (processed[child])
        _ijoints.push_back(joint);
      else
      {
        link_queue.push(child);
        _ejoints.push_back(joint);
      }
    }

    // indicate that the link has been processed
    processed[link] = true;
  }

  // recalculate the explicit joint degrees-of-freedom of this body
  _n_joint_DOF_explicit = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    _n_joint_DOF_explicit += _ejoints[i]->num_dof();

  // mark joints as the correct type
  for (unsigned i=0; i< _ejoints.size(); i++)
    _ejoints[i]->set_constraint_type(Joint::eExplicit);
  for (unsigned i=0; i< _ijoints.size(); i++)
    _ijoints[i]->set_constraint_type(Joint::eImplicit);

  // check to see whether user's numbering scheme is acceptable
  for (unsigned i=1; i< links.size(); i++)
  {
    // look for an unknown constraint
    BOOST_FOREACH(JointPtr joint, links[i]->get_inner_joints())
      if (joint->get_constraint_type() == Joint::eUnknown)
        throw std::runtime_error("Unknown constraint type found!");

    // no unknown constraint; look for an explicit constraint
    if (!links[i]->get_inner_joint_explicit())
      throw std::runtime_error("Nonzero link does not have an inner explicit joint!");
  }

  // look whether it's a floating base
  _floating_base = base->is_enabled();

  // call the parent method to update the link indices, etc.
  ArticulatedBody::set_links_and_joints(links, joints);
}

/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as
 * well as the base momentum; therefore, the derivative of the state vector is
 * composed of the joint-space accelerations and base forces (and torques).
 */
/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as
 * well as the base momentum; therefore, the derivative of the state vector is
 * composed of the joint-space accelerations and base forces (and torques).
 */
SharedVectorNd& RCArticulatedBody::get_generalized_acceleration(SharedVectorNd& ga)
{
  get_generalized_acceleration_generic(ga);
  return ga;
}

/// Updates the transforms of the links based on the current joint positions
/**
 * \note this doesn't actually calculate other than the joint positions; all
 *       links are defined with respect to the joints, which are defined
 *       with respect to their inner link
 */
void RCArticulatedBody::update_link_poses()
{
  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_poses() entered" << std::endl;

  // indicate factorized inertia matrix is no longer valid
  _position_invalidated = true;

  // update all joint poses
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->get_induced_pose();

  // print all link poses and joint poses
  if (LOGGING(LOG_DYNAMICS))
  {
    for (unsigned i=0; i< _links.size(); i++)
    {
      Transform3d Tx = Pose3d::calc_relative_pose(_links[i]->get_pose(), GLOBAL);
      FILE_LOG(LOG_DYNAMICS) << "  link " << _links[i]->id << " pose (relative to global frame): " << Tx.x << " " << AAngled(Tx.q) << std::endl;
    }
    for (unsigned i=0; i< _joints.size(); i++)
    {
      Transform3d Tx = Pose3d::calc_relative_pose(_joints[i]->get_pose(), GLOBAL);
      FILE_LOG(LOG_DYNAMICS) << "  joint " << _joints[i]->id << " pose (relative to global frame): " << Tx.x << " " << AAngled(Tx.q) << std::endl;
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_poses() exited" << std::endl;
}

/// Updates the link velocities
void RCArticulatedBody::update_link_velocities()
{
  queue<RigidBodyPtr> link_queue;
  vector<SVelocityd> sprime;

  // look for easy exit
  if (_links.empty() || _joints.empty())
    return;

  // get the base link
  RigidBodyPtr base = _links.front();

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_velocities() entered" << std::endl;

  // add all children of the base to the link queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links));
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);

  // reset processed vector
  for (unsigned i=0; i< _links.size(); i++)
    _processed[i] = false;
  _processed[base->get_index()] = true;

  // propagate link velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr outboard = link_queue.front();
    shared_ptr<const Pose3d> opose = outboard->get_computation_frame();
    link_queue.pop();
    unsigned i = outboard->get_index();

    // and push all children of the link onto the queue
    list<RigidBodyPtr> child_links;
    outboard->get_child_links(std::back_inserter(child_links));
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      if (!_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the inboard link
    JointPtr joint(outboard->get_inner_joint_explicit());
    RigidBodyPtr inboard(joint->get_inboard_link());
    shared_ptr<const Pose3d> ipose = inboard->get_computation_frame();
    unsigned h = inboard->get_index();

    // set this link's velocity to the parent's link velocity
    outboard->set_velocity(Pose3d::transform(opose, inboard->get_velocity()));

    // get the (transformed) link spatial axis
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    Pose3d::transform(opose, s, sprime);

    // determine the link velocity due to the parent velocity + joint velocity
    if (sprime.empty())
      outboard->set_velocity(outboard->get_velocity());
    else
      outboard->set_velocity(outboard->get_velocity() + mult(sprime, joint->qd));

    // indicate that the link has been processed
    _processed[i] = true;

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << outboard->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent velocity: " << inboard->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- link velocity : " << outboard->get_velocity()  << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_velocities() exited" << std::endl;
}

/// Calculates the column of a Jacobian matrix for a floating base with respect to a given point
/**
 * \param point the point in 3D space which the Jacobian is calculated against
 * \return a pointer to a 6x6 matrix; top three dimensions will be linear
 *         velocity and bottom three dimensions will be angular velocity
 */
/*
MatrixNd& RCArticulatedBody::calc_jacobian_floating_base(const Point3d& point, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC vector<SVelocityd> sbase, sbase_prime;
  SAFESTATIC shared_ptr<Pose3d> P(new Pose3d);

  // get the base link and the base pose
  RigidBodyPtr base = get_base_link();
  shared_ptr<const Pose3d> baseP = base->get_gc_pose();

  // setup the twists at the base link - first three vectors are linear motion
  // and second three are angular
  sbase.clear();
  sbase.resize(SPATIAL_DIM, SVelocityd::zero(baseP));
  sbase[3][0] = 1.0;  sbase[4][1] = 1.0;  sbase[5][2] = 1.0;
  sbase[0][3] = 1.0;  sbase[1][4] = 1.0;  sbase[2][5] = 1.0;

  // convert the poses to the point frame
  P->x = Origin3d(point);
  Pose3d::transform(P, sbase, sbase_prime);

  // init the base Jacobian
  J.resize(SPATIAL_DIM, SPATIAL_DIM);
  for (unsigned i=0; i< SPATIAL_DIM; i++)
  {
    SharedVectorNd Ji = J.column(i);
    sbase_prime[i].transpose_to_vector(Ji);
  }

  return J;
}

/// Calculates the Jacobian for the current robot configuration at a given point and with respect to a given link
MatrixNd& RCArticulatedBody::calc_jacobian(const Point3d& p, RigidBodyPtr link, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC MatrixNd Jsub;

  // resize the Jacobian
  J.set_zero(SPATIAL_DIM, num_generalized_coordinates(DynamicBody::eSpatial));

  // get the base link
  RigidBodyPtr base = get_base_link();

  if (is_floating_base())
  {
shared_ptr<Pose3d> frame(new Pose3d);
frame->rpose = p.pose;
frame->q.set_identity();
frame->x = Origin3d(p);

    // construct the spatial transform
    Pose3d::spatial_transform_to_matrix2(base->get_gc_pose(), frame, Jsub);

    // setup the floating base
    J.set_sub_mat(0, num_joint_dof_explicit(),Jsub);
  }

  // calculate all relevant columns
  while (link != base)
  {
    JointPtr joint = link->get_inner_joint_explicit();
    calc_jacobian_column(joint, p, Jsub);
    J.set_sub_mat(0, joint->get_coord_index(), Jsub);
    link = link->get_parent_link();
  }

  return J;
}
*/

/// Calculates the column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \param base_pose the transform to use for the base
 * \param q a mapping from joints to joint positions to use; joints without
 *        positions specified will be set to their current values
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 * \note this method works by changing the joint values, recomputing all link
 *       transforms, computing the Jacobian, restoring the old joint values and
 *       link transforms; thus, it will be slower than a special purpose method
 */
/*
MatrixNd& RCArticulatedBody::calc_jacobian(const Point3d& point, const Pose3d& base_pose, const map<JointPtr, VectorNd>& q, RigidBodyPtr link, MatrixNd& Jc)
{
  // store current joint values
  map<JointPtr, VectorNd> currentQ;
  for (unsigned i=0; i< _ejoints.size(); i++)
    currentQ[_ejoints[i]] = _ejoints[i]->q;

  // overwrite current joint values
  for (map<JointPtr, VectorNd>::const_iterator i = q.begin(); i != q.end(); i++)
    i->first->q = i->second;

  // save the current base pose and set it to that desired
  RigidBodyPtr base = get_base_link();
  Pose3d saved_base_pose = *base->get_pose();
  base->set_pose(base_pose);

  // update link transforms
  update_link_poses();

  // compute and store the Jacobian
  calc_jacobian(point, link, Jc);

  // restore joint values
  for (map<JointPtr, VectorNd>::const_iterator i = currentQ.begin(); i != currentQ.end(); i++)
    i->first->q = i->second;

  // restore the base pose and transforms
  base->set_pose(saved_base_pose);
  update_link_poses();

  return Jc;
}
*/

/// Calculates column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 */
MatrixNd& RCArticulatedBody::calc_jacobian_column(JointPtr joint, const Point3d& point, MatrixNd& Jc)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC shared_ptr<Pose3d> target(new Pose3d);
  SAFESTATIC vector<SVelocityd> sprime;

  // NOTE: spatial algebra provides us with a simple means to compute the
  // Jacobian of a joint with respect to a point.  The spatial axis of the
  // joint transforms joint velocity to spatial velocity. The spatial
  // transform is then used to transform the spatial velocity into a
  // convenient frame: the orientation of that frame is aligned with the
  // global frame, and the origin of the frame is aligned with the desired
  // point.  Applying the spatial transform from the spatial axis (either in
  // link or global frame) to the new frame will give us the desired vector.

  // get the spatial axis of the joint
  const vector<SVelocityd>& s = joint->get_spatial_axes();
  Jc.resize(SPATIAL_DIM, s.size());

  // compute a spatial transformation using the point as the target frame
  target->rpose = point.pose;
  target->x = Origin3d(point);
  target->q.set_identity();
  Pose3d::transform(target, s, sprime);

  // calculate the Jacobian column
  for (unsigned i=0; i< sprime.size(); i++)
  {
    Vector3d top = sprime[i].get_upper();
    Vector3d bot = sprime[i].get_lower();
    Jc(0,i) = bot[0];
    Jc(1,i) = bot[1];
    Jc(2,i) = bot[2];
    Jc(3,i) = top[0];
    Jc(4,i) = top[1];
    Jc(5,i) = top[2];
  }

  return Jc;
}

/// Resets the force and torque accumulators for all links and joints in the rigid body
void RCArticulatedBody::reset_accumulators()
{
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->reset_accumulators();

  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->reset_force();
}

/// The signum function
double RCArticulatedBody::sgn(double x)
{
  if (x < -NEAR_ZERO)
    return (double) -1.0;
  else if (x > NEAR_ZERO)
    return (double) 1.0;
  else
    return (double) 0.0;
}

/// Computes the forward dynamics
/**
 * Given the joint positions and velocities, joint forces, and external
 * forces on the links, determines the joint and link accelerations as well as
 * floating base accelerations (if applicable).  The joint
 * accelerations are stored in the individual joints and the link accelerations
 * are stored in the individual links.
 * \note only computes the forward dynamics if the state-derivative is no longer valid
 */
void RCArticulatedBody::calc_fwd_dyn()
{
  SAFESTATIC VectorNd ff;

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() entered" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  computing forward dynamics in ";
  if (get_computation_frame_type() == eGlobal)
    FILE_LOG(LOG_DYNAMICS) << "global ";
  else if (get_computation_frame_type() == eLink)
    FILE_LOG(LOG_DYNAMICS) << "link ";
  else if (get_computation_frame_type() == eLinkCOM)
    FILE_LOG(LOG_DYNAMICS) << "link c.o.m. ";
  else
    FILE_LOG(LOG_DYNAMICS) << "joint ";
  FILE_LOG(LOG_DYNAMICS) << "coordinate system" << std::endl;

  // do forward dynamics with old Coulomb model
  for (unsigned i=0; i< _joints.size(); i++)
  {
    ff.resize(_joints[i]->num_dof());
    for (unsigned j=0; j< _joints[i]->num_dof(); j++)
    {
      ff[j] = (double) -1.0;
      if (_joints[i]->qd[j] < (double) 0.0)
        ff[j] = (double) 1.0;
      ff[j] *= _joints[i]->mu_fc;
      ff[j] -= _joints[i]->qd[j]*_joints[i]->mu_fv;
    }
    _joints[i]->add_force(ff);
  }

  // if there are implicit joints, we must do the model with loops
  if (!_ijoints.empty())
    calc_fwd_dyn_loops();
  else
  {
    // use the proper dynamics algorithm
    switch (algorithm_type)
    {
      case eFeatherstone:
        if (!_position_invalidated)
          _fsab.calc_fwd_dyn_special();
        else
          _fsab.calc_fwd_dyn();
        break;

      case eCRB:
        if (!_position_invalidated)
          _crb.calc_fwd_dyn_special();
        else
          _crb.calc_fwd_dyn();
        break;

      default:
        assert(false);
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() exited" << std::endl;
}

/// Computes the forward dynamics with loops
void RCArticulatedBody::calc_fwd_dyn_loops()
{
  double Cx[6];
  SAFESTATIC MatrixNd Jx_iM_JxT, U, V;
  SAFESTATIC MatrixNd Jx_dot, iM_JxT;
  SAFESTATIC VectorNd v, fext, C, alpha_x, beta_x, Dx_v, Jx_v, Jx_dot_v, workv;
  SAFESTATIC VectorNd iM_fext, a, S;

  // get the generalized velocity, generalized forces, and inverse generalized
  // inertia matrix
  get_generalized_velocity(eSpatial, v);
  get_generalized_forces(fext);

  // determine how many implicit constraint equations
  unsigned N_EXPLICIT_CONSTRAINT_EQNS = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    N_EXPLICIT_CONSTRAINT_EQNS = _ijoints[i]->num_constraint_eqns();

  // evaluate implicit constraints
  C.resize(N_EXPLICIT_CONSTRAINT_EQNS);
  for (unsigned i=0, r=0, s=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->evaluate_constraints(Cx);
    for (unsigned j=0; j< _ijoints[i]->num_constraint_eqns(); j++)
      C[r++] = Cx[j];
  }

  // get the implicit constraint Jacobian and its time derivative
  determine_implicit_constraint_jacobian(_Jx);
  determine_implicit_constraint_jacobian_dot(Jx_dot);

  // compute Jx * v and Jx_dot * v
  _Jx.mult(v, Jx_v);
  Jx_dot.mult(v, Jx_dot_v);

  // get movement Jacobian for implicit constraints and compute velocities
  determine_implicit_constraint_movement_jacobian(_Dx);
  _Dx.mult(v, Dx_v);
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    Dx_v.get_sub_vec(k, k+_ijoints[i]->num_dof(), _ijoints[i]->qd);
    k += _ijoints[i]->num_dof();
  }

  // add in implicit actuator forces
  beta_x.resize(_Dx.rows());
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->get_scaled_force(workv);
    for (unsigned j=0; j< _ijoints[i]->num_dof(); j++)
      beta_x[k++] = workv[j];
  }
  _Dx.transpose_mult(beta_x, workv);
  fext += workv;

  // compute the constraint forces
  DynamicBody::solve_generalized_inertia(fext, iM_fext);
  _Jx.mult(iM_fext, alpha_x) += Jx_dot_v;
  _Jx.mult(v, workv) *= ((double) 2.0 * b_alpha);
  alpha_x += workv;
  C *= (b_beta*b_beta);
  alpha_x += C;
  DynamicBody::transpose_solve_generalized_inertia(_Jx, iM_JxT);
  _Jx.mult(iM_JxT, Jx_iM_JxT);
  _LA->svd(Jx_iM_JxT, U, S, V);
  _LA->solve_LS_fast(U, S, V, alpha_x);

  // compute generalized acceleration
  fext -= _Jx.transpose_mult(alpha_x, workv);
  DynamicBody::solve_generalized_inertia(fext, a);
  set_generalized_acceleration(a);
}

/// Determines the ndof x ngc Jacobian for implicit constraint movement (ndof is the number of degrees of freedom of the implicit constraints)
// TODO: fix this
void RCArticulatedBody::determine_implicit_constraint_jacobians(const UnilateralConstraintProblemData& q, MatrixNd& Jx, MatrixNd& Dx) const
{
/*
  SAFESTATIC vector<SVelocityd> so;
  SAFESTATIC MatrixNd sub, pinv_si;
  SAFESTATIC VectorNd vsub;
  double Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // determine the total number of implicit constraint DOF
  unsigned NDOF = 0, NEQ = 0;
  for (unsigned i=0; i< q.constraint_events.size(); i++)
  {
    JointPtr joint = q.constraint_events[i]->constraint_joint;
    if (joint->get_constraint_type() == Joint::eImplicit)
    {
      NDOF += joint->num_dof();
      NEQ += joint->num_constraint_eqns();
    }
  }

  // resize Jx and Dx
  Jx.set_zero(NEQ, NGC);
  Dx.set_zero(NDOF, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the Jacobian for all implicit constraints
  for (unsigned i=0, ii=0, jj=0; i< q.constraint_events.size(); i++)
  {
    // get the joint
    JointPtr joint = q.constraint_events[i]->constraint_joint;

    // if it's not implicit, skip it
    if (joint->get_constraint_type() != Joint::eImplicit)
      continue;

    // it is implicit; see whether it belongs to this body
    if (joint->get_articulated_body() != get_this())
    {
      ii += joint->num_dof();
      jj += joint->num_constraint_eqns();
      continue;
    }

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = joint->get_inboard_link();
    RigidBodyPtr rbo = joint->get_outboard_link();

    // get the outboard body's transform
    shared_ptr<const Pose3d> To = rbo->get_pose();

    // get the spatial movement (axes) of this joint (rbo frame)
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    pinv_si = si;
    try
    {
      LinAlg::pseudo_inverse(pinv_si, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      pinv_si = si;
      LinAlg::pseudo_inverse(pinv_si, LinAlg::svd2);
    }

    // special case: floating base
//    if (_floating_base)
//    {
//      SpatialTransform X(To, IDENTITY_3x3, base->get_position());
//      MatrixNd::transpose(pinv_si, sx);
//      X.transform(sx, so);
//      D.set_sub_mat(r, 0, so, true);
//    }


    // now, loop over explicit joints
    for (unsigned k=0; k< _ejoints.size(); k++)
    {
      // get the links
      RigidBodyPtr rba = _ejoints[k]->get_inboard_link();
      RigidBodyPtr rbb = _ejoints[k]->get_outboard_link();

      // get the coordinate index of the joint
      unsigned cidx = _ejoints[k]->get_coord_index();

      // get transform for the outboard link
      const Matrix4& rbbT = rbb->get_pose();

      // get the spatial axes for the joint (in rbb's frame)
      const vector<SVelocityd>& sk = _ejoints[k]->get_spatial_axes(eLink);

      // compute the spatial transformation to the outboard link
      SpatialTransform X(rbbT, To);
      X.transform(sk, so);
      pinv_si.mult(so, sub);

      // update Dx
      Dx.set_sub_mat(ii,cidx,sub);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, so);
//        so.transpose_mult(svo, dot);
//        sub -= dot;

      // iterate over the constraint equations
      for (unsigned j=0; j< joint->num_constraint_eqns(); j++, jj++)
      {
        // get constraint equations for inner link
        joint->calc_constraint_jacobian(DynamicBody::eSpatial, rbi, j, Cqi);

        // setup spatial vectors in rbi-centered frame and rbo-centered frame
        SVector6 svi(Cqi);

        // special case: floating base
        if (_floating_base)
        {
          // setup spatial vector in rbo-centered frame
          joint->calc_constraint_jacobian(DynamicBody::eSpatial,rbo, j, Cqo);
          SVector6 svo(Cqo);

          // setup spatial transform for svi/svo
          SpatialTransform X;
          X.r = base->get_position() - rbi->get_position();
          SVector6 svix = X.transform(svi);
          X.r = base->get_position() - rbo->get_position();
          SVector6 svox = X.transform(svo);

          // update Jx corresponding to the base
          for (unsigned k=0; k< 6; k++)
            Jx(jj,k) += svix[k] - svox[k];
        }

        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(sk, so);
        so.transpose_mult(svi, vsub);

        // update Jx
        Jx.set_sub_mat(jj,cidx,vsub,true);
      }
    }

    // update ii
    ii += joint->num_dof();
  }
*/
}

/// Determines the ndof x ngc Jacobian for implicit constraint movement (ndof is the number of degrees of freedom of the implicit constraints)
// TODO: fix this
void RCArticulatedBody::determine_implicit_constraint_movement_jacobian(MatrixNd& D)
{
/*
  SAFESTATIC vector<SVelocityd> so;
  SAFESTATIC MatrixNd sub, pinv_si;
  double Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // determine the number of implicit constraint DOF
  unsigned NDOF = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NDOF += _ijoints[i]->num_dof();

  // resize D
  D.set_zero(NDOF, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the Jacobian for all implicit constraints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the coordinate index of the joint
    unsigned r = _ijoints[i]->get_coord_index();

    // get the outboard body of this joint and its transform
    RigidBodyPtr rbo = _ijoints[i]->get_outboard_link();
    const Matrix4& To = rbo->get_pose();

    // get the spatial movement (axes) of this joint (rbo frame)
    const vector<SVelocityd>& si = _ijoints[i]->get_spatial_axes(eLink);
    pinv_si = si;
    try
    {
      LinAlg::pseudo_inverse(pinv_si, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      pinv_si = si;
      LinAlg::pseudo_inverse(pinv_si, LinAlg::svd2);
    }

    // special case: floating base
//    if (_floating_base)
//    {
//      SpatialTransform X(To, IDENTITY_3x3, base->get_position());
//      MatrixNd::transpose(pinv_si, sx);
//      X.transform(sx, so);
//      D.set_sub_mat(r, 0, so, true);
//    }

    // now, loop over implicit joints
    for (unsigned k=0; k< _ejoints.size(); k++)
    {
      // get the outboard link
      RigidBodyPtr rbb = _ejoints[k]->get_outboard_link();
      unsigned idx = _ejoints[k]->get_coord_index();
      const Matrix4& Tb = rbb->get_pose();

      // get the spatial axes for the joint (rbb frame)
      const vector<SVelocityd>& sk = _ejoints[k]->get_spatial_axes(eLink);

      // compute the spatial transformation
      SpatialTransform X(Tb, To);
      X.transform(sk, so);
      pinv_si.mult(so, sub);

      // update D
      D.set_sub_mat(r,idx,sub);
    }
  }
*/
}

/// Determines the constraint Jacobian for implicit constraints via impulse application (considerably slower than alternative method)
/*
void RCArticulatedBody::determine_implicit_constraint_jacobian(MatrixNd& J)
{
  double Cx[6];

  // determine the number of implicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NEQ += _ijoints[i]->num_constraint_eqns();
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // resize J
  J.set_zero(NEQ, NGC);
  MatrixNdN q(NGC);

  // init gj and dC
  VectorNd gj = VectorNd::zero(NGC);
  MatrixNd dC(NEQ, NGC);
  VectorNd gv(NGC);

  // loop
  for (unsigned i=0; i< NGC; i++)
  {
    // apply the generalized impulse
    if (i > 0)
      gj[i-1] = (double) 0.0;
    gj[i] = (double) 1.0;
    apply_generalized_impulse(DynamicBody::eSpatial, gj);

    // get the velocity
    get_generalized_velocity(DynamicBody::eSpatial, gv);
    q.set_column(i, gv);

    // compute the constraint dot
    for (unsigned r=0; r< _ijoints.size(); r++)
    {
      unsigned cidx = _ijoints[r]->get_constraint_index();
      _ijoints[r]->evaluate_constraints_dot(Cx);
      for (unsigned j=0; j< _ijoints[r]->num_constraint_eqns(); j++)
        dC(j+cidx, i) = Cx[j];
    }
  }

  // compute J
  MatrixNdN qinv;
  LinAlg::pseudo_inverse(q, qinv);
  dC.mult(qinv, J);

  // remove all impulses
  gj.set_one().negate();
  apply_generalized_impulse(DynamicBody::eSpatial, gj);
}
*/

/// Determines the constraint Jacobian for implicit constraints
// TODO: fix this
void RCArticulatedBody::determine_implicit_constraint_jacobian(MatrixNd& J)
{
/*
  SAFESTATIC vector<SVelocityd> so;
  SAFESTATIC VectorNd sub;
  double Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // determine the number of implicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NEQ += _ijoints[i]->num_constraint_eqns();

  // resize J
  J.set_zero(NEQ, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the constraint Jacobian for all implicit constraints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the constraint index for this joint
    unsigned ridx = _ijoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ijoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ijoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ijoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ijoints[i]->calc_constraint_jacobian(DynamicBody::eSpatial, rbi, j, Cqi);
      _ijoints[i]->calc_constraint_jacobian(DynamicBody::eSpatial, rbo, j, Cqo);

      // setup spatial vectors in rbi-centered frame and rbo-centered frame
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++)
          J(ridx,k) += svix[k] - svox[k];
      }

      // now, loop over explicit joints
      for (unsigned k=0; k< _ejoints.size(); k++)
      {
        // get the links
        RigidBodyPtr rba = _ejoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ejoints[k]->get_outboard_link();

        // get the coordinate index of the joint
        unsigned cidx = _ejoints[k]->get_coord_index();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_pose();

        // get the spatial axes for the joint (in rbb's frame)
        const vector<SVelocityd>& si = _ejoints[k]->get_spatial_axes(eLink);

        // setup components corresponding to rbi
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si, so);
        so.transpose_mult(svi, sub);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, so);
//        so.transpose_mult(svo, dot);
//        sub -= dot;

        // update J
        J.set_sub_mat(ridx,cidx,sub,true);
      }
    }
  }
*/
}

/// Determines the time derivative of the constraint Jacobian for implicit constraints
/**
 * Because this matrix is the product of two matrices, each of which is a
 * function of time, we have to add the results together.
 */
// TODO: fix this
void RCArticulatedBody::determine_implicit_constraint_jacobian_dot(MatrixNd& J) const
{
/*
  SAFESTATIC vector<SVelocityd> so;
  SAFESTATIC VectorNd sub;
  double Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // determine the number of implicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NEQ += _ijoints[i]->num_constraint_eqns();

  // resize J
  J.set_zero(NEQ, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // Jx * \dot{Jy}
  // compute the constraint Jacobian for all implicit constraints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the starting constraint index for this joint
    unsigned ridx = _ijoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ijoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ijoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ijoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ijoints[i]->calc_constraint_jacobian(DynamicBody::eSpatial, rbi, j, Cqi);
      _ijoints[i]->calc_constraint_jacobian(DynamicBody::eSpatial, rbo, j, Cqo);

      // setup spatial vectors
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++)
          J(ridx,k) += svix[k] - svox[k];
      }

      // now, loop over explicit joints
      for (unsigned k=0; k< _ejoints.size(); k++)
      {
        // get the coordinate index for the joint
        unsigned cidx = _ejoints[k]->get_coord_index();

        // get the links
        RigidBodyPtr rba = _ejoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ejoints[k]->get_outboard_link();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_pose();

        // get the time derivative of the spatial axes for the joint
        const vector<SVelocityd>& si_dot = _ejoints[k]->get_spatial_axes_dot(eLink);

        // setup components corresponding to rbi
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si_dot, so);
        so.transpose_mult(svi, sub);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si_dot, so);
//        so.transpose_mult(svo, dot);
//        sub -= dot;

        // update J
        J.set_sub_mat(ridx,cidx,sub,true);
      }
    }
  }

  // \dot{Jx} * Jy
  // compute the constraint Jacobian for all implicit constraints
  for (unsigned i=0, r=0; i< _ijoints.size(); i++)
  {
    // get the starting constraint index for this joint
    unsigned ridx = _ijoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ijoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ijoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ijoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ijoints[i]->calc_constraint_jacobian_dot(DynamicBody::eSpatial, rbi, j, Cqi);
      _ijoints[i]->calc_constraint_jacobian_dot(DynamicBody::eSpatial, rbo, j, Cqo);

      // setup spatial vectors
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++)
          J(ridx,k) += svix[k] - svox[k];
      }

      // now, loop over explicit joints
      for (unsigned k=0; k< _ejoints.size(); k++)
      {
        // get the coordinate index
        unsigned cidx = _ejoints[k]->get_coord_index();

        // get the links
        RigidBodyPtr rba = _ejoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ejoints[k]->get_outboard_link();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_pose();

        // get the spatial axes for the joint
        const vector<SVelocityd>& si = _ejoints[k]->get_spatial_axes(eLink);

        // setup components corresponding to rbi
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si, so);
        so.transpose_mult(svi, sub);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, so);
//        so.transpose_mult(svo, dot);
//        sub -= dot;

        // update J
        for (unsigned b=0; b< sub.size(); b++)
            J(ridx,cidx+b) += sub[b];
      }
    }
  }
*/
}

/// Sets the generalized acceleration for this body
void RCArticulatedBody::set_generalized_acceleration(const VectorNd& a)
{
  SAFESTATIC VectorNd base_a;

  if (_floating_base)
  {
    a.get_sub_vec(num_joint_dof_explicit(), a.size(), base_a);
    RigidBodyPtr base = _links.front();
    SAcceld xdd;
    xdd.set_linear(Vector3d(base_a[0], base_a[1], base_a[2]));
    xdd.set_angular(Vector3d(base_a[3], base_a[4], base_a[5]));
    base->set_accel(xdd);
  }

  // set joint accelerations
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    a.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->qdd);
  }
}

/// Applies an impulsive force at the specified point and propagates it up the articulated body chain
/**
 * \param w the impulse
 * \param link the link that the impulse is applied to
 * \param fsab_algo if non-null, already computed quanities from dynamics
 *         algorithm will be used; otherwise, forward dynamics algorithm will
 *         be called using FSAB algorithm
 */
void RCArticulatedBody::apply_impulse(const SMomentumd& w, RigidBodyPtr link)
{
  // compute the forward dynamics, given the algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      _fsab.apply_impulse(w, link);
      break;

    case eCRB:
      _crb.apply_impulse(w, link);
      break;

    default:
      assert(false);
  }
}

/// Synchronizes the visualization for all links and joints of this body to their corresponding transforms
void RCArticulatedBody::update_visualization()
{
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->update_visualization();
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->update_visualization();
}

/// Gets the generalized coordinates of this body
SharedVectorNd& RCArticulatedBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, SharedVectorNd& gc)
{
  get_generalized_coordinates_generic(gctype, gc);
  return gc;
}

/// Sets the generalized position of this body
void RCArticulatedBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const SharedVectorNd& gc)
{
  set_generalized_coordinates_generic(gctype, gc);
}

/// Sets the generalized velocity of this body
void RCArticulatedBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const SharedVectorNd& gv)
{
  set_generalized_velocity_generic(gctype, gv);
}

/// Gets the generalized velocity of this body
SharedVectorNd& RCArticulatedBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, SharedVectorNd& gv)
{
  get_generalized_velocity_generic(gctype, gv);
  return gv;
}

/// Gets the generalized inertia of this body
SharedMatrixNd& RCArticulatedBody::get_generalized_inertia(SharedMatrixNd& M)
{
  // calculate the generalized inertia matrix
  _crb.calc_generalized_inertia(M);

  return M;
}

/// Gets the number of degrees of freedom for implicit constraints
unsigned RCArticulatedBody::num_joint_dof_implicit() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    ndof += _ijoints[i]->num_dof();
  return ndof;
}

/// Gets the generalized forces on this body
/**
 * \note does not add forces from implicit joint constraints!
 */
SharedVectorNd& RCArticulatedBody::get_generalized_forces(SharedVectorNd& f)
{
  const unsigned SPATIAL_DIM = 6, X = 0, Y = 1, Z = 2, A = 3, B = 4, C = 5;

  // compute the generalized forces
  SForced f0;
  VectorNd CmQ;
  _crb.calc_generalized_forces(f0, CmQ);

  // determine the vector of joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    for (unsigned m=0; m< _ejoints[i]->num_dof(); m++, idx++)
      f[idx] = _ejoints[i]->force[m] - CmQ[idx];
  }

  // setup joint space part of f
  if (_floating_base)
  {
    // determine external and inertial forces on base
    RigidBodyPtr base = _links.front();
    unsigned idx = CmQ.size();
    SForced fb = Pose3d::transform(base->get_gc_pose(), f0);
    f[idx++] = -fb[X];  f[idx++] = -fb[Y];  f[idx++] = -fb[Z];
    f[idx++] = -fb[A];  f[idx++] = -fb[B];  f[idx++] = -fb[C];
  }

  return f;
}

/// Converts a force to a generalized force
SharedVectorNd& RCArticulatedBody::convert_to_generalized_force(SingleBodyPtr body, const SForced& w, SharedVectorNd& gf)
{
  const unsigned SPATIAL_DIM = 6;
  static vector<SVelocityd> J;
  static vector<SVelocityd> sprime;

  // get the body as a rigid body
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  assert(link);

  // get the gc frame
  shared_ptr<const Pose3d> P = _links.front()->get_gc_pose();

  // get w in the mixed frame
  SForced wP = Pose3d::transform(P, w);

  // clear the Jacobian
  J.resize(num_joint_dof_explicit());
  for (unsigned i=0; i< J.size(); i++)
    J[i] = SVelocityd::zero(P);

  // compute the Jacobian in w's frame
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // if link is not a descendent of this joint's outboard, keep looping
    RigidBodyPtr outboard = _ejoints[i]->get_outboard_link();
    if (!outboard->is_descendant_link(link))
      continue;

    // transform the Jacobian
    const vector<SVelocityd>& s = _ejoints[i]->get_spatial_axes();
    Pose3d::transform(P, s, sprime);
    for (unsigned j=0, k=_ejoints[i]->get_coord_index(); j < s.size(); j++, k++)
      J[k] = sprime[j];
  }

  // resize gf
  gf.resize(num_generalized_coordinates(DynamicBody::eSpatial));

  // get the torque on the joints
  SharedVectorNd jf = gf.segment(0, J.size());
  transpose_mult(J, wP, jf);

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::convert_to_generalized_force() - converting " << std::endl << "  " << w << " to generalized force" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  -- joint torques: " << gf.segment(0, J.size()) << std::endl;

  // determine the generalized force on the base, if the base is floating
  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    SharedVectorNd gfbase = gf.segment(J.size(), gf.size());
    wP.to_vector(gfbase);
  }

  // return the generalized force vector
  return gf;
}

/// Determines whether a joint supports a link
bool RCArticulatedBody::supports(JointPtr joint, RigidBodyPtr link)
{
  // save the original link
  RigidBodyPtr l = link;

  do
  {
    JointPtr j = l->get_inner_joint_explicit();
    // check for l is base
    if (!j)
      return false;
    if (j == joint)
      return true;

    // proceed up the chain
    l = j->get_inboard_link();
  }
  while (link != l);

  return false;
}

/// Implements Base::load_from_xml()
/**
 * \pre all links and joints must be loaded using their respective
 *      serialization methods before this method is called
 */
void RCArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // load the parent data
  ArticulatedBody::load_from_xml(node, id_map);

  // don't verify the node name -- this class has derived classes
//  assert(strcasecmp(node->name.c_str(), "RCArticulatedBody") == 0);

  // get whether the body has a floating base
  XMLAttrib* fb_attr = node->get_attrib("floating-base");
  if (fb_attr)
    _floating_base = fb_attr->get_bool_value();

  // read the pointer to the forward dynamics algorithm, if provided
  XMLAttrib* fdyn_algo_attr = node->get_attrib("fdyn-algorithm");
  if (fdyn_algo_attr)
  {
    // get the ID
    string algo = fdyn_algo_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = algo.find_first_not_of(" \t\n\r");
    size_t last_nws_index = algo.find_last_not_of(" \t\n\r");
    algo = algo.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the algorithm type
    if (strcasecmp(algo.c_str(), "fsab") == 0)
      algorithm_type = eFeatherstone;
    else if (strcasecmp(algo.c_str(), "crb") == 0)
      algorithm_type = eCRB;
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "forward dynamics algorithm " << std::endl << "  type '";
      std::cerr << algo << "' -- valid types are 'fsab' and 'crb'";
      std::cerr << std::endl;
    }
  }

  // read the forward dynamics algorithm computation frame, if provided
  XMLAttrib* fdyn_frame_attr = node->get_attrib("fdyn-algorithm-frame");
  if (fdyn_frame_attr)
  {
    // get the computation reference frame as a string
    string frame = fdyn_frame_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = frame.find_first_not_of(" \t\n\r");
    size_t last_nws_index = frame.find_last_not_of(" \t\n\r");
    frame = frame.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the frame type
    if (strcasecmp(frame.c_str(), "global") == 0)
      set_computation_frame_type(eGlobal);
    else if (strcasecmp(frame.c_str(), "link") == 0)
      set_computation_frame_type(eLink);
    else if (strcasecmp(frame.c_str(), "linkinertia") == 0 || strcasecmp(frame.c_str(), "link-inertia") == 0)
      set_computation_frame_type(eLinkInertia);
    else if (strcasecmp(frame.c_str(), "linkcom") == 0 || strcasecmp(frame.c_str(), "link-com") == 0)
      set_computation_frame_type(eLinkCOM);
    else if (strcasecmp(frame.c_str(), "joint") == 0)
      set_computation_frame_type(eJoint);
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "computation reference " << std::endl << "  frame type '";
      std::cerr << frame << "' -- valid types are 'link' and 'global'";
      std::cerr << std::endl;
    }
  }

  // get baumgarte paramters
  XMLAttrib* balpha_attr = node->get_attrib("baumgarte-alpha");
  if (balpha_attr)
    b_alpha = balpha_attr->get_real_value();
  XMLAttrib* bbeta_attr = node->get_attrib("baumgarte-beta");
  if (bbeta_attr)
    b_beta = bbeta_attr->get_real_value();

  // compile everything once again, for safe measure
  compile();

  // transform the body, if desired
  XMLAttrib* xlat_attr = node->get_attrib("translate");
  XMLAttrib* rotate_attr = node->get_attrib("rotate");

  if (rotate_attr)
  {
    RigidBodyPtr base = get_base_link();
    base->rotate(rotate_attr->get_rpy_value());
    update_link_poses();
  }
  if (xlat_attr)
  {
    RigidBodyPtr base = get_base_link();
    base->translate(xlat_attr->get_origin_value());
    update_link_poses();
  }
}

/// Implements Base::save_to_xml()
void RCArticulatedBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent method first
  ArticulatedBody::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "RCArticulatedBody";

  // write whether body has a floating base
  node->attribs.insert(XMLAttrib("floating-base", _floating_base));

  // write the forward dynamics algorithm type
  if (algorithm_type == eFeatherstone)
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("fsab")));
  else
  {
    assert(algorithm_type == eCRB);
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("crb")));
  }

  // write the forward dynamics algorithm frame -- note that the string()
  // is necessary on the second argument to XMLAttrib b/c the compiler
  // interprets a constant string as a bool, rather than as an string,
  // given a choice
  if (get_computation_frame_type() == eGlobal)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("global")));
  else if (get_computation_frame_type() == eLink)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("link")));
  else if (get_computation_frame_type() == eLinkCOM)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("linkcom")));
  else if (get_computation_frame_type() == eLinkInertia)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("linkinertia")));
  else
  {
    assert(get_computation_frame_type() == eJoint);
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("joint")));
  }

  // save baumgarte parameters
  node->attribs.insert(XMLAttrib("baumgarte-alpha", b_alpha));
  node->attribs.insert(XMLAttrib("baumgarte-beta", b_beta));
}

