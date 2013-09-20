/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
 
  // set default algorithm to FSAB and computation frame to global
  algorithm_type = eFeatherstone;
  set_computation_frame_type(eGlobal);

  // setup baumgarte parameters
  b_alpha = (double) 0.0;
  b_beta = (double) 0.0;

  // mark inverse inertia as invalid
  _fM_valid = false;
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
  // look for trivial cases
  if (_links.empty() || _joints.empty())
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
  if (_fM_valid)
    return;

  // we do need to update; mark M as non rank-deficient for now
  _M_rankdef = false;

  // indicate factorized inertia is valid    
  _fM_valid = true;
}

/// Solves using a generalized inertia matrix
VectorNd& RCArticulatedBody::solve_generalized_inertia(const VectorNd& v, VectorNd& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia();

  // make x/b one vector
  result = v;

  if (algorithm_type == eFeatherstone || _M_rankdef)
    _LA->solve_LS_fast(_uM, _sM, _vM, result);
  else
  {
    assert(algorithm_type == eCRB);
    _crb.M_solve(result);
  }
  
  return result;
}

/// Solves the transpose using a generalized inertia matrix
MatrixNd& RCArticulatedBody::transpose_solve_generalized_inertia(const MatrixNd& m, MatrixNd& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia();

  // setup the result
  MatrixNd::transpose(m, result);

  if (algorithm_type == eFeatherstone || _M_rankdef)
    _LA->solve_LS_fast(_uM, _sM, _vM, result);
  else
  {
    assert(algorithm_type == eCRB);
    _crb.M_solve(result);
  }
  
  return result;
}

/// Solves using a generalized inertia matrix
MatrixNd& RCArticulatedBody::solve_generalized_inertia(const MatrixNd& m, MatrixNd& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia();

  // setup the result
  result = m;

  if (algorithm_type == eFeatherstone || _M_rankdef)
    _LA->solve_LS_fast(_uM, _sM, _vM, result);
  else
  {
    assert(algorithm_type == eCRB);
    _crb.M_solve(result);
  }
  
  return result;
}

/// Solves using the transpose of the generalized inertia matrix
MatrixNd& RCArticulatedBody::solve_generalized_inertia_transpose(const MatrixNd& m, MatrixNd& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia();

  if (algorithm_type == eFeatherstone || _M_rankdef)
  {
    MatrixNd::transpose(m, result);
    _LA->solve_LS_fast(_uM, _sM, _vM, result);
  }
  else
  {
    assert(algorithm_type == eCRB);
    MatrixNd::transpose(m, result);
    _crb.M_solve(result);
  }
  
  return result;
}

/// Applies a generalized impulse to the articulated body
void RCArticulatedBody::apply_generalized_impulse(const VectorNd& gj)
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
    solve_generalized_inertia(gj, gv_delta);

    // apply the change in generalized velocity
    gv += gv_delta;
    set_generalized_velocity(DynamicBody::eSpatial, gv);
  }

  // reset the force and torque accumulators
  reset_accumulators();
}

/// Adds a generalized force to the articulated body
void RCArticulatedBody::add_generalized_force(const VectorNd& gf)
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
    base->add_force(fx);
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
VectorNd& RCArticulatedBody::get_generalized_acceleration(VectorNd& ga)
{
  const unsigned GC_AA_DIM = 6, GC_ROD_DIM = 7;
  SAFESTATIC VectorNd base_ga;

  // resize the state-derivative vector
  ga.resize(num_generalized_coordinates(DynamicBody::eSpatial));
  
  // setup the generalized acceleration for the base (if any) 
  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    base->get_generalized_acceleration_single(base_ga);
    ga.set_sub_vec(num_joint_dof_explicit(), base_ga);
  }

  // setup the state for the joints
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    ga.set_sub_vec(idx, _ejoints[i]->qdd);
  }

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
MatrixNd& RCArticulatedBody::calc_jacobian_floating_base(const Point3d& point, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC vector<SVelocityd> sbase, sbase_prime;
  SAFESTATIC shared_ptr<Pose3d> P(new Pose3d);

  // get the base link and the base pose
  RigidBodyPtr base = get_base_link();
  shared_ptr<const Pose3d> baseP = base->get_pose();

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

  if (is_floating_base())
  {
    // calculate the floating base
    calc_jacobian_floating_base(p, Jsub);

    // setup the floating base
    J.set_sub_mat(0,0,Jsub);
  }

  // calculate all relevant columns
  RigidBodyPtr base = get_base_link();
  while (link != base)
  {
    JointPtr joint = link->get_inner_joint_explicit();
    calc_jacobian_column(joint, p, Jsub); 
    J.set_sub_mat(0,joint->get_coord_index(), Jsub);
    link = link->get_parent_link();
  }

  return J;
}

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
  {
    calc_fwd_dyn_loops();
    return;
  } 

  // use the proper dynamics algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      _fsab.calc_fwd_dyn();
      break;

    case eCRB:
      _crb.calc_fwd_dyn();
      break;

    default:
      assert(false);
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
  solve_generalized_inertia(fext, iM_fext);
  _Jx.mult(iM_fext, alpha_x) += Jx_dot_v;
  _Jx.mult(v, workv) *= ((double) 2.0 * b_alpha);
  alpha_x += workv;
  C *= (b_beta*b_beta);
  alpha_x += C;
  solve_generalized_inertia_transpose(_Jx, iM_JxT);
  _Jx.mult(iM_JxT, Jx_iM_JxT);
  _LA->svd(Jx_iM_JxT, U, S, V);
  _LA->solve_LS_fast(U, S, V, alpha_x);

  // compute generalized acceleration
  fext -= _Jx.transpose_mult(alpha_x, workv);
  solve_generalized_inertia(fext, a);
  set_generalized_acceleration(a);
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
void RCArticulatedBody::calc_fwd_dyn_advanced_friction(double dt)
{
/*
  double Cx[6];
  const double QP_TOL = std::sqrt(NEAR_ZERO);  // tolerance for solving the QP
  SAFESTATIC MatrixNd iM, X, Jx_X_JxT, Dx_X_DxT, Jx_iM_JxT;
  SAFESTATIC MatrixNd Jx_dot, Y, X_JxT, X_DxT, Dx_X_JxT, Jx_iM, Jx_Y, A, RG;
  SAFESTATIC MatrixNd Jx_iM_DxT;
  SAFESTATIC VectorNd v, fext, C, alpha_x, beta_x, Dx_v, Jx_v, Jx_dot_v, workv;
  SAFESTATIC VectorNd vddt_plus_Y_fext, Y_fext, x, ff, delta, iM_fext, a;

  // get the generalized velocity, generalized forces, and inverse generalized 
  // inertia matrix
  get_generalized_velocity(eSpatial, v);
  get_generalized_forces(eSpatial, fext);
  get_generalized_inertia_inverse(eSpatial, iM);

  // setup X and Y
  const unsigned SPATIAL_DIM = 6;
  const unsigned START_GC = (_floating_base) ? SPATIAL_DIM : 0;
  iM.get_sub_mat(START_GC, iM.rows(), START_GC, iM.columns(), X);
  iM.get_sub_mat(0, iM.rows(), START_GC, iM.columns(), Y);

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

  // determine the number of kinematic loops and number of implicit joint DOF
  const unsigned N_LOOPS = _ijoints.size();
  unsigned N_EXPLICIT_DOF = 0, N_IMPLICIT_DOF = 0;
  for (unsigned i=0; i< _joints.size(); i++)
  {
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      N_EXPLICIT_DOF += _joints[i]->num_dof();
    else if (_joints[i]->get_constraint_type() == Joint::eExplicit)
      N_IMPLICIT_DOF += _joints[i]->num_dof();
    else
      assert(false);
  }

  // setup some constants
  const unsigned N_JOINT_DOF = N_IMPLICIT_DOF + N_EXPLICIT_DOF;
  const unsigned FF_START = 0;
  const unsigned ALPHAX_START = N_IMPLICIT_DOF;
  const unsigned BETA_START = N_IMPLICIT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;

  // precompute some things
  X.mult_transpose(_Jx, X_JxT);
  X.mult_transpose(_Dx, X_DxT);
  _Jx.mult(X_JxT, Jx_X_JxT);
  _Dx.mult(X_JxT, Dx_X_JxT);
  _Dx.mult(X_DxT, Dx_X_DxT);

  // setup data needed for convex optimization
  ABFwdDynOptData copt_data;
  copt_data.N_IMPLICIT_DOF = N_IMPLICIT_DOF;
  copt_data.N_EXPLICIT_CONSTRAINT_EQNS = N_EXPLICIT_CONSTRAINT_EQNS;
  copt_data.N_JOINT_DOF = N_JOINT_DOF;
  copt_data.N_LOOPS = N_LOOPS;
  copt_data.Dx = _Dx;
  copt_data.fext = fext;

  // first compute the nullspace and the homogeneous solution
  const unsigned NN = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS + N_LOOPS;
  MatrixNd& R = copt_data.R;
  VectorNd& z = copt_data.z;
  iM.mult(fext, iM_fext);
  if (N_EXPLICIT_CONSTRAINT_EQNS > 0)
  {
    // setup the right hand side
    _Jx.mult(iM_fext, z) *= dt;
    z += Jx_dot_v;
    _Jx.mult(v, workv) *= ((double) 2.0 * b_alpha);
    z += workv;
    C *= (b_beta*b_beta);
    z += C;

    // compute the nullspace
    _Jx.mult(Y, Jx_Y);
    _Jx.mult(iM, Jx_iM);
    Jx_iM.mult_transpose(_Jx, Jx_iM_JxT);
    Jx_iM.mult_transpose(_Dx, Jx_iM_DxT);
    A.resize(N_EXPLICIT_CONSTRAINT_EQNS, N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS);
    A.set_sub_mat(0, 0, Jx_Y);
    A.set_sub_mat(0, N_IMPLICIT_DOF, Jx_iM_JxT);
    A.set_sub_mat(0, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, Jx_iM_DxT);
    A *= dt;
    LinAlg::nullspace(A, RG);

    // now add in parts for delta (unrelated to nullspace)
    R.set_zero(NN, RG.columns()+N_LOOPS);
    R.set_sub_mat(0,0,RG);
    for (unsigned i=0; i< N_LOOPS; i++)
      R(N_JOINT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+i,RG.columns()+i) = (double) 1.0;

    // solve for implicit constraint forces 
    A = Jx_iM_JxT;
    A *= dt;
    workv = z;
    _LA->svd(A, U, S, V);
    _LA->solve_LS_fast(U, S, V, workv);
    z.set_zero(NN);
    z.set_sub_vec(ALPHAX_START, workv);
  }
  else
  {
    z.set_zero(N_JOINT_DOF);
    R.set_zero(N_JOINT_DOF, N_JOINT_DOF);
    for (unsigned i=0; i< N_JOINT_DOF; i++)
      R(i,i) = (double) 1.0;
  }

  // setup components of z and R to make things faster for gradient and Hessian
  // calculations
  z.get_sub_vec(FF_START, FF_START+N_IMPLICIT_DOF, copt_data.zff);
  z.get_sub_vec(BETA_START, BETA_START+N_EXPLICIT_DOF, workv);
  _Dx.transpose_mult(workv, copt_data.zbetax);
  R.get_sub_mat(FF_START, FF_START+N_IMPLICIT_DOF, 0, R.columns(), copt_data.Rff);
  R.get_sub_mat(BETA_START, BETA_START+N_EXPLICIT_DOF, 0, R.columns(), A);
  _Dx.transpose_mult(A, copt_data.DxTRbetax);

  // setup the data for the optimization problem
  const unsigned N = R.columns(); 
  const unsigned M = N_JOINT_DOF + N_LOOPS*2;
  OptParams copt_params(N, M, 0, calc_fwd_dyn_f0, calc_fwd_dyn_fx, NULL, calc_fwd_dyn_grad0, calc_fwd_dyn_cJac, NULL, calc_fwd_dyn_hess);
  copt_params.max_iterations = N*N*N;
  copt_params.data = (void*) &copt_data;

  // determine loops and compute Z matrices 
  vector<unsigned>& loop_indices = copt_data.loop_indices;
  vector<vector<unsigned> > loop_links;
  ArticulatedBody::find_loops(loop_indices, loop_links); 
  ArticulatedBody::compute_Z_matrices(loop_indices, loop_links, copt_data.Zd, copt_data.Z1d, copt_data.Z);

  // verify that joint indices are as we expect
  #ifndef NDEBUG
  for (unsigned i=0; i< _joints.size(); i++)
    assert(_joints[i]->get_index() == i);
  #endif

  // setup true indices: converts from an [explicit implicit] DOF to the
  // joint's index that contains that DOF
  vector<unsigned>& true_indices = copt_data.true_indices;
  true_indices.resize(N_JOINT_DOF);
  for (unsigned i=0, ki=0, ke=0; i< _joints.size(); i++)
  {
    if (_joints[i]->get_constraint_type() == Joint::eExplicit)
     for (unsigned j=0; j< _joints[i]->num_dof(); j++)
       true_indices[ki++] = i;
    else
    {
      assert(_joints[i]->get_constraint_type() == Joint::eImplicit);
      for (unsigned j=0; j< _joints[i]->num_dof(); j++)
        true_indices[N_IMPLICIT_DOF + ke++] = i;
    }
  }

  // initialize mu_c and viscous force vectors
  vector<double>& mu_c = copt_data.mu_c;
  vector<double>& visc = copt_data.visc;
  mu_c.resize(N_JOINT_DOF);
  visc.resize(N_JOINT_DOF);

  // setup mu_c and viscous force vectors for explicit joints
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
    for (unsigned j=0; j< _ejoints[i]->num_dof(); j++, k++)
    {
      mu_c[k] = _ejoints[i]->mu_fc*_ejoints[i]->mu_fc;
      double tmp = _ejoints[i]->mu_fv * std::fabs(_ejoints[i]->qd[j]);
      visc[k] = tmp*tmp;
    }

  // setup mu_c and viscous force vectors for implicit joints
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
    for (unsigned j=0; j< _ijoints[i]->num_dof(); j++, k++)
    {
      mu_c[k+N_IMPLICIT_DOF] = _ijoints[i]->mu_fc * _ijoints[i]->mu_fc;
      double tmp = _ijoints[i]->mu_fv * std::fabs(_ijoints[i]->qd[j]);
      visc[k+N_IMPLICIT_DOF] = tmp*tmp; 
    }

  // original (w/o nullspace) optimization problem is:
  // 0.5 * (v + h*inv(M)*tau)' * M * (v + h*inv(M)*tau)
  // ==> 0.5 * h^2 tau' inv(M) tau + v' * h * tau + 0.5 * v' * M * v'
  // ==> 0.5 v' M v +
  //     0.5 (h^2 fext' inv(M) fext + h^2 fext' inv(M) ff + 
  //     h^2 fext' inv(M) Jx'ax + h^2 fext' inv(M) Dx'bx) + v' * h * fext +
  //     0.5 (h^2 ff' inv(M) fext + h^2 ff' inv(M) ff + h^2 ff' inv(M) Jx' ax +
  //          h^2 ff' inv(M) Dx' bx) + v' * h * ff +
  //     0.5 (h^2 ax'Jx inv(M) fext + h^2 ax'Jx inv(M) ff + 
  //          h^2 ax'Jx inv(M) Jx' ax + h^2 ax'Jx inv(M) Dx' bx) + 
  //          v' * h * Jx' * ax +
  //     0.5 (h^2 bx'Dx inv(M) fext + h^2 bx'Dx inv(M) ff + 
  //          h^2 bx'Dx inv(M) Jx' ax + h^2 bx'Dx inv(M) Dx' bx) + 
  //          v' * h * Dx' * bx
  // G = | inv(M)     inv(M) Jx'     inv(M) Dx'    0 |
  //     | Jx inv(M)  Jx inv(M) Jx'  Jx inv(M) Dx' 0 |
  //     | Dx inv(M)  Dx inv(M) Jx'  Dx inv(M) Dx' 0 |
  //     | 0          0              0             0 |
  // c = | inv(M) fext + v/h       |
  //     | Jx inv(M) fext + Jx*v/h |
  //     | Dx inv(M) fext + Dx*v/h |
  //     | 0                       |

  // compute the quadratic matrix
  MatrixNd& G = copt_data.G;
  G.set_zero(NN, NN);
  G.set_sub_mat(0,0,X); 
  G.set_sub_mat(0,N_IMPLICIT_DOF,X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF,0,X_JxT,true);
  G.set_sub_mat(0,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,X_DxT);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,0,X_DxT,true);
  G.set_sub_mat(N_IMPLICIT_DOF,N_IMPLICIT_DOF,Jx_X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,Dx_X_JxT, true);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,N_IMPLICIT_DOF,Dx_X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,Dx_X_DxT);
  FILE_LOG(LOG_DYNAMICS) << "G (before nullspace): " << std::endl << G;
  R.transpose_mult(G, RG);
  RG.mult(R, G);

  // (z + Ry)G(z + Ry) + (z' + y'R')c = 0.5 y'R'G'Ry + y'(R'Gz + R'c)
  // compute the linear optimization vector
  VectorNd& c = copt_data.c;
  c.set_zero(NN);
  Y.mult(fext, Y_fext);
  (vddt_plus_Y_fext = v) /= dt;
  vddt_plus_Y_fext += Y_fext;
  c.set_sub_vec(0, vddt_plus_Y_fext);
  _Jx.mult(vddt_plus_Y_fext, workv);
  c.set_sub_vec(N_IMPLICIT_DOF, workv);
  _Dx.mult(vddt_plus_Y_fext, workv);
  c.set_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, workv);
  FILE_LOG(LOG_DYNAMICS) << "c (before nullspace): " << c << std::endl;
  R.transpose_mult(c, workv);
  c = workv;
  RG.mult(z, workv);
  c += workv;
  FILE_LOG(LOG_DYNAMICS) << "G (after nullspace): " << std::endl << G;
  FILE_LOG(LOG_DYNAMICS) << "c (after nullspace): " << std::endl << c;
  
  // setup an optimization vector
  x.set_zero(N);  

  // set delta to a middle value
  for (unsigned i=0,j=R.columns()-1; i< N_LOOPS; i++, j--)
    x[j] = (double) 0.5;

  // optimize
  if (G.norm_inf() > std::sqrt(NEAR_ZERO))
  {
//  Optimization::optimize_convex_pd(copt_params, x);
    Optimization::sqp(copt_params, x);
  }

  // output debugging data
  FILE_LOG(LOG_DYNAMICS) << "R: " << std::endl << R;
  FILE_LOG(LOG_DYNAMICS) << "z: " << z << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "x: " << x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "v: " << v << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "dt: " << dt << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "fext: " << fext << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "inv(M): " << std::endl << iM;
  FILE_LOG(LOG_DYNAMICS) << "Jx: " << std::endl << _Jx;
  FILE_LOG(LOG_DYNAMICS) << "\\dot{Jx}: " << std::endl << Jx_dot;
  FILE_LOG(LOG_DYNAMICS) << "Dx: " << std::endl << _Dx;

  // now retrieve the necessary variables
  R.mult(x, workv);
  z += workv;
  z.get_sub_vec(0, N_IMPLICIT_DOF, ff);
  z.get_sub_vec(N_IMPLICIT_DOF, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, alpha_x);
  z.get_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+N_EXPLICIT_DOF, beta_x);
  z.get_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+N_EXPLICIT_DOF, z.size(), delta);

  // output results
  FILE_LOG(LOG_DYNAMICS) << " explicit friction forces: " << ff << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " implicit constraint forces: " << alpha_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " implicit friction forces: " << beta_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " delta: " << delta << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " constraint evaluations: " << C << std::endl;

  // setup explicit joint friction forces
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    ff.get_sub_vec(k, k+_ejoints[i]->num_dof(), _ejoints[i]->ff);
    k += _ejoints[i]->num_dof();
  }

  // setup implicit joint friction forces
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    beta_x.get_sub_vec(k, k+_ijoints[i]->num_dof(), _ijoints[i]->ff);
    k += _ijoints[i]->num_dof();
  }

  // compute joint constraint forces
  fext += ff;
  fext += _Dx.transpose_mult(beta_x, workv);
  ArticulatedBody::calc_joint_constraint_forces(loop_indices, delta, copt_data.Zd, copt_data.Z1d, copt_data.Z, fext);

  // compute generalized acceleration
  fext -= _Jx.transpose_mult(alpha_x, workv);
  iM.mult(fext, a);
  set_generalized_acceleration(DynamicBody::eSpatial, a);
  if (LOGGING(LOG_DYNAMICS))
    for (unsigned i=0; i< _joints.size(); i++)
    {
      FILE_LOG(LOG_DYNAMICS) << "lambda for joint " << _joints[i]->id << ": " << _joints[i]->lambda << std::endl;
      FILE_LOG(LOG_DYNAMICS) << "Z: " << std::endl << copt_data.Z[i];
      FILE_LOG(LOG_DYNAMICS) << "Zd: " << std::endl << copt_data.Zd[i];
      FILE_LOG(LOG_DYNAMICS) << "Z1d: " << std::endl << copt_data.Z1d[i];
    }
  FILE_LOG(LOG_DYNAMICS) << " generalized acceleration: " << a << std::endl;

  // determine predicted energy
  if (LOGGING(LOG_DYNAMICS))
  {
    MatrixNd M;
    get_generalized_inertia(eSpatial, M);
    M.mult(v, workv);
    FILE_LOG(LOG_DYNAMICS) << " current energy: " << (0.5 * v.dot(workv)) << std::endl;
    a *= dt;
    v += a;
    M.mult(v, workv);
    FILE_LOG(LOG_DYNAMICS) << " new energy: " << (0.5 * v.dot(workv)) << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() exited" << std::endl;
*/
}

/// Determines the constraint Jacobian for implicit constraints via impulse application (considerably slower than alternative method)
/*
void RCArticulatedBody::determine_implicit_constraint_movement_jacobian(MatrixNd& D)
{
  // determine the number of implicit DOFs 
  unsigned NDOF = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NDOF += _ijoints[i]->num_dof();
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // resize J 
  D.set_zero(NDOF, NGC);
  MatrixNdN q(NGC);

  // init gj and dC
  VectorNd gj = VectorNd::zero(NGC);
  MatrixNd qx(NDOF, NGC);
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
    
    // compute the joint velocity
    for (unsigned r=0; r< _ijoints.size(); r++)
    {
      // get the coordinate index
      unsigned cidx = _ijoints[r]->get_coord_index();

      // get the spatial axis
      const vector<SVelocityd>& si = _ijoints[r]->get_spatial_axes(eGlobal);
      RigidBodyPtr inboard = _ijoints[r]->get_inboard_link();
      RigidBodyPtr outboard = _ijoints[r]->get_outboard_link();
      SVelocityd vh = inboard->get_spatial_velocity(eGlobal);
      SVelocityd vi = outboard->get_spatial_velocity(eGlobal);
// vi = vh + si*qdot, except these are all loops
      SVelocityd dv = vh - vi;
      MatrixNd si_inv = LinAlg::pseudo_inverse(si);
      VectorNd qdot = si_inv.mult(dv);
      for (unsigned j=0; j< _ijoints[r]->num_dof(); j++)
        qx(cidx + j, i) = qdot[j];
    }
  }

  // compute D 
  MatrixNdN qinv;
  LinAlg::pseudo_inverse(q, qinv);
  qx.mult(qinv, D);

  // remove all impulses
  gj.set_one().negate();
  apply_generalized_impulse(DynamicBody::eSpatial, gj);
}
*/

/// Determines the ndof x ngc Jacobian for implicit constraint movement (ndof is the number of degrees of freedom of the implicit constraints)
// TODO: fix this
void RCArticulatedBody::determine_implicit_constraint_jacobians(const EventProblemData& q, MatrixNd& Jx, MatrixNd& Dx) const
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
VectorNd& RCArticulatedBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gc)
{
  SAFESTATIC VectorNd base_gc;

  // resize gc
  gc.resize(num_generalized_coordinates(gctype));

  // get the joint positions of all explicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gc.set_sub_vec(idx, _ejoints[i]->q);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return gc;

  // get the generalized coordinates for the base -- NOTE: we put the base
  // acceleration in the base frame
  assert(!_links.empty());
  RigidBodyPtr base = _links.front();
  base->get_generalized_coordinates_single(gctype, base_gc);
  gc.set_sub_vec(num_joint_dof_explicit(), base_gc);
    
  return gc;
}

/// Sets the generalized position of this body
void RCArticulatedBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gc)
{
  SAFESTATIC VectorNd base_gc;
  assert(num_generalized_coordinates(gctype) == gc.size());

  // set the generalized coordinates for the implicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gc.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->q);
  }

  // update base gc, if necessary
  if (_floating_base)
  {
    assert(!_links.empty());
    RigidBodyPtr base = _links.front();
    gc.get_sub_vec(num_joint_dof_explicit(), gc.size(), base_gc);
    base->set_generalized_coordinates_single(gctype, base_gc);
  }

  // link transforms must now be updated
  update_link_poses();
}

/// Sets the generalized velocity of this body
void RCArticulatedBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gv)
{
  SAFESTATIC VectorNd Dx_qd, base_gv;
  assert(num_generalized_coordinates(gctype) == gv.size());

  // set the generalized velocities for the explicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gv.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->qd);
  }

  // see whether the body has a floating base  
  if (_floating_base)
  {
    assert(!_links.empty());
    RigidBodyPtr base = _links.front();
    gv.get_sub_vec(num_joint_dof_explicit(), gv.size(), base_gv);
    base->set_generalized_velocity(gctype, base_gv);
  }

  // compute implicit constraint velocities
  if (!_ijoints.empty()) 
  {
    determine_implicit_constraint_movement_jacobian(_Dx);
    _Dx.mult(gv, Dx_qd);
    for (unsigned i=0; i< _ijoints.size(); i++)
    {
      unsigned idx = _ijoints[i]->get_coord_index();
      Dx_qd.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _ijoints[i]->qd);
    }
  }

  // link velocities must now be updated
  update_link_velocities();
}

/// Gets the generalized velocity of this body
VectorNd& RCArticulatedBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gv)
{
  SAFESTATIC VectorNd base_gv;

  // resize gv
  gv.resize(num_generalized_coordinates(gctype));

  // get the joint velocities of all joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gv.set_sub_vec(idx, _ejoints[i]->qd);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return gv;

  // get the generalized velocity for the base
  assert(!_links.empty());
  RigidBodyPtr base = _links.front();
  base->get_generalized_velocity_single(gctype, base_gv);
  gv.set_sub_vec(num_joint_dof_explicit(), base_gv);

  return gv;
}

/// Gets the generalized inertia of this body
MatrixNd& RCArticulatedBody::get_generalized_inertia(MatrixNd& M) 
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

/*
/// Determines contact Jacobians
void RCArticulatedBody::determine_contact_jacobians(const EventProblemData& q, 
const VectorNd& v, const MatrixNd& M, MatrixNd& Jc, MatrixNd& Dc)
{
  SAFESTATIC MatrixNd workM, workM2, R;

  // resize Jc and Dc
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  Jc.set_zero(q.N_CONTACTS, NGC);
  Dc.set_zero(q.N_CONTACTS*2, NGC);

  // determine the Jacobian over all contact events 
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+= 2)
  {
    // don't process if contact event is inactive
    if (!q.contact_working_set[i])
      continue;

    // get the articulated bodies of the contacts
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
    
    // see whether we can skip contact event
    if (ab1.get() != this && ab2.get() != this)
      continue;

    // get the rigid bodies corresponding to sb1 and sb2
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(sb1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(sb2);

    // get the contact point and orthonormal basis 
    const Vector3d& p = q.contact_events[i]->contact_point;
    const Vector3d& normal = q.contact_events[i]->contact_normal;
    const Vector3d& tan1 = q.contact_events[i]->contact_tan1;
    const Vector3d& tan2 = q.contact_events[i]->contact_tan2;

    // compute R
    R.resize(3,3);
    R.set_row(0, normal);
    R.set_row(1, tan1);
    R.set_row(2, tan2); 

    // compute the Jacobian at the contact point
    if (ab1 == get_this())
      calc_jacobian(p, rb1, workM);
    else
    {
      calc_jacobian(p, rb2, workM);
      R.negate();
    }

    // transform the Jacobian
    workM.get_sub_mat(0,3,0,workM.columns(), workM2);
    R.mult(workM2, workM); 

    // set the appropriate blocks in the Jacobians
    workM.get_sub_mat(0,1,0,workM.columns(), workM2);
    Jc.set_sub_mat(i, 0, workM2);
    workM.get_sub_mat(1,3,0,workM.columns(), workM2);
    Dc.set_sub_mat(ii, 0, workM2);
  }

  // setup iM_JcT and iM_DcT
  solve_generalized_inertia_transpose(DynamicBody::eSpatial, Jc, _iM_JcT);
  solve_generalized_inertia_transpose(DynamicBody::eSpatial, Dc, _iM_DcT);

  FILE_LOG(LOG_EVENT) << "Jc:" << std::endl << Jc;
  FILE_LOG(LOG_EVENT) << "Dc:" << std::endl << Dc;
}

/// Updates the event data
void RCArticulatedBody::update_event_data(EventProblemData& q)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC MatrixNd tmpM;
  SAFESTATIC MatrixNd M;
  SAFESTATIC VectorNd workv, v;

  // get the generalized velocity (axis angle)
  get_generalized_velocity(DynamicBody::eSpatial, v);

  // get the generalized inertia matrix
  _crb.calc_generalized_inertia(DynamicBody::eSpatial, M);

  // determine contact normal and tangent Jacobians
  determine_contact_jacobians(q, v, M, _Jc, _Dc);

  // setup Jx (neqx x ngc) and Dx (nedof x ngc)
  determine_implicit_constraint_jacobians(q, _Jx, _Dx);

  // setup Jl (nl x ngc)
  const unsigned NGC = v.size();
  _Jl.set_zero(q.N_LIMITS, NGC);
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    // get the limit joint
    JointPtr joint = q.limit_events[i]->limit_joint;

    // see whether to skip this joint
    if (joint->get_articulated_body() != get_this())
      continue;

    // determine the degree-of-freedom
    unsigned dof = joint->get_coord_index() + q.limit_events[i]->limit_dof;

    // set the component of Jl appropriately
    _Jl(i, dof) = (q.limit_events[i]->limit_upper) ? (double) -1.0 : (double) 1.0;
  }

  // setup Dt (nt x ngc)
  _Dt.set_zero(q.N_CONSTRAINT_DOF_IMP, NGC);
  if (use_advanced_friction_model)
  {
    for (unsigned i=0, ii=0; i< q.constraint_events.size(); i++)
    {
      // get the joint
      JointPtr joint = q.constraint_events[i]->constraint_joint;

      // look for a fast skip
      if (joint->get_constraint_type() != Joint::eExplicit)
        continue;

      // determine whether the joint belongs to this body
      if (joint->get_articulated_body() != get_this())
      {
        ii += joint->num_dof();
        continue;
      }
 
      // get the starting index for this joint
      const unsigned ST_IDX = joint->get_coord_index();

      // setup entries of Dt
      for (unsigned j=0; j< joint->num_dof(); j++, ii++)
        _Dt(ii, ST_IDX+j) = (double) 1.0;
    }
  }

  // setup some temporary matrices
  solve_generalized_inertia_transpose(eSpatial, _Jl, _iM_JlT);
  solve_generalized_inertia_transpose(eSpatial, _Jx, _iM_JxT);
  solve_generalized_inertia_transpose(eSpatial, _Dx, _iM_DxT);
  solve_generalized_inertia_transpose(eSpatial, _Dt, _iM_DtT);

  // update all matrices
  q.Jc_iM_JcT += _Jc.mult(_iM_JcT, tmpM);
  q.Jc_iM_DcT += _Jc.mult(_iM_DcT, tmpM);
  q.Jc_iM_JlT += _Jc.mult(_iM_JlT, tmpM);
  q.Jc_iM_DtT += _Jc.mult(_iM_DtT, tmpM);
  q.Jc_iM_JxT += _Jc.mult(_iM_JxT, tmpM);
  q.Jc_iM_DxT += _Jc.mult(_iM_DxT, tmpM);
  q.Dc_iM_DcT += _Dc.mult(_iM_DcT, tmpM);
  q.Dc_iM_JlT += _Dc.mult(_iM_JlT, tmpM);
  q.Dc_iM_DtT += _Dc.mult(_iM_DtT, tmpM);
  q.Dc_iM_JxT += _Dc.mult(_iM_JxT, tmpM);
  q.Dc_iM_DxT += _Dc.mult(_iM_DxT, tmpM);
  q.Jl_iM_JlT += _Jl.mult(_iM_JlT, tmpM);
  q.Jl_iM_DtT += _Jl.mult(_iM_DtT, tmpM);
  q.Jl_iM_JxT += _Jl.mult(_iM_JxT, tmpM);
  q.Jl_iM_DxT += _Jl.mult(_iM_DxT, tmpM);
  q.Dt_iM_DtT += _Dt.mult(_iM_DtT, tmpM);
  q.Dt_iM_JxT += _Dt.mult(_iM_JxT, tmpM);
  q.Dt_iM_DxT += _Dt.mult(_iM_DxT, tmpM);
  q.Jx_iM_JxT += _Jx.mult(_iM_JxT, tmpM);
  q.Jx_iM_DxT += _Jx.mult(_iM_DxT, tmpM);
  q.Dx_iM_DxT += _Dx.mult(_iM_DxT, tmpM);

  // update velocity vectors
  q.Jc_v += _Jc.mult(v, workv);
  q.Dc_v += _Dc.mult(v, workv);
  q.Jl_v += _Jl.mult(v, workv);
  q.Jx_v += _Jx.mult(v, workv);
  q.Dx_v += _Dx.mult(v, workv);
}

/// Updates the body velocities
void RCArticulatedBody::update_velocity(const EventProblemData& q)
{
  SAFESTATIC VectorNd workv, v;
  SAFESTATIC MatrixNd workM;
  SAFESTATIC vector<unsigned> normal_contact_indices, tangent_contact_indices;

  // get the current spatial velocity
  get_generalized_velocity(DynamicBody::eSpatial, v);

  if (LOGGING(LOG_EVENT))
  {
    MatrixNd M;
    FILE_LOG(LOG_EVENT) << "alpha_c: " << q.alpha_c << std::endl;
    FILE_LOG(LOG_EVENT) << "beta_c: " << q.beta_c << std::endl;
    FILE_LOG(LOG_EVENT) << "alpha_l: " << q.alpha_l << std::endl;
    FILE_LOG(LOG_EVENT) << "beta_t: " << q.beta_t << std::endl;
    FILE_LOG(LOG_EVENT) << "velocity before: " << v << std::endl;
    FILE_LOG(LOG_EVENT) << "kinetic energy before: " << calc_kinetic_energy() << std::endl;
    FILE_LOG(LOG_EVENT) << "kinetic energy before (calculation 2): " << v.dot(get_generalized_inertia(DynamicBody::eSpatial, M).mult(v*0.5)) << std::endl;
  }

  // check whether we are using a subset of the contacts
  if (std::find(q.contact_working_set.begin(), q.contact_working_set.end(), false) != q.contact_working_set.end())
  {
    // setup normal and tangent contact indices
    normal_contact_indices.clear();
    tangent_contact_indices.clear();
    for (unsigned i=0, j=0; i< q.contact_working_set.size(); i++, j+=2)
      if (q.contact_working_set[i])
      {
        normal_contact_indices.push_back(i);
        tangent_contact_indices.push_back(j);
        tangent_contact_indices.push_back(j+1);
      }

    // compute change in velocities
    _iM_JcT.select_columns(normal_contact_indices.begin(), normal_contact_indices.end(), workM);
    v += workM.mult(q.alpha_c, workv);
    _iM_DcT.select_columns(tangent_contact_indices.begin(), tangent_contact_indices.end(), workM);
    v += workM.mult(q.beta_c, workv);
  }
  else
  {
    // compute change in velocities
    v += _iM_JcT.mult(q.alpha_c, workv);
    v += _iM_DcT.mult(q.beta_c, workv);
  }
  v += _iM_JlT.mult(q.alpha_l, workv);
  v += _iM_JxT.mult(q.alpha_x, workv);
  v += _iM_DxT.mult(q.beta_x, workv);
  v += _iM_DtT.mult(q.beta_t, workv);

  // set the spatial velocity
  set_generalized_velocity(DynamicBody::eSpatial, v);

  FILE_LOG(LOG_EVENT) << "velocity after: " << v << std::endl;
  FILE_LOG(LOG_EVENT) << "kinetic energy after: " << calc_kinetic_energy() << std::endl;
}
*/

/// Gets the generalized forces on this body
/**
 * \note does not add forces from implicit joint constraints!
 */
VectorNd& RCArticulatedBody::get_generalized_forces(VectorNd& f) 
{
  const unsigned SPATIAL_DIM = 6, X = 0, Y = 1, Z = 2, A = 3, B = 4, C = 5;

  // resize f
  f.resize(num_generalized_coordinates(DynamicBody::eSpatial));

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
VectorNd& RCArticulatedBody::convert_to_generalized_force(SingleBodyPtr body, const SForced& w, const Point3d& p, VectorNd& gf)
{
  const unsigned SPATIAL_DIM = 6;
  static vector<SVelocityd> J;
  static vector<SVelocityd> sprime;

  // get the body as a rigid body
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  assert(link);

  // get the gc frame
  shared_ptr<const Pose3d> P = _links.front()->get_gc_pose();

  // get w in P's frame
  SForced wP = Pose3d::transform(P, w);

  // compute the Jacobian in w's frame
  J.resize(num_joint_dof_explicit(), SVelocityd::zero(P));  
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // if link is not a descent of this joint's inboard, keep looping
    RigidBodyPtr inboard = _ejoints[i]->get_inboard_link();
    if (!inboard->is_descendant_link(link))
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
    else if (strcasecmp(frame.c_str(), "linkinertia") == 0)
      set_computation_frame_type(eLinkInertia);
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

  if (xlat_attr)
  {
    RigidBodyPtr base = get_base_link();
    base->translate(xlat_attr->get_origin_value());
    update_link_poses();
  }
  else if (rotate_attr)
  {
    RigidBodyPtr base = get_base_link();
    base->rotate(rotate_attr->get_rpy_value());
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

