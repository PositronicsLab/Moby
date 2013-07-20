/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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

  // setup inertial pose
  _jF = shared_ptr<Pose3d>(new Pose3d(Pose3d::identity()));
  _jF->rpose = _F;

  // setup visualization pose
  _vF->rpose = _F;

  // set everything else
  _J.set_zero();
  _xd.set_zero();
  _xdd.set_zero();
  _force.set_zero();
  _enabled = true;
  _link_idx = std::numeric_limits<unsigned>::max();
  viscous_coeff = VectorNd::zero(SPATIAL_DIM);
}

/// Gets the frame in which kinematics and dynamics computations occur
shared_ptr<const Pose3d> RigidBody::get_computation_frame() const
{
  switch (_rftype)
  {
    case eLink:
      return _F;

    case eGlobal:
      return shared_ptr<const Pose3d>();
    
    case eJoint:
      return (_abody.expired() || is_base()) ? _F : get_inner_joint_implicit()->get_pose();

    default:
      assert(false);
  }

  return shared_ptr<const Pose3d>();
}

/// Rotates the rigid body
void RigidBody::rotate(const Quatd& q)
{
  // update the rotation 
  _F->q *= q;
}

/// Computes the Jacobian
vector<SVelocityd>& RigidBody::calc_jacobian(shared_ptr<const Pose3d> frame, DynamicBodyPtr body, vector<SVelocityd>& J)
{
  const unsigned SPATIAL_DIM = 6;

  if (body != shared_from_this())
    throw std::runtime_error("RigidBody::calc_jacobian() called with wrong body!");

  // J will be a 6-dimensional vector
  J.resize(SPATIAL_DIM);

  // setup Jacobian
  J[0] = SVelocityd(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, _F);
  J[1] = SVelocityd(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, _F);
  J[2] = SVelocityd(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, _F);
  J[3] = SVelocityd(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, _F);
  J[4] = SVelocityd(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, _F);
  J[5] = SVelocityd(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, _F);

  // transform J to given pose 
  for (unsigned i=0; i< SPATIAL_DIM; i++) 
    J[i] = Pose3d::transform(_F, frame, J[i]);

  return J;
}

/// Translates the rigid body
void RigidBody::translate(const Origin3d& x)
{
  // update the translation
  _F->x += x;
}

/// (Re)sets the computation frame
void RigidBody::set_computation_frame_type(ReferenceFrameType rftype)
{
  // correct rftype if necessary
  if (_abody.expired() && rftype == eJoint)
    rftype = eLink;

  // check whether we need to do anything at all
  if (rftype == _rftype)
    return;

  // get the target reference frame and the necessary transform
  shared_ptr<const Pose3d> target;
  switch (rftype)
  {
    case eLink:   
      target = _F;     
      break;
 
   case eGlobal: 
      target = GLOBAL; 
      break;

    case eJoint:
      target = (is_base()) ? _F : get_inner_joint_implicit()->get_pose(); 
      break;

    default: 
      assert(false);
  }

  // otherwise, transform everything
  Transform3d T = Pose3d::calc_relative_pose(_J.pose, target);

  // convert all quantities
  _xd = T.transform(_xd);
  _xdd = T.transform(_xdd);
  _J = T.transform(_J);
  _force = T.transform(_force);

  // store the new reference frame type
  _rftype = rftype;
}

/// Integrates the body forward in time
void RigidBody::integrate(double t, double h, shared_ptr<Integrator> integrator)
{
  // don't attempt to integrate disabled bodies
  if (!is_enabled())
    return;

  // if this body is a link in an articulated body, don't attempt to integrate;
  // integration must occur at the articulated body
  if (!_abody.expired())
    return;

  // call parent method if still here
  DynamicBody::integrate(t, h, integrator);

  FILE_LOG(LOG_DYNAMICS) << "RigidBody::integrate()" << endl;
  FILE_LOG(LOG_DYNAMICS) << "  new transform: " << endl << _F;
  FILE_LOG(LOG_DYNAMICS) << "  new velocity: " << _xd << endl;
}

/// Computes the acceleration due to inertial forces on the body
SForced RigidBody::calc_inertial_forces() const
{
  return _xd.cross(get_inertia() * _xd); 
}

/// Computes the forward dynamics for this body
void RigidBody::calc_fwd_dyn(double dt)
{
  // if the body is free, just compute linear and angular acceleration via
  // Newton's and Euler's laws
  if (_abody.expired())
    _xdd = _J.inverse_mult(sum_forces() - calc_inertial_forces()); 
  else
  {
    // otherwise, need to call forward dynamics on the articulated body
    ArticulatedBodyPtr abody(_abody);
  
    // calculate forward dynamics on it
    abody->calc_fwd_dyn(dt);
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

  // if disabled, then zero the velocities
  if (!_enabled)
  {
    _xd.set_zero();
    _xdd.set_zero();
  }
}

/// Sets the linear velocity of this body
SVelocityd& RigidBody::velocity() 
{ 
  // set the velocity  
  return _xd; 
}

/// Sets the inertia tensor for this body
void RigidBody::set_inertia(const SpatialRBInertiad& inertia)
{
  // store the matrix and its inverse
  _J = inertia;
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

  *_F = P;
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

  // only need to update quantities if they're not in link frame
  // (or equivalently, set to joint frame and there is no inner joint)
  if (_rftype == eLink || (_rftype == eJoint && is_base()))
  {
    // just update the pose
    *_F = p;
  }
  else
  {
    // update pose now
    *_F = p;

    // get current frame for quantities
    // this will either be the joint frame or the global frame
    shared_ptr<const Pose3d> SOURCE = _J.pose;
 
    // compute transform for moving quantities to the new body frame
    Transform3d pTs = Pose3d::calc_relative_pose(SOURCE->rpose, _F);

    // compute target transform
    shared_ptr<const Pose3d> target = (_rftype == eJoint) ? get_inner_joint_implicit()->get_pose() : GLOBAL;
    Transform3d tTp = Pose3d::calc_relative_pose(_F, target);

    // combine the two transforms
    Transform3d T = tTp * pTs;

    // convert all quantities
    _xd = T.transform(_xd);
    _xdd = T.transform(_xdd);
    _J = T.transform(_J);
    _force = T.transform(_force);
  }
}

/// Gets the desired child link
RigidBodyPtr RigidBody::get_child_link(JointPtr j) const
{
  assert(_outer_joints.find(j) != _outer_joints.end());
  return j->get_outboard_link();
}

/// Adds a force to the body
void RigidBody::add_force(const SForced& w)
{
  // do not add forces to disabled bodies
  if (!_enabled)
    return;
  
  // update the force 
  _force += w;
}

/// Calculates the velocity of a point on this rigid body in the body frame
Vector3d RigidBody::calc_point_vel(const Point3d& point) const
{
  // if the body is disabled, point velocity is zero
  if (!_enabled)
    return Vector3d::zero(_F);

  // convert point to a vector in the body frame
  Vector3d r = Pose3d::transform(point.pose, _F, point); 

  // get the velocity in the body frame
  SVelocityd xd = Pose3d::transform(_xd.pose, _F, _xd);

  // compute the point velocity - in the body frame
  return xd.get_linear() + Vector3d::cross(xd.get_angular(), r);
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
 
  // read the viscous dampening coefficient, if provided
  const XMLAttrib* viscous_coeff_attr = node->get_attrib("viscous-dampening-coeff");
  if (viscous_coeff_attr)
    viscous_coeff_attr->get_vector_value(viscous_coeff);
 
  // read whether the body is enabled, if provided
  const XMLAttrib* enabled_attr = node->get_attrib("enabled");
  if (enabled_attr)
    _enabled = enabled_attr->get_bool_value();

  // read the mass, if provided
  const XMLAttrib* mass_attr = node->get_attrib("mass");
  if (mass_attr)
  {
    SpatialRBInertiad J = get_inertia();
    J.m = mass_attr->get_real_value();
    set_inertia(J);
  }

  // read the inertia matrix, if provided
  const XMLAttrib* inertia_attr = node->get_attrib("inertia");
  if (inertia_attr)
  {
    SpatialRBInertiad J = get_inertia();
    inertia_attr->get_matrix_value(J.J);
    set_inertia(J);
  }

  // read the position and orientation, if provided
  const XMLAttrib* position_attr = node->get_attrib("position");
  const XMLAttrib* rpy_attr = node->get_attrib("rpy");
  const XMLAttrib* quat_attr = node->get_attrib("quat");
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
  const XMLAttrib* com_attr = node->get_attrib("inertial-relative-com");
  const XMLAttrib* J_rpy_attr = node->get_attrib("inertial-relative-rpy");
  const XMLAttrib* J_quat_attr = node->get_attrib("inertial-relative-quat");
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
    SpatialRBInertiad J;

    // loop over all InertiaFromPrimitive nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = ifp_nodes.begin(); i != ifp_nodes.end(); i++)
    {
      // make sure the child node has the ID
      const XMLAttrib* pid_attr = (*i)->get_attrib("primitive-id");
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

      // get the inertia and inertial frame from the primitive
      const SpatialRBInertiad& Jx = primitive->get_inertia();
      shared_ptr<const Pose3d> Fx = primitive->get_inertial_pose();

      // set the relative transformation initially to identity for this primitive
      shared_ptr<Pose3d> rTR(new Pose3d); 

      // read the relative transformation, if specified
      const XMLAttrib* rel_origin_attr = (*i)->get_attrib("relative-origin");
      const XMLAttrib* rel_rpy_attr = (*i)->get_attrib("relative-rpy");
      if (rel_origin_attr)
        rTR->x = rel_origin_attr->get_origin_value();
      if (rel_rpy_attr)
        rTR->q = rel_rpy_attr->get_rpy_value();
      rTR->rpose = Fx;

      // transform the inertia and update the inertia for this
      J += Pose3d::transform(Fx, rTR, Jx); 
    }

    // set the mass and inertia of the RigidBody additively
    set_inertia(J);
  }

  // read the linear and/or velocity of the body, if provided
  const XMLAttrib* lvel_attr = node->get_attrib("linear-velocity");
  const XMLAttrib* avel_attr = node->get_attrib("angular-velocity");
  if (lvel_attr || avel_attr)
  {
    Vector3d lv = Vector3d::zero(), av = Vector3d::zero();
    shared_ptr<const Pose3d> TARGET(new Pose3d(Quatd::identity(), _F->x)); 
    SVelocityd v(TARGET);
    if (lvel_attr) lvel_attr->get_vector_value(lv);
    if (avel_attr) avel_attr->get_vector_value(av);
    v.set_linear(lv);
    v.set_angular(av);
    velocity() = Pose3d::transform(TARGET, get_computation_frame(), v);
  }

/*
  // read in the vector from the inner joint to the com in link coordinates
  const XMLAttrib* d_attr = node->get_attrib("inner-joint-to-com-vector-link");
  if (d_attr)
  {
    Vector3 d;
    d_attr->get_vector_value(d);
    set_inner_joint_to_com_vector_link(d);
  }
*/
  // read the articulated body, if given
  const XMLAttrib* ab_attr = node->get_attrib("articulated-body-id");
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

  // save the mass
  node->attribs.insert(XMLAttrib("mass", _J.m));

  // save the inertia
  node->attribs.insert(XMLAttrib("inertia", _J.J));

  // convert the current pose to be with respect to global coordinates
  Pose3d F0 = *_F;
  F0.update_relative_pose(GLOBAL);
  node->attribs.insert(XMLAttrib("position", F0.x));
  node->attribs.insert(XMLAttrib("quat", F0.q));

  // save the inertial frame
  node->attribs.insert(XMLAttrib("inertial-relative-com", _jF->x));
  node->attribs.insert(XMLAttrib("inertial-relative-quat", _jF->q));

  // save the linear and angular velocities
  shared_ptr<const Pose3d> TARGET(new Pose3d(Quatd::identity(), _F->x)); 
  SVelocityd v = Pose3d::transform(_xd.pose, TARGET, _xd); 
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

    // update linear and angular velocities 
    _xd += _J.inverse_mult(w);

    // reset the force and torque accumulators for this body
    _force.set_zero();
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
 
/// Adds a generalized force to this rigid body
void RigidBody::add_generalized_force(const VectorNd& gf)
{
  assert(gf.size() == num_generalized_coordinates(DynamicBody::eSpatial));
  SForced w;

  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  // set the pose for w
  w.pose = GLOBAL;

  // get the force and torque
  w.set_force(Vector3d(gf[0], gf[1], gf[2]));
  w.set_torque(Vector3d(gf[3], gf[4], gf[5]));

  // add the force to the sum of forces
  _force += Pose3d::transform(GLOBAL, _force.pose, w);
}

/// Applies a generalized impulse to this rigid body
void RigidBody::apply_generalized_impulse(const VectorNd& gj)
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
void RigidBody::apply_generalized_impulse_single(const VectorNd& gj)
{
  SMomentumd w;

  // don't do anything if this body is disabled
  if (!_enabled)
    return;

  // simple error check...
  assert(gj.size() == num_generalized_coordinates(DynamicBody::eSpatial));

  // clear the force and torque accumulators
  _force.set_zero();

  // get the impulses
  w.set_linear(Vector3d(gj[0], gj[1], gj[2]));
  w.set_angular(Vector3d(gj[3], gj[4], gj[5]));
  w.pose = GLOBAL;

  // determine the change in linear velocity
  _xd += _J.inverse_mult(Pose3d::transform(GLOBAL, _J.pose, w));
}

/// Solves using the generalized inertia matrix
MatrixNd& RigidBody::transpose_solve_generalized_inertia(const MatrixNd& B, MatrixNd& X)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->transpose_solve_generalized_inertia(B, X);
  }
  else
    return transpose_solve_generalized_inertia_single(B, X);
}

/// Solves using the generalized inertia matrix
MatrixNd& RigidBody::solve_generalized_inertia(const MatrixNd& B, MatrixNd& X)
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
MatrixNd& RigidBody::transpose_solve_generalized_inertia_single(const MatrixNd& B, MatrixNd& X)
{
  // get proper generalized inertia matrix
  MatrixNd M;
  get_generalized_inertia(M);
  MatrixNd::transpose(B, X);
  _LA.solve_fast(M, X);

  return X;
}

/// Solves using the generalized inertia matrix (does not call articulated body version)
MatrixNd& RigidBody::solve_generalized_inertia_single(const MatrixNd& B, MatrixNd& X)
{
  // get proper generalized inertia matrix
  MatrixNd M;
  get_generalized_inertia(M);
  X = B;
  _LA.solve_fast(M, X);

  return X;
}

/// Solves using the generalized inertia matrix
VectorNd& RigidBody::solve_generalized_inertia(const VectorNd& b, VectorNd& x)
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
VectorNd& RigidBody::solve_generalized_inertia_single(const VectorNd& b, VectorNd& x)
{
  // get proper generalized inertia matrix
  MatrixNd M;
  get_generalized_inertia(M);
  x = b;
  _LA.solve_fast(M, x);

  return x;
}

/// Gets the generalized position of this rigid body (does not call articulated body version)
VectorNd& RigidBody::get_generalized_coordinates_single(GeneralizedCoordinateType gctype, VectorNd& gc) 
{
  // special case: disabled body
  if (!_enabled)
    return gc.resize(0);

  // resize vector
  gc.resize(num_generalized_coordinates(gctype));

  // convert current pose to global frame
  Pose3d P = *_F;
  P.update_relative_pose(GLOBAL);

  // get linear components
  gc[0] = P.x[0];
  gc[1] = P.x[1];
  gc[2] = P.x[2];

  // get angular components 
  if (gctype == DynamicBody::eSpatial)
    P.q.to_rpy(gc[3], gc[4], gc[5]);
  else
  {
    // return the generalized position using Euler parameters
    assert(gctype == DynamicBody::eEuler);
    gc[3] = P.q.w;
    gc[4] = P.q.x;
    gc[5] = P.q.y;
    gc[6] = P.q.z;
  }

  return gc; 
}

/// Gets the generalized position of this rigid body
VectorNd& RigidBody::get_generalized_coordinates(GeneralizedCoordinateType gctype, VectorNd& gc) 
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->get_generalized_coordinates(gctype, gc);
  }
  else
    return get_generalized_coordinates_single(gctype, gc);
}

/// Sets the generalized coordinates of this rigid body (does not call articulated body)
void RigidBody::set_generalized_coordinates_single(GeneralizedCoordinateType gctype, const VectorNd& gc)
{
  assert(gc.size() == num_generalized_coordinates(gctype));

  // special case: disabled body
  if (!_enabled)
    return;

  // do easiest case first 
  if (gctype == DynamicBody::eSpatial)
  {
    // note: generalized coordinates in eSpatial are always set with regard to 
    // the global frame
    Origin3d x(gc[0], gc[1], gc[2]);
    Quatd q = Quatd::rpy(gc[3], gc[4], gc[5]); 

    // convert the pose to the correct relative frame
    Pose3d P(q, x);
    P.update_relative_pose(_F->rpose);

    // set the transform
    set_pose(P);
  }
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // get the position
    Origin3d x(gc[0], gc[1], gc[2]);

    // get the unit quaternion
    Quatd q;
    q.w = gc[3];
    q.x = gc[4];
    q.y = gc[5];
    q.z = gc[6];

    // normalize the unit quaternion, just in case
    q.normalize();

    // coordinates are in the global frame; must convert them to the
    // relative frame
    Pose3d P(q, x);
    P.update_relative_pose(_F->rpose); 

    // set the transform
    set_pose(P);
  }
}

/// Sets the generalized coordinates of this rigid body
void RigidBody::set_generalized_coordinates(GeneralizedCoordinateType gctype, const VectorNd& gc)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->set_generalized_coordinates(gctype, gc);
  }
  else
    set_generalized_coordinates_single(gctype, gc);
}

/// Sets the generalized velocity of this rigid body
void RigidBody::set_generalized_velocity(GeneralizedCoordinateType gctype, const VectorNd& gv)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    ab->set_generalized_velocity(gctype, gv);
  }
  else
    set_generalized_velocity_single(gctype, gv);
}

/// Sets the generalized velocity of this rigid body (does not call articulated body version)
void RigidBody::set_generalized_velocity_single(GeneralizedCoordinateType gctype, const VectorNd& gv)
{
  SVelocityd xd;

  // special case: disabled body
  if (!_enabled)
    return;
  
  // simplest case: spatial coordinates
  if (gctype == DynamicBody::eSpatial)
  {
    xd.set_linear(Vector3d(gv[0], gv[1], gv[2]));
    xd.set_angular(Vector3d(gv[3], gv[4], gv[5]));
    _xd = Pose3d::transform(GLOBAL, _xd.pose, xd);
  }
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // get the quaternion derivatives
    Quatd qd;
    qd.w = gv[3] * 2.0;
    qd.x = gv[4] * 2.0;
    qd.y = gv[5] * 2.0;
    qd.z = gv[6] * 2.0;

    // go ahead and set the linear contribution to xd
    xd.set_linear(Vector3d(gv[0], gv[1], gv[2]));

    // computation frame must be taken into account
    switch (get_computation_frame_type())
    {
      case eGlobal:
        xd.pose = GLOBAL;
        xd.set_angular(_F->q.G_mult(qd.w, qd.x, qd.y, qd.z));
        _xd = Pose3d::transform(GLOBAL, _xd.pose, xd);
        break;

      case eLink:
        xd.pose = _F;
        xd.set_angular(_F->q.L_mult(qd.w, qd.x, qd.y, qd.z));
        _xd = Pose3d::transform(_F, _xd.pose, xd);
        break;

      case eJoint:
        if (!(_abody.expired() || is_base()))
          throw std::runtime_error("Euler coordinates computed in joint frame not supported");
        break;
    }
  }
}

/// Gets the generalized velocity of this rigid body
VectorNd& RigidBody::get_generalized_velocity(GeneralizedCoordinateType gctype, VectorNd& gv) 
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->get_generalized_velocity(gctype, gv);
  }
  else
    return get_generalized_velocity_single(gctype, gv);
}

/// Gets the generalized velocity of this rigid body (does not call articulated body version)
VectorNd& RigidBody::get_generalized_velocity_single(GeneralizedCoordinateType gctype, VectorNd& gv) 
{
  SVelocityd xd;

  // special case: disabled body
  if (!_enabled)
    return gv.resize(0);

  // resize the generalized velocity
  gv.resize(num_generalized_coordinates(gctype));

  // determine the proper generalized coordinate type
  if (gctype == DynamicBody::eSpatial)
  {
    // transform velocity to global frame 
    xd = Pose3d::transform(_xd.pose, GLOBAL, _xd);

    // get linear and angular components of velocity
    Vector3d lv = xd.get_linear();
    Vector3d av = xd.get_angular();

    // setup the components
    gv[0] = lv[0];
    gv[1] = lv[1];
    gv[2] = lv[2];
    gv[3] = av[0];
    gv[4] = av[1];
    gv[5] = av[2];
  }
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // going to need Euler coordinate derivatives
    Quatd qd;

    // computation frame must be taken into account
    switch (get_computation_frame_type())
    {
      case eGlobal:
        xd = Pose3d::transform(_xd.pose, GLOBAL, _xd);
        qd = _F->q.G_transpose_mult(xd.get_angular()) * 0.5;
        break;

      case eLink:
        xd = Pose3d::transform(_xd.pose, _F, _xd);
        qd = _F->q.L_transpose_mult(xd.get_angular()) * 0.5;
        break;

      case eJoint:
        if (!(_abody.expired() || is_base()))
          throw std::runtime_error("Euler coordinates computed in joint frame not supported");
        break;
    }

    // setup the generalized velocity
    Vector3d lv = xd.get_linear();
    gv[0] = lv[0];
    gv[1] = lv[1];
    gv[2] = lv[2];
    gv[3] = qd.w;
    gv[4] = qd.x;
    gv[5] = qd.y;
    gv[6] = qd.z; 
  }

  return gv;
}

/// Gets the generalized acceleration of this body
VectorNd& RigidBody::get_generalized_acceleration(VectorNd& ga)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->get_generalized_acceleration(ga);
  }
  else
    return get_generalized_acceleration_single(ga);
}

/// Gets the generalized acceleration of this body (does not call articulated body version)
VectorNd& RigidBody::get_generalized_acceleration_single(VectorNd& ga)
{
  // special case: body is disabled
  if (!_enabled)
    return ga.resize(0);

  // setup the linear components
  ga.resize(num_generalized_coordinates(DynamicBody::eSpatial));

  // transform acceleration to global frame 
  SAcceld xdd = Pose3d::transform(_xdd.pose, GLOBAL, _xdd);

  // get linear and angular components
  Vector3d la = xdd.get_linear();
  Vector3d aa = xdd.get_angular();

  // set linear components
  ga[0] = la[0];
  ga[1] = la[1];
  ga[2] = la[2];
  ga[3] = aa[0];
  ga[4] = aa[1];
  ga[5] = aa[2];
 
  return ga;
}

/// Gets the generalized inertia of this rigid body
MatrixNd& RigidBody::get_generalized_inertia(MatrixNd& M)
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
MatrixNd& RigidBody::get_generalized_inertia_single(MatrixNd& M) 
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;

  // special case: disabled body
  if (!_enabled)
    return M.resize(0,0);

  // get the rigid body inertia as a matrix in the global frame
  SpatialRBInertiad J = Pose3d::transform(_J.pose, GLOBAL, _J);

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

/// Gets the generalized inertia of this rigid body
VectorNd& RigidBody::get_generalized_forces(VectorNd& gf)
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
VectorNd& RigidBody::get_generalized_forces_single(VectorNd& gf) 
{
  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // resize the generalized forces vector
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial); 
  gf.resize(NGC);

  // compute external forces and inertial forces in global frame
  SForced w = Pose3d::transform(_force.pose, GLOBAL, _force - calc_inertial_forces());

  // get force and torque
  Vector3d f = w.get_force();
  Vector3d t = w.get_torque();

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
VectorNd& RigidBody::convert_to_generalized_force(SingleBodyPtr body, const SForced& w, const Point3d& p, VectorNd& gf) 
{
  // if this belongs to an articulated body, call the articulated body method
  if (!_abody.expired())
  {
    ArticulatedBodyPtr ab(_abody);
    return ab->convert_to_generalized_force(body, w, p, gf); 
  }
  else
    return convert_to_generalized_force_single(body, w, gf);
}

/// Converts a force to a generalized force (does not call articulated body version)
VectorNd& RigidBody::convert_to_generalized_force_single(SingleBodyPtr body, const SForced& w, VectorNd& gf) 
{
  // verify that body == this
  assert(body.get() == this);

  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // transform w to computation frame
  SForced wt = Pose3d::transform(w.pose, GLOBAL, w);

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
double RigidBody::calc_kinetic_energy() const
{
  if (!_enabled)
    return (double) 0.0;

  // convert J to world frame
  return _xd.dot(_J * _xd) * 0.5;
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

/// Gets the implicit inner joint of this link; returns NULL if there is no implicit inner joint
/**
 * Throws an exception if this link has multiple implicit inner joints
 */
JointPtr RigidBody::get_inner_joint_implicit() const
{
  JointPtr ij;
  BOOST_FOREACH(JointPtr j, _inner_joints)
  {
    if (j->get_constraint_type() == Joint::eImplicit)
    {
      if (ij)
        throw std::runtime_error("Multiple implicit joints detected for a single link!"); 
      else
        ij = j;
    }
  }

  return ij; 
}

// TODO: fix this (if necessary)
/*
/// Updates the velocity of this body using impulses computed via event data
void RigidBody::update_velocity(const EventProblemData& q)
{
  // check for easy exit
  if (q.N_CONTACTS == 0 || !_enabled)
    return;

  // setup summed impulses
  Vector3 j .set_zero()

, k .set_zero()

;

  FILE_LOG(LOG_EVENT) << "RigidBody::update_velocity() entered" << std::endl;

  // update velocity using contact impulses
  for (unsigned i=0, s=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1.get() != this && sb2.get() != this)
    {
      s += 2;
      continue;
    }

    FILE_LOG(LOG_EVENT) << "  contact impulse magnitudes: " << q.alpha_c[i] << " " << q.beta_c[s] << " " << q.beta_c[s+1] << std::endl;

    // check whether to negate normal
    bool negate = (sb2.get() == this);

    // prepare for computation
    Vector3 n = q.contact_events[i]->contact_normal;
    if (negate)
      n = -n;
    const Vector3& p = q.contact_events[i]->contact_point;
    Vector3 r = p - _x;
    Vector3 cross = Vector3::cross(r, n);

    // update j and k
    j += n*q.alpha_c[i];
    k += cross*q.alpha_c[i];

    // now do tangent impulses
    n = q.contact_events[i]->contact_tan1;
    if (negate)
      n = -n;
    cross = Vector3::cross(r, n);
    j += n*q.beta_c[s];
    k += cross*q.beta_c[s];
    s++;
    n = q.contact_events[i]->contact_tan2;
    if (negate)
      n = -n;
    cross = Vector3::cross(r, n);
    j += n*q.beta_c[s];
    k += cross*q.beta_c[s];
    s++;
  }

  // look for whether we update
  if (j.norm() < NEAR_ZERO && k.norm() < NEAR_ZERO)
    return;
  else
    apply_impulse(j, k, _x);

  FILE_LOG(LOG_EVENT) << "  applying impulses: " << j << " / " << k << std::endl;
  FILE_LOG(LOG_EVENT) << "  new linear velocity for " << id << ": " << _xd << std::endl;
  FILE_LOG(LOG_EVENT) << "  new angular velocity for " << id << ": " << _omega << std::endl;
  for (unsigned i=0; i< q.contact_events.size(); i++)
    FILE_LOG(LOG_EVENT) << "  contact velocity at " << q.contact_events[i]->contact_point << " along normal " << q.contact_events[i]->contact_normal << ": " << q.contact_events[i]->contact_normal.dot(calc_point_vel(q.contact_events[i]->contact_point)) << std::endl;
  FILE_LOG(LOG_EVENT) << "RigidBody::update_velocity() exited" << std::endl;
}
*/

// TODO: fix this (if necessary)
/*
/// Adds contributions to the event matrices
void RigidBody::update_event_data(EventProblemData& q) 
{
  if (q.N_CONTACTS == 0 || !_enabled)
    return;

  // get inertia matrix in global frame
  Matrix3 R(&_q);
  Matrix3 invJ = R * _invJ * Matrix3::transpose(R);

  // NOTE: b/c this is an individual rigid body, we don't touch the constraint
  // or limit matrices or Ji

  // 1. update Jc_iM_JcT and Jc_v
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    // verify that it is the proper type
    assert(q.contact_events[i]->event_type == Event::eContact);

    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate normal
    bool negate1 = (sb2 == get_this());

    // prepare for computation
    Vector3 n1 = q.contact_events[i]->contact_normal;
    const Vector3& p1 = q.contact_events[i]->contact_point;
    Vector3 r1 = p1 - _x;
    Vector3 cross1 = Vector3::cross(r1, n1);

    // update Jc_v
    if (!negate1)
      q.Jc_v[i] += n1.dot(_xd + Vector3::cross(_omega, r1));
    else
      q.Jc_v[i] -= n1.dot(_xd + Vector3::cross(_omega, r1));

    // scale n1 and cross1
    n1 *= _inv_mass;
    cross1 = invJ * cross1;

    // loop again, note: the matrices are symmetric
    for (unsigned j=i; j< q.contact_events.size(); j++)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());

      // prepare for computation
      const Vector3& n2 = q.contact_events[j]->contact_normal;
      const Vector3& p2 = q.contact_events[j]->contact_point;
      Vector3 r2 = p2 - _x;
      Vector3 cross2 = Vector3::cross(r2, n2);

      // compute entry of Jc_iM_JcT
      double sum = n1.dot(n2) + cross1.dot(cross2);
      if ((negate1 && !negate2) || (negate2 && !negate1))
        sum = -sum;
      q.Jc_iM_JcT(i, j) += sum; 
      q.Jc_iM_JcT(j, i) = q.Jc_iM_JcT(i, j);
    }
  }

  // 2. update Dc_v
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
    {
      ii += 2;
      continue;
    }

    // check whether to negate normal
    bool negate = (sb2 == get_this());

    // prepare for computation
    const Vector3& p = q.contact_events[i]->contact_point;
    Vector3 r = p - _x;

    // get the contact tangents
    const Vector3& d1 = q.contact_events[i]->contact_tan1;
    const Vector3& d2 = q.contact_events[i]->contact_tan2;

    // determine velocity vector
    Vector3 vec =  _xd + Vector3::cross(_omega, r);

    // update Dc_v
    if (!negate)
    {
      q.Dc_v[ii++] += d1.dot(vec);
      q.Dc_v[ii++] += d2.dot(vec);
    }
    else
    {
      q.Dc_v[ii++] -= d1.dot(vec);
      q.Dc_v[ii++] -= d2.dot(vec);
    }
  } 

  // 3. update Dc_iM_JcT
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate normal
    bool negate1 = (sb2 == get_this());

    // prepare for computation
    Vector3 n1 = q.contact_events[i]->contact_normal;
    const Vector3& p1 = q.contact_events[i]->contact_point;
    Vector3 r1 = p1 - _x;
    Vector3 cross1 = Vector3::cross(r1, n1);

    // scale n1 and cross1
    n1 *= _inv_mass;
    cross1 = invJ * cross1;

    // loop over all contacts 
    for (unsigned j=0, jj=0; j< q.contact_events.size(); j++)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
      {
        jj += 2;
        continue;
      }

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());
      bool negate = ((negate1 && !negate2) || (negate2 && !negate1));

      // prepare for computation
      const Vector3& p2 = q.contact_events[j]->contact_point;
      Vector3 r2 = p2 - _x;

      // compute cross products for both tangent directions
      const Vector3& d21 = q.contact_events[j]->contact_tan1;
      const Vector3& d22 = q.contact_events[j]->contact_tan2;
      Vector3 cross21 = Vector3::cross(r2, d21);
      Vector3 cross22 = Vector3::cross(r2, d22);

      // compute entries of Jc_iM_DcT
      double sum1 = n1.dot(d21) + cross1.dot(cross21);
      double sum2 = n1.dot(d22) + cross1.dot(cross22);
      if (!negate)
      {
        q.Jc_iM_DcT(i,jj++) += sum1;
        q.Jc_iM_DcT(i,jj++) += sum2;
      }
      else
      {
        q.Jc_iM_DcT(i,jj++) -= sum1;
        q.Jc_iM_DcT(i,jj++) -= sum2;
      }
    }
  }
 
  // 4. update Dc_iM_DcT
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+=2)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate
    bool negate1 = (sb2 == get_this());

    // get the moment arm 
    Vector3 r1 = q.contact_events[i]->contact_point - _x;

    // get the contact tangents
    Vector3 d1a = q.contact_events[i]->contact_tan1;
    Vector3 d1b = q.contact_events[i]->contact_tan2;

    // compute the cross products
    Vector3 cross1a = Vector3::cross(r1, d1a);
    Vector3 cross1b = Vector3::cross(r1, d1b);

    // scale
    d1a *= _inv_mass;
    d1b *= _inv_mass;
    cross1a = invJ * cross1a;
    cross1b = invJ * cross1b;

    // loop over the remaining contacts 
    for (unsigned j=i, jj=ii; j< q.contact_events.size(); j++, jj+= 2)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());
      bool negate = ((negate1 && !negate2) || (negate2 && !negate1));

      // compute the second moment arm
      Vector3 r2 = q.contact_events[j]->contact_point - _x;

      // get the two tangent vectors 
      const Vector3& d2a = q.contact_events[j]->contact_tan1;
      const Vector3& d2b = q.contact_events[j]->contact_tan2;

      // compute the cross products
      Vector3 cross2a = Vector3::cross(r2, d2a);
      Vector3 cross2b = Vector3::cross(r2, d2b);

      // compute entries of Dc_iM_DcT
      double sum1 = d1a.dot(d2a) + cross1a.dot(cross2a);
      double sum2 = d1a.dot(d2b) + cross1a.dot(cross2b);
      double sum3 = d1b.dot(d2a) + cross1b.dot(cross2a);
      double sum4 = d1b.dot(d2b) + cross1b.dot(cross2b);

      if (!negate)
      {
        q.Dc_iM_DcT(ii, jj)     += sum1; 
        q.Dc_iM_DcT(ii, jj+1)   += sum2; 
        q.Dc_iM_DcT(ii+1, jj)   += sum3; 
        q.Dc_iM_DcT(ii+1, jj+1) += sum4;
      }
      else
      {
        q.Dc_iM_DcT(ii, jj)     -= sum1; 
        q.Dc_iM_DcT(ii, jj+1)   -= sum2; 
        q.Dc_iM_DcT(ii+1, jj)   -= sum3; 
        q.Dc_iM_DcT(ii+1, jj+1) -= sum4;
      }
    }
  }

  // ensure symmetry on Dc_iM_DcT
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+=2)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // loop over the remaining contacts 
    for (unsigned j=i, jj=ii; j< q.contact_events.size(); j++, jj+= 2)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // enforce symmetry
      q.Dc_iM_DcT(jj,ii) = q.Dc_iM_DcT(ii,jj);
      q.Dc_iM_DcT(jj+1,ii) = q.Dc_iM_DcT(ii,jj+1);
      q.Dc_iM_DcT(jj,ii+1) = q.Dc_iM_DcT(ii+1,jj);
      q.Dc_iM_DcT(jj+1,ii+1) = q.Dc_iM_DcT(ii+1,jj+1);
    }
  } 
}
*/

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
    // check whether inner implicit joints are present (if none are present, 
    // this is a base link)
    bool is_base = true;
    BOOST_FOREACH(JointPtr j, _inner_joints)
    {
      if (j->get_constraint_type() == Joint::eImplicit)
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

  // check whether no implicit joints are present
  BOOST_FOREACH(JointPtr j, _inner_joints)
  {
    if (j->get_constraint_type() == Joint::eImplicit)
      return false;
  }

  // no implicit joints... it's the base
  return true;
}

/// Outputs the object state to the specified stream
/**
 * This method outputs all of the low-level details to the stream
 */
std::ostream& Moby::operator<<(std::ostream& out, const Moby::RigidBody& rb)
{
  // write the ID of the object
  out << "rigid body ID: " << rb.id << endl;

  // indicate whether the body is enabled
  out << "  enabled? " << rb.is_enabled() << endl;

  // write the computation frame
  out << "  computation frame: "; 
  switch (rb.get_computation_frame_type())
  {
    case eGlobal: out << "global" << endl; break;
    case eLink:   out << "link" << endl; break;
    case eJoint:  out << "joint" << endl; break;
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
  out << "  velocity (twist): " << rb.velocity() << endl;

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

