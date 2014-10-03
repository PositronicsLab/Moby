/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iostream>
#include <limits>
#include <queue>
#include <Moby/RigidBody.h>
#include <Moby/Spatial.h>
#include <Moby/Joint.h>
#include <Moby/XMLTree.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/NumericalException.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

// shared (static) linear algebra object
LinAlgd Joint::_LA;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
Joint::Joint() : Jointd()
{
  // make the constraint type unknown
  _constraint_type = eUnknown;

  // initialize friction coefficients
  mu_fc = mu_fv = (double) 0.0;

  // initialize restitution coefficient
  limit_restitution = (double) 0.0;

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

  // indicate that q tare does not need to be determined
  _determine_q_tare = false;

  // setup the visualization pose
  _vF->set_identity();
  _vF->rpose = _F;
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * \note does not set the inner joint for the outboard link or add the outboard
 *       as a child of the inboard link
 */
Joint::Joint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Jointd()
{
  // make the constraint type unknown
  _constraint_type = eUnknown;

  // initialize friction coefficients
  mu_fc = mu_fv = (double) 0.0;

  // initialize restitution coefficient
  limit_restitution = (double) 0.0;

  // set the inboard and outboard links
  _inboard_link = inboard;
  _outboard_link = outboard;

  // set poses
  RigidBodyPtr ib(inboard);
  RigidBodyPtr ob(outboard);
  set_inboard_pose(ib->get_pose(), false);
  set_outboard_pose(ob->_F, false);

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

  // setup the visualization pose
  _vF->rpose = _F;
}  

/// (Relatively slow) method for determining the joint velocity from current link velocities
void Joint::determine_q_dot()
{
  MatrixNd m, U, V;
  VectorNd S;

  // get the spatial axes
  const vector<SVelocityd>& s = get_spatial_axes();

  // convert to a matrix
  to_matrix(s, m);

  // compute the SVD
  _LA.svd(m, U, S, V);

  // get the velocities in computation frames
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  const SVelocityd& vi = inboard->get_velocity();
  const SVelocityd& vo = outboard->get_velocity();

  // get velocities in s's frame
  shared_ptr<const Pose3d> spose = get_pose();
  SVelocityd svi = Pose3d::transform(spose, vi);
  SVelocityd svo = Pose3d::transform(spose, vo);

  // compute the change in velocity
  m.mult(svo - svi, this->qd);
}

/// Sets the pointer to the inboard link for this joint (and updates the spatial axes, if the outboard link has been set)
void Joint::set_inboard_link(RigidBodyPtr inboard, bool update_pose)
{
  _inboard_link = inboard;
  if (!inboard)
    return;

  // add this joint to the outer joints
  inboard->_outer_joints.insert(get_this());

  // setup F's pose relative to the inboard
  set_inboard_pose(inboard->get_pose(), update_pose);

  // update articulated body pointers, if possible
  if (!inboard->get_articulated_body() && !_abody.expired())
    inboard->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (inboard->get_articulated_body() && _abody.expired())
    set_articulated_body(ArticulatedBodyPtr(inboard->get_articulated_body()));

  // the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1(inboard->get_articulated_body());
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
}

/// Sets the pointer to the outboard link for this joint
/**
 * \note also points the outboard link to this joint
 */
void Joint::set_outboard_link(RigidBodyPtr outboard, bool update_pose)
{
  _outboard_link = outboard;
  if (!outboard)
    return;

  // add this joint to the outer joints
  outboard->_inner_joints.insert(get_this());

  // get the outboard pose
  if (outboard->_F->rpose)
    throw std::runtime_error("Joint::set_inboard_link() - relative pose on inboard link already set");

  // setup Fb's pose relative to the outboard 
  set_outboard_pose(outboard->_F, update_pose);

  // setup the frame
  outboard->_xdj.pose = get_pose();
  outboard->_xddj.pose = get_pose();
  outboard->_Jj.pose = get_pose();
  outboard->_forcej.pose = get_pose();

  // use one articulated body pointer to set the other, if possible
  if (!outboard->get_articulated_body() && !_abody.expired())
    outboard->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (outboard->get_articulated_body() && _abody.expired())
    set_articulated_body(ArticulatedBodyPtr(outboard->get_articulated_body()));

  // the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1(outboard->get_articulated_body());
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
}

/// Sets the number of degrees-of-freedom for this joint
/**
 * \note resets all joint values (q, qd, qdd) to zero, limits to -/+ infinity,
 *       actuator forces to zero, and actuator limit to infinity.
 */
void Joint::init_data()
{
  // call parent method first
  Jointd::init_data();

  const unsigned NDOF = num_dof();
  const unsigned NEQ = num_constraint_eqns();

  qdd.set_zero(NDOF);
  force.set_zero(NDOF);
  maxforce.set_one(NDOF) *= std::numeric_limits<double>::max();
  lolimit.set_one(NDOF) *= -std::numeric_limits<double>::max();
  hilimit.set_one(NDOF) *= std::numeric_limits<double>::max();
  ff.set_zero(NDOF);
  lambda.set_zero(NEQ);
  _s.resize(NDOF);
}

/// Resets the force / torque produced by this joint's actuator and by the joint friction forces
void Joint::reset_force()
{
  force.set_zero(num_dof());
  ff.set_zero(num_dof());
}

/// Adds to the force / torque produced by this joint's actuator
void Joint::add_force(const VectorNd& force)
{
  this->force += force;
}

/// Sets the location of this joint with specified inboard and outboard links
void Joint::set_location(const Point3d& point, RigidBodyPtr inboard, RigidBodyPtr outboard) 
{
  assert(inboard && outboard);

  // convert p to the inboard and outboard links' frames
  Point3d pi = Pose3d::transform_point(inboard->get_pose(), point);
  Point3d po = Pose3d::transform_point(outboard->get_pose(), point);

  // set _F's and Fb's origins
  _F->x = Origin3d(pi);
  _Fb->x = Origin3d(po);

  // invalidate all outboard pose vectors
  outboard->invalidate_pose_vectors();

  // set inboard and outboard links
  set_inboard_link(inboard, false);
  set_outboard_link(outboard, false);

  // setup joint pointers
  if (inboard) inboard->add_outer_joint(get_this());
  if (outboard) outboard->add_inner_joint(get_this());
}

/// Gets the scaled and limited actuator forces
/**
 * \note relevant only to reduced-coordinate articulated bodies
 */
VectorNd& Joint::get_scaled_force(VectorNd& f)
{
  // get the (limited) applied force
  f = force;
  for (unsigned i=0; i< f.size(); i++)
    if (f[i] > maxforce[i])
      f[i] = maxforce[i];
    else if (f[i] < -maxforce[i])
      f[i] = -maxforce[i];

  return f;
}

/// Implements Base::load_from_xml()
void Joint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  std::map<std::string, BasePtr>::const_iterator id_iter;
  RigidBodyPtr inboard, outboard;

  // ***********************************************************************
  // don't verify that the node is correct, b/c Joint will 
  // be subclassed
  // ***********************************************************************
  
  // load the parent data
  Visualizable::load_from_xml(node, id_map);

  // read the lower limits, if given
  XMLAttrib* lolimit_attr = node->get_attrib("lower-limits");
  if (lolimit_attr)
  {
    lolimit_attr->get_vector_value(lolimit);
    if (lolimit.size() != num_dof())
      throw std::runtime_error("lower-limits read from XML does not match joint DOF");
  }

  // read the upper limits, if given
  XMLAttrib* hilimit_attr = node->get_attrib("upper-limits");
  if (hilimit_attr)
  {
    hilimit_attr->get_vector_value(hilimit);
    if (hilimit.size() != num_dof())
      throw std::runtime_error("upper-limits read from XML does not match joint DOF");
  }

  // read the maximum actuator force, if given
  XMLAttrib* maxforce_attr = node->get_attrib("max-forces");
  if (maxforce_attr)
  {
    maxforce_attr->get_vector_value(maxforce);
    if (maxforce.size() != num_dof())
      throw std::runtime_error("max-forces read from XML does not match joint DOF");
  }

  // read the joint positions, if given
  XMLAttrib* q_attr = node->get_attrib("q");
  if (q_attr)
  {
    q_attr->get_vector_value(q);
    if (q.size() != num_dof())
      throw std::runtime_error("q read from XML does not match joint DOF");
  }
  else
    q.set_zero(num_dof());

  // read the joint velocities, if given
  XMLAttrib* qd_attr = node->get_attrib("qd");
  if (qd_attr)
  {
    qd_attr->get_vector_value(qd);
    if (qd.size() != num_dof())
      throw std::runtime_error("qd read from XML does not match joint DOF");
  }

  // read the joint positions, if given
  XMLAttrib* q_init_attr = node->get_attrib("q-tare");
  if (q_init_attr)
  {
    q_init_attr->get_vector_value(_q_tare);
    if (_q_tare.size() != num_dof())
      throw std::runtime_error("q-tare read from XML does not match joint DOF");
    _determine_q_tare = false;
  }
  else
  {
    _q_tare.set_zero(num_dof());
    _determine_q_tare = true;
  }

  // read the Coulomb friction coefficient, if given
  XMLAttrib* fc_attr = node->get_attrib("coulomb-friction-coeff");
  if (fc_attr)
    mu_fc = fc_attr->get_real_value();

  // read the viscous friction coefficient, if given
  XMLAttrib* fv_attr = node->get_attrib("viscous-friction-coeff");
  if (fv_attr)
    mu_fv = fv_attr->get_real_value();

  // read the restitution coefficient, if given
  XMLAttrib* resti_attr = node->get_attrib("restitution-coeff");
  if (resti_attr)
    limit_restitution = resti_attr->get_real_value();

  // read the articulated body, if given
  XMLAttrib* ab_attr = node->get_attrib("articulated-body-id");
  if (ab_attr)
  {
    // get the ID
    const std::string& ID = ab_attr->get_string_value();

    // look for the ID -- only warn if it is not found
    if ((id_iter = id_map.find(ID)) == id_map.end())
    {
      #ifdef _DEBUG_XML_
      std::cout << "Joint::load_from_xml() warning - ";
      std::cout << "articulated body" << std::endl << "  '" << ID << "' not ";
      std::cout << "found" << std::endl << "  ** This warning could result ";
      std::cout << "from joints being constructed before articulated bodies ";
      std::cout << std::endl << "    and may not be serious..." << std::endl; 
      std::cout << "  offending node: " << std::endl << *node;
      #endif
    }
    else
      set_articulated_body(boost::dynamic_pointer_cast<RCArticulatedBody>(id_iter->second));
  }

  // read the inboard link id, if given
  XMLAttrib* inboard_attr = node->get_attrib("inboard-link-id");
  if (inboard_attr)
  {
    // get the ID of the inboard link
    const std::string& id = inboard_attr->get_string_value();

    // complain if the link not found but don't treat it as an error-- there 
    // are circular dependencies
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      #ifdef _DEBUG_XML_
      std::cout << "Joint::load_from_xml() warning - link ";
      std::cout << id << " not found" << std::endl;
      std::cout << "  ** This warning could result from joints being ";
      std::cout << " constructed before links" << std::endl << "     and may ";
      std::cout << "not be serious" << std::endl;
      std::cout << "  offending node: " << std::endl << *node;
      #endif
    }
    else
      inboard = boost::dynamic_pointer_cast<RigidBody>(id_iter->second);
  }

  // read the outboard link id, if given
  XMLAttrib* outboard_attr = node->get_attrib("outboard-link-id");
  if (outboard_attr)
  {
    // get the ID of the outboard link
    const std::string& id = outboard_attr->get_string_value();

    // complain if the link not found but don't treat it as an error-- there 
    // are circular dependencies
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      #ifdef _DEBUG_XML_
      std::cout << "Joint::load_from_xml() warning - link ";
      std::cout << id << " not found" << std::endl;
      std::cout << "  ** This warning could result from joints being ";
      std::cout << " constructed before links" << std::endl << "     and may ";
      std::cout << "not be serious" << std::endl;
      std::cout << "  offending node: " << std::endl << *node;
      #endif
    }
    else
      outboard = boost::dynamic_pointer_cast<RigidBody>(id_iter->second);
  }

  // get the global position of the joint, if possible
  XMLAttrib* pos_attr = node->get_attrib("location");
  if (pos_attr)
  {
    // get the position of the joint
    Point3d position;
    pos_attr->get_vector_value(position);

    // make sure that both inboard and outboard links have been set
    if (!inboard || !outboard)
    {
      std::cerr << "Joint::load_from_xml() - global position";
      std::cerr << " specified w/o " << std::endl << "  inboard and/or";
      std::cerr << " outboard links set!" << std::endl;
      return;
    }

    // set the joint location
    set_location(position, inboard, outboard);
  }
  else // no location specified
  {
    // set the inboard and outboard links, as specified
    if (inboard) set_inboard_link(inboard, true);
    if (outboard) set_outboard_link(outboard, true);

    // add/replace this as an inner joint
    if (inboard) inboard->add_outer_joint(get_this());
    if (outboard) outboard->add_inner_joint(get_this());
  }
}

/// Implements Base::save_to_xml()
void Joint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data 
  Visualizable::save_to_xml(node, shared_objects);

  // set the name (even though the joint will be subclassed)
  node->name = "Joint";

  // save the joint limits
  node->attribs.insert(XMLAttrib("lower-limits", lolimit));
  node->attribs.insert(XMLAttrib("upper-limits", hilimit));

  // save the maximum applicable forces by the joint's actuators
  node->attribs.insert(XMLAttrib("max-forces", maxforce));

  // save the joint position and velocity
  node->attribs.insert(XMLAttrib("q", q));
  node->attribs.insert(XMLAttrib("qd", qd));
  node->attribs.insert(XMLAttrib("q-tare", _q_tare));

  // save the Coulomb and viscous friction coefficients
  node->attribs.insert(XMLAttrib("coulomb-friction-coeff", mu_fc));
  node->attribs.insert(XMLAttrib("viscous-friction-coeff", mu_fv));

  // save the restitution coefficient
  node->attribs.insert(XMLAttrib("restitution-coeff", limit_restitution));

  // save the force currently applied to this joint
  node->attribs.insert(XMLAttrib("force", force));

  // save the inboard link, if any
  if (!_inboard_link.expired())
  {
    RigidBodyPtr inboard(_inboard_link);
    node->attribs.insert(XMLAttrib("inboard-link-id", inboard->id));
  }  

  // save the outboard link, if any
  if (!_outboard_link.expired())
  {
    RigidBodyPtr outboard(_outboard_link);
    node->attribs.insert(XMLAttrib("outboard-link-id", outboard->id));
  }  

  // if both inboard and outboard links are set, set the global position
  if (!_inboard_link.expired() && !_outboard_link.expired())
    node->attribs.insert(XMLAttrib("location", get_location()));

  // save the ID articulated body (if any)
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody(_abody);
    node->attribs.insert(XMLAttrib("articulated-body-id", abody->id));
  }
}

