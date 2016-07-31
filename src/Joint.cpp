/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iostream>
#include <limits>
#include <queue>
#include <Ravelin/SpatialArithmeticd.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/XMLTree.h>
#include <Moby/RCArticulatedBody.h>

using std::vector;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
Joint::Joint() : Jointd()
{
  // initialize the implicit constraint stiffness
  implicit_constraint_stiffness = 0.2;

  // initialize boolean for determining to use inverse dynamics
  qd_des_set = false;

  // initialize spring joint bool to false
  spring_joint = false;

  // initialize friction coefficients
  mu_fc = mu_fv = (double) 0.0;

  // initialize restitution coefficient
  limit_restitution = (double) 0.0;

  // initialize the compliant layer depth
  compliant_layer_depth = 0.0;

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

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
  // initialize the implicit constraint stiffness
  implicit_constraint_stiffness = 0.2;

  // initialize boolean for determining to use inverse dynamics
  qd_des_set = false;

  // initialize spring joint bool to false
  spring_joint = false;

  // initialize friction coefficients
  mu_fc = mu_fv = (double) 0.0;

  // initialize restitution coefficient
  limit_restitution = (double) 0.0;

  // initialize the compliant layer depth
  compliant_layer_depth = 0.0;

  // set the inboard and outboard links
  _inboard_link = RigidBodyPtr(inboard);
  _outboard_link = RigidBodyPtr(outboard);

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

/// Applies impulse 

/// (Relatively slow) method for determining the joint velocity from current link velocities
void Joint::determine_q_dot()
{
  static LinAlgd LA;
  MatrixNd m, U, V;
  VectorNd S;

  // get the spatial axes
  const vector<SVelocityd>& s = get_spatial_axes();

  // convert to a matrix
  SpArithd::to_matrix(s, m);

  // get the velocities in computation frames
  shared_ptr<RigidBodyd> inboard = get_inboard_link();
  shared_ptr<RigidBodyd> outboard = get_outboard_link();
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
  if (!inboard)
    throw std::runtime_error("Inboard link is null!");

  // set the inboard link pointer
  _inboard_link = inboard;

  // setup F's pose relative to the inboard
  set_inboard_pose(inboard->_F, update_pose);
}

/// Sets the pointer to the outboard link for this joint
/**
 * \note also points the outboard link to this joint
 */
void Joint::set_outboard_link(RigidBodyPtr outboard, bool update_pose)
{
  if (!outboard)
    throw std::runtime_error("Outboard link is null!");

  // set the pointer
  _outboard_link = outboard;

  // set the outboard pose, if necessary 
  set_outboard_pose(outboard->_F, update_pose);
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

  // initialize desired state
  qd_des.set_zero(num_dof());

  // initialize max torque to infinity
  max_torque.set_one(num_dof());
  max_torque *= std::numeric_limits<double>::max();

  // initialize set point of spring
  set_point.set_zero(num_dof());

  // initialize stiffness of the spring
  stiffness.set_zero(num_dof());

  // initialize the damping of the spring
  damping.set_zero(num_dof());

  lolimit.set_one(NDOF) *= -std::numeric_limits<double>::max();
  hilimit.set_one(NDOF) *= std::numeric_limits<double>::max();
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
}

/// Gets the scaled actuator forces
/**
 * \note relevant only to reduced-coordinate articulated bodies
 */
VectorNd& Joint::get_scaled_force(VectorNd& f)
{
  // get the applied force
  f = force;

  return f;
}

/// Applies impulse using the constraint equations: J'*lambda
void Joint::apply_impulse(const VectorNd& lambda)
{
  // verify that the impulse is of the right dimension
  if (lambda.size() != num_constraint_eqns())
    throw std::runtime_error("Wrong size impulse applied");

  // get the inboard and outboard links
  shared_ptr<RigidBodyd> inboard = get_inboard_link();
  shared_ptr<RigidBodyd> outboard = get_outboard_link();

  // each Jacobian is of size joint equations x number of spatial coordinates
  // apply the impulse to the inboard body
  if (inboard->is_enabled())
  {
    MatrixNd J;
    VectorNd gf;
    calc_constraint_jacobian(true, J);
    J.transpose_mult(lambda, gf);
    inboard->apply_generalized_impulse(gf);
  }

  // apply the impulse to the inboard body
  if (outboard->is_enabled())
  {
    MatrixNd J;
    VectorNd gf;
    calc_constraint_jacobian(false, J);
    J.transpose_mult(lambda, gf);
    outboard->apply_generalized_impulse(gf);
  } 
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

  // save the id
  joint_id = id;

  // load the compliant layer depth, if specified
  XMLAttrib* cld_attr = node->get_attrib("compliant-layer-depth");
  if (cld_attr)
    compliant_layer_depth = cld_attr->get_real_value();

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
  }
  else
    _q_tare.set_zero(num_dof());

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
  }
}

/// Implements Base::save_to_xml()
void Joint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data 
  Visualizable::save_to_xml(node, shared_objects);

  // set the name (even though the joint will be subclassed)
  node->name = "Joint";

  // save the compliant layer depth
  node->attribs.insert(XMLAttrib("compliant-layer-depth", compliant_layer_depth));

  // save the joint limits
  node->attribs.insert(XMLAttrib("lower-limits", lolimit));
  node->attribs.insert(XMLAttrib("upper-limits", hilimit));

  // save the joint position and velocity
  node->attribs.insert(XMLAttrib("q", q));
  node->attribs.insert(XMLAttrib("qd", qd));
  node->attribs.insert(XMLAttrib("q-tare", _q_tare));

  // save the Coulomb and viscous friction coefficients
  node->attribs.insert(XMLAttrib("coulomb-friction-coeff", mu_fc));
  node->attribs.insert(XMLAttrib("viscous-friction-coeff", mu_fv));

  // save the restitution coefficient
  node->attribs.insert(XMLAttrib("restitution-coeff", limit_restitution));

  // save the inboard link, if any
  if (!_inboard_link.expired())
  {
    shared_ptr<RigidBodyd> inboard(_inboard_link);
    node->attribs.insert(XMLAttrib("inboard-link-id", inboard->body_id));
  }  

  // save the outboard link, if any
  if (!_outboard_link.expired())
  {
    shared_ptr<RigidBodyd> outboard(_outboard_link);
    node->attribs.insert(XMLAttrib("outboard-link-id", outboard->body_id));
  }  

  // if both inboard and outboard links are set, set the global position
  if (!_inboard_link.expired() && !_outboard_link.expired())
  {
    Origin3d loc(Pose3d::transform_point(GLOBAL, get_location()));
    node->attribs.insert(XMLAttrib("location", loc));
  }
}

