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
Joint::Joint()
{
  // make the constraint type unknown
  _constraint_type = eUnknown;

  // initialize friction coefficients
  mu_fc = mu_fv = (double) 0.0;

  // initialize restitution coefficient
  limit_restitution = (double) 0.0;

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

  // initialize _q_tare
  _q_tare.resize(0);

  // indicate that q tare does not need to be determined
  _determine_q_tare = false;

  // initialize the two frames
  _F = shared_ptr<Pose3d>(new Pose3d);
  _Fb = shared_ptr<Pose3d>(new Pose3d);
  _Fprime = shared_ptr<Pose3d>(new Pose3d);
  _Fprime->rpose = _F;

  // setup the visualization pose
  _vF->set_identity();
  _vF->rpose = _F;
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * \note does not set the inner joint for the outboard link or add the outboard
 *       as a child of the inboard link
 */
Joint::Joint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard)
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

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

  // initialize the two frames
  _F = shared_ptr<Pose3d>(new Pose3d);
  _Fb = shared_ptr<Pose3d>(new Pose3d);
  _Fprime = shared_ptr<Pose3d>(new Pose3d);
  _Fprime->rpose = _F;

  // setup the visualization pose
  _vF->rpose = _F;
}  

/// Sets the pose for the joint (this generally shouldn't be set by the user)
void Joint::set_pose(shared_ptr<const Pose3d> P)
{
  // verify that the relative pose is correct 
  if (_F->rpose != P->rpose)
    throw std::runtime_error("Joint::set_pose() - passed pose is not defined in correct frame");

  // update the frame
  *_F = *P;

  // copy the relative pose from _Fb (next operation will overwrite it)
  shared_ptr<const Pose3d> fb_rel = _Fb->rpose;

  // update _Fb
  *_Fb = *P;

  // reset _Fb's relative pose
  _Fb->update_relative_pose(fb_rel);
}

/// Determines q tare
void Joint::determine_q_tare()
{
  // determine q tare
  determine_q(_q_tare);
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

/// Evaluates the time derivative of the constraint
void Joint::evaluate_constraints_dot(double C[6])
{
  double Cx[6];

  // get the inboard and outboard links
  RigidBodyPtr in = get_inboard_link();
  RigidBodyPtr out = get_outboard_link();

  // get the linear angular velocities
  const SVelocityd& inv = in->get_velocity();
  const SVelocityd& outv = out->get_velocity();
  Vector3d lvi = inv.get_linear();
  Vector3d lvo = outv.get_linear();
  Vector3d avi = inv.get_angular();
  Vector3d avo = outv.get_angular();

  // compute
  const unsigned NEQ = num_constraint_eqns();
  for (unsigned i=0; i< NEQ; i++)
  {
    // TODO: fix this to do frame calculations
    calc_constraint_jacobian(DynamicBody::eSpatial, in, i, Cx);
    Vector3d lv(Cx[0], Cx[1], Cx[2]);
    Vector3d av(Cx[3], Cx[4], Cx[5]);
    C[i] = lv.dot(lvi) + av.dot(avi);
    calc_constraint_jacobian(DynamicBody::eSpatial, out, i, Cx);
    lv = Vector3d(Cx[0], Cx[1], Cx[2]);
    av = Vector3d(Cx[3], Cx[4], Cx[5]);
    C[i] += -lv.dot(lvo) - av.dot(avo);
  }
}

/// Abstract method to update the local spatial axes
/**
 * Only applicable for reduced-coordinate articulated bodies
 */
void Joint::update_spatial_axes()
{
  // setup the spatial axis vector and frame
  _s.resize(num_dof());
  for (unsigned i=0; i< _s.size(); i++)
    _s[i].pose = get_pose();
}

/// Sets s bar from si
void Joint::calc_s_bar_from_s()
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned NDOF = num_dof();
  MatrixNd sx, ns;

  // get s
  const vector<SVelocityd>& s = get_spatial_axes(); 

  // setup ns - it's the standard (i.e., non-spatial) transpose of s
  assert(s.size() == NDOF);
  to_matrix(s, sx);
  MatrixNd::transpose(sx, ns);

  // compute the nullspace
  _LA.nullspace(ns, sx);

  // setup s_bar
  from_matrix(sx, _s_bar);
  for (unsigned i=0; i< _s_bar.size(); i++)
    _s_bar[i].pose = get_pose();
}

/// Sets the pointer to the inboard link for this joint (and updates the spatial axes, if the outboard link has been set)
void Joint::set_inboard_link(RigidBodyPtr inboard)
{
  _inboard_link = inboard;
  if (!inboard)
    return;

  // add this joint to the outer joints
  inboard->_outer_joints.insert(get_this());

  // setup F's pose relative to the inboard
  _F->rpose = inboard->get_pose();

  // update spatial axes if both links are set
  if (!_outboard_link.expired() && !_inboard_link.expired())
    update_spatial_axes();

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
void Joint::set_outboard_link(RigidBodyPtr outboard)
{
  _outboard_link = outboard;
  if (!outboard)
    return;

  // add this joint to the outer joints
  outboard->_inner_joints.insert(get_this());

  // get the outboard pose
  shared_ptr<Pose3d> outboardF = outboard->_F; 
  assert(!outboardF->rpose);

  // setup Fb's pose relative to the outboard 
  _Fb->rpose = outboardF;
  outboardF->update_relative_pose(_Fprime);

  // setup the frame
  outboard->_xdj.pose = get_pose();
  outboard->_xddj.pose = get_pose();
  outboard->_Jj.pose = get_pose();
  outboard->_forcej.pose = get_pose();

  //update spatial axes if both links are set
  if (!_outboard_link.expired() && !_inboard_link.expired())
    update_spatial_axes();

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
  const unsigned NDOF = num_dof();
  const unsigned NEQ = num_constraint_eqns();

  q.set_zero(NDOF);
  _q_tare.set_zero(NDOF);
  qd.set_zero(NDOF);
  qdd.set_zero(NDOF);
  force.set_zero(NDOF);
  maxforce.set_one(NDOF) *= std::numeric_limits<double>::max();
  lolimit.set_one(NDOF) *= -std::numeric_limits<double>::max();
  hilimit.set_one(NDOF) *= std::numeric_limits<double>::max();
  ff.set_zero(NDOF);
  lambda.set_zero(NEQ);
  _s.resize(NDOF);
  _s_bar.resize(6-NDOF);
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

/// Sets the location of this joint
void Joint::set_location(const Point3d& point) 
{
  // verify inboard and outboard links are set
  if (_inboard_link.expired())
    throw std::runtime_error("Joint::set_location() called and inboard link not set");
  // verify inboard link is set
  if (_outboard_link.expired())
    throw std::runtime_error("Joint::set_location() called and outboard link not set");

  // get the inboard and outboard links
  RigidBodyPtr inboard(_inboard_link);
  RigidBodyPtr outboard(_outboard_link);

  // convert p to the inboard and outboard links' frames
  Point3d pi = Pose3d::transform_point(inboard->get_pose(), point);
  Point3d po = Pose3d::transform_point(outboard->get_pose(), point);

  // set _F's and Fb's origins
  _F->x = Origin3d(pi);
  _Fb->x = Origin3d(po);

  // invalidate all outboard pose vectors
  outboard->invalidate_pose_vectors();
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
  set_inboard_link(inboard);
  set_outboard_link(outboard);

  // setup joint pointers
  if (inboard) inboard->add_outer_joint(get_this());
  if (outboard) outboard->add_inner_joint(get_this());
}

/// Gets the location of this joint
/**
 * \param use_outboard if <b>true</b> then the joint position is calculated 
 *        using the outboard link rather than inboard link; the position will
 *        not be identical if the joint constraint is violated (therefore,
 *        this method will behave identically for reduced-coordinate 
 *        articulated bodies)
 */
Point3d Joint::get_location(bool use_outboard) const
{
  // get the outboard link
  RigidBodyPtr outboard(_outboard_link);

  // compute the global position
  if (!use_outboard)
  {
    // joint is defined with respect to inboard frame
    Point3d p(_F);
    p.set_zero();
    return p;
  }
  else
  {
    Point3d p(_Fb);
    p.set_zero();
    return p;
  }
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

/// Gets the spatial axes for this joint
/**
 * Spatial axes describe the motion of the joint. Note that for rftype = eLink,
 * spatial axes are given in outboard link's frame. 
 */
const vector<SVelocityd>& Joint::get_spatial_axes()
{
  return _s;
}

/// Gets the complement of the spatial axes for this joint
/**
 * Spatial axes describe the motion of the joint. Spatial axes complement are
 * given in identity-oriented frame located at the origin. 
 */
const vector<SVelocityd>& Joint::get_spatial_axes_complement()
{
  calc_s_bar_from_s();
  return _s_bar;
}

/*
/// Gets the spatial constraints for this joint
vector<SVelocityd>& Joint::get_spatial_constraints(ReferenceFrameType rftype, vector<SVelocityd>& s)
{
  const unsigned X = 0, Y = 1, Z = 2;
  double Cq[7];

  // get the outboard link and its orientation quaternion
  RigidBodyPtr outboard = get_outboard_link();
  assert(outboard);
  const Quat& q = outboard->get_orientation();

  // resize the spatial constraint matrix
  s.resize(6, num_constraint_eqns());

  // calculate the constraint Jacobian in relation to the outboard link
  for (unsigned i=0; i< num_constraint_eqns(); i++)
  {
    // calculate the constraint Jacobian
    calc_constraint_jacobian(outboard, i, Cq);

    // convert the differential quaternion constraints to an angular velocity
    // representation
    Vector3d omega = q.G_mult(Cq[3], Cq[4], Cq[5], Cq[6]) * (double) 0.5;

    // setup the column of the constraint Jacobian
    s(0,i) = omega[X];
    s(1,i) = omega[Y];
    s(2,i) = omega[Z];
    s(3,i) = Cq[X];
    s(4,i) = Cq[Y];
    s(5,i) = Cq[Z];
  }

  // TODO: this should be in link-global frame -- fix this!
  assert(false);

  // convert to link frame if necessary (constraints computed in global frame)
  if (rftype == eLink)
  {
    SpatialTransform Xi0 = outboard->get_spatial_transform_global_to_link();
    for (unsigned i=0; i< num_constraint_eqns(); i++)
    {
      SVector6 scol = s.get_column(i);
      scol = Xi0.transform(scol);
      s.set_column(i, scol);
    }
  }

  return s;
}
*/

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
    lolimit_attr->get_vector_value(lolimit);

  // read the upper limits, if given
  XMLAttrib* hilimit_attr = node->get_attrib("upper-limits");
  if (hilimit_attr)
    hilimit_attr->get_vector_value(hilimit);

  // read the maximum actuator force, if given
  XMLAttrib* maxforce_attr = node->get_attrib("max-forces");
  if (maxforce_attr)
    maxforce_attr->get_vector_value(maxforce);

  // read the joint positions, if given
  XMLAttrib* q_attr = node->get_attrib("q");
  if (q_attr)
    q_attr->get_vector_value(q);
  else
    q.set_zero(num_dof());

  // read the joint velocities, if given
  XMLAttrib* qd_attr = node->get_attrib("qd");
  if (qd_attr)
    qd_attr->get_vector_value(qd);

  // read the joint positions, if given
  XMLAttrib* q_init_attr = node->get_attrib("q-tare");
  if (q_init_attr)
  {
    q_init_attr->get_vector_value(_q_tare);
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
    if (inboard) set_inboard_link(inboard);
    if (outboard) set_outboard_link(outboard);

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

