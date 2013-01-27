/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <limits>
#include <queue>
#include <Moby/AAngle.h>
#include <Moby/Matrix4.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/RCArticulatedBody.h>

using namespace Moby;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
Joint::Joint()
{
  // make the constraint type unknown
  _constraint_type = eUnknown;

  // initialize friction coefficients
  mu_fc = mu_fv = (Real) 0.0;

  // initialize restitution coefficient
  limit_restitution = (Real) 0.0;

  // indicate that s0 is not valid
  _s0_valid = false;

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();

  // initialize _q_tare
  _q_tare.resize(0);
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
  mu_fc = mu_fv = (Real) 0.0;

  // initialize restitution coefficient
  limit_restitution = (Real) 0.0;

  // indicate that s0 is not valid
  _s0_valid = false;

  // set the inboard and outboard links
  _inboard_link = inboard;
  _outboard_link = outboard;

  // mark the indices as invalid initially
  _coord_idx = _joint_idx = _constraint_idx = std::numeric_limits<unsigned>::max();
}  

/// Determines q tare (if necessary)
void Joint::determine_q_tare()
{
  // see whether we need to determine q tare
  if (!_determine_q_tare)
    return;

  // determine q tare
  determine_q(_q_tare);

  // save current q
  VectorN q_save((const VectorN&) q);

  // we want the joint to induce identity transform, so set q to zero 
  q.set_zero();

  // get the global position - use the inboard link
  Vector3 position = get_position_global(false);

  // get the joint transform
  const Matrix4& Tj = get_transform(); 

  // get the inboard and outboard links
  RigidBodyPtr inboard(get_inboard_link());
  RigidBodyPtr outboard(get_outboard_link());

  // get the transforms for the two bodies
  const Matrix4& Ti = inboard->get_transform();
  const Matrix4& To = outboard->get_transform();

  // determine the vector from the inboard link to the joint (link coords)
  Vector3 inboard_to_joint = Ti.inverse_mult_point(position);

  // compute the vector from the joint to the outboard link in joint frame
  Vector3 joint_to_outboard_jf = Tj.inverse_mult_point((Ti.inverse_transform() * To).get_translation() - inboard_to_joint);    

  // update vector in joint frame
  RigidBody::InnerJointData& ijd = outboard->get_inner_joint_data(inboard);
  ijd.joint_to_com_vec_jf = joint_to_outboard_jf;

  // reset q
  q.copy_from(q_save);
}

/// (Relatively slow) method for determining the joint velocity from current link velocities
void Joint::determine_q_dot()
{
  // get the pseudo-inverse of the spatial axes
  MatrixN s;
  s.copy_from(get_spatial_axes(eGlobal));
  try
  {
    LinAlg::pseudo_inverse(s, LinAlg::svd1); 
  }
  catch (NumericalException e)
  {
    s.copy_from(get_spatial_axes(eGlobal));
    LinAlg::pseudo_inverse(s, LinAlg::svd2); 
  }

  // get the change in velocity
  RigidBodyPtr inboard = get_inboard_link();
  RigidBodyPtr outboard = get_outboard_link();
  SVector6 vi = inboard->get_spatial_velocity(eGlobal);
  SVector6 vo = outboard->get_spatial_velocity(eGlobal);
  s.mult(vo - vi, this->qd);
}

/// Evaluates the time derivative of the constraint
void Joint::evaluate_constraints_dot(Real C[6])
{
  Real Cx[6];

  // get the inboard and outboard links
  RigidBodyPtr in = get_inboard_link();
  RigidBodyPtr out = get_outboard_link();

  // get the linear angular velocities
  const Vector3& lvi = in->get_lvel();
  const Vector3& lvo = out->get_lvel();
  const Vector3& avi = in->get_avel();
  const Vector3& avo = out->get_avel();

  // compute
  const unsigned NEQ = num_constraint_eqns();
  for (unsigned i=0; i< NEQ; i++)
  {
    calc_constraint_jacobian(DynamicBody::eAxisAngle, in, i, Cx);
    Vector3 lv(Cx[0], Cx[1], Cx[2]);
    Vector3 av(Cx[3], Cx[4], Cx[5]);
    C[i] = lv.dot(lvi) + av.dot(avi);
    calc_constraint_jacobian(DynamicBody::eAxisAngle, out, i, Cx);
    lv = Vector3(Cx[0], Cx[1], Cx[2]);
    av = Vector3(Cx[3], Cx[4], Cx[5]);
    C[i] += -lv.dot(lvo) - av.dot(avo);
  }
}

/// Abstract method to update the local spatial axes
/**
 * Only applicable for reduced-coordinate articulated bodies
 */
void Joint::update_spatial_axes()
{
  // mark s0 as invalid
  _s0_valid = false;
}

/// Sets s bar from si
void Joint::calc_s_bar_from_si()
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned NDOF = num_dof();
  SAFESTATIC MatrixN ns;
  SAFESTATIC SMatrix6N sx;

  // transform sx to frame located at joint
  RigidBodyPtr outboard = get_outboard_link();
  if (!outboard)
    return;
  const Matrix4& To = outboard->get_transform();
  Vector3 x = get_position_global();
  SpatialTransform(To, IDENTITY_3x3, x).transform(_si, sx);

  // setup ns - it's the standard (i.e., non-spatial) transpose of sx
  assert(sx.columns() == NDOF);
  ns.resize(NDOF, SPATIAL_DIM);
  for (unsigned i=0; i< NDOF; i++)
    for (unsigned j=0; j< SPATIAL_DIM; j++)
      ns(i,j) = sx(j,i);

  // compute the nullspace
  LinAlg::nullspace(ns, _s_bar);
}

/// Sets the pointer to the inboard link for this joint (and updates the spatial axes, if the outboard link has been set)
void Joint::set_inboard_link(RigidBodyPtr inboard)
{
  _inboard_link = inboard;
  if (!inboard)
    return;

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
  maxforce.set_one(NDOF) *= std::numeric_limits<Real>::max();
  lolimit.set_one(NDOF) *= -std::numeric_limits<Real>::max();
  hilimit.set_one(NDOF) *= std::numeric_limits<Real>::max();
  ff.set_zero(NDOF);
  lambda.set_zero(NEQ);
  _si.resize(6, NDOF);
  _s0.resize(6, NDOF);
  _s_bar.resize(6, 6-NDOF);
}

/// Resets the force / torque produced by this joint's actuator and by the joint friction forces
void Joint::reset_force()
{
  force.set_zero(num_dof());
  ff.set_zero(num_dof());
}

/// Adds to the force / torque produced by this joint's actuator
void Joint::add_force(const VectorN& force)
{
  this->force += force;
}

/// Gets the global position of this joint
/**
 * \param use_outboard if <b>true</b> then the joint position is calculated 
 *        using the outboard link rather than inboard link; the position will
 *        not be identical if the joint constraint is violated (therefore,
 *        this method will behave identically for reduced-coordinate 
 *        articulated bodies)
 */
Vector3 Joint::get_position_global(bool use_outboard) const
{
  // get the inboard and outboard links
  RigidBodyPtr inboard(_inboard_link);
  RigidBodyPtr outboard(_outboard_link);

  // compute the global position
  if (!use_outboard)
  {
    const Matrix4& T = inboard->get_transform();
    const RigidBody::OuterJointData& o = inboard->get_outer_joint_data(outboard);
    return T.mult_point(o.com_to_joint_vec);
  }
  else
  {
    const Matrix4& T = outboard->get_transform();
    const RigidBody::InnerJointData& i = outboard->get_inner_joint_data(inboard);
    return T.mult_point(-i.joint_to_com_vec_of); 
  }
}

/// Resets the dynamics computation for this joint
/**
 * \note relevant only to reduced-coordinate articulated bodies
 */
void Joint::reset_spatial_axis()
{
  _s0_valid = false;
}

/// Gets the scaled and limited actuator forces
/**
 * \note relevant only to reduced-coordinate articulated bodies
 */
VectorN& Joint::get_scaled_force(VectorN& f)
{
  // get the (limited) applied force
  f.copy_from(force);
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
const SMatrix6N& Joint::get_spatial_axes(ReferenceFrameType rftype)
{
  if (rftype == eLink)
    return _si;
  else
  {
    if (!_s0_valid)
    {
      // NOTE: joint frame is linked to outboard link frame (see Mirtich
      //       joint and link numbering scheme 
      RigidBodyPtr outboard(_outboard_link);
      SpatialTransform X0i = outboard->get_spatial_transform_link_to_global();
      X0i.transform(_si, _s0);
      _s0_valid = true;
    }

    return _s0;
  }
}

/// Gets the complement of the spatial axes for this joint
/**
 * Spatial axes describe the motion of the joint. Spatial axes complement are
 * given in identity-oriented frame located at the origin. 
 */
const SMatrix6N& Joint::get_spatial_axes_complement()
{
  calc_s_bar_from_si();
  return _s_bar;
}

/// Gets the visualization transform for this joint
const Matrix4* Joint::get_visualization_transform()
{
  // make sure that there is an inboard link
  if (!get_inboard_link())
    return NULL;

  // get the inboard link
  RigidBodyPtr inboard(get_inboard_link());

  // get the transform for the inboard link
  _vtransform = inboard->get_transform();

  // set the translation for the joint visualization data to be its global position
  _vtransform.set_translation(get_position_global());

  return &_vtransform;
}

/// Gets the spatial constraints for this joint
SMatrix6N& Joint::get_spatial_constraints(ReferenceFrameType rftype, SMatrix6N& s)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Real Cq[7];

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
    calc_constraint_jacobian_rodrigues(outboard, i, Cq);

    // convert the differential quaternion constraints to an angular velocity
    // representation
    Vector3 omega = q.G_mult(Cq[3], Cq[4], Cq[5], Cq[6]) * (Real) 0.5;

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

/// Implements Base::load_from_xml()
void Joint::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  std::map<std::string, BasePtr>::const_iterator id_iter;

  // ***********************************************************************
  // don't verify that the node is correct, b/c Joint will 
  // be subclassed
  // ***********************************************************************
  
  // load the parent data
  Visualizable::load_from_xml(node, id_map);

  // read the lower limits, if given
  const XMLAttrib* lolimit_attr = node->get_attrib("lower-limits");
  if (lolimit_attr)
    lolimit_attr->get_vector_value(lolimit);

  // read the upper limits, if given
  const XMLAttrib* hilimit_attr = node->get_attrib("upper-limits");
  if (hilimit_attr)
    hilimit_attr->get_vector_value(hilimit);

  // read the maximum actuator force, if given
  const XMLAttrib* maxforce_attr = node->get_attrib("max-forces");
  if (maxforce_attr)
    maxforce_attr->get_vector_value(maxforce);

  // read the joint positions, if given
  const XMLAttrib* q_attr = node->get_attrib("q");
  if (q_attr)
    q_attr->get_vector_value(q);
  else
    q.set_zero(num_dof());

  // read the joint velocities, if given
  const XMLAttrib* qd_attr = node->get_attrib("qd");
  if (qd_attr)
    qd_attr->get_vector_value(qd);

  // read the joint positions, if given
  const XMLAttrib* q_init_attr = node->get_attrib("q-tare");
  if (q_init_attr)
  {
    _determine_q_tare = false;
    q_init_attr->get_vector_value(_q_tare);
  }
  else
    _determine_q_tare = true;

  // read the Coulomb friction coefficient, if given
  const XMLAttrib* fc_attr = node->get_attrib("coulomb-friction-coeff");
  if (fc_attr)
    mu_fc = fc_attr->get_real_value();

  // read the viscous friction coefficient, if given
  const XMLAttrib* fv_attr = node->get_attrib("viscous-friction-coeff");
  if (fv_attr)
    mu_fv = fv_attr->get_real_value();

  // read the restitution coefficient, if given
  const XMLAttrib* resti_attr = node->get_attrib("restitution-coeff");
  if (resti_attr)
    limit_restitution = resti_attr->get_real_value();

  // read the articulated body, if given
  const XMLAttrib* ab_attr = node->get_attrib("articulated-body-id");
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
  const XMLAttrib* inboard_attr = node->get_attrib("inboard-link-id");
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
    {
      RigidBodyPtr inboard(boost::dynamic_pointer_cast<RigidBody>(id_iter->second));
      set_inboard_link(inboard);
    }
  }

  // read the outboard link id, if given
  const XMLAttrib* outboard_attr = node->get_attrib("outboard-link-id");
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
    {
      RigidBodyPtr outboard(boost::dynamic_pointer_cast<RigidBody>(id_iter->second));
      set_outboard_link(outboard);
    }
  }

  // get the global position of the joint, if possible
  const XMLAttrib* pos_attr = node->get_attrib("global-position");
  if (pos_attr)
  {
    // get the position of the joint
    Vector3 position;
    pos_attr->get_vector_value(position);

    // make sure that both inboard and outboard links have been set
    if (!get_inboard_link() || !get_outboard_link())
    {
      std::cerr << "Joint::load_from_xml() - global position";
      std::cerr << " specified w/o " << std::endl << "  inboard and/or";
      std::cerr << " outboard links set!" << std::endl;
      return;
    }

    // get the inboard and outboard links
    RigidBodyPtr inboard(get_inboard_link());
    RigidBodyPtr outboard(get_outboard_link());

    // get the transforms for the two bodies
    const Matrix4& Ti = inboard->get_transform();
    const Matrix4& To = outboard->get_transform();

    // determine the vector from the inboard link to the joint (link coords)
    Vector3 inboard_to_joint = Ti.inverse_mult_point(position);

    // determine the vector from the joint to the outboard link (link coords)
    Vector3 joint_to_outboard_lf = -To.inverse_mult_point(position);

    // NOTE: the calculation immediately below assumes that the induced
    //       transform (i.e., the transform that the joint applies) is initally
    //       identity
    // compute the vector from the joint to the outboard link in joint frame
    Vector3 joint_to_outboard_jf = (Ti.inverse_transform() * To).get_translation() - inboard_to_joint;    

    // add/replace this as an inner joint
    inboard->add_outer_joint(outboard, get_this(), inboard_to_joint);
    outboard->add_inner_joint(inboard, get_this(), joint_to_outboard_jf, joint_to_outboard_lf);
  }

  // get the spatial axis in link coordinates; note that this must be done
  // after the link IDs are set
  const XMLAttrib* sa_link_attr = node->get_attrib("spatial-axis-link");
  if (sa_link_attr)
  {
    SMatrix6N si;
    sa_link_attr->get_matrix_value(si);
    if (si.columns() != num_dof())
      throw std::runtime_error("Incorrect spatial matrix size reading XML attribute spatial-axis-link");
    _si = si;
  }
}

/// Implements Base::save_to_xml()
void Joint::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
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
    node->attribs.insert(XMLAttrib("global-position", get_position_global()));

  // save the spatial axis in link coordinates
  node->attribs.insert(XMLAttrib("spatial-axis-link", _si));

  // save the ID articulated body (if any)
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody(_abody);
    node->attribs.insert(XMLAttrib("articulated-body-id", abody->id));
  }
}

