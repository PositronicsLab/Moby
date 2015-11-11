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
#include <Moby/Simulator.h>
#include <Moby/RigidBody.h>
#ifdef USE_OSG
#include "Color.h"
#include <osg/PositionAttitudeTransform>
#include <osg/PolygonMode>
#endif

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using boost::const_pointer_cast;
using std::vector;
using std::cerr;
using std::endl;
using std::map;
using std::list;
using std::queue;

Ravelin::LinAlgd LA_;

/// Default constructor
/**
 * Constructs a rigid body with zero mass, zero inertia tensor, and center
 * of mass at [0,0,0] with position at [0,0,0], identity orientation, and zero
 * linear and angular velocity.  Body is enabled by default.
 */
RigidBody::RigidBody()
{
  // setup visualization pose
  _vF->rpose = _F;
}

/// Applies a generalized impulse to the rigid body (calls the simulator)
void RigidBody::apply_generalized_impulse(const SharedVectorNd& gj)
{
  shared_ptr<Simulator> s(simulator);
  s->apply_impulse(RigidBodyd::get_this(), gj);
}

/// Sets the rigid body inertia for this body
void RigidBody::set_inertia(const SpatialRBInertiad& inertia)
{
  RigidBodyd::set_inertia(inertia);

  #ifdef VISUALIZE_INERTIA
    /// rigid body moment of inertia matrix (inertia aligned link COM frame)
    Ravelin::Origin3d _Jdiag;

    /// orientation of the inertia ellipsoid relative to inertial frame
    Ravelin::Matrix3d r_Jdiag;

    // determine the eigenvectors (principal inertial axes) and eigenvalues
    // ellipsoidal radii
    LA_.eig_symm_plus( r_Jdiag = inertia.J, _Jdiag);
    osg::Group* this_group = _vizdata->get_group();

    osg::Sphere* ellipse = new osg::Sphere( osg::Vec3(0,0,0), 1.0f);
    osg::ShapeDrawable* ellipse_draw = new osg::ShapeDrawable(ellipse);
    osg::Geode* ellipse_geode = new osg::Geode();
    ellipse_geode->addDrawable(ellipse_draw);

    osg::Node* n = ellipse_geode;
    CcolorVisitor  newColor;
    newColor.setColor( 1,0,1,0.1 );
    n->accept( newColor );

    // this code will show axes
//    n = x_axis_geode;
//    newColor.setColor( 1,0,0,0.25 );
//    n->accept( newColor );

//    n = y_axis_geode;
//    newColor.setColor( 0,1,0,0.25 );
//    n->accept( newColor );

//    n = z_axis_geode;
//    newColor.setColor( 0,0,1,0.25 );
//    n->accept( newColor );

    // Set to always wireframe
    osg::StateSet* stateset = new osg::StateSet;
    osg::PolygonMode* polymode = new osg::PolygonMode;
    polymode->setMode(osg::PolygonMode::FRONT_AND_BACK,osg::PolygonMode::LINE);
    stateset->setAttributeAndModes(polymode,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
    ellipse_geode->setStateSet(stateset);

    // setup the scaling, orientation, and translation
    osg::PositionAttitudeTransform *Transf = new osg::PositionAttitudeTransform();
    Transf->setScale(osg::Vec3(_Jdiag[0],_Jdiag[1],_Jdiag[2]));

    Ravelin::Quatd q_Jdiag(r_Jdiag);
    Ravelin::Transform3d iPose
        = Ravelin::Pose3d::calc_relative_pose(get_inertial_pose(),get_pose());
    q_Jdiag *= iPose.q;
    Transf->setAttitude(osg::Quat(q_Jdiag.x,q_Jdiag.y,q_Jdiag.z,q_Jdiag.w));
    Transf->setPosition(osg::Vec3(iPose.x[0],iPose.x[1],iPose.x[2]));
    Transf->addChild(ellipse_geode);

    // this code will show axes
//    Transf->addChild(x_axis_geode);
//    Transf->addChild(y_axis_geode);
//    Transf->addChild(z_axis_geode);

    // add to the groups
    this_group->removeChild(inertia_viz);
    inertia_viz = Transf;
    this_group->addChild(inertia_viz);
    #endif
}

/// Implements Base::load_from_xml()
void RigidBody::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  const unsigned X = 0, Y = 1, Z = 2;
  map<std::string, BasePtr>::const_iterator id_iter;

  // load parent data
  ControlledBody::load_from_xml(node, id_map);

  // save the id
  body_id = id;

  // ***********************************************************************
  // don't verify that the node is correct, b/c RigidBody can be subclassed
  // ***********************************************************************

  #ifdef USE_OSG
  /// Color to add to the rigid body when rendered
  Ravelin::VectorNd color_rgba;

  // Chnage the RGBA color of the link if provided
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

  // read whether the body is enabled, if provided
  XMLAttrib* enabled_attr = node->get_attrib("enabled");
  if (enabled_attr)
    set_enabled(enabled_attr->get_bool_value());

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
  XMLAttrib* aangle_attr = node->get_attrib("aangle");
  if (position_attr || rpy_attr || quat_attr || aangle_attr)
  {
    Pose3d T;
    if (position_attr)
      T.x = position_attr->get_origin_value();
    if (quat_attr)
      T.q = quat_attr->get_quat_value();
    else if (rpy_attr)
      T.q = rpy_attr->get_rpy_value();
    else if (aangle_attr)
    {
      VectorNd aa_vec;
      aangle_attr->get_vector_value(aa_vec);
      T.q = AAngled(aa_vec[0], aa_vec[1], aa_vec[2], aa_vec[3]);
    }
    set_pose(T);
  }

  // read the inertial frame here...
  XMLAttrib* com_attr = node->get_attrib("inertial-relative-com");
  XMLAttrib* J_rpy_attr = node->get_attrib("inertial-relative-rpy");
  XMLAttrib* J_aangle_attr = node->get_attrib("inertial-relative-aangle");
  XMLAttrib* J_quat_attr = node->get_attrib("inertial-relative-quat");
  if (com_attr || J_rpy_attr || J_aangle_attr || J_quat_attr)
  {
    shared_ptr<Pose3d> newjF(new Pose3d);
    newjF->rpose = _jF;
    
    // read the com
    if (com_attr)
      newjF->x = com_attr->get_origin_value();
    if (J_quat_attr)
      newjF->q = J_quat_attr->get_quat_value();
    else if (J_rpy_attr)
      newjF->q = J_rpy_attr->get_rpy_value();
    else if (J_aangle_attr)
    {
      VectorNd aa_vec;
      aangle_attr->get_vector_value(aa_vec);
      newjF->q = AAngled(aa_vec[0], aa_vec[1], aa_vec[2], aa_vec[3]);
    }

    // update newjF to refer to _jF's pose
    newjF->update_relative_pose(_jF->rpose);
    *_jF = *newjF;

    // update the mixed pose
    update_mixed_pose();
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
      cg->set_single_body(dynamic_pointer_cast<SingleBodyd>(SingleBodyd::shared_from_this()));

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
      XMLAttrib* rel_aangle_attr = (*i)->get_attrib("relative-aangle");
      XMLAttrib* rel_rpy_attr = (*i)->get_attrib("relative-rpy");
      if (rel_origin_attr)
        rTR->x = rel_origin_attr->get_origin_value();
      if (rel_rpy_attr)
        rTR->q = rel_rpy_attr->get_rpy_value();
      else if (rel_aangle_attr)
      {
        VectorNd aa_vec;
        aangle_attr->get_vector_value(aa_vec);
        rTR->q = AAngled(aa_vec[0], aa_vec[1], aa_vec[2], aa_vec[3]);
      }
      rTR->rpose = Fxx;
      Jx.pose = rTR;

      // transform the inertia and update the inertia for this
      J += Pose3d::transform(_jF, Jx);
    }

    // get the offset of J and subtract it from _jF
    shared_ptr<Pose3d> newjF(new Pose3d);
    newjF->x = J.h;
    newjF->rpose = _jF;
    SpatialRBInertiad Jnew = Pose3d::transform(newjF, J);

    // set _jF to be newJF;
    newjF->update_relative_pose(_jF->rpose);
    *_jF = *newjF;

    // update the mixed pose
    update_mixed_pose();

    // set the mass and inertia of the RigidBody additively
    set_inertia(Jnew);
  }

  // read the linear and/or velocity of the body, if provided
  XMLAttrib* lvel_attr = node->get_attrib("linear-velocity");
  XMLAttrib* avel_attr = node->get_attrib("angular-velocity");
  if (lvel_attr || avel_attr)
  {
    Vector3d lv = Vector3d::zero(), av = Vector3d::zero();
    SVelocityd v;
    v.pose = _F2;
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
  ControlledBody::save_to_xml(node, shared_objects);

  // get this as non-const for simplicity
  RigidBodyd* rb = (RigidBodyd*) this;

  // rename the node
  node->name = "RigidBody";

  // save whether the body is enabled
  node->attribs.insert(XMLAttrib("enabled", is_enabled()));

  // save whether the body is compliant
  node->attribs.insert(XMLAttrib("compliant", compliance == eCompliant));

  // save the mass
  node->attribs.insert(XMLAttrib("mass", rb->get_inertia().m));

   // save the inertia
  node->attribs.insert(XMLAttrib("inertia", rb->get_inertia().J));

  // convert the current pose to be with respect to global coordinates
  Pose3d F0 = *_F;
  F0.update_relative_pose(GLOBAL);
  node->attribs.insert(XMLAttrib("position", F0.x));
  node->attribs.insert(XMLAttrib("quat", F0.q));

  // save the inertial frame
  node->attribs.insert(XMLAttrib("inertial-relative-com", rb->get_inertial_pose()->x));
  node->attribs.insert(XMLAttrib("inertial-relative-quat", rb->get_inertial_pose()->q));

  // save the linear and angular velocities
  shared_ptr<Pose3d> TARGET(new Pose3d);
  TARGET->rpose = _F;
  TARGET->q = Quatd::invert(_F->q);
    shared_ptr<const Pose3d> TARGET_const = boost::const_pointer_cast<const Pose3d>(TARGET);
  SVelocityd v = Pose3d::transform(TARGET_const, rb->get_velocity());
  node->attribs.insert(XMLAttrib("linear-velocity", v.get_linear()));
  node->attribs.insert(XMLAttrib("angular-velocity", v.get_angular()));

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
    shared_ptr<ArticulatedBodyd> abody(_abody);
    node->attribs.insert(XMLAttrib("articulated-body-id", abody->body_id));
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

/// Determines whether the given link is a child link of this
bool RigidBody::is_child_link(shared_ptr<const RigidBody> query) const
{
  return RigidBodyd::is_child_link(dynamic_pointer_cast<const RigidBodyd>(query));
}

/// Determines whether the given link is a descendant of this
/**
 * \note returns <b>true</b> if query == this
 */
bool RigidBody::is_descendant_link(shared_ptr<const RigidBody> query) const
{
  return RigidBodyd::is_descendant_link(dynamic_pointer_cast<const RigidBodyd>(query));
}

/// Gets the first parent link of this link; returns NULL if there is no parent
RigidBodyPtr RigidBody::get_parent_link() const
{
  if (_inner_joints.size() > 1)
    throw std::runtime_error("Called RigidBody::get_parent_link() when multiple parent links present! It's not reasonable to call this method for links in maximal-coordinate articulated bodies.");

  return dynamic_pointer_cast<RigidBody>(RigidBodyd::get_parent_link());
}

/// Gets the explicit inner joint of this link; returns NULL if there is no explicit inner joint
/**
 * Throws an exception if this link has multiple explicit inner joints
 */
JointPtr RigidBody::get_inner_joint_explicit() const
{
  return dynamic_pointer_cast<Joint>(RigidBodyd::get_inner_joint_explicit());
}

/// Prepares to compute the ODE
void RigidBody::prepare_to_calc_ode(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the number of generalized coordinates
  const unsigned NGC_EUL = num_generalized_coordinates(DynamicBodyd::eEuler);

  // get the shared pointer to this
  RigidBodyPtr shared_this = dynamic_pointer_cast<RigidBody>(Base::shared_from_this());

  // get the generalized coordinates and velocity
  const SharedVectorNd gc = x.segment(0, NGC_EUL).get();
  const SharedVectorNd gv = x.segment(NGC_EUL, x.size()).get();

  // set the state and velocity
  set_generalized_coordinates_euler(gc);
  set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
  {
    VectorNd tmp;

    // get the clone as a dynamic body
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(clone);

    // update the clone
    get_generalized_coordinates_euler(tmp);
    db->set_generalized_coordinates_euler(tmp);    
    get_generalized_velocity(DynamicBodyd::eSpatial, tmp);
    db->set_generalized_velocity(DynamicBodyd::eSpatial, tmp);

    // get the generalized forces
    (*controller)(tmp, t, controller_arg);

    FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << id << std::endl;

    // apply the generalized forces
    add_generalized_force(tmp);
  }
}

/// Computes the ODE
void RigidBody::ode(double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the number of generalized coordinates
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the shared pointer to this
  RigidBodyPtr shared_this = dynamic_pointer_cast<RigidBody>(Base::shared_from_this());

  // get the derivative of generalized coordinates and velocity
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, dx.size());

  // we need the generalized velocity as Rodrigues coordinates
  get_generalized_velocity(DynamicBodyd::eEuler, dgc);

  // get the generalized acceleration
  get_generalized_acceleration(dgv);
}

/// Sets the articulated body corresponding to this body
/**
 * \param body a pointer to the articulated body or NULL if this body is
 *        not a link in an articulated body
 */
 void RigidBody::set_articulated_body(shared_ptr<ArticulatedBody> body) 
 { 
   _abody = dynamic_pointer_cast<ArticulatedBodyd>(body); 
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
  shared_ptr<ArticulatedBodyd> ab = rb.get_articulated_body();
  out << "  articulated body: " << ab->body_id << endl;

  // write all collision geometries
  out << "  collision geometries: ";
  BOOST_FOREACH(CollisionGeometryPtr g, rb.geometries)
    out << "    " << g << endl;

  return out;
}


