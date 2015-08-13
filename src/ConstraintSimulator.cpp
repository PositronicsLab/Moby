/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <unistd.h>
#include <boost/tuple/tuple.hpp>
#include <Moby/XMLTree.h>
#include <Moby/Dissipation.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/ControlledBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/SustainedUnilateralConstraintSolveFailException.h>
#include <Moby/InvalidStateException.h>
#include <Moby/InvalidVelocityException.h>
#include <Moby/ConstraintStabilization.h>
#include <Moby/ConstraintSimulator.h>

#ifdef USE_OSG
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/PositionAttitudeTransform>
#include <osg/Quat>
#endif // USE_OSG

using std::endl;
using std::set;
using std::list;
using std::vector;
using std::map;
using std::make_pair;
using std::multimap;
using std::pair;
using boost::tuple;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Default constructor
ConstraintSimulator::ConstraintSimulator()
{
  constraint_callback_fn = NULL;
  constraint_post_callback_fn = NULL;
  post_mini_step_callback_fn = NULL;
  get_contact_parameters_callback_fn = NULL;
  render_contact_points = false;

  // setup contact distance thresholds
  contact_dist_thresh = 1e-6;

  // setup the collision detector
  _coldet = shared_ptr<CollisionDetection>(new CCD);
}

/// Gets the contact data between a pair of geometries (if any)
/**
 * This method looks for contact data not only between the pair of geometries, but also
 * the rigid bodies that the geometries belong to, and any articulated bodies as well.
 * The search proceeds in the following manner: <br />
 * <ol>
 *  <li>two collision geometries</li>
 *  <li>one collision geometry, one rigid body</li>
 *  <li>two rigid bodies</li>
 *  <li>one collision geometry, one articulated body</li>
 *  <li>one rigid body, one articulated body</li>
 *  <li>two articulated bodies</li>
 * </ol>
 * The search order allows for multiple granularities; for example, a collision can easily
 * be specified between two geometries of two of a robot's links (i.e., representing different
 * surfaces on the links), between two links, or between two robots.
 * \param g1 the first collision geometry
 * \param g2 the second collision geometry
 * \return a pointer to the contact data, if any, found
 */
shared_ptr<ContactParameters> ConstraintSimulator::get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const
{
  map<sorted_pair<BasePtr>, shared_ptr<ContactParameters> >::const_iterator iter;

  // first see whether a user function for contact parameters is defined,
  // and if it is defined, attempt to get contact parameters from it
  if (get_contact_parameters_callback_fn)
  {
    shared_ptr<ContactParameters> cp = get_contact_parameters_callback_fn(geom1, geom2);
    if (cp)
      return cp;
  }

  // search for the two contact geometries first
  if ((iter = contact_params.find(make_sorted_pair(geom1, geom2))) != contact_params.end())
    return iter->second;

  // get the geometries as base pointers
  BasePtr g1(geom1);
  BasePtr g2(geom2);

  // get the two single bodies
  assert(geom1->get_single_body());
  assert(geom2->get_single_body());
  RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(geom1->get_single_body());
  RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(geom2->get_single_body());
  BasePtr sb1 = rb1;
  BasePtr sb2 = rb2;

  // get the articulated bodies, if any
  BasePtr ab1, ab2;
  if (rb1)
    ab1 = dynamic_pointer_cast<ArticulatedBody>(rb1->get_articulated_body());
  if (rb2)
    ab2 = dynamic_pointer_cast<ArticulatedBody>(rb2->get_articulated_body());

  // check collision geometry 2 and rigid body 2 against articulated body 1
  if (ab1)
  {
    if ((iter = contact_params.find(make_sorted_pair(g2, ab1))) != contact_params.end())
      return iter->second;
    if ((iter = contact_params.find(make_sorted_pair(sb2, ab1))) != contact_params.end())
      return iter->second;
  }

  // check collision geometry 1 and rigid body 1 against articulated body 2
  if (ab2)
  {
    if ((iter = contact_params.find(make_sorted_pair(g1, ab2))) != contact_params.end())
      return iter->second;
    if ((iter = contact_params.find(make_sorted_pair(sb1, ab2))) != contact_params.end())
      return iter->second;
  }

  // check the two articulated bodies against articulated body 2
  if (ab1 && ab2)
    if ((iter = contact_params.find(make_sorted_pair(ab1, ab2))) != contact_params.end())
      return iter->second;

  // search for contact geometry 1 and rigid body 2
  if ((iter = contact_params.find(make_sorted_pair(g1, sb2))) != contact_params.end())
    return iter->second;

  // search for contact geometry 2 and rigid body 1
  if ((iter = contact_params.find(make_sorted_pair(g2, sb1))) != contact_params.end())
    return iter->second;

  // search for both rigid bodies
  if ((iter = contact_params.find(make_sorted_pair(sb1, sb2))) != contact_params.end())
    return iter->second;

  // still here?  no contact data found
  return shared_ptr<ContactParameters>();
}

/// Draws a ray directed from a contact point along the contact normal
void ConstraintSimulator::visualize_contact( UnilateralConstraint& constraint ) 
{
  #ifdef USE_OSG

  // random color for this contact visualization
  double r = (double) rand() / (double) RAND_MAX;
  double g = (double) rand() / (double) RAND_MAX;
  double b = (double) rand() / (double) RAND_MAX;
  osg::Vec4 color = osg::Vec4( r, g, b, 1.0 );

  // knobs for tweaking
  const double point_radius = 0.75;
  const double point_scale = 0.01;
  const double line_length = 5.0;
  const double line_radius = 0.1;
  const double head_radius = 0.5;
  const double head_height = 2.0;

  // the osg node this constraint visualization will attach to
  osg::Group* contact_root = new osg::Group();

  // turn off lighting for this node
  osg::StateSet *contact_state = contact_root->getOrCreateStateSet();
  contact_state->setMode( GL_LIGHTING, osg::StateAttribute::PROTECTED | osg::StateAttribute::OFF );

  // a geode for the visualization geometry
  osg::Geode* contact_geode = new osg::Geode();

  // add some hints to reduce the polygonal complexity of the visualization
  osg::TessellationHints *hints = new osg::TessellationHints();
  hints->setTessellationMode( osg::TessellationHints::USE_TARGET_NUM_FACES );
  hints->setCreateNormals( true );
  hints->setDetailRatio( 0.2 );

  // add the contact point as a sphere at the origin of the geode's frame
  osg::Sphere* point_geometry = new osg::Sphere( osg::Vec3( 0, 0, 0 ), point_radius );
  osg::ShapeDrawable* point_shape = new osg::ShapeDrawable( point_geometry, hints );
  point_shape->setColor( color );
  contact_geode->addDrawable( point_shape );

  // add the contact normal as a cylinder in the geode's frame
  osg::Cylinder* line_geometry = new osg::Cylinder( osg::Vec3( 0.0, 0.0, line_length / 2 ), line_radius, line_length );
  osg::ShapeDrawable* line_shape = new osg::ShapeDrawable( line_geometry, hints );
  line_shape->setColor( color );
  contact_geode->addDrawable( line_shape );

  // add the arrow head as a cone in the geode's frame
  osg::Cone* head_geometry = new osg::Cone( osg::Vec3( 0, 0, line_length ), head_radius, head_height );
  osg::ShapeDrawable* head_shape = new osg::ShapeDrawable( head_geometry, hints );
  head_shape->setColor( color );
  contact_geode->addDrawable( head_shape );

  // calculate the orientation based upon the direction of the normal vector.
  // Note: the default orientation of the osg model is along the z-axis
  double theta;
  Vector3d z = Vector3d( 0.0, 0.0, 1.0 );
  Vector3d axis = Vector3d::cross( constraint.contact_normal, z );
  if( axis.norm_inf() < NEAR_ZERO) {
    // z and normal are parallel, axis ill defined
    if( constraint.contact_normal[2] > 0 ) {
      // normal is z
      axis = Vector3d( 0.0, 1.0, 0.0 );
      theta = 0.0;
    } else {
      // normal is -z
      axis = Vector3d( 0.0, 1.0, 0.0 );
      theta = osg::PI;
    }
  } else {
    // axis is well defined
    axis = Vector3d::normalize(axis);
    theta = -std::acos( Vector3d::dot( constraint.contact_normal, z ) );
    // Note: theta calculation works but is not robust, could be rotated in opposite direction
  }
  osg::Quat q = osg::Quat( axis[0]*std::sin(theta/2), axis[1]*std::sin(theta/2), axis[2]*std::sin(theta/2), std::cos(theta/2) );

  // create the visualization transform
  osg::PositionAttitudeTransform* contact_transform = new osg::PositionAttitudeTransform();
  contact_transform->setPosition( osg::Vec3( constraint.contact_point[0], constraint.contact_point[1], constraint.contact_point[2] ) );
  contact_transform->setScale( osg::Vec3( point_scale, point_scale, point_scale ) );
  contact_transform->setAttitude( q );

  // add the geode to the transform
  contact_transform->addChild( contact_geode );

  // add the transform to the root
  contact_root->addChild( contact_transform );

  // add the root to the transient data scene graph
  add_transient_vdata( contact_root );

  // JRT : remove validator once theta 100% proven
  // -----------------------------------------
  // Rotational Validator
  // -----------------------------------------

  // Validator is a simple sphere translated along the normal
  // such that the visualization above should point at the center
  // of the validator.  If it doesn't, then the calculation of
  // theta in the rotational code above needs correction for that case

  // knobs for tweaking
  const double validator_scale = point_scale / 3;
  const double validator_ray_length = line_length * 2.5;

  // a root for the validator
  osg::Group* validator_root = new osg::Group();

  // turn off lighting for this node
  osg::StateSet *validator_state = validator_root->getOrCreateStateSet();
  validator_state->setMode( GL_LIGHTING, osg::StateAttribute::PROTECTED | osg::StateAttribute::OFF );

  // colocate the validator position to the contact point
  osg::PositionAttitudeTransform* validator_transform = new osg::PositionAttitudeTransform();
  validator_transform->setPosition( osg::Vec3( constraint.contact_point[0], constraint.contact_point[1], constraint.contact_point[2] ) );
  validator_transform->setScale( osg::Vec3( validator_scale, validator_scale, validator_scale ) );
  validator_root->addChild( validator_transform );

  // validator geometry
  osg::Sphere* validator_geometry = new osg::Sphere( osg::Vec3( 0, 0, 0 ), 1.0 );
  osg::ShapeDrawable* validator_shape = new osg::ShapeDrawable( validator_geometry, hints );
  validator_shape->setColor( color );

  // validator transform follows the normal out to a distance of validator_ray_length
  // Note: the validator is not rotated at all.  It is translated from the point along the normal
  osg::PositionAttitudeTransform* validator_end_transform = new osg::PositionAttitudeTransform();
  validator_end_transform->setPosition( osg::Vec3( constraint.contact_normal[0] * validator_ray_length, constraint.contact_normal[1] * validator_ray_length, constraint.contact_normal[2] * validator_ray_length ) );
  validator_transform->addChild( validator_end_transform );

  // add all validator constituents to the group
  osg::Geode* validator_geode = new osg::Geode();
  validator_transform->addChild( validator_end_transform );
  validator_end_transform->addChild( validator_geode );
  validator_geode->addDrawable( validator_shape );
  add_transient_vdata( validator_root );

  #endif // USE_OSG
}

/// Handles constraints
void ConstraintSimulator::calc_impacting_unilateral_constraint_forces(double dt)
{
  // if there are no constraints, quit now
  if (_rigid_constraints.empty())
    return;

  // call the callback function, if any
  if (constraint_callback_fn)
    (*constraint_callback_fn)(_rigid_constraints, constraint_callback_data);

  // preprocess constraints
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    preprocess_constraint(_rigid_constraints[i]);

  // look for the case where there are no impacting constraints
  bool none_impacting = true;
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    if (_rigid_constraints[i].determine_constraint_class() == UnilateralConstraint::eNegative)
    {
      none_impacting = false;
      break;
    }

  // if there are no impacts, return
  if (none_impacting)
    return;

  // if the setting is enabled, draw all contact constraints
  if( render_contact_points ) {
    for ( std::vector<UnilateralConstraint>::iterator it = _rigid_constraints.begin(); it < _rigid_constraints.end(); it++ ) {
      UnilateralConstraint& constraint = *it;
      if( constraint.constraint_type != UnilateralConstraint::eContact ) continue;
      visualize_contact( constraint );
    }
  }

  // compute impulses here...
  try
  {
    _impact_constraint_handler.process_constraints(_rigid_constraints);
  }
  catch (ImpactToleranceException e)
  {
    #ifndef NDEBUG
    std::cerr << "warning: impacting constraint tolerances exceeded" << std::endl;
    #endif
  }

  // call the post application callback, if any
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_rigid_constraints, constraint_post_callback_data);
}

/// Computes compliant contact forces 
void ConstraintSimulator::calc_compliant_unilateral_constraint_forces()
{
  // if there are no compliant constraints, quit now
  if (_compliant_constraints.empty())
    return;

  // call the callback function, if any
  if (constraint_callback_fn)
    (*constraint_callback_fn)(_compliant_constraints, constraint_callback_data);

  // preprocess constraints
  for (unsigned i=0; i< _compliant_constraints.size(); i++)
    preprocess_constraint(_compliant_constraints[i]);

  // if the setting is enabled, draw all contact constraints
  if( render_contact_points ) {
    for ( std::vector<UnilateralConstraint>::iterator it = _compliant_constraints.begin(); it < _compliant_constraints.end(); it++ ) {
      UnilateralConstraint& constraint = *it;
      if( constraint.constraint_type != UnilateralConstraint::eContact ) continue;
      visualize_contact( constraint );
    }
  }

  // compute contact constraint penalty forces here...
  _penalty_constraint_handler.process_constraints(_compliant_constraints);

  // call the post application callback, if any 
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_compliant_constraints, constraint_post_callback_data);
}

/// Performs necessary preprocessing on an constraint
void ConstraintSimulator::preprocess_constraint(UnilateralConstraint& e)
{
  // no pre-processing for limit constraints currently...
  if (e.constraint_type == UnilateralConstraint::eLimit)
    return;

  // no pre-processing for (none) constraints
  if (e.constraint_type == UnilateralConstraint::eNone)
    return;

  // get the contact parameters
  assert(e.constraint_type == UnilateralConstraint::eContact);
  shared_ptr<ContactParameters> cparams = get_contact_parameters(e.contact_geom1, e.contact_geom2);
  if (cparams)
    e.set_contact_parameters(*cparams);
  else
  {
    shared_ptr<SingleBodyd> sb1(e.contact_geom1->get_single_body());
    shared_ptr<SingleBodyd> sb2(e.contact_geom2->get_single_body());
    std::cerr << "ConstraintSimulator::preprocess_constraint() warning- no contact ";
    std::cerr << "data for contact" << std::endl;
    std::cerr << "  between " << e.contact_geom1->id << " (body ";
    std::cerr << sb1->body_id << ") and " << e.contact_geom2->id;
    std::cerr << " (body " << sb2->body_id << ")" << std::endl;
    std::cerr << "  ... ignoring" << std::endl;
  }
}

/// Sets up the list of collision geometries
void ConstraintSimulator::determine_geometries()
{
  // clear the list at first
  _geometries.clear();

  // determine all geometries
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
    if (rb)
      _geometries.insert(_geometries.end(), rb->geometries.begin(), rb->geometries.end());
    else
    {
      ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
      BOOST_FOREACH(shared_ptr<RigidBodyd> rbd, ab->get_links())
      {
        RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(rbd);
        _geometries.insert(_geometries.end(), rb->geometries.begin(), rb->geometries.end());
      }
    }
  }

  // sort and remove duplicates
  std::sort(_geometries.begin(), _geometries.end());
  _geometries.erase(std::unique(_geometries.begin(), _geometries.end()), _geometries.end());
}

/// Computes pairwise distances of geometries at their current poses, using broad phase results to determine which pairs should be checked
/**
 * \param pairwise_distances on return, contains the pairwise distances
 */
void ConstraintSimulator::calc_pairwise_distances()
{
  // clear the vector
  _pairwise_distances.clear();

  for (unsigned i=0; i< _pairs_to_check.size(); i++)
  {
    PairwiseDistInfo pdi;
    pdi.a = _pairs_to_check[i].first;
    pdi.b = _pairs_to_check[i].second;
    pdi.dist = _coldet->calc_signed_dist(pdi.a, pdi.b, pdi.pa, pdi.pb);
    Pose3d poseA(*pdi.a->get_pose());
    Pose3d poseB(*pdi.b->get_pose());
    poseA.update_relative_pose(GLOBAL);
    poseB.update_relative_pose(GLOBAL);
    FILE_LOG(LOG_SIMULATOR) << "ConstraintSimulator::calc_pairwise_distances() - signed distance between " << pdi.a->get_single_body()->body_id << " and " << pdi.b->get_single_body()->body_id << ": " << pdi.dist << std::endl;
    _pairwise_distances.push_back(pdi);
  }
}

/// Does broad phase collision detection, identifying which pairs of geometries may come into contact over time step of dt
void ConstraintSimulator::broad_phase(double dt)
{
  // call the broad phase
  _coldet->broad_phase(dt, _bodies, _pairs_to_check);

  // remove pairs that are unchecked
  for (unsigned i=0; i< _pairs_to_check.size(); )
    if (unchecked_pairs.find(make_sorted_pair(_pairs_to_check[i].first, _pairs_to_check[i].second)) != unchecked_pairs.end())
    {
      _pairs_to_check[i] = _pairs_to_check.back();
      _pairs_to_check.pop_back();
    }
    else
      i++;
}

/// Updates constraint violation at beginning of integration step
void ConstraintSimulator::update_constraint_violations(const vector<PairwiseDistInfo>& pairwise_distances)
{
  FILE_LOG(LOG_SIMULATOR) << "ConstraintSimulator::update_constraint_violations() entered" << std::endl;

  // set possible constraint violation
  BOOST_FOREACH(const PairwiseDistInfo& pdi, pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

    _ip_tolerances[make_sorted_pair(pdi.a, pdi.b)] = std::min(pdi.dist, -NEAR_ZERO);
  }

  // update joint constraint interpenetration
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
    if (ab)
      ab->update_joint_constraint_violations();
  }

  FILE_LOG(LOG_SIMULATOR) << "ConstraintSimulator::update_constraint_violations() exited" << std::endl;
}

/// Computes forward dynamics for all bodies
void ConstraintSimulator::calc_fwd_dyn()
{
  // clear force accumulators, then add all recurrent and compliant
  // constraint forces
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    // get the body as a Ravelin dynamic body
    shared_ptr<DynamicBodyd> rdb = dynamic_pointer_cast<DynamicBodyd>(db);

    // clear the force accumulators on the body
    rdb->reset_accumulators();

    // add all recurrent forces on the body
    const list<RecurrentForcePtr>& rfs = db->get_recurrent_forces();
    BOOST_FOREACH(RecurrentForcePtr rf, rfs)
      rf->add_force(rdb);
    
    // call the body's controller
    if (db->controller)
    {
      FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << db->id << std::endl;
      (*db->controller)(db, current_time, db->controller_arg);
    }
  }

  // calculate compliant constraint forces
  calc_compliant_unilateral_constraint_forces();

  // compute controller forces and call forward dynamics
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    // get the body as a Ravelin dynamic body
    shared_ptr<DynamicBodyd> rdb = dynamic_pointer_cast<DynamicBodyd>(db);

    // calculate forward dynamics at state x
    rdb->calc_fwd_dyn();
  }
}

/// Finds the set of unilateral constraints
void ConstraintSimulator::find_unilateral_constraints(double contact_dist_thresh)
{
  FILE_LOG(LOG_SIMULATOR) << "ConstraintSimulator::find_unilateral_constraints() entered" << std::endl;

  // clear the vectors of constraints
  _rigid_constraints.clear();
  _compliant_constraints.clear();

  // process each articulated body, getting joint constraints
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // see whether the i'th body is articulated
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(_bodies[i]);
    if (!ab)
      continue;

    // if the body is kinematically controlled, do nothing
    if (ab->get_kinematic())
      continue;

    // get limit constraints
    ab->find_limit_constraints(std::back_inserter(_rigid_constraints));
  }

  // find contact constraints
  BOOST_FOREACH(const PairwiseDistInfo& pdi, _pairwise_distances)
    if (pdi.dist < contact_dist_thresh)
    {
      // see whether one of the bodies is compliant
      RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
      RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
      if (rba->compliance == RigidBody::eCompliant || 
          rbb->compliance == RigidBody::eCompliant)
        _coldet->find_contacts(pdi.a, pdi.b, _compliant_constraints);
      else        
        _coldet->find_contacts(pdi.a, pdi.b, _rigid_constraints, contact_dist_thresh);
    }

  // set constraints to proper type
  for (unsigned i=0; i< _compliant_constraints.size(); i++)
    _compliant_constraints[i].compliance = UnilateralConstraint::eCompliant;
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    _rigid_constraints[i].compliance = UnilateralConstraint::eRigid;

  if (LOGGING(LOG_SIMULATOR))
  {
    for (unsigned i=0; i< _rigid_constraints.size(); i++)
    FILE_LOG(LOG_SIMULATOR) << _rigid_constraints[i] << std::endl;
    for (unsigned i=0; i< _compliant_constraints.size(); i++)
    FILE_LOG(LOG_SIMULATOR) << _compliant_constraints[i] << std::endl;
  }

  FILE_LOG(LOG_SIMULATOR) << "ConstraintSimulator::find_unilateral_constraints() exited" << std::endl;
}

/// Implements Base::load_from_xml()
void ConstraintSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // do not verify node name b/c this is abstract class
  // assert(strcasecmp(node->name.c_str(), "ConstraintSimulator") == 0);

  // first, load all data specified to the Simulator object
  Simulator::load_from_xml(node, id_map);

  // read the dissipator, if any
  XMLAttrib* dissipator_attrib = node->get_attrib("dissipator-id");
  if (dissipator_attrib)
  {
    const std::string dissipator_id = dissipator_attrib->get_string_value();
    _dissipator = dynamic_pointer_cast<Dissipation>(id_map[dissipator_id]);
    if (!_dissipator)
      throw std::runtime_error("Unable to load dissipator");
  }

  // read the collision detection plugin, if any
  XMLAttrib* coldet_plugin_attrib = node->get_attrib("collision-detection-plugin");
  if (coldet_plugin_attrib)
  {
    const std::string coldet_plugin_id = coldet_plugin_attrib->get_string_value();
    shared_ptr<CollisionDetection> coldet = dynamic_pointer_cast<CollisionDetection>(id_map[coldet_plugin_id]);
    if (!coldet)
      throw std::runtime_error("Unable to load collision detection plugin");
    _coldet = coldet;
    shared_ptr<ConstraintSimulator> shared_this = dynamic_pointer_cast<ConstraintSimulator>(shared_from_this());
    _coldet->set_simulator(shared_this);
  }

  // read the contact distance threshold, if any
  XMLAttrib* contact_dist_thresh_attrib = node->get_attrib("contact-dist-thresh");
  if (contact_dist_thresh_attrib)
    contact_dist_thresh = contact_dist_thresh_attrib->get_real_value();
  
  // read in any ContactParameters
  child_nodes = node->find_child_nodes("ContactParameters");
  if (!child_nodes.empty())
    contact_params.clear();
  for (list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    boost::shared_ptr<ContactParameters> cd(new ContactParameters);
    cd->load_from_xml(*i, id_map);
    contact_params[cd->objects] = cd;
  }

  // read all disabled pairs
  child_nodes = node->find_child_nodes("DisabledPair");
  for (std::list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    // get the two ID attributes
    XMLAttrib* id1_attrib = (*i)->get_attrib("object1-id");
    XMLAttrib* id2_attrib = (*i)->get_attrib("object2-id");

    // make sure that they were read
    if (!id1_attrib || !id2_attrib)
    {
      std::cerr << "ConstraintSimulator::load_from_xml() - did not find ";
      std::cerr << "object1-id and/or object2-id" << std::endl;
      std::cerr << "  in offending node: " << std::endl << *node;
      continue;
    }

    // get the two IDs
    const std::string& ID1 = id1_attrib->get_string_value();
    const std::string& ID2 = id2_attrib->get_string_value();

    // setup pairs of geometries to disable
    std::list<CollisionGeometryPtr> disabled1, disabled2;

    // find the first object
    if ((id_iter = id_map.find(ID1)) == id_map.end())
    {
      std::cerr << "ConstraintSimulator::load_from_xml() - could not find ";
      std::cerr << "object with object1-id" << std::endl;
      std::cerr << "  '" << ID1 << "' in offending node: " << std::endl << *node;
      continue;
    }
    BasePtr o1 = id_iter->second;
    CollisionGeometryPtr g1 = dynamic_pointer_cast<CollisionGeometry>(o1);
    if (g1)
      disabled1.push_back(g1);
    else
    {
      RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(o1);
      if (rb1)
        disabled1 = rb1->geometries;
      else
      {
        ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(o1);
        if (ab1)
        {
          BOOST_FOREACH(shared_ptr<RigidBodyd> rbd, ab1->get_links())
          {
            RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(rbd);
            disabled1.insert(disabled1.end(), rb->geometries.begin(), rb->geometries.end());
          }
        }
        else
        {
          std::cerr << "ConstraintSimulator::load_from_xml() - object with object1-id is not a usable type!" << std::endl;
          continue;
        }
      }
    }

    // find the second object
    if ((id_iter = id_map.find(ID2)) == id_map.end())
    {
      std::cerr << "ConstraintSimulator::load_from_xml() - could not find ";
      std::cerr << "object with object2-id" << std::endl;
      std::cerr << "  '" << ID2 << "' in offending node: " << std::endl << *node;
      continue;
    }
    BasePtr o2 = id_iter->second;
    CollisionGeometryPtr g2 = dynamic_pointer_cast<CollisionGeometry>(o2);
    if (g2)
      disabled2.push_back(g2);
    else
    {
      RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(o2);
      if (rb2)
        disabled2 = rb2->geometries;
      else
      {
        ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(o2);
        if (ab2)
        {
          BOOST_FOREACH(shared_ptr<RigidBodyd> rbd, ab2->get_links())
          {
            RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(rbd);
            disabled2.insert(disabled2.end(), rb->geometries.begin(), rb->geometries.end());
          }
        }
        else
        {
          std::cerr << "ConstraintSimulator::load_from_xml() - object with object2-id is not a usable type!" << std::endl;
          continue;
        }
      }
    }


   // add the pairs to the unchecked pairs list
   BOOST_FOREACH(CollisionGeometryPtr cg1, disabled1)
     BOOST_FOREACH(CollisionGeometryPtr cg2, disabled2)
       if (cg1 != cg2 && cg1->get_single_body() != cg2->get_single_body())
         unchecked_pairs.insert(make_sorted_pair(cg1, cg2));
  }
}

/// Implements Base::save_to_xml()
void ConstraintSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call Simulator's save method first
  Simulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "ConstraintSimulator";

  // save the dissipation mechanism
  if (_dissipator)
  {
    node->attribs.insert(XMLAttrib("dissipator-id", _dissipator->id));
    shared_objects.push_back(_dissipator);
  }

  // save any collision detection plugins
  if (!dynamic_pointer_cast<CCD>(_coldet))
  {
    node->attribs.insert(XMLAttrib("collision-detection-plugin", _coldet->id));
    shared_objects.push_back(_coldet);
  }

  // save the distance thresholds
  node->attribs.insert(XMLAttrib("contact-dist-thesh", contact_dist_thresh));

  // save all ContactParameters
  for (map<sorted_pair<BasePtr>, shared_ptr<ContactParameters> >::const_iterator i = contact_params.begin(); i != contact_params.end(); i++)
  {
    XMLTreePtr new_node(new XMLTree("ContactParameters"));
    node->add_child(new_node);
    i->second->save_to_xml(new_node, shared_objects);
  }

  // save all disabled pairs
  for (std::set<sorted_pair<CollisionGeometryPtr> >::const_iterator i = unchecked_pairs.begin(); i != unchecked_pairs.end(); i++)
  {
    XMLTreePtr child_node(new XMLTree("DisabledPair"));
    child_node->attribs.insert(XMLAttrib("object1-id", i->first->id));
    child_node->attribs.insert(XMLAttrib("object2-id", i->second->id));
    node->add_child(child_node);
  }
}


