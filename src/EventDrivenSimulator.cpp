/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <unistd.h>
#include <boost/tuple/tuple.hpp>
#include <Moby/XMLTree.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/DynamicBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/VariableStepIntegrator.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/SustainedUnilateralConstraintSolveFailException.h>
#include <Moby/InvalidStateException.h>
#include <Moby/InvalidVelocityException.h>
#include <Moby/EventDrivenSimulator.h>

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
EventDrivenSimulator::EventDrivenSimulator()
{
  constraint_callback_fn = NULL;
  constraint_post_callback_fn = NULL;
  post_mini_step_callback_fn = NULL;
  get_contact_parameters_callback_fn = NULL;
  render_contact_points = false;

  // setup the maximum constraint processing time
  max_constraint_time = std::numeric_limits<double>::max();

  // setup the minimum advancement
  min_advance = 1e-6;

  // setup the standard Euler step
  euler_step = 1e-3;

  // setup contact distance thresholds
  impacting_contact_dist_thresh = 1e-6;
  sustained_contact_dist_thresh = 1e-6;

  // setup absolute and relative error tolerances
  rel_err_tol = NEAR_ZERO;
  abs_err_tol = NEAR_ZERO;
  minimum_step = 1e-5;
}

/// Compares two unilateral constraints for purposes of mapping velocity tolerances
bool EventDrivenSimulator::UnilateralConstraintCmp::operator()(const UnilateralConstraint& e1, const UnilateralConstraint& e2)
 const
{
  if (e1.constraint_type == UnilateralConstraint::eLimit)
  {
    // we'll place limit constraints before contact constraints
    if (e2.constraint_type == UnilateralConstraint::eContact)
      return true;

    // if here, both are limit constraints
    unsigned lj1 = e1.limit_joint->get_coord_index() + e1.limit_dof;
    unsigned lj2 = e2.limit_joint->get_coord_index() + e2.limit_dof;
    if (lj1 < lj2)
      return true;
    else
    {
      assert(lj1 != lj2 || e1.limit_upper == e2.limit_upper);
      return false;
    }
  }
  else
  {
    // first constraint is contact; check to see whether the second constraint is a contact
    if (e2.constraint_type == UnilateralConstraint::eContact)
    {
      long cg11 = (long) e1.contact_geom1.get();
      long cg12 = (long) e1.contact_geom2.get();
      long cg21 = (long) e2.contact_geom1.get();
      long cg22 = (long) e2.contact_geom2.get();
      if (cg11+cg12 < cg21+cg22)
        return true;
      else
      {
        assert(cg11+cg12 != cg21+cg22 ||
               ((e1.contact_geom1 == e2.contact_geom1 &&
                 e1.contact_geom2 == e2.contact_geom2) ||
                (e1.contact_geom1 == e2.contact_geom2 &&
                 e1.contact_geom2 == e2.contact_geom1)));
        return false;
      }
    }
    else
      return false; // limits returned before contacts
  }
}

/// Handles sustained rigid unilateral constraints 
void EventDrivenSimulator::calc_rigid_sustained_unilateral_constraint_forces()
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

  // begin timing for constraint handling
  clock_t start = clock();

  // compute forces here...
  _rigid_unilateral_constraint_handler.process_constraints(_rigid_constraints);

  // tabulate times for constraint handling
  tms cstop;
  clock_t stop = times(&cstop);
  constraint_time += (double) (stop-start)/sysconf(_SC_CLK_TCK);

  // call the post-force application callback, if any
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_rigid_constraints, constraint_post_callback_data);
}

/// Computes the ODE for systems with sustained unilateral constraints
VectorNd& EventDrivenSimulator::ode_sustained_constraints(const VectorNd& x, double t, double dt, void* data, VectorNd& dx)
{
  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::ode_sustained_constraints(t=" << t << ") entered" << std::endl;

  // get the simulator
  shared_ptr<EventDrivenSimulator>& s = *((shared_ptr<EventDrivenSimulator>*) data);

  // see whether t=current time and the derivative has already been computed
  if (t == s->current_time && s->_current_accel_dx.size() > 0)
  {
    dx = s->_current_accel_dx;
    return dx;
  }

  // initialize the ODE index
  unsigned idx = 0;

  // resize dx
  dx.resize(x.size());

  // loop through all bodies, preparing to compute the ODE
  BOOST_FOREACH(DynamicBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // get the number of generalized coordinates and velocities
    const unsigned NGC = db->num_generalized_coordinates(DynamicBody::eEuler);
    const unsigned NGV = db->num_generalized_coordinates(DynamicBody::eSpatial);

    // get x for the body
    SharedConstVectorNd xsub = x.segment(idx, idx+NGC+NGV);

    FILE_LOG(LOG_SIMULATOR) << "evaluating derivative for body " << db->id << " at state " << xsub << std::endl;

    // compute the ODE
    db->prepare_to_calc_ode_sustained_constraints(xsub, t, dt, &db);

    // update idx
    idx += NGC+NGV;
  }

  // update the velocity bounds
  s->update_bounds();

  // check pairwise constraint violations
  double min_dist = s->check_pairwise_constraint_violations(t);

  // find unilateral constraints
  s->find_unilateral_constraints(s->sustained_contact_dist_thresh);
  if (s->_rigid_constraints.empty())
    FILE_LOG(LOG_SIMULATOR) << " *** constraints vector is unexpectedly empty! ***" << std::endl;

  // check velocity violations for constraints
  s->check_constraint_velocity_violations(t);

  // convert rigid constraints to acceleration constraints
  for (unsigned i=0; i< s->_rigid_constraints.size(); i++)
    s->_rigid_constraints[i].deriv_type = UnilateralConstraint::eAccel;

  // loop through all bodies, computing forward dynamics
  BOOST_FOREACH(DynamicBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    db->calc_fwd_dyn();
  }

  // compute compliant contact forces
  s->calc_compliant_unilateral_constraint_forces();

  // compute rigid unilateral constraint forces
  s->calc_rigid_sustained_unilateral_constraint_forces();

  // get super bodies for constraints
  vector<DynamicBodyPtr> bodies;
  BOOST_FOREACH(const UnilateralConstraint& e, s->_rigid_constraints)
    e.get_super_bodies(std::back_inserter(bodies));
  BOOST_FOREACH(const UnilateralConstraint& e, s->_compliant_constraints)
    e.get_super_bodies(std::back_inserter(bodies));
  std::sort(bodies.begin(), bodies.end());
  bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());

  // recompute forward dynamics for bodies in constraints
   BOOST_FOREACH(DynamicBodyPtr body, bodies)
     body->calc_fwd_dyn();

  // report accelerations
  if (LOGGING(LOG_CONSTRAINT))
  {
    FILE_LOG(LOG_CONSTRAINT) << " -- post-contact acceleration (all contacts): " << std::endl;
    BOOST_FOREACH(const UnilateralConstraint& e, s->_rigid_constraints)
      FILE_LOG(LOG_CONSTRAINT) << e;
  }

// TODO: remove this
static double last_t = -1.0;
static std::vector<double> last_vels; 
std::vector<double> this_vels(s->_rigid_constraints.size());
for (unsigned i=0; i< s->_rigid_constraints.size(); i++)
  this_vels[i] = s->_rigid_constraints[i].calc_constraint_vel();
if (last_vels.size() == this_vels.size())
{
  double h = t - last_t;
  for (unsigned i=0; i< this_vels.size(); i++)
  {
    FILE_LOG(LOG_CONSTRAINT) << "Velocity at " << last_t << ": " << last_vels[i] << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "Velocity at " << t << ": " << this_vels[i] << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "Numerically computed acceleration: " << (this_vels[i] - last_vels[i])/h << std::endl;
  }
}
last_t = t;
last_vels = this_vels;

  // reset idx
  idx = 0;

  // loop through all bodies, computing the ODE
  BOOST_FOREACH(DynamicBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // get the number of generalized coordinates and velocities
    const unsigned NGC = db->num_generalized_coordinates(DynamicBody::eEuler);
    const unsigned NGV = db->num_generalized_coordinates(DynamicBody::eSpatial);

    // get dx for the body
    SharedVectorNd dxsub = dx.segment(idx, idx+NGC+NGV);

    // compute the ODE
    db->ode(t, dt, &db, dxsub);

    FILE_LOG(LOG_SIMULATOR) << " ODE evaluation for body " << db->id << ": " << dxsub << std::endl;

    // update idx
    idx += NGC+NGV;
  }

  // see whether to set current time derivative
  if (t == s->current_time)
    s->_current_accel_dx = dx;

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::ode_sustained_constraints(t=" << t << ") exited" << std::endl;

  // return the ODE
  return dx;
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
shared_ptr<ContactParameters> EventDrivenSimulator::get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const
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
  SingleBodyPtr singlebody1 = geom1->get_single_body();
  SingleBodyPtr singlebody2 = geom2->get_single_body();
  BasePtr sb1 = singlebody1;
  BasePtr sb2 = singlebody2;

  // search for contact geometry 1 and rigid body 2
  if ((iter = contact_params.find(make_sorted_pair(g1, sb2))) != contact_params.end())
    return iter->second;

  // search for contact geometry 2 and rigid body 1
  if ((iter = contact_params.find(make_sorted_pair(g2, sb1))) != contact_params.end())
    return iter->second;

  // search for both rigid bodies
  if ((iter = contact_params.find(make_sorted_pair(sb1, sb2))) != contact_params.end())
    return iter->second;

  // get the articulated bodies, if any
  RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(singlebody1);
  RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(singlebody2);
  BasePtr ab1, ab2;
  if (rb1)
    ab1 = rb1->get_articulated_body();
  if (rb2)
    ab2 = rb2->get_articulated_body();

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

  // still here?  no contact data found
  return shared_ptr<ContactParameters>();
}

/// Draws a ray directed from a contact point along the contact normal
void EventDrivenSimulator::visualize_contact( UnilateralConstraint& constraint ) 
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
void EventDrivenSimulator::calc_impacting_unilateral_constraint_forces(double dt)
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

  // if the setting is enabled, draw all contact constraints
  if( render_contact_points ) {
    for ( std::vector<UnilateralConstraint>::iterator it = _rigid_constraints.begin(); it < _rigid_constraints.end(); it++ ) {
      UnilateralConstraint& constraint = *it;
      if( constraint.constraint_type != UnilateralConstraint::eContact ) continue;
      visualize_contact( constraint );
    }
  }

  // begin timing for constraint handling
  clock_t start = clock();

  // compute impulses here...
  try
  {
    if (dt <= 0.0)
      _impact_constraint_handler.process_constraints(_rigid_constraints, max_constraint_time);
    else
      _impact_constraint_handler.process_constraints(_rigid_constraints, max_constraint_time, 0.1/dt);
  }
  catch (ImpactToleranceException e)
  {
    std::cerr << "warning: impacting constraint tolerances exceeded" << std::endl;
  }

  // tabulate times for constraint handling
  tms cstop;
  clock_t stop = times(&cstop);
  constraint_time += (double) (stop-start)/sysconf(_SC_CLK_TCK);

  // call the post application callback, if any
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_rigid_constraints, constraint_post_callback_data);
}

/// Computes compliant contact forces 
void EventDrivenSimulator::calc_compliant_unilateral_constraint_forces()
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
void EventDrivenSimulator::preprocess_constraint(UnilateralConstraint& e)
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
    SingleBodyPtr sb1(e.contact_geom1->get_single_body());
    SingleBodyPtr sb2(e.contact_geom2->get_single_body());
    std::cerr << "EventDrivenSimulator::preprocess_constraint() warning- no contact ";
    std::cerr << "data for contact" << std::endl;
    std::cerr << "  between " << e.contact_geom1->id << " (body ";
    std::cerr << sb1->id << ") and " << e.contact_geom2->id;
    std::cerr << " (body " << sb2->id << ")" << std::endl;
    std::cerr << "  ... ignoring" << std::endl;
  }
}

/// Sets up the list of collision geometries
void EventDrivenSimulator::determine_geometries()
{
  // clear the list at first
  _geometries.clear();

  // determine all geometries
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
    if (rb)
      _geometries.insert(_geometries.end(), rb->geometries.begin(), rb->geometries.end());
    else
    {
      ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
      BOOST_FOREACH(RigidBodyPtr rb, ab->get_links())
        _geometries.insert(_geometries.end(), rb->geometries.begin(), rb->geometries.end());
    }
  }

  // sort and remove duplicates
  std::sort(_geometries.begin(), _geometries.end());
  _geometries.erase(std::unique(_geometries.begin(), _geometries.end()), _geometries.end());
}

/// Steps the simulator forward by the given step size
double EventDrivenSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

  // clear timings
  dynamics_time = (double) 0.0;
  constraint_time = (double) 0.0;
  coldet_time = (double) 0.0;

  // clear statistics and step times
  std::fill(step_times, step_times+6, 0.0);
  std::fill(step_stats, step_stats+6, 0);
  int_min_step_stat = INF;
  int_max_step_stat = 0.0;
  int_mean_step_stat = 0.0;

  // setup timer
  clock_t start = clock();

  // determine the set of collision geometries
  determine_geometries();

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif

  // setup the time stepped
  double h = 0.0;

  // step until the requisite time has elapsed
  while (h < step_size)
  {
    // clear stored derivatives
    _current_dx.resize(0);
    _current_accel_dx.resize(0);

    FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << " by " << (step_size - h) << std::endl;
    if (LOGGING(LOG_SIMULATOR))
    {
      VectorNd q;
      BOOST_FOREACH(DynamicBodyPtr db, _bodies)
      {
        db->get_generalized_coordinates(DynamicBody::eEuler, q);
        FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " coordinates (before): " << q << std::endl;
      }
    }

    // start with initial estimates
    reset_limit_estimates();

    // get amount remaining to step
    double dt = step_size - h;

    // do broad phase collision detection (must be done before any Euler steps)
    broad_phase(dt);

    // compute pairwise distances at the current configuration
    calc_pairwise_distances();

    // if dt is sufficiently small, do Euler step
    if (dt <= euler_step)
    {
      FILE_LOG(LOG_SIMULATOR) << "  - step size really small; doing semi-implicit step" << std::endl;

      // do the Euler step
      step_si_Euler(dt);
      h += dt;

      // finish timing and record the stat
      step_times[0] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      step_stats[0]++;

       // call the mini-callback
      if (post_mini_step_callback_fn)
        post_mini_step_callback_fn(this);

      // break out of the while loop (could just call continue)
      break;
    }

    // update constraint violations
    update_constraint_violations(_pairwise_distances);

    // if there are any impacts at the current time, handle them
    FILE_LOG(LOG_SIMULATOR) << "  - preparing to handle any impacts at the current time" << std::endl;
    find_unilateral_constraints(impacting_contact_dist_thresh);
    calc_impacting_unilateral_constraint_forces(-1.0);

restart_with_new_limits:

    // mark velocity estimates as valid
    validate_limit_estimates();

    // determine the maximum step according to conservative advancement
    double safe_dt = std::min(calc_CA_step(), dt);
    FILE_LOG(LOG_SIMULATOR) << "  - conservative advancement step: " << safe_dt << std::endl;

    // get next possible constraint time
    double next_event_time = calc_next_CA_step(sustained_contact_dist_thresh);
    FILE_LOG(LOG_SIMULATOR) << "  - *next* conservative advancement step: " << next_event_time << std::endl;

    // we know that the current time is safe, so if the conservative
    // advancement step is zero, set it to the next event time
    if (safe_dt == 0.0)
      safe_dt = std::min(dt, next_event_time);

    // if (the distance between bodies is small (next_event_time will be
    // < INF if there is an event at the current time) and the next time of
    // contact is greater than an Euler step) or (safe_dt is greater than the
    // Euler step), we can *try* generic integration
    if ((next_event_time < INF && next_event_time > euler_step) ||
        safe_dt > euler_step)
    {
      // try the integration
      IntegrationResult istat = integrate_generic(safe_dt, start);

      // if limits were exceeded during integration, compute the conservative
      // advancement step and redo
      if (istat == eVelocityLimitExceeded)
        goto restart_with_new_limits;
      else if (istat == eMinStepReached)
      {
        // do Euler step; if we're here, we know euler_step < dt
        step_si_Euler(euler_step);
        h += euler_step;

        // finish timing and record the stat
        step_times[0] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
        step_stats[0]++;
        start = clock();
      }
      else
      {
        assert(istat == eIntegrationSuccessful);
        h += safe_dt;
      }
    }
    else
    {
      // do Euler step
      step_si_Euler(dt);
      h += dt;

      // finish timing, record the stat, and restart the clock
      step_times[0] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      step_stats[0]++;
      start = clock();
    }

    // call the mini-callback
    if (post_mini_step_callback_fn)
      post_mini_step_callback_fn(this);
  }

  // call the callback
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  return step_size;
}

/// Finds the next event time assuming constant velocity
/**
 * This method returns the next possible time of contact, discarding current
 * contacts from consideration. 
 * \note proper operation of this function is critical. If the function
 *       improperly designates an event as not occuring at the current time,
 *       calc_next_CA_Euler_step(.) will return a small value and prevent large
 *       integration steps from being taken. If the function improperly
 *       designates an event as occuring at the current time, constraint
 *       violation could occur.
 */
double EventDrivenSimulator::calc_next_CA_Euler_step(double contact_dist_thresh) const
{
  const double INF = std::numeric_limits<double>::max();
  double next_event_time = INF;

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::calc_next_CA_Euler_step entered" << std::endl; 

  // process each articulated body, looking for next joint events
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // see whether the i'th body is articulated
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(_bodies[i]);
    if (!ab)
      continue;

    // if the body is kinematically controlled, do nothing
    if (ab->get_kinematic())
      continue;

    // get limit events in [t, t+dt] (if any)
    const vector<JointPtr>& joints = ab->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
      for (unsigned j=0; j< joints[i]->num_dof(); j++)
      {
        if (joints[i]->q[j] < joints[i]->hilimit[j] && joints[i]->qd[j] > 0.0)
        {
          double t = (joints[i]->hilimit[j] - joints[i]->q[j])/joints[i]->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
        if (joints[i]->q[j] > joints[i]->lolimit[j] && joints[i]->qd[j] < 0.0)
        {
          double t = (joints[i]->lolimit[j] - joints[i]->q[j])/joints[i]->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
      }
  }

  // if the distance between any pair of bodies is sufficiently small
  // get next possible event time
  BOOST_FOREACH(const PairwiseDistInfo& pdi, _pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

    // TODO: this code needs to be extended such that all points on a rigid
    //       body that are within a single topological "hop" from the contact 
    //       point are evaluated for their next TOC with the other body
    //       We do not presently consider further contacts between a pair if
    //       there exists a contact at the current configuration (which means
    //       that contacts between polyhedral shapes will not be handled
    //       correctly.) 

    // if the distance is below the threshold, we have found a current event
    // (we want to skip current events)
    if (pdi.dist > contact_dist_thresh)
    {
      // compute an upper bound on the event time
      double event_time = _ccd.calc_CA_Euler_step(pdi);

      FILE_LOG(LOG_SIMULATOR) << "Next contact time between " << pdi.a->get_single_body()->id << " and " << pdi.b->get_single_body()->id << ": " << event_time << std::endl;

      // not a current event, find when it could become active
      next_event_time = std::min(next_event_time, event_time);
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::calc_next_CA_Euler_step exited" << std::endl; 

  return next_event_time;
}


/// Finds the next event time
/**
 * This method assumes an event is occurring at the current time. If an event
 * is not occuring at the current time, this method will return INF.
 * \note proper operation of this function is critical. If the function
 *       improperly designates an event as not occuring at the current time,
 *       calc_next_CA_step(.) will return a small value and prevent large
 *       integration steps from being taken. If the function improperly
 *       designates an event as occuring at the current time, constraint
 *       violation could occur.
 */
double EventDrivenSimulator::calc_next_CA_step(double contact_dist_thresh) const
{
  const double INF = std::numeric_limits<double>::max();
  bool found_one = false;
  double next_event_time = INF;

  // process each articulated body, looking for next joint events
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // see whether the i'th body is articulated
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(_bodies[i]);
    if (!ab)
      continue;

    // if the body is kinematically controlled, do nothing
    if (ab->get_kinematic())
      continue;

    // TODO: we need to modify this so we can see whether any joint limit
    // events are active. We should be able to remove
    // find_next_joint_limit_time(.) and add a CA step for each joint

    // get limit events in [t, t+dt] (if any)
    next_event_time = std::min(next_event_time, ab->find_next_joint_limit_time());
  }

  // if the distance between any pair of bodies is sufficiently small
  // get next possible event time
  BOOST_FOREACH(const PairwiseDistInfo& pdi, _pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

    // if the distance is below the threshold, we have found a current event
    if (pdi.dist < contact_dist_thresh)
      found_one = true;
    else
      // not a current event, find when it could become active
      next_event_time = std::min(next_event_time, _ccd.calc_CA_step(pdi));
  }

  if (!found_one)
    return INF;
  else
    return next_event_time;
}

/// Attempts to do generic integration
EventDrivenSimulator::IntegrationResult EventDrivenSimulator::integrate_generic(double dt, clock_t& start)
{
  do
  {
    // save the state of the system
    save_state();

    try
    {
      // do "smart" integration (watching for state violation)
      integrate_with_sustained_constraints(dt);

      FILE_LOG(LOG_SIMULATOR) << "Integration with sustained constraints successful" << std::endl;

      // check whether velocity estimates have been exceeded
      BOOST_FOREACH(DynamicBodyPtr db, _bodies)
      {
        if (db->limit_estimates_exceeded())
        {
          FILE_LOG(LOG_SIMULATOR) << " ** limit estimates exceeded for body " << db->id << "; retrying with new estimates" << std::endl;

          // reset the state of all bodies
          restore_state();

          // update the stats and the clock
          step_stats[4]++;
          step_times[4] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
          start = clock();

          return eVelocityLimitExceeded;
        }
      }

      // update the time
      current_time += dt;

      // update the stats and the clock
      step_stats[5]++;
      step_times[5] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      start = clock();

      return eIntegrationSuccessful;
    }
    catch (InvalidStateException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** attempted to evaluate derivative at invalid state; halving acceleration step size to " << (dt*0.5) << std::endl;

      // update the stats and the clock
      step_stats[1]++;
      step_times[1] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      start = clock();
    }
    catch (InvalidVelocityException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** attempted to evaluate derivative at invalid velocity; halfing acceleration step size to " << (dt*0.5) << std::endl;

      // update the stats and the clock
      step_stats[2]++;
      step_times[2] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      start = clock();
    }
    catch (SustainedUnilateralConstraintSolveFailException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** failed to solve an LCP; halving step size" << std::endl;

      // update the stats and the clock
      step_stats[3]++;
      step_times[3] += (double) (clock() - start)/sysconf(_SC_CLK_TCK);
      start = clock();
    }

    // ********************************************************
    // if we're still here, we caught an exception
    // ********************************************************

    // restore the state of the system (generalized coords/velocities)
    restore_state();

    // half the step size
    dt *= 0.5;
  }
  while (dt > euler_step);

  // if we're still here, step size has become too small
  return eMinStepReached;
}

/// Validates limit estimates
void EventDrivenSimulator::validate_limit_estimates()
{
  BOOST_FOREACH(DynamicBodyPtr body, _bodies)
    body->validate_limit_estimates();
}

/// Saves the state of the system (all dynamic bodies) at the current time
void EventDrivenSimulator::save_state()
{
  // resize the vector if necessary
  _qsave.resize(_bodies.size());
  _qdsave.resize(_bodies.size());

  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, _qsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, _qdsave[i]);
  }
}

/// Restores the state of the 'system' (all dynamic bodies)
void EventDrivenSimulator::restore_state()
{
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _qsave[i]);
    _bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, _qdsave[i]);
  }
}

/// Computes pairwise distances of geometries at their current poses, using broad phase results to determine which pairs should be checked
/**
 * \param pairwise_distances on return, contains the pairwise distances
 */
void EventDrivenSimulator::calc_pairwise_distances()
{
  // clear the vector
  _pairwise_distances.clear();

  for (unsigned i=0; i< _pairs_to_check.size(); i++)
  {
    PairwiseDistInfo pdi;
    pdi.a = _pairs_to_check[i].first;
    pdi.b = _pairs_to_check[i].second;
    FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::calc_pairwise_distances() - computing signed distance between " << pdi.a->get_single_body()->id << " and " << pdi.b->get_single_body()->id << std::endl;
    pdi.dist = CollisionGeometry::calc_signed_dist(pdi.a, pdi.b, pdi.pa, pdi.pb);
    _pairwise_distances.push_back(pdi);
  }
}

/// Does broad phase collision detection, identifying which pairs of geometries may come into contact over time step of dt
void EventDrivenSimulator::broad_phase(double dt)
{
  // call the broad phase
  _ccd.broad_phase(dt, _bodies, _pairs_to_check);

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

/// Checks whether bodies violate contact constraint velocity tolerances
void EventDrivenSimulator::check_constraint_velocity_violations(double t)
{
  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::check_constraint_velocity_violations() entered" << std::endl;

  // loop over all constraints
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
  {
    // get the constraint velocity
    double ev = _rigid_constraints[i].calc_constraint_vel();

    // look for the constraint in the mapping
    std::map<UnilateralConstraint, double, UnilateralConstraintCmp>::const_iterator zv_tol = _zero_velocity_tolerances.find(_rigid_constraints[i]);
    if (zv_tol == _zero_velocity_tolerances.end())
    {
      _zero_velocity_tolerances[_rigid_constraints[i]] = NEAR_ZERO;
      zv_tol = _zero_velocity_tolerances.find(_rigid_constraints[i]);
    }

    FILE_LOG(LOG_SIMULATOR) << " -- constraint velocity: " << ev << std::endl;

    // check whether it is larger than allowed
    if (ev < -zv_tol->second - NEAR_ZERO)
    {
      FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::check_constraint_velocity_violations() about to throw exception..." << std::endl;
//      throw InvalidVelocityException(t);
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::check_constraint_velocity_violations() exiting" << std::endl;
}

/// Checks whether bodies violate interpenetration constraints
double EventDrivenSimulator::check_pairwise_constraint_violations(double t)
{
  double min_dist = std::numeric_limits<double>::max();

  // update constraint violation due to increasing interpenetration
  // loop over all pairs of geometries
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
  {
    // get the two collision geometries
    CollisionGeometryPtr cg1 = _pairwise_distances[i].a;
    CollisionGeometryPtr cg2 = _pairwise_distances[i].b;

    // get the two rigid bodies
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(cg1->get_single_body());
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(cg2->get_single_body());

    // only process if neither of the bodies is compliant
    if (rb1->compliance == RigidBody::eCompliant || 
        rb2->compliance == RigidBody::eCompliant)
      continue; 

    // compute the distance between the two bodies
    Point3d p1, p2;
    double d = CollisionGeometry::calc_signed_dist(cg1, cg2, p1, p2);
    if (d <= _ip_tolerances[make_sorted_pair(cg1, cg2)] - NEAR_ZERO)
    {
      FILE_LOG(LOG_SIMULATOR) << "Interpenetration detected between " << cg1->get_single_body()->id << " and " << cg2->get_single_body()->id << ": " << d << std::endl;
      throw InvalidStateException();
    }

    // update the minimum distance
    min_dist = std::min(d, min_dist);
  }

  return std::max(0.0, min_dist);
}

/// Updates constraint violation at beginning of integration step
void EventDrivenSimulator::update_constraint_violations(const vector<PairwiseDistInfo>& pairwise_distances)
{
  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::update_constraint_violations() entered" << std::endl;

  // set possible constraint violation
  BOOST_FOREACH(const PairwiseDistInfo& pdi, pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

    if (pdi.dist <= 0)
      _ip_tolerances[make_sorted_pair(pdi.a, pdi.b)] = pdi.dist;
    else
      _ip_tolerances[make_sorted_pair(pdi.a, pdi.b)] = 0.0;
  }

  // update joint constraint interpenetration
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
    if (ab)
      ab->update_joint_constraint_violations();
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::update_constraint_violations() exited" << std::endl;
}

/// Computes a conservative advancement step
double EventDrivenSimulator::calc_CA_step()
{
  // setup safe amount to step
  double dt = std::numeric_limits<double>::max();

  // do joint limit CA step first (it's faster)
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    // try to get it as an articulated body
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
    if (!ab)
      continue;

    // compute best dt
    dt = std::min(dt, ab->calc_CA_time_for_joints());
    if (dt <= 0.0)
      return dt;
  }

  // do narrow-phase collision detection here
  BOOST_FOREACH(const PairwiseDistInfo& pdi, _pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

    // call conservative advancement method
    double step = _ccd.calc_CA_step(pdi);
    dt = std::min(dt, step);
    if (dt <= 0.0)
      return dt;
  }

  return dt;
}

void EventDrivenSimulator::reset_limit_estimates() const
{
  // now compute the bounds
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    // first, reset the limit estimates
    db->reset_limit_estimates();
  }
}

/// Computes forward dynamics for all bodies
void EventDrivenSimulator::calc_fwd_dyn()
{
  // clear force accumulators, then add all recurrent and compliant
  // constraint forces
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    // clear the force accumulators on the body
    db->reset_accumulators();

    // add all recurrent forces on the body
    const list<RecurrentForcePtr>& rfs = db->get_recurrent_forces();
    BOOST_FOREACH(RecurrentForcePtr rf, rfs)
      rf->add_force(db);
    
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
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    // calculate forward dynamics at state x
    db->calc_fwd_dyn();
  }
}

/// Integrates bodies' velocities forward by dt using Euler integration
void EventDrivenSimulator::integrate_velocities_Euler(double dt)
{
  VectorNd qd, qdd;

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::integrate_velocities_Euler() entered " << std::endl;

  // if forward dynamics are already computed, get the accelerations
  if (_current_dx.size() > 0)
  {
    // setup a coordinate index
    unsigned idx = 0;

    // loop through all bodies, computing the ODE
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      if (db->get_kinematic())
        continue;

      // get the number of generalized coordinates and velocities
      const unsigned NGC = db->num_generalized_coordinates(DynamicBody::eEuler);
      const unsigned NGV = db->num_generalized_coordinates(DynamicBody::eSpatial);

      // get the acceleration for the body and multiply by dt
      qdd = _current_dx.segment(idx+NGC, idx+NGC+NGV);
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " acceleration: " << qdd << std::endl;
      qdd *= dt;

      // update the generalized velocity
      db->get_generalized_velocity(DynamicBody::eSpatial, qd);
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " velocity: " << qd << std::endl;
      qd += qdd;
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " new velocity: " << qd << std::endl;
      db->set_generalized_velocity(DynamicBody::eSpatial, qd);

      // update idx
      idx += NGC+NGV;
    }
  }
  else
  {
    // first compute forward dynamics for all bodies
    calc_fwd_dyn();

    // now update all velocities
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      // get the generalized acceleration
      db->get_generalized_acceleration(qdd);
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " acceleration: " << qdd << std::endl;
      qdd *= dt;

      // update the generalized velocity
      db->get_generalized_velocity(DynamicBody::eSpatial, qd);
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " velocity: " << qd << std::endl;
      qd += qdd;
      FILE_LOG(LOG_SIMULATOR) << "body " << db->id << " new velocity: " << qd << std::endl;
      db->set_generalized_velocity(DynamicBody::eSpatial, qd);
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::integrate_velocities_Euler() exited " << std::endl;
}

/// Integrates bodies' positions forward by dt using Euler integration
void EventDrivenSimulator::integrate_positions_Euler(double dt)
{
  VectorNd q, qd;

  // update all positions
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    db->get_generalized_velocity(DynamicBody::eEuler, qd);
    qd *= dt;
    db->get_generalized_coordinates(DynamicBody::eEuler, q);
    q += qd;
    db->set_generalized_coordinates(DynamicBody::eEuler, q);
  }
}

/// Finds the set of unilateral constraints
void EventDrivenSimulator::find_unilateral_constraints(double contact_dist_thresh)
{
  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::find_unilateral_constraints() entered" << std::endl;

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
        _ccd.find_contacts(pdi.a, pdi.b, std::back_inserter(_compliant_constraints));
      else        
        _ccd.find_contacts(pdi.a, pdi.b, std::back_inserter(_rigid_constraints));
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

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::find_unilateral_constraints() exited" << std::endl;
}

/// Does a semi-implicit step
void EventDrivenSimulator::step_si_Euler(double dt)
{
  FILE_LOG(LOG_SIMULATOR) << "-- doing semi-implicit Euler step" << std::endl;

  // integrate bodies' velocities forward by dt
  integrate_velocities_Euler(dt);
  FILE_LOG(LOG_SIMULATOR) << "   integrating velocities forward by " << dt << std::endl;

  // setup target time
  double target_time = current_time + dt;

  // keep looping until we break out
  while (true)
  {
    // determine constraints (contacts, limits) that are currently active
    FILE_LOG(LOG_SIMULATOR) << "   finding constraints" << std::endl;
    find_unilateral_constraints(impacting_contact_dist_thresh);

    // solve constraints to yield new velocities
    FILE_LOG(LOG_SIMULATOR) << "   handling constraints" << std::endl;
    calc_impacting_unilateral_constraint_forces(-1.0);

    if (LOGGING(LOG_SIMULATOR))
    {
      VectorNd qd;
      BOOST_FOREACH(DynamicBodyPtr db, _bodies)
      {
        db->get_generalized_velocity(DynamicBody::eSpatial, qd);
        FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " velocity (after constraint treatment): " << qd << std::endl;
      }
    }

    // get the time of the next event(s)
    double h = std::min(calc_next_CA_Euler_step(impacting_contact_dist_thresh), target_time - current_time);
    FILE_LOG(LOG_SIMULATOR) << "   position integration: " << h << std::endl;

    // integrate bodies' positions forward by that time using new velocities
    integrate_positions_Euler(h);
    if (LOGGING(LOG_SIMULATOR))
    {
      VectorNd q;
      BOOST_FOREACH(DynamicBodyPtr db, _bodies)
      {
        db->get_generalized_coordinates(DynamicBody::eEuler, q);
        FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " position (after integration): " << q << std::endl;
      }
    }

    // update s and the current time
    current_time += h;

    // see whether to update pairwise distances
    if (current_time < target_time)
    {
      broad_phase(target_time - current_time);
      calc_pairwise_distances();
    }
    else
      break;
  }

  FILE_LOG(LOG_SIMULATOR) << "-- semi-implicit Euler step completed" << std::endl;
}

/// Implements Base::load_from_xml()
void EventDrivenSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify node name b/c this is abstract class
  assert(strcasecmp(node->name.c_str(), "EventDrivenSimulator") == 0);

  // first, load all data specified to the Simulator object
  Simulator::load_from_xml(node, id_map);

  // read the maximum time to process constraints, if any
  XMLAttrib* max_constraint_time_attrib = node->get_attrib("max-constraint-time");
  if (max_constraint_time_attrib)
    max_constraint_time = max_constraint_time_attrib->get_real_value();

  // read the maximum Euler step
  XMLAttrib* Euler_step_attrib = node->get_attrib("Euler-step");
  if (Euler_step_attrib)
    euler_step = Euler_step_attrib->get_real_value();

  // read the minimum advancement
  XMLAttrib* min_advance_attrib = node->get_attrib("min-advance");
  if (min_advance_attrib)
    min_advance = min_advance_attrib->get_real_value();

  // read the error tolerances
  XMLAttrib* rel_tol_attrib = node->get_attrib("rel-err-tol");
  XMLAttrib* abs_tol_attrib = node->get_attrib("abs-err-tol");
  if (rel_tol_attrib)
    rel_err_tol = rel_tol_attrib->get_real_value();
  if (abs_tol_attrib)
    abs_err_tol = abs_tol_attrib->get_real_value();

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
      std::cerr << "EventDrivenSimulator::load_from_xml() - did not find ";
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
      std::cerr << "EventDrivenSimulator::load_from_xml() - could not find ";
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
          BOOST_FOREACH(RigidBodyPtr rb, ab1->get_links())
            disabled1.insert(disabled1.end(), rb->geometries.begin(), rb->geometries.end());
        }
        else
        {
          std::cerr << "EventDrivenSimulator::load_from_xml() - object with object1-id is not a usable type!" << std::endl;
          continue;
        }
      }
    }

    // find the second object
    if ((id_iter = id_map.find(ID2)) == id_map.end())
    {
      std::cerr << "EventDrivenSimulator::load_from_xml() - could not find ";
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
          BOOST_FOREACH(RigidBodyPtr rb, ab2->get_links())
            disabled2.insert(disabled2.end(), rb->geometries.begin(), rb->geometries.end());
        }
        else
        {
          std::cerr << "EventDrivenSimulator::load_from_xml() - object with object2-id is not a usable type!" << std::endl;
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
void EventDrivenSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call Simulator's save method first
  Simulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "EventDrivenSimulator";

  // save the maximum constraint time
  node->attribs.insert(XMLAttrib("max-constraint-time", max_constraint_time));

  // save the maximum Euler step
  node->attribs.insert(XMLAttrib("Euler-step", euler_step));

  // save the minimum advancement step
  node->attribs.insert(XMLAttrib("min-advance", min_advance));

  // save the error tolerances
  node->attribs.insert(XMLAttrib("rel-err-tol", rel_err_tol));
  node->attribs.insert(XMLAttrib("abs-err-tol", abs_err_tol));

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


