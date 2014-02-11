/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

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
#include <Moby/InvalidStateException.h>
#include <Moby/EventDrivenSimulator.h>

#ifdef USE_OSG
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/PositionAttitudeTransform>
#include <osg/Quat>
#endif // USE_OSG

using std::endl;
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
  event_callback_fn = NULL;
  event_post_impulse_callback_fn = NULL;
  post_mini_step_callback_fn = NULL;
  get_contact_parameters_callback_fn = NULL;
  _simulation_violated = false;
  render_contact_points = false;

  // setup the maximum event processing time
  max_event_time = std::numeric_limits<double>::max();

  // setup the maximum Euler step
  max_Euler_step = 1e-3;

  // setup absolute and relative error tolerances
  rel_err_tol = NEAR_ZERO;
  abs_err_tol = NEAR_ZERO;
  minimum_step = 1e-5;

  // clear timings
  dynamics_time = (double) 0.0;
  event_time = (double) 0.0;
  coldet_time = (double) 0.0;
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
void EventDrivenSimulator::visualize_contact( Event& event ) {

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

  // the osg node this event visualization will attach to 
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
  Vector3d axis = Vector3d::cross( event.contact_normal, z );
  if( axis.norm_inf() < NEAR_ZERO) {
    // z and normal are parallel, axis ill defined
    if( event.contact_normal[2] > 0 ) {
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
    theta = -std::acos( Vector3d::dot( event.contact_normal, z ) );
    // Note: theta calculation works but is not robust, could be rotated in opposite direction
  }
  osg::Quat q = osg::Quat( axis[0]*std::sin(theta/2), axis[1]*std::sin(theta/2), axis[2]*std::sin(theta/2), std::cos(theta/2) );

  // create the visualization transform
  osg::PositionAttitudeTransform* contact_transform = new osg::PositionAttitudeTransform();
  contact_transform->setPosition( osg::Vec3( event.contact_point[0], event.contact_point[1], event.contact_point[2] ) );
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
  validator_transform->setPosition( osg::Vec3( event.contact_point[0], event.contact_point[1], event.contact_point[2] ) );
  validator_transform->setScale( osg::Vec3( validator_scale, validator_scale, validator_scale ) );
  validator_root->addChild( validator_transform );

  // validator geometry
  osg::Sphere* validator_geometry = new osg::Sphere( osg::Vec3( 0, 0, 0 ), 1.0 );
  osg::ShapeDrawable* validator_shape = new osg::ShapeDrawable( validator_geometry, hints );
  validator_shape->setColor( color );

  // validator transform follows the normal out to a distance of validator_ray_length
  // Note: the validator is not rotated at all.  It is translated from the point along the normal
  osg::PositionAttitudeTransform* validator_end_transform = new osg::PositionAttitudeTransform();
  validator_end_transform->setPosition( osg::Vec3( event.contact_normal[0] * validator_ray_length, event.contact_normal[1] * validator_ray_length, event.contact_normal[2] * validator_ray_length ) );
  validator_transform->addChild( validator_end_transform );

  // add all validator constituents to the group
  osg::Geode* validator_geode = new osg::Geode();
  validator_transform->addChild( validator_end_transform );
  validator_end_transform->addChild( validator_geode );
  validator_geode->addDrawable( validator_shape );
  add_transient_vdata( validator_root );

  #endif // USE_OSG
}

/// Handles events
void EventDrivenSimulator::handle_events()
{
  // call the callback function, if any
  if (event_callback_fn)
    (*event_callback_fn)(_events, event_callback_data);

  // preprocess events
  for (unsigned i=0; i< _events.size(); i++)
    preprocess_event(_events[i]);

  // if the setting is enabled, draw all contact events
  if( render_contact_points ) {
    for ( std::vector<Event>::iterator it = _events.begin(); it < _events.end(); it++ ) {
      Event event = *it;
      if( event.event_type != Event::eContact ) continue;
      visualize_contact( event );
    }
  }

  // begin timing for event handling 
  tms cstart;  
  clock_t start = times(&cstart);

  // compute impulses here...
  try
  {
    _impact_event_handler.process_events(_events, max_event_time);
  }
  catch (ImpactToleranceException e)
  {
    // process events, updating tolerances
    BOOST_FOREACH(Event* ev, e.events)
    {
      double event_v = ev->calc_event_vel();
      _event_tolerances[*ev] = std::fabs(event_v) + std::numeric_limits<double>::epsilon();  
    }
  }

  // tabulate times for event handling 
  tms cstop;  
  clock_t stop = times(&cstop);
  event_time += (double) (stop-start)/CLOCKS_PER_SEC;

  // call the post-impulse application callback, if any 
  if (event_post_impulse_callback_fn)
    (*event_post_impulse_callback_fn)(_events, event_post_impulse_callback_data);
}

/// Performs necessary preprocessing on an event
void EventDrivenSimulator::preprocess_event(Event& e) 
{
  // no pre-processing for limit events currently...
  if (e.event_type == Event::eLimit)
    return;

  // no pre-processing for (none) events
  if (e.event_type == Event::eNone)
    return;

  // get the contact parameters 
  assert(e.event_type == Event::eContact);
  shared_ptr<ContactParameters> cparams = get_contact_parameters(e.contact_geom1, e.contact_geom2);
  if (cparams)
    e.set_contact_parameters(*cparams);
  else
  {
    SingleBodyPtr sb1(e.contact_geom1->get_single_body());
    SingleBodyPtr sb2(e.contact_geom2->get_single_body());
    std::cerr << "EventDrivenSimulator::preprocess_event() warning- no contact ";
    std::cerr << "data for contact" << std::endl;
    std::cerr << "  between " << e.contact_geom1->id << " (body ";
    std::cerr << sb1->id << ") and " << e.contact_geom2->id;
    std::cerr << " (body " << sb2->id << ")" << std::endl;
    std::cerr << "  ... ignoring" << std::endl;
  }
}

/// Saves the coords of all bodies
void EventDrivenSimulator::get_coords(vector<VectorNd>& q) const
{
  // resize the vector if necessary
  q.resize(_bodies.size());

  for (unsigned i=0; i< _bodies.size(); i++)
    _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q[i]);
}

/// Saves the velocities of all bodies
void EventDrivenSimulator::get_velocities(vector<VectorNd>& qd) const
{
  // resize the vector if necessary
  qd.resize(_bodies.size());

  for (unsigned i=0; i< _bodies.size(); i++)
    _bodies[i]->get_generalized_velocity(DynamicBody::eEuler, qd[i]);
}

/// Sets the generalized coordinates of all bodies
void EventDrivenSimulator::set_coords(const vector<VectorNd>& q) const
{
  for (unsigned i=0; i< _bodies.size(); i++)
    _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q[i]);
}

/// Sets the generalized coordinates of all bodies using an interpolated value
void EventDrivenSimulator::set_coords(double t) 
{
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // do linear interpolation
    _workV = _qf[i];
    _workV -= _q0[i];
    _workV *= t;
    _workV += _q0[i];
    _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _workV);
  }
}

/// Sets the generalized velocities of all bodies
void EventDrivenSimulator::set_velocities(const vector<VectorNd>& qd) const
{
  for (unsigned i=0; i< _bodies.size(); i++)
    if (!_bodies[i]->get_kinematic())
      _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, qd[i]);
}

/// Sets the generalized velocities of all bodies using an interpolated value
void EventDrivenSimulator::set_velocities(double t) 
{
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // don't do this if the body is kinematically controlled
    if (_bodies[i]->get_kinematic())
      continue;

    // do linear interpolation
    _workV = _qdf[i];
    _workV -= _qd0[i];
    _workV *= t;
    _workV += _qd0[i];
    _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, _workV);
  }
}

/// Does a semi-implicit Euler integration
void EventDrivenSimulator::integrate_si_Euler(double step_size)
{
  VectorNd q, qd, x, dx;

  // begin timing dynamics
  tms cstart;  
  clock_t start = times(&cstart);

  // get the state-derivative for each dynamic body
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // if the body is kinematically updated, do not integrate it
    if (_bodies[i]->get_kinematic())
    {
      if (_bodies[i]->controller)
        _bodies[i]->controller(_bodies[i], current_time, _bodies[i]->controller_arg);
      continue;
    }

    // integrate the body
    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (before): " << _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (before): " << _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
    }

    // compute the velocity 
    _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q);
    _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
    x.resize(q.size()+qd.size());
    x.set_sub_vec(0, q);
    x.set_sub_vec(q.size(), qd);
    _bodies[i]->ode_both(x, current_time, step_size, &_bodies[i], dx);

    // update the velocity and position
    dx.segment(q.size(), dx.size()) *= step_size;
    qd += dx.segment(q.size(), dx.size());
    _bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd);
    _bodies[i]->get_generalized_velocity(DynamicBody::eEuler, qd);
    qd *= step_size;
    q += qd; 
    _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);

    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (after): " << _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (after): " << _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
    }
  }

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;
}

/// Steps the simulator forward by the given step size
double EventDrivenSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

  // clear timings
  dynamics_time = (double) 0.0;
  event_time = (double) 0.0;
  coldet_time = (double) 0.0;

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif
  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << std::endl;

  // TODO: should we clear the limits on each step? Compute a running
  //       mean? Compute a mean with a forgetting factor?

  // setup the time stepped
  double h = 0.0;

  // step until the requisite time has elapsed
  while (h < step_size)
  {
    // determine the maximum step according to conservative advancement
    double safe_dt = calc_CA_step(step_size-h);

    // called on integration restart
    restart: 

    // if safe_dt is effectively zero, use a semi-implicit Euler step to
    // solve events, etc.
    if (safe_dt < NEAR_ZERO)
    {
      double dt = std::min(step_size-h, max_Euler_step);
      step_si_Euler(dt);
      h += dt;

      // call the mini-callback
      if (post_mini_step_callback_fn)
        post_mini_step_callback_fn(this);

      // continue integrating
      continue;
    }

    // save the state of the system
    save_state();

    // attempt to integrate forward by safe_dt
    try
    {
      // do "smart" integration (watching for state violation) 
      integrate_DAE(safe_dt);
    }
    catch (InvalidStateException e)
    {
      // couldn't integrate that far; restart the integration with a smaller
      // step size
      safe_dt *= 0.5;
      goto restart;
    }

    // see whether there were any force or acceleration limits exceeded
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
      if (rb)
      {
        // check whether force limits were exceeded
        if (rb->force_limit_exceeded())
        {
          // reset the state of all bodies 
          restore_state();

          // attempt to integrate again using new CA info
          continue;
        }
      }
      else
      {
        ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
        if (ab->joint_accel_limit_exceeded())
        {
          // reset the state of all bodies
          restore_state();

          // attempt to integrate again using new CA info
          continue;
        }
      }
    }

    // no issues integrating; update h and call the mini-callback
    h += safe_dt;
    if (post_mini_step_callback_fn)
      post_mini_step_callback_fn(this);
  }

  // call the callback 
  if (post_step_callback_fn)
    post_step_callback_fn(this);
  
  return step_size;
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

/// Does effectively a integration of the underlying DAE's
void EventDrivenSimulator::integrate_DAE(double dt)
{
  // TODO: verify constraint violation will not affect semi-implicit step 

  // do normal integration
  integrate(dt);

  // update constraint violation after integration
  update_constraint_violations();
}

/// Updates constraint violation after integration
void EventDrivenSimulator::update_constraint_violations()
{
  // TODO: update constraint violation due to increasing interpenetration

  // update joint constraint interpenetration
  BOOST_FOREACH(DynamicBodyPtr db, _bodies)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
    if (ab)
      ab->update_joint_constraint_violations();
  }  
}

/// Computes a conservative advancement step
double EventDrivenSimulator::calc_CA_step(double dt) const
{
  const double INF = std::numeric_limits<double>::max();

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

  // TODO: do collision detection here

  return dt;
}

/// Steps the simulator forward
/*
double EventDrivenSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

  // setup the amount remaining to step
  double dt = step_size;

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif
  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << std::endl;

  // store the current generalized coordintes 
  get_coords(_q0);

  // integrate the systems forward by dt
  integrate_si_Euler(dt);

  // save the new phase coordinates 
  get_coords(_qf);
  get_velocities(_qdf);

  // methods below assume that coords/velocities of the bodies may be modified,
  // so we need to take precautions to save/restore them as necessary
  while (dt > (double) 0.0)
  {
    // look for events
    double t = find_and_handle_si_events(dt);
    if (t > dt)
      break; // no event.. finish up

    // events have been handled already; reduce dt and keep integrating
    dt -= t;

    // call the mini-callback
    if (post_mini_step_callback_fn)
      post_mini_step_callback_fn(this);

    // get the new velocities
    get_velocities(_qdf);

    // update the coordinates using the new velocities
    for (unsigned i=0; i< _q0.size(); i++)
    {
      // don't do this for kinematically updated bodies
      if (_bodies[i]->get_kinematic())
        continue;

      // update using semi-implicit integration
      _qf[i] = _qdf[i];
      _qf[i] *= dt;
      _qf[i] += _q0[i];
    }
  }

  // call the callback 
  if (post_step_callback_fn)
    post_step_callback_fn(this);
  
  return step_size;
}
*/

/// Computes forward dynamics for all bodies
void EventDrivenSimulator::calc_fwd_dyn() const
{
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
      (*db->controller)(db, current_time, db->controller_arg);

    // calculate forward dynamics at state x
    db->calc_fwd_dyn();
  }
}

/// Does a semi-implicit step 
void EventDrivenSimulator::step_si_Euler(double dt)
{
  // store the current generalized coordintes 
  get_coords(_q0);

  // integrate the systems forward by dt
  integrate_si_Euler(dt);

  // save the new phase coordinates 
  get_coords(_qf);
  get_velocities(_qdf);

  // methods below assume that coords/velocities of the bodies may be modified,
  // so we need to take precautions to save/restore them as necessary
  while (dt > (double) 0.0)
  {
    // look for events
    double t = find_and_handle_si_events(dt);
    if (t > dt)
      break; // no event.. finish up

    // events have been handled already; reduce dt and keep integrating
    dt -= t;

    // get the new velocities
    get_velocities(_qdf);

    // update the coordinates using the new velocities
    for (unsigned i=0; i< _q0.size(); i++)
    {
      _qf[i] = _qdf[i];
      _qf[i] *= dt;
      _qf[i] += _q0[i];
    }
  }
}

/// Finds and handles first impacting event(s) in [0,dt]; returns time t in [0,dt] of first impacting event(s) and advances bodies' dynamics to time t
double EventDrivenSimulator::find_events(double dt)
{
  vector<Event> cd_events, limit_events;
  typedef map<Event, double, EventCompare>::const_iterator EtolIter;

  // only for debugging purposes: verify that bodies aren't already interpenetrating
  #ifndef NDEBUG
  if (!_simulation_violated)
    check_violation();
  #endif

  // clear events 
  _events.clear();

  // begin timing for collision detection
  tms cstart;
  clock_t start = times(&cstart);

  FILE_LOG(LOG_SIMULATOR) << "-- checking for event in interval [" << (this->current_time) << ", " << (this->current_time+dt) << "] (dt=" << dt << ")" << std::endl;

  // make sure that dt is non-negative
  assert(dt >= (double) 0.0);

  // setup x0, x1
  if (!collision_detectors.empty())
  {
    _x0.resize(_q0.size());
    _x1.resize(_q0.size());
    for (unsigned i=0; i< _bodies.size(); i++)
    {
      _x0[i].first = _x1[i].first = _bodies[i];
      _x0[i].second = _q0[i];
      _x1[i].second = _qf[i];
    }
  }

  // call each collision detector
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
  {
    // indicate this is event driven
    cd->return_all_contacts = true;

    // do the collision detection routine
    cd_events.clear();
    cd->is_contact(dt, _x0, _x1, cd_events);

    // add to events
    _events.insert(_events.end(), cd_events.begin(), cd_events.end());
  }

  // tabulate times for collision detection 
  tms cstop;  
  clock_t stop = times(&cstop);
  coldet_time += (double) (stop-start)/CLOCKS_PER_SEC;

  // check each articulated body for a joint limit event
  limit_events.clear();
  find_limit_events(dt, limit_events);
  _events.insert(_events.end(), limit_events.begin(), limit_events.end());

  // sort the set of events
  std::sort(_events.begin(), _events.end()); 

  // set the "real" time for the events
  for (unsigned i=0; i< _events.size(); i++)
  {
    _events[i].t_true = current_time + _events[i].t * dt;
    EtolIter j = _event_tolerances.find(_events[i]);
    if (j != _event_tolerances.end())
      _events[i].tol = j->second;
  }

  // each group of events can be handled in x ways:
  // 1. one or more events is impacting; all events need to be handled with
  //    an impact method
  // 2. all events are separating at the velocity level; these events do not 
  //    need to be handled
  // 3. all events are resting at the velocity level; these events need to
  //    be checked at the acceleration level

  // step to first event time
  if (!_events.empty())
  {
    while (true)
    {
      // set the coordinates and velocities
      set_coords(_events.front().t);
      set_velocities(_events.front().t);

      // if all contacts at the current time are separating, remove those 
      // contacts and step to the next set of contacts
      bool all_separating = true;
      for (unsigned i=0; i< _events.size(); i++)
      {
        // look to see whether we can stop examining events
        if (std::fabs(_events[i].t - _events[0].t) > NEAR_ZERO)
          break;
        else if (!_events[i].is_separating())
        {
          all_separating = false;
          break;
        }
      }

      // if not all are separating, break out now
      if (!all_separating)
        break;

      // otherwise remove contacts not occurring simultaneously with
      // the first contacts
      vector<Event>::iterator i = _events.begin();
      while (++i != _events.end())
        if (std::fabs(i->t - _events[0].t) > NEAR_ZERO)
          break;
      _events.erase(_events.begin(), i);

      // if there are no events remaining, indicate no events
      if (_events.empty())
        return 1.0;
    } 
  }

  // check whether any events are at current time
  for (unsigned i=0; i< _events.size(); i++)
  {
    if (_events[i].t > NEAR_ZERO)
      break;

    // set event type as velocity initially
    _events[i].deriv_type = Event::eVel;

    // check whether we can encode the event as an acceleration event
    if (_events[i].determine_event_class() == Event::eZero)
      _events[i].deriv_type = Event::eAccel;
  }

  // output the events
  if (LOGGING(LOG_EVENT))
  {
    FILE_LOG(LOG_EVENT) << "Events to be processed:" << std::endl;
    for (unsigned i=0; i< _events.size(); i++)
      FILE_LOG(LOG_EVENT) << _events[i] << std::endl;
  }

  // if there are no events remaining, return now 
  if (_events.empty())
    return 1.0;

  // find the first TOI 
  return _events.front().t;
}

/// Removes events after time 0
void EventDrivenSimulator::remove_next_events()
{
  for (unsigned i=0; i< _events.size(); i++)
    if (_events[i].t > NEAR_ZERO)
    {
      _events.erase(_events.begin()+i, _events.end());
      return;
    }
}

/// Finds the next event time (after 0)
double EventDrivenSimulator::find_next_event_time() const
{
  const double INF = std::numeric_limits<double>::max();

  for (unsigned i=0; i< _events.size(); i++)
    if (_events[i].t > NEAR_ZERO)
      return _events[i].t;

  // still here? no next event time
  return INF;
}

/// Finds and handles first impacting event(s) in [0,dt]; returns time t in [0,dt] of first impacting event(s) and advances bodies' dynamics to time t
double EventDrivenSimulator::find_and_handle_si_events(double dt)
{
  vector<Event> cd_events, limit_events;
  typedef map<Event, double, EventCompare>::const_iterator EtolIter;

  // only for debugging purposes: verify that bodies aren't already interpenetrating
  #ifndef NDEBUG
  if (!_simulation_violated)
    check_violation();
  #endif

  // clear events 
  _events.clear();

  // begin timing for collision detection
  tms cstart;
  clock_t start = times(&cstart);

  FILE_LOG(LOG_SIMULATOR) << "-- checking for event in interval [" << (this->current_time) << ", " << (this->current_time+dt) << "] (dt=" << dt << ")" << std::endl;

  // make sure that dt is non-negative
  assert(dt >= (double) 0.0);

  // setup x0, x1
  if (!collision_detectors.empty())
  {
    _x0.clear();
    _x1.clear();
    for (unsigned i=0; i< _bodies.size(); i++)
    {
      if (!_bodies[i]->get_kinematic())
      {
        _x0.push_back(std::make_pair(_bodies[i], _q0[i]));
        _x1.push_back(std::make_pair(_bodies[i], _qf[i]));
      }
    }
  }

  // call each collision detector
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
  {
    // indicate this is event driven
    cd->return_all_contacts = true;

    // do the collision detection routine
    cd_events.clear();
    cd->is_contact(dt, _x0, _x1, cd_events);

    // add to events
    _events.insert(_events.end(), cd_events.begin(), cd_events.end());
  }

  // tabulate times for collision detection 
  tms cstop;  
  clock_t stop = times(&cstop);
  coldet_time += (double) (stop-start)/CLOCKS_PER_SEC;

  // check each articulated body for a joint limit event
  limit_events.clear();
  find_limit_events(dt, limit_events);
  _events.insert(_events.end(), limit_events.begin(), limit_events.end());

  // sort the set of events
  std::sort(_events.begin(), _events.end()); 

  // output the events
  if (LOGGING(LOG_EVENT))
  {
    FILE_LOG(LOG_EVENT) << "Events to be processed:" << std::endl;
    for (unsigned i=0; i< _events.size(); i++)
      FILE_LOG(LOG_EVENT) << _events[i] << std::endl;
  }

  // set the "real" time for the events
  for (unsigned i=0; i< _events.size(); i++)
  {
    _events[i].t_true = current_time + _events[i].t * dt;
    EtolIter j = _event_tolerances.find(_events[i]);
    if (j != _event_tolerances.end())
      _events[i].tol = j->second;
  }

  // find and integrate body positions to the time-of-impact
  double h = integrate_to_TOI(dt);

  // handle the events
  if (h < dt)
  {
    // if h = 0, revalidate all positions
    if (h < NEAR_ZERO)
    {
      BOOST_FOREACH(DynamicBodyPtr db, _bodies)
      {
        db->validate_position_variables();
      }
    }

    // handle the events
    handle_events();
  }

  return h;  
}

/// Finds joint limit events
void EventDrivenSimulator::find_limit_events(double dt, vector<Event>& events)
{
  // clear the vector of events
  events.clear();

  // process each articulated body, looking for joint events
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
    ab->find_limit_events(_q0[i], _qf[i], dt, std::back_inserter(events));
  }
}

/// Finds the next time-of-impact out of a set of events
double EventDrivenSimulator::integrate_to_TOI(double dt)
{
  const double INF = std::numeric_limits<double>::max();

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::integrate_to_TOI() entered with dt=" << dt << endl;

  // get the iterator start
  vector<Event>::iterator citer = _events.begin();

  // setup integration performed 
  double h = (double) 0.0;

  // loop while the iterator does not point to the end -- may need several
  // iterations b/c there may be no impacting events in a group 
  while (citer != _events.end())
  {
    // set tmin
    double tmin = citer->t*dt;
    FILE_LOG(LOG_SIMULATOR) << "  -- integrate_to_TOI() while loop, current time=" << current_time << " tmin=" << tmin << endl;

    // check for exit
    if (tmin > dt)
    {
      FILE_LOG(LOG_SIMULATOR) << "    " << tmin << " > " << dt << " --> exiting now w/o events" << endl;
      FILE_LOG(LOG_SIMULATOR) << "    .... but first, integrating bodies forward by " << (dt-h) << std::endl;

      // events vector no longer valid; clear it
      _events.clear();

      // set the coordinates
      for (unsigned i=0; i< _bodies.size(); i++)
      {
        if (_bodies[i]->get_kinematic())
          continue;
        _qf[i] = _qdf[i];
        _qf[i] *= dt;
        _qf[i] += _q0[i];
        _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _qf[i]);
        _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, _qdf[i]);
      }

      // update current_time
      current_time += dt;

      return INF;
    }

    // "integrate" starting coordinates to tmin
    h += tmin;
    for (unsigned i=0; i< _q0.size(); i++)
    {
      if (_bodies[i]->get_kinematic())
        continue;
      _qf[i] = _qdf[i];
      _qf[i] *= h;
      _qf[i] += _q0[i];
      _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _qf[i]);
      _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, _qdf[i]);
    }
    FILE_LOG(LOG_SIMULATOR) << "    current time is " << current_time << endl;
    FILE_LOG(LOG_SIMULATOR) << "    tmin (time to next event): " << tmin << endl;
    FILE_LOG(LOG_SIMULATOR) << "    moving forward by " << h << endl;

    // check for impacting event
    bool impacting = citer->determine_event_class() == Event::eNegative;

    // find all events at the same time as the event we are examining
    for (citer++; citer != _events.end(); citer++)
    {
      // see whether we are done
      if (citer->t*dt > tmin + std::numeric_limits<double>::epsilon())
        break;

      // see whether this event is impacting (if we don't yet have an
      // impacting event)
      if (!impacting)
        impacting = citer->determine_event_class() == Event::eNegative; 
    }

    // see whether we are done
    if (impacting)
    {
      // remove remainder of events
      _events.erase(citer, _events.end());

      // step positions to h (note that we we'll no longer need current value
      // of _qdf) 
      for (unsigned i=0; i< _q0.size(); i++)
      {
        if (_bodies[i]->get_kinematic())
          continue;
        _qf[i] = _qdf[i];
        _qf[i] *= h;
        _qf[i] += _q0[i];
        _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _qf[i]);
        _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, _qdf[i]);
      }

      // update current time
      current_time += h;

      return h;
    }
    else
      citer = _events.erase(citer, _events.end());
  }

  // contact map is empty, no contacts
  FILE_LOG(LOG_SIMULATOR) << "-- integrate_to_TOI(): no impacts detected; integrating forward by " << dt << endl;

  // events vector is no longer valid; clear it
  _events.clear();

  // set the coordinates (velocities are already set)
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    if (_bodies[i]->get_kinematic())
      continue;
    _qf[i] = _qdf[i];
    _qf[i] *= dt;
    _qf[i] += _q0[i];
    _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, _qf[i]);
    _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, _qdf[i]);
  }

  // update current_time
  current_time += dt;

  return INF;
}

/// Checks the simulator for a (contact/joint limit) violation
void EventDrivenSimulator::check_violation()
{
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
  {
    // do the collision detection routine
    if (cd->is_collision((double) 0.0))
    {
      if (!_simulation_violated)
      {
        std::cerr << "EventDrivenSimulator::is_contact() warning: detected interpenetrating geometries!" << endl;
        std::cerr << "  -- current time: " << current_time << "  NOTE: fidelity of simulation is no longer assured." << endl;
      }
      _simulation_violated = true;

/*
      // detailed contact information
      BOOST_FOREACH(const CollidingTriPair& ctp, cd->colliding_tris)
      {
        std::cerr << "    interpenetrating pair: " << endl;
        std::cerr << "      -- " << ctp.geom1->id << " (from " << ctp.geom1->get_single_body()->id << ")" << endl;
        std::cerr << "      -- " << ctp.geom2->id << " (from " << ctp.geom2->get_single_body()->id << ")" << endl;

        // get the triangles
        Triangle t1 = Triangle::transform(ctp.mesh1->get_triangle(ctp.tri1), ctp.geom1->get_transform());
        Triangle t2 = Triangle::transform(ctp.mesh2->get_triangle(ctp.tri2), ctp.geom2->get_transform());
        list<Vector3> isects;
        CompGeom::intersect_tris(t1, t2, std::back_inserter(isects));
        std::cerr << "      t1: " << t1 << std::endl;
        std::cerr << "      t2: " << t2 << std::endl;
        BOOST_FOREACH(const Vector3& point, isects)
        {
          std::cerr << "        isect: " << point << std::endl;
        }
      }
*/
      // standard contact information
      BOOST_FOREACH(sorted_pair<CollisionGeometryPtr> cg_pair, cd->colliding_pairs)
      {
        std::cerr << "    interpenetrating pair: " << endl;
        std::cerr << "      -- " << cg_pair.first->id << " (from " << cg_pair.first->get_single_body()->id << ")" << endl;
        std::cerr << "      -- " << cg_pair.second->id << " (from " << cg_pair.second->get_single_body()->id << ")" << endl;
      }
    }
  }
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

  // clear list of collision detectors
  collision_detectors.clear();

  // read the maximum time to process events, if any
  XMLAttrib* max_event_time_attrib = node->get_attrib("max-event-time");
  if (max_event_time_attrib)
    max_event_time = max_event_time_attrib->get_real_value(); 

  // read the maximum Euler step
  XMLAttrib* max_Euler_step_attrib = node->get_attrib("max-Euler-step");
  if (max_Euler_step_attrib)
    max_Euler_step = max_Euler_step_attrib->get_real_value();

  // read the error tolerances
  XMLAttrib* rel_tol_attrib = node->get_attrib("rel-err-tol");
  XMLAttrib* abs_tol_attrib = node->get_attrib("abs-err-tol");
  if (rel_tol_attrib)
    rel_err_tol = rel_tol_attrib->get_real_value();
  if (abs_tol_attrib)
    abs_err_tol = abs_tol_attrib->get_real_value();

  // get the collision detector, if specified
  XMLAttrib* coldet_attrib = node->get_attrib("collision-detector-id");
  if (coldet_attrib)
  {
    // get the ID of the collision detector
    const std::string& id = coldet_attrib->get_string_value(); 

    // find a collision detector
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "EventDrivenSimulator::load_from_xml() - could not find";
      std::cerr << std::endl << "  collision detector w/ID: " << id;
      std::cerr << " from offending node: " << std::endl << *node;
    }
    else
    {
      // make sure that it is castable to a collision detector before we
      // save the pointer
      shared_ptr<CollisionDetection> coldet = dynamic_pointer_cast<CollisionDetection>(id_iter->second);
      if (coldet)
      {
        collision_detectors.push_back(coldet);
        coldet->simulator = get_this();
      }
    }
  }

  // read in any CollisionDetection nodes
  child_nodes = node->find_child_nodes("CollisionDetector");
  BOOST_FOREACH(shared_ptr<const XMLTree> child_node, child_nodes)
  {
    XMLAttrib* id_attrib = child_node->get_attrib("id");
    if (!id_attrib)
      continue;

    // get the ID of the collision detector
    const std::string& id = id_attrib->get_string_value(); 

    // find a collision detector
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "EventDrivenSimulator::load_from_xml() - could not find";
      std::cerr << std::endl << "  collision detector w/ID: " << id;
      std::cerr << " from offending node: " << std::endl << *child_node;
    }
    else
    {
      // make sure that it is castable to a collision detector before we
      // save the pointer
      shared_ptr<CollisionDetection> coldet = dynamic_pointer_cast<CollisionDetection>(id_iter->second);
      if (coldet)
      {
        collision_detectors.push_back(coldet);
        coldet->simulator = get_this();
      }
    }
  }

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
}

/// Implements Base::save_to_xml()
void EventDrivenSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call Simulator's save method first
  Simulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "EventDrivenSimulator";

  // save the maximum event time
  node->attribs.insert(XMLAttrib("max-event-time", max_event_time));

  // save the maximum Euler step
  node->attribs.insert(XMLAttrib("max-Euler-step", max_Euler_step));

  // save the error tolerances
  node->attribs.insert(XMLAttrib("rel-err-tol", rel_err_tol));
  node->attribs.insert(XMLAttrib("abs-err-tol", abs_err_tol));

  // save the IDs of the collision detectors, if any 
  BOOST_FOREACH(shared_ptr<CollisionDetection> c, collision_detectors)
  {
    XMLTreePtr new_node(new XMLTree("CollisionDetector"));
    new_node->attribs.insert(XMLAttrib("id", c->id));
    node->add_child(new_node);
    shared_objects.push_back(c);
  }

  // save all ContactParameters
  for (map<sorted_pair<BasePtr>, shared_ptr<ContactParameters> >::const_iterator i = contact_params.begin(); i != contact_params.end(); i++)
  {
    XMLTreePtr new_node(new XMLTree("ContactParameters"));
    node->add_child(new_node);
    i->second->save_to_xml(new_node, shared_objects);
  }
}

/// Outputs this class data to the stream
/**
 * This method outputs all of the low-level details to the stream; if
 * serialization is desired, use save_to_xml() instead.
 * \sa save_to_xml()
 */
void EventDrivenSimulator::output_object_state(std::ostream& out) const
{
  // indicate the object type
  out << "EventDrivenSimulator object" << std::endl; 

  // output contact parameters
  out << "  contact parameters: " << std::endl;
  for  (map<sorted_pair<BasePtr>, shared_ptr<ContactParameters> >::const_iterator i = contact_params.begin(); i != contact_params.end(); i++)
  {
    out << "   object1: " << i->first.first << "  object2: ";
    out << i->first.second << "  parameters: " << i->second << std::endl;
  }

  // output collision detection pointers
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
    out << "  collision detector: " << cd << std::endl;

  // output event impulse callback function
   out << "  event post impulse callback fn: " << event_post_impulse_callback_fn << std::endl;

  // output event impulse callback data
   out << "  event post impulse callback data: " << event_post_impulse_callback_data << std::endl;

  // output event callback function
   out << "  event callback fn: " << event_callback_fn << std::endl;

  // output event callback data
   out << "  event callback data: " << event_callback_data << std::endl;
}

