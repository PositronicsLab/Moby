/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <sys/times.h>
#include <Moby/XMLTree.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/DynamicBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/VariableStepIntegrator.h>
#include <Moby/EventDrivenSimulator.h>

#define DO_TIMING

using namespace Moby;
using std::endl;
using std::list;
using std::vector;
using std::map;
using std::make_pair;
using std::multimap;
using std::pair;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// Default constructor
EventDrivenSimulator::EventDrivenSimulator()
{
  constraint_violation_tolerance = std::pow(std::numeric_limits<Real>::epsilon(), (Real) 0.25);
  toi_tolerance = NEAR_ZERO;
  max_Zeno_step = std::numeric_limits<Real>::max();
  event_callback_fn = NULL;
  event_post_impulse_callback_fn = NULL;
  _simulation_violated = false;
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

/// Handles events
void EventDrivenSimulator::handle_events()
{
  // call the callback function, if any
  if (event_callback_fn)
    (*event_callback_fn)(_events, event_callback_data);

  // preprocess events
  for (unsigned i=0; i< _events.size(); i++)
    preprocess_event(_events[i]);

  // compute impulses here...
  _impact_event_handler.process_events(_events, constraint_violation_tolerance);

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
void EventDrivenSimulator::get_coords_and_velocities(vector<pair<VectorN, VectorN> >& q) const
{
  // resize the vector if necessary
  q.resize(_bodies.size());

  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->get_generalized_coordinates(DynamicBody::eRodrigues, q[i].first);
    _bodies[i]->get_generalized_velocity(DynamicBody::eRodrigues, q[i].second);
  }
}

/// Sets the coords and velocities of all bodies
/**
 * Note that we set velocities in this manner (rather than using current
 * velocities) because our event finding methods assume constant velocities
 * over a time interval. Linearly interpolating the positions and velocities
 * (as we do) will result in warped dynamics.
 * \param q0 the coords at time 0
 * \param q1 the coords at time dt
 * \param t the mixture [0,1] to set the coordinates and velocities (t=0 is
 *        equivalent to time 0; t=1 equivalent to time dt)
 */
void EventDrivenSimulator::set_coords_and_velocities(const vector<pair<VectorN, VectorN> >& q0, const vector<pair<VectorN, VectorN> >& q1, Real t) const
{
  SAFESTATIC VectorN qx, qy;

  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // compute and set the generalized coordinates
    qx.copy_from(q0[i].first) *= ((Real) 1.0 - t);
    qy.copy_from(q1[i].first) *= t;
    qx += qy;
    _bodies[i]->set_generalized_coordinates(DynamicBody::eRodrigues, qx);
    FILE_LOG(LOG_SIMULATOR) << " -- set_coords_and_velocities() setting body " << _bodies[i]->id << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - state at q0: " << q0[i].first << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - state at q1: " << q1[i].first << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - generalized coords: " << qx << endl;

    // compute the generalized velocities using finite differencing
    qx.copy_from(q0[i].second) *= ((Real) 1.0 - t);
    qy.copy_from(q1[i].second) *= t;
    qx += qy;
    FILE_LOG(LOG_SIMULATOR) << "    - velocity at q0: " << q0[i].second << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - velocity at q1: " << q1[i].second << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - generalized vels: " << qx << endl;
    _bodies[i]->set_generalized_velocity(DynamicBody::eRodrigues, qx);
  }
}

/// Determines whether the given event will impact using only the generalized coordinates
/**
 * \pre Generalized coordinates of the bodies involved in the event are at the
 *      time of the event
 */
bool EventDrivenSimulator::will_impact(Event& e, const vector<pair<VectorN, VectorN> >& q0, const vector<pair<VectorN, VectorN> >& q1, Real dt) const
{
  SAFESTATIC VectorN qx, qy;

  // get the bodies of the event
  vector<DynamicBodyPtr> bodies;
  e.get_super_bodies(std::back_inserter(bodies));
  std::sort(bodies.begin(), bodies.end());

  // set the velocities of the bodies involved in the event
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // if the body is not involved in the event, skip it
    if (!std::binary_search(bodies.begin(), bodies.end(), _bodies[i]))
      continue;

    // set the velocity of the body
    qx.copy_from(q1[i].second) /= dt;
    qy.copy_from(q0[i].second) /= dt;
    qx -= qy;
    _bodies[i]->set_generalized_velocity(DynamicBody::eRodrigues, qx);
  }

  // determine whether the event reports as impacting
  bool impacting = e.is_impacting(constraint_violation_tolerance);

  // reset the velocities of the bodies
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // if the body is not involved in the event, skip it
    if (!std::binary_search(bodies.begin(), bodies.end(), _bodies[i]))
      continue;

    // set the (interpolated) velocity of the body
    qx.copy_from(q1[i].second) *= e.t;
    qy.copy_from(q0[i].second) *= ((Real) 1.0 - e.t);
    qx += qy;
    _bodies[i]->set_generalized_velocity(DynamicBody::eRodrigues, qx);
  }

  return impacting;
}

/// Sets the coords and velocities of all bodies
/**
 * \param q the coords and velocities
 */
void EventDrivenSimulator::set_coords_and_velocities(const vector<pair<VectorN, VectorN> >& q) const
{
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // set the generalized coordinates
    _bodies[i]->set_generalized_coordinates(DynamicBody::eRodrigues, q[i].first);
    _bodies[i]->set_generalized_velocity(DynamicBody::eRodrigues, q[i].second);
    FILE_LOG(LOG_SIMULATOR) << " -- set_coords_and_velocities() setting body " << _bodies[i]->id << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - generalized coords: " << q[i].first << endl;
    FILE_LOG(LOG_SIMULATOR) << "    - generalized velocities: " << q[i].second << endl;
  }
}

/// Copies one vector of generalized coordinates and velocities to another
void EventDrivenSimulator::copy(const vector<pair<VectorN, VectorN> >& source, vector<pair<VectorN, VectorN> >& dest)
{
  dest.resize(source.size());
  for (unsigned i=0; i< source.size(); i++)
  {
    dest[i].first.copy_from(source[i].first);
    dest[i].second.copy_from(source[i].second);
  }
}

/// Steps the simulator forward
Real EventDrivenSimulator::step(Real step_size)
{
  SAFESTATIC vector<shared_ptr<void> > x, xplus;
  SAFESTATIC vector<pair<VectorN, VectorN> > q0, q1, qstar;
  const Real INF = std::numeric_limits<Real>::max();

  // setup the amount remaining to step
  Real dt = step_size;

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif
  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << std::endl;

  // methods below assume that coords/velocities of the bodies may be modified,
  // so we need to take precautions to save/restore them as necessary
  while (dt > (Real) 0.0)
  {
    // get the current generalized coordinates and velocities
    get_coords_and_velocities(q0);

    // integrate the systems forward by dt
    integrate(dt);

    // save the current generalized coordinates and velocities
    get_coords_and_velocities(q1);

    // look for events in [0, dt], advance all bodies to the time of event,
    // and handle the event(s)
    bool Zeno = false;
    Real t = find_and_handle_events(dt, q0, q1, Zeno);
    if (t > dt)
      break; // no event.. finish up
    else if (Zeno) 
    {
      // move to time of designated Zeno point
      Real h = std::min(max_Zeno_step, dt);
      handle_Zeno_point(h, q0, q1);
      t += h;
      if (t > dt)    // don't want to accidentally step clock too far
        t = dt;
    }

    // events have been handled already; reduce dt and keep integrating
    dt -= t;
    current_time += t;
  }

  // update the current time
  current_time += dt;

  return step_size;
}

/// Handles Zeno points in the simulator
/**
 * The theory behind this is that, when a Zeno point is encountered, some force
 * components and some constraints "cancel" (for example, a box resting on a
 * plane experiences a downward force from gravity, which is immediately
 * counteracted by the plane constraint force). Event handling shapes the
 * evolution of the system; going back to the box example, if a horizontal
 * force is applied at the same time as the gravitational force, contact 
 * handling will counteract the gravitational force but may leave the 
 * horizontal force unaltered (if there is no friction present, for example).
 * So we can update the position of the system a small step forward using the 
 * velocity determined via the event handling method, thus reproducing the
 * desired behavior of the system (if the step is sufficiently small).
 */
void EventDrivenSimulator::handle_Zeno_point(Real dt, const vector<pair<VectorN, VectorN> >& q0, vector<pair<VectorN, VectorN> >& q1)
{
  // NOTE: _events must be composed strictly of Zeno point events 
  // (this is ensured by find_TOI())
  FILE_LOG(LOG_SIMULATOR) << "Zeno point detected! handling it..." << endl;

  // determine treated bodies for all events
  vector<DynamicBodyPtr> treated_bodies;
  for (unsigned i=0; i< _events.size(); i++)
    _events[i].get_super_bodies(std::back_inserter(treated_bodies));

  // remove duplicates
  std::sort(treated_bodies.begin(), treated_bodies.end());
  treated_bodies.erase(std::unique(treated_bodies.begin(), treated_bodies.end()), treated_bodies.end());

  // for each body in the simulator
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // each Zeno body will already have the proper velocity; use that to 
    // update q1; non-Zeno bodies will be updated to their already determined
    // q1
    if (std::binary_search(treated_bodies.begin(), treated_bodies.end(), _bodies[i]))
    {
      // get the body's current velocity (it was just treated) and use that to
      // update q1
      _bodies[i]->get_generalized_velocity(DynamicBody::eRodrigues, q1[i].second);
      q1[i].first.copy_from(q1[i].second) *= dt;
      q1[i].first += q0[i].first;
    }
    _bodies[i]->set_generalized_coordinates(DynamicBody::eRodrigues, q1[i].first);
  }
}

/// Determines (super) bodies treated in events
void EventDrivenSimulator::determine_treated_bodies(list<list<Event*> >& groups, vector<DynamicBodyPtr>& bodies)
{
  // clear the vector of bodies
  bodies.clear();

  // get the super bodies
  BOOST_FOREACH(list<Event*>& l, groups)
    BOOST_FOREACH(Event* e, l)
      e->get_super_bodies(std::back_inserter(bodies));

  // sort the vector of bodies and remove duplicates
  std::sort(bodies.begin(), bodies.end());
  bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());
}

/// Finds and handles first impacting event(s) in [0,dt]; returns time t in [0,dt] of first impacting event(s) and advances bodies' dynamics to time t
/**
 * Note that the velocities at both endpoints of the interval are input, but
 * our event checking mechanisms currently assume that velocities over the
 * interval are constant. The event checking mechanisms linearly interpolate
 * the generalized coordinates over the time interval. The dynamics may be 
 * "distorted" as a result b/c the generalized coordinates should change due to 
 * the velocity (rather than via linear interpolation). Two things to note,
 * however: 1) events will not be missed as long as the event finders can
 * search over [q0,q1] and 2) this allows us to handle non-explicit 
 * integration. 
 */
Real EventDrivenSimulator::find_and_handle_events(Real dt, const vector<pair<VectorN, VectorN> >& q0, const vector<pair<VectorN, VectorN> >& q1, bool& Zeno)
{
  SAFESTATIC VectorN qx;
  vector<Event> cd_events, limit_events;
  SAFESTATIC vector<pair<DynamicBodyPtr, VectorN> > x0, x1;

  FILE_LOG(LOG_SIMULATOR) << "-- checking for event in interval [" << this->current_time << ", " << (this->current_time+dt) << "] (dt=" << dt << ")" << std::endl;

  // make sure that dt is non-negative
  assert(dt >= (Real) 0.0);

  // only for debugging purposes: verify that bodies aren't already interpenetrating
  #ifndef NDEBUG
  if (!_simulation_violated)
    check_violation();
  #endif

  // clear events 
  _events.clear();

  #ifdef DO_TIMING
  tms start;
  times(&start);
  #endif

  // setup x0, x1
  if (!collision_detectors.empty())
  {
    x0.resize(q0.size());
    x1.resize(q0.size());
    for (unsigned i=0; i< _bodies.size(); i++)
    {
      x0[i].first = x1[i].first = _bodies[i];
      x0[i].second.copy_from(q0[i].first);
      x1[i].second.copy_from(q1[i].first);
    }
  }

  // call each collision detector
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
  {
    // indicate this is event driven
    cd->return_all_contacts = true;

    // do the collision detection routine
    cd_events.clear();
    cd->is_contact(dt, x0, x1, cd_events);

    // add to events
    _events.insert(_events.end(), cd_events.begin(), cd_events.end());
  }

  // check each articulated body for a joint limit event
  limit_events.clear();
  find_limit_events(q0, q1, dt, limit_events);
  _events.insert(_events.end(), limit_events.begin(), limit_events.end());

  // sort the set of events
  std::sort(_events.begin(), _events.end()); 

  // set the "real" time for the events
  for (unsigned i=0; i< _events.size(); i++)
    _events[i].t_true = current_time + _events[i].t * dt;

  #ifdef DO_TIMING
  tms stop;  
  times(&stop);
  Real utime = (Real) (stop.tms_utime-start.tms_utime)/CLOCKS_PER_SEC;
  Real stime = (Real) (stop.tms_stime-start.tms_stime)/CLOCKS_PER_SEC;
  std::ofstream out("coldet-timing", std::ofstream::app);
  out << current_time << " " << utime << " " << stime << endl;
  #endif

  // find and "integrate" to the time-of-impact
  Real TOI = find_TOI(dt, q0, q1);

  // check for Zeno point
  if (TOI < std::numeric_limits<Real>::epsilon())
  {
    // if all events are resting or separating, we have a Zeno point
    Zeno = true;
    for (unsigned i=0; i< _events.size(); i++)
      if (_events[i].is_impacting(constraint_violation_tolerance))
      {
        Zeno = false;
        break;
      }
  }
  else
    Zeno = true;

  // if it's a Zeno point, update the velocities and redetermine the event classes
  if (Zeno)
  {
    // determine bodies in the events
    vector<DynamicBodyPtr> treated_bodies;
    for (unsigned i=0; i< _events.size(); i++)
      _events[i].get_super_bodies(std::back_inserter(treated_bodies));

    // remove duplicates
    std::sort(treated_bodies.begin(), treated_bodies.end());
    treated_bodies.erase(std::unique(treated_bodies.begin(), treated_bodies.end()), treated_bodies.end());
 
    // set velocities for bodies in events  
    for (unsigned i=0; i< _bodies.size(); i++)
      if (std::binary_search(treated_bodies.begin(), treated_bodies.end(), _bodies[i]))
        _bodies[i]->set_generalized_velocity(DynamicBody::eRodrigues, q1[i].second);
  }

  // finally, handle the events
  if (TOI <= dt)
    handle_events();

  return TOI;
}

/// Finds joint limit events
void EventDrivenSimulator::find_limit_events(const vector<pair<VectorN, VectorN> >& q0, const vector<pair<VectorN, VectorN> >& q1, Real dt, vector<Event>& events)
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
    
    // get limit events in [t, t+dt] (if any)
    ab->find_limit_events(q0[i].first, q1[i].first, dt, std::back_inserter(events));
  }
}

/// Finds the next time-of-impact out of a set of events
Real EventDrivenSimulator::find_TOI(Real dt, const vector<pair<VectorN, VectorN> >& q0, const vector<pair<VectorN, VectorN> >& q1)
{
  const Real INF = std::numeric_limits<Real>::max();

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::find_TOI() entered with dt=" << dt << endl;

  // get the iterator start
  vector<Event>::iterator citer = _events.begin();

  // setup integration performed 
  Real h = (Real) 0.0;

  // loop while the iterator does not point to the end -- may need several
  // iterations b/c there may be no impacting events in a group 
  while (citer != _events.end())
  {
    // set tmin
    Real tmin = citer->t*dt;

    FILE_LOG(LOG_SIMULATOR) << "  -- find_TOI() while loop, current time=" << current_time << " tmin=" << tmin << endl;

    // check for exit
    if (tmin > dt)
    {
      FILE_LOG(LOG_SIMULATOR) << "    " << tmin << " > " << dt << " --> exiting now w/o events" << endl;
      FILE_LOG(LOG_SIMULATOR) << "    .... but first, integrating bodies forward by " << (dt-h) << std::endl;

      // events vector no longer valid; clear it
      _events.clear();

      // set the coordinates and velocities
      set_coords_and_velocities(q1);

      return INF;
    }

    // integrate to tmin
    h += tmin;
    set_coords_and_velocities(q0, q1, h/dt);
    FILE_LOG(LOG_SIMULATOR) << "    current time is " << current_time << endl;
    FILE_LOG(LOG_SIMULATOR) << "    tmin (time to next event): " << tmin << endl;
    FILE_LOG(LOG_SIMULATOR) << "    moving forward by " << h << endl;

    // check for impacting event
    bool impacting = (citer->is_impacting(constraint_violation_tolerance) || will_impact(*citer, q0, q1, dt));

    // find all events at the same time as the event we are examining
    for (citer++; citer != _events.end(); citer++)
    {
      // see whether we are done
      if (citer->t*dt > tmin + toi_tolerance)
        break;

      // see whether this event is impacting (if we don't yet have an
      // impacting event)
      if (!impacting)
        impacting = (citer->is_impacting(constraint_violation_tolerance) || will_impact(*citer, q0, q1, dt)); 
    }

    // see whether we are done
    if (impacting)
    {
      // remove remainder of events
      _events.erase(citer, _events.end());
      return h;
    }
    else
      citer = _events.erase(citer, _events.end());
  }

  // contact map is empty, no contacts
  FILE_LOG(LOG_SIMULATOR) << "-- find_TOI(): no impacts detected; integrating forward by " << dt << endl;

  // events vector is no longer valid; clear it
  _events.clear();

  // set the coordinates and velocities
  set_coords_and_velocities(q1);

  return INF;
}

/// Checks the simulator for a (contact/joint limit) violation
void EventDrivenSimulator::check_violation()
{
  BOOST_FOREACH(shared_ptr<CollisionDetection> cd, collision_detectors)
  {
    // do the collision detection routine
    if (cd->is_collision((Real) 0.0))
    {
      if (!_simulation_violated)
      {
        std::cerr << "EventDrivenSimulator::is_contact() warning: detected interpenetrating geometries!" << endl;
        std::cerr << "  -- current time: " << current_time << "  NOTE: fidelity of simulation is no longer assured." << endl;
      }
      _simulation_violated = true;

      // detailed contact information
/*
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
void EventDrivenSimulator::load_from_xml(XMLTreeConstPtr node, map<std::string, BasePtr>& id_map)
{
  list<XMLTreeConstPtr> child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify node name b/c this is abstract class
  assert(strcasecmp(node->name.c_str(), "EventDrivenSimulator") == 0);

  // first, load all data specified to the Simulator object
  Simulator::load_from_xml(node, id_map);

  // clear list of collision detectors
  collision_detectors.clear();

  // get the value of the TOI tolerance
  const XMLAttrib* toitol_attrib = node->get_attrib("TOI-tolerance");
  if (toitol_attrib)
    toi_tolerance = toitol_attrib->get_real_value();

  // get the maximum Zeno step 
  const XMLAttrib* maxZeno_attrib = node->get_attrib("max-Zeno-step");
  if (maxZeno_attrib)
    max_Zeno_step = maxZeno_attrib->get_real_value();

  // get the value of the Zeno velocity tolerance
  const XMLAttrib* zenotol_attrib = node->get_attrib("Zeno-vel-tolerance");
  if (zenotol_attrib)
    Zeno_vel_tolerance = zenotol_attrib->get_real_value();

  // get the value of collision tolerance, if specified
  const XMLAttrib* ivtol_attrib = node->get_attrib("constraint-violation-tolerance");
  if (ivtol_attrib)
    constraint_violation_tolerance = ivtol_attrib->get_real_value();

  // get the collision detector, if specified
  const XMLAttrib* coldet_attrib = node->get_attrib("collision-detector-id");
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
  BOOST_FOREACH(XMLTreeConstPtr child_node, child_nodes)
  {
    const XMLAttrib* id_attrib = child_node->get_attrib("id");
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
  for (list<XMLTreeConstPtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    boost::shared_ptr<ContactParameters> cd(new ContactParameters);
    cd->load_from_xml(*i, id_map);
    contact_params[cd->objects] = cd;
  }
}

/// Implements Base::save_to_xml()
void EventDrivenSimulator::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // call Simulator's save method first
  Simulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "EventDrivenSimulator";

  // save the TOI tolerance
  node->attribs.insert(XMLAttrib("TOI-tolerance", toi_tolerance));

  // save the maximum Zeno step 
  node->attribs.insert(XMLAttrib("max-Zeno-step", max_Zeno_step));

  // save the Zeno velocity tolerance
  node->attribs.insert(XMLAttrib("Zeno-vel-tolerance", Zeno_vel_tolerance));

  // save the value of impact velocity
  node->attribs.insert(XMLAttrib("constraint-violation-tolerance", constraint_violation_tolerance));

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

  // output the TOI tolerance
  out << "  time of impact tolerance: " << toi_tolerance << std::endl;

  // output scalars
  out << "  impact velocity tolerance: " << constraint_violation_tolerance << std::endl;

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

