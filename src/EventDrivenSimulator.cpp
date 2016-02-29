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
#include <Moby/ControlledBody.h>
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
  // setup the standard Euler step
  euler_step = 1e-3;
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
    return (lj1 < lj2);
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
      return (cg11+cg12 < cg21+cg22);
    }
    else
      return false; // limits returned before contacts
  }
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

    _ip_tolerances[make_sorted_pair(pdi.a, pdi.b)] = std::min(pdi.dist, -NEAR_ZERO);
  }

  // update joint constraint interpenetration
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
    if (ab)
      ab->update_joint_constraint_violations();
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::update_constraint_violations() exited" << std::endl;
}

/// Handles sustained rigid unilateral constraints 
void EventDrivenSimulator::calc_rigid_sustained_unilateral_constraint_forces()
{
  // if there are no constraints, quit now
  if (_rigid_constraints.empty())
    return;

  // remove velocity-level constraints
  std::vector<UnilateralConstraint> velocity_constraints;
  for (unsigned i=0; i< _rigid_constraints.size(); )
  {
    if (_rigid_constraints[i].deriv_type == UnilateralConstraint::eVel)
    {
      velocity_constraints.push_back(_rigid_constraints[i]);
      _rigid_constraints[i] = _rigid_constraints.back();
      _rigid_constraints.pop_back();
    }
    else
      i++;  
  }

  // call the callback function, if any
  if (constraint_callback_fn)
    (*constraint_callback_fn)(_rigid_constraints, constraint_callback_data);

  // preprocess constraints
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    preprocess_constraint(_rigid_constraints[i]);

  // compute forces here...
  _rigid_unilateral_constraint_handler.process_constraints(_rigid_constraints);

  // call the post-force application callback, if any
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_rigid_constraints, constraint_post_callback_data);

  // add velocity level constraints back in
  _rigid_constraints.insert(_rigid_constraints.end(), velocity_constraints.begin(), velocity_constraints.end());
}

/// Computes the ODE for systems with sustained unilateral constraints
VectorNd& EventDrivenSimulator::ode_sustained_constraints(const VectorNd& x, double t, double dt, void* data, VectorNd& dx)
{
  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::ode_sustained_constraints(t=" << t << ") entered" << std::endl;

  // get the simulator
  shared_ptr<EventDrivenSimulator>& s = *((shared_ptr<EventDrivenSimulator>*) data);

  // initialize the ODE index
  unsigned idx = 0;

  // resize dx
  dx.resize(x.size());

  // loop through all bodies, preparing to compute the ODE
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // get the number of generalized coordinates and velocities
    const unsigned NGC = db->num_generalized_coordinates(DynamicBodyd::eEuler);
    const unsigned NGV = db->num_generalized_coordinates(DynamicBodyd::eSpatial);

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
  s->calc_pairwise_distances();
  double min_dist = s->check_pairwise_constraint_violations(t);

  // find unilateral constraints
  s->find_unilateral_constraints(s->contact_dist_thresh);
  if (s->_rigid_constraints.empty())
  {
    if (LOGGING(LOG_SIMULATOR))
    {
      FILE_LOG(LOG_SIMULATOR) << " *** constraints vector is empty " << std::endl;

      // find contact constraints
      BOOST_FOREACH(const PairwiseDistInfo& pdi, s->_pairwise_distances)
        FILE_LOG(LOG_SIMULATOR) << " -- signed distance between " << pdi.a->get_single_body()->id << " and " << pdi.b->get_single_body()->id << ": " << pdi.dist << std::endl;
    }
  }

  // check velocity violations for constraints
  s->check_constraint_velocity_violations(t);

  // convert rigid constraints to acceleration constraints
  for (unsigned i=0; i< s->_rigid_constraints.size(); i++)
    if (std::fabs(s->_rigid_constraints[i].calc_constraint_vel()) < s->_rigid_constraints[i].sustained_tol)
      s->_rigid_constraints[i].deriv_type = UnilateralConstraint::eAccel;

  // loop through all bodies, computing forward dynamics
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
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
  vector<ControlledBodyPtr> bodies;
  BOOST_FOREACH(const UnilateralConstraint& e, s->_rigid_constraints)
    e.get_super_bodies(std::back_inserter(bodies));
  BOOST_FOREACH(const UnilateralConstraint& e, s->_compliant_constraints)
    e.get_super_bodies(std::back_inserter(bodies));
  std::sort(bodies.begin(), bodies.end());
  bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());

   // recompute forward dynamics for bodies in constraints
   BOOST_FOREACH(ControlledBodyPtr body, bodies)
     body->calc_fwd_dyn();

  // report accelerations
  if (LOGGING(LOG_CONSTRAINT))
  {
    FILE_LOG(LOG_CONSTRAINT) << " -- post-contact acceleration (all contacts): " << std::endl;
    BOOST_FOREACH(const UnilateralConstraint& e, s->_rigid_constraints)
      FILE_LOG(LOG_CONSTRAINT) << e;
  }

  // debugging code for checking numerical acceleration 
  #ifndef NDEBUG
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
  #endif

  // reset idx
  idx = 0;

  // loop through all bodies, computing the ODE
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // get the number of generalized coordinates and velocities
    const unsigned NGC = db->num_generalized_coordinates(DynamicBodyd::eEuler);
    const unsigned NGV = db->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // get dx for the body
    SharedVectorNd dxsub = dx.segment(idx, idx+NGC+NGV);

    // compute the ODE
    db->ode(t, dt, &db, dxsub);

    FILE_LOG(LOG_SIMULATOR) << " ODE evaluation for body " << db->id << ": " << dxsub << std::endl;

    // update idx
    idx += NGC+NGV;
  }

  FILE_LOG(LOG_SIMULATOR) << "EventDrivenSimulator::ode_sustained_constraints(t=" << t << ") exited" << std::endl;

  // return the ODE
  return dx;
}

/// Steps the simulator forward by the given step size
double EventDrivenSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

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

    FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << " by " << (step_size - h) << std::endl;
    if (LOGGING(LOG_SIMULATOR))
    {
      VectorNd q, qd;
      BOOST_FOREACH(ControlledBodyPtr db, _bodies)
      {
        db->get_generalized_coordinates(DynamicBodyd::eEuler, q);
        db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
        FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " Euler coordinates (before): " << q << std::endl;
        FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " spatial velocity (before): " << qd << std::endl;
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
    find_unilateral_constraints(contact_dist_thresh);
    calc_impacting_unilateral_constraint_forces(-1.0);

restart_with_new_limits:

    // mark velocity estimates as valid
    validate_limit_estimates();

    // determine the maximum step according to conservative advancement
    double safe_dt = std::min(calc_CA_step(), dt);
    FILE_LOG(LOG_SIMULATOR) << "  - conservative advancement step: " << safe_dt << std::endl;

    // get next possible constraint time
    double next_event_time = calc_next_CA_step(contact_dist_thresh);
    FILE_LOG(LOG_SIMULATOR) << "  - *next* conservative advancement step: " << next_event_time << std::endl;

    // we know that the current time is safe, so if the conservative
    // advancement step is zero, set it to the next event time
    if (safe_dt == 0.0)
      safe_dt = std::min(dt, next_event_time);

    // if (the distance between bodies is small [next_event_time will be
    // < INF if there is an event at the current time] and the next time of
    // contact is greater than an Euler step) or (safe_dt is greater than the
    // Euler step), we can *try* generic integration
    if ((next_event_time < INF && next_event_time > euler_step) ||
        safe_dt > euler_step)
    {
      // try the integration
      IntegrationResult istat = integrate_generic(safe_dt);

      // if limits were exceeded during integration, compute the conservative
      // advancement step and redo
      if (istat == eVelocityLimitExceeded)
        goto restart_with_new_limits;
      else if (istat == eMinStepReached)
      {
        // do Euler step; if we're here, we know euler_step < dt
        step_si_Euler(euler_step);
        h += euler_step;

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
      next_event_time = std::min(next_event_time, _coldet->calc_CA_step(pdi));
  }

  if (!found_one)
    return INF;
  else
    return next_event_time;
}

/// Attempts to do "generic" integration (i.e., trying not to pay attention to state violations)
/**
 * \param dt the time step to attempt to integrate, on output the step taken
 */
EventDrivenSimulator::IntegrationResult EventDrivenSimulator::integrate_generic(double& dt)
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
      BOOST_FOREACH(ControlledBodyPtr db, _bodies)
      {
        if (db->limit_estimates_exceeded())
        {
          FILE_LOG(LOG_SIMULATOR) << " ** limit estimates exceeded for body " << db->id << "; retrying with new estimates" << std::endl;

          // reset the state of all bodies
          restore_state();

          return eVelocityLimitExceeded;
        }
      }

      // update the time
      current_time += dt;

      return eIntegrationSuccessful;
    }
    catch (InvalidStateException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** attempted to evaluate derivative at invalid state; halving acceleration step size to " << (dt*0.5) << std::endl;
    }
    catch (InvalidVelocityException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** attempted to evaluate derivative at invalid velocity; halfing acceleration step size to " << (dt*0.5) << std::endl;
    }
    catch (SustainedUnilateralConstraintSolveFailException e)
    {
      FILE_LOG(LOG_SIMULATOR) << " ** failed to solve an LCP; halving step size" << std::endl;
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
  BOOST_FOREACH(ControlledBodyPtr body, _bodies)
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
    _bodies[i]->get_generalized_coordinates(DynamicBodyd::eEuler, _qsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, _qdsave[i]);
  }
}

/// Restores the state of the 'system' (all dynamic bodies)
void EventDrivenSimulator::restore_state()
{
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->set_generalized_coordinates(DynamicBodyd::eEuler, _qsave[i]);
    _bodies[i]->set_generalized_velocity(DynamicBodyd::eSpatial, _qdsave[i]);
  }
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
      throw InvalidVelocityException(t);
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
    double d = _coldet->calc_signed_dist(cg1, cg2, p1, p2);

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

/// Computes a conservative advancement step
double EventDrivenSimulator::calc_CA_step()
{
  // setup safe amount to step
  double dt = std::numeric_limits<double>::max();

  // do joint limit CA step first (it's faster)
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
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

    // if the distance between the two is zero, two are in contact already
    if (pdi.dist < contact_dist_thresh)
      return 0.0;

    // call conservative advancement method
    double step = _coldet->calc_CA_step(pdi);
    dt = std::min(dt, step);
    if (dt <= 0.0)
      return dt;
  }

  return dt;
}

void EventDrivenSimulator::reset_limit_estimates() const
{
  // now compute the bounds
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    // first, reset the limit estimates
    db->reset_limit_estimates();
  }
}

/// Implements Base::load_from_xml()
void EventDrivenSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify node name b/c this is not abstract class
  assert(strcasecmp(node->name.c_str(), "EventDrivenSimulator") == 0);

  // first, load all data specified to the parent object
  TimeSteppingSimulator::load_from_xml(node, id_map);
 
  // read the maximum Euler step
  XMLAttrib* Euler_step_attrib = node->get_attrib("Euler-step");
  if (Euler_step_attrib)
    euler_step = Euler_step_attrib->get_real_value();
}

/// Implements Base::save_to_xml()
void EventDrivenSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call parent's save method first
  TimeSteppingSimulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "EventDrivenSimulator";

  // save the maximum Euler step
  node->attribs.insert(XMLAttrib("Euler-step", euler_step));
}


