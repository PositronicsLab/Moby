/****************************************************************************
 * Copyright 2015 Evan Drumwright
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
#include <Moby/TimeSteppingSimulator.h>

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
TimeSteppingSimulator::TimeSteppingSimulator()
{
}

/// Steps the simulator forward by the given step size
double TimeSteppingSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

  // determine the set of collision geometries
  determine_geometries();

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif

  // clear stored derivatives
  _current_dx.resize(0);

  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << " by " << step_size << std::endl;
  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q, qd;
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      db->get_generalized_coordinates(DynamicBody::eEuler, q);
      db->get_generalized_velocity(DynamicBody::eSpatial, qd);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " Euler coordinates (before): " << q << std::endl;
      FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " spatial velocity (before): " << qd << std::endl;
    }
  }

  // do broad phase collision detection (must be done before any Euler steps)
  broad_phase(step_size);

  // compute pairwise distances at the current configuration
  calc_pairwise_distances();

  // do the Euler step
  step_si_Euler(step_size);

  // update constraint violations
  update_constraint_violations(_pairwise_distances);

  // call the mini-callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn(this);

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
double TimeSteppingSimulator::calc_next_CA_Euler_step(double contact_dist_thresh) const
{
  const double INF = std::numeric_limits<double>::max();
  double next_event_time = INF;

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::calc_next_CA_Euler_step entered" << std::endl; 

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
      double event_time = _coldet->calc_CA_Euler_step(pdi);

      FILE_LOG(LOG_SIMULATOR) << "Next contact time between " << pdi.a->get_single_body()->id << " and " << pdi.b->get_single_body()->id << ": " << event_time << std::endl;

      // not a current event, find when it could become active
      next_event_time = std::min(next_event_time, event_time);
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::calc_next_CA_Euler_step exited" << std::endl; 

  return next_event_time;
}

/// Integrates bodies' velocities forward by dt using Euler integration
void TimeSteppingSimulator::integrate_velocities_Euler(double dt)
{
  VectorNd qd, qdd;

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::integrate_velocities_Euler() entered " << std::endl;

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

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::integrate_velocities_Euler() exited " << std::endl;
}

/// Integrates bodies' positions forward by dt using Euler integration
void TimeSteppingSimulator::integrate_positions_Euler(double dt)
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

/// Does a semi-implicit step
void TimeSteppingSimulator::step_si_Euler(double dt)
{
  FILE_LOG(LOG_SIMULATOR) << "-- doing semi-implicit Euler step" << std::endl;
  const double INF = std::numeric_limits<double>::max();

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
    find_unilateral_constraints(contact_dist_thresh);

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

    // get the time of the next event(s), skipping events at current step
    double h = std::min(calc_next_CA_Euler_step(contact_dist_thresh), target_time - current_time);
    FILE_LOG(LOG_SIMULATOR) << "   position integration: " << h << std::endl;

    // look for small position integration events
    if (h < dt*dt*dt)
    {
      std::cerr << "TimeSteppingSimulator::step_si_Euler() warning: small position integration" << std::endl;
      std::cerr << "  timestep (" << h << ") taken, relative to nominal step (" << dt << ") " << std::endl;
    } 

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
void TimeSteppingSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // do not verify node name b/c this may be a parent class
  // assert(strcasecmp(node->name.c_str(), "TimeSteppingSimulator") == 0);

  // first, load all data specified to the ConstraintSimulator object
  ConstraintSimulator::load_from_xml(node, id_map);
}

/// Implements Base::save_to_xml()
void TimeSteppingSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call ConstraintSimulator's save method 
  ConstraintSimulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "TimeSteppingSimulator";
}


