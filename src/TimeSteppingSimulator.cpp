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
#include <Moby/Dissipation.h>
#include <Moby/ControlledBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
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
  min_step_size = NEAR_ZERO;
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
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(cb);
      db->get_generalized_coordinates(DynamicBodyd::eEuler, q);
      db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " Euler coordinates (before): " << q << std::endl;
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " spatial velocity (before): " << qd << std::endl;
    }
  }

  // do broad phase collision detection (must be done before any Euler steps)
  broad_phase(step_size);

  // compute pairwise distances at the current configuration
  calc_pairwise_distances();

  // do the Euler step
  step_si_Euler(step_size);

  // call the callback
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  std::ofstream out("cvio.dat", std::ostream::app);
  double d = std::numeric_limits<double>::max();
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
    d = std::min(d, _pairwise_distances[i].dist);
  out << d << std::endl;
  out.close();

  // do constraint stabilization
  shared_ptr<ConstraintSimulator> simulator = dynamic_pointer_cast<ConstraintSimulator>(shared_from_this());
  _cstab.stabilize(simulator);

  return step_size;
}

/// Computes impacting unilateral constraint forces, but does *not* do post-processing 
void TimeSteppingSimulator::calc_impacting_unilateral_constraint_forces2(double dt)
{
  // if there are no constraints, quit now
  if (_rigid_constraints.empty() && implicit_joints.empty())
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
}

/// Does a full integration cycle (but not necessarily a full step)
double TimeSteppingSimulator::do_mini_step(double dt)
{
  VectorNd q, qd, qdd;
  std::vector<VectorNd> qsave;

  // init qsave to proper size
  qsave.resize(_bodies.size());

  // save generalized coordinates for all bodies
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->get_generalized_coordinates(DynamicBodyd::eEuler, qsave[i]);
  }

  // set the amount stepped
  double h = 0.0;

  // integrate positions until a new event is detected
  while (h < dt)
  {
    // do broad phase collision detection
    broad_phase(dt-h);

    // compute pairwise distances
    calc_pairwise_distances();

    // get the conservative advancement step
    double tc = std::max(min_step_size, calc_next_CA_Euler_step(contact_dist_thresh));
    FILE_LOG(LOG_SIMULATOR) << "Conservative advancement step: " << tc << std::endl;

    // don't take too large a step
    tc = std::min(dt-h, tc); 

    // integrate the bodies' positions by h + conservative advancement step
    for (unsigned i=0; i< _bodies.size(); i++)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
      db->set_generalized_coordinates(DynamicBodyd::eEuler, qsave[i]);
      db->get_generalized_velocity(DynamicBodyd::eEuler, q);
      q *= (h + tc);
      q += qsave[i];
      db->set_generalized_coordinates(DynamicBodyd::eEuler, q);
    }

    // update h
    h += tc;
  }

  FILE_LOG(LOG_SIMULATOR) << "Position integration ended w/h = " << h << std::endl;

  // prepare to calculate forward dynamics
  precalc_fwd_dyn();

  // apply compliant unilateral constraint forces
  calc_compliant_unilateral_constraint_forces();

  // compute forward dynamics
  calc_fwd_dyn(h);

  // integrate the bodies' velocities forward by h
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->get_generalized_acceleration(qdd);
    qdd *= h;
    db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
    qd += qdd;
    db->set_generalized_velocity(DynamicBodyd::eSpatial, qd);
  }

  // dissipate some energy
  if (_dissipator)
  {
    vector<shared_ptr<DynamicBodyd> > bodies;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
      bodies.push_back(dynamic_pointer_cast<DynamicBodyd>(cb));
    _dissipator->apply(bodies);
  }

  FILE_LOG(LOG_SIMULATOR) << "Integrated velocity by " << h << std::endl;

  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints(contact_dist_thresh);

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(-1.0);

  // update the time
  current_time += h;

  // do a mini-step callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn((ConstraintSimulator*) this);

  return h;
/*
  VectorNd q, qd, qdd, v, vd;
  std::vector<VectorNd> qsave, qdsave, vsave, deltaqd, deltav;
  std::list<sorted_pair<CollisionGeometryPtr> > geom_diff;
  const double MIN_STEP_SIZE = 1e-10;

  // init qsave to proper size
  qsave.resize(_bodies.size());
  vsave.resize(_bodies.size());
  qdsave.resize(_bodies.size());
  deltaqd.resize(_bodies.size());
  deltav.resize(_bodies.size());


  // save all current body configurations and velocities
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->get_generalized_coordinates(DynamicBodyd::eEuler, qsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, vsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBodyd::eEuler, qdsave[i]);
  }

  // integrate acceleration in 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // get the current velocity  
    v = vsave[i]; 
   
    // integrate the acceleration into the velocity with a nominal step of 1.0
    _bodies[i]->get_generalized_acceleration(vd);
    v += vd;
    _bodies[i]->set_generalized_velocity(DynamicBodyd::eSpatial, v);
  }

  // compute impacts (if any), getting the new velocity
  calc_impacting_unilateral_constraint_forces2(-1.0);

  // get the change in velocity 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // get the generalized velocity
    _bodies[i]->get_generalized_velocity(DynamicBodyd::eEuler, deltaqd[i]);
    _bodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, deltav[i]);

    // get the old generalized velocity
    deltaqd[i] -= qdsave[i];
    deltav[i] -= vsave[i];

    // reset the old generalized velocity
    _bodies[i]->set_generalized_velocity(DynamicBodyd::eEuler, qdsave[i]);

    // set the new generalized acceleration
    _bodies[i]->set_generalized_acceleration(deltav[i]);
  }

  // setup amount remaining to step
  double dt_remaining = dt;  

  // get contacts between geometries
  find_unilateral_constraints(contact_dist_thresh);
  std::set<sorted_pair<CollisionGeometryPtr> > contact_geoms = get_current_contact_geoms(); 

  // copy the pairwise distance vector
  vector<PairwiseDistInfo> current_pairwise_distances = _pairwise_distances;

  // setup backtracking constant
  const double BETA = 0.9;
  FILE_LOG(LOG_SIMULATOR) << "About to check whether backtracking necessary using dt=" << dt_remaining << std::endl;

  // determine maximum step
  step_forward(dt_remaining, qsave, qdsave, deltaqd);
  calc_pairwise_distances();
  if (!constraints_met(current_pairwise_distances))
  {
    while (dt_remaining > MIN_STEP_SIZE)
    {
      dt_remaining *= BETA;
      step_forward(dt_remaining, qsave, qdsave, deltaqd);
      calc_pairwise_distances();
      if (constraints_met(current_pairwise_distances))
        break; 
    }
  }
  
  FILE_LOG(LOG_SIMULATOR) << "Found maximum dt=" << dt_remaining << std::endl;

  // setup accumulator for h
  double h_accum = 0.0;

  // get contacts between geometries
  step_forward(0.0, qsave, qdsave, deltaqd);
  calc_pairwise_distances();

  // loop until we have a new contact made
  while (h_accum < dt_remaining)
  {
    // get the time of the next event(s), skipping events at current step
    double h = std::min(calc_next_CA_Euler_step(contact_dist_thresh), dt_remaining);
    FILE_LOG(LOG_SIMULATOR) << "stepping bodies tentatively forward by " << h << std::endl;
    if (h < MIN_STEP_SIZE)
      h = std::min(MIN_STEP_SIZE, dt_remaining);

    // step forward by h 
    step_forward(h_accum+h, qsave, qdsave, deltaqd);

    // recompute pairwise distance
    calc_pairwise_distances();

    // see whether any new contacts were made
    geom_diff.clear();
    find_unilateral_constraints(contact_dist_thresh);
    std::set<sorted_pair<CollisionGeometryPtr> > new_contact_geoms = get_current_contact_geoms(); 
    std::set_difference(new_contact_geoms.begin(), new_contact_geoms.end(), contact_geoms.begin(), contact_geoms.end(), std::back_inserter(geom_diff));
    if (!geom_diff.empty())
    {
      step_forward(h_accum+h, qsave, qdsave, deltaqd);
      calc_pairwise_distances();
      if (LOGGING(LOG_SIMULATOR))
      {
        for (unsigned i=0; i< _pairwise_distances.size(); i++)
          if (_pairwise_distances[i].dist < 0.0)
            FILE_LOG(LOG_SIMULATOR) << "minimum distance: " << _pairwise_distances[i].dist << std::endl;
      }
      FILE_LOG(LOG_SIMULATOR) << "new contact reported: quitting" << std::endl;
      break;
    }

    // update h_accum
    h_accum += h;
    dt_remaining -= h;
  }

  // update rigid constraints
  for (unsigned i=0; i< _rigid_constraints.size(); i++)
  {
    _rigid_constraints[i].contact_impulse *= h_accum;
    _rigid_constraints[i].limit_impulse *= h_accum;
  }

  // call the post application callback, if any
  if (constraint_post_callback_fn)
    (*constraint_post_callback_fn)(_rigid_constraints, constraint_post_callback_data);

  if (LOGGING(LOG_SIMULATOR))
  {
    for (unsigned i=0; i< _pairwise_distances.size(); i++)
      if (_pairwise_distances[i].dist < 0.0)
        FILE_LOG(LOG_SIMULATOR) << "minimum distance: " << _pairwise_distances[i].dist << std::endl;
  }

  return h_accum; 
*/
}

/// Checks to see whether all constraints are met
bool TimeSteppingSimulator::constraints_met(const std::vector<PairwiseDistInfo>& current_pairwise_distances)
{
  // see whether constraints are met to specified tolerance
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
  {
    const PairwiseDistInfo& pdi = _pairwise_distances[i];
    if (current_pairwise_distances[i].dist < 0.0 &&
        pdi.dist < current_pairwise_distances[i].dist - NEAR_ZERO)
    {
      // check whether one of the bodies is compliant
      RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
      RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
      if (rba->compliance == RigidBody::eCompliant || 
          rbb->compliance == RigidBody::eCompliant)
        continue;

      FILE_LOG(LOG_SIMULATOR) << "signed distance between " << rba->id << " and " << rbb->id << "(" << pdi.dist << ") below tolerance: " << (pdi.dist - current_pairwise_distances[i].dist) << std::endl;

      return false;
    }
  }

  return true;
}

/// Gets the current set of contact geometries
std::set<sorted_pair<CollisionGeometryPtr> > TimeSteppingSimulator::get_current_contact_geoms() const
{
  std::set<sorted_pair<CollisionGeometryPtr> > contact_geoms;

  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    contact_geoms.insert(make_sorted_pair(_rigid_constraints[i].contact_geom1, _rigid_constraints[i].contact_geom2));

  return contact_geoms; 
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
    const vector<shared_ptr<Jointd> >& joints = ab->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
    {
      JointPtr joint = dynamic_pointer_cast<Joint>(joints[i]);
      for (unsigned j=0; j< joint->num_dof(); j++)
      {
        if (joint->q[j] < joint->hilimit[j] && joint->qd[j] > 0.0)
        {
          double t = (joint->hilimit[j] - joint->q[j])/joint->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
        if (joint->q[j] > joint->lolimit[j] && joint->qd[j] < 0.0)
        {
          double t = (joint->lolimit[j] - joint->q[j])/joint->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
      }
    }
  }

  // if the distance between any pair of bodies is sufficiently small
  // get next possible event time
  BOOST_FOREACH(PairwiseDistInfo pdi, _pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

      // compute an upper bound on the event time
      double event_time = _coldet->calc_CA_Euler_step(pdi);

      FILE_LOG(LOG_SIMULATOR) << "Next contact time between " << pdi.a->get_single_body()->body_id << " and " << pdi.b->get_single_body()->body_id << ": " << event_time << std::endl;

      // not a current event, find when it could become active
      next_event_time = std::min(next_event_time, event_time);
    }

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::calc_next_CA_Euler_step exited" << std::endl; 

  return next_event_time;
}

/// Does a semi-implicit step
void TimeSteppingSimulator::step_si_Euler(double dt)
{
  FILE_LOG(LOG_SIMULATOR) << "-- doing semi-implicit Euler step" << std::endl;
  const double INF = std::numeric_limits<double>::max();

  // do a number of mini-steps until integrated forward fully
  double h = 0.0;
  while (h < dt)
    h += do_mini_step(dt-h);

  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(cb);
      db->get_generalized_coordinates(DynamicBodyd::eEuler, q);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " position (after integration): " << q << std::endl;
    }
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

  // read the minimum step size
  XMLAttrib* min_step_attrib = node->get_attrib("min-step-size");
  if (min_step_attrib)
    min_step_size = min_step_attrib->get_real_value();
}

/// Implements Base::save_to_xml()
void TimeSteppingSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call ConstraintSimulator's save method 
  ConstraintSimulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "TimeSteppingSimulator";

  // save the minimum step size
  node->attribs.insert(XMLAttrib("min-step-size", min_step_size));
}


