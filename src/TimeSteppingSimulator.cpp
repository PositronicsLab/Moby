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
#include <Moby/RecurrentForce.h>
#include <Moby/DynamicBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/VariableStepIntegrator.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/SustainedUnilateralConstraintSolveFailException.h>
#include <Moby/InvalidStateException.h>
#include <Moby/InvalidVelocityException.h>
#include <Moby/Dissipation.h>
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

  // call the callback
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  std::ofstream out("cvio.dat", std::ostream::app);
  double d = std::numeric_limits<double>::max();
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
    d = std::min(d, _pairwise_distances[i].dist);
  out << d << std::endl;
  out.close();

  return step_size;
}

/// Computes impacting unilateral constraint forces, but does *not* do post-processing 
void TimeSteppingSimulator::calc_impacting_unilateral_constraint_forces2(double dt)
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
}

#define NEW_MINISTEP
#ifdef NEW_MINISTEP

// utility function for do_mini_step(.)
static void reset_positions_and_velocities(const vector<DynamicBodyPtr>& bodies, const vector<VectorNd>& qsave, const vector<VectorNd>& vssave)
{
  for (unsigned i=0; i< bodies.size(); i++)
  {
    bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
    bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vssave[i]);
  }
}

// utility function for do_mini_step(.)
static void integrate_positions(const vector<DynamicBodyPtr>& bodies, double h)
{
  VectorNd q, ve;
 
  for (unsigned i=0; i< bodies.size(); i++)
  {
    // get the generalized coordinates and velocity for the body
    bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q);
    bodies[i]->get_generalized_velocity(DynamicBody::eEuler, ve);
    
    // integrate position forward 
    ve *= h;
    q += ve;
    bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
  }
}

// utility function for do_mini_step(.)
static void integrate_velocities(const vector<DynamicBodyPtr>& bodies, double h)
{
  VectorNd qd, qdd;

   // integrate the bodies' velocities forward by h
   for (unsigned i=0; i< bodies.size(); i++)
   {
     bodies[i]->get_generalized_acceleration(qdd);
     bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
     qdd *= h;
     qd += qdd;
     bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd);
   }
}

// utility function for do_mini_step(.)
static void get_coordinates(const vector<DynamicBodyPtr>& bodies, vector<VectorNd>& q)
{
  for (unsigned i=0; i< bodies.size(); i++)
    bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q[i]);
}

// utility function for do_mini_step(.)
static void get_velocities(const vector<DynamicBodyPtr>& bodies, vector<VectorNd>& v)
{
  for (unsigned i=0; i< bodies.size(); i++)
    bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, v[i]);
}

//#define SEMI_IMPLICIT
// integrates bodies forward by h
void TimeSteppingSimulator::integrate_bodies(const vector<DynamicBodyPtr>& bodies, double h)
{ 
  // compute forward dynamics using old positions if fully explicit 
  #ifndef SEMI_IMPLICIT
  calc_fwd_dyn();
  #endif

  // integrate positions forward by h
  integrate_positions(bodies, h);

  // compute forward dynamics using new positions if semi implicit 
  #ifdef SEMI_IMPLICIT
  calc_fwd_dyn();
  #endif

  // integrate velocities forward by h
  integrate_velocities(bodies, h);

  // Handle any impacts 
  broad_phase(h);
  calc_pairwise_distances();
  find_unilateral_constraints(contact_dist_thresh);
  calc_impacting_unilateral_constraint_forces(-1.0);

  // dissipate energy
  if (_dissipator)
    _dissipator->apply(bodies);
}

// computes the kinetic energy for a rigid body in the global frame
static double calc_KE(RigidBodyPtr rb)
{
  SpatialRBInertiad M0; 
  SVelocityd v0;

  // save the computation frame type
  ReferenceFrameType rftype = rb->get_computation_frame_type();

  // set the computation frame to global
  rb->set_computation_frame_type(Moby::eGlobal);

  // get the rigid body inertia
  M0 = rb->get_inertia();

  // get the rigid body velocity
  v0 = rb->get_velocity();

  // compute the kinetic energy
  double ke = v0.dot(M0 * v0); 

  // reset the reference frame
  rb->set_computation_frame_type(rftype);

  return ke;
}

// stores the state for a rigid body
void TimeSteppingSimulator::store_state(RigidBodyPtr rb)
{
  VectorNd q, v;
  rb->get_generalized_coordinates(DynamicBody::eEuler, q);
  rb->get_generalized_velocity(DynamicBody::eSpatial, v);
  rb->abs_pos_err = q;
  rb->abs_vel_err = v;
}

// stores the state for a dynamic body
void TimeSteppingSimulator::store_state(DynamicBodyPtr db)
{
  ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
  if (ab)
  {
    for (unsigned i=0; i< ab->get_links().size(); i++)
      store_state(ab->get_links()[i]);
    return;
  }

  // it's a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
  store_state(rb);
}

// updates the state for a rigid body
void TimeSteppingSimulator::update_state(RigidBodyPtr rb, double h)
{
  VectorNd q, v;
  rb->get_generalized_coordinates(DynamicBody::eEuler, q);
  rb->get_generalized_velocity(DynamicBody::eSpatial, v);
  (rb->abs_pos_err -= q) /= h;
  (rb->abs_vel_err -= v) /= h;
}

// updates the state for a dynamic body
void TimeSteppingSimulator::update_state(DynamicBodyPtr db, double h)
{
  ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
  if (ab)
  {
    for (unsigned i=0; i< ab->get_links().size(); i++)
      update_state(ab->get_links()[i], h);
    return;
  }

  // it's a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
  update_state(rb, h);
}

// computes the kinetic energy for a body
static double calc_KE(DynamicBodyPtr body)
{
  ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);
  if (ab)
  {
    double T = 0.0;
    for (unsigned i=0; i< ab->get_links().size(); i++)
      T += calc_KE(ab->get_links()[i]);
    return T;
  }

  // it's a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(body);
  return calc_KE(rb);
}

// records the error for a rigid body
double TimeSteppingSimulator::record_error(RigidBodyPtr rb)
{
  double min_k = std::numeric_limits<double>::max();

  for(int j=0;j<6;j++){
    // record rel err
    // if rel err tol exceeded
    double k = 1.0;
    // Calc step scale
    if(std::fabs(rb->abs_pos_err[j]) > NEAR_ZERO)
      k = rb->abs_pos_err_tol[j]/std::fabs(rb->abs_pos_err[j]);
    min_k = std::min(min_k, k);
    if(std::fabs(rb->abs_vel_err[j]) > NEAR_ZERO)
      k = rb->abs_vel_err_tol[j]/std::fabs(rb->abs_vel_err[j]);
    min_k = std::min(min_k, k);
  }
  std::cerr << "absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
  //std::cerr << "relative position error (" <<  rb->id << ") " << rb->rel_pos_err << std::endl;
  std::cerr << "absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
  //std::cerr << "relative velocity error (" <<  rb->id << ") " << rb->rel_vel_err << std::endl;
  return min_k;
}

// records the error for a dynamic body
double TimeSteppingSimulator::record_error(DynamicBodyPtr db)
{
  ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(db);
  if (ab)
  {
    double min_k = std::numeric_limits<double>::max();
    for (unsigned i=0; i< ab->get_links().size(); i++)
      min_k = std::min(record_error(ab->get_links()[i]), min_k);
    return min_k;
  }

  // it's a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(db);
  return record_error(rb);
}

// computes integration error
bool TimeSteppingSimulator::calc_integration_error(const vector<DynamicBodyPtr>& bodies, double h, const vector<VectorNd>& qe_large, const vector<VectorNd>& v_large, const vector<VectorNd>& qe_small, const vector<VectorNd>& v_small, double& min_k , std::string flag)
{
  // determine if this integration is successful
  // NOTE: Relative error sucks
  bool error_exceeded = false;
  VectorNd energy_error(bodies.size());
  VectorNd energy_relative_error(bodies.size());
  for (unsigned i=0; i< bodies.size(); i++)
  {
    // set the coordinates and velocities
    reset_positions_and_velocities(bodies, qe_large, v_large);

    // calculate kinetic energy for the large step
    double Tlarge = calc_KE(bodies[i]);

    // store position and velocities
    store_state(bodies[i]);

    // set the coordinates and velocities
    reset_positions_and_velocities(bodies, qe_small, v_small);

    // calculate kinetic energy for the small step
    double Tsmall = calc_KE(bodies[i]);

    // use stored position and velocities to estimate integration errors
    update_state(bodies[i], h);

    energy_error[i] = Tlarge-Tsmall;
    energy_relative_error[i] = energy_error[i]/Tsmall;
    //std::cerr << "Velocity -- 1 steps (" <<  bodies[i]->id << ") " << v1_g << std::endl;
    //std::cerr << "Velocity -- 2 steps (" <<  bodies[i]->id << ") " << v2_g << std::endl;
    //std::cerr << "Inertia -- 1 steps (" <<  bodies[i]->id << ") " << M1_g << std::endl;
    //std::cerr << "Inertia -- 2 steps (" <<  bodies[i]->id << ") " << M2_g << std::endl;

    std::cerr << flag << "Energy -- 1 step (" <<  bodies[i]->id << ") " << Tlarge << std::endl;
    std::cerr << flag << "Energy -- 2 steps (" <<  bodies[i]->id << ") " << Tsmall << std::endl;
    std::cerr << flag << "Energy error (" <<  bodies[i]->id << ") " <<  energy_error[i] << std::endl;
    std::cerr << flag << "Energy relative error (" <<  bodies[i]->id << ") " <<  energy_relative_error[i] << std::endl;

    // record the error for the body
    min_k = record_error(bodies[i]);

    const double LOW_ENERGY = 1.0e-2;
    const double ENERGY_ABS_TOL = 1.0e-2;
    // NOTE: LINKS TOLERANCE 
     const double ENERGY_REL_TOL = 5e-4;
    // NOTE: BOX TOLERANCE
    //const double ENERGY_REL_TOL = 1e-2;
    if(Tsmall > LOW_ENERGY && Tlarge > LOW_ENERGY ){
      if(energy_relative_error[i] > ENERGY_REL_TOL){
        // Check if this integration step is accurate enough
        std::cerr << flag << "Limited by relative error: " << energy_relative_error[i] << " of " << Tsmall << std::endl;
        error_exceeded = true;  
      }
    } else if (std::fabs(energy_error[i]) > ENERGY_ABS_TOL){
        std::cerr << flag << "Limited by absolute error: " << energy_error[i] << " of " << Tsmall << std::endl;
      error_exceeded = true;
    }

  }
  return error_exceeded;
}

/// Updates position and velocity using Richardson extrapolation
static void do_richardson(const vector<DynamicBodyPtr>& bodies, const vector<VectorNd>& qe_small, const vector<VectorNd>& qe_large, const vector<VectorNd>& vs_small, const vector<VectorNd>& vs_large)
{
  VectorNd q, v;

  for (unsigned i=0; i< bodies.size(); i++)
  {
    // Position
    (q = qe_small[i]) *= 2.0;
    q -= qe_large[i];
    bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);

    // Velocity
    (v = vs_small[i]) *= 2.0;
    v -= vs_large[i];
    bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, v);
  }
}

// get active bodies
static void get_active_bodies(const vector<DynamicBodyPtr>& _bodies, vector<DynamicBodyPtr>& bodies)
{
  bodies.clear();
  for (unsigned i=0; i< _bodies.size(); i++){
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(_bodies[i]);
    if(rb)
      if(!rb->is_enabled())
        continue;
    bodies.push_back(_bodies[i]);
  }
}

/// Does a full integration cycle (but not necessarily a full step)
double TimeSteppingSimulator::do_mini_step(double dt)
{
  // only process dynamic bodies
  static std::vector<DynamicBodyPtr> bodies;
  get_active_bodies(_bodies, bodies);
  
  std::cerr << "begin do_mini_step() " << current_time << std::endl;
  VectorNd q, qd, qdd, eps;
  std::vector<VectorNd> qsave, vssave, vesave;

  // init qsave to proper size
  qsave.resize(bodies.size());
  vssave.resize(bodies.size());
  vesave.resize(bodies.size());

  // save generalized coordinates and velocities for all bodies
  for (unsigned i=0; i< bodies.size(); i++){
    bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
    bodies[i]->get_generalized_velocity(DynamicBody::eEuler, vesave[i]);
    bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, vssave[i]);
  }

  // Set rigid body accuracy tolerances
  const double INF = std::numeric_limits<double>::max();
  static SVector6d lin_ang_pos_tol(0.01,0.01,0.01,0.1,0.1,0.1);
  static SVector6d lin_ang_vel_tol(INF,INF,INF,INF,INF,INF);
  
  // see whether to initialize tolerances
  static bool tols_inited = false;
  if(!tols_inited){
    tols_inited = true;

    // init tolerances
    for (unsigned i=0; i< bodies.size(); i++){
      RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(bodies[i]);
      ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
      if(rb){
        lin_ang_pos_tol.to_vector(rb->abs_pos_err_tol ) ;
        lin_ang_vel_tol.to_vector(rb->abs_vel_err_tol ) ;
      } else if (ab) {
        BOOST_FOREACH(rb, ab->get_links()){
          lin_ang_pos_tol.to_vector(rb->abs_pos_err_tol);
          lin_ang_vel_tol.to_vector(rb->abs_vel_err_tol);
        }
      }
    }
  }
  
  // set the amount stepped
  double h = dt;
 
  // TODO: read this in properly
  double min_CA_step_size = 1.0;

  // get the conservative advancement step
  //double tc = std::max(min_CA_step_size, calc_next_CA_Euler_step(contact_dist_thresh));
  //std::cerr << "Conservative advancement step (initial): " << tc << std::endl;

  //h = std::min(h,tc);

  // variables for coordinates and velocities after big and small steps 
  std::vector<VectorNd> qe_large, qs_large, v_large;
  std::vector<VectorNd> qe_small, qs_small, v_small;
  qs_large.resize(bodies.size());
  qs_small.resize(bodies.size());
  qe_large.resize(bodies.size());
  qe_small.resize(bodies.size());
  v_large.resize(bodies.size());
  v_small.resize(bodies.size());
  
  while (true) { // Loop to collisison free time
    while (true) { // loop to accurate integration
      //////////////////////////////////////
      //// INTEGRATION one large and two small steps 
      //////////////////////////////////////
      
      std::cerr << "Finding maximum safe (accurate) integration step size: h =" << h << std::endl;
 
    
      //////////////////////////////////////
      //// INTEGRATION: one large step 
      // integrate the bodies' positions by h + conservative advancement step
      std::cerr << "one large step: h = " << h << std::endl;

      // reset positions and velocities first
      reset_positions_and_velocities(bodies, qsave, vssave);

      // integrate bodies forward by h
      integrate_bodies(bodies, h);

      // store the new coordinates and velocities
      get_coordinates(bodies, qe_large);
      get_velocities(bodies, v_large);

      //////////////////////////////////////
      //// INTEGRATION: two small steps 

      std::cerr << "two small steps: h/2 = " << h/2 << std::endl;
      // reset positions and velocities first
      reset_positions_and_velocities(bodies, qsave, vssave);

      // integrate bodies by h/2 twice 
      integrate_bodies(bodies, h*0.5);
      integrate_bodies(bodies, h*0.5);

      // store the new coordinates and velocities
      get_coordinates(bodies, qe_small);
      get_velocities(bodies, v_small);

      //////////////////////////////////////
      ////  CHECK FOR INTEGRATION ACCURACY  
      //////////////////////////////////////
      
      // Set by absolute error calculation, determines the size of the smaller step
      double min_k;
      std::cerr << "calc_integration_error: h = " << h << std::endl;
      bool error_exceeded = calc_integration_error(bodies, h, qe_large, v_large, qe_small, v_small, min_k); 
      
      // If we exceed relative error bounds
      if(error_exceeded && h > min_step_size){
        std::cerr << "Error too high at this step size!: " << h << std::endl;
        double new_h = 0.9*h*min_k;
        if(new_h > h)
          h = 0.9*h;
        else
          h = new_h;
        
        if(h < min_step_size)
          h = min_step_size;
        
        continue;
      } 
#if 0
      for (unsigned i=0; i< bodies.size(); i++)
      {
        std::cerr << "Final Energy error (" <<  bodies[i]->id << ") " <<  energy_error[i] << std::endl;
        std::cerr << "Final Energy relative error (" <<  bodies[i]->id << ") " <<  energy_relative_error[i] << std::endl;
        ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
        RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(bodies[i]);
        // if Articulated body, make the rigid body the base
        if(ab)
          rb = ab->get_base_link();
        if(rb){
          std::cerr << "Final absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
          std::cerr << "Final absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
        } else if (ab) {
          BOOST_FOREACH(RigidBodyPtr rb, ab->get_links()){
            std::cerr << "Final absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
            std::cerr << "Final absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
          }
        }
      }
#endif

      break;
    }
    std::cerr << "Found maximum safe (accurate) integration step size (or at min-step-size): h =" << h << std::endl;
    std::cerr << "Check if a collision occurs in this interval given error bounds" << std::endl;

    double min_k;
    calc_integration_error(bodies, h, qe_large, v_large, qe_small, v_small, min_k,"Final "); 
    // Set third order optimal position (Richardson extrapolation)
    // for conservative advancement (do not update velocity)
    do_richardson(bodies, qe_small, qe_large, v_small, v_large);
    break;
  }

// TODO: let's disable this for now since we're not using for experiment
/*
    // recompute conservative step 
    broad_phase(h);
    calc_pairwise_distances();
    double tc = min_step_size;
    //try {
    tc = std::max(min_CA_step_size, calc_next_CA_Euler_step(contact_dist_thresh));
    std::cerr << "Conservative advancement step (with error bounds): " << tc << std::endl;
    //} catch (...){}

    // If there is a possible impact in the interval.
    // set step size to time of collision (conservative estimate)
    // then check accuracy over that interval
    if(tc < h)
    {
      std::cerr << "Conservative advancement predicts collision in interval h: " << h << std::endl;
      h = tc;
      // And then restart process
      continue;
    } 
    // Else use current integration values to accurately integrate state
    // Then pass mini-step function on to rest of stepping function
    else {
      // Set third order optimal velocity (Richardson extrapolation)
      for (unsigned i=0; i< bodies.size(); i++)
      {
        // Velocity
        qd = vsave_small[i];
        qd *= 2;
        qd -= vsave_large[i];
        bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd);
      }
      // End step size adjustment
      break;
    }
*/

  std::cerr << "Safe integration ended w/ h = " << h << std::endl;

  // update the time
  current_time += h;

  // do a mini-step callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn((ConstraintSimulator*) this);

  std::cerr << "end do_mini_step() " << current_time << std::endl;
  return h;
}
#else
/// Does a full integration cycle (but not necessarily a full step)
double TimeSteppingSimulator::do_mini_step(double dt)
{
  VectorNd q, qd, qdd;
  std::vector<VectorNd> qsave;

  // init qsave to proper size
  qsave.resize(_bodies.size());

  // save generalized coordinates for all bodies
  for (unsigned i=0; i< _bodies.size(); i++)
    _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, qsave[i]);

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
      _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
      _bodies[i]->get_generalized_velocity(DynamicBody::eEuler, q);
      q *= (h + tc);
      q += qsave[i];
       _bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
    }

    // update h
    h += tc;
  }

  FILE_LOG(LOG_SIMULATOR) << "Position integration ended w/h = " << h << std::endl;

  // compute forward dynamics
  calc_fwd_dyn();

  // integrate the bodies' velocities forward by h
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    _bodies[i]->get_generalized_acceleration(qdd);
    qdd *= h;
    _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
    qd += qdd;
    _bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd);
  }

  FILE_LOG(LOG_SIMULATOR) << "Integrated velocity by " << h << std::endl;

  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints(contact_dist_thresh);

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(-1.0);

  // dissipate some energy
  if (_dissipator)
    _dissipator->apply(_bodies);

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
    _bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, vsave[i]);
    _bodies[i]->get_generalized_velocity(DynamicBody::eEuler, qdsave[i]);
  }

  // integrate acceleration in 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // get the current velocity  
    v = vsave[i]; 
   
    // integrate the acceleration into the velocity with a nominal step of 1.0
    _bodies[i]->get_generalized_acceleration(vd);
    v += vd;
    _bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, v);
  }

  // compute impacts (if any), getting the new velocity
  calc_impacting_unilateral_constraint_forces2(-1.0);

  // get the change in velocity 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // get the generalized velocity
    _bodies[i]->get_generalized_velocity(DynamicBody::eEuler, deltaqd[i]);
    _bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, deltav[i]);

    // get the old generalized velocity
    deltaqd[i] -= qdsave[i];
    deltav[i] -= vsave[i];

    // reset the old generalized velocity
    _bodies[i]->set_generalized_velocity(DynamicBody::eEuler, qdsave[i]);

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
#endif

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

      FILE_LOG(LOG_SIMULATOR) << "Next contact time between " << pdi.a->get_single_body()->id << " and " << pdi.b->get_single_body()->id << ": " << event_time << std::endl;

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
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      db->get_generalized_coordinates(DynamicBody::eEuler, q);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " position (after integration): " << q << std::endl;
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


