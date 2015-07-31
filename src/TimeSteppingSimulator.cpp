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
/// Does a full integration cycle (but not necessarily a full step)
double TimeSteppingSimulator::do_mini_step(double dt)
{
  // only process dynamic bodies
  static std::vector<DynamicBodyPtr> bodies;
  if(bodies.empty()){
    for (unsigned i=0; i< _bodies.size(); i++){
      RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(_bodies[i]);
      if(rb)
        if(!rb->is_enabled())
          continue;
      bodies.push_back(_bodies[i]);
    }
  }
  FILE_LOG(LOG_SIMULATOR) << "begin do_mini_step() " << current_time << std::endl;
  const double INF = std::numeric_limits<double>::max();
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
 
  // see whether to initialize tolerances
  static bool tols_inited = false;
  static SVector6d lin_ang_pos_tol(INF,INF,INF,INF,INF,INF);
  static SVector6d lin_ang_vel_tol(INF,INF,INF,INF,INF,INF);
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
 
  // EMD: why are these two calls disabled?
  // initially do CA calculation to save some time
  // do broad phase collision detection
  //broad_phase(h);

  // compute pairwise distances
  //calc_pairwise_distances();

  // TODO: read this in properly
  double min_CA_step_size = 0.005;

  // get the conservative advancement step
  //double tc = std::max(min_CA_step_size, calc_next_CA_Euler_step(contact_dist_thresh));
  //FILE_LOG(LOG_SIMULATOR) << "Conservative advancement step (initial): " << tc << std::endl;

  //h = std::min(h,tc);

  //   
  std::vector<VectorNd> qesave_large, qssave_large, vsave_large;
  std::vector<VectorNd> qesave_small, qssave_small, vsave_small;
  qssave_large.resize(bodies.size());
  qssave_small.resize(bodies.size());
  qesave_large.resize(bodies.size());
  qesave_small.resize(bodies.size());
  vsave_large.resize(bodies.size());
  vsave_small.resize(bodies.size());
  
  while (true) { // Loop to collisison free time
    while (true) { // loop to accurate integration
      //////////////////////////////////////
      //// INTEGRATION one large and two small steps 
      //////////////////////////////////////
      
      FILE_LOG(LOG_SIMULATOR) << "Finding maximum safe (accurate) integration step size: h =" << h << std::endl;
      
      //////////////////////////////////////
      //// INTEGRATION: one large step 
      // integrate the bodies' positions by h + conservative advancement step
      for (unsigned i=0; i< bodies.size(); i++)
      {
        bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
        bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vssave[i]);
        (q = vesave[i]) *= (h);
        q += qsave[i];
        bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
        // EMD: below line is dangerous b/c of likely rpy bugs
        bodies[i]->get_generalized_coordinates(DynamicBody::eSpatial, qssave_large[i]);
        qesave_large[i] = q;
      }

      // compute forward dynamics using new positions
      calc_fwd_dyn();

      // recompute pairwise distances
      calc_pairwise_distances();

      // find unilateral constraints
      find_unilateral_constraints(contact_dist_thresh);

      // handle any impacts
      calc_impacting_unilateral_constraint_forces(-1.0);

      // dissipate some energy
      if (_dissipator)
        _dissipator->apply(bodies);
  
      // integrate the bodies' velocities forward by h
      for (unsigned i=0; i< bodies.size(); i++)
      {
        bodies[i]->get_generalized_acceleration(qdd);
        bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
        FILE_LOG(LOG_SIMULATOR) << "qd1" << qd << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "qdd1" << qdd << std::endl;
        qdd *= h;
        qd += qdd;
        vsave_large[i] = qd;
      }
 
      //////////////////////////////////////
      //// INTEGRATION: 2 Small steps 
      // 1st
      for (unsigned i=0; i< bodies.size(); i++)
      {
        bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qsave[i]);
        bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vssave[i]);
        (q = vesave[i]) *= (h/2.0);
        q += qsave[i];
        bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
      }
      // compute forward dynamics
      calc_fwd_dyn();
  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints(contact_dist_thresh);

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(-1.0);

  // dissipate some energy
  if (_dissipator)
    _dissipator->apply(bodies);
  
      // integrate the bodies' velocities forward by h
      for (unsigned i=0; i< bodies.size(); i++)
      {
        // SI Euler ingegration for step 1 velocity
        bodies[i]->get_generalized_acceleration(qdd);
        bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
        qdd *= (h/2.0);
        qd += qdd;
        bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd);
        FILE_LOG(LOG_SIMULATOR) << "qd2" << qd << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "qdd2" << qdd << std::endl;
        bodies[i]->get_generalized_velocity(DynamicBody::eEuler, qd);
        // Integration position
        bodies[i]->get_generalized_coordinates(DynamicBody::eEuler, q);
        qd *= (h/2.0);
        q += qd;
        bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
        bodies[i]->get_generalized_coordinates(DynamicBody::eSpatial, qssave_small[i]);
        qesave_small[i] = q;
      }
      // compute forward dynamics
      calc_fwd_dyn();
  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints(contact_dist_thresh);

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(-1.0);

  // dissipate some energy
  if (_dissipator)
    _dissipator->apply(bodies);
  
      // integrate the bodies' velocities forward by h
      for (unsigned i=0; i< bodies.size(); i++)
      {
        // SI Euler ingegration for step 2 velocity
        bodies[i]->get_generalized_acceleration(qdd);
        qdd *= (h/2.0);
        bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd);
        qd += qdd;
        vsave_small[i] = qd;
        FILE_LOG(LOG_SIMULATOR) << "qd3" << qd << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "qdd3" << qdd << std::endl;
      }

      //////////////////////////////////////
      ////  CHECK FOR INTEGRATION ACCURACY  
      //////////////////////////////////////
      
      // Set by absolute error calculation, determines the size of the smaller step
      double min_k = INF;
      // determine if this integration is successful
      // NOTE: Relative error sucks
      double error_exceeded = false;
      VectorNd energy_error(bodies.size());
      VectorNd energy_relative_error(bodies.size());
      for (unsigned i=0; i< bodies.size(); i++)
      {
        VectorNd momentum, v1_g,v2_g, v0;
        MatrixNd M1_g(6,6), M2_g(6,6),I0;
        int N = M1_g.rows();
        ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
        RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(bodies[i]);
        // if Articulated body, make the rigid body the base
        if(ab)
          rb = ab->get_base_link();
        ReferenceFrameType rftype = rb->get_computation_frame_type();
        rb->set_computation_frame_type(Moby::eGlobal);
        {
          bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qesave_large[i]);
          bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vsave_large[i]);
          
          // Get the articulted body inertia
          if(ab)
            bodies[i]->get_generalized_inertia(M1_g);
         
          N = M1_g.rows();
         
          // Correct the velocity to global
          v1_g = vsave_large[i];
          const SVelocityd& _xd0 = rb->get_velocity();
          _xd0.transpose_to_vector(v0);
          v1_g.segment(N-6,N) = v0;
          
          // Correct the inertial to global
          const SpatialRBInertiad& _J0 = rb->get_inertia();
          _J0.to_matrix(I0);
          M1_g.block(N-6,N,N-6,N) = I0;
        
          // While we're here, record the errors into the rigid body structure
          rb->abs_pos_err = (qssave_large[i]);
          rb->abs_vel_err = (v0);
          if (ab) {
            BOOST_FOREACH(RigidBodyPtr rb, ab->get_links()){
              if(rb->is_base())
                continue;
              ReferenceFrameType rftype = rb->get_computation_frame_type();
              rb->set_computation_frame_type(Moby::eGlobal);
              const SVelocityd& _xd0 = rb->get_velocity();
              _xd0.transpose_to_vector(v0);
              rb->abs_vel_err = (v0);
              
              rb->abs_pos_err = (qssave_large[i]);
              rb->set_computation_frame_type(rftype);
            }
          }
        }
        
        // Calc energy for this system
        double energy1 = M1_g.mult(v1_g,momentum).dot(v1_g);
        
        {
          // --- Now do the same for the 2 step case ---
          bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, qesave_small[i]);
          bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vsave_small[i]);
          
          // Get the articulted body inertia
          if(ab)
            bodies[i]->get_generalized_inertia(M2_g);
         
          // Correct the velocity to global
          v2_g = vsave_small[i];
          const SVelocityd& _xd0 = rb->get_velocity();
          rb->_xd0.transpose_to_vector(v0);
          v2_g.segment(N-6,N) = v0;
          
          // Correct the inertial to global
          const SpatialRBInertiad& _J0 = rb->get_inertia();
          _J0.to_matrix(I0);
          M2_g.block(N-6,N,N-6,N) = I0;
          // While we're here, record the errors into the rigid body structure
          (rb->abs_pos_err -= qssave_small[i]) /=h;
          (rb->abs_vel_err -= v0) /=h;
          if (ab) {
            BOOST_FOREACH(RigidBodyPtr rb, ab->get_links()){
              if(rb->is_base())
                continue;
              ReferenceFrameType rftype = rb->get_computation_frame_type();
              rb->set_computation_frame_type(Moby::eGlobal);
              const SVelocityd& _xd0 = rb->get_velocity();
              _xd0.transpose_to_vector(v0);
              (rb->abs_vel_err -= v0) /=h;
              
              (rb->abs_pos_err -= qssave_small[i]) /=h;
              rb->set_computation_frame_type(rftype);
            }
          }
        } 
        rb->set_computation_frame_type(rftype);
        
        // Calc energy for this system
        double energy2 = M2_g.mult(v2_g,momentum).dot(v2_g);

        energy_error[i] = energy1-energy2;
        energy_relative_error[i] = energy_error[i]/energy2;
        //FILE_LOG(LOG_SIMULATOR) << "Velocity -- 1 steps (" <<  bodies[i]->id << ") " << v1_g << std::endl;
        //FILE_LOG(LOG_SIMULATOR) << "Velocity -- 2 steps (" <<  bodies[i]->id << ") " << v2_g << std::endl;
        //FILE_LOG(LOG_SIMULATOR) << "Inertia -- 1 steps (" <<  bodies[i]->id << ") " << M1_g << std::endl;
        //FILE_LOG(LOG_SIMULATOR) << "Inertia -- 2 steps (" <<  bodies[i]->id << ") " << M2_g << std::endl;

        FILE_LOG(LOG_SIMULATOR) << "Energy -- 1 step (" <<  bodies[i]->id << ") " << energy1 << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "Energy -- 2 steps (" <<  bodies[i]->id << ") " << energy2 << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "Energy error (" <<  bodies[i]->id << ") " <<  energy_error[i] << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "Energy relative error (" <<  bodies[i]->id << ") " <<  energy_relative_error[i] << std::endl;

        // Position
        ((bodies[i]->abs_pos_err = qssave_small[i]) -= qssave_large[i]) /= h;
        // Velocity
        ((bodies[i]->abs_vel_err = vsave_small[i]) -= vsave_large[i]) /= h;
        
        if(rb){
          for(int j=0;j<6;j++){
            // record rel err
            // if rel err tol exceeded
            double k = 1.0;
            // Calc step scale
            if(fabs(rb->abs_pos_err[j]) > NEAR_ZERO)
              k = rb->abs_pos_err_tol[j]/fabs(rb->abs_pos_err[j]);
            min_k = std::min(min_k, k);
            if(fabs(rb->abs_vel_err[j]) > NEAR_ZERO)
              k = rb->abs_vel_err_tol[j]/fabs(rb->abs_vel_err[j]);
            min_k = std::min(min_k, k);
          }
          FILE_LOG(LOG_SIMULATOR) << "absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
          //FILE_LOG(LOG_SIMULATOR) << "relative position error (" <<  rb->id << ") " << rb->rel_pos_err << std::endl;
          FILE_LOG(LOG_SIMULATOR) << "absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
          //FILE_LOG(LOG_SIMULATOR) << "relative velocity error (" <<  rb->id << ") " << rb->rel_vel_err << std::endl;
        } else if (ab) {
          BOOST_FOREACH(RigidBodyPtr rb, ab->get_links()){
            if(rb->is_base())
              continue;
            for(int j=0;j<6;j++){
              // record rel err
              // if rel err tol exceeded
              double k = 1.0;
              // Calc step scale
              if(fabs(rb->abs_pos_err[j]) > NEAR_ZERO)
                k = rb->abs_pos_err_tol[j]/fabs(rb->abs_pos_err[j]);
              min_k = std::min(min_k, k);
              if(fabs(rb->abs_vel_err[j]) > NEAR_ZERO)
                k = rb->abs_vel_err_tol[j]/fabs(rb->abs_vel_err[j]);
              min_k = std::min(min_k, k);
            }
            FILE_LOG(LOG_SIMULATOR) << "absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
            //FILE_LOG(LOG_SIMULATOR) << "relative position error (" <<  rb->id << ") " << rb->rel_pos_err << std::endl;
            FILE_LOG(LOG_SIMULATOR) << "absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
            //FILE_LOG(LOG_SIMULATOR) << "relative velocity error (" <<  rb->id << ") " << rb->rel_vel_err << std::endl;
          }
        }

        const double LOW_ENERGY = 1.0e-4;
        const double ENERGY_ABS_TOL = 1.0e-4;
        const double ENERGY_REL_TOL = 1.0e-2;
        if(energy1 > LOW_ENERGY || energy2 > LOW_ENERGY){
          if(fabs(energy_relative_error[i]) > ENERGY_REL_TOL){
            // Check if this integration step is accurate enough
            error_exceeded = true;  
          }
        } else if (fabs(energy_error[i]) > ENERGY_ABS_TOL){
          error_exceeded = true;
        }

      }
      
      // If we exceed relative error bounds
      if(error_exceeded && h > min_step_size){
        FILE_LOG(LOG_SIMULATOR) << "Error too high at this step size!: " << h << std::endl;
        double new_h = 0.9*h*min_k;
        if(new_h > h)
          h = 0.9*h;
        else
          h = new_h;
        
        if(h < min_step_size)
          h = min_step_size;
        
        continue;
      }
#ifndef NDEBUG
      for (unsigned i=0; i< bodies.size(); i++)
      {
        FILE_LOG(LOG_SIMULATOR) << "Final Energy error (" <<  bodies[i]->id << ") " <<  energy_error[i] << std::endl;
        FILE_LOG(LOG_SIMULATOR) << "Final Energy relative error (" <<  bodies[i]->id << ") " <<  energy_relative_error[i] << std::endl;
        ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
        RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(bodies[i]);
        // if Articulated body, make the rigid body the base
        if(ab)
          rb = ab->get_base_link();
        if(rb){
          FILE_LOG(LOG_SIMULATOR) << "Final absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
          FILE_LOG(LOG_SIMULATOR) << "Final absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
        } else if (ab) {
          BOOST_FOREACH(RigidBodyPtr rb, ab->get_links()){
            FILE_LOG(LOG_SIMULATOR) << "Final absolute position error (" <<  rb->id << ") " << rb->abs_pos_err << std::endl;
            FILE_LOG(LOG_SIMULATOR) << "Final absolute velocity error (" <<  rb->id << ") " << rb->abs_vel_err << std::endl;
          }
        }
      }
#endif
      break;
    }
    FILE_LOG(LOG_SIMULATOR) << "Found maximum safe (accurate) integration step size (or at min-step-size): h =" << h << std::endl;
    FILE_LOG(LOG_SIMULATOR) << "Check if a collision occurs in this interval given error bounds" << std::endl;

    // Set third order optimal position (Richardson extrapolation)
    // for conservative advancement (do not update velocity)
    for (unsigned i=0; i< bodies.size(); i++)
    {
      // Position
      q = qesave_small[i];
      q *= 2;
      q -= qesave_large[i];
      int N = q.rows();
      
      Quatd q1(qesave_small[i][N-4],qesave_small[i][N-3],qesave_small[i][N-2],qesave_small[i][N-1]), 
            q2(qesave_large[i][N-4],qesave_large[i][N-3],qesave_large[i][N-2],qesave_large[i][N-1]);
      q1.slerp(q2,0.5);
      q1.normalize();
      q[N-4] = q1.x;
      q[N-3] = q1.y;
      q[N-2] = q1.z;
      q[N-1] = q1.w;

      bodies[i]->set_generalized_coordinates(DynamicBody::eEuler, q);
      // Velocity
      bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, vssave[i]);
    }

    // Check if collisions occur in the interval
    // do broad phase collision detection
    broad_phase(h);

    // compute pairwise distances
    calc_pairwise_distances();

    // get the conservative advancement step
    // NOTE: This calulation is un necessary for our experiments, 
    // we can be _less_ conservative with contact timing if we want
    //double tc = std::max(min_step_size, calc_next_CA_Euler_step_error(contact_dist_thresh, epsq, epsv));
    double tc = min_step_size;
    //try {
    tc = std::max(min_CA_step_size, calc_next_CA_Euler_step(contact_dist_thresh));
    FILE_LOG(LOG_SIMULATOR) << "Conservative advancement step (with error bounds): " << tc << std::endl;
    //} catch (...){}
    // If there is a possible collision in the interval.
    // set step size to time of collision (conservative estimate)
    // then check accuracy over that interval
    if(tc < h){
      FILE_LOG(LOG_SIMULATOR) << "Conservative advancement predicts collision in interval h: " << h << std::endl;
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
  }

  FILE_LOG(LOG_SIMULATOR) << "Safe integration ended w/ h = " << h << std::endl;

  // update the time
  current_time += h;

  // do a mini-step callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn((ConstraintSimulator*) this);

  FILE_LOG(LOG_SIMULATOR) << "end do_mini_step() " << current_time << std::endl;
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


