/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Integrates both position and velocity of rigid _bodies
/**
 * \return the size of step taken
 */
template <class ForwardIterator>
double Simulator::integrate(double step_size, ForwardIterator begin, ForwardIterator end)
{
  // begin timing dynamics
  tms start;  
  times(&start);

  // get the state-derivative for each dynamic body
  for (ForwardIterator i = begin; i != end; i++)
  {
    // integrate the body
    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (before): " << (*i)->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (before): " << (*i)->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
      (*i)->calc_fwd_dyn(0.0);
      FILE_LOG(LOG_SIMULATOR) << "  generalized acceleration: " << (*i)->get_generalized_acceleration(q) << std::endl;
      RigidBodyPtr rb = boost::dynamic_pointer_cast<RigidBody>(*i);
      Ravelin::SVelocityd sv = Ravelin::Pose3d::transform(rb->get_pose(), rb->get_velocity());
      FILE_LOG(LOG_SIMULATOR) << "  spatial velocity (before, body frame): " << sv << std::endl;
    }
    (*i)->integrate(current_time, step_size, integrator);
    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (after): " << (*i)->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (after): " << (*i)->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
      RigidBodyPtr rb = boost::dynamic_pointer_cast<RigidBody>(*i);
      Ravelin::SVelocityd sv = Ravelin::Pose3d::transform(rb->get_pose(), rb->get_velocity());
      FILE_LOG(LOG_SIMULATOR) << "  spatial velocity (after, body frame): " << sv << std::endl;
    }
  }

  // tabulate dynamics computation
  tms stop;  
  times(&stop);
  dynamics_utime += (double) (stop.tms_utime-start.tms_utime)/CLOCKS_PER_SEC;
  dynamics_stime += (double) (stop.tms_stime-start.tms_stime)/CLOCKS_PER_SEC;

  return step_size;
}

