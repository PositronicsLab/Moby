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
  tms cstart;  
  clock_t start = times(&cstart);

  // get the state-derivative for each dynamic body
  for (ForwardIterator i = begin; i != end; i++)
  {
    // integrate the body
    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (before): " << (*i)->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (before): " << (*i)->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
    }
    (*i)->integrate(current_time, step_size, integrator);
    if (LOGGING(LOG_SIMULATOR))
    {
      Ravelin::VectorNd q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (after): " << (*i)->get_generalized_coordinates(DynamicBody::eEuler, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (after): " << (*i)->get_generalized_velocity(DynamicBody::eSpatial, q) << std::endl;
    }
  }

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;

  return step_size;
}

