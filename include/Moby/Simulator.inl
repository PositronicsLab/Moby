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
Real Simulator::integrate(Real step_size, ForwardIterator begin, ForwardIterator end)
{
  // get the state-derivative for each dynamic body
  for (ForwardIterator i = begin; i != end; i++)
  {
    // integrate the body
    if (LOGGING(LOG_SIMULATOR))
    {
      VectorN q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (before): " << (*i)->get_generalized_coordinates(DynamicBody::eRodrigues, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (before): " << (*i)->get_generalized_velocity(DynamicBody::eAxisAngle, q) << std::endl;
    }
    (*i)->integrate(current_time, step_size, integrator);
    if (LOGGING(LOG_SIMULATOR))
    {
      VectorN q;
      FILE_LOG(LOG_SIMULATOR) << "  generalized coordinates (after): " << (*i)->get_generalized_coordinates(DynamicBody::eRodrigues, q) << std::endl;
      FILE_LOG(LOG_SIMULATOR) << "  generalized velocities (after): " << (*i)->get_generalized_velocity(DynamicBody::eAxisAngle, q) << std::endl;
    }
  }

  return step_size;
}

