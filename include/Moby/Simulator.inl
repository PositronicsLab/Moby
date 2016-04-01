/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Integrates both position and velocity of rigid _bodies
/**
 * \return the size of step taken
 */
template <class ForwardIterator>
double Simulator::integrate(double dt, ForwardIterator begin, ForwardIterator end)
{
  Ravelin::VectorNd gc, gv, ga, gve;
  std::vector<JointPtr> island_ijoints;

  // begin timing dynamics
  tms cstart;  
  clock_t start = times(&cstart);

  // pre-calculate dynamics
  precalc_fwd_dyn();

  // calculate forward dynamics
  calc_fwd_dyn(dt);

  // update generalized velocities and use new generalized 
  // velocities to calculate new generalized coordinates 
  for (unsigned j=0; j< _bodies.size(); j++)
  {
    // NOTE: coordinates must be set first to ensure that frame information
    // is correct for velocities: mixed pose will be incorrect if
    // generalized_velocity is set first
    boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(_bodies[j]); 
    db->get_generalized_acceleration(ga);
    db->get_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);
    db->get_generalized_velocity(Ravelin::DynamicBodyd::eEuler, gve);
    ga *= dt;
    gv += ga;

   // dissipate some energy
   if (dissipator)
   {
     std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > bodies;
     BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
       bodies.push_back(boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(cb));
     dissipator->apply(bodies);
   }

    // integrate the generalized position forward
    gve *= dt;
    db->get_generalized_coordinates_euler(gc);
    gc += gve;
    db->set_generalized_coordinates_euler(gc);
    db->set_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);
  } 

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;

  return dt;
}


