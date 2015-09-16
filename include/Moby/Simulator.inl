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
  Ravelin::VectorNd gc, gv, ga;
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
    boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(_bodies[j]); 
    db->get_generalized_acceleration(ga);
    db->get_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);
    ga *= dt;
    gv += ga;
    db->set_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);

    // integrate the generalized position forward
    db->get_generalized_velocity(Ravelin::DynamicBodyd::eEuler, gv);
    db->get_generalized_coordinates(Ravelin::DynamicBodyd::eEuler, gc);
    gv *= dt;
    gc += gv;
    db->set_generalized_coordinates(Ravelin::DynamicBodyd::eEuler, gc);
  } 

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;

  return dt;
}


