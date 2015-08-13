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

  // begin timing dynamics
  tms cstart;  
  clock_t start = times(&cstart);

  // get the simulator pointer
  boost::shared_ptr<Simulator> shared_this = boost::dynamic_pointer_cast<Simulator>(shared_from_this());

  // integrate each body
  for (ForwardIterator i = begin; i != end; i++)
  {
    // cast the body as a Ravelin dynamic body
    boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(*i); 

    // if the body is kinematically updated, call its controller and otherwise
    // ignore it
    if ((*i)->get_kinematic() && (*i)->controller)
    {
      (*(*i)->controller)(*i, current_time, (*i)->controller_arg);
      continue;
    }

    // calculate forward dynamics
    db->calc_fwd_dyn();

    // integrate the generalized velocity forward
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
    db->set_generalized_coordinates(Ravelin::DynamicBodyd::eSpatial, gc);
  }

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;

  return dt;
}

