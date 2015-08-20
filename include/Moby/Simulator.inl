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
  Ravelin::VectorNd gc, gv, ga, f, dv, lambda;
  std::vector<JointPtr> island_ijoints;

  // begin timing dynamics
  tms cstart;  
  clock_t start = times(&cstart);

  // get the simulator pointer
  boost::shared_ptr<Simulator> shared_this = boost::dynamic_pointer_cast<Simulator>(shared_from_this());

  // find islands
  std::vector<std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > > islands; 
  find_islands(islands);

  // integrate each island 
  for (unsigned i=0; i< islands.size(); i++)
  {
    // get the i'th island
    std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& island = islands[i];

    // sort the island so that we can search it
    std::sort(island.begin(), island.end());

    // clear the set of implicit joints 
    island_ijoints.clear();

    // get the implicit joints in the island
    for (unsigned j=0; j< implicit_joints.size(); j++)
    {
      // get the inboard and outboard links for the joint
      boost::shared_ptr<Ravelin::RigidBodyd> ib = implicit_joints[j]->get_inboard_link();
      boost::shared_ptr<Ravelin::RigidBodyd> ob = implicit_joints[j]->get_outboard_link();

      // get the super bodies 
      boost::shared_ptr<Ravelin::DynamicBodyd> ib_super = ib->get_super_body(); 
      boost::shared_ptr<Ravelin::DynamicBodyd> ob_super = ib->get_super_body(); 

      if (std::binary_search(island.begin(), island.end(), ib_super) ||
          std::binary_search(island.begin(), island.end(), ob_super))
        island_ijoints.push_back(implicit_joints[j]);
    }

    // get all implicit joints from articulated bodies in the island
    for (unsigned j=0; j< island.size(); j++)
    {
      // see whether the body is articulated
      boost::shared_ptr<Ravelin::ArticulatedBodyd> ab = boost::dynamic_pointer_cast<Ravelin::ArticulatedBodyd>(island[j]);
      if (!ab)
        continue;

      // get the implicit joints for this body
      const std::vector<boost::shared_ptr<Ravelin::Jointd> >& ijoints = ab->get_implicit_joints();

      // add the joints
      for (unsigned k=0; k< ijoints.size(); k++)
        island_ijoints.push_back(boost::dynamic_pointer_cast<Joint>(ijoints[k]));
    }

    // get number of implicit constraints
    const unsigned N_IMPLICIT = island_ijoints.size();
  
    // process each body in the island
    for (unsigned j=0; j< island.size(); j++)
    {
      // get the body
      boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(island[j]); 
      ControlledBodyPtr cb = boost::dynamic_pointer_cast<ControlledBody>(db); 

      // reset accumulators
      db->reset_accumulators();
      
      // add recurrent forces
      const std::list<RecurrentForcePtr>& rfs = cb->get_recurrent_forces();
      BOOST_FOREACH(RecurrentForcePtr rf, rfs)
        rf->add_force(db);

      // call the body's controller
      if (cb->controller)
      {
        FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << cb->id << std::endl;
        (*cb->controller)(cb, current_time, cb->controller_arg);
      }

      // no implicit constraints? just calculate forward dynamics for the body
      if (N_IMPLICIT == 0)
        db->calc_fwd_dyn();
    }

    // if there are implicit constraints, we'll need to go through the solve
    // process
    if (N_IMPLICIT > 0)
    {
      // get the total number of generalized coordinates for the island
      const unsigned NGC_TOTAL = num_generalized_coordinates(island);

      // setup f
      f.resize(NGC_TOTAL);
      for (unsigned i=0, gc_index = 0; i< island.size(); i++)
      {
        const unsigned NGC = island[i]->num_generalized_coordinates(Ravelin::DynamicBodyd::eSpatial);
        Ravelin::SharedVectorNd f_sub = f.segment(gc_index, gc_index + NGC);
        island[i]->get_generalized_forces(f_sub);
        gc_index += NGC;
      }

      // compute change in velocity and constraint forces
      solve(island, island_ijoints, f, dv, lambda);

      // set accelerations
      dv /= dt;
      for (unsigned i=0, gc_index = 0; i< island.size(); i++)
      {
        const unsigned NGC = island[i]->num_generalized_coordinates(Ravelin::DynamicBodyd::eSpatial);
        Ravelin::SharedConstVectorNd dv_sub = dv.segment(gc_index, gc_index + NGC);
        island[i]->set_generalized_acceleration(dv_sub);
        gc_index += NGC;
      }

      // populate constraint forces
      for (unsigned i=0, c_index = 0; i< island_ijoints.size(); i++)
      {
        const unsigned NEQ = island_ijoints[i]->num_constraint_eqns();
        Ravelin::SharedConstVectorNd lambda_sub = lambda.segment(c_index, c_index + NEQ);
        island_ijoints[i]->lambda = lambda_sub;
        c_index += NEQ;
      }
    }

    // update generalized velocities and use new generalized 
    // velocities to calculate new generalized coordinates 
    for (unsigned j=0; j< island.size(); j++)
    {
      boost::shared_ptr<Ravelin::DynamicBodyd> db = island[j]; 
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
  } 

  // tabulate dynamics computation
  tms cstop;  
  clock_t stop = times(&cstop);
  dynamics_time += (double) (stop-start)/CLOCKS_PER_SEC;

  return dt;
}


