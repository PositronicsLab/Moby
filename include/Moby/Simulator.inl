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
  Ravelin::VectorNd iMfext_dt, fext;
  std::vector<Ravelin::VectorNd> forces;
  std::vector<Ravelin::MatrixNd> inertias;
  std::vector<Ravelin::MatrixNd> body_jacobians;
  std::vector<unsigned> gc_indices;

  // begin timing dynamics
  tms cstart;  
  clock_t start = times(&cstart);

  // get the simulator pointer
  boost::shared_ptr<Simulator> shared_this = boost::dynamic_pointer_cast<Simulator>(shared_from_this());

  // find islands
  std::vector<std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > > islands; 
//  find_islands(islands);

  // integrate each island 
  for (unsigned i=0; i< islands.size(); i++)
  {
    // get the i'th island
    const std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& island = islands[i];

    // TODO: clear the island quantities 
//    body_jacobians.clear(island.size());

    // setup the GC index
    gc_indices.clear();
    gc_indices.push_back(0); 

    // TODO: look whether there are any implicit constraints
    unsigned n_implicit = 0;
  
    // TODO: get the total number of implicit constraints 

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

      // no implicit constraints? just calculate forward dynamics
      if (n_implicit == 0)
        db->calc_fwd_dyn();
/*
      else
      {
        // TODO: get the constraint Jacobian for any implicit constraints for the body

        // TODO: set Jacobian indices 

        // update gc indices
        gc_indices.push_back(db->num_generalized_coordinates(Ravelin::DynamicBodyd::eSpatial));
      }

      // only do this if at least one implicit constraint
      if (n_implicit > 0)
      {
        // get the Jacobians for the simulator constraints

        // compute inv(M)*J'

        // compute inv(M)*fext*dt
        fext.resize(gc_indices.back());
        iMfext_dt.resize(gc_indices.back());
        for (unsigned j=0; j< island.size(); j++)
        {
          // get the body
          boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(island[j]); 

          // get the part of the force vector
          Ravelin::SharedVectorNd f = fext.segment(gc_indices[j], gc_indices[j+1]);

          // compute fext
          db->get_generalized_forces(f);

          // get the subvector
          Ravelin::SharedVectorNd iMfext_dt_sub = iMfext_dt.segment(gc_indices[j], gc_indices[j+1]);
          
          // solve
          db->solve_generalized_inertia(f, iMfext_dt_sub);

          // scale by dt
          iMfext_dt_sub *= dt;
        } 

        // TODO: compute J*inv(M)*J'

        // TODO: compute J*inv(M)*fext*dt

        // prepare to compute tau = J*inv(M)*J' \ J*inv(M)*fext*h
        double lambda = 1e-16;
        double expo = 10.0;
        JiMJT_copy = JiMJT;
        while (!Ravelin::LinAlgd::factor_chol(JiMJT_copy))
        {
          // update JiMJT_copy using lambda
          JiMJT_copy = JiMJT;
          Ravelin::RowIteratord iter = JiMJT_copy.row_iterator_begin();
          for (unsigned i=0; i< JiMJT.rows(); i++, iter += (JiMJT.rows()+1))
            *iter += lambda;

          // update lambda
          lambda *= expo;
        }

        // compute tau = J*inv(M)*J' \ J*inv(M)*fext*dt
        Ravelin::LinAlgd::solve_chol_fast(JiMJT_copy, JiMfext_dt);

        // TODO: populate tau into constraint forces

        // compute change in velocity
        iMJT.mult(JiMfexth, dv);
        dv *= dt;
        dv -= iMfext_dt;
        dv.negate();

        // compute acceleration
        dv /= dt;

        // store the generalized accelerations
        for (unsigned j=0; j< island.size(); j++)
        {
          boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(island[j]); 
          Ravelin::SharedVectorNd a_sub = dv.segment(gc_indices[j], gc_indices[j+1]);
          db->set_generalized_acceleration(a_sub);
        } 
      }

      // update generalized velocities and use new generalized 
      // velocities to calculate new generalized coordinates 
      for (unsigned j=0; j< island.size(); j++)
      {
        boost::shared_ptr<Ravelin::DynamicBodyd> db = boost::dynamic_pointer_cast<Ravelin::DynamicBodyd>(island[j]); 
        db->get_generalized_acceleration(ga);
        db->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
        ga *= dt;
        gv += ga;
        db->set_generalized_velocity(DynamicBodyd::eSpatial, gv);

        // integrate the generalized position forward
        db->get_generalized_velocity(DynamicBodyd::eEuler, gv);
        db->get_generalized_coordinates(DynamicBodyd::eEuler, gc);
        gv *= dt;
        gc += gv;
        db->set_generalized_coordinates(DynamicBodyd::eSpatial, gc);
*/
    }
  } 

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


