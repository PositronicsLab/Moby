/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <map>
#include <iostream>
#ifdef USE_OSG
#include <osg/Group>
#endif
#include <Moby/RecurrentForce.h>
#include <Moby/Dissipation.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/XMLTree.h>
#include <Moby/SparseJacobian.h>
#include <Moby/Simulator.h>

using std::list;
using std::set;
using std::multimap;
using std::queue;
using std::map;
using std::vector;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Sets up the simulator
/**
 * The simulator properties are set as follows:
 * <ul>
 * <li>simulator time = 0</li>
 * <li>no integrator</li>
 * </ul>
 */
Simulator::Simulator()
{
  this->current_time = 0;
  post_step_callback_fn = NULL;

  // clear dynamics timings
  dynamics_time = (double) 0.0;

  // setup the persistent and transient visualization data
  #ifdef USE_OSG
  _persistent_vdata = new osg::Group;
  _transient_vdata = new osg::Group;

  // add references to the visualization data
  _persistent_vdata->ref();
  _transient_vdata->ref();
  #endif
}

Simulator::~Simulator()
{
  #ifdef USE_OSG
  _persistent_vdata->unref();
  _transient_vdata->unref();
  #endif
}

/// Computes the ODE of the system
VectorNd& Simulator::ode(const VectorNd& x, double t, double dt, void* data, VectorNd& dx)
{
  // get the simulator
  shared_ptr<Simulator>& s = *((shared_ptr<Simulator>*) data);

  FILE_LOG(LOG_SIMULATOR) << "Simulator::ode(t=" << t << ") entered" << std::endl;

  // see whether t=current time and the derivative has already been computed
  if (t == s->current_time && s->_current_dx.size() > 0)
  {
    dx = s->_current_dx;
    return dx;
  }

  // initialize the ODE index
  unsigned idx = 0;

  // resize dx
  dx.resize(x.size());

  // loop through all bodies, preparing to compute the ODE
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // cast the body as a Ravelin dynamic body
    shared_ptr<DynamicBodyd> rdb = dynamic_pointer_cast<DynamicBodyd>(db);

    // get the number of generalized coordinates and velocities
    const unsigned NGC = rdb->num_generalized_coordinates(DynamicBodyd::eEuler);
    const unsigned NGV = rdb->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // get x for the body 
    SharedConstVectorNd xsub = x.segment(idx, idx+NGC+NGV);

    // compute the ODE
    db->prepare_to_calc_ode(xsub, t, dt, &db); 

    // update idx
    idx += NGC+NGV;
  }

  // check pairwise constraint violations
  s->check_pairwise_constraint_violations(t);

  // loop through all bodies, computing forward dynamics 
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
    if (!db->get_kinematic())
      dynamic_pointer_cast<DynamicBodyd>(db)->calc_fwd_dyn();

  // reset the index
  idx = 0;

  // loop through all bodies, computing the ODE
  BOOST_FOREACH(ControlledBodyPtr db, s->_bodies)
  {
    if (db->get_kinematic())
      continue;

    // cast the body as a Ravelin dynamic body
    shared_ptr<DynamicBodyd> rdb = dynamic_pointer_cast<DynamicBodyd>(db);

    // get the number of generalized coordinates and velocities
    const unsigned NGC = rdb->num_generalized_coordinates(DynamicBodyd::eEuler);
    const unsigned NGV = rdb->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // get dx for the body
    SharedVectorNd dxsub = dx.segment(idx, idx+NGC+NGV);

    // compute the ODE
    db->ode(t, dt, &db, dxsub); 

    // update idx
    idx += NGC+NGV;
  }

  FILE_LOG(LOG_SIMULATOR) << "Simulator::ode(t=" << t << ") exited" << std::endl;

  // see whether to set current time derivative
  if (t == s->current_time)
    s->_current_dx = dx;

  // return the ODE
  return dx;
}

/// Steps the Simulator forward in time without contact
/**
 * This pseudocode was inspired from [Baraff 1997] and [Mirtich 1996].
 * \param step_size the step size
 * \return step_size
 */
double Simulator::step(double step_size)
{
  #ifdef USE_OSG
  // clear one-step visualization data
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif

  // compute forward dynamics and integrate 
  current_time += integrate(step_size);

  // TODO: do any constraint stabilization
//  _cstab.stabilize(simulator);

  // call the callback
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  return step_size;
}

/// Finds the dynamic body in the simulator, if any
/**
 * Searches unarticulated bodies, articulated bodies, and links of
 * articulated bodies.
 */
ControlledBodyPtr Simulator::find_dynamic_body(const std::string& name) const
{
  BOOST_FOREACH(ControlledBodyPtr body, _bodies)
    if (body->id == name)
      return body;

  // failed, look through all links of articulated bodies
  BOOST_FOREACH(ControlledBodyPtr body, _bodies)
  {
    // try to cast the i'th DynamicBody as an ArticulatedBody
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);
    if (!ab)
      continue;
    
    // it was castable, get all links
    const vector<shared_ptr<RigidBodyd> >& links = ab->get_links();
    
    // look through all links for one matching the name
    BOOST_FOREACH(shared_ptr<RigidBodyd> rb, links)  
      if (rb->body_id == name)
        return dynamic_pointer_cast<RigidBody>(rb);
  }
    
  return ControlledBodyPtr();
}

/// Removes a dynamic body from the simulator
void Simulator::remove_dynamic_body(ControlledBodyPtr body)
{
  // remove the body from the list of bodies
  std::vector<ControlledBodyPtr>::iterator i = std::find(_bodies.begin(), _bodies.end(), body);
  if (i == _bodies.end())
    return;
  else
    _bodies.erase(i);

  #ifdef USE_OSG
  // see whether the body is articulated 
  ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(body);
  if (abody)
  {
    
    // remove visualization data for all links to the persistent visualization data
    const vector<shared_ptr<RigidBodyd> >& links = abody->get_links();
    for (unsigned i=0; i< links.size(); i++)
    {
      RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(links[i]);
      osg::Node* link_vdata = link->get_visualization_data();
      if (link_vdata)
        _persistent_vdata->removeChild(link_vdata);
    }
    
    // remove visualization data for all joints
    const vector<shared_ptr<Jointd> >& joints = abody->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
    {
      JointPtr joint = dynamic_pointer_cast<Joint>(joints[i]);
      osg::Node* joint_vdata = joint->get_visualization_data();
      if (joint_vdata)
        _persistent_vdata->removeChild(joint_vdata);
    }
  }
  else
  {
    // try rigid body
    RigidBodyPtr rigidbody = dynamic_pointer_cast<RigidBody>(body);
    assert(rigidbody);
    osg::Node* rb_vdata = rigidbody->get_visualization_data();
    if (rb_vdata)
      _persistent_vdata->removeChild(rb_vdata);
  }
  #endif
}

/// Adds a dynamic body to the simulator
/**
 * \pre list of bodies is sorted
 */
void Simulator::add_dynamic_body(ControlledBodyPtr body) 
{
  // if the body is already present in the simulator, skip it
  if (std::find(_bodies.begin(), _bodies.end(), body) != _bodies.end())
    return;

  #ifdef USE_OSG
  // see whether the body is articulated 
  ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(body);
  if (abody)
  {
    // add visualization data for all links to the persistent visualization data
    const vector<shared_ptr<RigidBodyd> >& links = abody->get_links();
    for (unsigned i=0; i< links.size(); i++)
    {
      RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(links[i]);
      osg::Node* link_vdata = link->get_visualization_data();
      if (link_vdata)
        _persistent_vdata->addChild(link_vdata);
    }
    
    // get visualization data for all joints
    const vector<shared_ptr<Jointd> >& joints = abody->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
    {
      JointPtr joint = dynamic_pointer_cast<Joint>(joints[i]);
      osg::Node* joint_vdata = joint->get_visualization_data();
      if (joint_vdata)
        _persistent_vdata->addChild(joint_vdata);
    }
  }
  else
  {
    // it must be a rigid body
    RigidBodyPtr rigidbody = dynamic_pointer_cast<RigidBody>(body);
    assert(rigidbody);

    // get the visualization data and add it to the simulator
    osg::Node* rb_vdata = rigidbody->get_visualization_data();
    if (rb_vdata)
      _persistent_vdata->addChild(rb_vdata);
  }
  #endif
  
  // add the body to the list of bodies and sort the list of bodies
  _bodies.push_back(body); 
  std::sort(_bodies.begin(), _bodies.end());
}

/// Updates all visualization under the simulator
void Simulator::update_visualization()
{
  BOOST_FOREACH(ControlledBodyPtr body, _bodies)
    body->update_visualization();
}

/// Adds transient visualization data to the simulator
void Simulator::add_transient_vdata(osg::Node* vdata)
{
  #ifdef USE_OSG
  _transient_vdata->addChild(vdata);
  #endif
}

/// Prepares to calculate forward dynamics for bodies
void Simulator::precalc_fwd_dyn()
{
  // clear force accumulators, then add all recurrent and compliant
  // constraint forces
  BOOST_FOREACH(ControlledBodyPtr db, _bodies)
  {
    // get the body as a Ravelin dynamic body
    shared_ptr<DynamicBodyd> rdb = dynamic_pointer_cast<DynamicBodyd>(db);

    // clear the force accumulators on the body
    rdb->reset_accumulators();

    // add all recurrent forces on the body
    const list<RecurrentForcePtr>& rfs = db->get_recurrent_forces();
    BOOST_FOREACH(RecurrentForcePtr rf, rfs)
      rf->add_force(rdb);
    
    // call the body's controller
    if (db->controller)
    {
      FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << db->id << std::endl;
      (*db->controller)(db, current_time, db->controller_arg);
    }
  }
}

/// Calculates forward dynamics for bodies (does not consider unilateral constraints)
void Simulator::calc_fwd_dyn(double dt)
{
  VectorNd v, a, lambda, f, tmpv;
  MatrixNd M;
  vector<JointPtr> island_ijoints; 

  // get the simulator pointer
  shared_ptr<Simulator> shared_this = dynamic_pointer_cast<Simulator>(shared_from_this());

  // find islands
  vector<vector<shared_ptr<DynamicBodyd> > > islands; 
  find_islands(islands);

  // calculate forward dynamics for each island 
  for (unsigned i=0; i< islands.size(); i++)
  {
    // get the i'th island
    vector<shared_ptr<DynamicBodyd> >& island = islands[i];

    // sort the island so that we can search it
    std::sort(island.begin(), island.end());

    // clear the set of implicit joints 
    island_ijoints.clear();

    // get the implicit joints in the island
    for (unsigned j=0; j< implicit_joints.size(); j++)
    {
      // get the inboard and outboard links for the joint
      shared_ptr<RigidBodyd> ib = implicit_joints[j]->get_inboard_link();
      shared_ptr<RigidBodyd> ob = implicit_joints[j]->get_outboard_link();

      // get the super bodies 
      shared_ptr<DynamicBodyd> ib_super = ib->get_super_body(); 
      shared_ptr<DynamicBodyd> ob_super = ob->get_super_body(); 

      if (std::binary_search(island.begin(), island.end(), ib_super) ||
          std::binary_search(island.begin(), island.end(), ob_super))
        island_ijoints.push_back(implicit_joints[j]);
    }

    // get all implicit joints from articulated bodies in the island
    for (unsigned j=0; j< island.size(); j++)
    {
      // see whether the body is articulated
      shared_ptr<ArticulatedBodyd> ab = dynamic_pointer_cast<ArticulatedBodyd>(island[j]);
      if (!ab)
        continue;

      // get the implicit joints for this body
      const vector<shared_ptr<Jointd> >& ijoints = ab->get_implicit_joints();

      // add the joints
      for (unsigned k=0; k< ijoints.size(); k++)
        island_ijoints.push_back(dynamic_pointer_cast<Joint>(ijoints[k]));
    }

    // get number of implicit constraints
    const unsigned N_IMPLICIT = island_ijoints.size();
  
    // if there are no implicit constraints, just call calc_fwd_dyn(.) on
    // each body
    if (N_IMPLICIT == 0)
    {
      for (unsigned j=0; j< island.size(); j++)
      {
        // get the body
        shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(island[j]); 

        // no implicit constraints? just calculate forward dynamics for the body
        db->calc_fwd_dyn();
      }
    }
    else // there are implicit constraints - must go through the solve process
    {
      // get the total number of generalized coordinates for the island
      const unsigned NGC_TOTAL = num_generalized_coordinates(island);

      // setup f
      f.resize(NGC_TOTAL);
      for (unsigned i=0, gc_index = 0; i< island.size(); i++)
      {
        const unsigned NGC = island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
        SharedVectorNd f_sub = f.segment(gc_index, gc_index + NGC);
        island[i]->get_generalized_forces(f_sub);
        gc_index += NGC;
      }

      // add in current momenta to f
      for (unsigned i=0, gc_index = 0; i< island.size(); i++)
      {
        const unsigned NGC = island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
        SharedVectorNd f_sub = f.segment(gc_index, gc_index + NGC);
        island[i]->get_generalized_velocity(DynamicBodyd::eSpatial, v);
        island[i]->get_generalized_inertia(M);
        M.mult(v, tmpv);
        f_sub += tmpv;
        gc_index += NGC;
      }

      // compute acceleration and constraint forces
      solve(island, island_ijoints, f, a, lambda);

      // set accelerations
      for (unsigned i=0, gc_index = 0; i< island.size(); i++)
      {
        const unsigned NGC = island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
        SharedConstVectorNd a_sub = a.segment(gc_index, gc_index + NGC);
        island[i]->set_generalized_acceleration(a_sub);
        gc_index += NGC;
      }

      // populate constraint forces
      for (unsigned i=0, c_index = 0; i< island_ijoints.size(); i++)
      {
        const unsigned NEQ = island_ijoints[i]->num_constraint_eqns();
        SharedConstVectorNd lambda_sub = lambda.segment(c_index, c_index + NEQ);
        island_ijoints[i]->lambda = lambda_sub;
        c_index += NEQ;
      }
    }
  }
}

// Solves | M   J' | | x      | = | f | 
//        | J   0  | | lambda |   | 0 |
// using M*x = -J'*lambda + f and
// J*inv(M)*J'lambda = J*inv(M)*f
void Simulator::solve(const vector<shared_ptr<DynamicBodyd> >& island, const vector<JointPtr>& island_ijoints, const VectorNd& f, VectorNd& x, VectorNd& lambda) const
{
  MatrixNd JiMJT_frr, JiM, iMJT, iMJT_frr, JiMJT, Jm, tmp;
  VectorNd JiMf_frr, JiMf, iMf, lambda_sub; 
  const unsigned N_SPATIAL = 6;
  map<shared_ptr<DynamicBodyd>, unsigned> gc_map;
  std::vector<MatrixBlock> inv_inertias;

  // get dynamic bodies in the island and total number of generalized coords
  unsigned NGC_TOTAL = num_generalized_coordinates(island);

  // get the number of implicit equations in the island
  unsigned n_implicit_eqns = 0;
  for (unsigned i=0; i< island_ijoints.size(); i++)
    n_implicit_eqns += island_ijoints[i]->num_constraint_eqns();

  // if there are no implicit equations, just solve with generalized inertia
  // matrices
  if (n_implicit_eqns == 0)
  {
    // resize lambda
    lambda.resize(0);

    for (unsigned i=0, gc_index = 0; i< island.size(); i++)
    {
      // get the body
      shared_ptr<DynamicBodyd> db = island[i];

      // get the components of f and x
      const unsigned NGC = db->num_generalized_coordinates(DynamicBodyd::eSpatial);
      SharedConstVectorNd f_sub = f.segment(gc_index, gc_index + NGC);
      SharedVectorNd x_sub = x.segment(gc_index, gc_index + NGC);

      // solve directly
      db->solve_generalized_inertia(f_sub, x_sub);

      // update gc_index
      gc_index += NGC;
    } 

    return;      
  }

  // setup the gc map
  for (unsigned i=0, gc_index = 0; i< island.size(); i++)
  {
    gc_map[island[i]] = gc_index;
    gc_index += island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
  }

  // compute iMf
  iMf.resize(NGC_TOTAL);
  for (unsigned i=0, gc_index = 0; i< island.size(); i++)
  {
    // add the new inverse inertia block
    inv_inertias.push_back(MatrixBlock());

    // compute the inverse of the generalized inertia matrix
    island[i]->get_generalized_inertia(inv_inertias.back().block);
    LinAlgd::inverse_SPD(inv_inertias.back().block);

    // get the number of generalized coordinates for this body
    const unsigned NGC = island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // get the appropriate parts of the vectors
    SharedVectorNd iMf_sub = iMf.segment(gc_index, gc_index + NGC);
    SharedConstVectorNd f_sub = f.segment(gc_index, gc_index + NGC);

    // do the multiplication
    inv_inertias.back().block.mult(f_sub, iMf_sub);

    // setup the starting row and column indices of the block
    inv_inertias.back().st_col_idx = inv_inertias.back().st_row_idx = gc_index;

    // update the generalized coordinate index 
    gc_index += NGC;
  }

  // form Jacobians here
  SparseJacobian J;
  J.rows = n_implicit_eqns;
  J.cols = NGC_TOTAL;
  for (unsigned i=0, eq_idx=0; i< island_ijoints.size(); i++)
  {
    // get the inboard and outboard links
    shared_ptr<RigidBodyd> inboard = island_ijoints[i]->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = island_ijoints[i]->get_outboard_link();

    // add the block to the Jacobian
    if (inboard->is_enabled())
    {
      J.blocks.push_back(MatrixBlock());
      island_ijoints[i]->calc_constraint_jacobian(true, J.blocks.back().block);
      J.blocks.back().st_row_idx = eq_idx;
      J.blocks.back().st_col_idx = gc_map[inboard];
    }

    // add the block to the Jacobian
    if (outboard->is_enabled())
    {
      J.blocks.push_back(MatrixBlock());
      island_ijoints[i]->calc_constraint_jacobian(false, J.blocks.back().block);
      J.blocks.back().st_row_idx = eq_idx;
      J.blocks.back().st_col_idx = gc_map[outboard];
    }

    // update the equation index
    eq_idx += island_ijoints[i]->num_constraint_eqns();
  } 

  if (LOGGING(LOG_DYNAMICS))
  {
    J.to_dense(tmp);
    FILE_LOG(LOG_DYNAMICS) << "dense J: " << std::endl << tmp;
  }

  // (J*inv(M)*J') * lambda = J*inv(M)*f
  J.mult(inv_inertias, NGC_TOTAL, JiM);
  MatrixNd::transpose(JiM, iMJT);
  J.mult(iMJT, JiMJT);
  JiM.mult(f, JiMf);

  // form the biggest full rank matrix
  vector<bool> indices(JiMJT.rows(), false);
  bool last_successful = false;
  for (unsigned i=0, n_active=0; i< indices.size(); i++)
  {
    // see whether the number of indices is maximized
    if (n_active == NGC_TOTAL)
      break;

    // update the number of active indices
    n_active++;

    // try to add this index
    indices[i] = true;

    // get the submatrix
    JiMJT.select_square(indices, JiMJT_frr);

    // attempt to factorize it
    if (!LinAlgd::factor_chol(JiMJT_frr))
    {
      n_active--;
      indices[i] = false;
      last_successful = false;
    }
    else
      last_successful = true;
  }

  // get the appropriate rows of J*inv(M)*f
  JiMf.select(indices, JiMf_frr);
  FILE_LOG(LOG_DYNAMICS) << "J*inv(M)*f: " << JiMf_frr << std::endl;

  // if the last factorization was not successful, refactor
  if (!last_successful)
  {
    JiMJT.select_square(indices, JiMJT_frr);
    LinAlgd::factor_chol(JiMJT_frr);
  }

  // solve for lambda_sub
  lambda_sub = JiMf_frr;
  LinAlgd::solve_chol_fast(JiMJT_frr, lambda_sub);

  // reform iMJT
  iMJT.select_columns(indices, iMJT_frr);

  // set lambda
  lambda.set_zero(n_implicit_eqns);
  lambda.set(indices, lambda_sub);
  FILE_LOG(LOG_DYNAMICS) << "lambda: " << lambda << std::endl;

  // now compute x using M*x = -J'*lambda + f and
  iMJT_frr.mult(lambda_sub, x);
  x.negate();
  x += iMf;

  FILE_LOG(LOG_DYNAMICS) << "dv: " << x << std::endl;
}

/// Gets the number of generalized coordinates in an island
unsigned Simulator::num_generalized_coordinates(const vector<shared_ptr<DynamicBodyd> >& island) const
{
  unsigned ngc_total = 0;  
  for (unsigned i=0; i< island.size(); i++)
    ngc_total += island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
  return ngc_total;
}

/// Implements Base::load_from_xml()
void Simulator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  std::list<shared_ptr<const XMLTree> > child_nodes;
  std::map<std::string, BasePtr>::const_iterator id_iter;

  // load parent data
  Base::load_from_xml(node, id_map);

  // ***********************************************************************
  // don't verify that the node is correct, b/c Simulator can be subclassed
  // ***********************************************************************
  
  // get the current time 
  XMLAttrib* time_attr = node->get_attrib("current-time");
  if (time_attr)
    this->current_time = time_attr->get_real_value();

  // get the dissipator, if any
  XMLAttrib* diss_attr = node->get_attrib("dissipator-id");
  if (diss_attr)
  {
    const std::string& id = diss_attr->get_string_value(); 
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "Simulator::load_from_xml() - could not find" << std::endl;
      std::cerr << "  dissipator w/ID: " << id << " from offending node: ";
      std::cerr << std::endl << *node;
    }
    else
      dissipator = dynamic_pointer_cast<Dissipation>(id_iter->second);
  }

  // get all dynamic bodies used in the simulator
  child_nodes = node->find_child_nodes("DynamicBody");
  if (!child_nodes.empty())
  {
    // safe to clear the vector of bodies
    _bodies.clear();

    // process all DynamicBody child nodes
    for (std::list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      // verify that the dynamic-body-id attribute exists
      XMLAttrib* id_attr = (*i)->get_attrib("dynamic-body-id");

      // make sure that the ID exists
      if (!id_attr)
      {
        std::cerr << "Simulator::load_from_xml() - no dynamic-body-id ";
        std::cerr << "attribute in " << std::endl << "  offending node: ";
        std::cerr << *i << std::endl;
        continue;
      }

      // look for the dynamic body with that ID
      const std::string& id = id_attr->get_string_value(); 
      if ((id_iter = id_map.find(id))== id_map.end())
      {
        std::cerr << "Simulator::load_from_xml() - could not find" << std::endl;
        std::cerr << "  dynamic body w/ID: '" << id << "' from offending node:";
        std::cerr << std::endl << *node;
      }
      else
        add_dynamic_body(dynamic_pointer_cast<ControlledBody>(id_iter->second));
    }
  }

  // clear the vector of implicit constraints
  implicit_joints.clear();

  // get all dynamic bodies used in the simulator
  child_nodes = node->find_child_nodes("ImplicitConstraint");
  if (!child_nodes.empty())
  {
    // process all ImplicitConstraint nodes
    for (std::list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      // verify that the dynamic-body-id attribute exists
      XMLAttrib* id_attr = (*i)->get_attrib("joint-id");

      // make sure that the ID exists
      if (!id_attr)
      {
        std::cerr << "Simulator::load_from_xml() - no joint-id ";
        std::cerr << "attribute in " << std::endl << "  offending node: ";
        std::cerr << *i << std::endl;
        continue;
      }

      // look for the joint with that ID
      const std::string& id = id_attr->get_string_value(); 
      if ((id_iter = id_map.find(id))== id_map.end())
      {
        std::cerr << "Simulator::load_from_xml() - could not find" << std::endl;
        std::cerr << "  dynamic body w/ID: '" << id << "' from offending node:";
        std::cerr << std::endl << *node;
      }
      else
        implicit_joints.push_back(dynamic_pointer_cast<Joint>(id_iter->second));
    }
  }

  // get all recurrent forces used in the simulator -- note: this must be done
  // *after* all bodies have been loaded
  child_nodes = node->find_child_nodes("RecurrentForce");
  if (!child_nodes.empty())
  {
    // process all child nodes
    for (std::list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      // verify that the dynamic-body-id attribute exists
      XMLAttrib* id_attr = (*i)->get_attrib("recurrent-force-id");

      // make sure that the ID exists
      if (!id_attr)
      {
        std::cerr << "Simulator::load_from_xml() - no recurrent-force-id ";
        std::cerr << "attribute in tag: " << node << std::endl;
        continue;
      }

      // look for the recurrent force with that ID
      const std::string& id = id_attr->get_string_value(); 
      if ((id_iter = id_map.find(id)) == id_map.end())
      {
        std::cerr << "Simulator::load_from_xml() - could not find" << std::endl;
        std::cerr << "  recurrent force w/ID: " << id << " from offending node: " << std::endl << *node;
      }
      else
      {
        RecurrentForcePtr rf = dynamic_pointer_cast<RecurrentForce>(id_iter->second);
        BOOST_FOREACH(ControlledBodyPtr db, _bodies)
          db->get_recurrent_forces().push_back(rf);
      }
    }
  }
}

/// Finds islands
void Simulator::find_islands(vector<vector<shared_ptr<DynamicBodyd> > >& islands)
{
  // clear the islands
  islands.clear();

  // The way that we'll determine the constraint islands is to treat each rigid
  // body present in the constraints as a node in a graph; nodes will be connected
  // to other nodes if (a) they are both present in constraint or (b) they are
  // part of the same articulated body.  Nodes will not be created for disabled
  // bodies.
  set<shared_ptr<DynamicBodyd> > nodes;
  multimap<shared_ptr<DynamicBodyd>, shared_ptr<DynamicBodyd> > edges;
  typedef multimap<shared_ptr<DynamicBodyd>, shared_ptr<DynamicBodyd> >::const_iterator EdgeIter;

  // first, add all nodes
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    shared_ptr<RigidBodyd> rb = dynamic_pointer_cast<RigidBodyd>(db);
    if (rb && !rb->is_enabled())
      continue;
    nodes.insert(db);
  }

  // loop through all implicit joints in the simulator
  for (unsigned i=0; i< implicit_joints.size(); i++)
  {
    // get the inboard link and set the dynamic body
    shared_ptr<RigidBodyd> inboard = implicit_joints[i]->get_inboard_link();
    if (!inboard->is_enabled())
      continue;
    shared_ptr<ArticulatedBodyd> abi = inboard->get_articulated_body();
    shared_ptr<DynamicBodyd> dbi;
    if (abi)
      dbi = abi;
    else
      dbi = inboard;

    // get the outboard link and set the dynamic body
    shared_ptr<RigidBodyd> outboard = implicit_joints[i]->get_outboard_link();
    if (!outboard->is_enabled())
      continue;
    shared_ptr<ArticulatedBodyd> abo = outboard->get_articulated_body();
    shared_ptr<DynamicBodyd> dbo;
    if (abo)
      dbo = abo;
    else
      dbo = outboard;

    // add an edge between the bodies
    edges.insert(std::make_pair(dbi, dbo));
    edges.insert(std::make_pair(dbo, dbi));
  }

  // remove nodes from the set until there are no more nodes
  while (!nodes.empty())
  {
    // create a new island
    islands.push_back(vector<shared_ptr<DynamicBodyd> >());

    // get the node from the front
    shared_ptr<DynamicBodyd> node = *nodes.begin();

    // create a node queue, with this node added
    queue<shared_ptr<DynamicBodyd> > node_q;
    node_q.push(node);

    // loop until the queue is empty
    while (!node_q.empty())
    {
      // get the node off of the front of the node queue
    // get the node off of the front of the node queue
      node = node_q.front();
      node_q.pop();

      // erase the node from the set of nodes
      nodes.erase(node);

      // add all neighbors of the node that have not been processed already 
      // to the node queue
      std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
      for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
        if (nodes.find(i->second) != nodes.end())
          node_q.push(i->second);

      // add the node to the island
      islands.back().push_back(node);
    }
  }
}

/// Implements Base::save_to_xml()
void Simulator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent save_to_xml() method
  Base::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "Simulator";

  // save the current time 
  node->attribs.insert(XMLAttrib("current-time", this->current_time));

  // save the ID of the dissipator
  if (dissipator)
  {
    node->attribs.insert(XMLAttrib("dissipator-id", dissipator->id));
    shared_objects.push_back(dissipator);
  }

  // save the IDs of all dynamic bodies in the simulator
  BOOST_FOREACH(ControlledBodyPtr body, _bodies)
  {
    XMLTreePtr child_node(new XMLTree("DynamicBody"));
    node->add_child(child_node);
    child_node->attribs.insert(XMLAttrib("dynamic-body-id", body->id));
    if (!body)
      throw std::runtime_error("dynamic-body-id does not belong to a dynamic body");
    shared_objects.push_back(body);
  }
}

