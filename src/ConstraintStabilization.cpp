/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/**
 * Thoughts about separating bodies:
 *
 * If we really want to separate bodies with minimum work, then we'll have to
 * solve an optimal motion planning problem (fun!)- not going to be fast enough 
 * for most simulations. Also, this strategy would not be able to account for 
 * change in potential energy.
 *
 * If we apply a force to the bodies at the "wrong" contact point (like pushing
 * the box out of the plane), this will not lead to an increase in signed
 * distance.
 *
 * If we separate bodies using the vector of minimum translational distance as
 * the normal and we put a contact point halfway through the two bodies
 * (this actually will have no effect on torque), then that will locally
 * minimize the work necessary to separate the bodies. A root finding process
 * can be used to determine when the bodies will be separated.
 *
 * The other concern is what to do when the bodies are separated. If two
 * bodies are separated, then then the root finding process can determine
 * when (if ever) they come into contact.
 */

#include <map>
#include <algorithm>
#include <Moby/Types.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/ConstraintStabilization.h>
#include <boost/algorithm/minmax_element.hpp>
#include <utility>

using namespace Ravelin;
using namespace Moby;
using boost::weak_ptr;
using boost::shared_ptr;
using boost::shared_array;
using boost::const_pointer_cast; 
using boost::dynamic_pointer_cast; 
using std::map;
using std::vector;
using std::make_pair;
using std::endl;
using std::list;

ConstraintStabilization::ConstraintStabilization()
{
  // set maximum iterations to ten, by default 
  max_iterations = 10;

  // set unilateral tolerance to negative NEAR_ZERO by default
  eps = -NEAR_ZERO;

  // set bilateral tolerance
  bilateral_eps = 1e-6;
}

/// Saves the velocities before constraint stabilization
void ConstraintStabilization::save_velocities(shared_ptr<ConstraintSimulator> sim, vector<VectorNd>& qd)
{
  qd.clear();
  for (unsigned i = 0; i< sim->_bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(sim->_bodies[i]);
    qd.push_back(VectorNd());
    db->get_generalized_velocity(DynamicBodyd::eSpatial, qd.back());
  }
}

/// Restores the velocities after constraint stabilization
void ConstraintStabilization::restore_velocities(shared_ptr<ConstraintSimulator> sim, const vector<VectorNd>& qd)
{
  for (unsigned i = 0; i< sim->_bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(sim->_bodies[i]);
    db->set_generalized_velocity(DynamicBodyd::eSpatial, qd[i]);
  }
}

/// Gets the minimum pairwise distance
double ConstraintStabilization::get_min_pairwise_dist(const vector<PairwiseDistInfo>& pdi)
{
  // set distance to infinite initially
  double dist = std::numeric_limits<double>::max();

  for (unsigned i=0; i< pdi.size(); i++)
    dist = std::min(dist, pdi[i].dist);

  return dist;
}

/// Evaluates implicit bilateral constraints and returns the most significant violation
double ConstraintStabilization::evaluate_implicit_constraints(shared_ptr<ConstraintSimulator> sim, vector<double>& C)
{
  const unsigned SPATIAL_DIM = 6;
  double eval[SPATIAL_DIM];

  // clear C - we'll build it as we go
  C.clear();

  // no need to continue if there are no implicit constraints
  if(sim->implicit_joints.size() == 0)
    return 0.0;

  // evaluate all constraints in the simulator
  for (unsigned i=0; i< sim->implicit_joints.size(); i++)
  {
    sim->implicit_joints[i]->evaluate_constraints(eval);
    std::copy(eval, eval+sim->implicit_joints[i]->num_constraint_eqns(), std::back_inserter(C));
  }

  if (LOGGING(LOG_CONSTRAINT))
  {
    std::stringstream C_str;
    for (unsigned i=0; i< C.size(); i++)
      C_str << " " << C[i]; 
    FILE_LOG(LOG_CONSTRAINT) << "Constraint evaluations:" << C_str.str() << std::endl;
  }

  // get the maximum value and minimum value of C
  std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> mmelm = boost::minmax_element(C.begin(), C.end());
  return (std::fabs(*mmelm.first) > std::fabs(*mmelm.second)) ? std::fabs(*mmelm.first) : std::fabs(*mmelm.second);
}

/// Stabilizes the constraints in the simulator
void ConstraintStabilization::stabilize(shared_ptr<ConstraintSimulator> sim)
{
  FILE_LOG(LOG_SIMULATOR)<< "======constraint stabilization start======"<<std::endl;
  VectorNd dq, q, v;
  std::vector<UnilateralConstraintProblemData> pd;
  std::vector<VectorNd> qd_save;
  std::vector<double> C;

  std::map<shared_ptr<DynamicBodyd>, unsigned> body_index_map;

  // save the generalized velocities
  save_velocities(sim, qd_save);

  // get the body configurations and setup the body index map
  get_body_configurations(q, sim);
  generate_body_index_map(body_index_map, sim);

  // evaluate the bilateral constraints
  double max_vio = evaluate_implicit_constraints(sim, C);

  // get the pairwise distances
  vector<PairwiseDistInfo>& pdi = sim->_pairwise_distances;

  // see whether any pairwise distances are below epsilon
  double min_dist = get_min_pairwise_dist(pdi);

  FILE_LOG(LOG_SIMULATOR) <<"minimum pairwise distance (before stabilization loop): "<<min_dist<<std::endl;
  FILE_LOG(LOG_SIMULATOR) <<"maximum constraint violation (before stabilization loop): "<< max_vio <<std::endl;
  unsigned iterations = 0;
  while (min_dist < eps || max_vio > bilateral_eps)
  {
    FILE_LOG(LOG_SIMULATOR) <<"minimum pairwise distance: "<<min_dist<<std::endl;
    FILE_LOG(LOG_SIMULATOR) <<"maximum constraint violation: "<< max_vio <<std::endl;

    // zero body velocities first (we only want to change positions based on
    // our updates)
    for (unsigned i=0; i< sim->_bodies.size(); i++)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(sim->_bodies[i]);
      db->get_generalized_velocity(DynamicBodyd::eSpatial, v);
      v.set_zero();
      db->set_generalized_velocity(DynamicBodyd::eSpatial, v);
    }

    // compute problem data (get M, N, alpha, etc.) 
    compute_problem_data(pd, sim);

    // determine dq's
    dq.set_zero(q.size());
    for (unsigned i=0; i< pd.size(); i++)
      determine_dq(pd[i], dq, body_index_map);
    FILE_LOG(LOG_SIMULATOR) << "dq: " << dq << std::endl;

    // determine s and update q; NOTE: update q computes the pairwise distances 
    if (!update_q(dq, q, sim))
    {
      FILE_LOG(LOG_SIMULATOR) << " -- failed to effectively finish the constraint stabilization process!" << std::endl;
      break;
    }

    // update minimum distance
    min_dist = get_min_pairwise_dist(pdi);  

    // update constraint violation
    max_vio = evaluate_implicit_constraints(sim, C);

    iterations++;
    if (iterations == max_iterations)
    {
      FILE_LOG(LOG_SIMULATOR) << " -- maximum number of iterations reached" << std::endl;
      FILE_LOG(LOG_SIMULATOR) << " -- failed to effectively finish the constraint stabilization process!" << std::endl;
      break;
    }
  }

  // restore the generalized velocities
  restore_velocities(sim, qd_save);

  // after stabilization, velocities may be in an impacting state; correct
//  sim->find_unilateral_constraints(sim->contact_dist_thresh);
//  sim->calc_impacting_unilateral_constraint_forces(-1.0); 

  FILE_LOG(LOG_SIMULATOR) << iterations << " iterations required" << std::endl;
  FILE_LOG(LOG_SIMULATOR) <<"=====constraint stabilization end ======" << std::endl;
}

/// Adds unilateral constraints for joint limits in an articulated body
void ConstraintStabilization::add_articulate_limit_constraint(std::vector<UnilateralConstraint>& constraints, shared_ptr<ArticulatedBodyd> abody)
{
  ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(abody);
  std::vector<UnilateralConstraint> limits;
  ab->find_limit_constraints(std::back_inserter(limits));
  constraints.insert(constraints.end(), limits.begin(), limits.end());
}

/// Adds contact constraints from a pair of bodies
void ConstraintStabilization::add_contact_constraints(std::vector<UnilateralConstraint>& constraints, CollisionGeometryPtr cg1, CollisionGeometryPtr cg2, shared_ptr<ConstraintSimulator> sim)
{
  boost::shared_ptr<const Pose3d> GLOBAL;
  Point3d p1, p2;

  // compute the signed distance between the geometries
  double dist = sim->get_collision_detection()->calc_signed_dist(cg1, cg2, p1, p2);

  // case 1: bodies are clearly separated      
  if (dist >= NEAR_ZERO)
  {
    boost::shared_ptr<const Pose3d> GLOBAL;
    boost::shared_ptr<const Pose3d> pose1(p1.pose);
    boost::shared_ptr<const Pose3d> pose2(p2.pose);
    Transform3d _1T2 = Ravelin::Pose3d::calc_relative_pose(pose1, pose2);
    Transform3d _1TG = Ravelin::Pose3d::calc_relative_pose(pose1, GLOBAL);
    Point3d p1_2 = _1T2.transform_point(p1);
    Point3d p1_g = _1TG.transform_point(p1);
    Ravelin::Vector3d normal = p2-p1_2;
    normal.normalize();
    UnilateralConstraint uc = CollisionDetection::create_contact(cg1, cg2, p1_g, normal, dist);
    FILE_LOG(LOG_CONSTRAINT) << "creating contact between separated bodies: " << cg1->get_single_body()->body_id << " and " << cg2->get_single_body()->body_id << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << uc << std::endl;
    constraints.insert(constraints.end(), uc); 
  }
  // case 2: bodies are contacting/interpenetrated 
  else// (std::fabs(dist) < NEAR_ZERO)
  {
    FILE_LOG(LOG_CONSTRAINT) << "creating contacts between interpenetrating/contacting bodies: " << cg1->get_single_body()->body_id << " and " << cg2->get_single_body()->body_id << std::endl;
    const unsigned OLD_CONSTRAINTS_SZ = constraints.size();
    sim->_coldet->find_contacts(cg1,cg2, constraints);
    if (LOGGING(LOG_CONSTRAINT))
    {
      for (unsigned i=OLD_CONSTRAINTS_SZ; i< constraints.size(); i++)
        FILE_LOG(LOG_CONSTRAINT) << constraints[i] << std::endl;
    }
  }
}

/// Computes the constraint data
void ConstraintStabilization::compute_problem_data(std::vector<UnilateralConstraintProblemData>& pd_vector, shared_ptr<ConstraintSimulator> sim)
{
  VectorNd tmpv;

  // clear the constraints vector
  constraints.clear();

  // clear the problem data vector 
  pd_vector.clear();

  // get all bodies
  const std::vector<shared_ptr<ControlledBody> >& bodies = sim->_bodies;

  // get the collision detection mechanism
  shared_ptr<CollisionDetection> coldet = sim->get_collision_detection();

  // do broad phase collision detection here
  coldet->broad_phase(0.0, bodies, _pairs_to_check);

  // remove pairs that are unchecked
  for (unsigned i=0; i< _pairs_to_check.size(); )
    if (sim->unchecked_pairs.find(make_sorted_pair(_pairs_to_check[i].first, _pairs_to_check[i].second)) != sim->unchecked_pairs.end())
    {
      _pairs_to_check[i] = _pairs_to_check.back();
      _pairs_to_check.pop_back();
    }
    else
      i++;

  // now add contact constraints for each checked pair
  typedef std::pair<CollisionGeometryPtr, CollisionGeometryPtr> CheckPair;
  BOOST_FOREACH(CheckPair& to_check, _pairs_to_check)
    add_contact_constraints(constraints, to_check.first, to_check.second, sim);

  //FILE_LOG(LOG_SIMULATOR) << "constraints added" << std::endl;
  // 2) for each articulated body, add as many UnilateralConstraint objects as
  //    there are joints at their limits


  // 3) for each pair of bodies in interpenetrating contact, add a single
  //    point of contact at the deepest interpenetrating point with normal
  //    in the direction of the signed distance function. 

  // find islands
  list<vector<shared_ptr<DynamicBodyd> > > remaining_islands;
  list<list<UnilateralConstraint*> > islands;
  UnilateralConstraint::determine_connected_constraints(constraints, sim->implicit_joints, islands, remaining_islands);

  // process unilateral constraint islands
  BOOST_FOREACH(list<UnilateralConstraint*>& island, islands)
  {
    // setup a UnilateralConstraintProblemData object
    pd_vector.push_back(UnilateralConstraintProblemData());
    UnilateralConstraintProblemData& pd = pd_vector.back();

    // set a pointer to the simulator
    pd.simulator = sim;

    // put the constraint into the appropriate place
    BOOST_FOREACH(UnilateralConstraint* c, island)
    { 
      if (c->constraint_type == UnilateralConstraint::eContact)
        pd.contact_constraints.push_back(c);
      else
        pd.limit_constraints.push_back(c);
    }

    // set number of contact and limit constraints
    pd.N_CONTACTS = pd.contact_constraints.size();
    pd.N_LIMITS = pd.limit_constraints.size(); 
    pd.N_K_TOTAL = 0;
    pd.N_TRUE_CONE = 0;

    // now set the unilateral constraint data
    set_unilateral_constraint_data(pd);
    
    // set the elements of Cn_v and L_v
    // L_v is always set to zero
    // Cn_v is set to the signed distance between the two bodies
    Point3d pa,pb;
    for (unsigned i = 0; i < pd.contact_constraints.size(); ++i)
    {
      pd.Cn_v[i] = pd.contact_constraints[i]->signed_violation;//CollisionGeometry::calc_signed_dist(pd.contact_constraints[i]->contact_geom1, pd.contact_constraints[i]->contact_geom2, pa, pb) - std::sqrt(std::fabs(eps));
    }

    // set Jx_v
    double C[6];
    tmpv.resize(pd.N_CONSTRAINT_EQNS_IMP);
    for (unsigned i=0, j=0; i< pd.island_ijoints.size(); i++)
    {
      pd.island_ijoints[i]->evaluate_constraints(C);
      for (unsigned k=0; k< pd.island_ijoints[i]->num_constraint_eqns(); k++)
        tmpv[j++] = C[k];
    }

    // select proper constraints
    tmpv.select(pd.active, pd.Jx_v);    

    // clear all impulses
    for (unsigned i=0; i< pd.N_CONTACTS; i++)
      pd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
    for (unsigned i=0; i< pd.N_LIMITS; i++)
      pd.limit_constraints[i]->limit_impulse = 0;
  }

  // process islands composed completely of bilateral constraints
  BOOST_FOREACH(vector<shared_ptr<DynamicBodyd> >& island, remaining_islands)
  {
    // setup a UnilateralConstraintProblemData object
    pd_vector.push_back(UnilateralConstraintProblemData());
    UnilateralConstraintProblemData& pd = pd_vector.back();

    // set a pointer to the simulator
    pd.simulator = sim;

    // now set the unilateral constraint data
    set_bilateral_only_constraint_data(pd, island);
    
    // set Jx_v
    double C[6];
    tmpv.resize(pd.N_CONSTRAINT_EQNS_IMP);
    for (unsigned i=0, j=0; i< pd.island_ijoints.size(); i++)
    {
      pd.island_ijoints[i]->evaluate_constraints(C);
      for (unsigned k=0; k< pd.island_ijoints[i]->num_constraint_eqns(); k++)
        tmpv[j++] = C[k];
    }

    // select proper constraints
    tmpv.select(pd.active, pd.Jx_v);    

    // clear all impulses
    pd.contact_constraints.clear();
    pd.limit_constraints.clear();
  }
}

/// Gets the super body (articulated if any)
shared_ptr<DynamicBodyd> ConstraintStabilization::get_super_body(shared_ptr<DynamicBodyd> sb)
{
  shared_ptr<ArticulatedBodyd> ab = dynamic_pointer_cast<ArticulatedBodyd>(sb);
  if (ab)
    return ab;
  else
  {
    shared_ptr<RigidBodyd> rb = dynamic_pointer_cast<RigidBodyd>(sb);
    if (rb)
    {
      if (!rb->is_enabled())
        return shared_ptr<DynamicBodyd>();
      ab = rb->get_articulated_body();
      if (ab)
        return ab;
    }
    else
      return shared_ptr<DynamicBodyd>();
  }
}

/// Gets the super body (articulated if any)
shared_ptr<DynamicBodyd> ConstraintStabilization::get_super_body_from_rigid_body(shared_ptr<RigidBodyd> sb)
{
  if (!sb->is_enabled())
    return shared_ptr<DynamicBodyd>();
  shared_ptr<ArticulatedBodyd> ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return sb;
}

/// Computes the data to the LCP / QP problems
void ConstraintStabilization::set_bilateral_only_constraint_data(UnilateralConstraintProblemData& q, const vector<shared_ptr<DynamicBodyd> >& island)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned N_SPATIAL = 6;
  MatrixNd X, tmp;
  VectorNd v;

  // determine set of "super" bodies 
  q.super_bodies.clear();
  for (unsigned i=0; i< island.size(); i++)
  {
    // get the super bodies
    shared_ptr<DynamicBodyd> sb = get_super_body(island[i]);
    if (sb)
      q.super_bodies.push_back(sb);
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // prepare to compute the number of implicit constraint equations from islands
  q.N_CONSTRAINT_EQNS_IMP = 0;

  // add island implicit joint constraints to q
  for (unsigned i=0; i< q.simulator->implicit_joints.size(); i++)
  {
    // see whether a body from the implicit constraint matches a body in the
    // island
    JointPtr j = q.simulator->implicit_joints[i];
    shared_ptr<DynamicBodyd> in = j->get_inboard_link()->get_super_body();
    shared_ptr<DynamicBodyd> out = j->get_outboard_link()->get_super_body();
    if (std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), in))
    {
      q.island_ijoints.push_back(j);
      q.N_CONSTRAINT_EQNS_IMP += j->num_constraint_eqns();
    }
    else if ( std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), out))
    {
      q.island_ijoints.push_back(j);
      q.N_CONSTRAINT_EQNS_IMP += j->num_constraint_eqns();
    }
  }

  // set total number of generalized coordinates
  q.N_GC = 0;
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    q.N_GC += q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // initialize constants and set easy to set constants
  q.N_CONTACTS = q.N_LIMITS = 0;

  // setup constants related to articulated bodies
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
    if (abody) {
      q.N_CONSTRAINT_EQNS_IMP += abody->num_constraint_eqns_implicit();
    }
  }

  // setup no friction constraints
  q.N_LIN_CONE = q.N_K_TOTAL = q.N_TRUE_CONE = 0; 

  // initialize the problem matrices / vectors
  q.Cn_X_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CsT.set_zero(q.N_CONTACTS, 0);
  q.Cn_X_CtT.set_zero(q.N_CONTACTS, 0);
  q.Cn_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cn_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Cs_X_CsT.set_zero(0, 0);
  q.Cs_X_CtT.set_zero(0, 0);
  q.Cs_X_LT.set_zero(0, q.N_LIMITS);
  q.Cs_X_JxT.set_zero(0, q.N_CONSTRAINT_EQNS_IMP);
  q.Ct_X_CtT.set_zero(0, 0);
  q.Ct_X_LT.set_zero(0, q.N_LIMITS);
  q.Ct_X_JxT.set_zero(0, q.N_CONSTRAINT_EQNS_IMP);
  q.L_X_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.L_X_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  q.Jx_X_JxT.set_zero(q.N_CONSTRAINT_EQNS_IMP, q.N_CONSTRAINT_EQNS_IMP);
  q.Cn_v.set_zero(q.N_CONTACTS);
  q.Cs_v.set_zero(0);
  q.Ct_v.set_zero(0);
  q.L_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(0);
  q.ct.set_zero(0);
  q.l.set_zero(q.N_LIMITS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX;
  q.NCS_IDX = q.CT_IDX;
  q.NCT_IDX = q.NCS_IDX;
  q.L_IDX = q.NCT_IDX;
  q.N_VARS = q.L_IDX + q.N_LIMITS;

  // get generalized velocity
  ImpactConstraintHandler::get_generalized_velocity(q, v);

  // setup the gc map
  map<shared_ptr<DynamicBodyd>, unsigned> gc_map;
  for (unsigned i=0, gc_index = 0; i< q.super_bodies.size(); i++)
  {
    gc_map[q.super_bodies[i]] = gc_index;
    unsigned NGC = q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
    gc_index += NGC;
  }

  // get the total number of implicit constraint equations
  unsigned n_implicit_eqns = 0;
  for (unsigned i=0; i< q.island_ijoints.size(); i++)
    n_implicit_eqns += q.island_ijoints[i]->num_constraint_eqns();

  // prepare to setup Jacobian
  q.Jfull.rows = n_implicit_eqns;
  q.Jfull.cols = q.N_GC;

  // determine implicit constraint Jacobian
  for (unsigned i=0, eq_idx=0; i< q.island_ijoints.size(); i++)
  {
    // resize the temporary matrix
    tmp.resize(q.island_ijoints[i]->num_constraint_eqns(), N_SPATIAL);

    // get the inboard and outboard links
    shared_ptr<RigidBodyd> inboard = q.island_ijoints[i]->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = q.island_ijoints[i]->get_outboard_link();

    // compute the Jacobian w.r.t. the inboard link
    if (inboard->is_enabled())
    {
      // add the block to the Jacobian
      q.Jfull.blocks.push_back(MatrixBlock());
      q.island_ijoints[i]->calc_constraint_jacobian(true, q.Jfull.blocks.back().block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[inboard];
    }
 
    if (outboard->is_enabled())
    {
      // add the block to the Jacobian
      q.Jfull.blocks.push_back(MatrixBlock());
      q.island_ijoints[i]->calc_constraint_jacobian(false, q.Jfull.blocks.back().block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[outboard];
    }

    // update the equation index
    eq_idx += q.island_ijoints[i]->num_constraint_eqns();
  } 

  // determine active set of implicit constraints
  ImpactConstraintHandler::get_full_rank_implicit_constraints(q.Jfull, q.active);

  // compute X
  ImpactConstraintHandler::compute_X(q, X);

  // setup Jacobian for Cn
  SparseJacobian Cn;
  Cn.cols = q.N_GC;

  // setup Jacobian for J
  q.J.mult(X, tmp); MatrixNd::transpose(tmp, q.X_JxT);

  // setup X*other Jacobians
  q.X_CnT.resize(q.N_GC, 0);
  q.X_CsT.resize(q.N_GC, 0);
  q.X_CtT.resize(q.N_GC, 0);
  q.X_LT.resize(q.N_GC, 0);
}

/// Computes the data to the LCP / QP problems
void ConstraintStabilization::set_unilateral_constraint_data(UnilateralConstraintProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned N_SPATIAL = 6;
  VectorNd v;
  MatrixNd X, tmp;

  // determine set of "super" bodies from contact constraints
  q.super_bodies.clear();
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // get the rigid bodies involved
    shared_ptr<RigidBodyd> rb1 = dynamic_pointer_cast<RigidBodyd>(q.contact_constraints[i]->contact_geom1->get_single_body());
    shared_ptr<RigidBodyd> rb2 = dynamic_pointer_cast<RigidBodyd>(q.contact_constraints[i]->contact_geom2->get_single_body());

    // get the super bodies
    shared_ptr<DynamicBodyd> sb1 = get_super_body_from_rigid_body(rb1);
    shared_ptr<DynamicBodyd> sb2 = get_super_body_from_rigid_body(rb2);

    // add the super bodies, if desired
    if (sb1)
      q.super_bodies.push_back(sb1);
    if (sb2)
      q.super_bodies.push_back(sb2);
  }

  // determine set of "super" bodies from limit constraints
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    RigidBodyPtr outboard = q.limit_constraints[i]->limit_joint->get_outboard_link();
    q.super_bodies.push_back(get_super_body_from_rigid_body(outboard));
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // prepare to compute the number of implicit constraint equations from islands
  q.N_CONSTRAINT_EQNS_IMP = 0;

  // add island implicit joint constraints to q
  for (unsigned i=0; i< q.simulator->implicit_joints.size(); i++)
  {
    // see whether a body from the implicit constraint matches a body in the
    // island
    JointPtr j = q.simulator->implicit_joints[i];
    shared_ptr<DynamicBodyd> in = j->get_inboard_link()->get_super_body();
    shared_ptr<DynamicBodyd> out = j->get_outboard_link()->get_super_body();
    if (std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), in))
    {
      q.island_ijoints.push_back(j);
      q.N_CONSTRAINT_EQNS_IMP += j->num_constraint_eqns();
    }
    else if ( std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), out))
    {
      q.island_ijoints.push_back(j);
      q.N_CONSTRAINT_EQNS_IMP += j->num_constraint_eqns();
    }
  }

  // set total number of generalized coordinates
  q.N_GC = 0;
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    q.N_GC += q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // initialize constants and set easy to set constants
  q.N_CONTACTS = q.contact_constraints.size();
  q.N_LIMITS = q.limit_constraints.size();

  // setup constants related to articulated bodies
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
    if (abody) {
      q.N_CONSTRAINT_EQNS_IMP += abody->num_constraint_eqns_implicit();
    }
  }

  // setup no friction constraints
  q.N_LIN_CONE = q.N_K_TOTAL = q.N_TRUE_CONE = 0; 

  // initialize the problem matrices / vectors
  q.Cn_X_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CsT.set_zero(q.N_CONTACTS, 0);
  q.Cn_X_CtT.set_zero(q.N_CONTACTS, 0);
  q.Cn_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cn_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Cs_X_CsT.set_zero(0, 0);
  q.Cs_X_CtT.set_zero(0, 0);
  q.Cs_X_LT.set_zero(0, q.N_LIMITS);
  q.Cs_X_JxT.set_zero(0, q.N_CONSTRAINT_EQNS_IMP);
  q.Ct_X_CtT.set_zero(0, 0);
  q.Ct_X_LT.set_zero(0, q.N_LIMITS);
  q.Ct_X_JxT.set_zero(0, q.N_CONSTRAINT_EQNS_IMP);
  q.L_X_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.L_X_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  q.Jx_X_JxT.set_zero(q.N_CONSTRAINT_EQNS_IMP, q.N_CONSTRAINT_EQNS_IMP);
  q.Cn_v.set_zero(q.N_CONTACTS);
  q.Cs_v.set_zero(0);
  q.Ct_v.set_zero(0);
  q.L_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(0);
  q.ct.set_zero(0);
  q.l.set_zero(q.N_LIMITS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX;
  q.NCS_IDX = q.CT_IDX;
  q.NCT_IDX = q.NCS_IDX;
  q.L_IDX = q.NCT_IDX;
  q.N_VARS = q.L_IDX + q.N_LIMITS;

  // get generalized velocity
  ImpactConstraintHandler::get_generalized_velocity(q, v);

  // setup the gc map
  map<shared_ptr<DynamicBodyd>, unsigned> gc_map;
  for (unsigned i=0, gc_index = 0; i< q.super_bodies.size(); i++)
  {
    gc_map[q.super_bodies[i]] = gc_index;
    unsigned NGC = q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
    gc_index += NGC;
  }

  // setup limit indices
  q.limit_indices.resize(q.N_LIMITS);
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    UnilateralConstraint& u = *q.limit_constraints[i];
    unsigned idx = u.limit_joint->get_coord_index() + u.limit_dof;
    shared_ptr<DynamicBodyd> ab = u.limit_joint->get_articulated_body();
    q.limit_indices[i] = idx + gc_map[ab];
  }

  // get the total number of implicit constraint equations
  unsigned n_implicit_eqns = 0;
  for (unsigned i=0; i< q.island_ijoints.size(); i++)
    n_implicit_eqns += q.island_ijoints[i]->num_constraint_eqns();

  // prepare to setup Jacobian
  q.Jfull.rows = n_implicit_eqns;
  q.Jfull.cols = q.N_GC;

  // determine implicit constraint Jacobian
  for (unsigned i=0, eq_idx=0; i< q.island_ijoints.size(); i++)
  {
    // resize the temporary matrix
    tmp.resize(q.island_ijoints[i]->num_constraint_eqns(), N_SPATIAL);

    // get the inboard and outboard links
    shared_ptr<RigidBodyd> inboard = q.island_ijoints[i]->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = q.island_ijoints[i]->get_outboard_link();

    // compute the Jacobian w.r.t. the inboard link
    if (inboard->is_enabled())
    {
      // add the block to the Jacobian
      q.Jfull.blocks.push_back(MatrixBlock());
      q.island_ijoints[i]->calc_constraint_jacobian(true, q.Jfull.blocks.back().block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[inboard];
    }
 
    if (outboard->is_enabled())
    {
      // add the block to the Jacobian
      q.Jfull.blocks.push_back(MatrixBlock());
      q.island_ijoints[i]->calc_constraint_jacobian(false, q.Jfull.blocks.back().block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[outboard];
    }

    // update the equation index
    eq_idx += q.island_ijoints[i]->num_constraint_eqns();
  } 

  // determine active set of implicit constraints
  ImpactConstraintHandler::get_full_rank_implicit_constraints(q.Jfull, q.active);

  // compute X
  ImpactConstraintHandler::compute_X(q, X);

  // setup Jacobian for Cn
  SparseJacobian Cn;

  // setup the number of columns in each Jacobian
  Cn.cols = q.N_GC;

  // process all contact constraints
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
    add_contact_to_Jacobian(*q.contact_constraints[i], Cn, gc_map, i); 

  // compute X_CnT
  Cn.mult(X, tmp);  MatrixNd::transpose(tmp, q.X_CnT);
  q.J.mult(X, tmp); MatrixNd::transpose(tmp, q.X_JxT);

  // setup X_CsT, X_CtT
  q.X_CsT.resize(q.N_GC, 0);
  q.X_CtT.resize(q.N_GC, 0);
  
  // compute limit components - must do this first
  ImpactConstraintHandler::compute_limit_components(X, q);

  // compute problem data for Cn rows
  Cn.mult(q.X_CnT, q.Cn_X_CnT); 
  Cn.mult(q.X_LT,  q.Cn_X_LT);  
  Cn.mult(q.X_JxT,  q.Cn_X_JxT);

  // compute problem data for limit rows
  q.L_X_JxT.resize(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  for (unsigned i=0; i< q.N_LIMITS; i++)
    q.L_X_JxT.row(i) = q.X_JxT.row(q.limit_indices[i]);
}

/// Adds a contact constraint to the contact Jacobians
void ConstraintStabilization::add_contact_to_Jacobian(const UnilateralConstraint& c, SparseJacobian& Cn, const std::map<shared_ptr<DynamicBodyd>, unsigned>& gc_map, unsigned contact_idx)
{
  MatrixNd tmp, tmp2, Jm;

  // get the two single bodies involved in the contact
  shared_ptr<SingleBodyd> b1 = c.contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> b2 = c.contact_geom2->get_single_body();

  // get the two rigid bodies
  shared_ptr<RigidBodyd> rb1 = dynamic_pointer_cast<RigidBodyd>(b1);
  shared_ptr<RigidBodyd> rb2 = dynamic_pointer_cast<RigidBodyd>(b2);

  // get the super bodies
  shared_ptr<ArticulatedBodyd> su1 = dynamic_pointer_cast<ArticulatedBodyd>(b1->get_super_body());
  shared_ptr<ArticulatedBodyd> su2 = dynamic_pointer_cast<ArticulatedBodyd>(b2->get_super_body());

  // add a row to the Jacobian
  Cn.rows++;

  // do this six times, one for each body and each direction
  ImpactConstraintHandler::add_contact_dir_to_Jacobian(rb1, su1, Cn, c.contact_point, c.contact_normal, gc_map, contact_idx);
  ImpactConstraintHandler::add_contact_dir_to_Jacobian(rb2, su2, Cn, c.contact_point, -c.contact_normal, gc_map, contact_idx);
}

/// Computes deltaq by solving a linear complementarity problem
void ConstraintStabilization::determine_dq(UnilateralConstraintProblemData& pd, VectorNd& dqm, const std::map<shared_ptr<DynamicBodyd>, unsigned>& body_index_map)
{
  VectorNd z, dq_sub;

  // initialize the LCP matrix and LCP vector
  MatrixNd MM(pd.N_CONTACTS + pd.N_LIMITS, pd.N_CONTACTS + pd.N_LIMITS);
  VectorNd qq(pd.N_CONTACTS + pd.N_LIMITS);

  // setup the LCP matrix and LCP vector
  MM.block(0, pd.N_CONTACTS, 0, pd.N_CONTACTS) = pd.Cn_X_CnT;
  MM.block(0, pd.N_CONTACTS, pd.N_CONTACTS, MM.columns()) = pd.Cn_X_LT;
  SharedMatrixNd L_X_CnT_block = MM.block(pd.N_CONTACTS, MM.rows(), 0, pd.N_CONTACTS);
  MatrixNd::transpose(pd.Cn_X_LT, L_X_CnT_block);
  MM.block(pd.N_CONTACTS, MM.rows(), pd.N_CONTACTS, MM.columns()) = pd.L_X_LT;
  qq.segment(0, pd.N_CONTACTS) = pd.Cn_v;
  qq.segment(pd.N_CONTACTS, qq.size()) = pd.L_v;

  FILE_LOG(LOG_SIMULATOR) << "# of constraints in determine_dq(): " << qq.size() << std::endl;
  FILE_LOG(LOG_SIMULATOR) << "MM: " << std::endl << MM;
  FILE_LOG(LOG_SIMULATOR) << "qq: " << qq << std::endl;

  // solve N*inv(M)*N'*dq = N*alpha for impulses 
  if (!_lcp.lcp_fast(MM, qq, z))
    _lcp.lcp_lemke_regularized(MM, qq, z);

  // update velocities
  ImpactConstraintHandler::update_from_stacked(pd, z);

  // populate dq 
  for(unsigned i = 0; i< pd.super_bodies.size(); i++)
  {
    assert(body_index_map.find(pd.super_bodies[i]) != body_index_map.end());
    unsigned start = (body_index_map.find(pd.super_bodies[i]))->second;
    unsigned coord_num = pd.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eEuler);
    pd.super_bodies[i]->get_generalized_velocity(DynamicBodyd::eEuler, dq_sub);  
    dqm.segment(start, start+coord_num) = dq_sub;
  }
}

/// Updates the body velocities based on problem data
void ConstraintStabilization::update_velocities(const UnilateralConstraintProblemData& pd)
{
  VectorNd v; 
  map<shared_ptr<DynamicBodyd>, VectorNd> gj;
  map<shared_ptr<DynamicBodyd>, VectorNd>::iterator gj_iter;

  // loop over all contact constraints first
  for (unsigned i=0; i< pd.contact_constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& e = *pd.contact_constraints[i];
    SForced w(e.contact_impulse);

    // get the two single bodies of the contact
    shared_ptr<SingleBodyd> sb1 = e.contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> b1 = sb1->get_super_body();
    shared_ptr<DynamicBodyd> b2 = sb2->get_super_body();

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, v);
      gj_iter->second += v;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, v);
      gj_iter->second += v;
    }
  }

  // loop over all limit constraints next
  for (unsigned i=0; i< pd.limit_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *pd.limit_constraints[i];
    shared_ptr<ArticulatedBodyd> ab = e.limit_joint->get_articulated_body();

    // get the iterator for the articulated body
    gj_iter = gj.find(ab);

    // apply limit impulses to bodies in independent coordinates
    if (dynamic_pointer_cast<RCArticulatedBody>(ab))
    {
      // get the index of the joint
      unsigned idx = e.limit_joint->get_coord_index() + e.limit_dof;

      // initialize the vector if necessary
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(ab);
      }

      // set the limit force
      gj_iter->second[idx] += e.limit_impulse;
    }
    else
    {
      assert(false);
    }
  }

  // apply all generalized impacts
  for (map<shared_ptr<DynamicBodyd>, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/// Simple squaring function
static double sqr(double x) { return x*x; }

/// Updates q doing a backtracking line search
/**
 * Two body version
 */
bool ConstraintStabilization::update_q(const VectorNd& dq, VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  VectorNd qstar, grad;
  vector<double> C, C_old;
  const double MIN_T = NEAR_ZERO;

  // evaluate the implicit constraints
  evaluate_implicit_constraints(sim, C_old);

  // compute old bilateral constraint violations
  double old_bilateral_cvio = 0.0;
  for (unsigned i=0; i< C_old.size(); i++)
    old_bilateral_cvio += sqr(C_old[i]);
  old_bilateral_cvio = std::sqrt(old_bilateral_cvio);

  // copy the pairwise distances
  vector<PairwiseDistInfo>& pdi = sim->_pairwise_distances;
  vector<PairwiseDistInfo> pdi_old = pdi;

  if (LOGGING(LOG_CONSTRAINT))
    sim->calc_pairwise_distances();

  FILE_LOG(LOG_CONSTRAINT) << "...about to compute unilateral brackets" << std::endl;
  // find the pairwise distances and implicit constraint evaluations at q + dq
  qstar = dq;
  qstar += q;
  update_body_configurations(qstar, sim);
  sim->calc_pairwise_distances();
  evaluate_implicit_constraints(sim, C);

  // evaluate new pdi
  vector<PairwiseDistInfo> pdi_new = pdi;

  // we may have to find roots for all pairwise distances
  vector<bool> unilateral_bracket;
  for (unsigned i=0; i< pdi.size(); i++)
  {
    // look for penetrating and then not penetrating
    if (pdi_old[i].dist < 0.0 && pdi_new[i].dist > 0.0)
      unilateral_bracket.push_back(true);
    else if (pdi_old[i].dist > 0.0 && pdi_new[i].dist < 0.0)
      unilateral_bracket.push_back(true);
    else
      unilateral_bracket.push_back(false); 
    FILE_LOG(LOG_CONSTRAINT) << "Unilateral bracket " << i << ": " << unilateral_bracket[i] << std::endl;
  }

  // add brackets for the implicit constraints
  vector<bool> bilateral_bracket;
  for (unsigned i=0; i< C.size(); i++)
  {
    // look for sign change
    if (C[i] < 0.0 && C_old[i] > 0.0)
      bilateral_bracket.push_back(true);
    else if (C[i] > 0.0 && C_old[i] < 0.0)
      bilateral_bracket.push_back(true);
    else
      bilateral_bracket.push_back(false);
  }

  // NOTE: if there is no sign change at the endpoints, there may still be
  //       a sign change at some point in between, which means that after the
  //       root finding process, it is possible that two bodies could be
  //       interpenetrating that we didn't expect to be interpenetrating
  //       or that two bodies could be separated by an unexpected amount 

  // setup initial t
  double t = 1.0;

  // iterate over all unilateral brackets
  for (unsigned i=0; i< unilateral_bracket.size(); i++)
  {
    // verify that we check this bracket 
    if (!unilateral_bracket[i])
      continue;

    // call Ridder's method to determine new t
    double root = ridders_unilateral(0, t, pdi_old[i].dist, pdi_new[i].dist, i, dq, q, sim);
    if (root > 0.0 && root < 1.0)
      t = std::min(root, t);
  }

  FILE_LOG(LOG_CONSTRAINT) << " t determined from Ridder's method: " << t << std::endl; 

  // iterate over all bilateral brackets
  for (unsigned i=0; i< bilateral_bracket.size(); i++)
  {
    // verify that we check this bracket 
    if (!bilateral_bracket[i])
      continue;

    // call Ridder's method to determine new t
    double root = ridders_bilateral(0, t, C_old[i], C[i], i, dq, q, sim);
    if (root > 0.0 && root < 1.0)
      t = std::min(root, t);
  }

  // determine new qstar 
  qstar = dq;
  qstar *= t;
  qstar += q;

  // re-evaluate signed distances and bilateral constraints
  update_body_configurations(qstar, sim);
  sim->calc_pairwise_distances();
  evaluate_implicit_constraints(sim, C);

  // for values that aren't bracketed, further decrease t until there is
  // no increase in error
  const double BETA = 0.6;
  while (true)
  {
    // indicate that we do *not* keep looping unless otherwise specified
    bool stop = true;

    // check unilateral brackets first
    for (unsigned i=0; i< unilateral_bracket.size(); i++)
    {
      if (!unilateral_bracket[i])
        FILE_LOG(LOG_CONSTRAINT) << " old distance " << pdi_old[i].dist << " new distance: " << pdi[i].dist << " stop? " << (!(pdi[i].dist < 0.0 && 
          pdi_old[i].dist > pdi[i].dist)) << std::endl;
      if (!unilateral_bracket[i] && pdi[i].dist < 0.0 && 
          pdi_old[i].dist > pdi[i].dist)
      {
        stop = false;
        break;
      }
    }

    // now check bilateral constraints 
    if (stop)
    {
      double bilateral_cvio = 0.0;
      for (unsigned i=0; i< C.size(); i++)
          bilateral_cvio += sqr(C[i]);
      bilateral_cvio = std::sqrt(bilateral_cvio);

      FILE_LOG(LOG_SIMULATOR) << "Old bilateral constraint violation: " << old_bilateral_cvio << " New: " << bilateral_cvio << std::endl;

      // see whether we can still stop
      if (bilateral_cvio < bilateral_eps || bilateral_cvio < old_bilateral_cvio)
        break;
    }

    // update t
    t *= BETA;

    // if t is too small, quit
    if (t < MIN_T)
      return false;

    // determine new qstar if necessary 
    qstar = dq;
    qstar *= t;
    qstar += q;

    // re-evaluate signed distances and bilateral constraints
    update_body_configurations(qstar, sim);
    sim->calc_pairwise_distances();
    evaluate_implicit_constraints(sim, C);
  }

  FILE_LOG(LOG_SIMULATOR) << "t: " << t << std::endl;
  FILE_LOG(LOG_SIMULATOR) << "new q (from update_q): " << q << std::endl;

  // update q
  q = qstar;

  return true;
}

/// Gets the body configurations, placing them into q 
void ConstraintStabilization::get_body_configurations(VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{  
  const std::vector<shared_ptr<ControlledBody> >& bodies = sim->_bodies;
  unsigned NGC = 0;

  // resize the vector appropriately
  BOOST_FOREACH(shared_ptr<ControlledBody> cb, bodies)
  {
    shared_ptr<DynamicBodyd> body = dynamic_pointer_cast<DynamicBodyd>(cb);
    NGC += body->num_generalized_coordinates(DynamicBodyd::eEuler);
  }
  q.resize(NGC);

  // set the appropriate part of the vector
  unsigned start = 0;
  BOOST_FOREACH(shared_ptr<ControlledBody> cb, bodies)
  {
    shared_ptr<DynamicBodyd> body = dynamic_pointer_cast<DynamicBodyd>(cb);
    SharedVectorNd body_gcs = q.segment(start, start + body->num_generalized_coordinates(DynamicBodyd::eEuler));
    body->get_generalized_coordinates_euler(body_gcs);
    start += body->num_generalized_coordinates(DynamicBodyd::eEuler);
  }
}

/// Computes mapping from bodies to generalized coordinate indices 
void ConstraintStabilization::generate_body_index_map(std::map<shared_ptr<DynamicBodyd>, unsigned>& body_index_map, shared_ptr<ConstraintSimulator> sim)
{
  const std::vector<shared_ptr<ControlledBody> >& bodies = sim->_bodies;
  unsigned cur_index = 0;

  BOOST_FOREACH(shared_ptr<ControlledBody> cb, bodies){
    shared_ptr<DynamicBodyd> body = dynamic_pointer_cast<DynamicBodyd>(cb);
    body_index_map.insert(std::make_pair(body, cur_index));
    cur_index += body->num_generalized_coordinates(DynamicBodyd::eEuler);
  }
}

/// Updates the body configurations given q
void ConstraintStabilization::update_body_configurations(const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  const std::vector<shared_ptr<ControlledBody> >& bodies = sim->_bodies;
  unsigned last = 0;
  BOOST_FOREACH(shared_ptr<ControlledBody> cb, bodies){
    shared_ptr<DynamicBodyd> body = dynamic_pointer_cast<DynamicBodyd>(cb);
    unsigned ngc = body->num_generalized_coordinates(DynamicBodyd::eEuler);
    Ravelin::SharedConstVectorNd gc_shared = q.segment(last,last+ngc);
    body->set_generalized_coordinates_euler(gc_shared);
    last += ngc;
  }
}

/// sign function
static double sign(double x, double y)
{
  if (y > 0.0)
    return std::fabs(x);
  else 
    return -std::fabs(x);
}

/// Evaluates the function for root finding
double ConstraintStabilization::eval_unilateral(double t, unsigned i, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  static VectorNd qstar;

  // setup qstar
  qstar = dq;
  qstar *= t;
  qstar += q; 

  // update body configurations
  update_body_configurations(qstar, sim);

  // compute new pairwise distance information
  sim->calc_pairwise_distances(); 

  return sim->_pairwise_distances[i].dist;
}

/// Evaluates the function for root finding
double ConstraintStabilization::eval_bilateral(double t, unsigned i, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  static VectorNd qstar;
  std::vector<double> C;

  // setup qstar
  qstar = dq;
  qstar *= t;
  qstar += q; 

  // update body configurations
  update_body_configurations(qstar, sim);

  // compute new constraint evaluations
  evaluate_implicit_constraints(sim, C); 

  return C[i];
}


/// Ridders method for root finding
double ConstraintStabilization::ridders_unilateral(double x1, double x2, double fl, double fh, unsigned idx, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  const unsigned MAX_ITERATIONS = 25;
  const double TOL = 1e-4;
  double ans=std::numeric_limits<double>::max(),fm,fnew,s,xh,xl,xm,xnew;

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) 
  {
    xl=x1;
    xh=x2;
    ans = std::numeric_limits<double>::max();
    for (unsigned j=0;j< MAX_ITERATIONS;j++) 
    { 
      xm=0.5*(xl+xh);
      fm=eval_unilateral(xm, idx, dq, q, sim); // First of two function evaluations per iteration
      s=std::sqrt(fm*fm-fl*fh); 
      if (s == 0.0) 
        return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); // Updating formula
      ans=xnew;
      fnew=eval_unilateral(ans, idx, dq, q, sim);
      if (std::fabs(fnew) < TOL && fnew >= 0.0)
        return xnew;
      if (sign(fm,fnew) != fm) 
      {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fl,fnew) != fl) 
      { 
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fh,fnew) != fh) 
      {
        xl=ans; // second of two function evaluations per iteration.
        fl=fnew; // Bookkeeping to keep the root bracketed on next iteration.
      } 
      else 
        assert(false);
    }

    // if here, maximum iterations have been exceeded
  }
  else 
  {
     if (fl == 0.0) 
       return x1;
     if (fh == 0.0) 
       return x2;

     throw std::runtime_error("Root was not bracketed");
  }

  return 0.0; // Never get here. 
}

/// Ridders method for root finding
double ConstraintStabilization::ridders_bilateral(double x1, double x2, double fl, double fh, unsigned idx, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  const unsigned MAX_ITERATIONS = 25;
  const double TOL = 1e-6;
  double ans=std::numeric_limits<double>::max(),fm,fnew,s,xh,xl,xm,xnew;

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) 
  {
    xl=x1;
    xh=x2;
    ans = std::numeric_limits<double>::max();
    for (unsigned j=0;j< MAX_ITERATIONS;j++) 
    { 
      xm=0.5*(xl+xh);
      fm=eval_bilateral(xm, idx, dq, q, sim); // First of two function evaluations per iteration
      s=std::sqrt(fm*fm-fl*fh); 
      if (s == 0.0) 
        return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); // Updating formula
      ans=xnew;
      fnew=eval_bilateral(ans, idx, dq, q, sim);
      if (std::fabs(fnew) < TOL)
        return ans; 
      if (sign(fm,fnew) != fm) 
      {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fl,fnew) != fl) 
      { 
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fh,fnew) != fh) 
      {
        xl=ans; // second of two function evaluations per iteration.
        fl=fnew; // Bookkeeping to keep the root bracketed on next iteration.
      } 
      else 
        assert(false);
    }

    // if here, maximum iterations have been exceeded
  }
  else 
  {
     if (fl == 0.0) 
       return x1;
     if (fh == 0.0) 
       return x2;

     throw std::runtime_error("Root was not bracketed");
  }

  return 0.0; // Never get here. 
}
