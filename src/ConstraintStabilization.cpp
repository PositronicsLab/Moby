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
#include <utility>
#include <boost/algorithm/minmax_element.hpp>
#include <Moby/Types.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/LCPSolverException.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/ConstraintStabilization.h>

#ifdef USE_QPOASES
#include <Moby/qpOASES.h>
#endif

using namespace Ravelin;
using namespace Moby;
using std::pair;
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
  // set maximum iterations to infinity, by default 
  max_iterations = std::numeric_limits<unsigned>::max();

  // set unilateral tolerance to negative NEAR_ZERO by default
  eps = NEAR_ZERO;

  // set bilateral tolerance
  bilateral_eps = 1e-6;

  // set the inequality tolerance
  inequality_tolerance = -1.0;
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

/// Get the maximum amount of unilateral constraint violation
double ConstraintStabilization::evaluate_unilateral_constraints(shared_ptr<ConstraintSimulator> sim, vector<double>& uC)
{
  // set violation to infinite initially
  double vio = std::numeric_limits<double>::max();

  // clear uC
  uC.clear();

  // calculate the pairwise distances
  sim->calc_pairwise_distances();

  // get the maximum amount of unilateral constraint violation
  const vector<PairwiseDistInfo>& pdi = sim->_pairwise_distances;
  for (unsigned i=0; i< pdi.size(); i++)
  {
    uC.push_back(pdi[i].dist);
    vio = std::min(vio, uC.back());
  }

  // look at all articulated bodies
  const vector<ControlledBodyPtr>& bodies = sim->get_dynamic_bodies();
  for (unsigned i=0; i< bodies.size(); i++)
  {
    shared_ptr<RCArticulatedBodyd> rcab = dynamic_pointer_cast<RCArticulatedBodyd>(bodies[i]);
    if (rcab)
    {
      const std::vector<shared_ptr<Jointd> >& joints = rcab->get_joints();
      for (unsigned j=0; j< joints.size(); j++)
      {
        shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]);
        const VectorNd& tare = joint->get_q_tare(); 
        for (unsigned k=0; k< joint->num_dof(); k++)
        {
          uC.push_back(joint->hilimit[k] - joint->q[k] - tare[k]);
          vio = std::min(vio, uC.back());
          uC.push_back(joint->q[k] + tare[k] - joint->lolimit[k]);
          vio = std::min(vio, uC.back());
        }
      }
    }
  }

  // also look at joint violations for all joints in the simulator
  if (!sim->implicit_joints.empty())
  {
    const std::vector<JointPtr>& joints = sim->implicit_joints;
    for (unsigned j=0; j< joints.size(); j++)
    {
      shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]);
      const VectorNd& tare = joint->get_q_tare(); 
      for (unsigned k=0; k< joint->num_dof(); k++)
      {
        uC.push_back(joint->hilimit[k] - joint->q[k] - tare[k]);
        vio = std::min(vio, uC.back());
        uC.push_back(joint->q[k] + tare[k] - joint->lolimit[k]);
        vio = std::min(vio, uC.back());
      }
    }
  }

  return vio;
}

/// Evaluates implicit bilateral constraints and returns the most significant violation
double ConstraintStabilization::evaluate_bilateral_constraints(shared_ptr<ConstraintSimulator> sim, vector<double>& C)
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
  std::vector<VectorNd> qd_save;
  std::vector<double> C, uC;
  std::map<shared_ptr<DynamicBodyd>, unsigned> body_index_map;

  // look for no constraint stabilization
  if (max_iterations == 0)
    return;

  // apply broad phase
  sim->broad_phase(0.0);
 
  // save the generalized velocities
  save_velocities(sim, qd_save);

  // get the body configurations and setup the body index map
  get_body_configurations(q, sim);
  generate_body_index_map(body_index_map, sim);

  // evaluate the bilateral constraints
  double max_bvio = evaluate_bilateral_constraints(sim, C);

  // see whether any pairwise distances are below epsilon
  double max_uvio = evaluate_unilateral_constraints(sim, uC);

  FILE_LOG(LOG_SIMULATOR) <<"maximum unilateral constraint violation (before stabilization loop): "<< max_uvio <<std::endl;
  FILE_LOG(LOG_SIMULATOR) <<"maximum bilateral constraint violation (before stabilization loop): "<< max_bvio <<std::endl;

  unsigned iterations = 0;
  while (max_uvio < eps || max_bvio > bilateral_eps)
  {
    // look for maximum iterations
    if (iterations == max_iterations)
    {
      FILE_LOG(LOG_SIMULATOR) << " -- maximum number of iterations reached" << std::endl;
      FILE_LOG(LOG_SIMULATOR) << " -- failed to effectively finish the constraint stabilization process!" << std::endl;
      break;
    }

    FILE_LOG(LOG_SIMULATOR) <<"maximum unilateral constraint violation: "<< max_uvio <<std::endl;
    FILE_LOG(LOG_SIMULATOR) <<"maximum bilateral constraint violation: "<< max_bvio <<std::endl;

    // zero body velocities first (we only want to change positions based on
    // our updates)
    for (unsigned i=0; i< sim->_bodies.size(); i++)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(sim->_bodies[i]);
      db->get_generalized_velocity(DynamicBodyd::eSpatial, v);
      v.set_zero();
      db->set_generalized_velocity(DynamicBodyd::eSpatial, v);
    }

    // determine islands and super bodies in the islands
    vector<vector<Constraint*> > pd;
    vector<vector<shared_ptr<DynamicBodyd> > > super_bodies;
    compute_problem_data(pd, super_bodies, sim);

    // determine dq's
    dq.set_zero(q.size());
    for (unsigned i=0; i< pd.size(); i++)
      determine_dq(pd[i], super_bodies[i], dq, body_index_map);
    FILE_LOG(LOG_SIMULATOR) << "dq: " << dq << std::endl;

    // TODO: compute maximum step size for each island

    // determine s and update q; NOTE: update q computes the pairwise distances 
    if (!update_q(dq, q, sim))
    {
      FILE_LOG(LOG_SIMULATOR) << " -- failed to effectively finish the constraint stabilization process!" << std::endl;
      break;
    }

    // save last violations
    double last_uvio = max_uvio;
    double last_bvio = max_bvio;

    // update minimum distance
    max_uvio = evaluate_unilateral_constraints(sim, uC);

    // update constraint violation
    max_bvio = evaluate_bilateral_constraints(sim, C);

    if (last_uvio >= max_uvio && max_bvio >= last_bvio)
    {
      FILE_LOG(LOG_SIMULATOR) << " -- could not stabilize constraints further" << std::endl;
      break;
    } 

    iterations++;
  }

  // restore the generalized velocities
  restore_velocities(sim, qd_save);

  // after stabilization, velocities may be in an impacting state; correct
//  sim->find_unilateral_constraints(sim->contact_dist_thresh);
//  sim->calc_impacting_unilateral_constraint_forces(-1.0); 

  FILE_LOG(LOG_SIMULATOR) << iterations << " iterations required" << std::endl;
  FILE_LOG(LOG_SIMULATOR) <<"=====constraint stabilization end ======" << std::endl;
}

/// Adds constraints for bilateral joints
void ConstraintStabilization::add_joint_constraints(const vector<ControlledBodyPtr>& bodies, std::vector<Constraint>& constraints, shared_ptr<ConstraintSimulator> sim)
{
  const double INF = std::numeric_limits<double>::max();

  // loop through all bodies looking for constraints
  for (unsigned i=0; i< bodies.size(); i++)
  {
    shared_ptr<RCArticulatedBodyd> rcab = dynamic_pointer_cast<RCArticulatedBodyd>(bodies[i]);
    if (rcab)
    {
      const std::vector<shared_ptr<Jointd> >& joints = rcab->get_joints();
      for (unsigned j=0; j< joints.size(); j++)
      {
        shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]); 
        if (joint->get_constraint_type() == Jointd::eImplicit)
        {
          Constraint e;
          e.constraint_type = Constraint::eImplicitJoint;
          e.implicit_joint = joint;
          constraints.push_back(e);
        }
      }
    }
  }

  // also look at joint violations for all implicit joints in the simulator
  const std::vector<JointPtr>& joints = sim->implicit_joints;
  for (unsigned j=0; j< joints.size(); j++)
  {
    shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]);
    Constraint e;
    e.constraint_type = Constraint::eImplicitJoint;
    e.implicit_joint = joint;
    constraints.push_back(e);
  }
}


/// Adds unilateral constraints for joint limits in an articulated body
void ConstraintStabilization::add_limit_constraints(const vector<ControlledBodyPtr>& bodies, std::vector<Constraint>& constraints, shared_ptr<ConstraintSimulator> sim)
{
  const double INF = std::numeric_limits<double>::max();

  // loop through all bodies looking for constraints
  for (unsigned i=0; i< bodies.size(); i++)
  {
    shared_ptr<RCArticulatedBodyd> rcab = dynamic_pointer_cast<RCArticulatedBodyd>(bodies[i]);
    if (rcab)
    {
      const std::vector<shared_ptr<Jointd> >& joints = rcab->get_joints();
      for (unsigned j=0; j< joints.size(); j++)
      {
        shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]); 
        for (unsigned k=0; k< joint->num_dof(); k++)
        {
          const VectorNd& tare = joint->get_q_tare(); 

          // only add an upper limit constraint if the upper limit < inf
          if (joint->hilimit[k] < INF)
          {
            Constraint e;
            e.constraint_type = Constraint::eLimit;
            e.limit_joint = joint;
            e.limit_dof = k;
            e.limit_epsilon = 0.0;
            e.limit_upper = true;
            e.limit_stiffness = 1.0;
            e.signed_distance = joint->hilimit[k] - joint->q[k] - tare[k];
            constraints.push_back(e);
          }

          // only add a lower limit constraint if the lower limit > -inf
          if (joint->lolimit[k] > -INF)
          {
            Constraint e;
            e.constraint_type = Constraint::eLimit;
            e.limit_joint = joint;
            e.limit_dof = k;
            e.limit_epsilon = 0.0;
            e.limit_upper = false;
            e.limit_stiffness = 1.0;
            e.signed_distance = joint->q[k] + tare[k] - joint->lolimit[k];
            constraints.push_back(e);
          }
        }
      }
    }
  }

  // also look at joint violations for all implicit joints in the simulator
  if (!sim->implicit_joints.empty())
  {
    const std::vector<JointPtr>& joints = sim->implicit_joints;
    for (unsigned j=0; j< joints.size(); j++)
    {
      shared_ptr<Joint> joint = dynamic_pointer_cast<Joint>(joints[j]);
      const VectorNd& tare = joint->get_q_tare(); 
      for (unsigned k=0; k< joint->num_dof(); k++)
      {
        // only add an upper limit constraint if the upper limit < inf
        if (joint->hilimit[k] < INF)
        {
          Constraint e;
          e.constraint_type = Constraint::eLimit;
          e.limit_joint = joint;
          e.limit_dof = k;
          e.limit_epsilon = 0.0;
          e.limit_upper = true;
          e.limit_stiffness = 1.0;
          e.signed_distance = joint->hilimit[k] - joint->q[k] - tare[k];
          constraints.push_back(e);
        }

        // only add a lower limit constraint if the lower limit > -inf
        if (joint->lolimit[k] > -INF)
        {
          Constraint e;
          e.constraint_type = Constraint::eLimit;
          e.limit_joint = joint;
          e.limit_dof = k;
          e.limit_epsilon = 0.0;
          e.limit_upper = false;
          e.limit_stiffness = 1.0;
          e.signed_distance = joint->q[k] + tare[k] - joint->lolimit[k];
          constraints.push_back(e);
        }
      }
    }
  }
}

/// Adds contact constraints from a pair of bodies
void ConstraintStabilization::add_contact_constraints(std::vector<Constraint>& constraints, CollisionGeometryPtr cg1, CollisionGeometryPtr cg2, shared_ptr<ConstraintSimulator> sim)
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
    Constraint uc = CollisionDetection::create_contact(cg1, cg2, p1_g, normal, dist);

    // make the contact frictionless
    uc.contact_mu_coulomb = 0.0;

    // maximize the contact stiffness 
    uc.contact_stiffness = 1.0;

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

    // make the contacts frictionless and rigid
    for (unsigned i=OLD_CONSTRAINTS_SZ; i< constraints.size(); i++)
    {
      constraints[i].contact_mu_coulomb = 0.0;
      constraints[i].contact_stiffness = 1.0;
      FILE_LOG(LOG_CONSTRAINT) << constraints[i] << std::endl;
    }
  }
}

/// Computes the constraint data
void ConstraintStabilization::compute_problem_data(vector<vector<Constraint*> >& pd, vector<vector<shared_ptr<DynamicBodyd> > >& super_bodies, shared_ptr<ConstraintSimulator> sim)
{
  VectorNd tmpv;

  // clear the constraints vector
  constraints.clear();

  // clear the problem data vector 
  pd.clear();

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

  // add limit constraints
  add_limit_constraints(bodies, constraints, sim);

  // add joint constraints
  add_joint_constraints(bodies, constraints, sim);

  //FILE_LOG(LOG_SIMULATOR) << "constraints added" << std::endl;
  // 2) for each articulated body, add as many UnilateralConstraint objects as
  //    there are joints at their limits


  // 3) for each pair of bodies in interpenetrating contact, add a single
  //    point of contact at the deepest interpenetrating point with normal
  //    in the direction of the signed distance function. 

  // find islands
  list<vector<shared_ptr<DynamicBodyd> > > remaining_islands;
  list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > > islands;
  Constraint::determine_connected_constraints(constraints, islands, remaining_islands);

  // construct problem data
  for (list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > >::const_iterator island_iter = islands.begin(); island_iter != islands.end(); island_iter++)
    pd.push_back(vector<Constraint*>(island_iter->first.begin(), island_iter->first.end())); 

  // get super bodies in the islands
  super_bodies.resize(pd.size());
  for (unsigned i=0; i< pd.size(); i++)
  {
    super_bodies[i].clear();
    ImpactConstraintHandler::get_super_bodies(pd[i].begin(), pd[i].end(), std::back_inserter(super_bodies[i]));
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
      else
        return sb;
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

/// Computes deltaq by solving a linear complementarity problem
void ConstraintStabilization::determine_dq(vector<Constraint*>& pd, const vector<shared_ptr<DynamicBodyd> >& super_bodies, VectorNd& dqm, const std::map<shared_ptr<DynamicBodyd>, unsigned>& body_index_map)
{
  VectorNd z, dq_sub;

  // determine the total number of variables, number of inequality constraints,
  // and number of equality constraints
  const unsigned N_VARS = ImpactConstraintHandler::num_variables(pd);
  const unsigned N_INEQ_CONSTRAINTS = ImpactConstraintHandler::num_inequality_constraints(pd);
  const unsigned N_EQ_CONSTRAINTS = ImpactConstraintHandler::num_equality_constraints(pd);

  // compute the linear term vector and quadratic matrix terms
  ImpactConstraintHandler::compute_quadratic_matrix(pd, N_VARS, _H);
  ImpactConstraintHandler::compute_linear_term(pd, N_VARS, 0.0, _c);

  // everything is happening instantaneously; nevertheless, we will set
  // inv_dt to 1.0 to incorporate constraint violations
  const double INV_DT = 1.0;
  const double INF = std::numeric_limits<double>::max();

  // compute constraint terms
  ImpactConstraintHandler::compute_inequality_terms(pd, _H, N_INEQ_CONSTRAINTS, INV_DT, _M, _q);
  ImpactConstraintHandler::compute_equality_terms(pd, _H, N_EQ_CONSTRAINTS, INV_DT, _A, _b);

  // modify q terms
  for (unsigned i=0; i< _q.size(); i++)
    _q[i] += std::fabs(eps) + NEAR_ZERO;

  // see whether an MLCP is necessary
  if (N_EQ_CONSTRAINTS == 0)
  {
    // no MLCP necessary
    // solve using an LCP formulation
    _MM.set_zero(N_VARS + N_INEQ_CONSTRAINTS, N_VARS + N_INEQ_CONSTRAINTS);
    _qq.resize(_MM.rows());

    // setup MM and qq
    _MM.block(0, N_VARS, 0, N_VARS) = _H;
    _MM.block(N_VARS, _MM.rows(), 0, N_VARS) = _M;
    SharedMatrixNd MT = _MM.block(0, N_VARS, N_VARS, _MM.rows());
    MatrixNd::transpose(_M, MT);
    MT.negate();
    _qq.segment(0, N_VARS) = _c;
    _qq.segment(N_VARS, _qq.size()) = _q;
    _qq.segment(N_VARS, _qq.size()).negate();

    // solve using the fast regularized solver first, then fall back to the
    // Lemke solver
    if (!_lcp.lcp_fast_regularized(_MM, _qq, z, -20, 4, -8))
    {
      // Lemke does not like warm starting
      z.set_zero();

      if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
        throw LCPSolverException();
    }

    FILE_LOG(LOG_SIMULATOR) << "# of constraints in determine_dq(): " << _qq.size() << std::endl;
    FILE_LOG(LOG_SIMULATOR) << "MM: " << std::endl << _MM;
    FILE_LOG(LOG_SIMULATOR) << "qq: " << _qq << std::endl;
    FILE_LOG(LOG_SIMULATOR) << "zz: " << z << std::endl;
  } 
  else
  {
    #ifdef USE_QPOASES
    static qpOASES qp;

    // setup the lower and upper bounds variables
    _lb.resize(N_VARS);
    _ub.resize(N_VARS);
    for (unsigned i=0, j=0; i< pd.size(); i++)
      for (unsigned k=0; k< pd[i]->num_variables(); k++)
      {
        _lb[j] = pd[i]->get_lower_variable_limit(k);
        _ub[j] = pd[i]->get_upper_variable_limit(k);
        j++;
      } 

    // add tolerance to inequality constraints
    double TOL = (inequality_tolerance >= 0.0) ? inequality_tolerance : (_H.rows() + _q.rows()) * (_H.norm_inf() + _M.norm_inf()) * std::numeric_limits<double>::epsilon();
    for (ColumnIteratord iter = _q.begin(); iter != _q.end(); iter++)
      *iter -= TOL;

    // add tolerance to lower bounds (if necessary)
    for (ColumnIteratord iter = _lb.begin(); iter != _lb.end(); iter++)
      if (*iter > -INF)
        *iter -= TOL;

    // add tolerance to upper bounds (if necessary)
    for (ColumnIteratord iter = _ub.begin(); iter != _ub.end(); iter++)
      if (*iter < INF)
        *iter += TOL;

    // no slackable equality constraints; try solving QP first w/o any
    // tolerance in the constraints
    bool result = qp.qp_activeset(_H, _c, _lb, _ub, _M, _q, _A, _b, z);
    if (!result)
    {
      FILE_LOG(LOG_CONSTRAINT) << "QP activeset solution failed without tolerance in the constraints" << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << " -- solving LP to find feasible solution" << std::endl;

      // determine new number of variables (one for inequality constraints,
      // two for equality constraints)
      unsigned NEW_VARS = N_VARS + 1 + N_EQ_CONSTRAINTS * 2;

      // augment M
      _Maug.resize(_M.rows(), NEW_VARS);
      _Maug.block(0, _M.rows(), 0, N_VARS) = _M;
      _Maug.block(0, _M.rows(), N_VARS, NEW_VARS).set_zero();
      _Maug.column(N_VARS).set_one();

      // augment A
      _Aaug.resize(_A.rows(), NEW_VARS);
      _Aaug.block(0, _A.rows(), 0, N_VARS) = _A;
      _Aaug.block(0, _A.rows(), N_VARS, NEW_VARS).set_zero();
      for (unsigned i=0, j=N_VARS+1; i< _A.rows(); i++)
      {
        _Aaug(i,j++) = 1.0;
        _Aaug(i,j++) = -1.0;
      }

      // augment lb
      _lbaug.resize(NEW_VARS);
      _lbaug.segment(0, _lb.size()) = _lb;
      _lbaug.segment(_lb.size(), NEW_VARS).set_zero();

      // augment ub
      _ubaug.resize(NEW_VARS);
      _ubaug.segment(0, _ub.size()) = _ub;
      _ubaug.segment(_ub.size(), NEW_VARS).set_one() *= INF;

      // create a c for the QP
      _caug.resize(NEW_VARS);
      _caug.segment(0, N_VARS).set_zero();
      _caug.segment(N_VARS, NEW_VARS).set_one();

      // solve the LP, using zero for z
      z.set_zero(NEW_VARS);
      result = LP::lp_simplex(_caug, _Maug, _q, _Aaug, _b, _lbaug, _ubaug, z);
      assert(result);
      FILE_LOG(LOG_CONSTRAINT) << " -- LP solution: " << z << std::endl;

      // modify c for re-solving the QP
      _caug.segment(0, N_VARS) = _c;
      _caug.segment(N_VARS, NEW_VARS).set_zero();

      // modify H for re-solving the QP
      _Haug.resize(NEW_VARS, NEW_VARS);
      _Haug.block(0, N_VARS, 0, N_VARS) = _H;
      _Haug.block(N_VARS, NEW_VARS, 0, N_VARS).set_zero();
      _Haug.block(0, N_VARS, N_VARS, NEW_VARS).set_zero();
      _Haug.block(N_VARS, NEW_VARS, N_VARS, NEW_VARS).set_zero();

      // modify lb and ub for the new variables
      for (unsigned i=N_VARS; i< NEW_VARS; i++)
        _lbaug[i] = _ubaug[i] = z[i];

      // resolve the QP
      result = qp.qp_activeset(_Haug, _caug, _lbaug, _ubaug, _Maug, _q, _Aaug, _b, z);
      FILE_LOG(LOG_CONSTRAINT) << "result of QP activeset with constraints: " << result << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << "robust QP activeset solution: " << z << std::endl;
    }
    else
      FILE_LOG(LOG_CONSTRAINT) << "QP activeset solution found" << std::endl;
    #else
    std::cerr << "Moby not built with qpOASES support; unable to solve QP" << std::endl;
    #endif
  }

  // update velocities
  ImpactConstraintHandler::apply_impulses(pd, z);

  // populate dq 
  for(unsigned i = 0; i< super_bodies.size(); i++)
  {
    assert(body_index_map.find(super_bodies[i]) != body_index_map.end());
    unsigned start = (body_index_map.find(super_bodies[i]))->second;
    unsigned coord_num = super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eEuler);
    super_bodies[i]->get_generalized_velocity(DynamicBodyd::eEuler, dq_sub);  
    dqm.segment(start, start+coord_num) = dq_sub;
  }
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
  vector<double> C, C_old, uC, uC_old;
  const double MIN_T = NEAR_ZERO;

  // evaluate the constraints
  evaluate_unilateral_constraints(sim, uC_old);
  evaluate_bilateral_constraints(sim, C_old);

  // compute old bilateral constraint violations
  double old_bilateral_cvio = 0.0;
  for (unsigned i=0; i< C_old.size(); i++)
    old_bilateral_cvio += sqr(C_old[i]);
  old_bilateral_cvio = std::sqrt(old_bilateral_cvio);

  FILE_LOG(LOG_CONSTRAINT) << "...about to compute unilateral brackets" << std::endl;
  // find the pairwise distances and implicit constraint evaluations at q + dq
  qstar = dq;
  qstar += q;
  update_body_configurations(qstar, sim);
  evaluate_unilateral_constraints(sim, uC); 
  evaluate_bilateral_constraints(sim, C);

  // we may have to find roots for all pairwise distances
  vector<bool> unilateral_bracket;
  for (unsigned i=0; i< uC.size(); i++)
  {
    // look for violating and then not violating 
    if (uC_old[i] < 0.0 && uC[i] > 0.0)
      unilateral_bracket.push_back(true);
    else if (uC_old[i] > 0.0 && uC[i] < 0.0)
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
    double root = ridders_unilateral(0, t, uC_old[i], uC[i], i, dq, q, sim);
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
  evaluate_unilateral_constraints(sim, uC); 
  evaluate_bilateral_constraints(sim, C);

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
        FILE_LOG(LOG_CONSTRAINT) << " old distance " << uC_old[i] << " new distance: " << uC[i] << " stop? " << (!(uC[i] < 0.0 && 
          uC_old[i] > uC[i])) << std::endl;
      if (!unilateral_bracket[i] && uC[i] < 0.0 && 
          uC_old[i] > uC[i])
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
    evaluate_bilateral_constraints(sim, C);
    evaluate_unilateral_constraints(sim, uC);
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

// TODO: update this
/// Evaluates the function for root finding
double ConstraintStabilization::eval_unilateral(double t, unsigned i, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  static VectorNd qstar;
  std::vector<double> uC;

  // setup qstar
  qstar = dq;
  qstar *= t;
  qstar += q; 

  // update body configurations
  update_body_configurations(qstar, sim);

  // compute new pairwise distance information
  evaluate_unilateral_constraints(sim, uC); 

  return uC[i];
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
  evaluate_bilateral_constraints(sim, C); 

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
