/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <algorithm>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/permute.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/SignedDistDot.h>
#ifdef HAVE_IPOPT
#include <Moby/NQP_IPOPT.h>
#include <Moby/LCP_IPOPT.h>
#endif

using namespace Ravelin;
using namespace Moby;
using std::list;
using boost::shared_ptr;
using std::vector;
using std::map;
using std::endl;
using std::cerr;
using std::pair;
using std::min_element;
using boost::dynamic_pointer_cast;

/// Sets up the default parameters for the impact event handler
ImpactConstraintHandler::ImpactConstraintHandler()
{
  ip_max_iterations = 100;
  ip_eps = 1e-6;
  use_ip_solver = false;

  // initialize IPOPT, if present
  #ifdef HAVE_IPOPT
  _app.Options()->SetNumericValue("tol", 1e-7);
  _app.Options()->SetStringValue("mu_strategy", "adaptive");
  _app.Options()->SetStringValue("output_file", "ipopt.out");
  #ifndef __APPLE__
  _app.RethrowNonIpoptException(true);
  #endif

  Ipopt::ApplicationReturnStatus status = _app.Initialize();
  if (status != Ipopt::Solve_Succeeded)
    throw std::runtime_error("Ipopt unable to initialize!");

  // setup the nonlinear IP solver
  _ipsolver = Ipopt::SmartPtr<NQP_IPOPT>(new NQP_IPOPT);
  _lcpsolver = Ipopt::SmartPtr<LCP_IPOPT>(new LCP_IPOPT);
  #endif
}

// Processes impacts
/**
 * \param constraints the vector of constraints
 */
void ImpactConstraintHandler::process_constraints(const vector<UnilateralConstraint>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // apply the method to all contacts
  apply_model(constraints);

  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of constraints
/**
 * \param constraints a set of constraints
 */
void ImpactConstraintHandler::apply_model(const vector<UnilateralConstraint>& constraints)
{
  const double INF = std::numeric_limits<double>::max();
  list<UnilateralConstraint*> impacting;
  VectorNd dv, v, f, lambda;

  // **********************************************************
  // determine sets of connected constraints
  // **********************************************************
  list<vector<shared_ptr<DynamicBodyd> > > remaining_islands;
  list<list<UnilateralConstraint*> > groups;
  UnilateralConstraint::determine_connected_constraints(constraints, _simulator->implicit_joints, groups, remaining_islands);
  UnilateralConstraint::remove_inactive_groups(groups);

  // **********************************************************
  // do method for each connected set
  // **********************************************************
  for (list<list<UnilateralConstraint*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
      // copy the list of constraints
      list<UnilateralConstraint*> rconstraints = *i;

      FILE_LOG(LOG_CONSTRAINT) << " -- pre-constraint velocity (all constraints): " << std::endl;
      for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;

      // look to see whether all contact constraints have zero or infinite Coulomb friction
      bool all_inf = true, all_frictionless = true;
      BOOST_FOREACH(UnilateralConstraint* e, rconstraints)
        if (e->constraint_type == UnilateralConstraint::eContact)
        {
          if (e->contact_mu_coulomb < 1e2)
            all_inf = false;
          if (e->contact_mu_coulomb > 0.0)
            all_frictionless = false;
        }

      // apply model to the reduced contacts
      if (all_inf)
        apply_no_slip_model_to_connected_constraints(rconstraints);
// TODO: fix viscous model- seems to be a bug in it
//      else if (all_frictionless)
//        apply_visc_friction_model_to_connected_constraints(rconstraints);
  #ifdef USE_AP_MODEL
      else {
        apply_ap_model_to_connected_constraints(rconstraints);
      }
  #else
      else
        apply_model_to_connected_constraints(rconstraints);
  #endif

      FILE_LOG(LOG_CONSTRAINT) << " -- post-constraint velocity (all constraints): " << std::endl;
      for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;
  }

  // determine whether there are any impacting constraints remaining
  for (list<list<UnilateralConstraint*> >::const_iterator i = groups.begin(); i != groups.end(); i++)
    for (list<UnilateralConstraint*>::const_iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->determine_constraint_class() == UnilateralConstraint::eNegative)
        impacting.push_back(*j);

  // if there are any constraints still impacting, throw an exception
  if (!impacting.empty())
    throw ImpactToleranceException(impacting);

// NOTE: this should already be handled for dynamics calculations
/*
  // process islands composed completely of bilateral constraints
  BOOST_FOREACH(vector<shared_ptr<DynamicBodyd> >& island, remaining_islands)
  {
    // sort the island so we can search it
    std::sort(island.begin(), island.end());

    // setup a set of implicit joints
    vector<JointPtr> island_ijoints;

    // get the implicit joints in the island
    const vector<JointPtr>& implicit_joints = _simulator->implicit_joints;
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

    // get the total number of generalized coordinates for the island
    const unsigned NGC_TOTAL = _simulator->num_generalized_coordinates(island);

    // setup f
    f.set_zero(NGC_TOTAL);

    // compute change in velocity and constraint forces
    _simulator->solve(island, island_ijoints, f, dv, lambda);

    // set new velocities 
    for (unsigned i=0, gc_index = 0; i< island.size(); i++)
    {
      const unsigned NGC = island[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
      SharedConstVectorNd dv_sub = dv.segment(gc_index, gc_index + NGC);
      island[i]->get_generalized_velocity(DynamicBodyd::eSpatial, v);
      v += dv_sub;
      island[i]->set_generalized_velocity(DynamicBodyd::eSpatial, v);
      gc_index += NGC;
      FILE_LOG(LOG_DYNAMICS) << "new velocity for " << island[i]->body_id << ": " << v << std::endl;
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
*/
}

/**
 * Applies purely viscous friction model to connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // set the simulator
  _epd.simulator = _simulator;

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd);

  // clear all impulses
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // solve the viscous friction model
  apply_visc_friction_model(_epd);

  // determine velocities due to impulse application
  update_constraint_velocities_from_impulses(_epd);

  // get the constraint violation before applying impulses
  double minv = calc_min_constraint_velocity(_epd);

  // apply restitution
  if (apply_restitution(_epd))
  {
    // determine velocities due to impulse application
    update_constraint_velocities_from_impulses(_epd);

    // check to see whether we need to solve another impact problem
    double minv_plus = calc_min_constraint_velocity(_epd);
    FILE_LOG(LOG_CONSTRAINT) << "Applying restitution" << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  compression v+ minimum: " << minv << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  restitution v+ minimum: " << minv_plus << std::endl;
    if (minv_plus < 0.0 && minv_plus < minv - NEAR_ZERO)
    {
      // need to solve another impact problem
      apply_visc_friction_model(_epd);
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints() exiting" << endl;
}

/**
 * Applies no slip friction model to connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_no_slip_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_no_slip_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // set the simulator
  _epd.simulator = _simulator;

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd);

  // clear all impulses
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // solve the no slip model
  apply_no_slip_model(_epd);

  // determine velocities due to impulse application
  update_constraint_velocities_from_impulses(_epd);

  // get the constraint violation before applying impulses
  double minv = calc_min_constraint_velocity(_epd);

  // apply restitution
  if (apply_restitution(_epd))
  {
    // determine velocities due to impulse application
    update_constraint_velocities_from_impulses(_epd);

    // check to see whether we need to solve another impact problem
    double minv_plus = calc_min_constraint_velocity(_epd);
    FILE_LOG(LOG_CONSTRAINT) << "Applying restitution" << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  compression v+ minimum: " << minv << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  restitution v+ minimum: " << minv_plus << std::endl;
    if (minv_plus < 0.0 && minv_plus < minv - NEAR_ZERO)
    {
      // need to solve another impact problem
      apply_no_slip_model(_epd);
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_no_slip_model_to_connected_constraints() exiting" << endl;
}

/// Updates determined impulses in UnilateralConstraintProblemData based on a QP/NQP solution
void ImpactConstraintHandler::update_from_stacked(UnilateralConstraintProblemData& q, const VectorNd& z)
{
  // save impulses in q
  if (use_qp_solver(q))
    q.update_from_stacked_qp(z);
  else
    q.update_from_stacked_nqp(z);

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save contact impulses
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = q.contact_constraints[i]->contact_normal * q.cn[i];

    // look whether friction is accounted for
    if (q.cs.size() > 0)
    {
      j += q.contact_constraints[i]->contact_tan1 * q.cs[i];
      j += q.contact_constraints[i]->contact_tan2 * q.ct[i];
    }

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += Pose3d::transform(GLOBAL, jx);
  }

  // save limit impulses
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    double limit_impulse = (q.limit_constraints[i]->limit_upper) ? -q.l[i] : q.l[i];
    q.limit_constraints[i]->limit_impulse += limit_impulse;
  }

  // get the change in velocity
  VectorNd dv, tmpv;
  q.X_CnT.mult(q.cn, dv);
  dv += q.X_CsT.mult(q.cs, tmpv);
  dv += q.X_CtT.mult(q.ct, tmpv);
  dv += q.X_LT.mult(q.l, tmpv);

  // compute lambda here using
  // | M  J' | | dv     | = | 0 | 
  // | J  0  | | lambda |   | -J*v |
  // M*dv + J'*lambda = 0
  // dv = -inv(M)*J'*lambda
  // -J*inv(M)*J'*lambda = -J*v
  q.lambda = q.Jx_v;
  FILE_LOG(LOG_CONSTRAINT) << "J*dv (constraint velocities): " << dv << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "bilateral constraint forces: " << q.lambda << std::endl;
MatrixNd tmp;
  FILE_LOG(LOG_CONSTRAINT) << "Jx " << std::endl << q.Jfull.to_dense(tmp);
  FILE_LOG(LOG_CONSTRAINT) << "inv(M)*Jx' " << std::endl << q.iM_JxT;
  FILE_LOG(LOG_CONSTRAINT) << "J*inv(M)*Jx' " << std::endl << q.Jx_iM_JxT;

  // solve for lambda
  MatrixNd tmpM = q.Jx_iM_JxT;
  LinAlgd::factor_chol(tmpM);
  LinAlgd::solve_chol_fast(tmpM, q.lambda);

  // update dv
  dv -= q.iM_JxT.mult(q.lambda, tmpv);
  FILE_LOG(LOG_CONSTRAINT) << "inv(M)*J'*lambda: " << tmpv << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "change in velocity after incorporating joint constraints: " << dv << std::endl;

  // compute J*dv
  if (LOGGING(LOG_CONSTRAINT))
  {
    VectorNd tmpv;
    q.Jfull.mult(dv, tmpv);
    FILE_LOG(LOG_CONSTRAINT) << "predicted constraint velocity: " << tmpv << std::endl; 
  }

  // update the bodies' velocities
  update_generalized_velocities(q, dv);

  // setup the full lambda
  VectorNd lambda;
  lambda.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  lambda.set(q.active, q.lambda);

  // push lambda back into individual joint lambdas
  for (unsigned i=0, j=0; i< q.island_ijoints.size(); i++)
  {
    const unsigned NEQ = q.island_ijoints[i]->num_constraint_eqns();
    q.island_ijoints[i]->lambda = lambda.segment(j, j+NEQ);
    j += NEQ;
  }
}

/// Gets the minimum constraint velocity
double ImpactConstraintHandler::calc_min_constraint_velocity(const UnilateralConstraintProblemData& q) const
{
  double minv = std::numeric_limits<double>::max();

  // see whether another QP must be solved
  if (q.Cn_v.size() > 0)
    minv = *min_element(q.Cn_v.column_iterator_begin(), q.Cn_v.column_iterator_end());
  if (q.L_v.size() > 0)
    minv = std::min(minv, *min_element(q.L_v.column_iterator_begin(), q.L_v.column_iterator_end()));

  return minv;
}

/// Updates post-impact velocities
void ImpactConstraintHandler::update_constraint_velocities_from_impulses(UnilateralConstraintProblemData& q)
{
  // update Cn_v
  q.Cn_v += q.Cn_X_CnT.mult(q.cn, _a);
  q.Cn_v += q.Cn_X_CsT.mult(q.cs, _a);
  q.Cn_v += q.Cn_X_CtT.mult(q.ct, _a);
  q.Cn_v += q.Cn_X_LT.mult(q.l, _a);

  // update Cs_v
  q.Cs_v += q.Cn_X_CsT.transpose_mult(q.cn, _a);
  q.Cs_v += q.Cs_X_CsT.mult(q.cs, _a);
  q.Cs_v += q.Cs_X_CtT.mult(q.ct, _a);
  q.Cs_v += q.Cs_X_LT.mult(q.l, _a);

  // update Ct_v
  q.Ct_v += q.Cn_X_CtT.transpose_mult(q.cn, _a);
  q.Ct_v += q.Cs_X_CtT.transpose_mult(q.cs, _a);
  q.Ct_v += q.Ct_X_CtT.mult(q.ct, _a);
  q.Ct_v += q.Ct_X_LT.mult(q.l, _a);

  // update L_v
  q.L_v += q.Cn_X_LT.transpose_mult(q.cn, _a);
  q.L_v += q.Cs_X_LT.transpose_mult(q.cs, _a);
  q.L_v += q.Ct_X_LT.transpose_mult(q.ct, _a);
  q.L_v += q.L_X_LT.mult(q.l, _a);

  // output results
  FILE_LOG(LOG_CONSTRAINT) << "results: " << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cn: " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs: " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct: " << q.ct << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "l: " << q.l << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "new Cn_v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "new Cs_v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "new Ct_v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "new L_v: " << q.L_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "new Jx_v: " << q.Jx_v << std::endl;
}

/// Applies restitution to impact problem
/**
 * \return false if no restitution was applied; true otherwise
 */
bool ImpactConstraintHandler::apply_restitution(const UnilateralConstraintProblemData& q, VectorNd& z) const
{
  bool changed = false;

  // apply (Poisson) restitution to contacts
  for (unsigned i=0, j=q.CN_IDX; i< q.N_CONTACTS; i++, j++)
  {
    z[j] *= q.contact_constraints[i]->contact_epsilon;
    if (!changed && z[j] > NEAR_ZERO)
      changed = true;
  }

  // apply (Poisson) restitution to limits
  for (unsigned i=0, j=q.L_IDX; i< q.N_LIMITS; i++, j++)
  {
    z[j] *= q.limit_constraints[i]->limit_epsilon;
    if (!changed && z[j] > NEAR_ZERO)
      changed = true;
  }

  return changed;
}

/// Applies restitution to impact problem
/**
 * \return false if no restitution was applied; true otherwise
 */
bool ImpactConstraintHandler::apply_restitution(UnilateralConstraintProblemData& q) const
{
  bool changed = false;

  // apply (Poisson) restitution to contacts
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    q.cn[i] *= q.contact_constraints[i]->contact_epsilon;
    if (!changed && q.cn[i] > NEAR_ZERO)
      changed = true;
  }

  // apply (Poisson) restitution to limits
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    q.l[i] *= q.limit_constraints[i]->limit_epsilon;
    if (!changed && q.l[i] > NEAR_ZERO)
      changed = true;
  }

  if (changed)
  {
    q.cs.set_zero();
    q.ct.set_zero();
  }

  return changed;
}

/**
 * Applies method of Drumwright and Shell to a set of connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints)
{
  double ke_minus = 0.0, ke_plus = 0.0;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // set the simulator
  _epd.simulator = _simulator;

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd);

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->body_id << " pre-constraint handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // use QP / NQP solver with warm starting to find the solution
  if (use_qp_solver(_epd))
    solve_qp(_z, _epd);
  else
    solve_nqp(_z, _epd);

  // update the impulses from z
  update_from_stacked(_epd, _z);

  // determine velocities due to impulse application
  update_constraint_velocities_from_impulses(_epd);

  // get the constraint violation before applying impulses
  double minv = calc_min_constraint_velocity(_epd);

  // apply restitution
  if (apply_restitution(_epd, _z))
  {
    // update the impulses from z
    update_from_stacked(_epd, _z);

    // determine velocities due to impulse application
    update_constraint_velocities_from_impulses(_epd);

    // check to see whether we need to solve another impact problem
    double minv_plus = calc_min_constraint_velocity(_epd);
    FILE_LOG(LOG_CONSTRAINT) << "Applying restitution" << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  compression v+ minimum: " << minv << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  restitution v+ minimum: " << minv_plus << std::endl;
    if (minv_plus < 0.0 && minv_plus < minv - NEAR_ZERO)
    {
      // use QP / NQP solver with warm starting to find the solution
      if (use_qp_solver(_epd))
        solve_qp(_z, _epd);
      else
        solve_nqp(_z, _epd);

      // update the impulses from z
      update_from_stacked(_epd, _z);
    }
  }

  // for debugging
/*
  if (LOGGING(LOG_CONSTRAINT))
  { 
    FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler debugging check (slow) " << std::endl;
    compute_problem_data(_epd);
    FILE_LOG(LOG_CONSTRAINT) << "new Cn_v (double check): " << _epd.Cn_v << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "new Cs_v (double check): " << _epd.Cs_v << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "new Ct_v (double check): " << _epd.Ct_v << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "new L_v (double check): " << _epd.L_v << std::endl;

    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->body_id << " post-constraint handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONSTRAINT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }
*/
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() exiting" << endl;
}

/// Determines whether we can use the QP solver
bool ImpactConstraintHandler::use_qp_solver(const UnilateralConstraintProblemData& epd)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // first, check whether any contact constraints use a true friction cone
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
    if (epd.contact_constraints[i]->contact_NK == UINF)
      return false;

  // still here? ok to use QP solver
  return true;
}

/// Propagates impulse data from ConstraintProblemData to the underlying unilateral constraints
void ImpactConstraintHandler::propagate_impulse_data(const UnilateralConstraintProblemData& q)
{
  shared_ptr<Pose3d> P(new Pose3d);

  // save normal contact impulses
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = q.contact_constraints[i]->contact_normal * q.cn[i];
    j += q.contact_constraints[i]->contact_tan1 * q.cs[i];
    j += q.contact_constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += Pose3d::transform(GLOBAL, jx);
  }

  // save normal contact impulses
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    q.limit_constraints[i]->limit_impulse += q.l[i];
  }
}

/// Applies impulses to bodies
void ImpactConstraintHandler::apply_impulses(const UnilateralConstraintProblemData& q)
{
  map<shared_ptr<DynamicBodyd>, VectorNd> gj;
  map<shared_ptr<DynamicBodyd>, VectorNd>::iterator gj_iter;

  // loop over all contact constraints first
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& e = *q.contact_constraints[i];
    SForced w(e.contact_impulse);

    // get the two single bodies of the contact
    shared_ptr<SingleBodyd> sb1 = e.contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> b1 = dynamic_pointer_cast<DynamicBodyd>(sb1->get_super_body());
    shared_ptr<DynamicBodyd> b2 = dynamic_pointer_cast<DynamicBodyd>(sb2->get_super_body());

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, _v);
      gj_iter->second += _v;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, _v);
      gj_iter->second += _v;
    }
  }

  // loop over all limit constraints next
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *q.limit_constraints[i];
    ArticulatedBodyPtr ab = e.limit_joint->get_articulated_body();

    // get the iterator for the articulated body
    gj_iter = gj.find(ab);

    // cast as an RCArticulatedBody
    shared_ptr<RCArticulatedBody> rcab = dynamic_pointer_cast<RCArticulatedBody>(ab);

    // apply limit impulses to bodies in independent coordinates
    if (rcab)
    {
      // get the index of the joint
      unsigned idx = e.limit_joint->get_coord_index() + e.limit_dof;

      // initialize the vector if necessary
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(rcab->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(ab);
      }

      // set the limit force
      gj_iter->second[idx] += e.limit_impulse;
    }
  }

  // apply all generalized impacts
  for (map<shared_ptr<DynamicBodyd>, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/*
/// Computes the data to the LCP / QP problems
void ImpactConstraintHandler::compute_problem_data(UnilateralConstraintProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // determine set of "super" bodies from contact constraints
  q.super_bodies.clear();
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom1->get_single_body()));
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom2->get_single_body()));
  }

  // determine set of "super" bodies from limit constraints
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    RigidBodyPtr outboard = q.limit_constraints[i]->limit_joint->get_outboard_link();
    q.super_bodies.push_back(get_super_body(outboard));
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

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

  // compute number of friction polygon edges
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    if (q.contact_constraints[i]->contact_NK < UINF)
    {
      q.N_K_TOTAL += q.contact_constraints[i]->contact_NK/2;
      q.N_LIN_CONE++;
    }
    else if (q.contact_constraints[i]->contact_NK == UINF)
      break;
  }

  // setup number of true cones
  q.N_TRUE_CONE = q.contact_constraints.size() - q.N_LIN_CONE;

  // verify contact constraints that use a true friction cone are at the end
  // of the contact vector
  #ifndef NDEBUG
  for (unsigned i=q.N_LIN_CONE; i< q.contact_constraints.size(); i++)
    assert(q.contact_constraints[i]->contact_NK == UINF);
  #endif

  // initialize the problem matrices / vectors
  q.Cn_X_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cn_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Cs_X_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cs_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Ct_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Ct_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Ct_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.L_X_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.L_X_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  q.Jx_X_JxT.set_zero(q.N_CONSTRAINT_EQNS_IMP, q.N_CONSTRAINT_EQNS_IMP);
  q.Cn_v.set_zero(q.N_CONTACTS);
  q.Cs_v.set_zero(q.N_CONTACTS);
  q.Ct_v.set_zero(q.N_CONTACTS);
  q.L_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  q.l.set_zero(q.N_LIMITS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
  q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
  q.NCT_IDX = q.NCS_IDX + q.N_LIN_CONE;
  q.L_IDX = q.NCT_IDX + q.N_LIN_CONE;
  q.ALPHA_X_IDX = q.L_IDX + q.N_LIMITS;
  q.N_VARS = q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_IMP;

  // get iterators to the proper matrices
  RowIteratord CnCn = q.Cn_X_CnT.row_iterator_begin();
  RowIteratord CnCs = q.Cn_X_CsT.row_iterator_begin();
  RowIteratord CnCt = q.Cn_X_CtT.row_iterator_begin();
  RowIteratord CsCs = q.Cs_X_CsT.row_iterator_begin();
  RowIteratord CsCt = q.Cs_X_CtT.row_iterator_begin();
  RowIteratord CtCt = q.Ct_X_CtT.row_iterator_begin();

  // process contact constraints, setting up matrices
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // compute cross constraint data for contact constraints
    for (unsigned j=0; j< q.contact_constraints.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 3);

      // check whether i==j (single contact constraint)
      if (i == j)
      {
        // compute matrix / vector for contact constraint i
        _v.set_zero(3);
        q.contact_constraints[i]->compute_constraint_data(_MM, _v);

        // setup appropriate parts of contact inertia matrices
        RowIteratord_const data = _MM.row_iterator_begin();
        *CnCn = *data++;
        *CnCs = *data++;
        *CnCt = *data; data += 2; // advance past Cs_X_CnT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_X_CtT
        *CtCt = *data;

        // setup appropriate parts of contact velocities
        data = _v.row_iterator_begin();
        q.Cn_v[i] = *data++;
        q.Cs_v[i] = *data++;
        q.Ct_v[i] = *data;
      }
      else
      {
        // compute matrix for cross constraint
        q.contact_constraints[i]->compute_cross_constraint_data(*q.contact_constraints[j], _MM);

        // setup appropriate parts of contact inertia matrices
        RowIteratord_const data = _MM.row_iterator_begin();
        *CnCn = *data++;
        *CnCs = *data++;
        *CnCt = *data; data += 2; // advance to Cs_X_CsT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_X_CtT
        *CtCt = *data;
      }

      // advance the iterators
      CnCn++;
      CnCs++;
      CnCt++;
      CsCs++;
      CsCt++;
      CtCt++;
    }

    // compute cross constraint data for contact/limit constraints
    for (unsigned j=0; j< q.limit_constraints.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 1);

      // compute matrix for cross constraint
      q.contact_constraints[i]->compute_cross_constraint_data(*q.limit_constraints[j], _MM);

      // setup appropriate parts of contact / limit inertia matrices
      ColumnIteratord_const data = _MM.column_iterator_begin();
      q.Cn_X_LT(i,j) = *data++;
      q.Cs_X_LT(i,j) = *data++;
      q.Ct_X_LT(i,j) = *data;
    }
  }

  // process limit constraints, setting up matrices
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    // compute matrix / vector for contact constraint i
    q.limit_constraints[i]->compute_constraint_data(_MM, _v);

    // setup appropriate entry of limit inertia matrix and limit velocity
    q.L_X_LT(i,i) = _MM.data()[0];
    q.L_v[i] = _v.data()[0];

    // compute cross/cross limit constraint data
    for (unsigned j=i+1; j< q.limit_constraints.size(); j++)
    {
      // reset _MM
      _MM.resize(1,1);

      // compute matrix for cross constraint
      q.limit_constraints[i]->compute_cross_constraint_data(*q.limit_constraints[j], _MM);

      // setup appropriate part of limit / limit inertia matrix
      q.L_X_LT(i,j) = q.L_iM_LT(j,i) = _MM.data()[0];
    }

    // NOTE: cross data has already been computed for contact/limit constraints
  }
}
*/

/// Solves the viscous friction LCP
void ImpactConstraintHandler::apply_visc_friction_model(UnilateralConstraintProblemData& q)
{
  // compute the (Coulomb) frictionless LCP
  VectorNd z;
  solve_frictionless_lcp(q, z);

  // setup impulses
  q.cn = z.segment(q.CN_IDX, q.N_CONTACTS);
  q.l = z.segment(q.L_IDX, q.L_IDX+q.N_LIMITS);
  q.cs = _cs_visc;
  q.ct = _ct_visc;
  q.cs.negate();
  q.ct.negate();

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save contact impulses
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = q.contact_constraints[i]->contact_normal * q.cn[i];
    j += q.contact_constraints[i]->contact_tan1 * q.cs[i];
    j += q.contact_constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += Pose3d::transform(GLOBAL, jx);
  }

  // save limit impulses
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    double limit_impulse = (q.limit_constraints[i]->limit_upper) ? -q.l[i] : q.l[i];
    q.limit_constraints[i]->limit_impulse += limit_impulse;
  }

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model() exited" << std::endl;
}

/// Solves the no-slip model LCP
void ImpactConstraintHandler::apply_no_slip_model(UnilateralConstraintProblemData& q)
{
  std::vector<unsigned> S_indices, T_indices;
  const unsigned NCONTACTS = q.N_CONTACTS;
  const unsigned NLIMITS = q.N_LIMITS;
  const unsigned N_IDX = 0;
  const unsigned L_IDX = N_IDX + NCONTACTS;
  VectorNd lb, ub, b;
  MatrixNd A;
  double ke_plus = 0.0, ke_minus = 0.0;

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->body_id << " pre-constraint handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_X_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cs': " << std::endl << q.Cn_X_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Ct': " << std::endl << q.Cn_X_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Cs': " << std::endl << q.Cs_X_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Ct': " << std::endl << q.Cs_X_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * L': " << std::endl << q.Cs_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Ct': " << std::endl << q.Ct_X_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * L': " << std::endl << q.Ct_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * L': " << std::endl << q.L_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * v: " << q.Cs_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * v: " << q.Ct_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << q.L_v << std::endl;

  // we do this by solving the MLCP:
  // |  A  C  | | u | + | a | = | 0 |
  // |  D  B  | | v |   | b |   | r |

  // a = [-M*v'; 0]
  // b = 0
  // B = 0
  // C = [ -N' -L' ]
  // D = -C'
  // u = [ v^+; cs; ct; alphax ]
  // v = [ cn; l ]
  // r = [ Cn*v+; L*v+ ]

  // Assuming that C is of full row rank (no dependent joint constraints)
  // A is invertible; then we just need to solve the LCP:

  // | B - D*inv(A)*C | | v | + | b - D*inv(A)*a | = | w |
  // and use the result to solve for u:
  // u = -inv(A)*(a + Cv)

  // A is the matrix | M X'|
  //                 | X 0 |  where X is [ S; T; J ]
  // blockwise inversion yields inv(A) =
  // inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y
  // Y*X*inv(M)                    -Y
  // where Y = inv(X*inv(M)*X')

  // defining Q = [Cn; L; 0] and using the result above yields following LCP:
  // matrix: Q*inv(A)*Q' = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  // vector: -Q*inv(A)*a  = -Q*v + Q*inv(M)*X'*Y*X*v

  // NOTE: we can check whether joint constraints are linearly dependent
  // (linearly dependent constraints should be discarded) if J*inv(M)*J'
  // is singular (using a Cholesky factorization). We can use the same
  // principle to determine whether a contact direction should be discarded

  // ********************************************************
  // find largest non-singular set of S and T indices
  // ********************************************************

  // loop through contacts, forming matrix below and checking its condition
  // | S*inv(M)*S'  S*inv(M)*T' |
  // | T*inv(M)*S'  T*inv(M)*T' |
  for (unsigned i=0; i< NCONTACTS; i++)
  {
    // update S indices
    S_indices.push_back(i);

    // setup indices
    unsigned S_IDX = 0;
    unsigned T_IDX = S_indices.size();
    unsigned J_IDX = T_IDX + T_indices.size();
    _Y.resize(J_IDX, J_IDX);

    // add S/S, T/T, components to 'check' matrix
    q.Cs_X_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, S_IDX, _MM);
    q.Ct_X_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, T_IDX, _MM);

    // add S/T components to 'check' matrix
    q.Cs_X_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, T_IDX, _MM);
    _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< _Y.rows(); j++)
      _Y(j,j) -= NEAR_ZERO;

    // see whether check matrix can be Cholesky factorized
    if (!_LA.factor_chol(_Y))
      S_indices.pop_back();

    // add index for T
    T_indices.push_back(i);

    // resize the check matrix
    T_IDX = S_indices.size();
    J_IDX = T_IDX + T_indices.size();
    _Y.resize(J_IDX, J_IDX);

    // add S/S, T/T components to 'check' matrix
    q.Cs_X_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, S_IDX, _MM);
    q.Ct_X_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, T_IDX, _MM);

    // add S/T components to 'check' matrix
    q.Cs_X_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, T_IDX, _MM);
    _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< _Y.rows(); j++)
      _Y(j,j) -= NEAR_ZERO;

    // see whether check matrix can be Cholesky factorized
    if (!_LA.factor_chol(_Y))
      T_indices.pop_back();
  }

  // output indices
  if (LOGGING(LOG_CONSTRAINT))
  {
    std::ostringstream oss;
    oss << "s indices:";
    for (unsigned i=0; i< S_indices.size(); i++)
      oss << " " << S_indices[i];
    oss << "  t indices:";
    for (unsigned i=0; i< T_indices.size(); i++)
      oss << " " << T_indices[i];
    FILE_LOG(LOG_CONSTRAINT) << oss.str() << std::endl;
  }

  // ********************************************************
  // reform Y (last matrix was de-regularized) 
  // ********************************************************

  // setup indices
  const unsigned S_IDX = 0;
  const unsigned T_IDX = S_indices.size();
  const unsigned J_IDX = T_IDX + T_indices.size();
  _Y.resize(J_IDX, J_IDX);

  // add S/S, T/T components to X
  q.Cs_X_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
  _Y.set_sub_mat(S_IDX, S_IDX, _MM);
  q.Ct_X_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
  _Y.set_sub_mat(T_IDX, T_IDX, _MM);

  // add S/T components to X
  q.Cs_X_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
  _Y.set_sub_mat(S_IDX, T_IDX, _MM);
  _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

  // do the Cholesky factorization (should not fail)
  bool success = _LA.factor_chol(_Y);
  assert(success);

  // defining Y = inv(X*inv(M)*X') and Q = [Cn; L; 0]
  // and using the result above yields following LCP:
  // matrix: Q*inv(A)*Q' = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  // vector: Q*inv(A)*a  = Q*v - Q*inv(M)*X'*Y*X*v

  // setup Q*inv(M)*Q'
  _MM.set_zero(q.N_CONTACTS + q.N_LIMITS, q.N_CONTACTS + q.N_LIMITS);
  _MM.set_sub_mat(N_IDX, N_IDX, q.Cn_X_CnT);
  _MM.set_sub_mat(N_IDX, L_IDX, q.Cn_X_LT);
  _MM.set_sub_mat(L_IDX, N_IDX, q.Cn_X_LT, Ravelin::eTranspose);
  _MM.set_sub_mat(L_IDX, L_IDX, q.L_X_LT);

  // setup Q*inv(M)*X'
  _Q_X_XT.resize(q.N_CONTACTS + q.N_LIMITS, S_indices.size() + T_indices.size());
  q.Cn_X_CsT.select_columns(S_indices.begin(), S_indices.end(), _workM);
  _Q_X_XT.set_sub_mat(N_IDX, S_IDX, _workM);
  q.Cn_X_CtT.select_columns(T_indices.begin(), T_indices.end(), _workM);
  _Q_X_XT.set_sub_mat(N_IDX, T_IDX, _workM);
  q.Cs_X_LT.select_rows(S_indices.begin(), S_indices.end(), _workM);
  _Q_X_XT.set_sub_mat(L_IDX, S_IDX, _workM, Ravelin::eTranspose);
  q.Ct_X_LT.select_rows(T_indices.begin(), T_indices.end(), _workM);
  _Q_X_XT.set_sub_mat(L_IDX, T_IDX, _workM, Ravelin::eTranspose);

  // compute Y*X*inv(M)*Q'
  MatrixNd::transpose(_Q_X_XT, _workM);
  _LA.solve_chol_fast(_Y, _workM);

  // compute Q*inv(M)*X'*Y*X*inv(M)*Q'
  _Q_X_XT.mult(_workM, _workM2);
  _MM -= _workM2;

  // setup -Q*v
  _qq.resize(q.N_CONTACTS + q.N_LIMITS);
  _qq.set_sub_vec(N_IDX, q.Cn_v);
  _qq.set_sub_vec(L_IDX, q.L_v);

  // setup X*v
  _Xv.resize(S_indices.size() + T_indices.size());
  q.Cs_v.select(S_indices.begin(), S_indices.end(), _workv);
  _Xv.set_sub_vec(S_IDX, _workv);
  q.Ct_v.select(T_indices.begin(), T_indices.end(), _workv);
  _Xv.set_sub_vec(T_IDX, _workv);

  // compute Y*X*v
  _YXv = _Xv;
  _LA.solve_chol_fast(_Y, _YXv);

  // compute Q*inv(M)*X' * Y*X*v
  _Q_X_XT.mult(_YXv, _workv);

  // setup remainder of LCP vector
  _qq -= _workv;

  // attempt to solve the LCP using the fast method
  if (!_lcp.lcp_fast(_MM, _qq, _v))
  {
    FILE_LOG(LOG_CONSTRAINT) << "Principal pivoting method LCP solver failed; falling back to slower solvers" << std::endl;

    #ifdef USE_QLCPD
    // solve didn't work; attempt to solve using QP solver
    (_workv = _qq) *= 0.5;
    lb.set_zero(_qq.size());
    ub.set_one(_qq.size()) *= 1e+29;
    A.set_zero(0, _qq.size());
    b.resize(0);
    (_workv2 = _qq).negate();
    if (!_qp.qp_activeset(_MM, _workv, lb, ub, _MM, _workv2, A, b, _v))
    {
      FILE_LOG(LOG_CONSTRAINT) << "QLCPD failed to find feasible point; finding closest feasible point" << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << "  old LCP q: " << _qq << std::endl;

      // QP solver didn't work; solve LP to find closest feasible solution
      if (!_qp.find_closest_feasible(lb, ub, _MM, _workv2, A, b, _v))
        throw std::runtime_error("Unable to solve constraint LCP!");

      // modify constraints
      _MM.mult(_v, _workv2) += _qq;
      for (unsigned i=0; i< _qq.size(); i++)
        if (_workv2[i] < 0.0)
          _qq[i] += (_workv2[i] - NEAR_ZERO);
      FILE_LOG(LOG_CONSTRAINT) << "  new LCP q: " << _qq << std::endl;

      // now try solving again
      (_workv2 = _qq).negate();
      if (!_qp.qp_activeset(_MM, _workv, lb, ub, _MM, _workv2, A, b, _v))
      {
        FILE_LOG(LOG_CONSTRAINT) << "QLCPD failed to find feasible point *twice*" << std::endl;
        throw std::runtime_error("Unable to solve constraint LCP!");
      }
    }
    #else
    if (!_lcp.lcp_lemke_regularized(_MM, _qq, _v))
      throw std::runtime_error("Unable to solve constraint LCP!");
    #endif
  }

  // compute the joint constraint forces and friction forces
  // u = -inv(A)*(a + Cv)
  // u = inv(A)*(Q'*[cn; l; 0])  [b/c we don't care about new velocity]
  // recalling that inv(A) =
  // | inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y | ngc x ngc,    ngc x sz(x)
  // | Y*X*inv(M)                    -Y          | sz(x) x ngc,  sz(x) x sz(x)
  // Q is nlcp x (ngc + sz(x))
  // [cs; ct; alphax] = -Y*X*v - Y*X*inv(M)*Q'*[cn; ct]
  _cs_ct_alphax = _YXv;
  _Q_X_XT.transpose_mult(_v, _workv);
  _LA.solve_chol_fast(_Y, _workv);
  _cs_ct_alphax += _workv;
  _cs_ct_alphax.negate();

  // setup impulses
  q.cn = _v.segment(0, q.N_CONTACTS);
  q.l = _v.segment(q.N_CONTACTS, _v.size());
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  SharedConstVectorNd cs_vec = _cs_ct_alphax.segment(S_IDX, T_IDX);
  SharedConstVectorNd ct_vec = _cs_ct_alphax.segment(T_IDX, J_IDX);
  q.cs.set(S_indices.begin(), S_indices.end(), cs_vec);
  q.ct.set(T_indices.begin(), T_indices.end(), ct_vec);

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save contact impulses
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = q.contact_constraints[i]->contact_normal * q.cn[i];
    j += q.contact_constraints[i]->contact_tan1 * q.cs[i];
    j += q.contact_constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += Pose3d::transform(GLOBAL, jx);
  }

  // save limit impulses
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    double limit_impulse = (q.limit_constraints[i]->limit_upper) ? -q.l[i] : q.l[i];
    q.limit_constraints[i]->limit_impulse += limit_impulse;
  }

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->body_id << " post-constraint handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONSTRAINT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }


  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_no_slip_lcp() exited" << std::endl;
}

/// Solves the (frictionless) LCP
void ImpactConstraintHandler::solve_frictionless_lcp(UnilateralConstraintProblemData& q, VectorNd& z)
{
  const unsigned NCONTACTS = q.N_CONTACTS;
  const unsigned NLIMITS = q.N_LIMITS;
  const unsigned NIMP = q.N_CONSTRAINT_EQNS_IMP;

  // we do this by solving the MLCP:
  // |  A  C  | | u | + | a | = | 0 |
  // |  D  B  | | v |   | b |   | r |

  // A is the matrix Jx*inv(M)*Jx', Jx is implicit joint constraint Jacobians
  // NOTE: we assume that Jx is of full row rank (no dependent constraints)

  // u = alphax
  // v = [ cn; l ]
  // r = [ Cn*v+; L*v+ ]
  // a = v - inv(M)*S'*muv*S*v - inv(M)*T'*muv*T*v
  // b = 0

  // Assuming that C is of full row rank (no dependent joint constraints)
  // A is invertible; then we just need to solve the LCP:

  // | B - D*inv(A)*C | | v | + | b - D*inv(A)*a | = | w |
  // and use the result to solve for u:
  // u = -inv(A)*(a + Cv)

  // compute SVD of Jx*inv(M)*Jx'
  _A = q.Jx_X_JxT;
  _LA.svd(_A, _AU, _AS, _AV);

  // setup the B matrix
  // B = [ Cn; L ]*inv(M)*[ Cn' L' ]
  _B.resize(NCONTACTS+NLIMITS, NCONTACTS+NLIMITS);
  _B.set_sub_mat(0, 0, q.Cn_X_CnT);
  _B.set_sub_mat(0, NCONTACTS, q.Cn_X_LT);
  _B.set_sub_mat(NCONTACTS, 0, q.Cn_X_LT, Ravelin::eTranspose);
  _B.set_sub_mat(NCONTACTS, NCONTACTS, q.L_X_LT);

  // setup the C matrix and compute inv(A)*C
  // C = Jx*inv(M)*[ Cn' L' ]; note: D = C'
  _C.resize(NIMP, NCONTACTS+NLIMITS);
  _C.set_sub_mat(0,0, q.Cn_X_JxT, Ravelin::eTranspose);
  _C.set_sub_mat(0,NCONTACTS, q.L_X_JxT, Ravelin::eTranspose);
  MatrixNd::transpose(_C, _D);
  _LA.solve_LS_fast(_AU, _AS, _AV, _C);

  // setup the a vector and compute inv(A)*a
  // a = [ Jx*v ]
  _a = q.Jx_v;
  _LA.solve_LS_fast(_AU, _AS, _AV, _a);

  // setup the b vector
  // b = [ Cn*v; L*v ]
  _b.resize(NLIMITS+NCONTACTS);
  _b.set_sub_vec(0, q.Cn_v);
  _b.set_sub_vec(NCONTACTS, q.L_v);

  // compute viscous friction terms
  _cs_visc = q.Cs_v;
  _ct_visc = q.Ct_v;
  RowIteratord cs_visc_iter = _cs_visc.row_iterator_begin();
  RowIteratord ct_visc_iter = _ct_visc.row_iterator_begin();
  for (unsigned i=0; i< NCONTACTS; i++, cs_visc_iter++, ct_visc_iter++)
  {
    (*cs_visc_iter) *= q.contact_constraints[i]->contact_mu_viscous;
    (*ct_visc_iter) *= q.contact_constraints[i]->contact_mu_viscous;
  }

  // compute viscous friction terms contributions in normal directions
  SharedVectorNd bsub = _b.segment(0, NCONTACTS);
  q.Cn_X_CsT.mult(_cs_visc, _workv);
  bsub -= _workv;
  q.Cn_X_CtT.mult(_ct_visc, _workv);
  bsub -= _workv;

  // setup the LCP matrix
  _D.mult(_C, _MM);
  _MM -= _B;
  _MM.negate();

  // setup the LCP vector
  _D.mult(_a, _qq);
  _qq -= _b;
  _qq.negate();

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_lcp() entered" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_X_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP vector: " << _qq << std::endl;

  // solve the LCP
  if (!_lcp.lcp_fast(_MM, _qq, _v) && !_lcp.lcp_lemke_regularized(_MM, _qq, _v))
    throw std::runtime_error("Unable to solve constraint LCP!");

  // determine the value of kappa
  SharedConstVectorNd cn = _v.segment(0, q.N_CONTACTS);
  SharedConstVectorNd l = _v.segment(q.N_CONTACTS, _v.size());
  q.Cn_X_CnT.mult(cn, _Cn_vplus) += q.Cn_v;

  // setup the homogeneous solution
  z.set_zero(q.N_VARS);
  z.set_sub_vec(q.CN_IDX, cn);
  z.set_sub_vec(q.L_IDX, l);

  FILE_LOG(LOG_CONSTRAINT) << "  LCP result: " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_lcp() exited" << std::endl;
}

/// Multiplies a block diagonal matrix formed of inertia matrices (or inverse inertia matrix) by a matrix
MatrixNd& ImpactConstraintHandler::mult(const std::vector<MatrixNd>& inertias, const MatrixNd& X, MatrixNd& B)
{
  MatrixNd tmp; 

  // resize B and set to zero
  B.set_zero(X.rows(), X.columns());

  // carry out the multiplication, one block at a time
  for (unsigned i=0, gc=0; i< inertias.size(); i++)
  {
    // get the number of rows from the i'th inertia
    const unsigned R = inertias[i].rows();

    // get the appropriate blocks from X and B
    SharedConstMatrixNd X_block = X.block(gc, gc+R, 0, X.columns());
    SharedMatrixNd B_block = B.block(gc, gc+R, 0, B.columns());

    // do the multiplication
    inertias[i].mult(X_block, tmp);
    B_block += tmp;

    // update gc 
    gc += R;
  }
}

/// Converts a sparse form of block diagonal matrices into a dense matrix
MatrixNd& ImpactConstraintHandler::to_dense(const vector<MatrixNd>& J, MatrixNd& dense)
{
  // get the number of rows in all matrices
  unsigned N = 0;
  for (unsigned i=0; i< J.size(); i++)
  {
    assert(J[i].rows() == J[i].columns());
    N += J[i].rows();
  }

  // resize the dense matrix
  dense.set_zero(N, N);

  // set the blocks
  for (unsigned i=0, j=0; i< J.size(); i++)
  {
    dense.block(j, j+J[i].rows(), j, j+J[i].rows()) = J[i];
    j += J[i].rows();
  }

  return dense;
}

/// Function for counting
static bool IsTrue(bool x) { return x; }

/// Computes the simplification used for X = C'*M*C
/**
 * C is defined as the upper left hand block of 
 *  | M J' |^-1  
 *  | J 0  |
 */
void ImpactConstraintHandler::compute_X(UnilateralConstraintProblemData& q, MatrixNd& X)
{
  MatrixNd JiM, iJ_iM_JT, HT, HTJiM, G, MG, GTMG;
  vector<MatrixNd> inertias, inv_inertias;

  // count number of active constraints
  const unsigned N_ACTIVE = std::count_if(q.active.begin(), q.active.end(), IsTrue);

  // form inertias and inverse inertia matrices
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    inertias.push_back(MatrixNd());
    shared_ptr<RigidBodyd> rb = dynamic_pointer_cast<RigidBodyd>(q.super_bodies[i]);
    if (!rb || rb->is_enabled())
    {
      q.super_bodies[i]->get_generalized_inertia(inertias.back());
      inv_inertias.push_back(inertias.back());
      LinAlgd::inverse_SPD(inv_inertias.back());
    }
    else
      inv_inertias.push_back(MatrixNd());
  } 

  // get the full row rank version of J
  q.J.blocks.clear();
  q.J.rows = N_ACTIVE;
  q.J.cols = q.Jfull.cols;
  for (unsigned i=0, j=0; i< q.Jfull.blocks.size(); i++)
  {
    unsigned row_start = q.Jfull.blocks[i].st_row_idx;
    unsigned row_end = row_start; 
    while (row_end < q.Jfull.blocks[i].st_row_idx + q.Jfull.blocks[i].rows())
    {
      // see whether the index is active for row_end
      if (!q.active[row_end])
      {
        // it's not active, put everything in a block up to this point
        if (row_end - row_start > 0)
        {
          const unsigned N_SUBTRACTED = std::count(q.active.begin(), q.active.begin()+q.Jfull.blocks[i].st_row_idx, false);
          q.J.blocks.push_back(MatrixBlock());
          q.J.blocks.back().st_col_idx = q.Jfull.blocks[i].st_col_idx;
          q.J.blocks.back().st_row_idx = q.Jfull.blocks[i].st_row_idx - N_SUBTRACTED;
          q.J.blocks.back().block = q.Jfull.blocks[i].block.block(row_start - q.Jfull.blocks[i].st_row_idx, row_end - q.Jfull.blocks[i].st_row_idx, 0, q.Jfull.blocks[i].block.columns());

          // update row_start and row_end
          row_start = row_end + 1;
          row_end = row_start;
        }
        else
        {
          row_start = row_end + 1;
          row_end = row_start;
        }
      }
      else
        row_end++; 
    }

    // add the last block if necessary
    if (row_end - row_start > 0)
    {
       const unsigned N_SUBTRACTED = std::count(q.active.begin(), q.active.begin()+q.Jfull.blocks[i].st_row_idx, false);
      q.J.blocks.push_back(MatrixBlock());
      q.J.blocks.back().st_col_idx = q.Jfull.blocks[i].st_col_idx;
      q.J.blocks.back().st_row_idx = q.Jfull.blocks[i].st_row_idx - N_SUBTRACTED;
      q.J.blocks.back().block = q.Jfull.blocks[i].block.block(row_start - q.Jfull.blocks[i].st_row_idx, row_end - q.Jfull.blocks[i].st_row_idx, 0, q.Jfull.blocks[i].block.columns());
    }
  }

  // compute J*inv(M)*J'
  q.J.mult(inv_inertias, JiM);
  MatrixNd::transpose(JiM, q.iM_JxT);
  q.J.mult(q.iM_JxT, q.Jx_iM_JxT);

  // compute H' = (inv(J*inv(M)*J')*J)'
  iJ_iM_JT = q.Jx_iM_JxT;
  LinAlgd::inverse_SPD(iJ_iM_JT);
  q.J.transpose_mult(iJ_iM_JT, HT);

  // compute G
  HT.mult_transpose(q.iM_JxT, HTJiM);
  mult(inv_inertias, HTJiM, G);

  // form inverse inertia as a single matrix
  to_dense(inv_inertias, X);

  // compute G'*M*G
  mult(inertias, G, MG);
  G.transpose_mult(MG, GTMG);

  // scale G by 2
  G *= 2.0;

  // compute X
  X -= G;
  X += GTMG;
/*
  for (unsigned i=0, j=0; i< q.super_bodies.size(); i++)
  {
    const unsigned NGC = q.super_bodies->num_generalized_coordinates(DynamicBodyd::eSpatial);
    SharedConstMatrixNd Gshared = G.block(j, j+NGC, j, j+NGC);
    SharedConstMatrixNd GTMT
    j += NGC;
*/
} 

/// Gets the full rank set of implicit constraints
void ImpactConstraintHandler::get_full_rank_implicit_constraints(const SparseJacobian& J, vector<bool>& active)
{
  MatrixNd JJT, JJT_sub;
  const double DEREGULARIZATION_FACTOR = NEAR_ZERO;

  // resize the set of the active constraints
  active.resize(J.rows, false);
  if (active.empty())
    return;

  // compute J*J'
  J.mult_transpose(J, JJT);

  // keep trying to add constraints 
  for (unsigned i=0, n_active = 0; i< active.size(); i++)
  {
    // see whether the number of indices is maximized
    if (n_active == J.cols)
      break;
  
    // tentatively update the number of active indices
    n_active++;

    // try to add this indices
    active[i] = true;

    // get the submatrix
    JJT.select_square(active, JJT_sub);

    // deregularize
    RowIteratord row_iter = JJT_sub.row_iterator_begin();
    for (unsigned j=0; j< JJT_sub.rows(); j++, row_iter += (JJT_sub.rows()+1))
      *row_iter -= DEREGULARIZATION_FACTOR;

    // attempt to factorize it
    if (!LinAlgd::factor_chol(JJT_sub))
    {
      n_active--; 
      active[i] = false;
    }
  }
}

/// Gets the super body (articulated if any)
shared_ptr<DynamicBodyd> ImpactConstraintHandler::get_super_body(shared_ptr<SingleBodyd> sb)
{
  shared_ptr<ArticulatedBodyd> ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return dynamic_pointer_cast<DynamicBodyd>(sb);
}

/// Computes limit components for the UnilateralConstraintProblemData
/**
 * Computes X_LT matrix and computes L_v vector
 */
void ImpactConstraintHandler::compute_limit_components(const MatrixNd& X, UnilateralConstraintProblemData& q)
{
  // compute X_LT
  q.X_LT.resize(q.N_GC, q.N_LIMITS);
  for (unsigned i=0; i< q.N_LIMITS; i++)
    q.X_LT.column(i) = X.row(q.limit_indices[i]);

  // compute L_X_LT
  q.L_X_LT.resize(q.N_LIMITS, q.N_LIMITS);
  for (unsigned i=0; i< q.N_LIMITS; i++)
    for (unsigned j=i; j< q.N_LIMITS; j++)
    {
      q.L_X_LT(i,j) = X.row(q.limit_indices[i])[q.limit_indices[j]];
      q.L_X_LT(j,i) = q.L_X_LT(i,j);
    }    

  // compute L_v
  q.L_v.resize(q.N_LIMITS);
  for (unsigned i=0; i< q.N_LIMITS; i++)
  {
    const UnilateralConstraint& u = *q.limit_constraints[i];
    q.L_v[i] = u.limit_joint->qd[u.limit_dof];
    if (u.limit_upper)
      q.L_v[i] = q.L_v[i];
  }
}

/// Updates the generalized velocity for all bodies in an island
void ImpactConstraintHandler::update_generalized_velocities(const UnilateralConstraintProblemData& q, const VectorNd& dv)
{
  VectorNd tmp;

  // loop through all super bodies, getting the generalized velocity
  for (unsigned i=0, j=0; i< q.super_bodies.size(); i++)
  {
    const unsigned N_GC = q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
    SharedConstVectorNd dv_sub = dv.segment(j, j+N_GC);
    q.super_bodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, tmp);
    tmp += dv_sub;
    q.super_bodies[i]->set_generalized_velocity(DynamicBodyd::eSpatial, tmp);
    j += N_GC; 
  }
}

/// Computes the generalized velocity for all bodies in an island
void ImpactConstraintHandler::get_generalized_velocity(const UnilateralConstraintProblemData& q, VectorNd& v)
{
  // resize v
  v.resize(q.N_GC);

  // loop through all super bodies, getting the generalized velocity
  for (unsigned i=0, j=0; i< q.super_bodies.size(); i++)
  {
    const unsigned N_GC = q.super_bodies[i]->num_generalized_coordinates(DynamicBodyd::eSpatial);
    SharedVectorNd v_sub = v.segment(j, j+N_GC);
    q.super_bodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, v_sub);
    j += N_GC; 
  }
}

/// Adds a contact constraint to the contact Jacobians
void ImpactConstraintHandler::add_contact_to_Jacobian(const UnilateralConstraint& c, SparseJacobian& Cn, SparseJacobian& Cs, SparseJacobian& Ct, const std::map<shared_ptr<DynamicBodyd>, unsigned>& gc_map, unsigned contact_idx)
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

  // add a row to each Jacobian
  Cn.rows++;
  Cs.rows++;
  Ct.rows++;

  // do this six times, one for each body and each direction
  add_contact_dir_to_Jacobian(rb1, su1, Cn, c.contact_point, c.contact_normal, gc_map, contact_idx);
  add_contact_dir_to_Jacobian(rb2, su2, Cn, c.contact_point, -c.contact_normal, gc_map, contact_idx);
  add_contact_dir_to_Jacobian(rb1, su1, Cs, c.contact_point, c.contact_tan1, gc_map, contact_idx);
  add_contact_dir_to_Jacobian(rb2, su2, Cs, c.contact_point, -c.contact_tan1, gc_map, contact_idx);
  add_contact_dir_to_Jacobian(rb1, su1, Ct, c.contact_point, c.contact_tan2, gc_map, contact_idx);
  add_contact_dir_to_Jacobian(rb2, su2, Ct, c.contact_point, -c.contact_tan2, gc_map, contact_idx);
} 

void ImpactConstraintHandler::add_contact_dir_to_Jacobian(shared_ptr<RigidBodyd> rb, shared_ptr<ArticulatedBodyd> ab, SparseJacobian& C, const Vector3d& contact_point, const Vector3d& d, const std::map<shared_ptr<DynamicBodyd>, unsigned>& gc_map, unsigned contact_index)
{
  const unsigned N_SPATIAL = 6;
  MatrixNd contact_wrench(1, N_SPATIAL), Jm, tmp;

  // check whether the body is enabled
  if (!rb->is_enabled())
    return;

  // get the vector from the center of mass to the contact point 
  Vector3d x0(Pose3d::calc_relative_pose(rb->get_inertial_pose(), GLOBAL).x, GLOBAL);
  Vector3d r = contact_point - x0; 

  // compute r x d
  Vector3d rxd = Vector3d::cross(r, d);

  // setup the contact wrench
  contact_wrench.row(0).set_sub_vec(0, d);
  contact_wrench.row(0).set_sub_vec(3, rxd);

  // put the Jacobian into independent coordinates if necessary
  if (ab)
  {
    // get the reduced coordinate body
    shared_ptr<RCArticulatedBodyd> rcab = dynamic_pointer_cast<RCArticulatedBodyd>(ab);
    if (rcab)
    {
      // get the Jacobian and carry out the multiplication
      rcab->calc_jacobian(rb->get_mixed_pose(), rb, Jm);
      contact_wrench.mult(Jm, tmp);
      contact_wrench = tmp;
    }
  }

  // add the block to the Jacobian
  C.blocks.push_back(MatrixBlock());
  C.blocks.back().block = contact_wrench;
  C.blocks.back().st_row_idx = contact_index;
  if (ab)
  {
    assert(gc_map.find(ab) != gc_map.end());
    C.blocks.back().st_col_idx = gc_map.find(ab)->second;
  }
  else
  {
    assert(gc_map.find(rb) != gc_map.end());
    C.blocks.back().st_col_idx = gc_map.find(rb)->second;
  }
}

/// Computes the data to the LCP / QP problems
void ImpactConstraintHandler::compute_problem_data(UnilateralConstraintProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned N_SPATIAL = 6;
  VectorNd v;
  MatrixNd X, tmp, Jm;
 
  // determine set of "super" bodies from contact constraints
  q.super_bodies.clear();
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom1->get_single_body()));
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom2->get_single_body()));
  }

  // make vector super bodies unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // get the signed distances
  q.signed_distances = _simulator->get_pairwise_distances();

  // remove signed distances that do not correspond to a contact- we do not
  // want to consider those (bodies may be well separated)
  for (unsigned i=0; i< q.signed_distances.size(); i++)
  {
    // get the two single bodies involved in contact
    shared_ptr<SingleBodyd> s1 = q.signed_distances[i].a->get_single_body();
    shared_ptr<SingleBodyd> s2 = q.signed_distances[i].b->get_single_body();

    // get the two super bodies involved in contact
    shared_ptr<DynamicBodyd> sb1 = ImpactConstraintHandler::get_super_body(s1); 
    shared_ptr<DynamicBodyd> sb2 = ImpactConstraintHandler::get_super_body(s2);

    // see whether we must consider this signed distance
    if (std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), sb1) ||        std::binary_search(q.super_bodies.begin(), q.super_bodies.end(), sb2))
      continue;

    // remove the signed distance
    q.signed_distances[i] = q.signed_distances.back();
    q.signed_distances.pop_back();
  }

  // determine set of "super" bodies from limit constraints
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    RigidBodyPtr outboard = q.limit_constraints[i]->limit_joint->get_outboard_link();
    q.super_bodies.push_back(get_super_body(outboard));
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

  // compute number of friction polygon edges
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    if (q.contact_constraints[i]->contact_NK < UINF)
    {
      q.N_K_TOTAL += q.contact_constraints[i]->contact_NK/2;
      q.N_LIN_CONE++;
    }
    else if (q.contact_constraints[i]->contact_NK == UINF)
      break;
  }

  // setup number of true cones
  q.N_TRUE_CONE = q.contact_constraints.size() - q.N_LIN_CONE;

  // verify contact constraints that use a true friction cone are at the end
  // of the contact vector
  #ifndef NDEBUG
  for (unsigned i=q.N_LIN_CONE; i< q.contact_constraints.size(); i++)
    assert(q.contact_constraints[i]->contact_NK == UINF);
  #endif

  // initialize the problem matrices / vectors
  q.Cn_X_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cn_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Cs_X_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cs_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Ct_X_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Ct_X_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Ct_X_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.L_X_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.L_X_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  q.Jx_X_JxT.set_zero(q.N_CONSTRAINT_EQNS_IMP, q.N_CONSTRAINT_EQNS_IMP);
  q.Cn_v.set_zero(q.N_CONTACTS);
  q.Cs_v.set_zero(q.N_CONTACTS);
  q.Ct_v.set_zero(q.N_CONTACTS);
  q.L_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  q.l.set_zero(q.N_LIMITS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
  q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
  q.NCT_IDX = q.NCS_IDX + q.N_LIN_CONE;
  q.L_IDX = q.NCT_IDX + q.N_LIN_CONE;
  q.N_VARS = q.L_IDX + q.N_LIMITS;

  // get generalized velocity
  get_generalized_velocity(q, v);

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
      MatrixNd& block = q.Jfull.blocks.back().block;
      q.island_ijoints[i]->calc_constraint_jacobian(true, block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[inboard];
    }
 
    if (outboard->is_enabled())
    {
      // add the block to the Jacobian
      q.Jfull.blocks.push_back(MatrixBlock());
      MatrixNd& block = q.Jfull.blocks.back().block;
      q.island_ijoints[i]->calc_constraint_jacobian(false, block);
      q.Jfull.blocks.back().st_row_idx = eq_idx;
      q.Jfull.blocks.back().st_col_idx = gc_map[outboard];
    }

    // update the equation index
    eq_idx += q.island_ijoints[i]->num_constraint_eqns();
  } 

  // determine active set of implicit constraints
  get_full_rank_implicit_constraints(q.Jfull, q.active);

  // compute X
  compute_X(q, X);

  // setup Jacobians for Cn, Cs, Ct
  SparseJacobian Cn, Cs, Ct;

  // setup the number of columns in each Jacobian
  Cn.cols = q.N_GC;
  Cs.cols = q.N_GC;
  Ct.cols = q.N_GC;

  // process all contact constraints
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
    add_contact_to_Jacobian(*q.contact_constraints[i], Cn, Cs, Ct, gc_map, i); 

  // compute X_CnT, X_CsT, and X_CtT
  Cn.mult(X, tmp);  MatrixNd::transpose(tmp, q.X_CnT);
  Cs.mult(X, tmp);  MatrixNd::transpose(tmp, q.X_CsT);
  Ct.mult(X, tmp);  MatrixNd::transpose(tmp, q.X_CtT);
  q.J.mult(X, tmp); MatrixNd::transpose(tmp, q.X_JxT);
  
  // compute limit components - must do this first
  compute_limit_components(X, q);

  // compute problem data for Cn rows
  Cn.mult(q.X_CnT, q.Cn_X_CnT); 
  Cn.mult(q.X_CsT, q.Cn_X_CsT);  
  Cn.mult(q.X_CtT, q.Cn_X_CtT);  
  Cn.mult(q.X_LT,  q.Cn_X_LT);  
  Cn.mult(q.X_JxT,  q.Cn_X_JxT);

  // compute problem data for Cs rows
  Cs.mult(q.X_CsT, q.Cs_X_CsT);  
  Cs.mult(q.X_CtT, q.Cs_X_CtT);  
  Cs.mult(q.X_LT,  q.Cs_X_LT);  
  Cs.mult(q.X_JxT,  q.Cs_X_JxT);  

  // compute problem data for Ct rows
  Ct.mult(q.X_CtT, q.Ct_X_CtT);  
  Ct.mult(q.X_LT,  q.Ct_X_LT);  
  Ct.mult(q.X_JxT,  q.Ct_X_JxT);  

  // compute problem data for limit rows
  q.L_X_JxT.resize(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  for (unsigned i=0; i< q.N_LIMITS; i++)
    q.L_X_JxT.row(i) = q.X_JxT.row(q.limit_indices[i]);

  // compute vectors
  Cn.mult(v, q.Cn_v);
  Cs.mult(v, q.Cs_v);
  Ct.mult(v, q.Ct_v);
  q.J.mult(v, q.Jx_v);

  // compute signed distance Jacobians
  #ifdef USE_SIGNED_DIST_CONSTRAINT
  SignedDistDot::compute_signed_dist_dot_Jacobians(q, q.Cdot_iM_CnT, q.Cdot_iM_CsT, q.Cdot_iM_CtT, q.Cdot_iM_LT, q.Cdot_v);
  #endif
}

