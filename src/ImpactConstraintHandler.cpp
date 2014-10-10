/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/permute.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactConstraintHandler.h>
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

  // setup variables to help warmstarting
  _last_contacts = _last_limits = _last_contact_constraints = 0;

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
 * \param max_time the maximum time to solve the constraints
 * \param inv_dt 1/dt (the time step) for correcting interpenetration
 */
void ImpactConstraintHandler::process_constraints(const vector<UnilateralConstraint>& constraints, double max_time, double inv_dt)
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // apply the method to all contacts
  if (!constraints.empty())
    apply_model(constraints, max_time, inv_dt);
  else
    FILE_LOG(LOG_CONSTRAINT) << " (no constraints?!)" << endl;
    
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of constraints 
/**
 * \param constraints a set of constraints
 * \param max_time the maximum time to solve the constraints
 * \param inv_dt 1/dt (the time step) for correcting interpenetration
 */
void ImpactConstraintHandler::apply_model(const vector<UnilateralConstraint>& constraints, double max_time, double inv_dt)
{
  const double INF = std::numeric_limits<double>::max();
  list<UnilateralConstraint*> impacting;

  // **********************************************************
  // determine sets of connected constraints 
  // **********************************************************
  list<list<UnilateralConstraint*> > groups;
  UnilateralConstraint::determine_connected_constraints(constraints, groups);
  UnilateralConstraint::remove_inactive_groups(groups);

  // **********************************************************
  // do method for each connected set 
  // **********************************************************
  for (list<list<UnilateralConstraint*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
    // determine contact tangents
    for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->constraint_type == UnilateralConstraint::eContact)
        (*j)->determine_contact_tangents();

      // copy the list of constraints
      list<UnilateralConstraint*> rconstraints = *i;

      FILE_LOG(LOG_CONSTRAINT) << " -- pre-constraint velocity (all constraints): " << std::endl;
      for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;

      // determine a reduced set of constraints
      UnilateralConstraint::determine_minimal_set(rconstraints);

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
        apply_no_slip_model_to_connected_constraints(rconstraints, inv_dt);
// TODO: fix viscous model- seems to be a bug in it
//      else if (all_frictionless)
//        apply_visc_friction_model_to_connected_constraints(rconstraints, inv_dt);

      if (max_time < INF)   
        apply_model_to_connected_constraints(rconstraints, max_time, inv_dt);
      else
        apply_model_to_connected_constraints(rconstraints, inv_dt);

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
}

/**
 * Applies method of Drumwright and Shell to a set of connected constraints.
 * "Anytime" version
 * \param constraints a set of connected constraints 
 */
void ImpactConstraintHandler::apply_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints, double max_time, double inv_dt)
{
  double ke_minus = 0.0, ke_plus = 0.0;
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd, inv_dt);

  // clear all impulses 
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->id << " pre-constraint handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // solve the no-slip linear complementarity problem to determine
  // the kappa constant
  VectorNd z;
  solve_frictionless_lcp(_epd, z);

  // update constraint problem data and z
  permute_problem(_epd, z);

  // determine N_ACT_K
  _epd.N_ACT_K = 0;
  for (unsigned i=0; i< _epd.N_ACT_CONTACTS; i++)
    if (_epd.contact_constraints[i]->contact_NK < UINF)
      _epd.N_ACT_K += _epd.contact_constraints[i]->contact_NK/2;

  // use QP / NQP solver with warm starting to find the solution
  if (use_qp_solver(_epd))
    solve_qp(z, _epd, max_time);
  else
    solve_nqp(z, _epd, max_time);

  // update the impulses from z 
  update_from_stacked(_epd, z);

  // determine velocities due to impulse application
  update_constraint_velocities_from_impulses(_epd);

  // get the constraint violation before applying impulses
  double minv = calc_min_constraint_velocity(_epd);

  // apply restitution
  if (apply_restitution(_epd, z))
  {
    // update the impulses from z 
    update_from_stacked(_epd, z);

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
      solve_frictionless_lcp(_epd, z);

      // update constraint problem data and z
      permute_problem(_epd, z);

      // determine N_ACT_K
      _epd.N_ACT_K = 0;
      for (unsigned i=0; i< _epd.N_ACT_CONTACTS; i++)
        if (_epd.contact_constraints[i]->contact_NK < UINF)
          _epd.N_ACT_K += _epd.contact_constraints[i]->contact_NK/2;

      // use QP / NQP solver with warm starting to find the solution
      if (use_qp_solver(_epd))
        solve_qp(z, _epd, max_time);
      else
        solve_nqp(z, _epd, max_time);

      // update the impulses from z 
      update_from_stacked(_epd, z);
    }  
  }

  // apply impulses 
  apply_impulses(_epd);

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->id << " post-constraint handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONSTRAINT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() exiting" << endl;
}

/**
 * Applies purely viscous friction model to connected constraints 
 * \param constraints a set of connected constraints 
 */
void ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints, double inv_dt)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd, inv_dt);

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

  // apply impulses 
  apply_impulses(_epd);

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model_to_connected_constraints() exiting" << endl;
}

/**
 * Applies no slip friction model to connected constraints 
 * \param constraints a set of connected constraints 
 */
void ImpactConstraintHandler::apply_no_slip_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints, double inv_dt)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_no_slip_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd, inv_dt);

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

  // apply impulses 
  apply_impulses(_epd);

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
  q.Cn_v += q.Cn_iM_CnT.mult(q.cn, _a);
  q.Cn_v += q.Cn_iM_CsT.mult(q.cs, _a);
  q.Cn_v += q.Cn_iM_CtT.mult(q.ct, _a);
  q.Cn_v += q.Cn_iM_LT.mult(q.l, _a);
  q.Cn_v += q.Cn_iM_JxT.mult(q.alpha_x, _a);

  // update Cs_v
  q.Cs_v += q.Cn_iM_CsT.transpose_mult(q.cn, _a);
  q.Cs_v += q.Cs_iM_CsT.mult(q.cs, _a);
  q.Cs_v += q.Cs_iM_CtT.mult(q.ct, _a);
  q.Cs_v += q.Cs_iM_LT.mult(q.l, _a);
  q.Cs_v += q.Cs_iM_JxT.mult(q.alpha_x, _a);

  // update Ct_v
  q.Ct_v += q.Cn_iM_CtT.transpose_mult(q.cn, _a);
  q.Ct_v += q.Cs_iM_CtT.transpose_mult(q.cs, _a);
  q.Ct_v += q.Ct_iM_CtT.mult(q.ct, _a);
  q.Ct_v += q.Ct_iM_LT.mult(q.l, _a);
  q.Ct_v += q.Ct_iM_JxT.mult(q.alpha_x, _a);

  // update L_v
  q.L_v += q.Cn_iM_LT.transpose_mult(q.cn, _a);
  q.L_v += q.Cs_iM_LT.transpose_mult(q.cs, _a);
  q.L_v += q.Ct_iM_LT.transpose_mult(q.ct, _a);
  q.L_v += q.L_iM_LT.mult(q.l, _a);
  q.L_v += q.L_iM_JxT.mult(q.alpha_x, _a);

  // update Jx_v
  q.Jx_v += q.Cn_iM_JxT.transpose_mult(q.cn, _a);
  q.Jx_v += q.Cs_iM_JxT.transpose_mult(q.cs, _a);
  q.Jx_v += q.Ct_iM_JxT.transpose_mult(q.ct, _a);
  q.Jx_v += q.L_iM_JxT.transpose_mult(q.l, _a);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, _a);

  // output results
  FILE_LOG(LOG_CONSTRAINT) << "results: " << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cn: " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs: " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct: " << q.ct << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "l: " << q.l << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "alpha_x: " << q.alpha_x << std::endl;
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
  for (unsigned i=0, j=q.CN_IDX; i< q.N_ACT_CONTACTS; i++, j++)
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

  return changed;
}

/// Permutes the problem to reflect active contact constraints
void ImpactConstraintHandler::permute_problem(UnilateralConstraintProblemData& epd, VectorNd& z)
{
  // determine mapping of old contact constraint indices to new contact constraint
  // indices
  std::vector<unsigned> mapping(epd.N_CONTACTS);

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::permute_problem() entered" << std::endl;

  // 1. compute active indices
  epd.N_ACT_CONTACTS = 0;
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
    if (z[i] > NEAR_ZERO)
    {
      FILE_LOG(LOG_CONSTRAINT) << " -- contact " << i << " is active" << std::endl;
      mapping[epd.N_ACT_CONTACTS++] = i;
    }

  // 2. compute inactive indices
  for (unsigned i=0, j=epd.N_ACT_CONTACTS; i< epd.N_CONTACTS; i++)
    if (z[i] < NEAR_ZERO)
      mapping[j++] = i;

  // permute inactive indices
  std::random_shuffle(mapping.begin()+epd.N_ACT_CONTACTS, mapping.end());  

  // set solution vector to reflect indices in active set
  for (unsigned i=0; i< epd.N_ACT_CONTACTS; i++)
    z[i] = z[mapping[i]];
  std::fill(z.row_iterator_begin()+epd.N_ACT_CONTACTS, z.row_iterator_begin()+epd.N_CONTACTS, 0.0);
  FILE_LOG(LOG_CONSTRAINT) << "-- permuted frictionless lcp solution: " << z << std::endl;

  // permute contact constraints
  std::vector<UnilateralConstraint*> new_contact_constraints(epd.contact_constraints.size());
  permute(mapping.begin(), mapping.end(), epd.contact_constraints.begin(), new_contact_constraints.begin());
  epd.contact_constraints = new_contact_constraints;

  // TODO: add constraint computation and cross computation methods to Joint

  // get iterators to the proper matrices
  RowIteratord CnCn = epd.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = epd.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = epd.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCs = epd.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = epd.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCt = epd.Ct_iM_CtT.row_iterator_begin();

  // process contact constraints, setting up matrices
  for (unsigned i=0; i< epd.contact_constraints.size(); i++) 
  {
    // compute cross constraint data for contact constraints
    for (unsigned j=0; j< epd.contact_constraints.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 3);

      // check whether i==j (single contact constraint)
      if (i == j)
      {
        // compute matrix / vector for contact constraint i
        _v.set_zero(3);
        epd.contact_constraints[i]->compute_constraint_data(_MM, _v);

        // setup appropriate parts of contact inertia matrices
        RowIteratord_const data = _MM.row_iterator_begin();
        *CnCn = *data++;
        *CnCs = *data++;
        *CnCt = *data; data += 2; // advance past Cs_iM_CnT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_iM_CtT
        *CtCt = *data;

        // setup appropriate parts of contact velocities
        data = _v.row_iterator_begin();
        epd.Cn_v[i] = *data++;
        epd.Cs_v[i] = *data++;
        epd.Ct_v[i] = *data;
      }
      else
      {
        // compute matrix for cross constraint
        epd.contact_constraints[i]->compute_cross_constraint_data(*epd.contact_constraints[j], _MM);

        // setup appropriate parts of contact inertia matrices
        RowIteratord_const data = _MM.row_iterator_begin();
        *CnCn = *data++;
        *CnCs = *data++;
        *CnCt = *data; data += 2; // advance to Cs_iM_CsT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_iM_CtT
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
    for (unsigned j=0; j< epd.limit_constraints.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 1);

      // compute matrix for cross constraint
      epd.contact_constraints[i]->compute_cross_constraint_data(*epd.limit_constraints[j], _MM);

      // setup appropriate parts of contact / limit inertia matrices
      ColumnIteratord_const data = _MM.column_iterator_begin();
      epd.Cn_iM_LT(i,j) = *data++;
      epd.Cs_iM_LT(i,j) = *data++;
      epd.Ct_iM_LT(i,j) = *data; 
    }
  }

  // NOTE: no need to update limit constraints of form L_iM_X 
  //       (cross data has already been computed for contact/limit constraints)
  
  // update cn, l, alpha and kappa 
  z.get_sub_vec(epd.CN_IDX, epd.CS_IDX, epd.cn);
  z.get_sub_vec(epd.L_IDX, epd.ALPHA_X_IDX, epd.l);
  z.get_sub_vec(epd.ALPHA_X_IDX, epd.N_VARS, epd.alpha_x);

  // mark active set as contact constraint set
  epd.N_CONTACT_CONSTRAINTS = epd.N_ACT_CONTACTS;
  epd.contact_constraint_set.resize(epd.N_CONTACTS);
  for (unsigned i=0; i< epd.N_ACT_CONTACTS; i++)
    epd.contact_constraint_set[i] = true;
  for (unsigned i=epd.N_ACT_CONTACTS; i< epd.N_CONTACTS; i++)
    epd.contact_constraint_set[i] = false;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::permute_problem() entered" << std::endl;
}

/**
 * Applies method of Drumwright and Shell to a set of connected constraints
 * \param constraints a set of connected constraints 
 */
void ImpactConstraintHandler::apply_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints, double inv_dt)
{
  double ke_minus = 0.0, ke_plus = 0.0;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd, inv_dt);

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->id << " pre-constraint handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // solve the (non-frictional) linear complementarity problem to determine
  // the kappa constant
  solve_frictionless_lcp(_epd, _z);

  // permute the problem
  permute_problem(_epd, _z);

  // mark all contacts as active
  _epd.N_ACT_CONTACTS = _epd.N_CONTACTS;
  _epd.N_ACT_K = _epd.N_K_TOTAL;

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
      // need to solve another impact problem 
      solve_frictionless_lcp(_epd, _z);

      // update constraint problem data and z
      permute_problem(_epd, _z);

      // determine N_ACT_K
      _epd.N_ACT_CONTACTS = _epd.N_CONTACTS;
      _epd.N_ACT_K = _epd.N_K_TOTAL;

      // use QP / NQP solver with warm starting to find the solution
      if (use_qp_solver(_epd))
        solve_qp(_z, _epd);
      else
        solve_nqp(_z, _epd);

      // update the impulses from z 
      update_from_stacked(_epd, _z);
    }  
  }

  // apply impulses 
  apply_impulses(_epd);

  // compute energy
  if (LOGGING(LOG_CONSTRAINT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << _epd.super_bodies[i]->id << " post-constraint handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONSTRAINT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

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

/// Applies impulses to bodies
void ImpactConstraintHandler::apply_impulses(const UnilateralConstraintProblemData& q)
{
  map<DynamicBodyPtr, VectorNd> gj;
  map<DynamicBodyPtr, VectorNd>::iterator gj_iter;

  // loop over all contact constraints first
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& e = *q.contact_constraints[i];
    SForced w(e.contact_impulse);

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

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

    // apply limit impulses to bodies in independent coordinates
    if (dynamic_pointer_cast<RCArticulatedBody>(ab))
    {
      // get the index of the joint
      unsigned idx = e.limit_joint->get_coord_index() + e.limit_dof;

      // initialize the vector if necessary
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBody::eSpatial));
        gj_iter = gj.find(ab);
      }

      // set the limit force
      gj_iter->second[idx] += e.limit_impulse;
    }
    else
    {
      // TODO: handle bodies in absolute coordinates here
      assert(false);
    }
  }

  // TODO: apply constraint impulses

  // apply all generalized impacts
  for (map<DynamicBodyPtr, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/// Computes the data to the LCP / QP problems
void ImpactConstraintHandler::compute_problem_data(UnilateralConstraintProblemData& q, double inv_dt)
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
    q.N_GC += q.super_bodies[i]->num_generalized_coordinates(DynamicBody::eSpatial);

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
  q.Cn_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_iM_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cn_iM_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Cs_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cs_iM_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cs_iM_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.Ct_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Ct_iM_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Ct_iM_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_IMP);
  q.L_iM_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.L_iM_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_IMP);
  q.Jx_iM_JxT.set_zero(q.N_CONSTRAINT_EQNS_IMP, q.N_CONSTRAINT_EQNS_IMP);
  q.Cn_v.set_zero(q.N_CONTACTS);
  q.Cs_v.set_zero(q.N_CONTACTS);
  q.Ct_v.set_zero(q.N_CONTACTS);
  q.L_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_IMP);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  q.l.set_zero(q.N_LIMITS);
  q.alpha_x.set_zero(q.N_CONSTRAINT_EQNS_IMP);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
  q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
  q.NCT_IDX = q.NCS_IDX + q.N_LIN_CONE;
  q.L_IDX = q.NCT_IDX + q.N_LIN_CONE;
  q.ALPHA_X_IDX = q.L_IDX + q.N_LIMITS;
  q.N_VARS = q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_IMP;

  // TODO: add constraint computation and cross computation methods to Joint

  // get iterators to the proper matrices
  RowIteratord CnCn = q.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = q.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = q.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCs = q.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = q.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCt = q.Ct_iM_CtT.row_iterator_begin();

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
        *CnCt = *data; data += 2; // advance past Cs_iM_CnT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_iM_CtT
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
        *CnCt = *data; data += 2; // advance to Cs_iM_CsT
        *CsCs = *data++;
        *CsCt = *data; data += 3; // advance to Ct_iM_CtT
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
      q.Cn_iM_LT(i,j) = *data++;
      q.Cs_iM_LT(i,j) = *data++;
      q.Ct_iM_LT(i,j) = *data; 
    }
  }

  // process limit constraints, setting up matrices
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    // compute matrix / vector for contact constraint i
    q.limit_constraints[i]->compute_constraint_data(_MM, _v);

    // setup appropriate entry of limit inertia matrix and limit velocity
    q.L_iM_LT(i,i) = _MM.data()[0];
    q.L_v[i] = _v.data()[0];

    // compute cross/cross limit constraint data
    for (unsigned j=i+1; j< q.limit_constraints.size(); j++)
    {
      // reset _MM
      _MM.resize(1,1);

      // compute matrix for cross constraint
      q.limit_constraints[i]->compute_cross_constraint_data(*q.limit_constraints[j], _MM);

      // setup appropriate part of limit / limit inertia matrix
      q.L_iM_LT(i,j) = q.L_iM_LT(j,i) = _MM.data()[0];
    }

    // NOTE: cross data has already been computed for contact/limit constraints
  }

  // correct constraint violations on contacts
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
    q.Cn_v[i] += q.contact_constraints[i]->signed_violation * inv_dt;

  // correct constraint violations on limits
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
    q.L_v[i] += q.limit_constraints[i]->signed_violation * inv_dt;
}

/// Solves the viscous friction LCP
void ImpactConstraintHandler::apply_visc_friction_model(UnilateralConstraintProblemData& q)
{
  // compute the (Coulomb) frictionless LCP
  VectorNd z;
  solve_frictionless_lcp(q, z);

  // setup impulses 
  q.cn = z.segment(q.CN_IDX, q.N_CONTACTS);
  q.l = z.segment(q.L_IDX, q.L_IDX+q.N_LIMITS);
  q.alpha_x = z.segment(q.ALPHA_X_IDX, q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_IMP);
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

  // TODO: setup joint constraint impulses here

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_visc_friction_model() exited" << std::endl;
}

/// Solves the no-slip model LCP
void ImpactConstraintHandler::apply_no_slip_model(UnilateralConstraintProblemData& q)
{
  std::vector<unsigned> J_indices, S_indices, T_indices;
  const unsigned NCONTACTS = q.N_CONTACTS;
  const unsigned NLIMITS = q.N_LIMITS;
  const unsigned NIMP = q.N_CONSTRAINT_EQNS_IMP;
  const unsigned N_IDX = 0;
  const unsigned L_IDX = N_IDX + NCONTACTS;
  VectorNd lb, ub, b;
  MatrixNd A;

  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cs': " << std::endl << q.Cn_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Ct': " << std::endl << q.Cn_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Cs': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * Ct': " << std::endl << q.Cs_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cs * inv(M) * L': " << std::endl << q.Cs_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * Ct': " << std::endl << q.Ct_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "  Ct * inv(M) * L': " << std::endl << q.Ct_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "  L * inv(M) * L': " << std::endl << q.L_iM_LT;
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
  // find largest non-singular set of J, S, and T indices 
  // ********************************************************

  // loop through joint constraints, forming J*inv(M)*J' and checking condition
  for (unsigned i=0; i< NIMP; i++)
  {
    // add the index tentatively to the set
    J_indices.push_back(i);

    // select the rows and columns
    q.Jx_iM_JxT.select_square(J_indices.begin(), J_indices.end(), _Y);

    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< J_indices.size(); j++)
      _Y(j,j) -= NEAR_ZERO;
    
    // attempt Cholesky factorization
    if (!_LA.factor_chol(_Y))
      J_indices.pop_back();
  }
  
  // get the reduced Jx*iM*Jx' matrix  
  q.Jx_iM_JxT.select_square(J_indices.begin(), J_indices.end(), _rJx_iM_JxT);

  // loop through contacts, forming matrix below and checking its condition 
  // | S*inv(M)*S'  S*inv(M)*T' S*inv(M)*J' |
  // | T*inv(M)*S'  T*inv(M)*T' T*inv(M)*J' |
  // | J*inv(M)*S'  J*inv(M)*T' J*inv(M)*J' |
  bool last_success = false;
  for (unsigned i=0; i< NCONTACTS; i++)
  {
    // update S indices
    S_indices.push_back(i);
    
    // setup indices
    unsigned S_IDX = 0;
    unsigned T_IDX = S_indices.size();
    unsigned J_IDX = T_IDX + T_indices.size();
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());

    // add S/S, T/T, J/J components to 'check' matrix
    q.Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, S_IDX, _MM);
    q.Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, T_IDX, _MM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);

    // add S/T components to 'check' matrix
    q.Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, T_IDX, _MM);
    _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

    // add S/J components to check matrix
    q.Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, S_IDX, _MM, Ravelin::eTranspose);
 
    // add T/J components to check matrix
    q.Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, T_IDX, _MM, Ravelin::eTranspose);
    
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
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());

    // add S/S, T/T, J/J components to 'check' matrix
    q.Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, S_IDX, _MM);
    q.Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, T_IDX, _MM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);

    // add S/T components to 'check' matrix
    q.Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, T_IDX, _MM);
    _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

    // add S/J components to check matrix
    q.Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, S_IDX, _MM, Ravelin::eTranspose);
 
    // add T/J components to check matrix
    q.Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, T_IDX, _MM, Ravelin::eTranspose);
 
    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< _Y.rows(); j++)
      _Y(j,j) -= NEAR_ZERO;

    // see whether check matrix can be Cholesky factorized
    last_success = _LA.factor_chol(_Y);
    if (!last_success)
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
  // reform Y if necessary 
  // ********************************************************

  // setup indices
  const unsigned S_IDX = 0;
  const unsigned T_IDX = S_indices.size();
  const unsigned J_IDX = T_IDX + T_indices.size();
  if (!last_success)
  {
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());

    // add S/S, T/T, J/J components to X 
    q.Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, S_IDX, _MM);
    q.Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, T_IDX, _MM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);

    // add S/T components to X 
    q.Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, T_IDX, _MM);
    _Y.set_sub_mat(T_IDX, S_IDX, _MM, Ravelin::eTranspose);

    // add S/J components to X 
    q.Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(S_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, S_IDX, _MM, Ravelin::eTranspose);
 
    // add T/J components to X 
    q.Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _MM);
    _Y.set_sub_mat(T_IDX, J_IDX, _MM);
    _Y.set_sub_mat(J_IDX, T_IDX, _MM, Ravelin::eTranspose);
  
    // do the Cholesky factorization (should not fail) 
    bool success = _LA.factor_chol(_Y);
    assert(success);
  }

  // defining Y = inv(X*inv(M)*X') and Q = [Cn; L; 0]
  // and using the result above yields following LCP:
  // matrix: Q*inv(A)*Q' = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  // vector: Q*inv(A)*a  = Q*v - Q*inv(M)*X'*Y*X*v

  // setup Q*inv(M)*Q' 
  _MM.set_zero(q.N_CONTACTS + q.N_LIMITS, q.N_CONTACTS + q.N_LIMITS);
  _MM.set_sub_mat(N_IDX, N_IDX, q.Cn_iM_CnT);
  _MM.set_sub_mat(N_IDX, L_IDX, q.Cn_iM_LT);
  _MM.set_sub_mat(L_IDX, N_IDX, q.Cn_iM_LT, Ravelin::eTranspose);
  _MM.set_sub_mat(L_IDX, L_IDX, q.L_iM_LT);

  // setup Q*inv(M)*X'
  _Q_iM_XT.resize(q.N_CONTACTS + q.N_LIMITS, S_indices.size() + T_indices.size() + J_indices.size());
  q.Cn_iM_CsT.select_columns(S_indices.begin(), S_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, S_IDX, _workM);
  q.Cn_iM_CtT.select_columns(T_indices.begin(), T_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, T_IDX, _workM);
  q.Cn_iM_JxT.select_columns(J_indices.begin(), J_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, J_IDX, _workM);
  q.Cs_iM_LT.select_rows(S_indices.begin(), S_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(L_IDX, S_IDX, _workM, Ravelin::eTranspose);
  q.Ct_iM_LT.select_rows(T_indices.begin(), T_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(L_IDX, T_IDX, _workM, Ravelin::eTranspose);
  q.L_iM_JxT.select_columns(J_indices.begin(), J_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(L_IDX, J_IDX, _workM);

  // compute Y*X*inv(M)*Q'
  MatrixNd::transpose(_Q_iM_XT, _workM);
  _LA.solve_chol_fast(_Y, _workM);

  // compute Q*inv(M)*X'*Y*X*inv(M)*Q'
  _Q_iM_XT.mult(_workM, _workM2);
  _MM -= _workM2;

  // setup -Q*v
  _qq.resize(q.N_CONTACTS + q.N_LIMITS);
  _qq.set_sub_vec(N_IDX, q.Cn_v);
  _qq.set_sub_vec(L_IDX, q.L_v);

  // setup X*v
  _Xv.resize(S_indices.size() + T_indices.size() + J_indices.size());
  q.Cs_v.select(S_indices.begin(), S_indices.end(), _workv);
  _Xv.set_sub_vec(S_IDX, _workv);
  q.Ct_v.select(T_indices.begin(), T_indices.end(), _workv);
  _Xv.set_sub_vec(T_IDX, _workv);
  q.Jx_v.select(J_indices.begin(), J_indices.end(), _workv);
  _Xv.set_sub_vec(J_IDX, _workv);

  // compute Y*X*v
  _YXv = _Xv;
  _LA.solve_chol_fast(_Y, _YXv);

  // compute Q*inv(M)*X' * Y*X*v
  _Q_iM_XT.mult(_YXv, _workv);

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
  _Q_iM_XT.transpose_mult(_v, _workv);
  _LA.solve_chol_fast(_Y, _workv);
  _cs_ct_alphax += _workv;
  _cs_ct_alphax.negate();

  // setup impulses 
  q.cn = _v.segment(0, q.N_CONTACTS);
  q.l = _v.segment(q.N_CONTACTS, _v.size());
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  q.alpha_x.set_zero(NIMP);
  SharedConstVectorNd cs_vec = _cs_ct_alphax.segment(S_IDX, T_IDX);
  SharedConstVectorNd ct_vec = _cs_ct_alphax.segment(T_IDX, J_IDX);
  SharedConstVectorNd alphax_vec = _cs_ct_alphax.segment(J_IDX, _cs_ct_alphax.size()); 
  q.cs.set(S_indices.begin(), S_indices.end(), cs_vec);
  q.ct.set(T_indices.begin(), T_indices.end(), ct_vec);
  q.alpha_x.set(J_indices.begin(), J_indices.end(), alphax_vec);

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

  // TODO: setup joint constraint impulses here

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
  _A = q.Jx_iM_JxT; 
  _LA.svd(_A, _AU, _AS, _AV);

  // setup the B matrix
  // B = [ Cn; L ]*inv(M)*[ Cn' L' ]
  _B.resize(NCONTACTS+NLIMITS, NCONTACTS+NLIMITS);
  _B.set_sub_mat(0, 0, q.Cn_iM_CnT);  
  _B.set_sub_mat(0, NCONTACTS, q.Cn_iM_LT);
  _B.set_sub_mat(NCONTACTS, 0, q.Cn_iM_LT, Ravelin::eTranspose);
  _B.set_sub_mat(NCONTACTS, NCONTACTS, q.L_iM_LT);

  // setup the C matrix and compute inv(A)*C
  // C = Jx*inv(M)*[ Cn' L' ]; note: D = C'
  _C.resize(NIMP, NCONTACTS+NLIMITS);
  _C.set_sub_mat(0,0, q.Cn_iM_JxT, Ravelin::eTranspose);
  _C.set_sub_mat(0,NCONTACTS, q.L_iM_JxT, Ravelin::eTranspose);
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
  q.Cn_iM_CsT.mult(_cs_visc, _workv);
  bsub -= _workv;
  q.Cn_iM_CtT.mult(_ct_visc, _workv);
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
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP vector: " << _qq << std::endl;

  // solve the LCP
  if (!_lcp.lcp_fast(_MM, _qq, _v) && !_lcp.lcp_lemke_regularized(_MM, _qq, _v))
    throw std::runtime_error("Unable to solve constraint LCP!");

  // compute alphax
  // u = -inv(A)*(a + Cv)
  _C.mult(_v, _alpha_x) += _a;
  _alpha_x.negate();   

  // determine the value of kappa
  SharedConstVectorNd cn = _v.segment(0, q.N_CONTACTS);
  SharedConstVectorNd l = _v.segment(q.N_CONTACTS, _v.size());
  q.Cn_iM_CnT.mult(cn, _Cn_vplus) += q.Cn_v;
  q.kappa = _Cn_vplus.norm1();

  // setup the homogeneous solution
  z.set_zero(q.N_VARS);
  z.set_sub_vec(q.CN_IDX, cn);
  z.set_sub_vec(q.L_IDX, l);
  z.set_sub_vec(q.ALPHA_X_IDX, _alpha_x);

  FILE_LOG(LOG_CONSTRAINT) << "  LCP result: " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  kappa: " << q.kappa << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_lcp() exited" << std::endl;
}

/// Gets the super body (articulated if any)
DynamicBodyPtr ImpactConstraintHandler::get_super_body(SingleBodyPtr sb)
{
  ArticulatedBodyPtr ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return sb;
}

