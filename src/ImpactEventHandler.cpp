/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/ArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/Event.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/LinAlg.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/Optimization.h>
#include <Moby/ImpactEventHandler.h>

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
ImpactEventHandler::ImpactEventHandler()
{
  ip_max_iterations = 100;
  ip_eps = 1e-6;
  use_ip_solver = false;
  poisson_eps = NEAR_ZERO;
}

// Processes impacts
void ImpactEventHandler::process_events(const vector<Event>& events, Real tol)
{
  FILE_LOG(LOG_CONTACT) << "*************************************************************";
  FILE_LOG(LOG_CONTACT) << endl;
  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::process_events() entered";
  FILE_LOG(LOG_CONTACT) << endl;
  FILE_LOG(LOG_CONTACT) << "*************************************************************";
  FILE_LOG(LOG_CONTACT) << endl;

  // apply the method to all contacts
  if (!events.empty())
  {
    vector<Event> nev;
    nev.push_back(events.front());
    apply_model(events, tol);
  }
  else
    FILE_LOG(LOG_CONTACT) << " (no events?!)" << endl;
    
  FILE_LOG(LOG_CONTACT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::process_events() exited" << endl;
  FILE_LOG(LOG_CONTACT) << "*************************************************************" << endl;
}

/// Applies the model to a set of events 
/**
 * \param events a set of events
 */
void ImpactEventHandler::apply_model(const vector<Event>& events, Real tol) const
{
  // **********************************************************
  // determine sets of connected events 
  // **********************************************************
  list<list<Event*> > groups;
  Event::determine_connected_events(events, groups);
  Event::remove_nonimpacting_groups(groups, tol);

  // **********************************************************
  // do method for each connected set 
  // **********************************************************
  for (list<list<Event*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
    list<Event*>& revents = *i;
    apply_model_to_connected_events(revents);
  }
}

/**
 * Applies method of Drumwright and Shell to a set of connected events
 * \param events a set of connected events 
 */
void ImpactEventHandler::apply_model_to_connected_events(const list<Event*>& events) const
{
  Real ke_minus = 0.0, ke_plus = 0.0;
  vector<Event> constraint_event_objects;
  SAFESTATIC EventProblemData epd;

  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::apply_model_to_connected_events() entered" << endl;

  // reset problem data
  epd.reset();

  // determine sets of contact and limit events
  partition_events(events, epd.contact_events, epd.limit_events);

  // add events for constraints for articulated bodies
  add_constraint_events(events, constraint_event_objects, epd.constraint_events);

  // compute all event cross-terms
  compute_problem_data(epd);

  // compute energy
  if (LOGGING(LOG_CONTACT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      Real ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONTACT) << "  body " << epd.super_bodies[i]->id << " pre-event handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

// NOTE: we disable this per Ruina's suggestion
/*
  // solve the (non-frictional) linear complementarity problem to determine
  // the kappa constant
  solve_lcp(epd);
*/
  epd.kappa = (Real) -std::numeric_limits<float>::max();

  // determine what type of QP solver to use
  if (use_qp_solver(epd))
    solve_qp(epd, poisson_eps);
  else
    solve_nqp(epd, poisson_eps);

  // set new generalized velocities 
  set_generalized_velocities(epd);

  // compute energy
  if (LOGGING(LOG_CONTACT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      Real ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONTACT) << "  body " << epd.super_bodies[i]->id << " post-event handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONTACT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::apply_model_to_connected_events() exiting" << endl;
}

/// Determines whether we can use the QP solver
bool ImpactEventHandler::use_qp_solver(const EventProblemData& epd)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // first, check whether any contact events use a true friction cone
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
    if (epd.contact_events[i]->contact_NK == UINF || 
        epd.contact_events[i]->contact_NK < 2)
      return false;

  // now, check whether any articulated bodies use the advanced friction
  // model
  for (unsigned i=0; i< epd.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(epd.super_bodies[i]);
    if (abody && abody->use_advanced_friction_model)
      return false;
  }

  // still here? ok to use QP solver
  return true;
}

/// Adds constraint events
void ImpactEventHandler::add_constraint_events(const list<Event*>& events, vector<Event>& constraint_events_objects, vector<Event*>& constraint_events)
{
  // clear the vectors
  constraint_events_objects.clear();
  constraint_events.clear();

  // determine the articulated bodies
  vector<ArticulatedBodyPtr> abodies;
  BOOST_FOREACH(const Event* e, events)
  {
    if (e->event_type == Event::eContact)
    {
      SingleBodyPtr sb1 = e->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = e->contact_geom2->get_single_body();
      ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
      ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
      if (ab1)
        abodies.push_back(ab1);
      if (ab2)
        abodies.push_back(ab2);
    }
    else if (e->event_type == Event::eLimit)
    {
      RigidBodyPtr rb = e->limit_joint->get_outboard_link();
      abodies.push_back(rb->get_articulated_body());
    }
    else
      assert(false);
  }

  // make the vector of articulated bodies unique
  std::sort(abodies.begin(), abodies.end());
  abodies.erase(std::unique(abodies.begin(), abodies.end()), abodies.end());

  // determine the constraint events
  for (unsigned i=0; i< abodies.size(); i++)
    abodies[i]->get_constraint_events(constraint_events_objects);

  // add the pointers to the constraint event objects to the constraint events
  constraint_events.resize(constraint_events_objects.size());
  for (unsigned i=0; i< constraint_events_objects.size(); i++)
    constraint_events[i] = &constraint_events_objects[i];
}

/// Partitions the events into contact and limit events
void ImpactEventHandler::partition_events(const list<Event*>& all, vector<Event*>& contacts, vector<Event*>& limits)
{
  contacts.clear();
  limits.clear();

  BOOST_FOREACH(Event* e, all)
  {
    if (e->event_type == Event::eContact)
      contacts.push_back(e);
    else
    {
      assert(e->event_type == Event::eLimit);
      limits.push_back(e);
    }
  }
}

/// Determines and sets the new generalized velocities
void ImpactEventHandler::set_generalized_velocities(const EventProblemData& q)
{
  // determine the change in generalized velocities
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    q.super_bodies[i]->update_velocity(q); 
}

/// Computes the data to the LCP / QP problems
void ImpactEventHandler::compute_problem_data(EventProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // first, setup all tangent directions for contact events
  for (unsigned i=0; i< q.contact_events.size(); i++)
    q.contact_events[i]->determine_contact_tangents();

  // determine set of "super" bodies from contact events
  q.super_bodies.clear();
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    q.super_bodies.push_back(get_super_body(q.contact_events[i]->contact_geom1->get_single_body()));
    q.super_bodies.push_back(get_super_body(q.contact_events[i]->contact_geom2->get_single_body()));
  }

  // determine set of "super" bodies from limit events
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    RigidBodyPtr outboard = q.limit_events[i]->limit_joint->get_outboard_link();
    q.super_bodies.push_back(get_super_body(outboard));
  }

  // determine set of "super" bodies from constraint events
  for (unsigned i=0; i< q.constraint_events.size(); i++)
  {
    RigidBodyPtr outboard = q.constraint_events[i]->constraint_joint->get_outboard_link();
    q.super_bodies.push_back(get_super_body(outboard));
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // initialize constants and set easy to set constants
  q.N_CONSTRAINTS = q.constraint_events.size();
  q.N_CONTACTS = q.contact_events.size();
  q.N_LIMITS = q.limit_events.size();

  // setup constants related to articulated bodies
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
    if (abody) {
      q.N_CONSTRAINT_EQNS_EXP += abody->num_constraint_eqns_explicit();
      if (abody->use_advanced_friction_model)
      {
        q.N_CONSTRAINT_DOF_IMP += abody->num_joint_dof_implicit();
        q.N_CONSTRAINT_DOF_EXP += abody->num_joint_dof_explicit();
      }
    }
  }

  // determine which contact constraints use a true friction cone
  for (unsigned i=0; i< q.contact_events.size(); i++)
    if (q.contact_events[i]->contact_NK < UINF &&
        q.contact_events[i]->contact_NK >= 2)
      q.N_K_TOTAL += q.contact_events[i]->contact_NK;

  // determine number of true friction constraints
  for (unsigned i=0; i< q.N_CONTACTS; i++)
    if (q.contact_events[i]->contact_NK == UINF ||
        q.contact_events[i]->contact_NK < 2)
      q.N_TRUE_CONE++;

  // determine number of linearized friction constraints
  q.N_LIN_CONE = q.N_CONTACTS - q.N_TRUE_CONE;

  // initialize the problem matrices / vectors
  q.Jc_iM_JcT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Jc_iM_DcT.set_zero(q.N_CONTACTS, q.N_CONTACTS*2);
  q.Jc_iM_JlT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Jc_iM_DtT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_DOF_IMP);
  q.Jc_iM_JxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_EQNS_EXP);
  q.Jc_iM_DxT.set_zero(q.N_CONTACTS, q.N_CONSTRAINT_DOF_EXP);
  q.Dc_iM_DcT.set_zero(q.N_CONTACTS*2, q.N_CONTACTS*2);
  q.Dc_iM_JlT.set_zero(q.N_CONTACTS*2, q.N_LIMITS);
  q.Dc_iM_DtT.set_zero(q.N_CONTACTS*2, q.N_CONSTRAINT_DOF_IMP);
  q.Dc_iM_JxT.set_zero(q.N_CONTACTS*2, q.N_CONSTRAINT_EQNS_EXP);
  q.Dc_iM_DxT.set_zero(q.N_CONTACTS*2, q.N_CONSTRAINT_DOF_EXP);
  q.Jl_iM_JlT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.Jl_iM_DtT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_DOF_IMP);
  q.Jl_iM_JxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_EQNS_EXP);
  q.Jl_iM_DxT.set_zero(q.N_LIMITS, q.N_CONSTRAINT_DOF_EXP);
  q.Dt_iM_DtT.set_zero(q.N_CONSTRAINT_DOF_IMP, q.N_CONSTRAINT_DOF_IMP);
  q.Dt_iM_JxT.set_zero(q.N_CONSTRAINT_DOF_IMP, q.N_CONSTRAINT_EQNS_EXP);
  q.Dt_iM_DxT.set_zero(q.N_CONSTRAINT_DOF_IMP, q.N_CONSTRAINT_DOF_EXP);
  q.Jx_iM_JxT.set_zero(q.N_CONSTRAINT_EQNS_EXP, q.N_CONSTRAINT_EQNS_EXP);
  q.Jx_iM_DxT.set_zero(q.N_CONSTRAINT_EQNS_EXP, q.N_CONSTRAINT_DOF_EXP);
  q.Dx_iM_DxT.set_zero(q.N_CONSTRAINT_DOF_EXP, q.N_CONSTRAINT_DOF_EXP);
  q.Jc_v.set_zero(q.N_CONTACTS);
  q.Dc_v.set_zero(q.N_CONTACTS*2);
  q.Jl_v.set_zero(q.N_LIMITS);
  q.Jx_v.set_zero(q.N_CONSTRAINT_EQNS_EXP);
  q.Dx_v.set_zero(q.N_CONSTRAINT_DOF_EXP);
  q.alpha_c.set_zero(q.N_CONTACTS);
  q.beta_c.set_zero(q.N_CONTACTS*2);
  q.alpha_l.set_zero(q.N_LIMITS);
  q.beta_t.set_zero(q.N_CONSTRAINT_DOF_IMP);
  q.alpha_x.set_zero(q.N_CONSTRAINT_EQNS_EXP);
  q.beta_x.set_zero(q.N_CONSTRAINT_DOF_EXP);

  // for each super body, update the problem data
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    q.super_bodies[i]->update_event_data(q);
}

/// Sets up optimization data for nonlinear QP
void ImpactEventHandler::set_optimization_data(EventProblemData& q, ImpactOptData& iopt)
{
  vector<vector<unsigned> > loop_links;

  // setup some needed constants
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned N_SUPERS = q.super_bodies.size();
  const unsigned BETA_C_IDX = q.N_CONTACTS;
  const unsigned NBETA_C_IDX = BETA_C_IDX + q.N_CONTACTS*2;
  const unsigned ALPHA_L_IDX = NBETA_C_IDX + q.N_LIN_CONE*2;
  const unsigned ALPHA_X_IDX = ALPHA_L_IDX + q.N_LIMITS;
  const unsigned BETA_T_IDX = ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_EXP;

  // setup iopt data
  iopt.loop_indices.resize(N_SUPERS);
  iopt.Z.resize(N_SUPERS);
  iopt.Zd.resize(N_SUPERS);
  iopt.Z1d.resize(N_SUPERS);
  iopt.cone_contacts.resize(q.N_TRUE_CONE);
  iopt.c_mu_c.resize(q.N_TRUE_CONE);
  iopt.c_visc.resize(q.N_TRUE_CONE);
  iopt.j_mu_c.clear();
  iopt.j_visc.clear();

  // clear vectors
  for (unsigned i=0; i< N_SUPERS; i++)
  {
    iopt.loop_indices[i].clear();
    iopt.Z[i].clear();
    iopt.Zd[i].clear();
    iopt.Z1d[i].clear();
  }

  // determine which contact constraints use a true friction cone
  for (unsigned i=0, j=0; i< q.contact_events.size(); i++)
    if (q.contact_events[i]->contact_NK == UINF ||
        q.contact_events[i]->contact_NK < 2)
    {
      iopt.cone_contacts[j] = i;
      iopt.c_mu_c[j] = sqr(q.contact_events[i]->contact_mu_coulomb);
      iopt.c_visc[j] = (sqr(q.Dc_v[i*2]) + sqr(q.Dc_v[i*2+1])) *
                       sqr(q.contact_events[i]->contact_mu_viscous);
      j++;
    }

  // init vectors
  iopt.body_indices.clear();
  iopt.delta_start.resize(N_SUPERS);
  iopt.joint_friction_start.resize(N_SUPERS);
  iopt.implicit_start.resize(N_SUPERS);
  iopt.explicit_start.resize(N_SUPERS);
  iopt.true_indices.resize(N_SUPERS);
  iopt.alpha_c_indices.resize(N_SUPERS);
  iopt.beta_nbeta_c_indices.resize(N_SUPERS);
  iopt.alpha_l_indices.resize(N_SUPERS);
  iopt.beta_t_indices.resize(N_SUPERS);
  iopt.beta_x_indices.resize(N_SUPERS);
  iopt.j_mu_c.clear();
  iopt.j_visc.clear();

  // setup the current joint indices
  unsigned st_idx = 0;

  // set current index starts
  unsigned implicit_start = 0, explicit_start = 0, joint_friction_start = 0;

  // determine the mapping from contacts with true friction cones to *all*
  // contact indices
  for (unsigned i=0; i< N_SUPERS; i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
    if (abody)
    {
      // update number of kinematic loops, Z variables, and numbers of joint
      // DOFs if using advanced friction model
      if (abody->use_advanced_friction_model)
      {
        // search for this articulated body in contact events to setup indices
        iopt.alpha_c_indices[i].clear();
        iopt.beta_nbeta_c_indices[i].clear();
        for (unsigned j=0, r=BETA_C_IDX+1, s=NBETA_C_IDX+1; j< q.N_CONTACTS; j++, r+=2)
        {
          DynamicBodyPtr d1 = get_super_body(q.contact_events[j]->contact_geom1->get_single_body());
          DynamicBodyPtr d2 = get_super_body(q.contact_events[j]->contact_geom2->get_single_body());
          if (d1 == abody)
          {
            iopt.alpha_c_indices[i].push_back((int) j+1);
            iopt.beta_nbeta_c_indices[i].push_back((int) r);
            iopt.beta_nbeta_c_indices[i].push_back((int) r+1);
            if (q.contact_events[j]->contact_NK < UINF && 
                q.contact_events[j]->contact_NK >= 2)
            {
              iopt.beta_nbeta_c_indices[i].push_back((int) s);
              iopt.beta_nbeta_c_indices[i].push_back((int) s+1);
            }
          }
          if (d2 == abody)
          {
            iopt.alpha_c_indices[i].push_back((int) -j-1);
            iopt.beta_nbeta_c_indices[i].push_back((int) -r);
            iopt.beta_nbeta_c_indices[i].push_back((int) -r-1);
            if (q.contact_events[j]->contact_NK < UINF &&
                q.contact_events[j]->contact_NK >= 2)
            {
              iopt.beta_nbeta_c_indices[i].push_back((int) -s);
              iopt.beta_nbeta_c_indices[i].push_back((int) -s-1);
            }
          }

          // see whether to advance s
          if (q.contact_events[j]->contact_NK < UINF &&
              q.contact_events[j]->contact_NK >= 2)
            s += 2;
        }

        // NOTE: beta indices *must* be sorted for contact_select() to work
        std::sort(iopt.beta_nbeta_c_indices[i].begin(), iopt.beta_nbeta_c_indices[i].end());
 
        // search for this articulated body in limit events to setup indices
        iopt.alpha_l_indices[i].clear();
        for (unsigned j=0, r=ALPHA_L_IDX; j< q.N_LIMITS; j++, r++)
          if (q.limit_events[j]->limit_joint->get_articulated_body() == abody)
            iopt.alpha_l_indices[i].push_back(r);

        // setup beta_t indices
        const unsigned N_IMPLICIT_DOF = abody->num_joint_dof_implicit();
        iopt.beta_t_indices[i].clear();
        for (unsigned j=0; j< N_IMPLICIT_DOF; j++)
          iopt.beta_t_indices[i].push_back(BETA_T_IDX + implicit_start + j);

        // setup beta_x indices - note: these will have to be updated later...
        const unsigned N_EXPLICIT_DOF = abody->num_joint_dof_explicit();
        iopt.beta_x_indices[i].clear();
        for (unsigned j=0; j< N_EXPLICIT_DOF; j++)
          iopt.beta_x_indices[i].push_back(BETA_T_IDX + explicit_start + j);

        // starting indices
        iopt.implicit_start[i] = implicit_start;
        iopt.explicit_start[i] = explicit_start;
        iopt.joint_friction_start[i] = joint_friction_start;

        // update starting indices
        implicit_start += N_IMPLICIT_DOF;
        explicit_start += N_EXPLICIT_DOF;
        joint_friction_start += N_IMPLICIT_DOF + N_EXPLICIT_DOF;

        // update joint friction related stuff
        abody->find_loops(iopt.loop_indices[i], loop_links);
        abody->compute_Z_matrices(iopt.loop_indices[i], loop_links, iopt.Zd[i], iopt.Z1d[i], iopt.Z[i]);
        iopt.delta_start[i] = q.N_LOOPS;
        q.N_LOOPS += loop_links.size();  
        loop_links.clear();

        // update numbers of degrees-of-freedom
        q.N_CONSTRAINT_DOF_IMP += abody->num_joint_dof_implicit();
        q.N_CONSTRAINT_DOF_EXP += abody->num_joint_dof_explicit();
        const unsigned N_JOINT_DOF = N_IMPLICIT_DOF + abody->num_joint_dof_explicit();

        // setup true indices [implicit; explicit] and joint friction coeffs
        const vector<JointPtr>& joints = abody->get_joints();
        vector<unsigned>& true_indices = iopt.true_indices[i];
        true_indices.resize(N_JOINT_DOF);
        for (unsigned j=0, ki=st_idx, ke=st_idx+N_IMPLICIT_DOF; j< joints.size(); j++)
        {
          if (joints[j]->get_constraint_type() == Joint::eImplicit)
          {
            for (unsigned k=0; k< joints[j]->num_dof(); k++, ki++)
            {
              true_indices[ki] = j;
              iopt.j_mu_c[ki] = sqr(joints[j]->mu_fc);
              iopt.j_visc[ki] = sqr(joints[j]->qd[k] * joints[j]->mu_fv);
            }
          }
          else
          {
            assert(joints[j]->get_constraint_type() == Joint::eExplicit);
            for (unsigned k=0; k< joints[j]->num_dof(); k++, ke++)
            {
              true_indices[ke] = j;
              iopt.j_mu_c[ke] = sqr(joints[j]->mu_fc);
              iopt.j_visc[ke] = sqr(joints[j]->qd[k] * joints[j]->mu_fv);
            }
          }
        }

        // setup body indices
        for (unsigned j=0; j< N_JOINT_DOF; j++)
          iopt.body_indices.push_back(i);

        // update st_idx
        st_idx += N_JOINT_DOF;
      }
    }
  }

  // finally, update explicit friction start
  for (unsigned i=0; i< N_SUPERS; i++)
    iopt.explicit_start[i] += implicit_start;
}

/// Solves the quadratic program (potentially solves two QPs, actually)
void ImpactEventHandler::solve_qp(EventProblemData& q, Real poisson_eps)
{
  SAFESTATIC VectorN z, tmp, tmp2;
  const Real TOL = poisson_eps;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_CONSTRAINT_DOF_IMP = q.N_CONSTRAINT_DOF_IMP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // solve the QP
  solve_qp_work(q, z);

  // apply (Poisson) restitution to contacts
  for (unsigned i=0; i< N_CONTACTS; i++)
    z[i] *= ((Real) 1.0 + q.contact_events[i]->contact_epsilon);

  // apply (Poisson) restitution to limits
  for (unsigned i=0; i< N_LIMITS; i++)
    z[N_CONTACTS*5+i] *= ((Real) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  update_impulses(q, z);

  // update Jc_v, Dc_v, Jl_v, and Jx_v
  q.Jc_v += q.Jc_iM_JcT.mult(q.alpha_c, tmp);
  q.Jc_v += q.Jc_iM_DcT.mult(q.beta_c, tmp);
  q.Jc_v += q.Jc_iM_JlT.mult(q.alpha_l, tmp);
  q.Jc_v += q.Jc_iM_JxT.mult(q.alpha_x, tmp);
  q.Dc_v += q.Jc_iM_DcT.transpose_mult(q.alpha_c, tmp);
  q.Dc_v += q.Dc_iM_DcT.mult(q.beta_c, tmp);
  q.Dc_v += q.Dc_iM_JlT.mult(q.alpha_l, tmp);
  q.Dc_v += q.Dc_iM_JxT.mult(q.alpha_x, tmp);
  q.Jl_v += q.Jc_iM_JlT.transpose_mult(q.alpha_c, tmp);
  q.Jl_v += q.Dc_iM_JlT.transpose_mult(q.beta_c, tmp);
  q.Jl_v += q.Jl_iM_JlT.mult(q.alpha_l, tmp);
  q.Jl_v += q.Jl_iM_JxT.mult(q.alpha_x, tmp);
  q.Jx_v += q.Jc_iM_JxT.transpose_mult(q.alpha_c, tmp);
  q.Jx_v += q.Dc_iM_JxT.transpose_mult(q.beta_c, tmp);
  q.Jx_v += q.Jl_iM_JxT.transpose_mult(q.alpha_l, tmp);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, tmp);

  // output results
  FILE_LOG(LOG_CONTACT) << "results: " << std::endl;
  FILE_LOG(LOG_CONTACT) << "new Jc_v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "new Dc_v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "new Jl_v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "new Jx_v: " << q.Jx_v << std::endl;

  // see whether another QP must be solved
  if (q.Jc_v.size() > 0 && *min_element(q.Jc_v.begin(), q.Jc_v.end()) < -TOL)
  {
    FILE_LOG(LOG_CONTACT) << "minimum Jc*v: " << *min_element(q.Jc_v.begin(), q.Jc_v.end()) << std::endl;
    FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
    solve_qp_work(q, z);
    update_impulses(q, z);
  }
  else 
    if (q.Jl_v.size() > 0 && *min_element(q.Jl_v.begin(), q.Jl_v.end()) < -TOL)
    {
      FILE_LOG(LOG_CONTACT) << "minimum Jl*v: " << *min_element(q.Jl_v.begin(), q.Jl_v.end()) << std::endl;
      FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      update_impulses(q, z);
    }
  else
  {
    pair<Real*, Real*> mm = boost::minmax_element(q.Jx_v.begin(), q.Jx_v.end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_CONTACT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_CONTACT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
      solve_qp_work(q, z);
      update_impulses(q, z);
    }
  }

  // save normal contact impulses
  for (unsigned i=0; i< N_CONTACTS; i++)
    q.contact_events[i]->contact_impulse = q.contact_events[i]->contact_normal * q.alpha_c[i];

  // save tangent contact impulses
  for (unsigned i=0, j=0; i< N_CONTACTS; i++)
  {
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan1 * q.beta_c[j++];
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan2 * q.beta_c[j++];
  }

  // save limit impulses
  for (unsigned i=0; i< N_LIMITS; i++)
    q.limit_events[i]->limit_impulse = q.alpha_l[i]; 
}

/// Solves the nonlinearly constrained quadratic program (potentially solves two nQPs, actually)
void ImpactEventHandler::solve_nqp(EventProblemData& q, Real poisson_eps)
{
  SAFESTATIC VectorN z, tmp, tmp2;
  const Real TOL = poisson_eps;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_CONSTRAINT_DOF_IMP = q.N_CONSTRAINT_DOF_IMP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;
  const unsigned ALPHA_L_IDX = N_CONTACTS*3 + q.N_LIN_CONE*1;

  // solve the nonlinearly constrained QP
  solve_nqp_work(q, z);

  // apply (Poisson) restitution to contacts
  for (unsigned i=0; i< N_CONTACTS; i++)
    z[i] *= ((Real) 1.0 + q.contact_events[i]->contact_epsilon);

  // apply (Poisson) restitution to limits
  for (unsigned i=0; i< N_LIMITS; i++)
    z[ALPHA_L_IDX+i] *= ((Real) 1.0 + q.limit_events[i]->limit_epsilon);

  // save impulses in q
  update_impulses(q, z);

  // update Jc_v, Dc_v, Jl_v, and Jx_v
  q.Jc_v += q.Jc_iM_JcT.mult(q.alpha_c, tmp);
  q.Jc_v += q.Jc_iM_DcT.mult(q.beta_c, tmp);
  q.Jc_v += q.Jc_iM_JlT.mult(q.alpha_l, tmp);
  q.Jc_v += q.Jc_iM_JxT.mult(q.alpha_x, tmp);
  q.Dc_v += q.Jc_iM_DcT.transpose_mult(q.alpha_c, tmp);
  q.Dc_v += q.Dc_iM_DcT.mult(q.beta_c, tmp);
  q.Dc_v += q.Dc_iM_JlT.mult(q.alpha_l, tmp);
  q.Dc_v += q.Dc_iM_JxT.mult(q.alpha_x, tmp);
  q.Jl_v += q.Jc_iM_JlT.transpose_mult(q.alpha_c, tmp);
  q.Jl_v += q.Dc_iM_JlT.transpose_mult(q.beta_c, tmp);
  q.Jl_v += q.Jl_iM_JlT.mult(q.alpha_l, tmp);
  q.Jl_v += q.Jl_iM_JxT.mult(q.alpha_x, tmp);
  q.Jx_v += q.Jc_iM_JxT.transpose_mult(q.alpha_c, tmp);
  q.Jx_v += q.Dc_iM_JxT.transpose_mult(q.beta_c, tmp);
  q.Jx_v += q.Jl_iM_JxT.transpose_mult(q.alpha_l, tmp);
  q.Jx_v += q.Jx_iM_JxT.mult(q.alpha_x, tmp);

  // see whether another QP must be solved
  if (q.Jc_v.size() > 0 && *min_element(q.Jc_v.begin(), q.Jc_v.end()) < -TOL)
  {
    FILE_LOG(LOG_CONTACT) << "minimum Jc*v: " << *min_element(q.Jc_v.begin(), q.Jc_v.end()) << std::endl;
    FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
    solve_nqp_work(q, z);
    update_impulses(q, z);
  }
  else 
    if (q.Jl_v.size() > 0 && *min_element(q.Jl_v.begin(), q.Jl_v.end()) < -TOL)
    {
      FILE_LOG(LOG_CONTACT) << "minimum Jl*v: " << *min_element(q.Jl_v.begin(), q.Jl_v.end()) << std::endl;
      FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
      solve_nqp_work(q, z);
      update_impulses(q, z);
    }
  else
  {
    pair<Real*, Real*> mm = boost::minmax_element(q.Jx_v.begin(), q.Jx_v.end());
    if (q.Jx_v.size() > 0 && (*mm.first < -TOL || *mm.second > TOL))
    {
      FILE_LOG(LOG_CONTACT) << "minimum J*v: " << *mm.first << std::endl;
      FILE_LOG(LOG_CONTACT) << "maximum J*v: " << *mm.second << std::endl;
      FILE_LOG(LOG_CONTACT) << " -- running another QP iteration..." << std::endl;
      solve_nqp_work(q, z);
      update_impulses(q, z);
    }
  }

  // save normal contact impulses
  for (unsigned i=0; i< N_CONTACTS; i++)
    q.contact_events[i]->contact_impulse = q.contact_events[i]->contact_normal * q.alpha_c[i];

  // save tangent contact impulses
  for (unsigned i=0, j=0; i< N_CONTACTS; i++)
  {
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan1 * q.beta_c[j++];
    q.contact_events[i]->contact_impulse += q.contact_events[i]->contact_tan2 * q.beta_c[j++];
  }

  // save limit impulses
  for (unsigned i=0; i< N_LIMITS; i++)
    q.limit_events[i]->limit_impulse = q.alpha_l[i]; 
}

/// Updates impulses in q using concatenated vector of impulses z
void ImpactEventHandler::update_impulses(EventProblemData& q, const VectorN& z)
{
  SAFESTATIC VectorN tmp;

  // get the number of different types of each event
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = BETA_C_IDX + N_CONTACTS*2;
  const unsigned ALPHA_L_IDX = NBETA_C_IDX + q.N_LIN_CONE*2;
  const unsigned BETA_T_IDX = ALPHA_L_IDX + N_LIMITS;
  const unsigned ALPHA_X_IDX = BETA_T_IDX + q.N_CONSTRAINT_DOF_IMP;
  const unsigned BETA_X_IDX = ALPHA_X_IDX + N_CONSTRAINT_EQNS;
  const unsigned NVARS = BETA_X_IDX + q.N_CONSTRAINT_DOF_EXP;
  assert(NVARS == z.size());

  // add to impulses in q
  q.alpha_c += z.get_sub_vec(ALPHA_C_IDX, BETA_C_IDX, tmp);
  q.beta_c += z.get_sub_vec(BETA_C_IDX, NBETA_C_IDX, tmp);
  q.alpha_l += z.get_sub_vec(ALPHA_L_IDX, BETA_T_IDX, tmp);
  q.alpha_x += z.get_sub_vec(ALPHA_X_IDX, BETA_X_IDX, tmp);

  // update components of beta_c with negative impulses
  for (unsigned i=0, j=0, k=NBETA_C_IDX; i< N_CONTACTS; i++, j+= 2)
    if (q.contact_events[i]->contact_NK < UINF &&
        q.contact_events[i]->contact_NK >= 2)
    {
      q.beta_c[j] -= z[k++];
      q.beta_c[j+1] -= z[k++];
    }
}

/// Solves the nonlinearly constrained quadratic program (does all of the work)
void ImpactEventHandler::solve_nqp_work(EventProblemData& q, VectorN& z)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  SAFESTATIC MatrixN sub, t1, t2, t3, t4, A, MR, RTH;
  SAFESTATIC VectorN tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_CONSTRAINT_DOF_EXP = q.N_CONSTRAINT_DOF_EXP;
  const unsigned N_CONSTRAINT_DOF_IMP = q.N_CONSTRAINT_DOF_IMP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;
  const unsigned N_LIN_CONE = q.N_LIN_CONE;
  const unsigned N_TRUE_CONE = q.N_TRUE_CONE;
  const unsigned N_LOOPS = q.N_LOOPS;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*3;
  const unsigned ALPHA_L_IDX = N_LIN_CONE*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned BETA_T_IDX = N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;
  const unsigned BETA_X_IDX = N_CONSTRAINT_DOF_IMP + BETA_T_IDX;
  const unsigned DELTA_IDX = BETA_X_IDX + N_CONSTRAINT_DOF_EXP;
  const unsigned NVARS = N_LOOPS + DELTA_IDX; 

  // setup the optimization data
  SAFESTATIC ImpactOptData opt_data;
  MatrixN& R = opt_data.R;
  MatrixNN& H = opt_data.H;
  VectorN& c = opt_data.c; 

  // init z
  z.set_zero(NVARS);

  // first, compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_EXP > 0)
  {
    // compute the homogeneous solution
    A.copy_from(q.Jx_iM_JxT);
    tmpv.copy_from(q.Jx_v).negate();
    LinAlg::solve_LS_fast(A, tmpv);
    z.set_sub_vec(ALPHA_X_IDX, tmpv);

    // compute the nullspace
    A.resize(N_CONSTRAINT_EQNS_EXP, NVARS-N_LOOPS);
    MatrixN::transpose(q.Jc_iM_JxT, t1);
    MatrixN::transpose(q.Dc_iM_JxT, t2);
    MatrixN::transpose(q.Jl_iM_JxT, t3);
    A.set_sub_row_block(0, &t1, &t2, &t3, &q.Jx_iM_JxT);
    LinAlg::nullspace(A, t2);

    // augment nullspace with loops
    R.set_zero(t2.rows()+N_LOOPS, t2.columns()+N_LOOPS);
    R.set_sub_mat(0,0,t2);
    for (unsigned i=0, j=R.columns(), k=R.rows(); i< N_LOOPS; i++, j++, k++)
      R(j,k) = (Real) 1.0;
  }
  else
  {
    // setup the nullspace
    R.set_zero(NVARS,NVARS);
    for (unsigned i=0; i< NVARS; i++) 
      R(i,i) = (Real) 1.0;
  }

  // get number of qp variables
  const unsigned N_PRIMAL = R.columns();

  // setup number of nonlinear inequality constraints
  const unsigned N_JF_DOF = N_CONSTRAINT_DOF_IMP + N_CONSTRAINT_DOF_EXP;
  const unsigned NONLIN_INEQUAL = N_TRUE_CONE + N_JF_DOF + N_LOOPS;

  // setup the optimization data
  set_optimization_data(q, opt_data);
  opt_data.z.copy_from(z);
  opt_data.epd = &q;

  // setup the optimization parameters
  SAFESTATIC OptParams oparams(N_PRIMAL, NONLIN_INEQUAL, 0, &sqp_f0, &sqp_fx, NULL, &sqp_grad0, &sqp_cJac, NULL, &sqp_hess);
  oparams.n = N_PRIMAL;
  oparams.m = NONLIN_INEQUAL;
  oparams.eps_feas = NEAR_ZERO;
  oparams.A.resize(0, N_PRIMAL);
  oparams.b.resize(0);
  oparams.data = (void*) &opt_data;
  MatrixN& M = oparams.M;
  VectorN& qq = oparams.q;

  // init the QP matrix and vector
  H.resize(N_PRIMAL);
  c.resize(H.size());

  // setup quadratic matrix 
  unsigned row = 0;

  // row (block) 1 -- Jc * iM * [Jc' Dc' Jl' Dt' Jx' Dx']
  H.set_sub_row_block(0,&q.Jc_iM_JcT, &q.Jc_iM_DcT, &q.Jc_iM_JlT, 
                        &q.Jc_iM_DtT, &q.Jc_iM_JxT, &q.Jc_iM_DxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Jc' Dc' Jl' Dt' Jx' Dx']
  MatrixN::transpose(q.Jc_iM_DcT, t1);
  H.set_sub_row_block(row, &t1,          &q.Dc_iM_DcT, &q.Dc_iM_JlT, 
                           &q.Dc_iM_DtT, &q.Dc_iM_JxT, &q.Dc_iM_DxT);
  row += N_CONTACTS*2;

  // row (block 3) -- Jl * iM * [Jc' Dc' Jl' Dt' Jx' Dx']
  MatrixN::transpose(q.Jc_iM_JlT, t1);
  MatrixN::transpose(q.Dc_iM_JlT, t2);
  H.set_sub_row_block(row, &t1,          &t2,          &q.Jl_iM_JlT, 
                           &q.Jl_iM_DtT, &q.Jl_iM_JxT, &q.Jl_iM_DxT);
  row += N_LIMITS;
  
  // row (block 4) -- Jx * iM * [Jc' Dc' Jl Dt' Jx' Dx']
  MatrixN::transpose(q.Jc_iM_JxT, t1);
  MatrixN::transpose(q.Dc_iM_JxT, t2);
  MatrixN::transpose(q.Jl_iM_JxT, t3);
  MatrixN::transpose(q.Dt_iM_JxT, t4);
  H.set_sub_row_block(row, &t1,          &t2,          &t3, 
                           &t4, &q.Jx_iM_JxT, &q.Jx_iM_DxT);

  // setup linear vector 
  c.set_sub_vec(ALPHA_C_IDX, q.Jc_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  c.set_sub_vec(ALPHA_L_IDX, q.Jl_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // ****** now setup linear inequality constraints ******
  // (ordered) components of z and corresponding length:
  // alpha_c / N_CONTACTS
  // beta_c  / N_LIN_CONE*2 + N_TRUE_CONE*2
  // nbeta_c / N_LIN_CONE*2
  // alpha_l / N_LIMITS

  // determine whether to use kappa constraint
  const unsigned KAPPA = (q.use_kappa) ? 1 : 0;

  // determine number of linear inequality constraints
  const unsigned N_INEQUAL = N_CONTACTS*2 + N_LIN_CONE*4 + N_K_TOTAL + N_LIMITS*2 + KAPPA;

  // resize matrices and vectors
  M.set_zero(N_INEQUAL, N_PRIMAL);
  qq.set_zero(N_INEQUAL);

  // setup the alpha_c >= 0 constraint
  row = 0; 
  for (unsigned i=0; i< N_CONTACTS; i++)
    M(row++, ALPHA_C_IDX+i) = (Real) 1.0;

  // setup the alpha_l >= 0 constraint
  for (unsigned i=0; i< N_LIMITS; i++)
    M(row++, ALPHA_L_IDX+i) = (Real) 1.0;

  // setup the beta_c >= 0 constraints (where necessary)
  for (unsigned i=0, r=0, s=0; i< N_CONTACTS; i++, r+= 2)
  {
    if (q.contact_events[i]->contact_NK == UINF ||
        q.contact_events[i]->contact_NK < 2)
      continue;

    // setup beta_c >= 0
    M(row++, BETA_C_IDX+r) = (Real) 1.0;
    M(row++, BETA_C_IDX+r+1) = (Real) 1.0;
    M(row++, NBETA_C_IDX+s) = (Real) 1.0;
    M(row++, NBETA_C_IDX+s+1) = (Real) 1.0;

    // update s
    s += 2;
  }

  // setup the Jc*v+ >= 0 constraint
  // Jc*(inv(M)*impulses + v) >= 0, Jc*inv(M)*impulses >= -Jc*v
  qq.set_sub_vec(row, q.Jc_v);
  H.get_sub_mat(ALPHA_C_IDX, ALPHA_C_IDX+N_CONTACTS, 0, H.columns(), sub);
  M.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the Jl*v+ >= 0 constraint
  qq.set_sub_vec(row, q.Jl_v);
  H.get_sub_mat(ALPHA_L_IDX, ALPHA_L_IDX+N_LIMITS, 0, H.columns(), sub);
  M.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the linearized contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, r=0, s=0; i< N_CONTACTS; i++, r+= 2)
  {
    // check whether we are using a true friction cone or not
    if (q.contact_events[i]->contact_NK == std::numeric_limits<unsigned>::max())
      continue;

    // determine contact velocity for viscous friction
    Real vel = std::sqrt(sqr(q.Dc_v[0]) + sqr(q.Dc_v[1]));

    // setup the Coulomb friction constraints 
    for (unsigned j=0; j< q.contact_events[i]->contact_NK; j++)
    {
      M(row, ALPHA_C_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      Real theta = (Real) j/(q.contact_events[i]->contact_NK-1) * M_PI_2;
      const Real ct = std::cos(theta);
      const Real st = std::sin(theta);
      M(row, BETA_C_IDX+r) = -ct;
      M(row, NBETA_C_IDX+s) = -ct;
      M(row, BETA_C_IDX+r+1) = -st;
      M(row, NBETA_C_IDX+s+1) = -st;

      // setup the viscous friction component
      qq[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }

    // update s
    s+= 2;
  }

  // setup the normal impulse constraint
  // 1'cn <= kappa (equiv. to -1'cn >= -kappa)
  if (q.use_kappa)
  {
    for (unsigned i=0; i< N_CONTACTS; i++)
      M(row, ALPHA_C_IDX+i) = (Real) -1.0;
    qq[row] = q.kappa;
  }

  // negate qq (all constraints above were for -qq)
  qq.negate();

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H);
  R.transpose_mult(c, tmpv);
  c.copy_from(tmpv);
  RTH.mult(z, tmpv);
  c += tmpv;

  // setup constraints M*x >= qq in nullspace, yielding
  // M(R*y + z) >= qq, yielding M*R*y >= qq - M*z
  M.mult(z, tmpv);
  qq -= tmpv;
  M.mult(R, MR);
  M.copy_from(MR);

  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::solve_nqp_work() entered" << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_CONTACT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_CONTACT) << "M matrix: " << std::endl << M;
  FILE_LOG(LOG_CONTACT) << "q vector: " << qq << std::endl;

  // setup the maximum numbers of iterations 
  const unsigned MAX_SQP_ITER = 10 + oparams.n + oparams.m + oparams.q.size(); 
  oparams.max_iterations = MAX_SQP_ITER;

  // solve the nonlinear QP using sequential QP algorithm
  y.set_zero(N_PRIMAL);
y[0] = .004905;
y[3] = y[0];
oparams.max_iterations = 10000;
  for (unsigned i=N_PRIMAL-N_LOOPS; i< N_PRIMAL; i++)
    y[i] = (Real) 0.5;
  Optimization::sqp(oparams, y);
//oparams.max_iterations = 1000;
//  Optimization::optimize_convex_pd(oparams, y);

  // compute the particular solution
  z += R.mult(y, tmpv);

  FILE_LOG(LOG_CONTACT) << "nonlinear QP solution: " << z << std::endl; 
  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::solve_nqp() exited" << std::endl;
}

/// Solves the quadratic program (does all of the work)
/**
 * \note this is the version without joint friction forces
 */
void ImpactEventHandler::solve_qp_work(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixN sub, t1, t2, t3, neg1, A, AR, R, RTH;
  SAFESTATIC MatrixNN H, MM;
  SAFESTATIC VectorN negv, c, qq, nb, tmpv, y;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*2 + BETA_C_IDX;
  const unsigned ALPHA_L_IDX = N_CONTACTS*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned NVARS = N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;

  // first, solve for impulses that satisfy explicit constraint equations
  // and compute the appropriate nullspace 
  if (N_CONSTRAINT_EQNS_EXP > 0)
  {
    // compute the homogeneous solution
    A.copy_from(q.Jx_iM_JxT);
    z.copy_from(q.Jx_v).negate();
    LinAlg::solve_LS_fast(A, z);

    // compute the nullspace
    A.resize(N_CONSTRAINT_EQNS_EXP, NVARS);
    MatrixN::transpose(q.Jc_iM_JxT, t1);
    MatrixN::transpose(q.Dc_iM_JxT, t2);
    MatrixN::transpose(q.Jl_iM_JxT, t3);
    neg1.copy_from(t2).negate();
    A.set_sub_row_block(0, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);
    LinAlg::nullspace(A, R);
  }
  else
  {
    R.set_zero(NVARS,NVARS);
    for (unsigned i=0; i< NVARS; i++) 
      R(i,i) = (Real) 1.0;
    z.set_zero(NVARS);
  }

  // get number of qp variables
  const unsigned N_PRIMAL = R.columns();

  // init the QP matrix and vector
  const unsigned KAPPA = (q.use_kappa) ? 1 : 0;
  const unsigned N_INEQUAL = N_CONTACTS + N_K_TOTAL + N_LIMITS + KAPPA;
  H.resize(N_PRIMAL);
  c.resize(H.size());
  A.set_zero(N_INEQUAL, N_PRIMAL);
  nb.set_zero(N_INEQUAL);
  MM.set_zero(N_PRIMAL + N_INEQUAL);
  qq.resize(MM.size());

  // setup [Q M'; -M 0]
  unsigned col = 0, row = 0;

  // row (block) 1 -- Jc * iM * [Jc' Dc' -Dc' Jl' Jx']
  neg1.copy_from(q.Jc_iM_DcT).negate();
  H.set_sub_row_block(0, &q.Jc_iM_JcT, &q.Jc_iM_DcT, &neg1, &q.Jc_iM_JlT, 
                      &q.Jc_iM_JxT);
  row += N_CONTACTS;
  
  // row (block) 2 -- Dc * iM * [Jc' Dc' -Dc' Jl' Jx']
  MatrixN::transpose(q.Jc_iM_DcT, t1);
  neg1.copy_from(q.Dc_iM_DcT).negate();
  H.set_sub_row_block(row, &t1, &q.Dc_iM_DcT, &neg1, &q.Dc_iM_JlT, 
                      &q.Dc_iM_JxT);

  // row (block 3) -- negated block 2
  H.get_sub_mat(row, row+N_CONTACTS*2, 0, H.columns(), sub);
  H.set_sub_mat(row+N_CONTACTS*2, 0, sub.negate());
  row += N_CONTACTS*4;

  // row (block 4) -- Jl * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JlT, t1);
  MatrixN::transpose(q.Dc_iM_JlT, t2);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &q.Jl_iM_JlT, &q.Jl_iM_JxT);
  row += N_LIMITS;
  
  // row (block 5) -- Jx * iM * [Jc' Dc' -Dc' Jl Jx']
  MatrixN::transpose(q.Jc_iM_JxT, t1);
  MatrixN::transpose(q.Dc_iM_JxT, t2);
  MatrixN::transpose(q.Jl_iM_JxT, t3);
  neg1.copy_from(t2).negate();
  H.set_sub_row_block(row, &t1, &t2, &neg1, &t3, &q.Jx_iM_JxT);

  // setup c 
  c.set_sub_vec(ALPHA_C_IDX, q.Jc_v);         
  c.set_sub_vec(BETA_C_IDX, q.Dc_v);         
  negv.copy_from(q.Dc_v).negate();
  c.set_sub_vec(NBETA_C_IDX, negv);           
  c.set_sub_vec(ALPHA_L_IDX, q.Jl_v);         
  c.set_sub_vec(ALPHA_X_IDX, q.Jx_v);         

  // setup the Jc*v+ >= 0 constraint
  // Jc*(inv(M)*impulses + v) >= 0, Jc*inv(M)*impulses >= -Jc*v
  row = 0; col = 0;
  nb.set_sub_vec(row, q.Jc_v);
  H.get_sub_mat(ALPHA_C_IDX, ALPHA_C_IDX+N_CONTACTS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_CONTACTS;
  
  // setup the Jl*v+ >= 0 constraint
  nb.set_sub_vec(row, q.Jl_v);
  H.get_sub_mat(ALPHA_L_IDX, ALPHA_L_IDX+N_LIMITS, 0, H.columns(), sub);
  A.set_sub_mat(row, 0, sub);
  row += N_LIMITS;

  // setup the contact friction constraints
  // mu_c*cn + mu_v*cvel >= beta
  for (unsigned i=0, k=0; i< N_CONTACTS; i++, k+= 2)
  {
    // initialize the contact velocity
    Real vel = std::sqrt(sqr(q.Dc_v[k]) + sqr(q.Dc_v[k+1]));

    // setup the Coulomb friction inequality constraints for this contact
    for (unsigned j=0; j< q.contact_events[i]->contact_NK; j++)
    {
      Real theta = (Real) j/(q.contact_events[i]->contact_NK-1) * M_PI_2;
      const Real ct = std::cos(theta);
      const Real st = std::sin(theta);
      A(row, ALPHA_C_IDX+i) = q.contact_events[i]->contact_mu_coulomb;
      A(row, BETA_C_IDX+k) = -ct;
      A(row, NBETA_C_IDX+k) = -ct;
      A(row, BETA_C_IDX+k+1) = -st;
      A(row, NBETA_C_IDX+k+1) = -st;

      // setup the viscous friction component
      nb[row++] = q.contact_events[i]->contact_mu_viscous * vel;
    }
  }

  // setup the normal impulse constraint
  // 1'cn <= kappa (equiv. to -1'cn >= -kappa)
  if (q.use_kappa)
  {
    for (unsigned i=0; i< N_CONTACTS; i++)
      A(row, ALPHA_C_IDX+i) = (Real) -1.0;
    nb[row] = q.kappa;
  }

  // setup optimizations in nullspace 
  R.transpose_mult(H, RTH);
  RTH.mult(R, H);
  R.transpose_mult(c, tmpv);
  c.copy_from(tmpv);
  RTH.mult(z, tmpv);
  c += tmpv;

  // setup constraints A*x >= b in nullspace, yielding
  // A(R*y + z) >= b, yielding R*y >= b - A*z
  A.mult(z, tmpv);
  nb += tmpv;
  A.mult(R, AR);

  // setup the LCP matrix
  MM.set_sub_mat(0, 0, H);
  MM.set_sub_mat(N_PRIMAL, 0, AR);
  MM.set_sub_mat(0, N_PRIMAL, AR.negate(), true);

  // setup the LCP vector
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N_PRIMAL, nb);

  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::solve_qp() entered" << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Dc': " << std::endl << q.Jc_iM_DcT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jl': " << std::endl << q.Jc_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jx': " << std::endl << q.Jc_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Dc': " << std::endl << q.Dc_iM_DcT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Jl': " << std::endl << q.Dc_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Dc * inv(M) * Jx': " << std::endl << q.Dc_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jl * inv(M) * Jl': " << std::endl << q.Jl_iM_JlT;
  FILE_LOG(LOG_CONTACT) << "  Jl * inv(M) * Jx': " << std::endl << q.Jl_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jx * inv(M) * Jx': " << std::endl << q.Jx_iM_JxT;
  FILE_LOG(LOG_CONTACT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Dc * v: " << q.Dc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jx * v: " << q.Jx_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "H matrix: " << std::endl << H;
  FILE_LOG(LOG_CONTACT) << "c vector: " << c << std::endl;
  FILE_LOG(LOG_CONTACT) << "A matrix: " << std::endl << A;
  FILE_LOG(LOG_CONTACT) << "b vector: " << (-nb) << std::endl;
  FILE_LOG(LOG_CONTACT) << "LCP matrix: " << std::endl << MM; 
  FILE_LOG(LOG_CONTACT) << "LCP vector: " << qq << std::endl; 

  // solve the LCP using Lemke's algorithm
  if (!Optimization::lcp_lemke_regularized(MM, qq, tmpv))
    throw std::runtime_error("Unable to solve event QP!");

  // get the nullspace solution out
  FILE_LOG(LOG_CONTACT) << "LCP solution: " << tmpv << std::endl; 
  tmpv.get_sub_vec(0, N_PRIMAL, y);
  R.mult(y, tmpv);
  z += tmpv;

  FILE_LOG(LOG_CONTACT) << "QP solution: " << z << std::endl; 
  FILE_LOG(LOG_CONTACT) << "ImpactEventHandler::solve_qp() exited" << std::endl;
}

/// Solves the (frictionless) LCP
void ImpactEventHandler::solve_lcp(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixNN UL, LR, MM;
  SAFESTATIC MatrixN UR, t2, iJx_iM_JxT;
  SAFESTATIC VectorN alpha_c, alpha_l, alpha_x, v1, v2, qq;

  // get the number of different types of each event
  const unsigned N_CONTACTS = q.N_CONTACTS;
  const unsigned N_LIMITS = q.N_LIMITS;
  const unsigned N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
  const unsigned N_K_TOTAL = q.N_K_TOTAL;
  const unsigned N_LIN_CONE = q.N_LIN_CONE;
  const unsigned N_TRUE_CONE = q.N_TRUE_CONE;
  const unsigned N_LOOPS = q.N_LOOPS;

  // setup variable indices
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = N_CONTACTS;
  const unsigned NBETA_C_IDX = N_CONTACTS*3;
  const unsigned ALPHA_L_IDX = N_LIN_CONE*2 + NBETA_C_IDX;
  const unsigned ALPHA_X_IDX = N_LIMITS + ALPHA_L_IDX;
  const unsigned BETA_T_IDX = q.N_CONSTRAINT_EQNS_EXP + ALPHA_X_IDX;
  const unsigned BETA_X_IDX = q.N_CONSTRAINT_DOF_IMP + BETA_T_IDX;
  const unsigned DELTA_IDX = BETA_X_IDX + q.N_CONSTRAINT_DOF_EXP;
  const unsigned NVARS = N_LOOPS + DELTA_IDX; 

  // setup sizes
  UL.resize(N_CONTACTS);
  UR.resize(N_CONTACTS, N_LIMITS);
  LR.resize(N_LIMITS);

  // setup primary terms -- first upper left hand block of matrix
  LinAlg::pseudo_inverse(q.Jx_iM_JxT, iJx_iM_JxT);
  q.Jc_iM_JxT.mult(iJx_iM_JxT, t2);
  t2.mult_transpose(q.Jc_iM_JxT, UL);
  // now do upper right hand block of matrix
  t2.mult_transpose(q.Jl_iM_JxT, UR);
  // now lower right hand block of matrix
  q.Jl_iM_JxT.mult(iJx_iM_JxT, t2);
  t2.mult_transpose(q.Jl_iM_JxT, LR);

  // subtract secondary terms
  UL -= q.Jc_iM_JcT;
  UR -= q.Jc_iM_JlT;
  LR -= q.Jl_iM_JlT;

  // now negate all terms
  UL.negate();
  UR.negate();
  LR.negate();

  // setup the LCP matrix
  MM.resize(N_CONTACTS + N_LIMITS);
  MM.set_sub_mat(0, 0, UL);
  MM.set_sub_mat(0, N_CONTACTS, UR);
  MM.set_sub_mat(N_CONTACTS, 0, UR, true);
  MM.set_sub_mat(N_CONTACTS, N_CONTACTS, LR);

  // setup the LCP vector
  qq.resize(MM.size());
  iJx_iM_JxT.mult(q.Jx_v, v2); 
  q.Jc_iM_JxT.mult(v2, v1);
  v1 -= q.Jc_v;
  qq.set_sub_vec(0, v1);
  q.Jl_iM_JxT.mult(v2, v1);
  v1 -= q.Jl_v;
  qq.set_sub_vec(N_CONTACTS, v1);
  qq.negate();

  FILE_LOG(LOG_CONTACT) << "ImpulseEventHandler::solve_lcp() entered" << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_CONTACT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_CONTACT) << "  LCP matrix: " << std::endl << MM;
  FILE_LOG(LOG_CONTACT) << "  LCP vector: " << qq << std::endl;

  // solve the LCP
  if (!Optimization::lcp_lemke_regularized(MM, qq, z))
    throw std::runtime_error("Unable to solve event LCP!");

  // determine the value of kappa
  q.kappa = (Real) 0.0;
  for (unsigned i=0; i< N_CONTACTS; i++)
    q.kappa += z[i];

  // get alpha_c and alpha_l
  z.get_sub_vec(0, N_CONTACTS, alpha_c);
  z.get_sub_vec(N_CONTACTS, z.size(), alpha_l);

  // Mv^* - Mv = Jc'*alpha_c + Jl'*alpha_l + Jx'*alpha_x

  // Mv^* - Mv^- = Jx'*alpha_x
  // Jx*v^*     = 0
  // v^* = v^- + inv(M)*Jx'*alpha_x
  // Jx*v^- + Jx*inv(M)*Jx'*alpha_x = 0

  // Jx*inv(M)*Jx'*alpha_x = -Jx*(v + inv(M)*Jc'*alpha_c + inv(M)*Jl'*alpha_l)

  // compute alpha_x 
  q.Jc_iM_JxT.transpose_mult(alpha_c, v1);
  q.Jl_iM_JxT.transpose_mult(alpha_l, v2); 
  v1 += v2;
  v1 += q.Jx_v;
  v1.negate();
  iJx_iM_JxT.mult(v1, alpha_x);

  // setup the homogeneous solution
  z.set_zero(NVARS);
  z.set_sub_vec(ALPHA_C_IDX, alpha_c);
  z.set_sub_vec(ALPHA_L_IDX, alpha_l);
  z.set_sub_vec(ALPHA_X_IDX, alpha_x);

  FILE_LOG(LOG_CONTACT) << "  LCP result: " << z << std::endl;
  FILE_LOG(LOG_CONTACT) << "  kappa: " << q.kappa << std::endl;
  FILE_LOG(LOG_CONTACT) << "ImpulseEventHandler::solve_lcp() exited" << std::endl;
}

/// Gets the super body (articulated if any)
DynamicBodyPtr ImpactEventHandler::get_super_body(SingleBodyPtr sb)
{
  ArticulatedBodyPtr ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return sb;
}

/// The Hessian
void ImpactEventHandler::sqp_hess(const VectorN& x, Real objscal, const VectorN& hlambda, const VectorN& nu, MatrixNN& H, void* data)
{
  SAFESTATIC VectorN Gx, w, tmpv, ff, lambda;
  SAFESTATIC MatrixN t1, t2;
  const Real INFEAS_TOL = 1e-8;

  // get the optimization data
  const ImpactOptData& opt_data = *(const ImpactOptData*) data;

  // get the event problem data
  const EventProblemData& epd = *opt_data.epd;

  // setup constants
  const unsigned N_CONTACTS = epd.N_CONTACTS;
  const unsigned N_LIMITS = epd.N_LIMITS;
  const unsigned N_TRUE_CONE = epd.N_TRUE_CONE;
  const unsigned N_LIN_CONE = epd.N_LIN_CONE;
  const unsigned N_LOOPS = epd.N_LOOPS;
  const unsigned N_CONSTRAINT_EQNS_EXP = epd.N_CONSTRAINT_EQNS_EXP; 
  const unsigned N_CONSTRAINT_DOF_IMP = epd.N_CONSTRAINT_DOF_IMP;
  const unsigned N_CONSTRAINT_DOF_EXP = epd.N_CONSTRAINT_DOF_EXP;
  const unsigned N_JOINT_DOF = N_CONSTRAINT_DOF_EXP + N_CONSTRAINT_DOF_IMP;
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = ALPHA_C_IDX + N_CONTACTS;
  const unsigned NBETA_C_IDX = BETA_C_IDX + N_CONTACTS*2;
  const unsigned ALPHA_L_IDX = NBETA_C_IDX + N_LIN_CONE*2;
  const unsigned ALPHA_X_IDX = ALPHA_L_IDX + N_LIMITS;
  const unsigned BETA_T_IDX = ALPHA_X_IDX + N_CONSTRAINT_EQNS_EXP;
  const unsigned BETA_X_IDX = BETA_T_IDX + N_CONSTRAINT_DOF_IMP;
  const unsigned DELTA_IDX = BETA_X_IDX + N_CONSTRAINT_DOF_EXP;

  // get necessary data
  const VectorN& z = opt_data.z;
  const MatrixN& R = opt_data.R;
  const MatrixNN& G = opt_data.H;
  const VectorN& c = opt_data.c;

  // objective function is quadratic
  H.copy_from(G) *= objscal;

  // setup the constraint index
  unsigned index = 0;

  // add in constraints for true friction cone
  for (unsigned i=0; i< N_TRUE_CONE; i++)
  {
    // get the contact event index
    const unsigned CIDX = opt_data.cone_contacts[i];
    const unsigned BETA_CX = BETA_C_IDX + CIDX*2;
    const unsigned BETA_CY = BETA_CX + 1;
    const unsigned ALPHA_C = ALPHA_C_IDX + CIDX;

    // compute contact friction
    R.get_row(BETA_CX, tmpv);
    (*VectorN::outer_prod(tmpv, tmpv, &t2));
    R.get_row(BETA_CY, tmpv);
    (*VectorN::outer_prod(tmpv, tmpv, &t1));
    t2 += t1;
    t2 *= (Real) 2.0;
    R.get_row(ALPHA_C, tmpv);
    (*VectorN::outer_prod(tmpv, tmpv, &t1)) *= (Real) (2.0 * opt_data.c_mu_c[CIDX]);
    t2 -= t1;

    // add in t2 to hessian
    t2 *= hlambda[index++];
    H += t2;
  }

/*
    // setup numerical Hessian
    unsigned n = x.size();
    const Real h = NEAR_ZERO;
    const Real INV_H2 = 1.0 / (h*2);
    VectorN xx = x;
    VectorN v1, v2;
    H.resize(n);
    for (unsigned i=0; i< n; i++)
    {
      xx[i] += h;
      sqp_grad(xx, m, v1, data);
      xx[i] -= 2*h;
      sqp_grad(xx, m, v2, data);
      xx[i] += h;
      v1 -= v2;
      v1 *= INV_H2;
      H.set_column(i, v1);
    }
    // average values of the Hessian
    for (unsigned i=0; i< n; i++)
      for (unsigned j=i+1; j< n; j++)
        H(i,j) = H(j,i) = 0.5*(H(i,j) + H(j,i));
    return true;
*/
}

/// The gradient
void ImpactEventHandler::sqp_cJac(const VectorN& x, MatrixN& J, void* data)
{
  SAFESTATIC VectorN Gx, tmpv, tmpv2, ff, fff, w, wdelta, walphac, wbetac, walphal;
  SAFESTATIC VectorN wbetat, wbetax, zalphac, zbetac, zalphal, zbetat, zbetax;
  SAFESTATIC VectorN Zf, Zdf, Z1df, grad;
  SAFESTATIC MatrixN dX, Rd, tmpM, tmpM2, JcTRalphac, DcTRbetac, JlTRalphal;
  SAFESTATIC MatrixN Rbetat, DxTRbetax;
  const Real INFEAS_TOL = 1e-8;

  // get the optimization data
  const ImpactOptData& opt_data = *(const ImpactOptData*) data;

  // get the event problem data
  const EventProblemData& epd = *opt_data.epd;

  // setup constants
  const unsigned N_CONTACTS = epd.N_CONTACTS;
  const unsigned N_LIMITS = epd.N_LIMITS;
  const unsigned N_TRUE_CONE = epd.N_TRUE_CONE;
  const unsigned N_LIN_CONE = epd.N_LIN_CONE;
  const unsigned N_LOOPS = epd.N_LOOPS;
  const unsigned N_CONSTRAINT_EQNS_EXP = epd.N_CONSTRAINT_EQNS_EXP; 
  const unsigned N_CONSTRAINT_DOF_IMP = epd.N_CONSTRAINT_DOF_IMP;
  const unsigned N_CONSTRAINT_DOF_EXP = epd.N_CONSTRAINT_DOF_EXP;
  const unsigned N_JOINT_DOF = N_CONSTRAINT_DOF_EXP + N_CONSTRAINT_DOF_IMP;
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = ALPHA_C_IDX + N_CONTACTS;
  const unsigned NBETA_C_IDX = BETA_C_IDX + N_CONTACTS*2;
  const unsigned ALPHA_L_IDX = NBETA_C_IDX + N_LIN_CONE*2;
  const unsigned ALPHA_X_IDX = ALPHA_L_IDX + N_LIMITS;
  const unsigned BETA_T_IDX = ALPHA_X_IDX + N_CONSTRAINT_EQNS_EXP;
  const unsigned BETA_X_IDX = BETA_T_IDX + N_CONSTRAINT_DOF_IMP;
  const unsigned DELTA_IDX = BETA_X_IDX + N_CONSTRAINT_DOF_EXP;

  // get necessary data
  const VectorN& z = opt_data.z;
  const MatrixN& R = opt_data.R;
  const MatrixNN& G = opt_data.H;
  const VectorN& c = opt_data.c;

  // resize J
  J.resize(N_TRUE_CONE + N_LOOPS*2 + N_JOINT_DOF, x.size());

  // setup constraint index
  unsigned index = 0;

  // compute w
  R.mult(x, w) += z;

  // compute gradients for true friction cone
  for (unsigned i=0; i< N_TRUE_CONE; i++)
  {
    // get the contact event index
    const unsigned CIDX = opt_data.cone_contacts[i];
    const unsigned BETA_CX = BETA_C_IDX + CIDX*2;
    const unsigned BETA_CY = BETA_CX + 1;
    const unsigned ALPHA_C = ALPHA_C_IDX + CIDX;

    // compute contact friction
    R.get_row(BETA_CX, grad) *= ((Real) 2.0 * w[BETA_CX]);
    R.get_row(BETA_CY, tmpv) *= ((Real) 2.0 * w[BETA_CY]);
    grad += tmpv;
    R.get_row(ALPHA_C, tmpv) *= ((Real) 2.0 * w[ALPHA_C] * opt_data.c_mu_c[CIDX]);
    grad -= tmpv;
    J.set_row(index++, grad); 
  }

  // compute gradients for delta >= 0 constraints
  for (unsigned i=0; i< N_LOOPS; i++) 
  {
    const unsigned DIDX = DELTA_IDX + i;
    R.get_row(DIDX, grad);
    grad.negate();
    J.set_row(index++, grad);
  }

  // compute gradients for delta <= 1 constraints
  for (unsigned i=0; i< N_LOOPS; i++) 
  {
    // this is delta <= 1
    const unsigned DIDX = DELTA_IDX + i; 
    R.get_row(DIDX, grad);
    J.set_row(index++, grad);
  }

  // compute gradients for joint friction constraints
  for (unsigned i=0; i< N_JOINT_DOF; i++)
  {
    // original equation is mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||ff||
    //                   or mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||beta_x||

    // get the body index and body
    const unsigned AIDX = opt_data.body_indices[i];
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(epd.super_bodies[AIDX]);    

    // get indices for this body
    const vector<int>& alpha_c_indices = opt_data.alpha_c_indices[AIDX];
    const vector<int>& beta_nbeta_c_indices = opt_data.beta_nbeta_c_indices[AIDX];
    const vector<unsigned>& alpha_l_indices = opt_data.alpha_l_indices[AIDX];
    const vector<unsigned>& beta_t_indices = opt_data.beta_t_indices[AIDX];
    const vector<unsigned>& beta_x_indices = opt_data.beta_x_indices[AIDX];

    // get the relative joint friction index for this articulated body
    const unsigned RIDX = i - opt_data.joint_friction_start[AIDX];

    // get the index of the joint in the articulated body
    const unsigned JIDX = opt_data.true_indices[AIDX][RIDX];

    // get the number of implicit DOFs for this body, so we can determine
    // whether the relative index is implicit or explicit
    const unsigned N_ABODY_IMPLICIT = abody->num_joint_dof_implicit();

    // get the joint frictional force in w
    const unsigned FIDX = (RIDX < N_ABODY_IMPLICIT) ? BETA_T_IDX + opt_data.implicit_start[AIDX] : BETA_X_IDX + opt_data.explicit_start[AIDX];

    // get the three Z's
    const MatrixN& Zd = opt_data.Zd[AIDX][JIDX];
    const MatrixN& Z1d = opt_data.Z1d[AIDX][JIDX];
    const MatrixN& Z = opt_data.Z[AIDX][JIDX];

    // get squared Coulomb joint friction coefficient
    const Real MUCSQ = opt_data.j_mu_c[i];

    // get components of R corresponding to alpha_c, beta_c, nbeta_c,
    // alpha_l, beta_t, beta_x, and delta for this body and
    // multiply components of R by corresponding Jacobian
    contact_select(alpha_c_indices, beta_nbeta_c_indices, R, tmpM, tmpM2);
    abody->transpose_Jc_mult(tmpM, JcTRalphac);
    abody->transpose_Dc_mult(tmpM2, DcTRbetac);
    R.select_rows(alpha_l_indices.begin(), alpha_l_indices.end(), tmpM);
    abody->transpose_Jl_mult(tmpM, JlTRalphal);
    R.select_rows(beta_x_indices.begin(), beta_x_indices.end(), tmpM);
    abody->transpose_Dx_mult(tmpM, DxTRbetax);

    // get components of z and multiply by the proper Jacobians
    contact_select(alpha_c_indices, beta_nbeta_c_indices, z, tmpv, tmpv2);
    abody->transpose_Jc_mult(tmpv, zalphac);
    abody->transpose_Dc_mult(tmpv2, zbetac);
    z.select(alpha_l_indices.begin(), alpha_l_indices.end(), tmpv);
    abody->transpose_Jl_mult(tmpv, zalphal);
    z.select(beta_t_indices.begin(), beta_t_indices.end(), zbetat);
    z.select(beta_x_indices.begin(), beta_x_indices.end(), tmpv);
    abody->transpose_Dx_mult(tmpv, zbetax);

    // setup impulse 
    JcTRalphac.mult(x, walphac) += zalphac;
    DcTRbetac.mult(x, wbetac) += zbetac;      
    Rbetat.mult(x, wbetat) += zbetat;
    DxTRbetax.mult(x, wbetax) += zbetax;
    fff.copy_from(walphac);
    fff += wbetac;
    fff += walphal;
    fff += wbetat;
    fff += wbetax;

    // compute gradient of impulse 
    dX.copy_from(JcTRalphac) += DcTRbetac;
    dX += JlTRalphal;
    dX += Rbetat;
    dX += DxTRbetax;

    // get loop indices for this articulated body
    const vector<unsigned>& loop_indices = opt_data.loop_indices[AIDX];

    // two cases: joint is part of a loop or not
    if (loop_indices[JIDX] != std::numeric_limits<unsigned>::max())
    {
      // setup delta stuff
      const unsigned DELTA_START = DELTA_IDX + opt_data.delta_start[AIDX];
      const unsigned LOOP_IDX = DELTA_START+loop_indices[JIDX];
      R.get_sub_mat(LOOP_IDX, LOOP_IDX+1, 0, R.columns(), Rd);
      const Real ZDELTA = opt_data.z[LOOP_IDX];
      Rd.mult(x, wdelta);
      const Real DELTA = wdelta[0] + ZDELTA;

      // compute first component of gradient (dX' * Z' * Z * f)
      Z.mult(fff, Zf);
      Z.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, grad);

      // compute second component of gradient (Rd' * f * Z' * Zd * f)
      Zd.mult(fff, Zdf);
      Z.transpose_mult(Zdf, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv));
      grad += tmpv2;    

      // compute third component of gradient (d * dX' * Z' * Zd * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= DELTA;
      grad += tmpv2;

      // compute fourth component of gradient (d * dX' * Zd' * Z * f)
      Z.mult(fff, Zf);
      Zd.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, tmpv2) *= DELTA;
      grad += tmpv2;

      // compute fifth component of gradient (-Rd' * f' * Z' * Z1d * f)
      Z1d.mult(fff, Z1df);
      Z.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (-fff.dot(tmpv));
      grad += tmpv2;

      // compute six and seventh components of gradient 
      // ((1-d) * dX' * [(Z' * Z1d * f) + (Z1d' * Z * f)]
      Z1d.transpose_mult(Zf, tmpv);
      Z.transpose_mult(Z1df, tmpv2);
      tmpv += tmpv2;
      dX.transpose_mult(tmpv, tmpv2) *= ((Real) 1.0 - DELTA);
      grad += tmpv2;

      // compute eight component of gradient (d * Rd' * f' * Zd' * Zd * f)
      Zd.transpose_mult(Zdf, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * DELTA);
      grad += tmpv2;

      // compute ninth component of gradient (d^2 * dX' * Zd' * Zd * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= (DELTA * DELTA);
      grad += tmpv2;

      // compute 10th and 11th components of gradient 
      // (-d * Rd' * f' * Zd' * Z1d * f) and ((1-d) * Rd' * f' * Zd' * Z1d * f)
      Zd.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * (1 - DELTA*2));
      grad += tmpv2;

      // compute 12th and 13th components of gradient
      // ((1-d) * dX * R' * Zd' * Z1d * f) and ((1-d) * d * dX' * Z1d' * Zd * f)
      Zd.transpose_mult(Z1df, tmpv);
      Z1d.transpose_mult(Zdf, tmpv2);
      tmpv += tmpv2;
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= (DELTA - DELTA*DELTA); 
      grad += tmpv2;

      // compute 14th component of gradient (-(1-d) * Rd' * f' * Z1d' * Z1d * f)
      Z1d.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * (DELTA - (Real) 1.0));
      grad += tmpv2;

      // compute 15th component of gradient ((1-d)^2 * dX' * Z1d' * Z1d * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= ((Real) 1.0 - (Real) 2.0*DELTA + DELTA*DELTA);
      grad += tmpv2;

      // scale gradient
      grad *= (-MUCSQ);

      // add in gradient due to applied frictional force
      R.get_row(FIDX, ff);
      ff *= (ff.dot(x) + z[FIDX]);
      grad += ff;
    }
    else
    {
      // compute gradient
      Z.mult(fff, Zf);
      Z.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, grad) *= (-MUCSQ);

      // add in gradient due to applied frictional force
      R.get_row(FIDX, ff);
      ff *= (ff.dot(x) + z[FIDX]);
      grad += ff;
    }

    // set appropriate row of the Jacobian
    J.set_row(index++, grad);
  }
/*
    // setup numerical gradient 
    unsigned n = x.size();
    const Real h = NEAR_ZERO;
    const Real INV_H2 = 1.0 / (h*2);
    VectorN xx = x;
    grad.resize(n);
    for (unsigned i=0; i< n; i++)
    {
      xx[i] += h;
      Real v1 = sqp_fx(xx, m, data);
      xx[i] -= 2*h;
      Real v2 = sqp_fx(xx, m, data);
      xx[i] += h;
      v1 -= v2;
      v1 *= INV_H2;
      grad[i] = v1;
    }
    return;
*/
}

/// The gradient
void ImpactEventHandler::sqp_grad0(const VectorN& x, VectorN& grad, void* data)
{
  // get the optimization data
  const ImpactOptData& opt_data = *(const ImpactOptData*) data;

  // get necessary data
  const MatrixNN& G = opt_data.H;
  const VectorN& c = opt_data.c;

  // objective function is quadratic
  G.mult(x, grad) += c;
}

/// Evaluates objective functions at x
Real ImpactEventHandler::sqp_f0(const VectorN& x, void* data)
{
  SAFESTATIC VectorN Gx;

  // get the optimization data
  const ImpactOptData& opt_data = *(const ImpactOptData*) data;

  // get the event problem data
  const EventProblemData& epd = *opt_data.epd;

  // get necessary data
  const MatrixNN& G = opt_data.H;
  const VectorN& c = opt_data.c;

  // objective function is quadratic
  G.mult(x, Gx) *= (Real) 0.5;
  Gx += c;
  return x.dot(Gx);
} 

/// Evaluates objective and constraint functions at x
void ImpactEventHandler::sqp_fx(const VectorN& x, VectorN& fc, void* data)
{
  SAFESTATIC VectorN Gx, w, tmpv, tmpv2, fff, lambda, walphac, wbetac, walphal;
  SAFESTATIC VectorN wbetat, wbetax;
  SAFESTATIC MatrixN tmpM, tmpM2;
  const Real INFEAS_TOL = 1e-8;

  // get the optimization data
  const ImpactOptData& opt_data = *(const ImpactOptData*) data;

  // get the event problem data
  const EventProblemData& epd = *opt_data.epd;

  // setup constants
  const unsigned N_CONTACTS = epd.N_CONTACTS;
  const unsigned N_LIMITS = epd.N_LIMITS;
  const unsigned N_TRUE_CONE = epd.N_TRUE_CONE;
  const unsigned N_LIN_CONE = epd.N_LIN_CONE;
  const unsigned N_LOOPS = epd.N_LOOPS;
  const unsigned N_CONSTRAINT_EQNS_EXP = epd.N_CONSTRAINT_EQNS_EXP; 
  const unsigned N_CONSTRAINT_DOF_IMP = epd.N_CONSTRAINT_DOF_IMP;
  const unsigned N_CONSTRAINT_DOF_EXP = epd.N_CONSTRAINT_DOF_EXP;
  const unsigned N_JOINT_DOF = N_CONSTRAINT_DOF_EXP + N_CONSTRAINT_DOF_IMP;
  const unsigned ALPHA_C_IDX = 0;
  const unsigned BETA_C_IDX = ALPHA_C_IDX + N_CONTACTS;
  const unsigned NBETA_C_IDX = BETA_C_IDX + N_CONTACTS*2;
  const unsigned ALPHA_L_IDX = NBETA_C_IDX + N_LIN_CONE*2;
  const unsigned ALPHA_X_IDX = ALPHA_L_IDX + N_LIMITS;
  const unsigned BETA_T_IDX = ALPHA_X_IDX + N_CONSTRAINT_EQNS_EXP;
  const unsigned BETA_X_IDX = BETA_T_IDX + N_CONSTRAINT_DOF_IMP;
  const unsigned DELTA_IDX = BETA_X_IDX + N_CONSTRAINT_DOF_EXP;

  // get necessary data
  const VectorN& z = opt_data.z;
  const MatrixN& R = opt_data.R;
  const MatrixNN& G = opt_data.H;
  const VectorN& c = opt_data.c;

  // resize fc
  fc.resize(N_TRUE_CONE + N_LOOPS*2 + N_JOINT_DOF);

  // setup the index
  unsigned index = 0;

  // compute preneeded things
   R.mult(x, w) += z;

  // true friction cone constraint 
  for (unsigned i=0; i< N_TRUE_CONE; i++)
  {
    // get the contact event index
    const unsigned CIDX = opt_data.cone_contacts[i];
    const unsigned BETA_CX = BETA_C_IDX + CIDX*2;
    const unsigned BETA_CY = BETA_CX + 1;
    const unsigned ALPHA_C = ALPHA_C_IDX + CIDX;

    // compute contact friction
    fc[index++] = sqr(w[BETA_CX]) + sqr(w[BETA_CY]) - opt_data.c_mu_c[CIDX]*sqr(w[ALPHA_C]) - opt_data.c_visc[CIDX] - INFEAS_TOL;
  }

  // delta >= 0 constraint
  for (unsigned i=0; i< N_LOOPS; i++)
  {
    const unsigned DIDX = DELTA_IDX + i;
    fc[index++] = -w[DIDX] - INFEAS_TOL;
  }

  // delta <= 1 constraint
  for (unsigned i=0; i< N_LOOPS; i++) 
  {
    // this is delta <= 1
    const unsigned DIDX = DELTA_IDX + i;
    fc[index++] = w[DIDX] - (Real) 1.0 - INFEAS_TOL;
  }

  for (unsigned i=0; i< N_JOINT_DOF; i++)
  {
    // original equation is mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||ff||
    //                   or mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||beta_x||

    // get the body index and body
    const unsigned AIDX = opt_data.body_indices[i];
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(epd.super_bodies[AIDX]);    

   // get indices for this body
    const vector<int>& alpha_c_indices = opt_data.alpha_c_indices[AIDX];
    const vector<int>& beta_nbeta_c_indices = opt_data.beta_nbeta_c_indices[AIDX];
    const vector<unsigned>& alpha_l_indices = opt_data.alpha_l_indices[AIDX];
    const vector<unsigned>& beta_t_indices = opt_data.beta_t_indices[AIDX];
    const vector<unsigned>& beta_x_indices = opt_data.beta_x_indices[AIDX];

    // get the relative joint friction index for this articulated body
    const unsigned RIDX = i - opt_data.joint_friction_start[AIDX];

    // get the index of the joint in the articulated body
    const unsigned JIDX = opt_data.true_indices[AIDX][RIDX];

    // get the number of implicit DOFs for this body, so we can determine
    // whether the relative index is implicit or explicit
    const unsigned N_ABODY_IMPLICIT = abody->num_joint_dof_implicit();

    // get the joint frictional force in w
    const unsigned FIDX = (RIDX < N_ABODY_IMPLICIT) ? BETA_T_IDX + opt_data.implicit_start[AIDX] : BETA_X_IDX + opt_data.explicit_start[AIDX];

    // get the three Z's
    const MatrixN& Zd = opt_data.Zd[AIDX][JIDX];
    const MatrixN& Z1d = opt_data.Z1d[AIDX][JIDX];
    const MatrixN& Z = opt_data.Z[AIDX][JIDX];

    // get components of w and multiply by the proper Jacobians
    contact_select(alpha_c_indices, beta_nbeta_c_indices, w, tmpv, tmpv2);
    abody->transpose_Jc_mult(tmpv, walphac);
    abody->transpose_Dc_mult(tmpv2, wbetac);
    w.select(alpha_l_indices.begin(), alpha_l_indices.end(), tmpv);
    abody->transpose_Jl_mult(tmpv, walphal);
    w.select(beta_t_indices.begin(), beta_t_indices.end(), wbetat);
    w.select(beta_x_indices.begin(), beta_x_indices.end(), tmpv);
    abody->transpose_Dx_mult(tmpv, wbetax);

    // setup impulse 
    fff.copy_from(walphac);
    fff += wbetac;
    fff += walphal;
    fff += wbetat;
    fff += wbetax;
   
    // get the loop indices
    const vector<unsigned>& loop_indices = opt_data.loop_indices[AIDX];

    // two cases: joint is part of a loop or not
    if (loop_indices[JIDX] != std::numeric_limits<unsigned>::max())
    {
      // determine delta stuff
      const unsigned DELTA_START = DELTA_IDX + opt_data.delta_start[AIDX];
      const Real DELTA = w[DELTA_START+loop_indices[JIDX]];
 
      // compute lambda
      Zd.mult(fff, lambda) *= DELTA;
      lambda += (Z1d.mult(fff, tmpv) *= ((Real) 1.0 - DELTA));
      lambda += Z.mult(fff, tmpv);
    }
    else
      Z.mult(fff, lambda);

    // evaluate the equation
    const Real MUCSQ = opt_data.j_mu_c[i];
    const Real VSC = opt_data.j_visc[i];
    fc[index++] = w[FIDX]*w[FIDX] - MUCSQ*lambda.norm_sq() - VSC - INFEAS_TOL;
  }
} 

/// Selects the contact related components from a vector
void ImpactEventHandler::contact_select(const vector<int>& alpha_c_indices, const vector<int>& beta_nbeta_c_indices, const VectorN& x, VectorN& alpha_c, VectorN& beta_c)
{
  // resize vectors
  const unsigned N = alpha_c_indices.size();
  const unsigned M = beta_nbeta_c_indices.size();
  const unsigned TWON = N*2;
  alpha_c.resize(N);
  beta_c.resize(TWON);

  // first get the alpha_c components
  for (unsigned i=0; i< N; i++)
  {
    if (alpha_c_indices[i] > 0)
      alpha_c[i] = x[(unsigned) alpha_c_indices[i]-1];
    else
      alpha_c[i] = -x[(unsigned) -alpha_c_indices[i]-1];
  }

  // now get the beta_c components
  for (unsigned i=0; i< TWON; i++)
  {
    if (beta_nbeta_c_indices[i] > 0)
      beta_c[i] = x[(unsigned) beta_nbeta_c_indices[i]-1];
    else
      beta_c[i] = -x[(unsigned) -beta_nbeta_c_indices[i]-1];
  }

  // now get the nbeta_c components
  for (unsigned i=TWON, j=0; i< M; i++, j++)
  {
    if (beta_nbeta_c_indices[i] > 0)
      beta_c[j] -= x[(unsigned) beta_nbeta_c_indices[i]-1];
    else
      beta_c[j] += x[(unsigned) -beta_nbeta_c_indices[i]-1];
  }
}

/// Selects the contact related rows from a matrix
void ImpactEventHandler::contact_select(const vector<int>& alpha_c_indices, const vector<int>& beta_nbeta_c_indices, const MatrixN& m, MatrixN& alpha_c_rows, MatrixN& beta_c_rows)
{
  SAFESTATIC VectorN tmpv;

  // resize matrices 
  const unsigned N = alpha_c_indices.size();
  const unsigned M = beta_nbeta_c_indices.size();
  const unsigned TWON = N*2;
  alpha_c_rows.resize(N,m.columns());
  beta_c_rows.resize(TWON,m.columns());

  // first get the alpha_c components
  for (unsigned i=0; i< N; i++)
  {
    if (alpha_c_indices[i] > 0)
      alpha_c_rows.set_row(i, m.get_row((unsigned) alpha_c_indices[i]-1, tmpv));
    else
      alpha_c_rows.set_row(i, m.get_row((unsigned) -alpha_c_indices[i]-1, tmpv).negate());
  }

  // now get the beta_c components
  for (unsigned i=0; i< TWON; i++)
  {
    if (beta_nbeta_c_indices[i] > 0)
      beta_c_rows.set_row(i, m.get_row((unsigned) beta_nbeta_c_indices[i]-1, tmpv));
    else
      beta_c_rows.set_row(i, m.get_row((unsigned) -beta_nbeta_c_indices[i]-1, tmpv).negate());
  }

  // now get the nbeta_c components
  for (unsigned i=TWON, j=0; i< M; i++, j++)
  {
    if (beta_nbeta_c_indices[i] > 0)
      beta_c_rows.set_row(j, m.get_row((unsigned) beta_nbeta_c_indices[i]-1, tmpv).negate());
    else
      beta_c_rows.set_row(j, m.get_row((unsigned) -beta_nbeta_c_indices[i]-1, tmpv));
  }
}

