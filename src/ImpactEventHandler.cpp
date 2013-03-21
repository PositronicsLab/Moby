/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
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
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
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
void ImpactEventHandler::process_events(const vector<Event>& events)
{
  FILE_LOG(LOG_EVENT) << "*************************************************************";
  FILE_LOG(LOG_EVENT) << endl;
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::process_events() entered";
  FILE_LOG(LOG_EVENT) << endl;
  FILE_LOG(LOG_EVENT) << "*************************************************************";
  FILE_LOG(LOG_EVENT) << endl;

  // apply the method to all contacts
  if (!events.empty())
  {
    vector<Event> nev;
    nev.push_back(events.front());
    apply_model(events);
  }
  else
    FILE_LOG(LOG_EVENT) << " (no events?!)" << endl;
    
  FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::process_events() exited" << endl;
  FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
}

/// Applies the model to a set of events 
/**
 * \param events a set of events
 */
void ImpactEventHandler::apply_model(const vector<Event>& events) const
{
  list<Event*> impacting;

  // **********************************************************
  // determine sets of connected events 
  // **********************************************************
  list<list<Event*> > groups;
  Event::determine_connected_events(events, groups);
  Event::remove_nonimpacting_groups(groups);

  // **********************************************************
  // do method for each connected set 
  // **********************************************************
  for (list<list<Event*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
    // determine contact tangents
    for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->event_type == Event::eContact)
        (*j)->determine_contact_tangents();

      // copy the list of events
      list<Event*> revents = *i;

      FILE_LOG(LOG_EVENT) << " -- pre-event velocity (all events): " << std::endl;
      for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_EVENT) << "    event: " << std::endl << **j;

      // determine a reduced set of events
      Event::determine_minimal_set(revents);

      // apply model to the reduced contacts   
      apply_model_to_connected_events(revents);

      FILE_LOG(LOG_EVENT) << " -- post-event velocity (all events): " << std::endl;
      for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_EVENT) << "    event: " << std::endl << **j;
  }

  // determine whether there are any impacting events remaining
  for (list<list<Event*> >::const_iterator i = groups.begin(); i != groups.end(); i++)
    for (list<Event*>::const_iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->is_impacting())

  // if there are any events still impacting, throw an exception 
  if (!impacting.empty())
    throw ImpactToleranceException(impacting);
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

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::apply_model_to_connected_events() entered" << endl;

  // reset problem data
  epd.reset();

  // determine sets of contact and limit events
  partition_events(events, epd.contact_events, epd.limit_events);

  // add events for constraints for articulated bodies
  add_constraint_events(events, constraint_event_objects, epd.constraint_events);

  // compute all event cross-terms
  compute_problem_data(epd);

  // compute energy
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      Real ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " pre-event handling KE: " << ke << endl;
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
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      Real ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " post-event handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_EVENT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::apply_model_to_connected_events() exiting" << endl;
}

/// Determines whether we can use the QP solver
bool ImpactEventHandler::use_qp_solver(const EventProblemData& epd)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // first, check whether any contact events use a true friction cone
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
    if (epd.contact_events[i]->contact_NK == UINF)
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
  const unsigned UINF = std::numeric_limits<unsigned>::max();
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

  // now, sort the contact events such that events that use a true friction
  // cone are at the end
  for (unsigned i=0, j=contacts.size()-1; i< j; )
  {
    if (contacts[i]->contact_NK == UINF)
    {
      std::swap(contacts[i], contacts[j]);
      j--;
    } 
    else
      i++;
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

  // setup contact working set
  q.contact_working_set.clear();
  q.contact_working_set.resize(q.N_CONTACTS, true);

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

  // compute number of friction polygon edges
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    if (q.contact_events[i]->contact_NK < UINF)
    {
      q.N_K_TOTAL += q.contact_events[i]->contact_NK/2;
      q.N_LIN_CONE++;
    }
    else if (q.contact_events[i]->contact_NK == UINF)
      break;
  }

  // setup number of true cones
  q.N_TRUE_CONE = q.contact_events.size() - q.N_LIN_CONE; 

  // verify contact constraints that use a true friction cone are at the end 
  // of the contact vector
  #ifndef NDEBUG
  for (unsigned i=q.N_LIN_CONE; i< q.contact_events.size(); i++)
    assert(q.contact_events[i]->contact_NK == UINF);
  #endif
   
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

  // setup indices
  q.ALPHA_C_IDX = 0;
  q.BETA_C_IDX = q.ALPHA_C_IDX + q.N_CONTACTS;
  q.NBETA_C_IDX = q.BETA_C_IDX + q.N_LIN_CONE*2;
  q.BETAU_C_IDX = q.NBETA_C_IDX + q.N_LIN_CONE*2;
  q.ALPHA_L_IDX = q.BETAU_C_IDX + q.N_TRUE_CONE;
  q.BETA_T_IDX = q.ALPHA_L_IDX + q.N_LIMITS;
  q.ALPHA_X_IDX = q.BETA_T_IDX + q.N_CONSTRAINT_DOF_IMP;
  q.BETA_X_IDX = q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_EXP;
  q.N_VARS = q.BETA_X_IDX + q.N_CONSTRAINT_DOF_EXP;

  // for each super body, update the problem data
  for (unsigned i=0; i< q.super_bodies.size(); i++)
    q.super_bodies[i]->update_event_data(q);
}

/// Solves the (frictionless) LCP
void ImpactEventHandler::solve_lcp(EventProblemData& q, VectorN& z)
{
  SAFESTATIC MatrixN UL, LR, MM;
  SAFESTATIC MatrixN UR, t2, iJx_iM_JxT;
  SAFESTATIC VectorN alpha_c, alpha_l, alpha_x, v1, v2, qq;

  // setup sizes
  UL.resize(q.N_CONTACTS, q.N_CONTACTS);
  UR.resize(q.N_CONTACTS, q.N_LIMITS);
  LR.resize(q.N_LIMITS, q.N_LIMITS);

  // setup primary terms -- first upper left hand block of matrix
  iJx_iM_JxT.copy_from(q.Jx_iM_JxT);
  try
  {
    LinAlg::pseudo_inverse(iJx_iM_JxT, LinAlg::svd1);
  }
  catch (NumericalException e)
  {
    iJx_iM_JxT.copy_from(q.Jx_iM_JxT);
    LinAlg::pseudo_inverse(iJx_iM_JxT, LinAlg::svd2);
  }
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
  MM.resize(q.N_CONTACTS + q.N_LIMITS, q.N_CONTACTS + q.N_LIMITS);
  MM.set_sub_mat(0, 0, UL);
  MM.set_sub_mat(0, q.N_CONTACTS, UR);
  MM.set_sub_mat(q.N_CONTACTS, 0, UR, true);
  MM.set_sub_mat(q.N_CONTACTS, q.N_CONTACTS, LR);

  // setup the LCP vector
  qq.resize(MM.rows());
  iJx_iM_JxT.mult(q.Jx_v, v2); 
  q.Jc_iM_JxT.mult(v2, v1);
  v1 -= q.Jc_v;
  qq.set_sub_vec(0, v1);
  q.Jl_iM_JxT.mult(v2, v1);
  v1 -= q.Jl_v;
  qq.set_sub_vec(q.N_CONTACTS, v1);
  qq.negate();

  FILE_LOG(LOG_EVENT) << "ImpulseEventHandler::solve_lcp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jc * inv(M) * Jc': " << std::endl << q.Jc_iM_JcT;
  FILE_LOG(LOG_EVENT) << "  Jc * v: " << q.Jc_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  Jl * v: " << q.Jl_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  LCP matrix: " << std::endl << MM;
  FILE_LOG(LOG_EVENT) << "  LCP vector: " << qq << std::endl;

  // solve the LCP
  if (!Optimization::lcp_lemke_regularized(MM, qq, z))
    throw std::runtime_error("Unable to solve event LCP!");

  // determine the value of kappa
  q.kappa = (Real) 0.0;
  for (unsigned i=0; i< q.N_CONTACTS; i++)
    q.kappa += z[i];

  // get alpha_c and alpha_l
  z.get_sub_vec(0, q.N_CONTACTS, alpha_c);
  z.get_sub_vec(q.N_CONTACTS, z.size(), alpha_l);

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
  z.set_zero();
  z.set_sub_vec(q.ALPHA_C_IDX, alpha_c);
  z.set_sub_vec(q.ALPHA_L_IDX, alpha_l);
  z.set_sub_vec(q.ALPHA_X_IDX, alpha_x);

  FILE_LOG(LOG_EVENT) << "  LCP result: " << z << std::endl;
  FILE_LOG(LOG_EVENT) << "  kappa: " << q.kappa << std::endl;
  FILE_LOG(LOG_EVENT) << "ImpulseEventHandler::solve_lcp() exited" << std::endl;
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

