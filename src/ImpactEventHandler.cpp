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
#include <Moby/permute.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/Event.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/NumericalException.h>
#include <Moby/ImpactEventHandler.h>
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
ImpactEventHandler::ImpactEventHandler()
{
  ip_max_iterations = 100;
  ip_eps = 1e-6;
  use_ip_solver = false;
  poisson_eps = NEAR_ZERO;

  // initialize IPOPT, if present
  #ifdef HAVE_IPOPT
  _app.Options()->SetNumericValue("tol", 1e-7);
  _app.Options()->SetStringValue("mu_strategy", "adaptive");
  _app.Options()->SetStringValue("output_file", "ipopt.out");
//  _app.RethrowNonIpoptException(true);

  Ipopt::ApplicationReturnStatus status = _app.Initialize();
  if (status != Ipopt::Solve_Succeeded)
    throw std::runtime_error("Ipopt unable to initialize!");

  // setup the nonlinear IP solver
  _ipsolver = Ipopt::SmartPtr<NQP_IPOPT>(new NQP_IPOPT);  
  _lcpsolver = Ipopt::SmartPtr<LCP_IPOPT>(new LCP_IPOPT);
  #endif
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
    apply_model(events);
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
void ImpactEventHandler::apply_model(const vector<Event>& events)
{
  list<Event*> impacting;

  // **********************************************************
  // determine sets of connected events 
  // **********************************************************
  list<list<Event*> > groups;
  Event::determine_connected_events(events, groups);
  Event::remove_inactive_groups(groups);

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
      apply_model_to_connected_events(revents, 1000.0);

      FILE_LOG(LOG_EVENT) << " -- post-event velocity (all events): " << std::endl;
      for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_EVENT) << "    event: " << std::endl << **j;
  }

  // determine whether there are any impacting events remaining
  for (list<list<Event*> >::const_iterator i = groups.begin(); i != groups.end(); i++)
    for (list<Event*>::const_iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->determine_event_class() == Event::eNegative)
        impacting.push_back(*j);

  // if there are any events still impacting, throw an exception 
  if (!impacting.empty())
    throw ImpactToleranceException(impacting);
}

/**
 * Applies method of Drumwright and Shell to a set of connected events.
 * "Anytime" version
 * \param events a set of connected events 
 */
void ImpactEventHandler::apply_model_to_connected_events(const list<Event*>& events, double max_time)
{
  double ke_minus = 0.0, ke_plus = 0.0;
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  SAFESTATIC EventProblemData epd;

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::apply_model_to_connected_events() entered" << endl;

  // reset problem data
  epd.reset();

  // save the events
  epd.events = vector<Event*>(events.begin(), events.end());

  // determine sets of contact and limit events
  epd.partition_events();

  // compute all event cross-terms
  compute_problem_data(epd);

  // compute energy
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      double ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " pre-event handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // solve the (non-frictional) linear complementarity problem to determine
  // the kappa constant
  VectorNd z;
  solve_lcp(epd, z);

  // update event problem data and z
  permute_problem(epd, z);

  // determine N_ACT_K
  epd.N_ACT_K = 0;
  for (unsigned i=0; i< epd.N_ACT_CONTACTS; i++)
    if (epd.contact_events[i]->contact_NK < UINF)
      epd.N_ACT_K += epd.contact_events[i]->contact_NK/2;

  // use QP / NQP solver with warm starting to find the solution
  if (use_qp_solver(_epd))
    solve_qp(z, epd, poisson_eps, max_time);
  else
    solve_nqp(z, epd, poisson_eps, max_time);

  // apply impulses 
  apply_impulses(epd);

  // compute energy
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< epd.super_bodies.size(); i++)
    {
      double ke = epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " post-event handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_EVENT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::apply_model_to_connected_events() exiting" << endl;

}

/// Permutes the problem to reflect active contact events
void ImpactEventHandler::permute_problem(EventProblemData& epd, VectorNd& z)
{
  // determine mapping of old contact event indices to new contact event
  // indices
  std::vector<unsigned> mapping(epd.N_CONTACTS);

  // 1. compute active indices
  epd.N_ACT_CONTACTS = epd.N_MIN_CONTACTS = 0;
  for (unsigned i=0; i< epd.N_CONTACTS; i++)
    if (z[i] > NEAR_ZERO)
      mapping[epd.N_MIN_CONTACTS++] = i;

  // 2. compute inactive indices
  for (unsigned i=0, j=epd.N_MIN_CONTACTS; i< epd.N_CONTACTS; i++)
    if (z[i] < NEAR_ZERO)
      mapping[j++] = i;

  // permute inactive indices
  std::random_shuffle(mapping.begin()+epd.N_MIN_CONTACTS, mapping.end());  

  // set solution vector to reflect indices in active set
  for (unsigned i=0; i< epd.N_MIN_CONTACTS; i++)
    z[i] = z[mapping[i]];
  std::fill(z.row_iterator_begin()+epd.N_MIN_CONTACTS, z.row_iterator_begin()+epd.N_CONTACTS, 0.0);

  // permute contact events
  std::vector<Event*> new_contact_events(epd.contact_events.size());
  permute(mapping.begin(), mapping.end(), epd.contact_events.begin(), new_contact_events.begin());
  epd.contact_events = new_contact_events;

  // TODO: add event computation and cross computation methods to Joint

  // get iterators to the proper matrices
  RowIteratord CnCn = epd.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = epd.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = epd.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCs = epd.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = epd.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCt = epd.Ct_iM_CtT.row_iterator_begin();

  // process contact events, setting up matrices
  for (unsigned i=0; i< epd.contact_events.size(); i++) 
  {
    // compute cross event data for contact events
    for (unsigned j=0; j< epd.contact_events.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 3);

      // check whether i==j (single contact event)
      if (i == j)
      {
        // compute matrix / vector for contact event i
        _v.set_zero(3);
        epd.contact_events[i]->compute_event_data(_MM, _v);

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
        // compute matrix for cross event
        epd.contact_events[i]->compute_cross_event_data(*epd.contact_events[j], _MM);

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

    // compute cross event data for contact/limit events 
    for (unsigned j=0; j< epd.limit_events.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 1);

      // compute matrix for cross event
      epd.contact_events[i]->compute_cross_event_data(*epd.limit_events[j], _MM);

      // setup appropriate parts of contact / limit inertia matrices
      ColumnIteratord_const data = _MM.column_iterator_begin();
      epd.Cn_iM_LT(i,j) = *data++;
      epd.Cs_iM_LT(i,j) = *data++;
      epd.Ct_iM_LT(i,j) = *data; 
    }
  }

  // NOTE: no need to update limit events of form L_iM_X 
  //       (cross data has already been computed for contact/limit events)
  
  // update cn, l, alpha and kappa 
  z.get_sub_vec(epd.CN_IDX, epd.CS_IDX, epd.cn);
  z.get_sub_vec(epd.L_IDX, epd.ALPHA_X_IDX, epd.l);
  z.get_sub_vec(epd.ALPHA_X_IDX, epd.N_VARS, epd.alpha_x);

  // set active contact indices
  epd.active_contacts.resize(epd.N_CONTACTS);
  for (unsigned i=0; i< epd.N_MIN_CONTACTS; i++)
    epd.active_contacts[i] = true;
  for (unsigned i=epd.N_MIN_CONTACTS; i< epd.N_CONTACTS; i++)
    epd.active_contacts[i] = false;

  // setup active contacts
  epd.N_ACT_CONTACTS = epd.N_MIN_CONTACTS;
}

/**
 * Applies method of Drumwright and Shell to a set of connected events
 * \param events a set of connected events 
 */
void ImpactEventHandler::apply_model_to_connected_events(const list<Event*>& events)
{
  double ke_minus = 0.0, ke_plus = 0.0;
  VectorNd z;

  FILE_LOG(LOG_EVENT) << "ImpactEventHandler::apply_model_to_connected_events() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the events
  _epd.events = vector<Event*>(events.begin(), events.end());

  // determine sets of contact and limit events
  _epd.partition_events();

  // compute all event cross-terms
  compute_problem_data(_epd);

  // compute energy
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << _epd.super_bodies[i]->id << " pre-event handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // mark all contacts as active
  _epd.N_ACT_CONTACTS = _epd.N_MIN_CONTACTS = _epd.N_CONTACTS;
  _epd.N_ACT_K = _epd.N_K_TOTAL;

  // determine what type of QP solver to use
  #ifdef HAVE_IPOPT
  if (use_qp_solver(_epd))
  {
    // solve the (non-frictional) linear complementarity problem to determine
    // the kappa constant
    solve_lcp(_epd, z);
    _epd.set_qp_indices();
    solve_qp(z, _epd, poisson_eps);
  }
  else
  {
    solve_lcp(_epd, z);
    _epd.set_nqp_indices();
    solve_nqp(z, _epd, poisson_eps);
  }
  #else
  solve_lcp(_epd, z);
  _epd.set_qp_indices();
  solve_qp(z, _epd, poisson_eps);
  #endif

  // apply impulses 
  apply_impulses(_epd);

  // compute energy
  if (LOGGING(LOG_EVENT))
  {
    for (unsigned i=0; i< _epd.super_bodies.size(); i++)
    {
      double ke = _epd.super_bodies[i]->calc_kinetic_energy();
      FILE_LOG(LOG_EVENT) << "  body " << _epd.super_bodies[i]->id << " post-event handling KE: " << ke << endl;
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

  // still here? ok to use QP solver
  return true;
}

/// Applies impulses to bodies
void ImpactEventHandler::apply_impulses(const EventProblemData& q)
{
  map<DynamicBodyPtr, VectorNd> gj;
  map<DynamicBodyPtr, VectorNd>::iterator gj_iter;

  // loop over all contact events first
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    // get the contact force
    const Event& e = *q.contact_events[i];
    SForced w(e.contact_impulse);
    const Point3d& p = e.contact_point;

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, p, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, p, _v);
      gj_iter->second += _v; 
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, p, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, p, _v);
      gj_iter->second += _v; 
    }
  }

  // loop over all limit events next
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    const Event& e = *q.limit_events[i];
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

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // initialize constants and set easy to set constants
  q.N_CONTACTS = q.contact_events.size();
  q.N_LIMITS = q.limit_events.size();

  // setup constants related to articulated bodies
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
    if (abody) {
      q.N_CONSTRAINT_EQNS_IMP += abody->num_constraint_eqns_implicit();
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

  // TODO: add event computation and cross computation methods to Joint

  // get iterators to the proper matrices
  RowIteratord CnCn = q.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = q.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = q.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCs = q.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = q.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCt = q.Ct_iM_CtT.row_iterator_begin();

  // process contact events, setting up matrices
  for (unsigned i=0; i< q.contact_events.size(); i++) 
  {
    // compute cross event data for contact events
    for (unsigned j=0; j< q.contact_events.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 3);

      // check whether i==j (single contact event)
      if (i == j)
      {
        // compute matrix / vector for contact event i
        _v.set_zero(3);
        q.contact_events[i]->compute_event_data(_MM, _v);

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
        // compute matrix for cross event
        q.contact_events[i]->compute_cross_event_data(*q.contact_events[j], _MM);

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

    // compute cross event data for contact/limit events 
    for (unsigned j=0; j< q.limit_events.size(); j++)
    {
      // reset _MM
      _MM.set_zero(3, 1);

      // compute matrix for cross event
      q.contact_events[i]->compute_cross_event_data(*q.limit_events[j], _MM);

      // setup appropriate parts of contact / limit inertia matrices
      ColumnIteratord_const data = _MM.column_iterator_begin();
      q.Cn_iM_LT(i,j) = *data++;
      q.Cs_iM_LT(i,j) = *data++;
      q.Ct_iM_LT(i,j) = *data; 
    }
  }

  // process limit events, setting up matrices
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    // compute matrix / vector for contact event i
    q.limit_events[i]->compute_event_data(_MM, _v);

    // setup appropriate entry of limit inertia matrix and limit velocity
    q.L_iM_LT(i,i) = _MM.data()[0];
    q.L_v[i] = _v.data()[0];

    // compute cross/cross limit event data
    for (unsigned j=i+1; j< q.limit_events.size(); j++)
    {
      // reset _MM
      _MM.resize(1,1);

      // compute matrix for cross event
      q.limit_events[i]->compute_cross_event_data(*q.limit_events[j], _MM);

      // setup appropriate part of limit / limit inertia matrix
      q.L_iM_LT(i,j) = q.L_iM_LT(j,i) = _MM.data()[0];
    }

    // NOTE: cross data has already been computed for contact/limit events
  }
}

/// Solves the (frictionless) LCP
void ImpactEventHandler::solve_lcp(EventProblemData& q, VectorNd& z)
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

  // setup the LCP matrix
  _D.mult(_C, _MM);
  _MM -= _B;
  _MM.negate();

  // setup the LCP vector
  _D.mult(_a, _qq);
  _qq -= _b;
  _qq.negate();

  FILE_LOG(LOG_EVENT) << "ImpulseEventHandler::solve_lcp() entered" << std::endl;
  FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_EVENT) << "  Cn * v: " << q.Cn_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  L * v: " << q.L_v << std::endl;
  FILE_LOG(LOG_EVENT) << "  LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_EVENT) << "  LCP vector: " << _qq << std::endl;

  // solve the LCP
  if (!_lcp.lcp_fast(_MM, _qq, _v) && !_lcp.lcp_lemke_regularized(_MM, _qq, _v))
    throw std::runtime_error("Unable to solve event LCP!");

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

