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
#include <Moby/Constraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/ConstraintSimulator.h>
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

// static linear algebra variable
LinAlgd ImpactConstraintHandler::_LA;

/// Sets up the default parameters for the impact event handler
ImpactConstraintHandler::ImpactConstraintHandler()
{
  ip_max_iterations = 100;
  ip_eps = 1e-6;
  use_ip_solver = false;

  // trigger automatic tolerance calculation
  inequality_tolerance = -1.0;

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
void ImpactConstraintHandler::process_constraints(const vector<Constraint>& constraints, double inv_dt)
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // store the set of constraints
  last_constraints = constraints;

  // apply the method to all contacts
  apply_model(constraints, inv_dt);

  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of constraints
/**
 * \param constraints a set of constraints
 */
void ImpactConstraintHandler::apply_model(const vector<Constraint>& constraints, double inv_dt)
{
  const double INF = std::numeric_limits<double>::max();
  list<Constraint*> impacting;
  VectorNd dv, v, f, lambda;

  // **********************************************************
  // determine sets of connected constraints
  // **********************************************************
  list<vector<shared_ptr<DynamicBodyd> > > remaining_islands;
  list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > > groups;
  Constraint::determine_connected_constraints(constraints, groups, remaining_islands);
  Constraint::remove_inactive_groups(groups);

  // **********************************************************
  // do method for each connected set
  // **********************************************************
  for (list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > >::iterator i = groups.begin(); i != groups.end(); i++)
  {
    // copy the lists
    pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > rconstraints = *i;

    // prepare to remove joint limit constraints if they also exist 
    // as an inverse dynamics constraint
    std::vector<JointPtr> id_constraints;
    for (list<Constraint*>::iterator j = i->first.begin(); j != i->first.end(); j++)
      if ((*j)->constraint_type == Constraint::eInverseDynamics)
        id_constraints.push_back((*j)->inv_dyn_joint);

    // sort the inverse dynamics constraints
    std::sort(id_constraints.begin(), id_constraints.end());

    // remove joint limit constraints if they also exist 
    // as an inverse dynamics constraint
    FILE_LOG(LOG_CONSTRAINT) << " -- pre-constraint velocity (all constraints): " << std::endl;
    for (list<Constraint*>::iterator j = i->first.begin(); j != i->first.end(); )
    {
      if ((*j)->constraint_type == Constraint::eLimit && std::binary_search(id_constraints.begin(), id_constraints.end(), (*j)->limit_joint))
        j = i->first.erase(j); 
      else
      {
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;
        j++;
      }
    }

    // apply the model to the constraints
    #ifdef USE_AP_MODEL
    apply_ap_model_to_connected_constraints(rconstraints.first, rconstraints.second, inv_dt);
    #else
    apply_model_to_connected_constraints(rconstraints.first, rconstraints.second, inv_dt);
    #endif

    FILE_LOG(LOG_CONSTRAINT) << " -- post-constraint velocity (all constraints): " << std::endl;
    if (LOGGING(LOG_CONSTRAINT))
    {
      for (list<Constraint*>::iterator j = i->first.begin(); j != i->first.end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;
    }
  }
}

/**
 * Applies method of Drumwright and Shell to a set of connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_model_to_connected_constraints(const list<Constraint*>& constraints, const list<shared_ptr<SingleBodyd> >& single_bodies, double inv_dt)
{
  double ke_minus = 0.0, ke_plus = 0.0;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // compute energy before
  if (LOGGING(LOG_CONSTRAINT))
  {
    // get the super bodies
    vector<shared_ptr<DynamicBodyd> > supers;
    get_super_bodies(constraints.begin(), constraints.end(), std::back_inserter(supers));
    for (unsigned i=0; i< supers.size(); i++)
    {
      double ke = supers[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << supers[i]->body_id << " pre-constraint handling KE: " << ke << endl;
      ke_minus += ke;
    }
  }

  // get the constraints as a vector
  vector<Constraint*> vconstraints(constraints.begin(), constraints.end());

  // apply the model
  apply_model_to_connected_constraints(vconstraints, inv_dt);

  // compute energy after 
  if (LOGGING(LOG_CONSTRAINT))
  {
    // get the super bodies
    vector<shared_ptr<DynamicBodyd> > supers;
    get_super_bodies(constraints.begin(), constraints.end(), std::back_inserter(supers));
    for (unsigned i=0; i< supers.size(); i++)
    {
      double ke = supers[i]->calc_kinetic_energy();
      FILE_LOG(LOG_CONSTRAINT) << "  body " << supers[i]->body_id << " post-constraint handling KE: " << ke << endl;
      ke_plus += ke;
    }
    if (ke_plus > ke_minus)
      FILE_LOG(LOG_CONSTRAINT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
  }

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_model_to_connected_constraints() exiting" << endl;
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

/// Get the total number of variables
unsigned ImpactConstraintHandler::num_variables(const vector<Constraint*>& constraints)
{
  unsigned n_vars = 0;
  for (unsigned m=0; m< constraints.size(); m++)
    n_vars += constraints[m]->num_variables(); 
  return n_vars;
}


