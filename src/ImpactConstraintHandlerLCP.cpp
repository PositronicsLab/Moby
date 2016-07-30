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
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/LCPSolverException.h>
#include <Moby/ImpactConstraintHandler.h>

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

/**
 * Applies Anitescu-Potra model to connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_ap_model_to_connected_constraints(const std::list<Constraint*>& constraints, const list<shared_ptr<SingleBodyd> >& single_bodies, double inv_dt)
{
  double ke_minus = 0.0, ke_plus = 0.0;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() entered" << endl;

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
  apply_ap_model_to_connected_constraints(vconstraints, inv_dt);

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

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() exiting" << endl;
}

/// Applies the AP model to connected constraints
void ImpactConstraintHandler::apply_ap_model_to_connected_constraints(const std::vector<Constraint*>& constraints, double inv_dt)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() entered" << std::endl;

  // mark starting time
  tms cstart;
  clock_t start = times(&cstart);

  // determine the total number of variables
  const unsigned N_VARS = num_variables(constraints);

  // indicate that the LCP matrix must be constructed
  _M.resize(0,0);

  // form and solve the problem
  VectorNd z;
  form_and_solve_ap(constraints, inv_dt, N_VARS, _M, z);

  // apply impulses
  apply_impulses(constraints, z);

  // get the constraint violation resulting from solving the LCP 
  double max_vio_pre = calc_max_constraint_violation(constraints);
  
  // apply restitution
  VectorNd zplus = z;
  if (apply_restitution(constraints, zplus))
  {
    apply_impulses(constraints, zplus);

    // update z
    z += zplus;
    FILE_LOG(LOG_CONSTRAINT) << "z (after restitution): " << z << std::endl;

    // see whether we need to solve another impact problem
    // -- if the constraint violation is negative and the constraint violation
    //    is larger than after solving the inelastic impact QP, another impact 
    //    problem must be solved 
    double max_vio_post = calc_max_constraint_violation(constraints);
    if (max_vio_post > max_vio_pre + NEAR_ZERO)
    {
      // form and solve the problem
      form_and_solve_ap(constraints, inv_dt, N_VARS, _M, z);

      // apply impulses
      apply_impulses(constraints, z);
      FILE_LOG(LOG_CONSTRAINT) << "z (from second impact phase): " << z << std::endl;
    }
  }

  // get the elapsed time
  const long TPS = sysconf(_SC_CLK_TCK);
  tms cstop;
  clock_t stop = times(&cstop);
  double elapsed = (double) (stop - start)/TPS;
  FILE_LOG(LOG_CONSTRAINT) << "Elapsed time: " << elapsed << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() exited" << std::endl;
}

/// Forms the mixed LCP vector
static void compute_mlcp_vector(const vector<Constraint*>& constraints, unsigned N_CONSTRAINTS, double inv_dt, VectorNd& q)
{
  q.resize(N_CONSTRAINTS);
  for (unsigned i=0, j=0; i< constraints.size(); i++)
  {
    assert(constraints[i]->num_constraint_equations() == constraints[i]->num_variables());
    for (unsigned k=0; k< constraints[i]->num_constraint_equations(); k++)
      q[j++] = -constraints[i]->get_constraint_rhs(k, inv_dt);
  }
}

/// Forms the mixed LCP matrix
void ImpactConstraintHandler::compute_mlcp_matrix(const vector<Constraint*>& constraints, unsigned N_CONSTRAINTS, double inv_dt, MatrixNd& M)
{
  // get all super bodies
  vector<shared_ptr<DynamicBodyd> > supers;
  get_super_bodies(constraints.begin(), constraints.end(), std::back_inserter(supers));
  
  // store current velocities
  vector<VectorNd> qd(supers.size());
  for (unsigned i=0; i< supers.size(); i++)
    supers[i]->get_generalized_velocity(DynamicBodyd::eSpatial, qd[i]);

  // initialize the MLCP matrix
  M.resize(N_CONSTRAINTS, N_CONSTRAINTS);

  // compute the MLCP matrix
  for (unsigned m=0, n=0; m< constraints.size(); m++)
  {
    // get the constraint equation
    constraints[m]->get_mlcp_rows(constraints, M, n); 

    // update n
    n += constraints[m]->num_constraint_equations();
  }

  // restore velocities
  for (unsigned i=0; i< supers.size(); i++)
    supers[i]->set_generalized_velocity(DynamicBodyd::eSpatial, qd[i]);
}

/// Forms and solves the impact problem
void ImpactConstraintHandler::form_and_solve_ap(const vector<Constraint*>& constraints, double inv_dt, unsigned N_VARS, MatrixNd& MM, VectorNd& z)
{
  const double INF = (double) std::numeric_limits<double>::max();

  // compute the LCP matrix if necessary
  if (MM.rows() == 0)
  {
    compute_mlcp_matrix(constraints, N_VARS, inv_dt, MM);

    // add in damping terms
    for (unsigned m=0, o=0; m< constraints.size(); m++)
      for (unsigned n=0; n< constraints[m]->num_constraint_equations(); n++, o++)
      {
        // get the damping term (if any)
        double damping = constraints[m]->get_damping(n);
        if (damping <= 0.0)
          continue;

        // update H
        MM(o,o) += damping;
      }
  }

  // compute the linear term vector
  compute_mlcp_vector(constraints, N_VARS, inv_dt, _qq);

  // setup the lower and upper bounds variables; simultaneously determine
  // whether the problem is a pure LCP
  bool pure_lcp = true;
  _lb.resize(N_VARS);
  _ub.resize(N_VARS);
  for (unsigned i=0, j=0; i< constraints.size(); i++)
    for (unsigned k=0; k< constraints[i]->num_variables(); k++)
    {
      _lb[j] = constraints[i]->get_lower_variable_limit(k);
      _ub[j] = constraints[i]->get_upper_variable_limit(k);
      if (std::fabs(_lb[j]) > NEAR_ZERO || _ub[j] < INF)
        pure_lcp = false;
      j++;
    } 

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::form_and_solve()" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "M: " << std::endl  << MM;
  FILE_LOG(LOG_CONSTRAINT) << "q: " << _qq << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "lb: " << _lb << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ub: " << _ub << std::endl;

  // attempt to solve using the fast regularized solver
  if (!_lcp.mlcp_fast_regularized(MM, _qq, _lb, _ub, z, -16, 4, -4))
  {
    // if the LCP is not truly mixed, we can attempt Lemke's algorithm
    if (pure_lcp)
    {
      // Lemke does not like warm starting
      z.set_zero();

      FILE_LOG(LOG_CONSTRAINT) << "-- fast LCP solver failed; dropping back to Lemke" << std::endl;
      if (!_lcp.lcp_lemke_regularized(MM, _qq, z))
          throw LCPSolverException();
    }
    else
      throw LCPSolverException();
  }

  FILE_LOG(LOG_CONSTRAINT) << "mixed LCP solution z: " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "mixed LCP pivots: " << _lcp.pivots << std::endl;
}


