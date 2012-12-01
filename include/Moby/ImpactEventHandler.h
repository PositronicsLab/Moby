/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _IMPACT_EVENT_HANDLER_H
#define _IMPACT_EVENT_HANDLER_H

#include <list>
#include <vector>
#include <map>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/Event.h>

namespace Moby {

/// Defines the mechanism for handling impact events 
class ImpactEventHandler
{
  private:
    struct ImpactOptData
    {
      /// Homogeneous solution
      VectorN z;

      /// Nullspace for optimization
      MatrixN R;

      /// Quadratic optimization matrix
      MatrixN H;

      /// Linear optimization vector
      VectorN c;

      /// pointer to the EventProblemData
      const EventProblemData* epd;

      /// maps a true cone contact (0..N_TRUE_CONE-1) to its contact event index
      std::vector<unsigned> cone_contacts;

      /// Squared Coulomb friction coefficients for contacts (size N_TRUE_CONE)
      std::vector<Real> c_mu_c;

      /// Squared viscous terms for contacts (size N_TRUE_CONE); this is squared viscous friction coefficient times squared tangential contact velocity
      std::vector<Real> c_visc;

      /// Squared Coulomb friction coefficients for joints (size N_JOINT_DOF)
      std::vector<Real> j_mu_c;

      /// Squared viscous terms for joints (size N_JOINT_DOF); this is squared viscous friction coefficient times squared joint velocity
      std::vector<Real> j_visc;

      /// Z matrices (for joint friction) for each joint of each super body (actually articulated body)
      std::vector<std::vector<MatrixN> > Z;

      /// Zd matrices (for joint friction) for each joint of each super body (actually articulated body)
      std::vector<std::vector<MatrixN> > Zd;

      /// Z1d matrices (for joint friction) for each joint of each super body (actually articulated body)
      std::vector<std::vector<MatrixN> > Z1d;

      /// loop indices (loop that each joint belongs to) for each articulated body
      /**
       * The top level vector is for each super body (actually articulated 
       * body); the second level vector maps from the [implicit; explicit]
       * joint for the articulated body to the loop index containing that 
       * joint. If the joint index is not a member of a loop, it is mapped
       * to std::numeric_limits<unsigned>::max(), aka, UINT_MAX.
       */
      std::vector<std::vector<unsigned> > loop_indices;

      /// maps joint friction indices to the containing articulated body in the EventProblemData super bodies vector
      std::vector<unsigned> body_indices;

      /// maps articulated body to the start of its delta indices (relative to the start of *all* delta indices)
      std::vector<unsigned> delta_start;

      /// maps articulated body to the start of its friction indices in the optimization problem (first articulated body will have index 0, second articulated body will have index n [n is the number of joint degrees-of-freedom in the first articulated body], etc.)
      std::vector<unsigned> joint_friction_start; 

      /// maps articulated body to the start of its implicit joint friction indices (beta_t, relative to the start of *all* beta_t indices)
      std::vector<unsigned> implicit_start;

      /// maps articulated body to the start of its explicit joint friction indices (beta_x, relative to the start of *all* beta_x indices)
      std::vector<unsigned> explicit_start; 

      /// maps articulated bodies to "true" indices
      /**
       * The top level vector is for each super body (actually articulated
       * body); the second level vector maps from the [implicit; explicit]
       * index for a joint to the true (get_index()) index for that joint in
       * the articulated body.
       */
      std::vector<std::vector<unsigned> > true_indices;

      /// gets alpha_c indices for articulated bodies
      std::vector<std::vector<int> > alpha_c_indices;

      /// gets beta_c indices for articulated bodies
      std::vector<std::vector<int> > beta_nbeta_c_indices;

      /// gets alpha_l indices for articulated bodies
      std::vector<std::vector<unsigned> > alpha_l_indices;

      /// gets beta_t indices for articulated bodies
      std::vector<std::vector<unsigned> > beta_t_indices;

      /// gets beta_x indices for articulated bodies
      std::vector<std::vector<unsigned> > beta_x_indices;
    };

  public:
    ImpactEventHandler();
    void process_events(const std::vector<Event>& events, Real tol = NEAR_ZERO);

    /// If set to true, uses the interior-point solver (default is false)
    bool use_ip_solver;

    /// The maximum number of iterations to use for the interior-point solver 
    unsigned ip_max_iterations;

    /// The tolerance for to the interior-point solver (default 1e-6)
    Real ip_eps;

    /// The velocity tolerance above which another iteration of the solver is run after applying Poisson restitution
    Real poisson_eps;

  private:
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);
    static bool use_qp_solver(const EventProblemData& epd);
    void apply_model(const std::vector<Event>& events, Real tol) const;
    void apply_model_to_connected_events(const std::list<Event*>& events) const;
    static void compute_problem_data(EventProblemData& epd);
    static void solve_lcp(EventProblemData& epd, VectorN& z);
    static void solve_qp(EventProblemData& epd, Real eps);
    static void solve_nqp(EventProblemData& epd, Real eps);
    static void solve_qp_work(EventProblemData& epd, VectorN& z);
    static void solve_nqp_work(EventProblemData& epd, VectorN& z);
    static void set_generalized_velocities(const EventProblemData& epd);
    static void partition_events(const std::list<Event*>& events, std::vector<Event*>& contacts, std::vector<Event*>& limits);
    static void add_constraint_events(const std::list<Event*>& events, std::vector<Event>& constraint_event_objects, std::vector<Event*>& constraint_events);
    static void update_impulses(EventProblemData& q, const VectorN& z);
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const VectorN& x, VectorN& alpha_c, VectorN& beta_c);
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const MatrixN& m, MatrixN& alpha_c_rows, MatrixN& beta_c_rows);
    static Real sqr(Real x) { return x*x; }
    static void sqp_hess(const VectorN& x, Real objscal, const VectorN& lambda, const VectorN& nu, MatrixN& H, void* data);
    static void sqp_grad0(const VectorN& x, VectorN& g, void* data);
    static void sqp_cJac(const VectorN& x, MatrixN& J, void* data);
    static Real sqp_f0(const VectorN& x, void* data);
    static void sqp_fx(const VectorN& x, VectorN& fc, void* data);
    static void set_optimization_data(EventProblemData& q, ImpactOptData& iopt);
}; // end class
} // end namespace

#endif

