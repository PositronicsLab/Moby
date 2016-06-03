/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _IMPACT_EVENT_HANDLER_H
#define _IMPACT_EVENT_HANDLER_H

#include <list>
#include <vector>
#include <map>
#ifdef USE_QLCPD
#include <Moby/QLCPD.h>
#endif
#ifdef HAVE_IPOPT
#include <coin/IpTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <Moby/NQP_IPOPT.h>
#include <Moby/LCP_IPOPT.h>
#endif
#include <Ravelin/LinAlgd.h>
#include <Ravelin/RigidBodyd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/SparseJacobian.h>
#include <Moby/Constraint.h>

namespace Moby {

class ConstraintSimulator;
class NQP_IPOPT;
class LCP_IPOPT;

/// Defines the mechanism for handling impact constraints
class ImpactConstraintHandler
{
  friend class ConstraintSimulator;
  friend class ConstraintStabilization;

  public:
    ImpactConstraintHandler();
    void process_constraints(const std::vector<Constraint>& constraints, double inv_dt);
    static boost::shared_ptr<Ravelin::DynamicBodyd> get_super_body(boost::shared_ptr<Ravelin::SingleBodyd> sb);

    /// If set to true, uses the interior-point solver (default is false)
    bool use_ip_solver;

    /// The maximum number of iterations to use for the interior-point solver
    unsigned ip_max_iterations;

    /// The tolerance for to the interior-point solver (default 1e-6)
    double ip_eps;

    /// Tolerance for which to solve QP problems (tolerance less than zero triggers automatic tolerance determination)
    double inequality_tolerance;

  private:
    static unsigned num_variables(const std::vector<Constraint*>& constraints);
    static unsigned num_inequality_constraints(const std::vector<Constraint*>& constraints);
    static unsigned num_slackable_constraints(const std::vector<Constraint*>& constraints);
    static unsigned num_equality_constraints(const std::vector<Constraint*>& constraints);
    static double calc_max_constraint_violation(const std::vector<Constraint*>& constraints);
    void form_and_solve(const std::vector<Constraint*>& constraints, double inv_dt, unsigned N_VARS, Ravelin::MatrixNd& H, Ravelin::VectorNd& z);
    void apply_model(const std::vector<Constraint>& constraints, double inv_dt);
    void apply_model_to_connected_constraints(const std::list<Constraint*>& constraints, const std::list<boost::shared_ptr<Ravelin::SingleBodyd> >& single_bodies, double inv_dt);
    void apply_model_to_connected_constraints(const std::vector<Constraint*>& constraints, double dt);
    bool apply_restitution(const std::vector<Constraint*>& constraints, Ravelin::VectorNd& x);
    static void apply_impulses(const std::vector<Constraint*>& constraints, const Ravelin::VectorNd& x);

    static double sqr(double x) { return x*x; }
    static void compute_quadratic_matrix(const std::vector<Constraint*>& constraints, unsigned N_VARS, Ravelin::MatrixNd& H);
    static void compute_linear_term(const std::vector<Constraint*>& constraints, unsigned N_VARS, double inv_dt, Ravelin::VectorNd& c);
    static void compute_equality_terms(const std::vector<Constraint*>& constraints, const Ravelin::MatrixNd& H, unsigned N_EQ_CONSTRAINTS, double inv_dt, Ravelin::MatrixNd& A, Ravelin::VectorNd& b);
    static void compute_inequality_terms(const std::vector<Constraint*>& constraints, const Ravelin::MatrixNd& H, unsigned N_INEQ_CONSTRAINTS, double inv_dt, Ravelin::MatrixNd& M, Ravelin::VectorNd& q);

    LCP _lcp;

    // a pointer to the simulator
    boost::shared_ptr<ConstraintSimulator> _simulator;

    // temporaries for apply_model_to_connected_constraints()
    Ravelin::MatrixNd _H, _M, _A, _MM, _Maug, _Aaug, _Haug, _workM;
    Ravelin::VectorNd _c, _q, _b, _qq, _lb, _ub, _zlast, _caug, _lbaug, _ubaug;
    Ravelin::VectorNd _workv;

    // interior-point solver "application"
    #ifdef HAVE_IPOPT
    Ipopt::IpoptApplication _app;
    #endif

    // temporaries for solve_nqp_work()
    #ifdef HAVE_IPOPT
    Ipopt::SmartPtr <NQP_IPOPT> _ipsolver;
    Ipopt::SmartPtr <LCP_IPOPT> _lcpsolver;
    #endif

    // QLCPD solver
    #ifdef USE_QLCPD
    QLCPD _qp;
    #endif

    // linear algebra variable (for reducing equality constraints)
    static Ravelin::LinAlgd _LA;

    #include "ImpactConstraintHandler.inl"
/*
    // temporaries for IPOPT
    boost::shared_ptr<double> _h_obj, _cJac;
    std::vector<boost::shared_ptr<double> > _h_con;
    unsigned _nnz_h_obj;
    std::vector<unsigned> _nnz_h_con;
    boost::shared_ptr<unsigned> _h_iRow, _h_jCol, _cJac_iRow, _cJac_jCol;
*/
}; // end class
} // end namespace

#endif

