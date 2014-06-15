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
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/Event.h>
#include <Moby/EventProblemData.h>

namespace Moby {

class NQP_IPOPT;
class LCP_IPOPT;

/// Defines the mechanism for handling impact events 
class ImpactEventHandler
{
  public:
    ImpactEventHandler();
    void process_events(const std::vector<Event>& events, double max_time);

    /// If set to true, uses the interior-point solver (default is false)
    bool use_ip_solver;

    /// The maximum number of iterations to use for the interior-point solver 
    unsigned ip_max_iterations;

    /// The tolerance for to the interior-point solver (default 1e-6)
    double ip_eps;

  private:
    void apply_visc_friction_model_to_connected_events(const std::list<Event*>& events);
    void apply_inf_friction_model_to_connected_events(const std::list<Event*>& events);
    void update_from_stacked(EventProblemData& q, const Ravelin::VectorNd& z);
    double calc_min_constraint_velocity(const EventProblemData& q) const;
    void update_event_velocities_from_impulses(EventProblemData& q);
    bool apply_restitution(const EventProblemData& q, Ravelin::VectorNd& z) const;
    bool apply_restitution(EventProblemData& q) const;
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);
    static bool use_qp_solver(const EventProblemData& epd);
    void apply_model(const std::vector<Event>& events, double max_time);
    void apply_model_to_connected_events(const std::list<Event*>& events);
    void apply_model_to_connected_events(const std::list<Event*>& events, double max_time);
    void compute_problem_data(EventProblemData& epd);
    void solve_frictionless_lcp(EventProblemData& epd, Ravelin::VectorNd& z);
    void apply_visc_friction_model(EventProblemData& epd);
    void apply_inf_friction_model(EventProblemData& epd);
    void solve_qp(Ravelin::VectorNd& z, EventProblemData& epd, double max_time = std::numeric_limits<double>::max());
    void solve_nqp(Ravelin::VectorNd& z, EventProblemData& epd, double max_time = std::numeric_limits<double>::max());
    void solve_qp_work(EventProblemData& epd, Ravelin::VectorNd& z);
    double calc_ke(EventProblemData& epd, const Ravelin::VectorNd& z);
    void update_problem(const EventProblemData& qorig, EventProblemData& qnew);
    void update_solution(const EventProblemData& q, const Ravelin::VectorNd& x, const std::vector<bool>& working_set, unsigned jidx, Ravelin::VectorNd& z);
    void solve_nqp_work(EventProblemData& epd, Ravelin::VectorNd& z);
    void apply_impulses(const EventProblemData& epd);
    static void contact_select(const std::vector<int>& cn_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::VectorNd& x, Ravelin::VectorNd& cn, Ravelin::VectorNd& beta_c);
    static void contact_select(const std::vector<int>& cn_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& cn_rows, Ravelin::MatrixNd& beta_c_rows);
    static double sqr(double x) { return x*x; }
    void permute_problem(EventProblemData& epd, Ravelin::VectorNd& z);
    void setup_QP(EventProblemData& epd, Ravelin::SharedMatrixNd& H, Ravelin::SharedVectorNd& c, Ravelin::SharedMatrixNd& M, Ravelin::SharedVectorNd& q, Ravelin::SharedMatrixNd& A, Ravelin::SharedVectorNd& b);

    Ravelin::LinAlgd _LA;
    LCP _lcp;

    // persistent event data
    EventProblemData _epd;

    // temporaries for compute_problem_data(), solve_qp_work(), solve_lcp(), and apply_impulses()
    Ravelin::MatrixNd _MM;
    Ravelin::VectorNd _zlast, _v, _zsuccess;

    // temporaries for solve_qp_work() and solve_nqp_work() 
    Ravelin::VectorNd _Cnstar_v, _workv, _new_Cn_v;
    Ravelin::MatrixNd _Cnstar_Cn, _Cnstar_Cs, _Cnstar_Ct, _Cnstar_L;

    // temporaries shared between solve_lcp(), solve_qp(), and solve_nqp()
    Ravelin::VectorNd _a, _b;

    // temporaries for solve_qp() and solve_nqp()
    Ravelin::VectorNd _z;

    // temporaries for solve_frictionless_lcp()
    Ravelin::VectorNd _cs_visc, _ct_visc;

    // temporaries for solve_inf_friction_lcp()
    Ravelin::MatrixNd _rJx_iM_JxT, _Y, _Q_iM_XT, _workM, _workM2;
    Ravelin::VectorNd _YXv, _Xv, _cs_ct_alphax;

    // interior-point solver "application"
    #ifdef HAVE_IPOPT
    Ipopt::IpoptApplication _app;
    #endif

    // temporaries shared between solve_nqp_work() and solve_lcp()
    Ravelin::MatrixNd _A;

    // temporaries for solve_nqp_work()
    #ifdef HAVE_IPOPT
    Ipopt::SmartPtr <NQP_IPOPT> _ipsolver;
    Ipopt::SmartPtr <LCP_IPOPT> _lcpsolver;
    #endif
    Ravelin::MatrixNd _RTH;
    Ravelin::VectorNd _w, _workv2, _x;

    // temporaries for solve_lcp()
    Ravelin::MatrixNd _AU, _AV, _B, _C, _D;
    Ravelin::VectorNd _AS, _alpha_x, _qq, _Cn_vplus;

    // last number of (active) contacts handled
    unsigned _last_contacts;

    // last number of limits handled
    unsigned _last_limits;

    // last number of friction directions
    unsigned _last_contact_nk;

    // last number of contact constraints handled
    unsigned _last_contact_constraints;

    // QLCPD solver
    #ifdef USE_QLCPD
    QLCPD _qp;
    #endif

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

