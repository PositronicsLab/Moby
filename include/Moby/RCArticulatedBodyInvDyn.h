/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RC_ARTICULATED_BODY_INV_DYN_ALGO_H
#define _RC_ARTICULATED_BODY_INV_DYN_ALGO_H

#include <map>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/RigidBodyd.h>
#include <Moby/Base.h>
#include <Moby/LCP.h>
#include <Moby/Constraint.h>

namespace Moby {

class RCArticulatedBody;

/// Data structure for passing data to inverse dynamics algorithm
class RCArticulatedBodyInvDynData
{
  public:

    /// Desired inner joint accelerations
    Ravelin::VectorNd qdd_des;

    /// Constraints applied to the body
    std::vector<Constraint> constraints;
};
  
/// Class for performing inverse dynamics computation on a reduced-coordinate articulated body
class RCArticulatedBodyInvDyn : public virtual Base
{
  public:
    RCArticulatedBodyInvDyn() { }
    virtual ~RCArticulatedBodyInvDyn() {}

    /**
     * Given the current base position and velocity, the joint positions and 
     * joint velocities, the external forces on all links (including the 
     * base), and the desired joint accelerations, the inverse 
     * dynamics method must compute the generalized forces that achieve these 
     * accelerations.
     */
    void calc_inv_dyn(boost::shared_ptr<RCArticulatedBody> body, const RCArticulatedBodyInvDynData& inv_dyn_data, double dt, Ravelin::VectorNd& u);

  private:
    void calc_contact_jacobians(const std::vector<Constraint>& c, Ravelin::MatrixNd& N, Ravelin::MatrixNd& S, Ravelin::MatrixNd& T, const std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& index_start, unsigned NDOFS);
    void jacobian_ST_to_D(const Ravelin::MatrixNd& S, const Ravelin::MatrixNd& T, Ravelin::MatrixNd& D);
    bool inverse_dynamics_comp(RCArticulatedBodyPtr body, const std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& start_indices, const Ravelin::VectorNd& vel, const Ravelin::VectorNd& qdd, const std::vector<Ravelin::MatrixNd>& M,const  Ravelin::MatrixNd& NT, const Ravelin::MatrixNd& D_, const Ravelin::VectorNd& fext, double dt, const std::vector<double>& mu, Ravelin::VectorNd& x, Ravelin::VectorNd& cf);
    bool inverse_dynamics_two_stage_simple_no_slip(const Ravelin::VectorNd& v, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& N, const Ravelin::MatrixNd& ST, const Ravelin::VectorNd& fext, double h, Ravelin::VectorNd& x, Ravelin::VectorNd& cf_final,double damping);
    bool inverse_dynamics_one_stage(const Ravelin::VectorNd& v, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& N, const Ravelin::MatrixNd& ST, const Ravelin::VectorNd& fext_, double h, const Ravelin::MatrixNd& MU, Ravelin::VectorNd& x, Ravelin::VectorNd& cf_final, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS);
    bool inverse_dynamics_no_contact(const Ravelin::VectorNd& vel, const Ravelin::VectorNd& qdd_des, const Ravelin::MatrixNd& M, const Ravelin::VectorNd& fext, double dt, Ravelin::VectorNd& x);
    bool inverse_dynamics_no_slip_fast(const Ravelin::VectorNd& vel, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& nT, const Ravelin::MatrixNd& D, const Ravelin::VectorNd& fext, double dt, Ravelin::VectorNd& x, Ravelin::VectorNd& cf, bool frictionless, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS);
    bool inverse_dynamics_two_stage_simple(const Ravelin::VectorNd& v, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& N, const Ravelin::MatrixNd& ST, const Ravelin::VectorNd& fext, double h, const Ravelin::MatrixNd& MU, Ravelin::VectorNd& x, Ravelin::VectorNd& cf_final, double damping);
    bool inverse_dynamics_two_stage(const Ravelin::VectorNd& v, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& N, const Ravelin::MatrixNd& ST, const Ravelin::VectorNd& fext_, double h, const Ravelin::MatrixNd& MU, Ravelin::VectorNd& x, Ravelin::VectorNd& cf_final, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS);
    bool inverse_dynamics_no_slip(const Ravelin::VectorNd& v, const Ravelin::VectorNd& qdd, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& N, const Ravelin::MatrixNd& ST, const Ravelin::VectorNd& fext_, double h, Ravelin::VectorNd& x, Ravelin::VectorNd& cf_final, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS);
    bool solve_qp_pos(const Ravelin::MatrixNd& Q, const Ravelin::VectorNd& c, Ravelin::VectorNd& x, bool warm_start);
    bool solve_qp(const Ravelin::MatrixNd& Q, const Ravelin::VectorNd& c, const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
    bool solve_qp_pos(const Ravelin::MatrixNd& Q, const Ravelin::VectorNd& c, const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, Ravelin::VectorNd& x, Ravelin::VectorNd& v, bool warm_start = false, bool regularize = true);
    bool lcp_fast(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, const std::vector<unsigned>& indices, Ravelin::VectorNd& z, double zero_tol);
    static unsigned get_body_indices(const std::vector<Constraint>& c, std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& index_start);
    static bool factor_chols(const std::vector<Ravelin::MatrixNd>& M, std::vector<Ravelin::MatrixNd>& M_chols);
    static void solve_chol_fast(const std::vector<Ravelin::MatrixNd>& M_chols, Ravelin::MatrixNd& M);
    static void solve_chol_fast(const std::vector<Ravelin::MatrixNd>& M_chols, Ravelin::VectorNd& v);
    static unsigned select_pivot(const Ravelin::VectorNd& znbas, const std::vector<unsigned>& nonbas, const std::vector<unsigned>& indices, double zero_tol = NEAR_ZERO);
    static unsigned rand_min(const Ravelin::VectorNd& v, double zero_tol);

    std::vector<Ravelin::MatrixNd> _M_chols;
    Moby::LCP _lcp;
    Ravelin::LinAlgd _LA;
    Ravelin::VectorNd _workv, _workv2, _qq, _Xv, _YXv, _Yvqstar, _cs_ct_tau;
    Ravelin::VectorNd _v_plus, _z, _qbas, _w;
    Ravelin::MatrixNd _workM, _workM2, _MM, _R, _iM, _M_chol, _D, _DT, _N, _S;
    Ravelin::MatrixNd _Q_iM_XT, _Y, _P, _X, _iX, _PT, _E, _ET, _F, _F_chol;
    Ravelin::MatrixNd _ST, _TT, _rJx_iM_JxT, _T, _B, _LL, _UR;
    Ravelin::MatrixNd _Cn_iM_CnT, _Cn_iM_CsT, _Cn_iM_CtT, _Cn_iM_CdT;
    Ravelin::MatrixNd _Cn_iM_JxT, _Jx_iM_JxT, _Msub, _Mmix; 
    Ravelin::MatrixNd _Cs_iM_CsT, _Cs_iM_CtT, _Ct_iM_CtT; 
    Ravelin::MatrixNd _Cd_iM_CdT, _Cd_iM_CnT, _Cd_iM_JxT;
    std::vector<bool> _represented;
    std::vector<unsigned> _nonbas, _bas;
 
}; // end class
} // end namespace

#endif

