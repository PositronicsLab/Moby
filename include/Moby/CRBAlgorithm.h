/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CRB_ALGORITHM_H
#define _CRB_ALGORITHM_H

#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/MatrixNd.h>

namespace Moby {

/// Computes forward dynamics using composite-rigid body method
class CRBAlgorithm
{
  friend class RCArticulatedBody;

  public:
    CRBAlgorithm();
    ~CRBAlgorithm() {}
    RCArticulatedBodyPtr get_body() const { return RCArticulatedBodyPtr(_body); }
    void set_body(RCArticulatedBodyPtr body) { _body = body; setup_parent_array(); }
    void calc_fwd_dyn();
    void apply_impulse(const Ravelin::SMomentumd& w, RigidBodyPtr link);
    void calc_generalized_inertia(Ravelin::MatrixNd& M);
    void calc_generalized_forces(Ravelin::SForced& f0, Ravelin::VectorNd& C);
    bool factorize_cholesky(Ravelin::MatrixNd& M);
    Ravelin::VectorNd& M_solve(Ravelin::VectorNd& xb);
    Ravelin::MatrixNd& M_solve(Ravelin::MatrixNd& XB);

  private:
    std::vector<unsigned> _lambda;
    void setup_parent_array();

    /// The body that this algorithm operates on
    boost::weak_ptr<RCArticulatedBody> _body;

    /// The spatial acceleration of the base computed on the last call to calc_fwd_dyn()
    Ravelin::SAcceld _a0;

    /// The vector of joint accelerations computed on the last call to calc_fwd_dyn()
    Ravelin::VectorNd _qdd;

    /// The joint space inertia matrix H (fixed base) or augmented matrix [I_0^c K; K^s H] (floating base, see [Featherstone 1987], p. 123) used to compute forward dynamics for floating bases
    Ravelin::MatrixNd _M;

    /// A factorization (or possibly inverse) of the matrix M; note that we compute this b/c we generally may need to solve multiple systems of linear equations using this matrix as a LHS at different times -- always in global frame
    Ravelin::MatrixNd _fM;

    /// Determines whether the system of equations for forward dynamics is rank-deficient
     bool _rank_deficient;

    void calc_joint_space_inertia(RCArticulatedBodyPtr body, Ravelin::MatrixNd& H, std::vector<Ravelin::SpatialRBInertiad>& Ic);
    void apply_coulomb_joint_friction(RCArticulatedBodyPtr body);
    void precalc(RCArticulatedBodyPtr body);
    void calc_generalized_inertia(RCArticulatedBodyPtr body);
    void calc_fwd_dyn_fixed_base(RCArticulatedBodyPtr body);
    void calc_fwd_dyn_floating_base(RCArticulatedBodyPtr body);
    void update_link_accelerations(RCArticulatedBodyPtr body);
    static void to_spatial7_inertia(const Ravelin::SpatialRBInertiad& I, const Ravelin::Quatd& q, Ravelin::MatrixNd& I7);
    Ravelin::VectorNd& M_solve_noprecalc(Ravelin::VectorNd& xb);
    Ravelin::MatrixNd& M_solve_noprecalc(Ravelin::MatrixNd& XB);
    void transform_and_mult(RigidBodyPtr link, const Ravelin::SpatialRBInertiad& I, const std::vector<Ravelin::SAxisd>& s, std::vector<Ravelin::SMomentumd>& Is);

  private:
    // temporaries for transform_and_transpose_mult() functions
    std::vector<Ravelin::SForced> _tandt_fx;
    std::vector<Ravelin::SMomentumd> _tandt_wx;
    std::vector<Ravelin::SAxisd> _tandt_tx;

    // temporary for calc_fwd_dyn() 
    std::vector<Ravelin::SAcceld> _a;

    // temporary for calc_generalized_forces() 
    std::vector<Ravelin::SForced> _w;

    // temporary spatial axes
    std::vector<Ravelin::SAxisd> _sprime;

    // temporaries for solving and linear algebra
    boost::shared_ptr<Ravelin::LinAlgd> _LA;
    Ravelin::MatrixNd _uM, _vM;
    Ravelin::VectorNd _sM;

    // temporaries for calc_generalized_inertia()
    Ravelin::MatrixNd _H;
    std::vector<Ravelin::SpatialRBInertiad> _Ic;
    std::vector<Ravelin::SMomentumd> _Is;

    // temporaries for calc_joint_space_inertia()
    Ravelin::MatrixNd _workM, _sub;
    std::vector<std::vector<bool> > _supports;
    std::vector<std::vector<Ravelin::SMomentumd> > _momenta;

    // temporaries for calc_fwd_dyn_fixed_base(), calc_fwd_dyn_floating_base()
    Ravelin::VectorNd _C, _Q, _Qi, _b, _augV;

    // temporaries for applying impulse
    Ravelin::VectorNd _workv;
    std::vector<Ravelin::SAxisd> _J;

    #include "CRBAlgorithm.inl"
};
}

#endif
