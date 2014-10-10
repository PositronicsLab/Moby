/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _FS_AB_ALGORITHM_H
#define _FS_AB_ALGORITHM_H

#include <queue>
#include <Ravelin/SpatialABInertiad.h>

namespace Moby {

/// Implements Featherstone's algorithm for forward dynamics
/**
 * Implements Featherstone's algorithm for articulated bodies.  Featherstone's 
 * algorithm runs in O(n) time [n = # of joints].  This implementation is based
 * on Brian Mirtich's Ph. D. thesis, and remains pretty consistent with it. 
 * There are a couple of changes, to produce a nice implementation.  The user 
 * need not be concerned with these issues, but they are useful to know for 
 * debugging.
 * <ol>
 * <li>Mirtich labels his links from 1..n, and considers the base to be link 0; 
 * the total number of links is considered to be n, rather than n+1.  I make 
 * the total number of links n+1 and treat the links the same as the base.  I 
 * do this so that the user thinks of the base as a link for purposes of link 
 * connectivity.</li>
 * <li>Mirtich labels his joints from 0..n-1</li>.  When labeling the link in 
 * Mirtich's style, link i and joint i match up (joint i is link i's inner 
 * joint).  When labeling the link in my style, joint i-1 is the corresponding 
 * joint for link i.</li>
 * </ol>
 * Note that one critical note for manipulator setup is that the base is the 
 * first link in the list of links.
 */
class FSABAlgorithm 
{
  friend class RCArticulatedBody;

  public:
    FSABAlgorithm();
    ~FSABAlgorithm() {}
    RCArticulatedBodyPtr get_body() const { return RCArticulatedBodyPtr(_body); }
    void set_body(RCArticulatedBodyPtr body) { _body = body; }
    void calc_fwd_dyn();
    void calc_inverse_generalized_inertia_noprecalc(Ravelin::MatrixNd& iM);
    void solve_generalized_inertia_noprecalc(Ravelin::SharedVectorNd& v);
    void solve_generalized_inertia_noprecalc(Ravelin::SharedMatrixNd& Y);
    void apply_generalized_impulse(const Ravelin::VectorNd& gj);
    void apply_impulse(const Ravelin::SMomentumd& j, RigidBodyPtr link);
    void calc_spatial_inertias(RCArticulatedBodyPtr body);

    /// The body that this algorithm operates on
    boost::weak_ptr<RCArticulatedBody> _body;

    /// The spatial accelerations
    std::vector<Ravelin::SAcceld> _a;

    /// The articulated body inertias
    std::vector<Ravelin::SpatialABInertiad> _I;

    /// The articulated body spatial zero accelerations
    std::vector<Ravelin::SForced> _Z;

    /// Vector of link velocity updates
    std::vector<Ravelin::SVelocityd> _dv;

    /// The spatial coriolis vectors
    std::vector<Ravelin::SAcceld> _c;

    /// The expressions I*s
    std::vector<std::vector<Ravelin::SMomentumd> > _Is;

    /// Cholesky factorizations sIs
    std::vector<Ravelin::MatrixNd> _sIs;

    /// SVDs of sIs
    std::vector<Ravelin::MatrixNd> _usIs, _vsIs;
    std::vector<Ravelin::VectorNd> _ssIs;

    /// Determines whether the equations for a joint are rank deficient 
    std::vector<bool> _rank_deficient;

    /// The temporary expression Q - I*s'*c - s'*Z
    std::vector<Ravelin::VectorNd> _mu;

  private:
    // pointer to the linear algebra routines
    boost::shared_ptr<Ravelin::LinAlgd> _LA;

    /// work variables 
    Ravelin::VectorNd _workv, _workv2, _sTY, _qd_delta, _sIsmu, _Qi, _Q;
    Ravelin::MatrixNd _sIss, _workM;
    std::vector<Ravelin::SMomentumd> _Y;

    /// processed vector
    std::vector<bool> _processed;

    void calc_fwd_dyn_special();
    static double sgn(double x);
    static void push_children(RigidBodyPtr link, std::queue<RigidBodyPtr>& q);
    void apply_coulomb_joint_friction(RCArticulatedBodyPtr body);
    void apply_generalized_impulse(unsigned index, Ravelin::VectorNd& vgj);
    void set_spatial_velocities(RCArticulatedBodyPtr body);
    void calc_spatial_accelerations(RCArticulatedBodyPtr body);
    void calc_spatial_zero_accelerations(RCArticulatedBodyPtr body);
    void calc_spatial_coriolis_vectors(RCArticulatedBodyPtr body);
    Ravelin::VectorNd& solve_sIs(unsigned idx, const Ravelin::VectorNd& v, Ravelin::VectorNd& result) const;
    Ravelin::MatrixNd& solve_sIs(unsigned idx, const Ravelin::MatrixNd& v, Ravelin::MatrixNd& result) const;
    Ravelin::MatrixNd& transpose_solve_sIs(unsigned idx, const std::vector<Ravelin::SVelocityd>& m, Ravelin::MatrixNd& result) const;
}; // end class
} // end namespace

#endif

