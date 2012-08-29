/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _FS_AB_ALGORITHM_H
#define _FS_AB_ALGORITHM_H

#include <queue>
#include <Moby/SpatialABInertia.h>

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
  public:
    FSABAlgorithm();
    ~FSABAlgorithm() {}
    RCArticulatedBodyPtr get_body() const { return RCArticulatedBodyPtr(_body); }
    void set_body(RCArticulatedBodyPtr body) { _body = body; }
    void calc_fwd_dyn();
    void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gj);
    void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& contact_point, RigidBodyPtr link);
    void calc_inverse_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixNN& iM);
    void invalidate_position_data() { _position_data_valid = false; }
    void invalidate_velocity_data() { _velocity_data_valid = false; }

    /// The body that this algorithm operates on
    boost::weak_ptr<RCArticulatedBody> _body;

    /// The reference frame all quantities are computed in
    ReferenceFrameType _rftype;

    /// The spatial velocities
    std::vector<SVector6> _v;

    /// The spatial accelerations
    std::vector<SVector6> _a;

    /// The spatial isolated inertias
    std::vector<SpatialRBInertia> _Iiso;

    /// The articulated body inertias
    std::vector<SpatialABInertia> _I;

    /// The articulated body spatial zero accelerations
    std::vector<SVector6> _Z;

    /// Vector of link velocity updates
    std::vector<SVector6> _dv;

    /// The spatial coriolis vectors
    std::vector<SVector6> _c;

    /// The expression I*s
    std::vector<SMatrix6N> _Is;

    /// The temporary factorizations inv(sIs)
    std::vector<MatrixNN> _sIs;

    /// Determines whether the equations for a joint are rank deficient 
    std::vector<bool> _rank_deficient;

    /// The temporary expression Q - I*s'*c - s'*Z
    std::vector<VectorN> _mu;

  private:
    /// Is positional data valid?
    bool _position_data_valid;

    /// Is velocity data valid?
    bool _velocity_data_valid;

    static Real sgn(Real x);
    static void push_children(RigidBodyPtr link, std::queue<RigidBodyPtr>& q);
    void apply_coulomb_joint_friction(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_impulse_dyn(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void apply_generalized_impulse(unsigned index, const std::vector<MatrixN>& sTI, VectorN& vgj) const;
    void set_spatial_iso_inertias(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void set_spatial_velocities(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_spatial_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_spatial_zero_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_spatial_inertias(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_spatial_coriolis_vectors(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    VectorN& solve_sIs(unsigned idx, const VectorN& v, VectorN& result) const;
    MatrixN& solve_sIs(unsigned idx, const MatrixN& v, MatrixN& result) const;
    MatrixN& transpose_solve_sIs(unsigned idx, const SMatrix6N& m, MatrixN& result) const;
}; // end class
} // end namespace

#endif

