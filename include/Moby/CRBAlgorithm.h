/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CRB_ALGORITHM_H
#define _CRB_ALGORITHM_H

#include <Moby/SpatialRBInertia.h>
#include <Moby/MatrixNN.h>

namespace Moby {

/// Computes forward dynamics using composite-rigid body method
class CRBAlgorithm
{
  public:
    CRBAlgorithm();
    ~CRBAlgorithm() {}
    RCArticulatedBodyPtr get_body() const { return RCArticulatedBodyPtr(_body); }
    void set_body(RCArticulatedBodyPtr body) { _body = body; setup_parent_array(); }
    void calc_fwd_dyn();
    void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& contact_point, RigidBodyPtr link);
    void calc_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixNN& M);
    void calc_generalized_forces(SVector6& f0, VectorN& C);
    void invalidate_position_data() { _position_data_valid = false; }
    void invalidate_velocity_data() { _velocity_data_valid = false; }
    bool factorize_cholesky(MatrixNN& M);
    VectorN& M_solve(const VectorN& v, VectorN& result);
    MatrixN& M_solve(const MatrixN& v, MatrixN& result);

  private:
    std::vector<unsigned> _lambda;
    void setup_parent_array();

    /// Is positional data valid?
    bool _position_data_valid;

    /// Is velocity data valid?
    bool _velocity_data_valid;

    /// The body that this algorithm operates on
    boost::weak_ptr<RCArticulatedBody> _body;

    /// The spatial acceleration of the base computed on the last call to calc_fwd_dyn()
    SVector6 _a0;

    /// The vector of joint accelerations computed on the last call to calc_fwd_dyn()
    VectorN _qdd;

    /// The joint space inertia matrix H (fixed base) or augmented matrix [I_0^c K; K^s H] (floating base, see [Featherstone 1987], p. 123) used to compute forward dynamics for floating bases
    MatrixNN _M;

    /// A factorization (or possibly inverse) of the matrix M; note that we compute this b/c we generally may need to solve multiple systems of linear equations using this matrix as a LHS at different times -- always in global frame
    MatrixNN _fM;

    /// Determines whether the system of equations for forward dynamics is rank-deficient
     bool _rank_deficient;

    /// The vector of spatial velocities determined in last call to calc_fwd_dyn()
     std::vector<SVector6> _velocities;

    /// The last reference frame used for computation 
     ReferenceFrameType _rftype;

    void calc_generalized_inertia_axisangle(MatrixNN& M) const;
    void calc_generalized_inertia_rodrigues(MatrixNN& M) const;
    void calc_joint_space_inertia(RCArticulatedBodyPtr body, ReferenceFrameType rftype, MatrixNN& H, std::vector<SpatialRBInertia>& Ic) const;
    void apply_coulomb_joint_friction(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void precalc(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void apply_impulse_fixed_base(RCArticulatedBodyPtr body, const Vector3& j, const Vector3& k, const Vector3& point, RigidBodyPtr link);
    void apply_impulse_floating_base(RCArticulatedBodyPtr body, const Vector3& j, const Vector3& k, const Vector3& point, RigidBodyPtr link);
    void set_spatial_velocities(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_generalized_inertia(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_fwd_dyn_fixed_base(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void calc_fwd_dyn_floating_base(RCArticulatedBodyPtr body, ReferenceFrameType rftype);
    void update_link_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype) const;
    static void to_spatial7_inertia(const SpatialRBInertia& I, const Quat& q, MatrixNN& I7);
    VectorN& M_solve_noprecalc(const VectorN& v, VectorN& result) const;
    MatrixN& M_solve_noprecalc(const MatrixN& v, MatrixN& result) const;
};
}

#endif
