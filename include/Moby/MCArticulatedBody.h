/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MC_ARTICULATED_BODY_H
#define _MC_ARTICULATED_BODY_H

#include <pthread.h>
#include <map>
#include <list>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/SpatialTransform.h>
#include <Moby/ArticulatedBody.h>

namespace Moby {

class RigidBody;

/// Defines an articulated body simulated using maximal coordinate methods
/**
 * Maximal-coordinate methods are generally inferior to reduced-coordinate
 * methods for simulating loop-free structures, due to increased computation
 * requirements as well as numerical and physical stability issues that
 * arise in correcting drift.  Regardless, a maximal-coordinate formulation
 * is often the best option for structures with kinematic loops.  Moby
 * does not currently provide means to simulate structures with kinematic
 * loops using reduced-coordinate methods, so a maximal-coordinate approach
 * is necessary.  For loop-free structures, both reduced and 
 * maximal-coordinate methods are feasible.
 */
class MCArticulatedBody : public ArticulatedBody
{
  public:
    MCArticulatedBody();
    virtual ~MCArticulatedBody() {}
    virtual MatrixNN& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixNN& M);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gj);
    virtual VectorN& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv);
    virtual VectorN& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, VectorN& ga); 
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gv);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gc);
    virtual VectorN& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gc);
    virtual VectorN& get_generalized_velocities(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv) { return get_generalized_velocity(gctype, gv); }
    virtual VectorN& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, VectorN& Qf);
    virtual void reset_accumulators();
    virtual void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& contact_point, RigidBodyPtr link);
    virtual void calc_fwd_dyn(Real dt);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    MCArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<MCArticulatedBody>(shared_from_this()); }
    MCArticulatedBodyConstPtr get_this() const { return boost::dynamic_pointer_cast<const MCArticulatedBody>(shared_from_this()); }
    virtual VectorN& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr link, const Vector3& f, const Vector3& t, VectorN& gf);
    virtual void update_event_data(EventProblemData& epd);
    virtual void update_velocity(const EventProblemData& epd);
    virtual void integrate(Real t, Real h, boost::shared_ptr<Integrator<VectorN> > integrator);
    virtual unsigned num_joint_dof_explicit() const { return num_joint_dof(); }
    virtual unsigned num_joint_dof_implicit() const { return 0; }
    virtual VectorN& transpose_Jc_mult(const VectorN& v, VectorN& result) { return mult_transpose_sparse(_Jc, v, result); } 
    virtual MatrixN& transpose_Jc_mult(const MatrixN& m, MatrixN& result) { return mult_transpose_sparse(_Jc, m, result); }
    virtual VectorN& transpose_Dc_mult(const VectorN& v, VectorN& result) { return mult_transpose_sparse(_Dc, v, result); }
    virtual MatrixN& transpose_Dc_mult(const MatrixN& m, MatrixN& result) { return mult_transpose_sparse(_Dc, m, result); }
    virtual VectorN& transpose_Jl_mult(const VectorN& v, VectorN& result) { return mult_transpose_sparse(_Jl, v, result); }
    virtual MatrixN& transpose_Jl_mult(const MatrixN& m, MatrixN& result) { return mult_transpose_sparse(_Jl, m, result); }
    virtual VectorN& transpose_Dx_mult(const VectorN& v, VectorN& result) { return mult_transpose_sparse(_Dx, v, result); }
    virtual MatrixN& transpose_Dx_mult(const MatrixN& m, MatrixN& result) { return mult_transpose_sparse(_Dx, m, result); }

    /// The Baumgarte stabilization constant alpha >= 0
    Real b_alpha;

    /// The Baumgarte stabilization constant beta >= 0
    Real b_beta;

  protected:
    virtual void compile();

    /// There is no visualization transform (returns NULL)
    virtual const Matrix4* get_visualization_transform() { return NULL; }

  private:
    enum JacobianType { eNone, eContactNormal, eContactTangent, eLimit,
                        eJointFriction, eConstraint };

    class SparseJacobian : public MatrixN
    {
      public:
        std::vector<std::vector<unsigned> > indices;
    };

    /// An inverse inertia of a rigid body matrix with frame referenced at the body's c.o.m. 
    struct InvInertia
    {
      Real inv_mass;      // The inverse of the rigid body mass
      Matrix3 inv_inertia;  // The inverse of the rigid body inertia 
    };

    static MatrixN& reverse_transform(const SpatialTransform& X, const MatrixN& pinv_s, MatrixN& sx);
    static SpatialTransform calc_special_spatial_transform(const SpatialTransform& X);
    static bool affects(RigidBodyPtr rb, Event* e);
    static unsigned num_sub_events(JacobianType jt, Event* e);
    static void get_event_data(JacobianType jt, Event* e, RigidBodyPtr rb, unsigned subidx, Vector3& tx, Vector3& rx);
    static const std::vector<Event*>& get_events_vector(const EventProblemData& q, JacobianType jt);
    void update_Jx_iM_JyT(EventProblemData& q, MatrixN& Jx_iM_JyT, JacobianType j1t, JacobianType j2t);
    void update_Jx_iM_JyT(RigidBodyPtr rb, EventProblemData& q, MatrixN& Jx_iM_JyT, JacobianType j1t, JacobianType j2t);
    void calc_joint_accelerations();
    static Real sgn(Real x); 
    MatrixN dense_J(const SparseJacobian& J) const;
    MatrixN& dense_J(const SparseJacobian& J, MatrixN& dJ) const;
    void precalc();
    void get_constraint_jacobian(SparseJacobian& J) const;
    void get_constraint_jacobian_dot(SparseJacobian& J) const;
    void get_constraint_jacobian_numerically(SparseJacobian& J) const;
    void get_mechanism_jacobian(SparseJacobian& J, SparseJacobian& J_dot) const;
    VectorN& get_constraint_evals(VectorN& C) const;
    VectorN& iM_mult(const VectorN& v, VectorN& result) const;
    static void transform(RigidBodyPtr rb, Real CJrb[7]);
    void form_Jm_iM_Km(const std::vector<unsigned>& Jm_indices, const std::vector<unsigned>& Km_indices, MatrixN& M);
    VectorN& scale_inverse_inertia(unsigned i, VectorN& v) const;
    void apply_joint_limit_impulses();
    void update_Jx_v(EventProblemData& q);
    void update_Jl_v(EventProblemData& q);
    void update_Dx_v(EventProblemData& q);
    void determine_inertias();
    void update_link_accelerations() const;
    void update_link_velocities() const;
    VectorN& mult_transpose_sparse(const SparseJacobian& J, const VectorN& v, VectorN& result) const;
    MatrixN& mult_transpose_sparse(const SparseJacobian& J, const MatrixN& v, MatrixN& result) const;
    VectorN& mult_sparse(const SparseJacobian& J, const VectorN& v, VectorN& result) const;
    void calc_Dx_iM_DxT(MatrixNN& Dx_iM_DxT) const;
    void calc_Dx_iM(SparseJacobian& Dx_iM) const;
    MatrixN& calc_Jx_iM_JyT(const SparseJacobian& Jx, const SparseJacobian& Jy, MatrixN& Jx_iM_JyT) const;
    static void get_sub_jacobian(const std::vector<unsigned>& rows, const SparseJacobian& J, SparseJacobian& Jx);
    static void increment_dof(RigidBodyPtr rb1, RigidBodyPtr rb2, unsigned k, Real h);
    VectorN& solve_Jx_iM_JxT(const VectorN& rhs, VectorN& x) const;

    /// The last-computed generalized velocity (in axis-angle representation)
    VectorN _xd;

    /// The acceleration vector (generalized coordinates in axis-angle representation) calculated by calc_fwd_dyn()
    VectorN _xdd;

    /// The block inverse inertia matrix
    std::vector<InvInertia> _iM;

    /// The matrix J*iM*J' 
    MatrixNN _Jx_iM_JxT;

    /// The factorization or regularized inverse of Jx*iM*Jx'
    MatrixNN _inv_Jx_iM_JxT;

    /// Indicates whether J*iM*J' is rank deficient
    bool _rank_def;

    /// Constraint Jacobian
    SparseJacobian _Jx;

    /// Time derivative of the constraint Jacobian
    SparseJacobian _Jx_dot;

    /// The joint limit Jacobian
    SparseJacobian _Jl;

    /// The mechanism Jacobian (used for joint friction)
    SparseJacobian _Dx;

    /// The contact normal Jacobian
    SparseJacobian _Jc;

    /// The contact transpose Jacobian
    SparseJacobian _Dc;

    /// The time derivative of the mechanism Jacobian
    SparseJacobian _Dx_dot;

    /// The gc indices: maps link indices to starting coordinate indices
    std::vector<unsigned> _gc_indices;
}; // end class

} // end namespace
#endif

