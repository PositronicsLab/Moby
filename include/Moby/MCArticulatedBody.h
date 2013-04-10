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
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Moby/Constants.h>
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
    virtual Ravelin::MatrixNd& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::MatrixNd& M);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gj);
    virtual Ravelin::VectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv);
    virtual Ravelin::VectorNd& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& ga); 
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_velocities(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv) { return get_generalized_velocity(gctype, gv); }
    virtual Ravelin::VectorNd& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& Qf);
    virtual void reset_accumulators();
    virtual void apply_impulse(const Ravelin::Wrenchd& j, RigidBodyPtr link);
    virtual void calc_fwd_dyn(double dt);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    MCArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<MCArticulatedBody>(shared_from_this()); }
    MCArticulatedBodyConstPtr get_this() const { return boost::dynamic_pointer_cast<const MCArticulatedBody>(shared_from_this()); }
    virtual Ravelin::VectorNd& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr link, const Ravelin::Wrenchd& f, Ravelin::VectorNd& gf);
    virtual void update_event_data(EventProblemData& epd);
    virtual void update_velocity(const EventProblemData& epd);
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);
    virtual unsigned num_joint_dof_explicit() const { return num_joint_dof(); }
    virtual unsigned num_joint_dof_implicit() const { return 0; }
    virtual Ravelin::VectorNd& transpose_Jc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return mult_transpose_sparse(_Jc, v, result); } 
    virtual Ravelin::MatrixNd& transpose_Jc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return mult_transpose_sparse(_Jc, m, result); }
    virtual Ravelin::VectorNd& transpose_Dc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return mult_transpose_sparse(_Dc, v, result); }
    virtual Ravelin::MatrixNd& transpose_Dc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return mult_transpose_sparse(_Dc, m, result); }
    virtual Ravelin::VectorNd& transpose_Jl_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return mult_transpose_sparse(_Jl, v, result); }
    virtual Ravelin::MatrixNd& transpose_Jl_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return mult_transpose_sparse(_Jl, m, result); }
    virtual Ravelin::VectorNd& transpose_Dx_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return mult_transpose_sparse(_Dx, v, result); }
    virtual Ravelin::MatrixNd& transpose_Dx_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return mult_transpose_sparse(_Dx, m, result); }

    /// The Baumgarte stabilization constant alpha >= 0
    double b_alpha;

    /// The Baumgarte stabilization constant beta >= 0
    double b_beta;

  protected:
    virtual void compile();

    /// There is no visualization transform (returns NULL)
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_visualization_transform() { return boost::shared_ptr<const Ravelin::Pose3d>(); }

  private:
    enum JacobianType { eNone, eContactNormal, eContactTangent, eLimit,
                        eJointFriction, eConstraint };

    class SparseJacobian : public Ravelin::MatrixNd
    {
      public:
        std::vector<std::vector<unsigned> > indices;
    };

    /// An inverse inertia of a rigid body matrix with frame referenced at the body's c.o.m. 
    struct InvInertia
    {
      double inv_mass;      // The inverse of the rigid body mass
      Ravelin::Matrix3d inv_inertia;  // The inverse of the rigid body inertia 
    };

    static bool affects(RigidBodyPtr rb, Event* e);
    static unsigned num_sub_events(JacobianType jt, Event* e);
    static void get_event_data(JacobianType jt, Event* e, RigidBodyPtr rb, unsigned subidx, Ravelin::Vector3d& tx, Ravelin::Vector3d& rx);
    static const std::vector<Event*>& get_events_vector(const EventProblemData& q, JacobianType jt);
    void update_Jx_iM_JyT(EventProblemData& q, Ravelin::MatrixNd& Jx_iM_JyT, JacobianType j1t, JacobianType j2t);
    void update_Jx_iM_JyT(RigidBodyPtr rb, EventProblemData& q, Ravelin::MatrixNd& Jx_iM_JyT, JacobianType j1t, JacobianType j2t);
    void calc_joint_accelerations();
    static double sgn(double x); 
    Ravelin::MatrixNd dense_J(const SparseJacobian& J) const;
    Ravelin::MatrixNd& dense_J(const SparseJacobian& J, Ravelin::MatrixNd& dJ) const;
    void precalc();
    void get_constraint_jacobian(SparseJacobian& J) const;
    void get_constraint_jacobian_dot(SparseJacobian& J) const;
    void get_constraint_jacobian_numerically(SparseJacobian& J) const;
    void get_mechanism_jacobian(SparseJacobian& J, SparseJacobian& J_dot) const;
    Ravelin::VectorNd& get_constraint_evals(Ravelin::VectorNd& C) const;
    Ravelin::VectorNd& iM_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) const;
    static void transform(RigidBodyPtr rb, double CJrb[7]);
    void form_Jm_iM_Km(const std::vector<unsigned>& Jm_indices, const std::vector<unsigned>& Km_indices, Ravelin::MatrixNd& M);
    Ravelin::VectorNd& scale_inverse_inertia(unsigned i, Ravelin::VectorNd& v) const;
    void apply_joint_limit_impulses();
    void update_Jx_v(EventProblemData& q);
    void update_Jl_v(EventProblemData& q);
    void update_Dx_v(EventProblemData& q);
    void determine_inertias();
    void update_link_accelerations() const;
    void update_link_velocities() const;
    Ravelin::VectorNd& mult_transpose_sparse(const SparseJacobian& J, const Ravelin::VectorNd& v, Ravelin::VectorNd& result) const;
    Ravelin::MatrixNd& mult_transpose_sparse(const SparseJacobian& J, const Ravelin::MatrixNd& v, Ravelin::MatrixNd& result) const;
    Ravelin::VectorNd& mult_sparse(const SparseJacobian& J, const Ravelin::VectorNd& v, Ravelin::VectorNd& result) const;
    void calc_Dx_iM_DxT(Ravelin::MatrixNd& Dx_iM_DxT) const;
    void calc_Dx_iM(SparseJacobian& Dx_iM) const;
    Ravelin::MatrixNd& calc_Jx_iM_JyT(const SparseJacobian& Jx, const SparseJacobian& Jy, Ravelin::MatrixNd& Jx_iM_JyT) const;
    static void get_sub_jacobian(const std::vector<unsigned>& rows, const SparseJacobian& J, SparseJacobian& Jx);
    static void increment_dof(RigidBodyPtr rb1, RigidBodyPtr rb2, unsigned k, double h);
    virtual Ravelin::VectorNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
    virtual Ravelin::MatrixNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    void select_sub_contact_Jacobians(const EventProblemData& q, SparseJacobian& Jc_sub, SparseJacobian& Dc_sub) const;
    Ravelin::VectorNd& solve_Jx_iM_JxT(const Ravelin::VectorNd& rhs, Ravelin::VectorNd& x) const;

    /// The last-computed generalized velocity (in axis-angle representation)
    Ravelin::VectorNd _xd;

    /// The acceleration vector (generalized coordinates in axis-angle representation) calculated by calc_fwd_dyn()
    Ravelin::VectorNd _xdd;

    /// The block inverse inertia matrix
    std::vector<InvInertia> _iM;

    /// The matrix J*iM*J' 
    Ravelin::MatrixNd _Jx_iM_JxT;

    /// The factorization or regularized inverse of Jx*iM*Jx'
    Ravelin::MatrixNd _inv_Jx_iM_JxT;

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

