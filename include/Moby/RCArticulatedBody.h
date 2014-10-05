/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RC_ARTICULATED_BODY_H
#define _RC_ARTICULATED_BODY_H

#include <pthread.h>
#include <map>
#include <list>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SForced.h>
#include <Moby/Constants.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/FSABAlgorithm.h>
#include <Moby/CRBAlgorithm.h>

namespace Moby {

class Joint;
class UnilateralConstraintProblemData;

/// Defines an articulated body for use with reduced-coordinate dynamics algorithms
/**
 * Reduced-coordinate articulated bodies cannot rely upon the integrator to automatically update
 * the states (i.e., positions, velocities) of the links, as is done with maximal-coordinate 
 * articulated bodies.  Rather, the integrator updates the joint positions and velocities; the
 * states are obtained from this reduced-coordinate representation.
 * Notes about concurrency: <br /><br />
 *
 * It is generally desirable to be able to run forward dynamics and inverse 
 * dynamics algorithms concurrently to simulate actual robotic systems.   In 
 * general, derived classes should not operate on state variables
 * (joint positions, velocities, accelerations and floating base positions, 
 * velocites, and accelerations) directly during execution of the algorithm.  
 * Rather, derived classes should operate on copies of the state
 * variables, updating the state variables on conclusion of the algorithms.  
 */
class RCArticulatedBody : public ArticulatedBody
{
  friend class CRBAlgorithm;
  friend class FSABAlgorithm;

  public:
    enum ForwardDynamicsAlgorithmType { eFeatherstone, eCRB }; 
    RCArticulatedBody();
    virtual ~RCArticulatedBody() {}
/*
    virtual Ravelin::MatrixNd& calc_jacobian(const Point3d& point, RigidBodyPtr link, Ravelin::MatrixNd& J);
    virtual Ravelin::MatrixNd& calc_jacobian(const Point3d& point, const Ravelin::Pose3d& base_pose, const std::map<JointPtr, Ravelin::VectorNd>& q, RigidBodyPtr link, Ravelin::MatrixNd& J);
*/
    virtual void reset_accumulators();
    virtual void update_link_poses();    
    virtual void update_link_velocities();
    virtual void apply_impulse(const Ravelin::SMomentumd& w, RigidBodyPtr link);
    virtual void calc_fwd_dyn();
    virtual void update_visualization();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    RCArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<RCArticulatedBody>(shared_from_this()); }
    boost::shared_ptr<const RCArticulatedBody> get_this() const { return boost::dynamic_pointer_cast<const RCArticulatedBody>(shared_from_this()); }
    virtual void set_generalized_forces(const Ravelin::SharedVectorNd& gf);
    virtual void set_generalized_forces(const Ravelin::VectorNd& gf) { DynamicBody::set_generalized_forces(gf); }
    virtual void add_generalized_force(const Ravelin::SharedVectorNd& gf);
    virtual void add_generalized_force(const Ravelin::VectorNd& gf) { DynamicBody::add_generalized_force(gf); }
    virtual void apply_generalized_impulse(const Ravelin::SharedVectorNd& gj);
    virtual void apply_generalized_impulse(const Ravelin::VectorNd& gj) { DynamicBody::apply_generalized_impulse(gj); }
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gc);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc) { DynamicBody::set_generalized_coordinates(gctype, gc); }
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gv);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv) { DynamicBody::set_generalized_velocity(gctype, gv); }
    virtual Ravelin::SharedMatrixNd& get_generalized_inertia(Ravelin::SharedMatrixNd& M);
    virtual Ravelin::MatrixNd& get_generalized_inertia(Ravelin::MatrixNd& M) { return DynamicBody::get_generalized_inertia(M); }
    virtual Ravelin::SharedVectorNd& get_generalized_forces(Ravelin::SharedVectorNd& f);
    virtual Ravelin::VectorNd& get_generalized_forces(Ravelin::VectorNd& f) { return DynamicBody::get_generalized_forces(f); }
    virtual Ravelin::SharedVectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::SharedVectorNd& gf);
    virtual Ravelin::VectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::VectorNd& gf) { return DynamicBody::convert_to_generalized_force(body, w, gf); }
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual void set_links_and_joints(const std::vector<RigidBodyPtr>& links, const std::vector<JointPtr>& joints);
    virtual unsigned num_joint_dof_implicit() const;
    virtual unsigned num_joint_dof_explicit() const { return _n_joint_DOF_explicit; }
    void set_floating_base(bool flag);
    virtual Ravelin::VectorNd& transpose_Jc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return _Jc.transpose_mult(v, result); } 
    virtual Ravelin::MatrixNd& transpose_Jc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return _Jc.transpose_mult(m, result); }
    virtual Ravelin::VectorNd& transpose_Dc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return _Dc.transpose_mult(v, result); }
    virtual Ravelin::MatrixNd& transpose_Dc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return _Dc.transpose_mult(m, result); }
    virtual Ravelin::VectorNd& transpose_Jl_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return _Jl.transpose_mult(v, result); }
    virtual Ravelin::MatrixNd& transpose_Jl_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return _Jl.transpose_mult(m, result); }
    virtual Ravelin::VectorNd& transpose_Dx_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) { return _Dx.transpose_mult(v, result); }
    virtual Ravelin::MatrixNd& transpose_Dx_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) { return _Dx.transpose_mult(m, result); }
    virtual void set_computation_frame_type(ReferenceFrameType rftype);
    virtual Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::VectorNd& x) { return DynamicBody::solve_generalized_inertia(b, x); }
    virtual Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);
    virtual Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& v, Ravelin::SharedVectorNd& result);
    virtual Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& m, Ravelin::SharedMatrixNd& result);
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_gc_pose() const; 
    virtual void validate_position_variables();
    virtual Ravelin::SharedVectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc) { return DynamicBody::get_generalized_coordinates(gctype, gc); }
    virtual Ravelin::SharedVectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gv);
    virtual Ravelin::VectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv) { return DynamicBody::get_generalized_velocity(gctype, gv); }
    virtual Ravelin::SharedVectorNd& get_generalized_acceleration(Ravelin::SharedVectorNd& ga);
    virtual Ravelin::VectorNd& get_generalized_acceleration(Ravelin::VectorNd& ga) { return DynamicBody::get_generalized_acceleration(ga); }

    template <class V>
    void get_generalized_acceleration_generic(V& ga);

    template <class V>
    void get_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv);

    template <class V>
    void get_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv);

    /// Gets whether the base of this body is fixed or "floating"
    virtual bool is_floating_base() const { return _floating_base; }

    /// Gets the number of DOF of the explicit joints in the body, not including floating base DOF
    virtual unsigned num_joint_dof() const { return _n_joint_DOF_explicit + num_joint_dof_implicit(); }

    /// Gets the base link
    virtual RigidBodyPtr get_base_link() const { return (!_links.empty()) ? _links.front() : RigidBodyPtr(); }

    /// The forward dynamics algorithm
    ForwardDynamicsAlgorithmType algorithm_type;

    /// Baumgarte alpha parameter >= 0
    double b_alpha;

    /// Baumgarte beta parameter >= 0
    double b_beta;

  protected:
    /// Whether this body uses a floating base
    bool _floating_base;
  
     virtual void compile();

    /// The number of DOF of the explicit joint constraints in the body (does not include floating base DOF!)
    unsigned _n_joint_DOF_explicit;

    /// Gets the vector of explicit joint constraints
    const std::vector<JointPtr>& get_explicit_joints() const { return _ejoints; }

  private:
    virtual Ravelin::MatrixNd& calc_jacobian_column(JointPtr joint, const Point3d& point, Ravelin::MatrixNd& Jc);
/*
    virtual Ravelin::MatrixNd& calc_jacobian_floating_base(const Point3d& point, Ravelin::MatrixNd& J);
*/
    bool all_children_processed(RigidBodyPtr link) const;
    void calc_fwd_dyn_loops();
    void calc_fwd_dyn_advanced_friction(double dt);

    /// The vector of explicit joint constraints
    std::vector<JointPtr> _ejoints;

    /// The vector of implicit joint constraints
    std::vector<JointPtr> _ijoints;

    /// Variables used for events
    Ravelin::MatrixNd _Jc, _Dc, _Jl, _Jx, _Dx, _Dt;
    Ravelin::MatrixNd _iM_JcT, _iM_DcT, _iM_JlT, _iM_JxT, _iM_DxT, _iM_DtT;

    /// Indicates when position data has been invalidated
    bool _position_invalidated;

    /// The CRB algorithm
    CRBAlgorithm _crb;

    /// The FSAB algorithm
    FSABAlgorithm _fsab;

    /// Linear algebra object
    boost::shared_ptr<Ravelin::LinAlgd> _LA;

    static double sgn(double x);
    bool treat_link_as_leaf(RigidBodyPtr link) const;
    void update_factorized_generalized_inertia();
    void determine_contact_jacobians(const UnilateralConstraintProblemData& q, const Ravelin::VectorNd& v, const Ravelin::MatrixNd& M, Ravelin::MatrixNd& Jc, Ravelin::MatrixNd& Dc);
    static bool supports(JointPtr joint, RigidBodyPtr link);
    void determine_generalized_forces(Ravelin::VectorNd& gf) const;
    void determine_generalized_accelerations(Ravelin::VectorNd& xdd) const;
    void determine_constraint_force_transform(Ravelin::MatrixNd& K) const;
    void set_generalized_acceleration(const Ravelin::VectorNd& a);
    void determine_implicit_constraint_movement_jacobian(Ravelin::MatrixNd& D);
    void determine_implicit_constraint_jacobians(const UnilateralConstraintProblemData& q, Ravelin::MatrixNd& Jx, Ravelin::MatrixNd& Dx) const;
    void determine_implicit_constraint_jacobian(Ravelin::MatrixNd& J);
    void determine_implicit_constraint_jacobian_dot(Ravelin::MatrixNd& J) const;
    void set_implicit_constraint_forces(const Ravelin::VectorNd& lambda);
}; // end class

#include "RCArticulatedBody.inl"

} // end namespace
#endif

