/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RC_ARTICULATED_BODY_H
#define _RC_ARTICULATED_BODY_H

#include <pthread.h>
#include <map>
#include <list>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/SVector6.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/FSABAlgorithm.h>
#include <Moby/CRBAlgorithm.h>

namespace Moby {

class Joint;

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
    MatrixN& generalized_inertia_mult(const MatrixN& M, MatrixN& result);
    virtual MatrixN calc_jacobian_column(JointPtr joint, const Vector3& point);
    virtual MatrixN& calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const std::map<JointPtr, VectorN>& q, MatrixN& Jc);
    virtual MatrixN& calc_jacobian_column(JointPtr joint, const Vector3& point, MatrixN& Jc);
    virtual MatrixN calc_jacobian_floating_base(const Vector3& point);
    virtual void reset_accumulators();
    virtual void update_link_transforms();    
    virtual void update_link_velocities();
    virtual void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& p, RigidBodyPtr link);
    virtual void calc_fwd_dyn(Real dt);
    virtual void update_visualization();
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    RCArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<RCArticulatedBody>(shared_from_this()); }
    RCArticulatedBodyConstPtr get_this() const { return boost::dynamic_pointer_cast<const RCArticulatedBody>(shared_from_this()); }
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gj);
    virtual VectorN& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gc);
    virtual VectorN& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv);
    virtual VectorN& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gv);
    virtual MatrixN& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M);
    virtual MatrixN& get_generalized_inertia_inverse(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M);
    virtual VectorN& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, VectorN& f);
    virtual VectorN& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& f, const Vector3& t, VectorN& gf);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual void set_links(const std::vector<RigidBodyPtr>& links);
    virtual void set_joints(const std::vector<JointPtr>& links);
    virtual void update_event_data(EventProblemData& epd);
    virtual void update_velocity(const EventProblemData& epd);
    virtual void invalidate_positions();
    virtual void invalidate_velocities();
    virtual unsigned num_joint_dof_explicit() const;
    virtual unsigned num_joint_dof_implicit() const { return _n_joint_DOF_implicit; }
    void set_floating_base(bool flag);
    virtual VectorN& transpose_Jc_mult(const VectorN& v, VectorN& result) { return _Jc.transpose_mult(v, result); } 
    virtual MatrixN& transpose_Jc_mult(const MatrixN& m, MatrixN& result) { return _Jc.transpose_mult(m, result); }
    virtual VectorN& transpose_Dc_mult(const VectorN& v, VectorN& result) { return _Dc.transpose_mult(v, result); }
    virtual MatrixN& transpose_Dc_mult(const MatrixN& m, MatrixN& result) { return _Dc.transpose_mult(m, result); }
    virtual VectorN& transpose_Jl_mult(const VectorN& v, VectorN& result) { return _Jl.transpose_mult(v, result); }
    virtual MatrixN& transpose_Jl_mult(const MatrixN& m, MatrixN& result) { return _Jl.transpose_mult(m, result); }
    virtual VectorN& transpose_Dx_mult(const VectorN& v, VectorN& result) { return _Dx.transpose_mult(v, result); }
    virtual MatrixN& transpose_Dx_mult(const MatrixN& m, MatrixN& result) { return _Dx.transpose_mult(m, result); }

    /// Gets whether the base of this body is fixed or "floating"
    bool is_floating_base() const { return _floating_base; }

    /// Gets the number of DOF of the implicit joints in the body, not including floating base DOF
    virtual unsigned num_joint_dof() const { return _n_joint_DOF_implicit + num_joint_dof_explicit(); }

    /// Gets the base link
    RigidBodyPtr get_base_link() const { return _links.front(); }

   /// The computation reference frame
    ReferenceFrameType computation_frame_type;

    /// The forward dynamics algorithm
    ForwardDynamicsAlgorithmType algorithm_type;

    /// Gets constraint events (currently not any)
    virtual void get_constraint_events(std::vector<Event>& events) const { }

    /// Baumgarte alpha parameter >= 0
    Real b_alpha;

    /// Baumgarte beta parameter >= 0
    Real b_beta;

  protected:
    /// Whether this body uses a floating base
    bool _floating_base;
  
     virtual void compile();

    /// There is no visualization transform
    virtual const Matrix4* get_visualization_transform() { return NULL; }    

    /// The number of DOF of the implicit joint constraints in the body (does not include floating base DOF!)
    unsigned _n_joint_DOF_implicit;

    /// Gets the vector of implicit joint constraints
    const std::vector<JointPtr>& get_implicit_joints() const { return _ijoints; }

  private:
    bool all_children_processed(RigidBodyPtr link) const;
    void calc_fwd_dyn_loops();
    void calc_fwd_dyn_advanced_friction(Real dt);
    void generalized_inertia_mult_fixed(const MatrixN& M, MatrixN& result);
    void generalized_inertia_mult_floating(const MatrixN& M, MatrixN& result);
    virtual MatrixN calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const std::map<JointPtr, VectorN>& q);

    /// The vector of implicit joint constraints
    std::vector<JointPtr> _ijoints;

    /// The vector of explicit joint constraints
    std::vector<JointPtr> _ejoints;

    /// Variables used for events
    MatrixN _Jc, _Dc, _Jl, _Jx, _Dx, _Dt;
    MatrixN _iM_JcT, _iM_DcT, _iM_JlT, _iM_JxT, _iM_DxT, _iM_DtT;

    /// The CRB algorithm
    CRBAlgorithm _crb;

    /// The FSAB algorithm
    FSABAlgorithm _fsab;

    /// The inverse or factorized generalized inertia matrix
    MatrixN _invM;

    /// Indicates whether the inverse generalized inertia needs to be recomputed
    bool _invM_valid;

    /// Indicates whether the inverse generalized inertia is rank deficient
    bool _invM_rankdef;

    /// Indicates the type of the inverse generalized inertia
    DynamicBody::GeneralizedCoordinateType _invM_type;

    static Real sgn(Real x);
    bool treat_link_as_leaf(RigidBodyPtr link) const;
    void update_factorized_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype);
    virtual VectorN& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& v, VectorN& result);
    virtual MatrixN& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const MatrixN& m, MatrixN& result);
    MatrixN& solve_generalized_inertia_transpose(DynamicBody::GeneralizedCoordinateType gctype, const MatrixN& m, MatrixN& result);
    void determine_contact_jacobians(const EventProblemData& q, const VectorN& v, const MatrixN& M, MatrixN& Jc, MatrixN& Dc);
    static bool supports(JointPtr joint, RigidBodyPtr link);
    void determine_generalized_forces(VectorN& gf) const;
    void determine_generalized_accelerations(VectorN& xdd) const;
    void determine_constraint_force_transform(MatrixN& K) const;
    void M_mult(const VectorN& v, VectorN& result) const;
    void set_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& a);
    void determine_explicit_constraint_movement_jacobian(MatrixN& D);
    void determine_explicit_constraint_jacobians(const EventProblemData& q, MatrixN& Jx, MatrixN& Dx) const;
    void determine_explicit_constraint_jacobian(MatrixN& J);
    void determine_explicit_constraint_jacobian_dot(MatrixN& J) const;
    void set_explicit_constraint_forces(const VectorN& lambda);
}; // end class

} // end namespace
#endif

