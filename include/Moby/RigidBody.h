/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RIGID_BODY_H
#define _RIGID_BODY_H

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/SForced.h>
#include <Ravelin/SVelocityd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>

namespace osg { class Node; }

namespace Moby {

class ArticulatedBody;
class CollisionGeometry;

/// Represents a single rigid body
/**
 *  Contains information needed to represent a rigid body, including position
 *  and velocity (both linear and angular), mass, center of
 *  mass, inertia matrix, collision data, and visualization data.  This class
 *  is used for both non-articulated and articulated rigid bodies, though not
 *  all member data may be used in the latter.
 * \todo implement rest matrix
 */
class RigidBody : public SingleBody
{
  friend class ArticulatedBody;
  friend class RCArticulatedBody;
  friend class MCArticulatedBody;
  friend class Joint;
  
  public:
    enum Compliance { eRigid, eCompliant};
    RigidBody();
    virtual ~RigidBody() {}
    void add_force(const Ravelin::SForced& w);
    void set_pose(const Ravelin::Pose3d& pose);
    void set_inertial_pose(const Ravelin::Pose3d& pose);
    void set_enabled(bool flag);
    void apply_impulse(const Ravelin::SMomentumd& w);
    virtual void rotate(const Ravelin::Quatd& q);
    virtual void translate(const Ravelin::Origin3d& o);
    virtual void calc_fwd_dyn();
    const Ravelin::SpatialRBInertiad& get_inertia();
    void set_inertia(const Ravelin::SpatialRBInertiad& J);
    boost::shared_ptr<const Ravelin::Pose3d> get_inertial_pose() const { return _jF; }

    virtual void set_visualization_data(osg::Node* vdata) { Visualizable::set_visualization_data(vdata); }

    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    bool is_child_link(boost::shared_ptr<const RigidBody> query) const;
    bool is_descendant_link(boost::shared_ptr<const RigidBody> query) const;
    const Ravelin::SVelocityd& get_velocity();
    void set_velocity(const Ravelin::SVelocityd& xd);
    void set_accel(const Ravelin::SAcceld& xdd);
    virtual const Ravelin::SAcceld& get_accel();
    void set_velocity(const Ravelin::SAcceld& xdd);
    virtual void set_generalized_forces(const Ravelin::SharedVectorNd& gf);
    virtual Ravelin::SharedMatrixNd& get_generalized_inertia(Ravelin::SharedMatrixNd& M);
    virtual Ravelin::SharedVectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gc);
    virtual Ravelin::SharedVectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gv);
    virtual Ravelin::SharedVectorNd& get_generalized_acceleration(Ravelin::SharedVectorNd& ga);
    virtual void add_generalized_force(const Ravelin::SharedVectorNd& gf);
    virtual void apply_generalized_impulse(const Ravelin::SharedVectorNd& gf);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gv);
    virtual Ravelin::SharedVectorNd& get_generalized_forces(Ravelin::SharedVectorNd& f);
    virtual Ravelin::SharedVectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::SharedVectorNd& gf);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);
    Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia_single(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);
    virtual Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);
    virtual Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::SharedVectorNd& x);
    RigidBodyPtr get_parent_link() const;
    JointPtr get_inner_joint_explicit() const;
    void add_inner_joint(JointPtr j);
    void add_outer_joint(JointPtr j);
    void remove_inner_joint(JointPtr joint);
    void remove_outer_joint(JointPtr joint);
    virtual double calc_kinetic_energy();
    virtual Ravelin::Vector3d calc_point_vel(const Point3d& p) const;
    bool is_base() const;
    bool is_ground() const;
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_computation_frame() const;
    virtual void set_computation_frame_type(ReferenceFrameType rftype);
    virtual Ravelin::MatrixNd& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J);
    virtual Ravelin::MatrixNd& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J);
    const Ravelin::SForced& sum_forces();
    void reset_accumulators();
    Ravelin::SForced calc_euler_torques();
    virtual void ode_noexcept(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) { prepare_to_calc_ode(x, t, dt, data); }
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void reset_limit_estimates();
    virtual bool limit_estimates_exceeded() const { return compliance == eRigid && _vel_limit_exceeded; }
    const Ravelin::SVelocityd& get_vel_upper_bounds() const { return _vel_limit_hi; }
    const Ravelin::SVelocityd& get_vel_lower_bounds() const { return _vel_limit_lo; }
    void update_vel_limits();

    template <class OutputIterator>
    OutputIterator get_parent_links(OutputIterator begin) const;

    template <class OutputIterator>
    OutputIterator get_child_links(OutputIterator begin) const;

    /// Collision geometries, if any, for this rigid body
    std::list<CollisionGeometryPtr> geometries;

    /// Gets the shared pointer for <b>this</b>
    RigidBodyPtr get_this() { return boost::dynamic_pointer_cast<RigidBody>(shared_from_this()); }
  
    /// Gets the shared const pointer for <b>this</b>
    boost::shared_ptr<const RigidBody> get_this() const { return boost::dynamic_pointer_cast<const RigidBody>(shared_from_this()); }

    /// Gets the current pose of this body
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _F; }

    // Gets the pose used in generalized coordinates calculations
    boost::shared_ptr<const Ravelin::Pose3d> get_gc_pose() const { return _F2; }


    // Gets the "mixed" pose of this body (pose origin at the body's reference point but pose aligned with global frame)
    boost::shared_ptr<const Ravelin::Pose3d> get_mixed_pose() const { return _F2; }

    /// Synonym for get_mass() (implements SingleBody::calc_mass())
    double calc_mass() const { return _Jm.m; }

    /// Gets the mass of this body
    virtual double get_mass() const { return _Jm.m; }
    
    /// Gets whether this body is enabled
    bool is_enabled() const { return _enabled; }

    /// Gets the articulated body corresponding to this body
    /**
     * \return a pointer to the articulated body, or NULL if this body is not 
     *         a link an articulated body
     */
    boost::shared_ptr<ArticulatedBody> get_articulated_body() const { return (_abody.expired()) ? boost::shared_ptr<ArticulatedBody>() : boost::shared_ptr<ArticulatedBody>(_abody); }

    /// Sets the articulated body corresponding to this body
    /**
     * \param body a pointer to the articulated body or NULL if this body is
     *        not a link in an articulated body
     */
    void set_articulated_body(boost::shared_ptr<ArticulatedBody> body) { _abody = body; }

    /// Gets the number of child links of this link
    unsigned num_child_links() const { return _outer_joints.size(); }

    /// Gets whether this body is an end-effector (i.e., the number of child links is zero) in an articulated body
    bool is_end_effector() const { assert (!_abody.expired()); return !is_base() && _outer_joints.empty(); }
 
    /// Removes all inner joints from this link
    void clear_inner_joints() { _inner_joints.clear(); }

    /// Removes all outer joints from this link
    void clear_outer_joints() { _outer_joints.clear(); }

    /// Gets the link index (returns std::numeric_limits<unsigned>::max() if not set)
    unsigned get_index() const { return _link_idx; }

    /// Sets the link index
    /**
     * This is set automatically by the articulated body.  Users should not
     * change this index or unknown behavior will result.
     */
    void set_index(unsigned index) { _link_idx = index; }

    /// Gets the set of inner joints for this link
    const std::set<JointPtr>& get_inner_joints() const { return _inner_joints; }
   
    /// Gets the list of outer joints for this link
    const std::set<JointPtr>& get_outer_joints() const { return _outer_joints; }

    /// Viscous coefficient for dampening the body motion
    Ravelin::VectorNd viscous_coeff;

    /// Validates the limit estimates
    virtual void validate_limit_estimates() { _vel_limit_exceeded = false; }

    /// Limit bound expansion scalar (default = 0.15 = 15%)
    double limit_bound_expansion;
    
    /// Compliance value, determines event type
    Compliance compliance;

  private:  
    template <class V>
    void get_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void get_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv);

    template <class V>
    void get_generalized_acceleration_generic(V& ga);

    template <class V>
    void set_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv);

    void set_force(const Ravelin::SForced& w);
    void invalidate_pose_vectors();
    void apply_generalized_impulse_single(const Ravelin::SharedVectorNd& gf);
    Ravelin::SharedMatrixNd& get_generalized_inertia_single(Ravelin::SharedMatrixNd& M);
    virtual Ravelin::SharedMatrixNd& get_generalized_inertia_inverse(Ravelin::SharedMatrixNd& M) const;
    Ravelin::SharedVectorNd& get_generalized_forces_single(Ravelin::SharedVectorNd& f);
    Ravelin::SharedVectorNd& convert_to_generalized_force_single(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::SharedVectorNd& gf);
    unsigned num_generalized_coordinates_single(DynamicBody::GeneralizedCoordinateType gctype) const;
    Ravelin::SharedMatrixNd& solve_generalized_inertia_single(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);
    Ravelin::SharedVectorNd& solve_generalized_inertia_single(const Ravelin::SharedVectorNd& b, Ravelin::SharedVectorNd& x);
    RigidBodyPtr get_parent_link(JointPtr j) const;
    RigidBodyPtr get_child_link(JointPtr j) const;
    void check_vel_limit_exceeded_and_update();

    /// Indicates whether link frame velocity is valid (up-to-date)
    bool _xdi_valid;

    /// Indicates whether inner joint frame velocity is valid (up-to-date)
    bool _xdj_valid;

    /// Indicates whether inertial frame velocity is valid (up-to-date)
    bool _xdm_valid;

    /// Indicates whether global frame velocity is valid (up-to-date)
    bool _xd0_valid;

    /// Indicates whether link frame acceleration is valid (up-to-date)
    bool _xddi_valid;

    /// Indicates whether global frame acceleration is valid (up-to-date)
    bool _xdd0_valid;

    /// Indicates whether inner joint frame acceleration is valid (up-to-date)
    bool _xddj_valid;

    /// Indicates whether inertial frame acceleration is valid (up-to-date)
    bool _xddm_valid;

    /// Indicates whether link frame force is valid (up-to-date)
    bool _forcei_valid;

    /// Indicates whether inner joint frame force is valid (up-to-date)
    bool _forcej_valid;

    /// Indicates whether inertial frame force is valid (up-to-date)
    bool _forcem_valid;

    /// Indicates whether global frame force is valid (up-to-date)
    bool _force0_valid;

    /// Indicates whether the global frame inertia matrix is valid
    bool _J0_valid;

    /// Indicates whether the link frame inertia matrix is valid
    bool _Ji_valid;

    /// Indicates whether the link com frame inertia matrix is valid
    bool _Jcom_valid;

    /// Indicates whether the inner joint frame inertia matix is valid 
    bool _Jj_valid;

    /// Spatial rigid body inertia matrix (global frame) 
    Ravelin::SpatialRBInertiad _J0;

    /// Velocity (global frame)
    Ravelin::SVelocityd _xd0;

    /// Acceleration (global frame)
    Ravelin::SAcceld _xdd0;

    /// Cumulative force on the body (global frame)
    Ravelin::SForced _force0;

    /// Spatial rigid body inertia matrix (inertial frame) 
    Ravelin::SpatialRBInertiad _Jm;

    /// Velocity (inertial frame)
    Ravelin::SVelocityd _xdm;

    /// Acceleration (inertial frame)
    Ravelin::SAcceld _xddm;

    /// Cumulative force on the body (inertial frame)
    Ravelin::SForced _forcem;

    /// Spatial rigid body inertia matrix (link frame) 
    Ravelin::SpatialRBInertiad _Ji;

    /// Velocity (link frame)
    Ravelin::SVelocityd _xdi;

    /// Acceleration (link frame)
    Ravelin::SAcceld _xddi;

    /// Cumulative force on the body (link frame)
    Ravelin::SForced _forcei;

    /// Spatial rigid body inertia matrix (inner joint frame) 
    Ravelin::SpatialRBInertiad _Jj;

    /// Velocity (inner joint frame)
    Ravelin::SVelocityd _xdj;

    /// Acceleration (inner joint frame)
    Ravelin::SAcceld _xddj;

    /// Cumulative force on the body (inner joint frame)
    Ravelin::SForced _forcej;

    /// Spatial rigid body inertia matrix (link COM frame) 
    Ravelin::SpatialRBInertiad _Jcom;

    /// Velocity (link com frame)
    Ravelin::SVelocityd _xdcom;

    /// Acceleration (link com frame)
    Ravelin::SAcceld _xddcom;

    /// Cumulative force on the body (com frame)
    Ravelin::SForced _forcecom;

    /// reference pose for this body
    boost::shared_ptr<Ravelin::Pose3d> _F;

    /// secondary pose for this body
    boost::shared_ptr<Ravelin::Pose3d> _F2;

    /// inertial pose for this body
    boost::shared_ptr<Ravelin::Pose3d> _jF;

    /// The link index (if a link in an articulated body)
    unsigned _link_idx;

    /// Flag for determining whether or not the body is physically enabled
    bool _enabled;

    /// Pointer to articulated body (if this body is a link)
    boost::weak_ptr<ArticulatedBody> _abody;

    /// Inner joints and associated data 
    std::set<JointPtr> _inner_joints;

    /// Outer joints and associated data 
    std::set<JointPtr> _outer_joints; 

    /// Lower velocity limits on the body 
    Ravelin::SVelocityd _vel_limit_lo;

    /// Upper velocity limits on the body
    Ravelin::SVelocityd _vel_limit_hi;

    /// Indicates whether the velocity limit has been exceeded
    bool _vel_limit_exceeded;

}; // end class

std::ostream& operator<<(std::ostream&, RigidBody&);

// incline inline functions
#include "RigidBody.inl"

} // end namespace

#endif

