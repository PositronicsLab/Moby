/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
  friend class RCArticulatedBody;
  friend class MCArticulatedBody;
  friend class Joint;

  public:
    RigidBody();
    virtual ~RigidBody() {}
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);
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
    virtual Ravelin::VectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv);
    virtual Ravelin::VectorNd& get_generalized_acceleration(Ravelin::VectorNd& ga);
    virtual void add_generalized_force(const Ravelin::VectorNd& gf);
    virtual void apply_generalized_impulse(const Ravelin::VectorNd& gf);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv);
    virtual Ravelin::MatrixNd& get_generalized_inertia(Ravelin::MatrixNd& M);
    virtual Ravelin::VectorNd& get_generalized_forces(Ravelin::VectorNd& f);
    virtual Ravelin::VectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, const Point3d& p, Ravelin::VectorNd& gf);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual Ravelin::MatrixNd& transpose_solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    Ravelin::MatrixNd& transpose_solve_generalized_inertia_single(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    virtual Ravelin::MatrixNd& solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    virtual Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
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
    virtual std::vector<Ravelin::SVelocityd>& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, std::vector<Ravelin::SVelocityd>& J);
    virtual std::vector<Ravelin::SVelocityd>& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, std::vector<Ravelin::SVelocityd>& J);
    const Ravelin::SForced& sum_forces();
    void reset_accumulators();

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

  private:  
    Ravelin::VectorNd& get_generalized_coordinates_single(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc);
    Ravelin::VectorNd& get_generalized_velocity_single(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv);
    Ravelin::VectorNd& get_generalized_acceleration_single(Ravelin::VectorNd& ga);
    void invalidate_pose_vectors();
    void apply_generalized_impulse_single(const Ravelin::VectorNd& gf);
    void set_generalized_coordinates_single(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc);
    void set_generalized_velocity_single(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv);
    Ravelin::MatrixNd& get_generalized_inertia_single(Ravelin::MatrixNd& M);
    Ravelin::VectorNd& get_generalized_forces_single(Ravelin::VectorNd& f);
    Ravelin::VectorNd& convert_to_generalized_force_single(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::VectorNd& gf);
    unsigned num_generalized_coordinates_single(DynamicBody::GeneralizedCoordinateType gctype) const;
    Ravelin::MatrixNd& solve_generalized_inertia_single(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    Ravelin::VectorNd& solve_generalized_inertia_single(const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
    RigidBodyPtr get_parent_link(JointPtr j) const;
    RigidBodyPtr get_child_link(JointPtr j) const;

    /// Indicates whether link frame velocity is valid (up-to-date)
    bool _xdi_valid;

    /// Indicates whether inner joint frame velocity is valid (up-to-date)
    bool _xdj_valid;

    /// Indicates whether inertial frame velocity is valid (up-to-date)
    bool _xdm_valid;

    /// Indicates whether link frame acceleration is valid (up-to-date)
    bool _xddi_valid;

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

    /// Indicates whether the global frame inertia matrix is valid
    bool _J0_valid;

    /// Indicates whether the link frame inertia matrix is valid
    bool _Ji_valid;

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

    Ravelin::LinAlgd _LA;
}; // end class

std::ostream& operator<<(std::ostream&, RigidBody&);

// incline inline functions
#include "RigidBody.inl"

} // end namespace

#endif

