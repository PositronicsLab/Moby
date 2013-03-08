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
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Quat.h>
#include <Moby/Vector3.h>
#include <Moby/SVector6.h>
#include <Moby/SpatialRBInertia.h>
#include <Moby/SpatialTransform.h>

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
  public:

    /// Data for determining joint hierarchy
    struct InnerJointData
    {
      // the inner joint in question
      boost::weak_ptr<Joint> inner_joint;

      // vector from inner joint to center-of-mass of this link (outboard link frame)
      Vector3 joint_to_com_vec_of;

      // vector from inner joint to center-of-mass of this link (inner joint frame)
      Vector3 joint_to_com_vec_jf;


      // the parent link
      boost::weak_ptr<RigidBody> parent;
    };

    /// Data for determining joint hierarchy
    struct OuterJointData
    {
      // the outer joint in question
      boost::weak_ptr<Joint> outer_joint;

      // vector from c.o.m. of this link to inner joint of child link (this link coordinates)
      Vector3 com_to_joint_vec;

      // the child link
      boost::weak_ptr<RigidBody> child;
    };

    RigidBody();
    virtual ~RigidBody() {}
    virtual void integrate(Real t, Real h, boost::shared_ptr<Integrator<VectorN> > integrator);
    void set_avel(const Vector3& avel);
    void set_lvel(const Vector3& lvel);
    void add_force(const Vector3& force);
    void add_torque(const Vector3& torque);
    void set_transform(const Matrix4& transform);
    void set_transform(const Quat& q, const Vector3& x);
    void set_inertia(const Matrix3& m);
    void set_enabled(bool flag);
    void set_orientation(const Quat& q);
    void set_position(const Vector3& pos);
    void add_force(const Vector3& f, const Vector3& p);
    void apply_impulse(const Vector3& j, const Vector3& p);
    void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& p);  
    void set_mass(Real mass);
    virtual void transform(const Matrix4& transform) { set_transform(transform * _F); }
    virtual void calc_fwd_dyn(Real dt);
    SpatialRBInertia get_spatial_iso_inertia(ReferenceFrameType rftype) const;

    virtual void set_visualization_data(osg::Node* vdata) { Visualizable::set_visualization_data(vdata); synchronize(); }

    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    Real calc_point_accel(const Vector3& point, const Vector3& dir, Real dt);
    static Real calc_sep_accel(RigidBody& rb1, RigidBody& rb2, const Vector3& point, const Vector3& dir, const Vector3& dir_dot, Real dt);
    bool is_child_link(RigidBodyConstPtr query) const;
    bool is_descendant_link(RigidBodyConstPtr query) const;
    SpatialTransform get_spatial_transform_forward() const;
    SpatialTransform get_spatial_transform_backward() const;
    SpatialTransform get_spatial_transform_link_to_global() const;
    SpatialTransform get_spatial_transform_global_to_link() const;
    SVector6 get_spatial_accel(ReferenceFrameType rftype);
    void set_spatial_accel(const SVector6& a, ReferenceFrameType rftype);
    SVector6 get_spatial_velocity(ReferenceFrameType rftype);
    void set_spatial_velocity(const SVector6& v, ReferenceFrameType rftype);
    boost::shared_ptr<const DynamicBody> get_dynamic_body() const;
    DynamicBodyPtr get_dynamic_body();
    virtual VectorN& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gc);
    virtual VectorN& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv);
    virtual VectorN& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, VectorN& ga);
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gv);
    virtual MatrixN& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M);
    virtual VectorN& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, VectorN& f);
    virtual VectorN& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& f, const Vector3& t, VectorN& gf);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    Vector3 calc_inertial_forces() const;
    const InnerJointData& get_inner_joint_data(RigidBodyPtr parent) const { return get_inner_joint_data(parent); }
    InnerJointData& get_inner_joint_data(RigidBodyPtr parent);
    const InnerJointData& get_inner_joint_data(JointPtr inner_joint) const { return get_inner_joint_data(inner_joint); }
    InnerJointData& get_inner_joint_data(JointPtr inner_joint);
    const OuterJointData& get_outer_joint_data(RigidBodyPtr child) const { return get_outer_joint_data(child); }
    OuterJointData& get_outer_joint_data(RigidBodyPtr child);
    const OuterJointData& get_outer_joint_data(JointPtr outer_joint) const { return get_outer_joint_data(outer_joint); }
    OuterJointData& get_outer_joint_data(JointPtr outer_joint);
    RigidBodyPtr get_parent_link() const;
    JointPtr get_inner_joint_implicit() const;
    void add_inner_joint(RigidBodyPtr parent, JointPtr j, const Vector3& joint_to_com_vec_joint, const Vector3& joint_to_com_vec_link);
    void add_outer_joint(RigidBodyPtr child, JointPtr j, const Vector3& com_to_joint_vec_link);
    bool remove_inner_joints(RigidBodyPtr parent);
    bool remove_inner_joint(JointPtr joint);
    bool remove_outer_joints(RigidBodyPtr child);
    bool remove_outer_joint(JointPtr joint);
    virtual Real calc_kinetic_energy() const;
    virtual Vector3 calc_point_vel(const Vector3& p) const;
    virtual void update_velocity(const EventProblemData& q);
    virtual void update_event_data(EventProblemData& q);
    RigidBodyPtr get_child_link(unsigned i) const;
    bool is_base() const;
    bool is_ground() const;

    template <class OutputIterator>
    OutputIterator get_parent_links(OutputIterator begin) const;

    template <class OutputIterator>
    OutputIterator get_child_links(OutputIterator begin) const;

    template <class OutputIterator>    
    OutputIterator get_all_collision_geometries(OutputIterator begin) const;

    /// Collision geometries, if any, for this rigid body
    std::list<CollisionGeometryPtr> geometries;

    /// Gets the shared pointer for <b>this</b>
    RigidBodyPtr get_this() { return boost::dynamic_pointer_cast<RigidBody>(shared_from_this()); }
  
    /// Gets the shared const pointer for <b>this</b>
    RigidBodyConstPtr get_this() const { return boost::dynamic_pointer_cast<const RigidBody>(shared_from_this()); }

    /// Gets the current transform of this body
    const Matrix4& get_transform() const { return _F; }

    /// Synonym for get_mass() (implements SingleBody::calc_mass())
    Real calc_mass() const { return _mass; }

    /// Gets the mass of this body
    Real get_mass() const { return _mass; }
    
    /// Gets the inverse mass of this body
    Real get_inv_mass() const { return _inv_mass; }
    
    /// Gets the 3x3 inertia tensor of this body
    const Matrix3& get_inertia() const { return _J; }

    /// Gets the inverse of the 3x3 inertia tensor of this body
    const Matrix3& get_inv_inertia() const { return _invJ; }

    /// Gets the position of this body
    const Vector3& get_position() const { return _x; }
    
    /// Gets the quaternion orientation of this body
    const Quat& get_orientation() const { return _q; }
    
    /// Gets the linear velocity of this body
    const Vector3& get_lvel() const { return _xd; }
    
    /// Gets the angular velocity of this body (world frame)
    const Vector3& get_avel() const { return _omega; }
    
    /// Gets the linear acceleration of this body
    /**
     * \note It is the user's responsibility to call calc_fwd_dyn() before
     * calling this method (if necessary)!
     */
    const Vector3& get_laccel() { return _xdd; }

    /// Gets the angular acceleration of this body
    /**
     * \note It is the user's responsibility to call calc_fwd_dyn() before
     * calling this method (if necessary)!
     */
    const Vector3& get_aaccel() { return _alpha; }

    /// Sets the angular acceleration of this link
    void set_aaccel(const Vector3& a) { _alpha = a; }

    /// Sets the angular acceleration of this link
    void set_laccel(const Vector3& a) { _xdd = a; }

    /// Resets the force and torque accumulators of this body
    void reset_accumulators() { _forces = ZEROS_3; _torques = ZEROS_3; }
    
    /// Gets the sum of forces on this body 
    /**
      * The frame of the sum of forces is centered at the center-of-mass of
      * the body and is aligned with the global frame.
      */
    const Vector3& sum_forces() const { return _forces; }
    
    /// Gets the sum of torques on this body
    /**
      * The frame of the sum of torques is centered at the center-of-mass of
      * the body and is aligned with the global frame.
      */
    const Vector3& sum_torques() const { return _torques; }

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

    /// Gets the list of inner joint data for this link
    const std::list<InnerJointData>& get_inner_joints_data() const { return _inner_joints; }
   
    /// Gets the list of outer joint data for this link
    const std::list<OuterJointData>& get_outer_joints_data() const { return _outer_joints; }

    /// Coulomb coefficient for dampening the body motion
    VectorN coulomb_coeff;

    /// Viscous coefficient for dampening the body motion
    VectorN viscous_coeff;

  protected:
    /// Gets the transform for visualization
    virtual const Matrix4* get_visualization_transform() { return &_F; }

  private:  
    void invalidate_position();
    void invalidate_velocity();
    void synchronize();
    static bool valid_transform(const MatrixN& T, Real tol);

    /// Mass of the rigid body
    Real _mass;

    /// Inverse of the mass of the rigid body
    Real _inv_mass;

    /// Inertia matrix for the rigid body
    Matrix3 _J;

    /// Inverse inertia matrix for the rigid body
    Matrix3 _invJ;

    /// Rigid body angular velocity (world frame)
    Vector3 _omega;

    /// Position of center of mass
    Vector3 _x;
   
    /// Velocity of center of mass
    Vector3 _xd;

    /// Body orientation
    Quat _q;
   
    /// Transform for this body
    Matrix4 _F;

    /// Cumulative force on the body
    Vector3 _forces;

    /// Cumulative torque on the body (world frame)
     Vector3 _torques;

    /// The link index (if a link in an articulated body)
    unsigned _link_idx;

    /// Linear acceleration
    Vector3 _xdd;

    /// Angular acceleration
    Vector3 _alpha;  

    /// Flag for determining whether or not the body is physically enabled
    bool _enabled;

    /// Pointer to articulated body (if this body is a link)
    boost::weak_ptr<ArticulatedBody> _abody;

    /// Inner joints and associated data 
    std::list<InnerJointData> _inner_joints;

    /// Outer joints and associated data 
    std::list<OuterJointData> _outer_joints; 

    static VectorN ode_p(const VectorN& x, Real t, void* data);
    static VectorN ode_v(const VectorN& x, Real t, void* data);
}; // end class

std::ostream& operator<<(std::ostream&, const RigidBody&);

// incline inline functions
#include "RigidBody.inl"

} // end namespace

#endif

