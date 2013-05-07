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
#include <Ravelin/Wrenchd.h>
#include <Ravelin/Twistd.h>
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
  public:

    /// Data for determining joint hierarchy
    struct InnerJointData
    {
      // the inner joint in question
      boost::weak_ptr<Joint> inner_joint;

      // vector from inner joint to center-of-mass of this link (outboard link frame)
      Ravelin::Vector3d joint_to_com_vec_of;

      // vector from inner joint to center-of-mass of this link (inner joint frame)
      Ravelin::Vector3d joint_to_com_vec_jf;


      // the parent link
      boost::weak_ptr<RigidBody> parent;
    };

    /// Data for determining joint hierarchy
    struct OuterJointData
    {
      // the outer joint in question
      boost::weak_ptr<Joint> outer_joint;

      // vector from c.o.m. of this link to inner joint of child link (this link coordinates)
      Ravelin::Vector3d com_to_joint_vec;

      // the child link
      boost::weak_ptr<RigidBody> child;
    };

    RigidBody();
    virtual ~RigidBody() {}
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);
    void add_wrench(const Ravelin::Wrenchd& w);
    void set_pose(const Ravelin::Pose3d& pose);
    void set_inertia(const Ravelin::SpatialRBInertiad& m);
    void set_enabled(bool flag);
    void apply_impulse(const Ravelin::Wrenchd& w);
    virtual void transform(const Ravelin::Pose3d& transform) { set_pose(transform * (*_F)); }
    virtual void calc_fwd_dyn(double dt);
    const Ravelin::SpatialRBInertiad& get_inertia() const;

    virtual void set_visualization_data(osg::Node* vdata) { Visualizable::set_visualization_data(vdata); synchronize(); }

    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    bool is_child_link(boost::shared_ptr<const RigidBody> query) const;
    bool is_descendant_link(boost::shared_ptr<const RigidBody> query) const;
    Ravelin::Twistd& accel() { return _xdd; }
    const Ravelin::Twistd& accel() const { return _xdd; } 
    Ravelin::Twistd& velocity();
    const Ravelin::Twistd& velocity() const { return _xd; }
    boost::shared_ptr<const DynamicBody> get_dynamic_body() const;
    DynamicBodyPtr get_dynamic_body();
    virtual Ravelin::VectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv);
    virtual Ravelin::VectorNd& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& ga);
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gf);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv);
    virtual Ravelin::MatrixNd& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::MatrixNd& M);
    virtual Ravelin::VectorNd& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& f);
    virtual Ravelin::VectorNd& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Ravelin::Wrenchd& w, Ravelin::VectorNd& gf);
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual Ravelin::MatrixNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gc, const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    virtual Ravelin::VectorNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gc, const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
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
    void add_inner_joint(RigidBodyPtr parent, JointPtr j, const Ravelin::Vector3d& joint_to_com_vec_joint, const Ravelin::Vector3d& joint_to_com_vec_link);
    void add_outer_joint(RigidBodyPtr child, JointPtr j, const Ravelin::Vector3d& com_to_joint_vec_link);
    bool remove_inner_joints(RigidBodyPtr parent);
    bool remove_inner_joint(JointPtr joint);
    bool remove_outer_joints(RigidBodyPtr child);
    bool remove_outer_joint(JointPtr joint);
    virtual double calc_kinetic_energy() const;
    virtual Ravelin::Vector3d calc_point_vel(const Ravelin::Point3d& p) const;
    virtual void update_velocity(const EventProblemData& q);
    virtual void update_event_data(EventProblemData& q);
    RigidBodyPtr get_child_link(unsigned i) const;
    bool is_base() const;
    bool is_ground() const;
    virtual Ravelin::Point3d get_position() const;
    boost::shared_ptr<const Ravelin::Pose3d> get_computation_frame() const;
    virtual void set_computation_frame_type(ReferenceFrameType rftype);

    virtual void calc_event_data(const Event& e, const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q);
    virtual void calc_event_data(const Event& e1, const Event& e2, const Ravelin::MatrixNd& M);
    Ravelin::Wrenchd calc_inertial_forces() const;

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
    boost::shared_ptr<const RigidBody> get_this() const { return boost::dynamic_pointer_cast<const RigidBody>(shared_from_this()); }

    /// Gets the current transform of this body
    boost::shared_ptr<Ravelin::Pose3d> get_pose() const { return _F; }

    /// Synonym for get_mass() (implements SingleBody::calc_mass())
    double calc_mass() const { return _J.m; }

    /// Gets the mass of this body
    double get_mass() const { return _J.m; }
    
    /// Resets the force and torque accumulators of this body
    void reset_accumulators() { _wrench.set_zero(); }
    
    /// Gets the external wrench on this body 
    const Ravelin::Wrenchd& sum_wrenches() const { return _wrench; }
    
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

    /// Viscous coefficient for dampening the body motion
    Ravelin::VectorNd viscous_coeff;

  protected:
    /// Gets the transform for visualization
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_visualization_pose() { return _F; }

  private:  
    void invalidate_position();
    void invalidate_velocity();
    void synchronize();

    /// Spatial rigid body inertia matrix (given computation frame) 
    Ravelin::SpatialRBInertiad _J;

    /// Velocity (given computation frame)
    Ravelin::Twistd _xd;

    /// reference pose for this body
    boost::shared_ptr<Ravelin::Pose3d> _F;

    /// inertial pose for this body
    boost::shared_ptr<Ravelin::Pose3d> _jF;

    /// Cumulative wrench on the body
    Ravelin::Wrenchd _wrench;

    /// The link index (if a link in an articulated body)
    unsigned _link_idx;

    /// Acceleration (given computation frame)
    Ravelin::Twistd _xdd;

    /// Flag for determining whether or not the body is physically enabled
    bool _enabled;

    /// Pointer to articulated body (if this body is a link)
    boost::weak_ptr<ArticulatedBody> _abody;

    /// Inner joints and associated data 
    std::list<InnerJointData> _inner_joints;

    /// Outer joints and associated data 
    std::list<OuterJointData> _outer_joints; 

    static Ravelin::VectorNd ode_p(const Ravelin::VectorNd& x, double t, void* data);
    static Ravelin::VectorNd ode_v(const Ravelin::VectorNd& x, double t, void* data);
    Ravelin::LinAlgd _LA;
}; // end class

std::ostream& operator<<(std::ostream&, const RigidBody&);

// incline inline functions
#include "RigidBody.inl"

} // end namespace

#endif

