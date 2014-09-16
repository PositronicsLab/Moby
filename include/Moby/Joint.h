/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _JOINT_H
#define _JOINT_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/Jointd.h>
#include <Moby/Base.h>
#include <Moby/DynamicBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Visualizable.h>

namespace Moby {

class VisualizationData;

/// Defines a joint used in articulated rigid body dynamics simulation
/**
 * \todo implement a rest position for q?
 */
class Joint : public Visualizable, public virtual Ravelin::Jointd
{
  public:
    /**  
     * The constraint type is implicit if constraint forces must be determined
     * to hold the inner and outer links together and explicit otherwise.
     */ 
    enum ConstraintType { eUnknown, eExplicit, eImplicit };
    enum DOFs { DOF_1=0, DOF_2=1, DOF_3=2, DOF_4=3, DOF_5=4, DOF_6=5 };

    Joint();
    Joint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual ~Joint() {}
    void add_force(const Ravelin::VectorNd& force);
    void reset_force();    
    void set_location(const Point3d& p, RigidBodyPtr inboard, RigidBodyPtr outboard);
    Ravelin::VectorNd& get_scaled_force(Ravelin::VectorNd& f);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_link(RigidBodyPtr link, bool update_pose);
    virtual void set_outboard_link(RigidBodyPtr link, bool update_pose);
    ConstraintType get_constraint_type() const { return _constraint_type; }
    void evaluate_constraints_dot(double C[6]);

    /// Sets whether this constraint is implicit or explicit (or unknown)
    void set_constraint_type(ConstraintType type) { _constraint_type = type; }

    /// Gets the shared pointer to this joint
    JointPtr get_this() { return boost::dynamic_pointer_cast<Joint>(shared_from_this()); }

    /// Gets the constant shared pointer to this joint
    boost::shared_ptr<const Joint> get_this() const { return boost::dynamic_pointer_cast<const Joint>(shared_from_this()); }

    /// Gets the inboard link for this joint
    RigidBodyPtr get_inboard_link() const { return (_inboard_link.expired()) ? RigidBodyPtr() : RigidBodyPtr(_inboard_link); }

    /// Gets the outboard link for this joint
    RigidBodyPtr get_outboard_link() const { return (_outboard_link.expired()) ? RigidBodyPtr() : RigidBodyPtr(_outboard_link); }

    /// Gets the articulated body corresponding to this body
    /**
     * \return a pointer to the articulated body, or NULL if this body is not 
     *         a link an articulated body
     */
    ArticulatedBodyPtr get_articulated_body() { return (_abody.expired()) ? ArticulatedBodyPtr() : ArticulatedBodyPtr(_abody); }

    /// Sets the articulated body corresponding to this body
    /**
     * \param body a pointer to the articulated body or NULL if this body is
     *        not a link in an articulated body
     */
    void set_articulated_body(ArticulatedBodyPtr abody) { _abody = abody; }

    /// The lower joint limit 
    Ravelin::VectorNd lolimit;

    /// The upper joint limit 
    Ravelin::VectorNd hilimit;

    /// The maximum force that can be applied to each DOF of the joint
    Ravelin::VectorNd maxforce;

    /// The acceleration of this joint
    Ravelin::VectorNd qdd;

    /// The coefficient of restitution applied when this joint reaches a limit
    double limit_restitution;

    /// The coulomb friction coefficient for this joint
    double mu_fc;

    /// The viscous friction coefficient for this joint
    double mu_fv;

    /// The total friction applied to this joint (computed by forward dynamics / impulse application)
    Ravelin::VectorNd ff;

    /// The actuator force (user/controller sets this)
    Ravelin::VectorNd force;

    /// Constraint forces calculated by forward dynamics
    Ravelin::VectorNd lambda;

    /// Gets the joint index (returns UINT_MAX if not set)
    unsigned get_index() const { return _joint_idx; }

    /// Sets the joint index
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_index(unsigned index) { _joint_idx = index; }

    /// Gets the constraint index for this joint
    unsigned get_constraint_index() const { return _constraint_idx; }

    /// Sets the constraint index for this joint
    void set_constraint_index(unsigned idx) { _constraint_idx = idx; } 

    /// Sets the coordinate index for this joint
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_coord_index(unsigned index) { _coord_idx = index; }

    /// Gets the starting coordinate index for this joint
    unsigned get_coord_index() const { return _coord_idx; }

  protected:
    void determine_q_dot();
    void invalidate_pose_vectors() { get_outboard_link()->invalidate_pose_vectors(); }

    /// Computes the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
//    virtual void calc_constraint_jacobian(RigidBodyPtr body, unsigned index, double Cq[]) = 0;
 
     /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
//    virtual void calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, double Cq[]) = 0;

    /// Method for initializing all variables in the joint
    /**
     * This method should be called at the beginning of all constructors of
     * all derived classes.
     */
    virtual void init_data();

    /// Set whether _q_tare needs to be determined
    bool _determine_q_tare;

  private:
    ConstraintType _constraint_type;
    unsigned _joint_idx;
    unsigned _coord_idx;
    unsigned _constraint_idx;
    boost::shared_ptr<Ravelin::Pose3d> _vtransform;
    boost::weak_ptr<RigidBody> _inboard_link;
    boost::weak_ptr<RigidBody> _outboard_link;
    boost::weak_ptr<ArticulatedBody> _abody;

    // shared linear algebra object
    static Ravelin::LinAlgd _LA;

}; // end class
} // end namespace

#endif

