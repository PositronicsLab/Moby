/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_JOINT_H
#define _MOBY_JOINT_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/Jointd.h>
#include <Moby/Base.h>
#include <Moby/ControlledBody.h>
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
  friend class ArticulatedBody;

  public:
    /**  
     * The constraint type is implicit if constraint forces must be determined
     * to hold the inner and outer links together and explicit otherwise.
     */ 
    enum ConstraintType { eUnknown, eExplicit, eImplicit };
    enum DOFs { DOF_1=0, DOF_2=1, DOF_3=2, DOF_4=3, DOF_5=4, DOF_6=5 };

    // the constraint stiffness for this joint, if it is implicit (0 = no stabilization, 1 = immediate stabilization); default = 0.2
    double implicit_constraint_stiffness;

    /// if set, use qd_des for inverse dynamics
    bool qd_des_set;  

    /// max torque soon initialized to infinity
    Ravelin::VectorNd max_torque;

    /// the set point of the spring
    Ravelin::VectorNd set_point;

    /// does this joint use spring and damper
    bool spring_joint;

    /// vectors for the stiffness and damping of joints
    Ravelin::VectorNd stiffness, damping;

    /// desired joint velocity (only used if qd_des_set is true)
    Ravelin::VectorNd qd_des;

    Joint();
    Joint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual ~Joint() {}
    void set_location(const Point3d& p, RigidBodyPtr inboard, RigidBodyPtr outboard);
    Ravelin::VectorNd& get_scaled_force(Ravelin::VectorNd& f);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_link(RigidBodyPtr link, bool update_pose);
    virtual void set_outboard_link(RigidBodyPtr link, bool update_pose);
    void apply_impulse(const Ravelin::VectorNd& lambda);

    /// Gets the inboard link for this joint
    RigidBodyPtr get_inboard_link() const { return (_inboard_link.expired()) ? RigidBodyPtr() : boost::dynamic_pointer_cast<RigidBody>(boost::shared_ptr<Ravelin::RigidBodyd>(_inboard_link)); }

    /// Gets the outboard link for this joint
    RigidBodyPtr get_outboard_link() const { return (_outboard_link.expired()) ? RigidBodyPtr() : boost::dynamic_pointer_cast<RigidBody>(boost::shared_ptr<Ravelin::RigidBodyd>(_outboard_link)); }

    /// The lower joint limit 
    Ravelin::VectorNd lolimit;

    /// The upper joint limit 
    Ravelin::VectorNd hilimit;

    /// The coefficient of restitution applied when this joint reaches a limit
    double limit_restitution;

    /// The coulomb friction coefficient for this joint
    double mu_fc;

    /// The viscous friction coefficient for this joint
    double mu_fv;

    /// The depth of the compliant layer around the joint limit
    double compliant_layer_depth;

  protected:
    void determine_q_dot();

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

  private:
    boost::shared_ptr<Ravelin::Pose3d> _vtransform;
}; // end class
} // end namespace

#endif

