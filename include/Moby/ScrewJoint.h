/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _SCREW_JOINT_H
#define _SCREW_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated screw (helical) joint
class ScrewJoint : public Joint
{
  public:
    ScrewJoint();
    virtual ~ScrewJoint() {}
    ScrewJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    const Ravelin::Vector3d& get_axis() const { return _u; }
    void set_axis(const Ravelin::Vector3d& axis);
    virtual void update_spatial_axes();    
    virtual void determine_q(VectorN& q);
    virtual boost::shared_ptr<const Ravelin::Pose3d>& get_induced_pose();
    virtual const std::vector<Ravelin::SVelocityd>& get_spatial_axes_dot();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_dof() const { return 1; }
    virtual void evaluate_constraints(double C[]);

    /// Gets the pitch of the joint
    double get_pitch() const { return _pitch; }

    /// Screw joint can never be in a singular configuration
    virtual bool is_sngular_config() const { return false; }

    /// Sets the pitch of the joint
    void set_pitch(double pitch) { _pitch = pitch; update_spatial_axes(); }

  private:
    virtual void calc_constraint_jacobian(RigidBodyPtr body, unsigned index, double Cq[7]);
    virtual void calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, double Cq[7]);

    /// The pitch of the joint
    double _pitch;

    /// The joint axis (defined in inner link coordinates)
    Ravelin::Vector3d _u;

    /// Two unit vectors that make a orthonormal basis with _u
    Ravelin::Vector3d _ui, _uj; 

    /// The derivative of the spatial axis -- should be zero vector 6x1
    std::vector<Ravelin::SVelocityd> _s_dot;

    /// The joint axis (defined in outer link coordinates) [used only for maximal coordinate formulations]
    Ravelin::Vector3d _v2;
}; // end class
} // end namespace

#endif

