/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _REVOLUTE_JOINT_H
#define _REVOLUTE_JOINT_H

#include <Ravelin/RevoluteJointd.h>
#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated revolute joint
class RevoluteJoint : public Joint, public Ravelin::RevoluteJointd
{
  public:
    RevoluteJoint();
    virtual ~RevoluteJoint() {}
    RevoluteJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual unsigned num_dof() const { return RevoluteJointd::num_dof(); }
    virtual bool is_singular_config() const { return RevoluteJointd::is_singular_config(); }
    virtual void evaluate_constraints(double C[]) { RevoluteJointd::evaluate_constraints(C); }
    virtual const std::vector<Ravelin::SVelocityd>& get_spatial_axes_dot() { return RevoluteJointd::get_spatial_axes_dot(); }
    virtual void determine_q(Ravelin::VectorNd& q) { RevoluteJointd::determine_q(q); }
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose() { return RevoluteJointd::get_induced_pose(); }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
}; // end class
} // end namespace

#endif
