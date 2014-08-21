/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _FIXED_JOINT_H
#define _FIXED_JOINT_H

#include <Ravelin/FixedJointd.h>
#include <Moby/Joint.h>

namespace Moby {

/// Defines a joint for fixing two bodies together or fixing one body to the ground
class FixedJoint : public Joint, public Ravelin::FixedJointd
{
  public:
    FixedJoint();
    virtual ~FixedJoint() {}
    FixedJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual unsigned num_dof() const { return FixedJointd::num_dof(); }
    virtual bool is_singular_config() const { return FixedJointd::is_singular_config(); }
    virtual void evaluate_constraints(double C[]) { FixedJointd::evaluate_constraints(C); }
    virtual const std::vector<Ravelin::SVelocityd>& get_spatial_axes_dot() { return FixedJointd::get_spatial_axes_dot(); }
    virtual void determine_q(Ravelin::VectorNd& q) { FixedJointd::determine_q(q); }
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose() { return FixedJointd::get_induced_pose(); }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
}; // end class
} // end namespace

#endif

