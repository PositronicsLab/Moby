/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _UNIVERSAL_JOINT_H
#define _UNIVERSAL_JOINT_H

#include <Ravelin/UniversalJointd.h>
#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated universal joint
class UniversalJoint : public Joint, public Ravelin::UniversalJointd
{
  public:
    UniversalJoint();
    virtual ~UniversalJoint() {}
    UniversalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual unsigned num_dof() const { return UniversalJointd::num_dof(); }
    virtual bool is_singular_config() const { return UniversalJointd::is_singular_config(); }
    virtual void evaluate_constraints(double C[]) { UniversalJointd::evaluate_constraints(C); }
    virtual const std::vector<Ravelin::SVelocityd>& get_spatial_axes_dot() { return UniversalJointd::get_spatial_axes_dot(); }
    virtual void determine_q(Ravelin::VectorNd& q) { UniversalJointd::determine_q(q); }
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose() { return UniversalJointd::get_induced_pose(); }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
}; // end class
} // end namespace

#endif

