/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _UNIVERSAL_JOINT_H
#define _UNIVERSAL_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated universal joint
class UniversalJoint : public Joint
{
  public:
    enum Axis { eAxis1, eAxis2 };

    UniversalJoint();
    virtual ~UniversalJoint() {}
    UniversalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    Ravelin::Vector3d get_axis(Axis a) const;
    void set_axis(const Ravelin::Vector3d& axis, Axis a);
    virtual void update_spatial_axes();    
    virtual void determine_q(Ravelin::VectorNd& q);
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes_dot();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_dof() const { return 2; }
    virtual void evaluate_constraints(double C[]);

    /// Universal joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }
  
  private:
    virtual void calc_constraint_jacobian_euler(RigidBodyPtr body, unsigned index, double Cq[7]);
    virtual void calc_constraint_jacobian_dot_euler(RigidBodyPtr body, unsigned index, double Cq[7]);
    static bool rel_equal(double x, double y);
    Ravelin::Matrix3d get_rotation() const;

    /// The joint axes (inner link frame)
    Ravelin::Vector3d _u[2];

    /// The derivative of the spatial axis
    std::vector<Ravelin::Twistd> _s_dot;

    /// The second joint axis in outer link frame (only used for maximal coordinate articulated bodies)
    Ravelin::Vector3d _h2;
}; // end class
} // end namespace

#endif

