/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _REVOLUTE_JOINT_H
#define _REVOLUTE_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated revolute joint
class RevoluteJoint : public Joint
{
  public:
    RevoluteJoint();
    virtual ~RevoluteJoint() {}
    RevoluteJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    Ravelin::Vector3d get_axis_global() const;
    void set_axis_global(const Ravelin::Vector3d& axis);
    void set_axis_local(const Ravelin::Vector3d& axis);
    virtual void update_spatial_axes();    
    virtual void determine_q(Ravelin::VectorNd& q);
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_pose();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes_dot();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_dof() const { return 1; }
    virtual void evaluate_constraints(double C[]);

    /// Gets the unit vector describing the local axis of rotation for this joint
    /**
     * The local axis for this joint does not take the orientation of the 
     * inboard link into account; thus, if the orientation of the inboard link 
     * changes, then the local axis remains constant.
     * \sa get_axis_global()
     * \sa set_axis_global()
     */
    const Ravelin::Vector3d& get_axis_local() const { return _u; }

    /// Revolute joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

  private:
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, double Cq[7]);

    /// The joint axis (defined in inner link coordinates)
    Ravelin::Vector3d _u;

    /// Two unit vectors that make a orthonormal basis with _u
    Ravelin::Vector3d _ui, _uj; 

    /// The transform induced by this joint
    boost::shared_ptr<Ravelin::Pose3d> _T;

    /// The derivative of the spatial axis -- should be zero vector 6x1
    std::vector<Ravelin::Twistd> _si_deriv;

    /// The joint axis (defined in outer link coordinates) [used only for maximal coordinate formulations]
    Ravelin::Vector3d _v2;
}; // end class
} // end namespace

#endif
