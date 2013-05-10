/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPHERICAL_JOINT_H
#define _SPHERICAL_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated spherical joint
class SphericalJoint : public Joint
{
  public:
    enum Axis { eAxis1, eAxis2, eAxis3 };

    SphericalJoint();
    virtual ~SphericalJoint() {}
    SphericalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    Ravelin::Vector3d get_axis(Axis a) const;
    void set_axis(const Ravelin::Vector3d& axis, Axis a);
    virtual void update_spatial_axes();    
    virtual void determine_q(Ravelin::VectorNd& q);
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes_dot();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_dof() const { return 3; }
    virtual unsigned num_constraint_eqns() const { return 3; }
    virtual void evaluate_constraints(double C[]);

    /// Spherical joint is singular if sin(q1) = 0 and cos(q2) = 0
    virtual bool is_sngular_config() const { return std::fabs(std::sin(q[DOF_1])) < SINGULAR_TOL && std::fabs(std::cos(q[DOF_2])) < SINGULAR_TOL; }

    /// The tolerance to which a joint configuration is considered singular
    /**
     * \note if this tolerance is too low, then dynamics computation may
     * become unstable; if the tolerance is too high, then dynamics computation
     * will be slower.  A safe value is 1e-2.
     */
    double SINGULAR_TOL;

  private:
    bool assign_axes();
    virtual void calc_constraint_jacobian_euler(RigidBodyPtr body, unsigned index, double Cq[7]);
    virtual void calc_constraint_jacobian_dot_euler(RigidBodyPtr body, unsigned index, double Cq[7]);
    static bool rel_equal(double x, double y);
    Ravelin::Matrix3d get_rotation() const;

    /// The local joint axes
    Ravelin::Vector3d _u[3];

    /// The derivative of the spatial axis
    std::vector<Ravelin::Twistd> _s_dot;
}; // end class
} // end namespace

#endif

