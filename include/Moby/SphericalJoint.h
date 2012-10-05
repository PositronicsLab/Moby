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
    Vector3 get_axis_global(Axis a) const;
    void set_axis_global(const Vector3& axis, Axis a);
    void set_axis_local(const Vector3& axis, Axis a);
    virtual void update_spatial_axes();    
    virtual void determine_q(VectorN& q);
    virtual const Matrix4& get_transform();
    virtual const SMatrix6N& get_spatial_axes(ReferenceFrameType type);
    virtual const SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_dof() const { return 3; }
    virtual unsigned num_constraint_eqns() const { return 3; }
    virtual void evaluate_constraints(Real C[]);

    /// Spherical joint is singular if sin(q1) = 0 and cos(q2) = 0
    virtual bool is_singular_config() const { return std::fabs(std::sin(q[DOF_1])) < SINGULAR_TOL && std::fabs(std::cos(q[DOF_2])) < SINGULAR_TOL; }

    /// Gets the unit vector describing the local axis of rotation for this joint
    /**
     * The local axis for this joint does not take the orientation of the 
     * inboard link into account; thus, if the orientation of the inboard link 
     * changes, then the local axis remains constant.
     * \sa get_axis_global()
     * \sa set_axis_global()
     */
    const Vector3& get_axis_local(Axis a) const { return _u[a]; }

    /// The tolerance to which a joint configuration is considered singular
    /**
     * \note if this tolerance is too low, then dynamics computation may
     * become unstable; if the tolerance is too high, then dynamics computation
     * will be slower.  A safe value is 1e-2.
     */
    Real SINGULAR_TOL;

  private:
    bool assign_axes();
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    static bool rel_equal(Real x, Real y);
    Matrix3 get_rotation() const;

    /// Rotation applied to axes to make computation simpler
    Matrix3 _R;

    /// The local joint axes
    Vector3 _u[3];

    /// The transform induced by this joint
    Matrix4 _T;

    /// The derivative of the spatial axis (link frame)
    SMatrix6N _si_dot;

    /// The derivative of the spatial axis (global frame)
    SMatrix6N _s0_dot;
}; // end class
} // end namespace

#endif

