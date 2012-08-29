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
    Vector3 get_axis_global(Axis a) const;
    void set_axis_global(const Vector3& axis, Axis a);
    void set_axis_local(const Vector3& axis, Axis a);
    virtual void update_spatial_axes();    
    virtual void determine_Q();
    virtual const Matrix4& get_transform();
    virtual const SMatrix6N& get_spatial_axes(ReferenceFrameType type);
    virtual const SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_dof() const { return 2; }
    virtual void evaluate_constraints(Real C[]);

    /// Universal joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }
  
    /// Gets the unit vector describing the local axis of rotation for this joint
    /**
     * The local axis for this joint does not take the orientation of the 
     * inboard link into account; thus, if the orientation of the inboard link 
     * changes, then the local axis remains constant.
     * \sa get_axis_global()
     * \sa set_axis_global()
     */
    const Vector3& get_axis_local(Axis a) const { return _u[a]; }

  private:
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    static bool rel_equal(Real x, Real y);
    Matrix3 get_rotation() const;

    /// Rotation applied to axes to make computation simpler
    Matrix3 _R;

    /// The joint axes (inner link frame)
    Vector3 _u[2];

    /// The transform induced by this joint
    Matrix4 _T;

    /// The derivative of the spatial axis (link frame)
    SMatrix6N _si_dot;

    /// The derivative of the spatial axis (global frame)
    SMatrix6N _s0_dot;

    /// The second joint axis in outer link frame (only used for maximal coordinate articulated bodies)
    Vector3 _h2;
}; // end class
} // end namespace

#endif

