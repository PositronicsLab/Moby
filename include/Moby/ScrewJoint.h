/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
    Vector3 get_axis_global() const;
    void set_axis_global(const Vector3& axis);
    void set_axis_local(const Vector3& axis);
    virtual void update_spatial_axes();    
    virtual void determine_q(VectorN& q);
    virtual const Matrix4& get_transform();
    virtual const SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_dof() const { return 1; }
    virtual void evaluate_constraints(Real C[]);

    /// Gets the pitch of the joint
    Real get_pitch() const { return _pitch; }

    /// Gets the unit vector describing the local axis of rotation/translation for this joint
    /**
     * The local axis for this joint does not take the orientation of the 
     * inboard link into account; thus, if the orientation of the inboard link 
     * changes, then the local axis remains constant.
     * \sa get_axis_global()
     * \sa set_axis_global()
     */
    const Vector3& get_axis_local() const { return _u; }

    /// Screw joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

    /// Sets the pitch of the joint
    void set_pitch(Real pitch) { _pitch = pitch; update_spatial_axes(); }

  private:
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);

    /// The pitch of the joint
    Real _pitch;

    /// The joint axis (defined in inner link coordinates)
    Vector3 _u;

    /// Two unit vectors that make a orthonormal basis with _u
    Vector3 _ui, _uj; 

    /// The transform induced by this joint
    Matrix4 _T;

    /// The derivative of the spatial axis -- should be zero vector 6x1
    SMatrix6N _si_deriv;

    /// The joint axis (defined in outer link coordinates) [used only for maximal coordinate formulations]
    Vector3 _v2;
}; // end class
} // end namespace

#endif

