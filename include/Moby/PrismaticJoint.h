/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _PRISMATIC_JOINT_H
#define _PRISMATIC_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines an actuated prismatic joint
/**
 * \todo implement a rest position for q?
 */
class PrismaticJoint : public Joint
{
  public:
    PrismaticJoint();
    virtual ~PrismaticJoint() {}
    PrismaticJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    Vector3 get_axis_global() const;
    void set_axis_global(const Vector3& axis);
    void set_axis_local(const Vector3& axis);
    virtual void update_spatial_axes();    
    virtual void determine_Q();
    virtual const Matrix4& get_transform();
    virtual SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_dof() const { return 1; } 
    virtual void evaluate_constraints(Real C[]);

    /// Gets the unit vector describing the local axis of translation for this joint
    /**
     * The local axis for this joint does not take the orientation of the 
     * inboard link into account; thus, if the orientation of the inboard link 
     * changes, then the local axis remains constant.
     * \sa get_axis_global()
     * \sa set_axis_global()
     */
    const Vector3& get_axis_local() const { return _u; }

    /// Prismatic joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }
    
  private:
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[7]);

    /// The axis of the joint (inner link coordinates)
    Vector3 _u;

    /// The derivative of the spatial axis matrix (used in reduced-coordinate articulated bodies only)
    SMatrix6N _si_deriv;

    /// The 4x4 homogeneous transform induced by this joint
    Matrix4 _T;

    /// Vector attached to inner link and initially orthogonal to joint axis (used in maximal coordinate articulated bodies only); vector specified in inner link frame
    Vector3 _ui;

    /// Vector attached to outer link and initially orthogonal to joint axis (used in maximal coordinate articulated bodies only); vector specified in outer link frame
    Vector3 _uj;

    /// The joint axis on the second body (used in maximal coordinate articulated bodies only); vector specified in outer link frame 
    Vector3 _v2;
}; // end class
} // end namespace

#endif

