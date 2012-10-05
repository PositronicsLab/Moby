/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _FIXED_JOINT_H
#define _FIXED_JOINT_H

#include <Moby/Joint.h>

namespace Moby {

/// Defines a joint for fixing two bodies together or fixing one body to the ground
class FixedJoint : public Joint
{
  public:
    FixedJoint();
    virtual ~FixedJoint() {}
    FixedJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual void update_spatial_axes();    
    virtual void set_inboard_link(RigidBodyPtr link);
    virtual void set_outboard_link(RigidBodyPtr link);
    virtual void determine_q(VectorN& q) { }
    virtual const Matrix4& get_transform();
    virtual const SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_dof() const { return 0; }
    virtual void evaluate_constraints(Real C[]);

    /// Fixed joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

  private:
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr, unsigned index, Real Cq[7]);
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr, unsigned index, Real Cq[7]);
    void setup_joint();

    /// The relative transform from the inboard link to the outboard link
    Matrix4 _T;

    /// The orientation constant that we attempt to maintain
    Vector3 _rconst;

    /// The vector from the inner link to the outer link in inner link frame
    Vector3 _ui;

    /// The time derivative of the spatial axis -- should be zero vector 6x1
    SMatrix6N _si_deriv;
}; // end class
} // end namespace

#endif
