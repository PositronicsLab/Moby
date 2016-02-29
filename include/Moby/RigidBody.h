/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RIGID_BODY_H
#define _RIGID_BODY_H

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/weak_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/SForced.h>
#include <Ravelin/SVelocityd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/ControlledBody.h>
#include <Ravelin/RigidBodyd.h>

namespace osg { class Node; }

namespace Moby {

class ArticulatedBody;
class CollisionGeometry;

/// Represents a single rigid body
/**
 *  Contains information needed to represent a rigid body, including position
 *  and velocity (both linear and angular), mass, center of
 *  mass, inertia matrix, collision data, and visualization data.  This class
 *  is used for both non-articulated and articulated rigid bodies, though not
 *  all member data may be used in the latter.
 * \todo implement rest matrix
 */
class RigidBody : public virtual Ravelin::RigidBodyd, public virtual ControlledBody 
{
  friend class ArticulatedBody;
  friend class RCArticulatedBody;
  friend class MCArticulatedBody;
  friend class TimeSteppingSimulator;
  friend class Joint;

  public:
    enum Compliance { eRigid, eCompliant};
    RigidBody();
    virtual ~RigidBody() {}

    virtual void set_visualization_data(osg::Node* vdata) { Visualizable::set_visualization_data(vdata); }

    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    bool is_child_link(boost::shared_ptr<const RigidBody> query) const;
    bool is_descendant_link(boost::shared_ptr<const RigidBody> query) const;
    RigidBodyPtr get_parent_link() const;
    JointPtr get_inner_joint_explicit() const;
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) { prepare_to_calc_ode(x, t, dt, data); }
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void set_articulated_body(boost::shared_ptr<ArticulatedBody> body);
    virtual void set_inertia(const Ravelin::SpatialRBInertiad& J);
    virtual void apply_generalized_impulse(const Ravelin::SharedVectorNd& gj);

    template <class OutputIterator>
    OutputIterator get_parent_links(OutputIterator begin) const;

    template <class OutputIterator>
    OutputIterator get_child_links(OutputIterator begin) const;

    /// Collision geometries, if any, for this rigid body
    std::list<CollisionGeometryPtr> geometries;

    /// Compliance value, determines event type
    Compliance compliance;

  private:
    RigidBodyPtr get_parent_link(boost::shared_ptr<Ravelin::Jointd> j) const;
    RigidBodyPtr get_child_link(boost::shared_ptr<Ravelin::Jointd> j) const;

#ifdef USE_OSG
    osg::Node * inertia_viz;
#endif
}; // end class

std::ostream& operator<<(std::ostream&, RigidBody&);

// incline inline functions
#include "RigidBody.inl"

} // end namespace

#endif

