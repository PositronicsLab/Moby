/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _COLLISION_DETECTION_H
#define _COLLISION_DETECTION_H

#include <map>
#include <set>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/sorted_pair>
#include <Moby/Base.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/Constraint.h>
#include <Moby/RigidBody.h>

namespace Moby {

class ConstraintSimulator;
class ArticulatedBody;
class CollisionGeometry;
class Triangle;
class IndexedTriArray;

/// Defines an abstract collision detection mechanism
/**
 * Contact finding and collision detection are two separate, but related, 
 * procedures. Collision detection quickly determines whether two bodies are 
 * geometrically intersecting at their current configurations. Contact finding
 * determines whether two bodies will come into contact (or are already in 
 * contact) in a given time interval; if the bodies do/will contact, contact
 * finding computes contact points and normals. 
 */
class CollisionDetection : public virtual Base
{
  public:
    CollisionDetection() {}
    virtual ~CollisionDetection() {}
    virtual void set_simulator(boost::shared_ptr<ConstraintSimulator> sim) {}
    virtual void broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi) = 0;
    virtual void find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<Constraint>& contacts, double TOL = NEAR_ZERO) = 0;

    /// Calculates the signed distance between two geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cg1, CollisionGeometryPtr cg2, Point3d& p1, Point3d& p2)
    {
      return CollisionGeometry::calc_signed_dist(cg1, cg2, p1, p2);
    }

    /// Get the shared pointer for this
    boost::shared_ptr<CollisionDetection> get_this() { return boost::dynamic_pointer_cast<CollisionDetection>(shared_from_this()); }

  protected:
    virtual double calc_next_CA_Euler_step(const PairwiseDistInfo& pdi) = 0;
    static Constraint create_contact(CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Ravelin::Vector3d& normal, double violation = 0.0);

  friend class ConstraintStabilization;
}; // end class

} // end namespace Moby

#endif

