/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/**
 * A bunch of typedefs to make things more readable.
 */

#ifndef _MOBY_TYPES_H
#define _MOBY_TYPES_H

#include <utility>
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Point3d.h>
#include <Ravelin/Point2d.h>

namespace Moby {

class Polyhedron;
class Simulator;
class RigidBody;
class DeformableBody;
class SingleBody;
class ArticulatedBody;
class RCArticulatedBody;
class MCArticulatedBody;
class Joint;
class CollisionGeometry;
class Contact;
class Primitive;
class Base;
class RecurrentForce;
class DynamicBody;
class Octree;
class XMLTree;
class OSGGroupWrapper;
class AABB;
class OBB;
class BV;

/// Typedef to make specifying line segments easier
typedef std::pair<Ravelin::Point3d, Ravelin::Point3d> LineSeg3;

/// Typedef to make specifying line segments easier
typedef std::pair<Ravelin::Point2d, Ravelin::Point2d> LineSeg2;

/// reference frame type for reduced-coordinate dynamics computations
enum ReferenceFrameType { eGlobal, eLink };

/// Polyhedron smart pointer
typedef boost::shared_ptr<Polyhedron> PolyhedronPtr;

/// Simulator smart pointer
typedef boost::shared_ptr<Simulator> SimulatorPtr;

/// Single body smart pointer
typedef boost::shared_ptr<SingleBody> SingleBodyPtr;

/// Deformable body smart pointer
typedef boost::shared_ptr<DeformableBody> DeformableBodyPtr;

/// Rigid body smart pointer
typedef boost::shared_ptr<RigidBody> RigidBodyPtr;

/// Rigid body constant smart pointer
typedef boost::shared_ptr<const RigidBody> RigidBodyConstPtr;

/// Articulated body smart pointer
typedef boost::shared_ptr<ArticulatedBody> ArticulatedBodyPtr;

/// Articulated body constant smart pointer
typedef boost::shared_ptr<const ArticulatedBody> ArticulatedBodyConstPtr;

/// Reduced-coordinate articulated body smart pointer
typedef boost::shared_ptr<RCArticulatedBody> RCArticulatedBodyPtr;

/// Reduced-coordinate articulated body const smart pointer
typedef boost::shared_ptr<const RCArticulatedBody> RCArticulatedBodyConstPtr;

/// Maximal-coordinate articulated body smart pointer
typedef boost::shared_ptr<MCArticulatedBody> MCArticulatedBodyPtr;

/// Maximal-coordinate articulated body const smart pointer
typedef boost::shared_ptr<const MCArticulatedBody> MCArticulatedBodyConstPtr;

/// Reduced-coordinate articulated body constant joint smart pointer
typedef boost::shared_ptr<const Joint> JointConstPtr;

/// Reduced-coordinate articulated body joint smart pointer
typedef boost::shared_ptr<Joint> JointPtr;

/// Collision geometry smart pointer
typedef boost::shared_ptr<CollisionGeometry> CollisionGeometryPtr;

/// Collision geometry constant smart pointer
typedef boost::shared_ptr<const CollisionGeometry> CollisionGeometryConstPtr;

/// Primitive smart pointer
typedef boost::shared_ptr<Primitive> PrimitivePtr;

/// constant Primitive smart pointer
typedef boost::shared_ptr<const Primitive> PrimitiveConstPtr;

/// Base smart pointer
typedef boost::shared_ptr<Base> BasePtr;

/// Constant Base smart pointer
typedef boost::shared_ptr<const Base> BaseConstPtr;

/// Recurrent force smart pointer
typedef boost::shared_ptr<RecurrentForce> RecurrentForcePtr;

/// Recurrent force constant smart pointer
typedef boost::shared_ptr<const RecurrentForce> RecurrentForceConstPtr;

/// Dynamic body smart pointer
typedef boost::shared_ptr<DynamicBody> DynamicBodyPtr;

/// Octree smart pointer
typedef boost::shared_ptr<Octree> OctreePtr;

/// Bounding volume (BV) smart pointer
typedef boost::shared_ptr<BV> BVPtr;

/// Axis-aligned bounding box (AABB) smart pointer
typedef boost::shared_ptr<AABB> AABBPtr; 

/// Oriented bounding box (OBB) smart pointer
typedef boost::shared_ptr<OBB> OBBPtr;

/// OSGGroupWrapper smart pointer
typedef boost::shared_ptr<OSGGroupWrapper> OSGGroupWrapperPtr;

/// OSGGroupWrapper constant smart pointer
typedef boost::shared_ptr<const OSGGroupWrapper> OSGGroupWrapperConstPtr;

/// XML tree smart pointer
typedef boost::shared_ptr<XMLTree> XMLTreePtr;

/// constant XML tree smart pointer
typedef boost::shared_ptr<const XMLTree> XMLTreeConstPtr;

} // end namespace

#endif

