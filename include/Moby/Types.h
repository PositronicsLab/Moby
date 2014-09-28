/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
#include <Ravelin/Vector2d.h>
#include <Ravelin/Vector3d.h>

namespace Moby {

class Polyhedron;
class Simulator;
class RigidBody;
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
class XMLTree;
class OSGGroupWrapper;
class AABB;
class OBB;
class BV;

/// Typedef to distinguish between a 2D vector and a point
typedef Ravelin::Vector2d Point2d;

/// Typedef to distinguish between a 3D vector and a point
typedef Ravelin::Vector3d Point3d;

/// Typedef to make specifying line segments easier
typedef std::pair<Point3d, Point3d> LineSeg3;

/// Typedef to make specifying line segments easier
typedef std::pair<Point2d, Point2d> LineSeg2;

/// reference frame type for reduced-coordinate dynamics computations
enum ReferenceFrameType { eGlobal, eLink, eLinkInertia, eLinkCOM, eJoint };

/// Polyhedron smart pointer
typedef boost::shared_ptr<Polyhedron> PolyhedronPtr;

/// Simulator smart pointer
typedef boost::shared_ptr<Simulator> SimulatorPtr;

/// Single body smart pointer
typedef boost::shared_ptr<SingleBody> SingleBodyPtr;

/// Rigid body smart pointer
typedef boost::shared_ptr<RigidBody> RigidBodyPtr;

/// Articulated body smart pointer
typedef boost::shared_ptr<ArticulatedBody> ArticulatedBodyPtr;

/// Reduced-coordinate articulated body smart pointer
typedef boost::shared_ptr<RCArticulatedBody> RCArticulatedBodyPtr;

/// Maximal-coordinate articulated body smart pointer
typedef boost::shared_ptr<MCArticulatedBody> MCArticulatedBodyPtr;

/// Reduced-coordinate articulated body joint smart pointer
typedef boost::shared_ptr<Joint> JointPtr;

/// Collision geometry smart pointer
typedef boost::shared_ptr<CollisionGeometry> CollisionGeometryPtr;

/// Primitive smart pointer
typedef boost::shared_ptr<Primitive> PrimitivePtr;

/// Base smart pointer
typedef boost::shared_ptr<Base> BasePtr;

/// Recurrent force smart pointer
typedef boost::shared_ptr<RecurrentForce> RecurrentForcePtr;

/// Dynamic body smart pointer
typedef boost::shared_ptr<DynamicBody> DynamicBodyPtr;

/// Bounding volume (BV) smart pointer
typedef boost::shared_ptr<BV> BVPtr;

/// Axis-aligned bounding box (AABB) smart pointer
typedef boost::shared_ptr<AABB> AABBPtr; 

/// Oriented bounding box (OBB) smart pointer
typedef boost::shared_ptr<OBB> OBBPtr;

/// OSGGroupWrapper smart pointer
typedef boost::shared_ptr<OSGGroupWrapper> OSGGroupWrapperPtr;

/// XML tree smart pointer
typedef boost::shared_ptr<XMLTree> XMLTreePtr;

} // end namespace

#endif

