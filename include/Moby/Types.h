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
#ifdef BUILD_ARBITRARY_PRECISION
#include <Moby/mpreal.h>
#endif

namespace Moby {

class Vector2;
class Vector3;
class VectorN;
class SVector6;
class Matrix3;
class Matrix4;
class MatrixN;
class MatrixNN;
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
class SMatrix6N;
class Octree;
class XMLTree;
class OSGGroupWrapper;
class AABB;
class OBB;
class BV;

#ifdef BUILD_SINGLE
/// default floating-point type
typedef float Real;
	
/// floating-point type to use when extra precision is necessary
typedef double LongReal;
#else
#ifdef BUILD_DOUBLE
/// default floating-point type
typedef double Real;
	
/// floating-point type to use when extra precision is necessary
typedef long double LongReal;
#else
#ifdef BUILD_ARBITRARY_PRECISION

typedef mpfr::mpreal Real;
typedef mpfr::mpreal LongReal;

#endif
#endif
#endif

/// Typedef to make specifying line segments easier
typedef std::pair<Vector3, Vector3> LineSeg3;

/// Typedef to make specifying line segments easier
typedef std::pair<Vector2, Vector2> LineSeg2;

/// reference frame type for reduced-coordinate dynamics computations
enum ReferenceFrameType { eGlobal, eLink };

/// Vector2 smart pointer
typedef boost::shared_ptr<Vector2> Vector2Ptr;

/// Vector2 const smart pointer
typedef boost::shared_ptr<const Vector2> Vector2ConstPtr;

/// Vector3 smart pointer
typedef boost::shared_ptr<Vector3> Vector3Ptr;

/// Vector3 const smart pointer
typedef boost::shared_ptr<const Vector3> Vector3ConstPtr;

/// VectorN smart pointer
typedef boost::shared_ptr<VectorN> VectorNPtr;

/// VectorN const smart pointer
typedef boost::shared_ptr<const VectorN> VectorNConstPtr;

/// SVector6 smart pointer
typedef boost::shared_ptr<SVector6> SVector6Ptr;

/// Constant Matrix3 smart pointer
typedef boost::shared_ptr<const Matrix3> Matrix3ConstPtr;

/// Matrix3 smart pointer
typedef boost::shared_ptr<Matrix3> Matrix3Ptr;

/// Matrix4 smart pointer
typedef boost::shared_ptr<Matrix4> Matrix4Ptr;

/// MatrixN smart pointer
typedef boost::shared_ptr<MatrixN> MatrixNPtr;

/// MatrixNN smart pointer
typedef boost::shared_ptr<MatrixN> MatrixNPtr;

/// MatrixNN smart pointer
typedef boost::shared_ptr<MatrixNN> MatrixNNPtr;

/// const MatrixNN smart pointer
typedef boost::shared_ptr<const MatrixNN> MatrixNNConstPtr;

/// constant MatrixN smart pointer
typedef boost::shared_ptr<const MatrixN> MatrixNConstPtr;

/// SMatrix6N smart pointer
typedef boost::shared_ptr<SMatrix6N> SMatrix6NPtr;

/// const SMatrix6N smart pointer
typedef boost::shared_ptr<const SMatrix6N> SMatrix6NConstPtr;

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

