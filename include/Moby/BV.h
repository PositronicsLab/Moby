/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_BV_H_
#define _MOBY_BV_H_

#include <stack>
#include <iostream>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/foreach.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>

namespace Moby {

class AABB;
class OBB;
class BoundingSphere;
class SSR;
class SSL;

/// An abstract bounding volume 
/**
 * \note the BV is generally constructed such that its frame is aligned with
 *       that of the underlying rigid body (or collision geometry).  
 *       Therefore, the center of the BV is computed relative to the 
 *       center-of-mass of the body (or the center of the geometry).  The
 *       orientation of the BV will always be identical to the orientation of
 *       the rigid body (or collision geometry).
 */
class BV : public boost::enable_shared_from_this<BV>
{
  public:
    virtual ~BV() {}

    /// Virtual function for outputting the bounding volume to VRML
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const = 0
;

    /// Determines whether a point is outside the bounding volume
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const = 0;

    /// Determines whether a line segment intersects the bounding volume
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const = 0;

    /// Gets the associated pose for this bounding volume
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const = 0;

    /// Virtual function for transforming the BV
    virtual void transform(const Ravelin::Transform3d& T, BV* result) const = 0;

    /// Virtual function that calculates a velocity-expanded BV
    /**
     * \param g the geometry that this bounding volume represents
     * \param v the "velocity" to sweep by
     * \return the velocity-expanded bounding volume
     */ 
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const = 0;

    /// Convenience method
    static bool intersects(BVPtr a, BVPtr b) { return intersects(a.get(), b.get()); }

    /// Convenience method
    static bool intersects(BVPtr a, BVPtr b, const Ravelin::Transform3d& T) { return intersects(a.get(), b.get(), T); }

    /// Convenience method
    static double calc_distance(BVPtr a, BVPtr b, Point3d& cp1, Point3d& cp2) { return calc_distance(a.get(), b.get(), cp1, cp2); }

    /// Convenience method
    static double calc_distance(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, Point3d& cp1, Point3d& cp2) { return calc_distance(a.get(), b.get(), aTb, cp1, cp2); }

    static bool intersects(const BV* a, const BV* b);
    static bool intersects(const BV* a, const BV* b, const Ravelin::Transform3d& T);
    static double calc_distance(const BV* a, const BV* b, Point3d& cp1, Point3d& cp2);
    static double calc_distance(const BV* a, const BV* b, const Ravelin::Transform3d& aTb, Point3d& cp1, Point3d& cp2);

    BVPtr get_this() { return boost::dynamic_pointer_cast<BV>(shared_from_this()); }
    boost::shared_ptr<const BV> get_this() const { return boost::dynamic_pointer_cast<const BV>(shared_from_this()); }
    bool is_leaf() const { return children.empty(); }

    template <class OutputIterator>
    OutputIterator get_all_BVs(OutputIterator begin) const;

    template <class OutputIterator>
    OutputIterator get_all_leafs(OutputIterator begin) const;

    template <class OutputIterator>
    static OutputIterator intersect_BV_trees(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, const Ravelin::Transform3d& bTa, OutputIterator output_begin);

    /// Userdata for the BV
    boost::shared_ptr<void> userdata;

    /// The collision geometry associated with this bounding volume
    CollisionGeometryPtr geom;

    /// The children of this BV
    std::list<BVPtr> children;

    /// Gets the volume for this bounding volume
    virtual double calc_volume() const = 0;

    /// Gets the lower bound on a AABB around the bounding volume when a transform of T is applied
    virtual Point3d get_lower_bounds() const = 0;

    /// Gets the upper bound on a AABB around the bounding volume when a transform of T is applied
    virtual Point3d get_upper_bounds() const = 0;

  private:

    static bool intersects(const OBB* O, const BoundingSphere* S);
    static bool intersects(const OBB* O, const BoundingSphere* S, const Ravelin::Transform3d& OTS);
    static bool intersects(const OBB* O, const AABB* A);
    static bool intersects(const OBB* O, const AABB* A, const Ravelin::Transform3d& OTA);
    static bool intersects(const OBB* O, const SSR* S);
    static bool intersects(const OBB* O, const SSR* S, const Ravelin::Transform3d& OTS);
    static bool intersects(const OBB* O, const SSL* S);
    static bool intersects(const OBB* O, const SSL* S, const Ravelin::Transform3d& OTS);
    static bool intersects(const AABB* A, const BoundingSphere* S);
    static bool intersects(const AABB* A, const BoundingSphere* S, const Ravelin::Transform3d& ATS);
    static bool intersects(const AABB* A, const SSL* S);
    static bool intersects(const AABB* A, const SSL* S, const Ravelin::Transform3d& ATS);
    static bool intersects(const SSR* S, const AABB* A);
    static bool intersects(const SSR* S, const AABB* A, const Ravelin::Transform3d& STA);
    static bool intersects(const SSR* S, const BoundingSphere* B);
    static bool intersects(const SSR* S, const BoundingSphere* B, const Ravelin::Transform3d& STB);
    static bool intersects(const SSR* S, const SSL* B);
    static bool intersects(const SSR* S, const SSL* B, const Ravelin::Transform3d& STB);
    static bool intersects(const SSL* S, const BoundingSphere* B);
    static bool intersects(const SSL* S, const BoundingSphere* B, const Ravelin::Transform3d& STB);
}; // end class

// include inline functions
#include "BV.inl"

} // end namespace

#endif

