/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_AABB_H_
#define _MOBY_AABB_H_

#include <Moby/OBB.h>
#include <Moby/BV.h>

namespace Moby {

/// An axis-aligned bounding box 
class AABB : public BV
{
  public:
    AABB() {}

    template <class InputIterator>
    AABB(InputIterator begin, InputIterator end);

    virtual void transform(const Ravelin::Transform3d& T, BV* result) const;
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const { return AABB::outside(*this, point, tol); }
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const { return AABB::intersects(*this, seg, tmin, tmax, q); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const;
    virtual double calc_volume() const;
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const { return minp.pose; }
    virtual Point3d get_lower_bounds() const;
    virtual Point3d get_upper_bounds() const;
    OBB get_OBB() const;
    static bool intersects(const AABB& a, const AABB& b);
    static bool intersects(const AABB& a, const AABB& b, const Ravelin::Transform3d& aTb);
    static bool outside(const AABB& a, const Point3d& point, double tol = NEAR_ZERO);
    static bool intersects(const AABB& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q);
    static void get_closest_point(const AABB& a, const Point3d& p, Point3d& closest);
    static double get_farthest_point(const AABB& a, const Point3d& p, Point3d& farthest);
    AABB& operator=(const AABB& a);

    /// The lower corner of the AABB
    Point3d minp;

    /// The upper corner of the AABB;
    Point3d maxp;
}; // end class

#include "AABB.inl"

} // end namespace

#endif

