/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_BSPHERE_H_
#define _MOBY_BSPHERE_H_

#include <Moby/BV.h>
#include <Moby/DummyBV.h>

namespace Moby {

/// A sphere used for bounding geometry
class BoundingSphere : public BV 
{
  public:
    BoundingSphere();
    BoundingSphere(const BoundingSphere& bsphere) { operator=(bsphere); }
    BoundingSphere(const Point3d& center, double radius);
    BoundingSphere& operator=(const BoundingSphere& bsphere);
    virtual void transform(const Ravelin::Transform3d& T, BV* result) const;
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const;
    static double calc_dist(const BoundingSphere& s1, const BoundingSphere& s2);
    static bool intersects(const BoundingSphere& a, const BoundingSphere& b);
    static bool intersects(const BoundingSphere& a, const BoundingSphere& b, const Ravelin::Transform3d& aTb);
    static bool intersects(const BoundingSphere& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q);
    static bool outside(const BoundingSphere& a, const Point3d& point, double tol = NEAR_ZERO);
    boost::shared_ptr<const BoundingSphere> get_this() const { return boost::dynamic_pointer_cast<const BoundingSphere>(shared_from_this()); }
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const { return BoundingSphere::outside(*this, point, tol); }
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const { return BoundingSphere::intersects(*this, seg, tmin, tmax, q); }

    template <class ForwardIterator>
    BoundingSphere(ForwardIterator begin, ForwardIterator end);

    /// Gets the pose for this BV
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const { return center.pose; }

    /// Gets the lower bounds on the bounding sphere
    virtual Point3d get_lower_bounds() const { return Point3d(center[0]-radius, center[1]-radius, center[2] - radius, center.pose); }

    /// Gets the upper bounds on the bounding sphere
    virtual Point3d get_upper_bounds() const { return Point3d(center[0]+radius, center[1]+radius, center[2] + radius, center.pose); }

    /// Center of the bounding box
    Point3d center;

    /// The radius of the bounding sphere (we use a float b/c accuracy here not so important)
    double radius;

    /// Calculates the volume of this bounding volume
    virtual double calc_volume() const { return M_PI * radius * radius * radius; }
}; // end class

// inline functions
#include "BoundingSphere.inl"

} // end namespace

#endif

