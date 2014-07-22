/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_SSL_H_
#define _MOBY_SSL_H_

#include <stack>
#include <iostream>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/AAngled.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>

namespace Moby {

/// A sphere-swept line (SSL) that optionally allows building an SSL tree
class SSL : public BV 
{
  public:
    SSL();
    SSL(const SSL& obb) { operator=(obb); }
    SSL(const Point3d& p1, const Point3d& p2, double radius);
    SSL(const SSL& s, const Ravelin::Vector3d& v);
    SSL& operator=(const SSL& s);
    virtual void transform(const Ravelin::Transform3d& T, BV* result) const;
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const;
    static double calc_dist(const SSL& o, const Point3d& p);
    static double calc_dist(const SSL& a, const SSL& b, Point3d& cpa, Point3d& cpb);
    static double calc_dist(const SSL& a, const SSL& b, const Ravelin::Transform3d& aTb, Point3d& cpa, Point3d& cpb);
    static bool intersects(const SSL& a, const SSL& b);
    static bool intersects(const SSL& a, const SSL& b, const Ravelin::Transform3d& aTb);
    static bool intersects(const SSL& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q);
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const { return SSL::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSL& a, const Point3d& point, double tol = NEAR_ZERO);
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const { return SSL::outside(*this, point, tol); }
    boost::shared_ptr<SSL> get_this() { return boost::dynamic_pointer_cast<SSL>(shared_from_this()); }
    boost::shared_ptr<const SSL> get_this() const { return boost::dynamic_pointer_cast<const SSL>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    unsigned calc_size() const;
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const { return p1.pose; }
    virtual Point3d get_lower_bounds() const;
    virtual Point3d get_upper_bounds() const;

    /// Calculates (approximate?) volume of the SSL 
    virtual double calc_volume() const { return (p1 - p2).norm() * M_PI * radius * radius * radius; }

    /// The first point of the line segment 
    Point3d p1;

    /// The second point of the line segment
    Point3d p2;

    /// Radius of the spherical addition 
    double radius;

  private:
    static double calc_sq_dist(const LineSeg3& seg, const Point3d& q, Point3d& cp);
    static double calc_sq_dist(const LineSeg3& seg0, const LineSeg3& seg1, Point3d& cp0, Point3d& cp1);
    static double calc_sq_dist(const Point3d& origin, const Ravelin::Vector3d& dir, const LineSeg3& seg, Point3d& cp_line, Point3d& cp_seg, double& line_param);

}; // end class

// include inline functions
#include "SSL.inl"

} // end namespace

#endif

