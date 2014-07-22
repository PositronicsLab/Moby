/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _THICK_TRIANGLE_H_
#define _THICK_TRIANGLE_H_

#include <iostream>
#include <list>
#include <cmath>
#include <boost/foreach.hpp>
#include <Moby/Triangle.h>
#include <Moby/Plane.h>

namespace Moby {

/// Thick triangle for intersection / query testing
class ThickTriangle
{
  public:
    ThickTriangle(const Triangle& t, double epsilon) { construct_from_triangle(t, epsilon); }
    void operator=(const ThickTriangle& t) { _planes = t._planes; _normal = t._normal; tri = t.tri; }
    bool intersect_seg(const LineSeg3& seg, double& tnear, Point3d& isect) const;
    bool intersect_seg(const LineSeg3& seg, double& tnear, Point3d& isect, Ravelin::Vector3d& normal) const;
    bool point_inside(const Point3d& point) const;
    Ravelin::Vector3d determine_normal(const Point3d& p) const;
    const std::list<Plane>& planes() const { return _planes; }
    void construct_from_triangle(const Triangle& t, double epsilon);
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const;

    /// The original triangle  
    Triangle tri;

  private:
    std::list<Plane> _planes;
    Ravelin::Vector3d _normal;  
}; // end class

inline std::ostream& operator<<(std::ostream& out, const ThickTriangle& t)
{
  BOOST_FOREACH(const Plane& p, t.planes())
    std::cout << "plane: " << p << std::endl;

  return out;
}

} // end namespace

#endif

