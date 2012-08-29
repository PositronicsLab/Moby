/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
    ThickTriangle(const Triangle& t, Real epsilon) { construct_from_triangle(t, epsilon); }
    void operator=(const ThickTriangle& t) { _planes = t._planes; _normal = t._normal; tri = t.tri; }
    bool intersect_seg(const LineSeg3& seg, Real& tnear, Vector3& isect) const;
    bool intersect_seg(const LineSeg3& seg, Real& tnear, Vector3& isect, Vector3& normal) const;
    bool point_inside(const Vector3& point) const;
    Vector3 determine_normal(const Vector3& p) const;
    const std::list<Plane>& planes() const { return _planes; }
    void construct_from_triangle(const Triangle& t, Real epsilon);

    /// The original triangle  
    Triangle tri;

  private:
    std::list<Plane> _planes;
    Vector3 _normal;  
}; // end class

inline std::ostream& operator<<(std::ostream& out, const ThickTriangle& t)
{
  BOOST_FOREACH(const Plane& p, t.planes())
    std::cout << "plane: " << p << std::endl;

  return out;
}

} // end namespace

#endif

