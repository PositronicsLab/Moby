/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_TRIANGLE_H
#define _MOBY_TRIANGLE_H

#include <ostream>
#include <Moby/Vector2.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix4.h>

namespace Moby {

/// Defines a triangle for use in collision detection and visualization geometries
class Triangle
{
  public:
    enum FeatureType { eNone, eVertexA, eVertexB, eVertexC, eEdgeAB, eEdgeBC, eEdgeAC, eFace };
    Triangle(const Vector3& a, const Vector3& b, const Vector3& c);
    Triangle() { }
    ~Triangle();
    Triangle(const Triangle& t) { operator=(t); }
    FeatureType determine_feature(const Vector3& point) const;
    FeatureType determine_feature(Real s, Real t) const;
    void operator=(const Triangle& t);
    static Real dist_sq(const Triangle& t1, const Triangle& t2);
    static Triangle transform(const Triangle& t, const Matrix4& m);
    const Vector3& get_vertex(unsigned i) const;
    Vector3 calc_normal() const;
    void to_vrml(std::ostream& o, Vector3 diffuse_color = Vector3(1,1,1), bool wireframe = false) const;
    Real calc_area() const;
    static Real calc_sq_dist(const Triangle& tri, const Vector3& point, Vector3& closest_point);
    static Real calc_sq_dist(const Triangle& tri, const Vector3& origin, const Vector3& dir, Vector3& cp_seg, Vector3& cp_tri, Real& t_line);
    static Real calc_sq_dist(const Triangle& tri, const LineSeg3& seg, Vector3& cp_tri, Vector3& cp_seg);
    static Real calc_sq_dist(const Triangle& t1, const Triangle& t2, Vector3& cp1, Vector3& cp2);
    std::pair<Vector3, Vector3> calc_AABB() const;
    static void determine_barycentric_coords(const Vector2 tri[3], const Vector2& v, Real& s, Real& t);
    void determine_barycentric_coords(const Vector3& v, Real& s, Real& t) const;
    Real calc_signed_dist(const Vector3& p) const;

    /// Calculates a point in barycentric coordinates
    Vector3 calc_point(Real s, Real t) const { return a*s + b*t + c*(1.0-s-t); }

    /// Calculates the centroid of the triangle
    Vector3 calc_centroid() const { return (a + b + c)*0.33333333; }

    /// Gets the value 'd' in the equation n[X]*x + n[Y]*y + n[Z]*z - d = 0, where n is the normal of the triangle, and [x y z] is a point on the triangle
    Real calc_offset(const Vector3& normal) const { return normal.dot(a); }

    /// The first vertex of the triangle
    Vector3 a;
    
    /// The second vertex of the triangle
    Vector3 b;
 
    /// The third vertex of the triangle
    Vector3 c;

  private:

    void determine_params(const Vector3& x, Real& s, Real& t) const;
    static Real calc_sq_dist(const Vector3& origin, const Vector3& dir, const LineSeg3& seg, Vector3& cp_line, Vector3& cp_seg, Real& t_line, Real& t_seg);
}; // end class


std::ostream& operator<<(std::ostream& o, const Triangle& t);

} // end namespace Moby

#endif
