/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_TRIANGLE_H
#define _MOBY_TRIANGLE_H

#include <ostream>
#include <Ravelin/Vector2d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/LinAlgd.h>
#include <Moby/Types.h>

namespace Moby {

/// Defines a triangle for use in collision detection and visualization geometries
class Triangle
{
  public:
    enum FeatureType { eNone, eVertexA, eVertexB, eVertexC, eEdgeAB, eEdgeBC, eEdgeAC, eFace };
    Triangle(const Point3d& a, const Point3d& b, const Point3d& c);
    Triangle() { }
    ~Triangle();
    Triangle(const Triangle& t) { operator=(t); }
    FeatureType determine_feature(const Point3d& point) const;
    FeatureType determine_feature(double s, double t) const;
    void operator=(const Triangle& t);
    static double dist_sq(const Triangle& t1, const Triangle& t2);
    static Triangle transform(const Triangle& t, const Ravelin::Transform3d& Tx);
    const Point3d& get_vertex(unsigned i) const;
    Ravelin::Vector3d calc_normal() const;
    void to_vrml(std::ostream& o, Point3d diffuse_color = Point3d(1,1,1), bool wireframe = false) const;
    double calc_area() const;
    static double calc_signed_dist(const Triangle& tri, const Point3d& point, Point3d& closest_point);
    double calc_signed_dist(const Point3d& point, Point3d& closest_point) const;
    static double calc_sq_dist(const Triangle& tri, const Point3d& point, Point3d& closest_point);
    static double calc_sq_dist(const Triangle& tri, const Point3d& point, double& s, double& t);
    static double calc_sq_dist(const Triangle& tri, const Point3d& origin, const Ravelin::Vector3d& dir, Point3d& cp_seg, Point3d& cp_tri, double& t_line);
    static double calc_sq_dist(const Triangle& tri, const LineSeg3& seg, Point3d& cp_tri, Point3d& cp_seg);
    static double calc_sq_dist(const Triangle& t1, const Triangle& t2, Point3d& cp1, Point3d& cp2);
    std::pair<Point3d, Point3d> calc_AABB() const;
    static void determine_barycentric_coords(const Point2d tri[3], const Point2d& v, double& s, double& t);
    void determine_barycentric_coords(const Point3d& v, double& s, double& t) const;
    double calc_signed_dist(const Point3d& p) const;

    /// Calculates a point in barycentric coordinates
    Point3d calc_point(double s, double t) const { return a*s + b*t + c*(1.0-s-t); }

    /// Calculates the centroid of the triangle
    Point3d calc_centroid() const { return (a + b + c)*0.33333333; }

    /// Gets the value 'd' in the equation n[X]*x + n[Y]*y + n[Z]*z - d = 0, where n is the normal of the triangle, and [x y z] is a point on the triangle
    double calc_offset(const Ravelin::Vector3d& normal) const { return normal.dot(a); }

    /// The first vertex of the triangle
    Point3d a;
    
    /// The second vertex of the triangle
    Point3d b;
 
    /// The third vertex of the triangle
    Point3d c;

  private:

    void determine_params(const Point3d& x, double& s, double& t) const;
    static double calc_sq_dist(const Point3d& origin, const Ravelin::Vector3d& dir, const LineSeg3& seg, Point3d& cp_line, Point3d& cp_seg, double& t_line, double& t_seg);
    static Ravelin::LinAlgd _LA;
}; // end class


std::ostream& operator<<(std::ostream& o, const Triangle& t);

} // end namespace Moby

#endif
