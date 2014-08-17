/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _COMPGEOM_H
#define _COMPGEOM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/Origin2d.h>
#include <Ravelin/Vector2d.h>
#include <Ravelin/Matrix2d.h>
#include <Ravelin/Origin3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Moby/Constants.h>
#include <Moby/Triangle.h>
#include <Moby/FastThreadable.h>
#include <Moby/Polyhedron.h>
#include <Moby/Log.h>
#include <Moby/NumericalException.h>
#include <Moby/Types.h>

#ifdef THREADSAFE
#include <pthread.h>
#endif

// needed for qhull
extern "C"
{
  #include <qhull/qhull_a.h>
}

namespace Moby {

class Polyhedron;  

/// Class for performing assorted computational geometry functions
class CompGeom
{
  template <class T, class V> friend class CompGeomSpecOne;
  template <class T, class U, class V> friend class CompGeomSpecTwo;

  public:
    enum SegSegIntersectType { eSegSegNoIntersect, eSegSegIntersect, eSegSegVertex, eSegSegEdge };
    enum LineLineIntersectType { eLineLineNoIntersect, eLineLineIntersect, eLineLineVertex, eLineLineEdge };
    enum SegTriIntersectType { eSegTriNoIntersect, eSegTriInclFace, eSegTriInclVertex, eSegTriInclEdge, eSegTriInside, eSegTriVertex, eSegTriEdge, eSegTriFace, eSegTriEdgeOverlap, eSegTriPlanarIntersect };
    enum SegPlaneIntersectType { eSegPlaneInPlane, eSegPlaneToSide, eSegPlaneFirst, eSegPlaneSecond, eSegPlaneOtherIntersect };
    enum PolygonLocationType { ePolygonInside, ePolygonOutside, ePolygonOnVertex, ePolygonOnEdge };
    enum SegLocationType { eSegInterior, eSegOrigin, eSegEndpoint, eSegOutside };
    enum OrientationType { eLeft, eOn, eRight };
    enum VisibilityType { eVisible, eInvisible, eCoplanar };

    static SegLocationType determine_seg_location(const LineSeg3& seg, double t);
    static SegLocationType determine_seg_location(const LineSeg2& seg, double t);
    static unsigned get_num_intersects(SegSegIntersectType t);
    static unsigned get_num_intersects(SegTriIntersectType t);
    static unsigned get_num_intersects(LineLineIntersectType t);
    static double determine_3D_to_2D_offset(const Ravelin::Origin3d& o, const Ravelin::Matrix3d& R);
    static Ravelin::Matrix3d calc_3D_to_2D_matrix(const Ravelin::Vector3d& normal);
    static Point3d generate_point_on_plane(const Ravelin::Vector3d& normal, double d);
    static bool intersect_noncoplanar_tris(const Triangle& t1, const Triangle& t2, Point3d& p1, Point3d& p2);
    static SegTriIntersectType intersect_seg_tri(const LineSeg3& seg, const Triangle& t, Point3d& isect, Point3d& isect2, double tol = NEAR_ZERO);
    static SegTriIntersectType intersect_seg_tri(const LineSeg2& seg, const Point2d tri[3], Point2d& isect, Point2d& isect2, double tol = NEAR_ZERO);
    static SegSegIntersectType intersect_segs(const LineSeg3& s1, const LineSeg3& s2, Point3d& isect, Point3d& isect2);
    static double intersect_seg_plane(const Ravelin::Vector3d& nrm, double d, const LineSeg3& seg);
    static SegPlaneIntersectType intersect_seg_plane(const Triangle& tri, const LineSeg3& seg, Point3d& isect, unsigned int& big_dim, double tol = NEAR_ZERO);
    static SegSegIntersectType intersect_segs(const LineSeg2& s1, const LineSeg2& s2, Point2d& isect, Point2d& isect2);
    static LineLineIntersectType intersect_lines(const Point2d& o1, const Ravelin::Vector2d& dir1, double min1, double max1, const Point2d& o2, const Ravelin::Vector2d& dir2, double min2, double max2, double& s, double& t);
    static bool collinear(const Point3d& a, const Point3d& b, const Point3d& c, double tol = NEAR_ZERO);
    static PolygonLocationType in_tri(const Point2d t[3], const Point2d& p, double tol = NEAR_ZERO);
    static bool collinear(const Point2d& a, const Point2d& b, const Point2d& c, double tol = NEAR_ZERO);
    static long double volume(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d);
    static long double volume(const Ravelin::Origin3d& a, const Ravelin::Origin3d& b, const Ravelin::Origin3d& c, const Ravelin::Origin3d& d);
    static VisibilityType volume_sign(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d, double tol = NEAR_ZERO);
    static VisibilityType volume_sign(const Ravelin::Origin3d& a, const Ravelin::Origin3d& b, const Ravelin::Origin3d& c, const Ravelin::Origin3d& d, double tol = NEAR_ZERO);
    static bool coplanar(const Triangle& t1, const Triangle& t2, double tol = NEAR_ZERO);
    static Ravelin::Origin3d to_3D(const Ravelin::Origin2d& p, const Ravelin::Matrix3d& RT);
    static Ravelin::Origin3d to_3D(const Ravelin::Origin2d& p, const Ravelin::Matrix3d& RT, double offset);
    static Ravelin::Origin2d to_2D(const Point3d& p, const Ravelin::Matrix3d& R);
    static bool point_in_tri(const Point3d& p, const Triangle& t, double tol = NEAR_ZERO);
    static bool coplanar(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d, double tol = NEAR_ZERO);
    static double calc_dist(const LineSeg3& line_seg, const Point3d& point, double& t, Point3d& closest);
    static double calc_closest_points(const LineSeg3& s1, const LineSeg3& s2, Point3d& p1, Point3d& p2); 
    static bool query_intersect_tri_tri(const Triangle& t1, const Triangle& t2);
    static PolygonLocationType in_tri(const Triangle& t, const Point3d& p, double tol = NEAR_ZERO);

    /// Returns a if b > 0, -a if b < 0, and 0 if b = 0
    static double sgn(double a, double b, double tol = NEAR_ZERO) { if (b > NEAR_ZERO) return a; else if (b < -NEAR_ZERO) return -a; else return 0.0; }

    /// Determines whether two numbers are relatively/absolutely equal to within some tolerance
    /**
     * If the numbers are both smaller than zero, an absolute tolerance is used;
     * otherwise, a relative tolerance is used.
     * \note rationale for why this test is useful: [Ericson, 2005] p. 443
     */ 
    static bool rel_equal(double x, double y, double tol = NEAR_ZERO) { return std::fabs(x-y) <= tol*std::max(std::fabs(x), std::max(std::fabs(y), (double) 1.0)); }

    template <class T, class U>
    static double determine_line_param(const T& p, const U& dir, const T& v);

    template <class OutputIterator>
    static OutputIterator intersect_seg_tri(const LineSeg2& seg, const Point2d tri[3], OutputIterator output_begin);

    template <class OutputIterator>
    static OutputIterator intersect_segs(const LineSeg2& s1, const LineSeg2& s2, OutputIterator output_begin);

     template <class OutputIterator>
     static OutputIterator intersect_tris(const Point2d t1[3], const Point2d t2[3], OutputIterator output_begin);

    template <class ForwardIterator>
    static void calc_min_area_bounding_rect(ForwardIterator begin, ForwardIterator end, Point2d& p1, Point2d& p2, Point2d& p3, Point2d& p4);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator intersect_polygons(ForwardIterator pbegin, ForwardIterator pend, ForwardIterator qbegin, ForwardIterator qend, const Ravelin::Vector3d& normal, OutputIterator isect_begin);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator intersect_polygons(ForwardIterator pbegin, ForwardIterator pend, ForwardIterator qbegin, ForwardIterator qend, OutputIterator isect_begin);

    template <class ForwardIterator>
    static bool ccw(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol = NEAR_ZERO);

    template <class ForwardIterator>
    static bool ccw(ForwardIterator begin, ForwardIterator end, double tol = NEAR_ZERO);

    template <class ForwardIterator>
    static bool is_convex_polygon(ForwardIterator begin, ForwardIterator end, double tol = NEAR_ZERO);

    template <class ForwardIterator>
    static bool is_convex_polygon(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol = NEAR_ZERO);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, const Ravelin::Vector3d& normal, OutputIterator target_begin);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin);

    template <class ForwardIterator>
    static PolygonLocationType polygon_location(ForwardIterator begin, ForwardIterator end, const Point2d& point);

    template <class ForwardIterator>
    static Point2d calc_centroid_2D(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static Point3d calc_centroid_2D(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal);

    template <class ForwardIterator>
    static double calc_polygon_area(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal);  

    template <class ForwardIterator>
    static double calc_polygon_area(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void project_plane(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double offset);

    template <class ForwardIterator>
    static double fit_plane(ForwardIterator begin, ForwardIterator end, Ravelin::Vector3d& normal, double& offset);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator triangulate_convex_polygon(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin);

    template <class OutputIterator>
    static OutputIterator intersect_tris(const Triangle& t1, const Triangle& t2, OutputIterator begin);

    template <class OutputIterator>
    static OutputIterator intersect_coplanar_tris(const Triangle& t1, const Triangle& t2, const Ravelin::Vector3d& normal, OutputIterator begin);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator to_2D(ForwardIterator begin, ForwardIterator end, OutputIterator begin_target, const Ravelin::Matrix3d& R);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator to_3D(ForwardIterator begin, ForwardIterator end, OutputIterator begin_target, const Ravelin::Matrix3d& RT, double offset);

    template <class ForwardIterator>
    static void determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point3d, Point3d>& endpoints);

    template <class ForwardIterator>
    static void determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point3d*, Point3d*>& endpoints);

    template <class ForwardIterator>
    static void determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point2d, Point2d>& endpoints);

    template <class ForwardIterator>
    static void determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point2d*, Point2d*>& endpoints);

    template <class ForwardIterator>
    static PolyhedronPtr calc_convex_hull(ForwardIterator first, ForwardIterator last);

    template <class ForwardIterator>
    static Point3d calc_centroid_3D(ForwardIterator first, ForwardIterator last);
    
    template <class ForwardIterator>
    static PolyhedronPtr calc_hs_intersection(ForwardIterator start, ForwardIterator end, const Ravelin::VectorNd& interior_point);

    template <class ForwardIterator>
    static double find_hs_interior_point(ForwardIterator start, ForwardIterator end, Point3d& point);

    template <class ForwardIterator>
    static unsigned calc_dimensionality(ForwardIterator first, ForwardIterator last, double tol = NEAR_ZERO);

    template <class ForwardIterator>
    static bool intersect_seg_convex_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, double& te, double& tl, double tol = NEAR_ZERO);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator intersect_seg_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, OutputIterator outbegin);

    template <class BidirectionalIterator, class OutputIterator>
    static OutputIterator triangulate_polygon_2D(BidirectionalIterator begin, BidirectionalIterator end, OutputIterator begintarget, double tol = NEAR_ZERO);

    /// Returns the point a*t + (1-t)*b on the line segment ab
    template <class T>
    static T point_on_seg(T a, T b, double t) { return b + (a - b)*t; }

  private:
    static void normalize_or_zero(Ravelin::Vector3d& v);
    static void update_box(const Point2d& lP, const Point2d& rP, const Point2d& bP, const Point2d& tP, const Ravelin::Vector2d& U, const Ravelin::Vector2d& V, double& min_area_div4, Point2d& center, Ravelin::Vector2d axis[2], double extent[2]);
    static void project_onto_axis(const Triangle& t, const Ravelin::Vector3d& v, double& rfmin, double& rfmax);
    static PolygonLocationType in_tri(OrientationType ori1, OrientationType ori2, OrientationType ori3);
    static bool intersect(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d, double tol = NEAR_ZERO);
    static bool intersect_prop(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d, double tol = NEAR_ZERO);
    static void clip_convex_polygon_against_line(const Ravelin::Vector2d& rkN, double fC, unsigned& riQuantity, Point2d akV[6]);

    template <class BidirectionalIterator>
    static bool in_cone(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol = NEAR_ZERO);

    template <class BidirectionalIterator>
    static bool diagonal(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol = NEAR_ZERO);

    template <class BidirectionalIterator>
    static bool diagonalie(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol = NEAR_ZERO);

    template <class OutputIterator>
    static unsigned advance(unsigned a, unsigned* aa, unsigned n, bool inside, const Point2d& p, OutputIterator& current);

    static bool same_side(const Point3d& p1, const Point3d& p2, const Point3d& a, const Point3d& b, double tol = std::sqrt(std::numeric_limits<double>::epsilon()));
    static bool compute_intervals_isectline(const Triangle& t1, double vv0, double vv1, double vv2, double d0, double d1, double d2, double d0d1, double d0d2, double& isect0, double& isect1, Point3d& isectpoint0, Point3d& isectpoint1);
    static void isect2X(const Point3d& vtx0, const Point3d& vtx1, const Point3d& vtx2, double vv0, double vv1, double vv2, double d0, double d1, double d2, double& isect0, double& isect1, Point3d& isectpoint0, Point3d& isectpoint1);
    static SegSegIntersectType get_parallel_intersect_type(const LineSeg2& s1, const LineSeg2& s2, Point2d& isect, Point2d& isect2);
    static bool between(const Point2d& a, const Point2d& b, const Point2d& c, double tol = NEAR_ZERO);
    static SegTriIntersectType seg_tri_cross(const LineSeg3& seg, const Triangle& t, double tol = NEAR_ZERO);
    static SegTriIntersectType intersect_seg_tri_in_plane(const LineSeg3& seg, const Triangle& t, Point3d& isect, Point3d& isect2, double tol = NEAR_ZERO);
    static PolygonLocationType in_tri(const Triangle& t, const Point3d& p, unsigned int skip_dim, double tol = NEAR_ZERO);
    static void find_plane_coeff(const Triangle& t, Ravelin::Vector3d& N, unsigned int& m, double& dot);
    static Ravelin::Vector3d normal_vec(const Point3d& a, const Point3d& b, const Point3d& c);
    static OrientationType area_sign(const Point2d& a, const Point2d& b, const Point2d& c, double tol = NEAR_ZERO);
    static long double area(const Point2d& a, const Point2d& b, const Point2d& c);
    static bool test_edge_edge(const Point3d& p, const Point3d& a, const Point3d& b, double Ax, double Ay, unsigned i0, unsigned i1);
    static bool test_edge_tri(const Point3d& a, const Point3d& b, const Triangle& t, unsigned i0, unsigned i1);
    static bool test_point_in_tri(const Point3d& p, const Triangle& t, unsigned i0, unsigned i1);
    static bool test_coplanar_tri_tri(const Ravelin::Vector3d& N, const Triangle& t1, const Triangle& t2);

    #ifdef THREADSAFE
    static pthread_mutex_t _qhull_mutex;
    #endif
}; // end class

/// Specialized one-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class V>
class CompGeomSpecOne
{
};

template <class T>
class CompGeomSpecOne<T, Point3d>
{
  friend class CompGeom;

  private:
    static unsigned calc_dimensionality(T begin, T end, double tol);
    static PolyhedronPtr calc_convex_hull(T begin, T end); 
    static bool is_convex_polygon(T begin, T end, const Ravelin::Vector3d& normal, double tol);
    static bool ccw(T begin, T end, const Ravelin::Vector3d& normal, double tol = NEAR_ZERO);
    static double fit_plane(T begin, T end, Ravelin::Vector3d& normal, double& offset);
};

template <class T>
class CompGeomSpecOne<T, Point3d*>
{
  friend class CompGeom;

  private:
    static unsigned calc_dimensionality(T begin, T end, double tol); 
    static PolyhedronPtr calc_convex_hull(T begin, T end);
    static bool is_convex_polygon(T begin, T end, const Ravelin::Vector3d& normal, double tol);
    static bool ccw(T begin, T end, const Ravelin::Vector3d& normal, double tol = NEAR_ZERO);
    static double fit_plane(T begin, T end, Ravelin::Vector3d& normal, double& offset);
};

template <class T>
class CompGeomSpecOne<T, Point2d>
{
  friend class CompGeom;

  private:
    static bool intersect_seg_convex_polygon(T begin, T end, const LineSeg2& seg, double& te, double& tl, double tol);
    static unsigned calc_dimensionality(T begin, T end, double tol) { }
    static bool is_convex_polygon(T begin, T end, double tol);
    static void calc_min_area_bounding_rect(T begin, T end, Point2d& p1, Point2d& p2, Point2d& p3, Point2d& p4); 
    static bool ccw(T begin, T end, double tol = NEAR_ZERO);
};

template <class T>
class CompGeomSpecOne<T, Point2d*>
{
  friend class CompGeom;

  private:
    static bool intersect_seg_convex_polygon(T begin, T end, const LineSeg2& seg, double& te, double& tl, double tol);
    static bool is_convex_polygon(T begin, T end, double tol);
    static void calc_min_area_bounding_rect(T begin, T end, Point2d& p1, Point2d& p2, Point2d& p3, Point2d& p4); 
    static bool ccw(T begin, T end, double tol = NEAR_ZERO);
};

/// Specialized two-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class U, class V>
class CompGeomSpecTwo
{
};

/// Specialized two-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class U>
class CompGeomSpecTwo<T, U, Point2d*>
{
  friend class CompGeom;

  private:
    static U calc_convex_hull(T source_begin, T source_end, U target_begin);
    static U intersect_seg_polygon(T begin, T end, const LineSeg2& seg, U outbegin);
    static U to_3D(T begin, T end, U begin_target, const Ravelin::Matrix3d& RT, double offset);
};

/// Specialized two-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class U>
class CompGeomSpecTwo<T, U, Point2d>
{
  friend class CompGeom;

  private:
    static U calc_convex_hull(T source_begin, T source_end, U target_begin);
    static U intersect_seg_polygon(T begin, T end, const LineSeg2& seg, U outbegin);
    static U to_3D(T begin, T end, U begin_target, const Ravelin::Matrix3d& RT, double offset);
};

/// Specialized two-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class U>
class CompGeomSpecTwo<T, U, Point3d>
{
  friend class CompGeom;

  private:
    static U calc_convex_hull(T source_begin, T source_end, U target_begin);
    static U calc_convex_hull(T source_begin, T source_end, const Ravelin::Vector3d& normal, U target_begin);
    static U to_2D(T begin_source, T end_source, U begin_target, const Ravelin::Matrix3d& R);
};

/// Specialized two-parameter computational geometry routines
/**
 * This hack is necessary because C++ function templates do not support
 * partial specialization.
 */
template <class T, class U>
class CompGeomSpecTwo<T, U, Point3d*>
{
  friend class CompGeom;

  private:
    static U calc_convex_hull(T source_begin, T source_end, U target_begin);
    static U calc_convex_hull(T source_begin, T source_end, const Ravelin::Vector3d& normal, U target_begin);
    static U to_2D(T begin_source, T end_source, U begin_target, const Ravelin::Matrix3d& R);
};

// include inline functions
#include "CompGeom.inl"

} // end namespace

#endif
