/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _COMPGEOM_H
#define _COMPGEOM_H

#include <pthread.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <Moby/Constants.h>
#include <Moby/Triangle.h>
#include <Moby/FastThreadable.h>
#include <Moby/Vector2.h>
#include <Moby/Matrix2.h>
#include <Moby/Optimization.h>
#include <Moby/Polyhedron.h>
#include <Moby/LinAlg.h>
#include <Moby/Log.h>
#include <Moby/NumericalException.h>
#include <Moby/Types.h>

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
  public:
    enum SegSegIntersectType { eSegSegNoIntersect, eSegSegIntersect, eSegSegVertex, eSegSegEdge };
    enum LineLineIntersectType { eLineLineNoIntersect, eLineLineIntersect, eLineLineVertex, eLineLineEdge };
    enum SegTriIntersectType { eSegTriNoIntersect, eSegTriInclFace, eSegTriInclVertex, eSegTriInclEdge, eSegTriInside, eSegTriVertex, eSegTriEdge, eSegTriFace, eSegTriEdgeOverlap, eSegTriPlanarIntersect };
    enum SegPlaneIntersectType { eSegPlaneInPlane, eSegPlaneToSide, eSegPlaneFirst, eSegPlaneSecond, eSegPlaneOtherIntersect };
    enum PolygonLocationType { ePolygonInside, ePolygonOutside, ePolygonOnVertex, ePolygonOnEdge };
    enum SegLocationType { eSegInterior, eSegOrigin, eSegEndpoint, eSegOutside };
    enum OrientationType { eLeft, eOn, eRight };
    enum VisibilityType { eVisible, eInvisible, eCoplanar };

    static SegLocationType determine_seg_location(const LineSeg3& seg, Real t);
    static SegLocationType determine_seg_location(const LineSeg2& seg, Real t);
    static unsigned get_num_intersects(SegSegIntersectType t);
    static unsigned get_num_intersects(SegTriIntersectType t);
    static unsigned get_num_intersects(LineLineIntersectType t);
    static Real determine_3D_to_2D_offset(const Vector3& point, const Matrix3& R);
    static Matrix3 calc_3D_to_2D_matrix(const Vector3& normal);
    static Vector3 generate_point_on_plane(const Vector3& normal, Real d);
    static bool intersect_noncoplanar_tris(const Triangle& t1, const Triangle& t2, Vector3& p1, Vector3& p2);
    static SegTriIntersectType intersect_seg_tri(const LineSeg3& seg, const Triangle& t, Vector3& isect, Vector3& isect2, Real tol = NEAR_ZERO);
    static SegTriIntersectType intersect_seg_tri(const LineSeg2& seg, const Vector2 tri[3], Vector2& isect, Vector2& isect2, Real tol = NEAR_ZERO);
    static SegSegIntersectType intersect_segs(const LineSeg3& s1, const LineSeg3& s2, Vector3& isect, Vector3& isect2);
    static Real intersect_seg_plane(const Vector3& nrm, Real d, const LineSeg3& seg);
    static SegPlaneIntersectType intersect_seg_plane(const Triangle& tri, const LineSeg3& seg, Vector3& isect, unsigned int& big_dim, Real tol = NEAR_ZERO);
    static SegSegIntersectType intersect_segs(const LineSeg2& s1, const LineSeg2& s2, Vector2& isect, Vector2& isect2);
    static LineLineIntersectType intersect_lines(const Vector2& o1, const Vector2& dir1, Real min1, Real max1, const Vector2& o2, const Vector2& dir2, Real min2, Real max2, Real& s, Real& t);
    static bool collinear(const Vector3& a, const Vector3& b, const Vector3& c, Real tol = NEAR_ZERO);
    static PolygonLocationType in_tri(const Vector2 t[3], const Vector2& p, Real tol = NEAR_ZERO);
    static bool collinear(const Vector2& a, const Vector2& b, const Vector2& c, Real tol = NEAR_ZERO);
    static LongReal volume(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d);
    static VisibilityType volume_sign(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, Real tol = NEAR_ZERO);
    static bool coplanar(const Triangle& t1, const Triangle& t2, Real tol = NEAR_ZERO);
    static Vector3 to_3D(const Vector2& v, const Matrix3& RT);
    static Vector3 to_3D(const Vector2& v, const Matrix3& RT, Real offset);
    static Vector2 to_2D(const Vector3& v, const Matrix3& R);
    static bool point_in_tri(const Vector3& p, const Triangle& t, Real tol = NEAR_ZERO);
    static bool coplanar(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, Real tol = NEAR_ZERO);
    static Real calc_dist(const LineSeg3& line_seg, const Vector3& point, Real& t);
    static Real calc_closest_points(const LineSeg3& s1, const LineSeg3& s2, Vector3& p1, Vector3& p2); 
    static bool query_intersect_tri_tri(const Triangle& t1, const Triangle& t2);
    static PolygonLocationType in_tri(const Triangle& t, const Vector3& p, Real tol = NEAR_ZERO);

    /// Returns a if b > 0, -a if b < 0, and 0 if b = 0
    static Real sgn(Real a, Real b, Real tol = NEAR_ZERO) { if (b > NEAR_ZERO) return a; else if (b < -NEAR_ZERO) return -a; else return 0.0; }

    /// Determines whether two numbers are relatively/absolutely equal to within some tolerance
    /**
     * If the numbers are both smaller than zero, an absolute tolerance is used;
     * otherwise, a relative tolerance is used.
     * \note rationale for why this test is useful: [Ericson, 2005] p. 443
     */ 
    static bool rel_equal(Real x, Real y, Real tol = NEAR_ZERO) { return std::fabs(x-y) <= tol*std::max(std::fabs(x), std::max(std::fabs(y), (Real) 1.0)); }

    template <class T>
    static Real determine_line_param(const T& p, const T& dir, const T& v);

    template <class OutputIterator>
    static OutputIterator intersect_seg_tri(const LineSeg2& seg, const Vector2 tri[3], OutputIterator output_begin);

    template <class OutputIterator>
    static OutputIterator intersect_segs(const LineSeg2& s1, const LineSeg2& s2, OutputIterator output_begin);

     template <class OutputIterator>
     static OutputIterator intersect_tris(const Vector2 t1[3], const Vector2 t2[3], OutputIterator output_begin);

    template <class InputIterator>
    static void calc_min_area_bounding_rect(InputIterator begin, InputIterator end, Vector2& p1, Vector2& p2, Vector2& p3, Vector2& p4);

    template <class InputIterator, class OutputIterator>
    static OutputIterator intersect_polygons(InputIterator pbegin, InputIterator pend, InputIterator qbegin, InputIterator qend, const Vector3& normal, OutputIterator isect_begin);

    template <class InputIterator, class OutputIterator>
    static OutputIterator intersect_polygons(InputIterator pbegin, InputIterator pend, InputIterator qbegin, InputIterator qend, OutputIterator isect_begin);

    template <class InputIterator>
    static bool ccw(InputIterator begin, InputIterator end, const Vector3& normal, Real tol = NEAR_ZERO);

    template <class InputIterator>
    static bool ccw(InputIterator begin, InputIterator end, Real tol = NEAR_ZERO);

    template <class InputIterator>
    static bool is_convex_polygon(InputIterator begin, InputIterator end, Real tol = NEAR_ZERO);

    template <class InputIterator>
    static bool is_convex_polygon(InputIterator begin, InputIterator end, const Vector3& normal, Real tol = NEAR_ZERO);

    template <class InputIterator, class OutputIterator>
    static OutputIterator calc_convex_hull_2D(InputIterator source_begin, InputIterator source_end, const Vector3& normal, OutputIterator target_begin);

    template <class InputIterator, class OutputIterator>
    static OutputIterator calc_convex_hull_2D(InputIterator source_begin, InputIterator source_end, OutputIterator target_begin);

    template <class InputIterator>
    static PolygonLocationType polygon_location(InputIterator begin, InputIterator end, const Vector2& point);

    template <class ForwardIterator>
    static Vector2 calc_centroid_2D(ForwardIterator begin, ForwardIterator end);

    template <class InputIterator>
    static Vector3 calc_centroid_2D(InputIterator begin, InputIterator end, const Vector3& normal);

    template <class InputIterator>
    static Real calc_polygon_area(InputIterator begin, InputIterator end, const Vector3& normal);  

    template <class InputIterator>
    static Real calc_polygon_area(InputIterator begin, InputIterator end);

    template <class InputIterator>
    static void project_plane(InputIterator begin, InputIterator end, const Vector3& normal, Real offset);

    template <class InputIterator>
    static Real fit_plane(InputIterator begin, InputIterator end, Vector3& normal, Real& offset);

    template <class InputIterator, class OutputIterator>
    static OutputIterator triangulate_convex_polygon(InputIterator source_begin, InputIterator source_end, OutputIterator target_begin);

    template <class OutputIterator>
    static OutputIterator intersect_tris(const Triangle& t1, const Triangle& t2, OutputIterator begin);

    template <class OutputIterator>
    static OutputIterator intersect_coplanar_tris(const Triangle& t1, const Triangle& t2, const Vector3& normal, OutputIterator begin);

    template <class InputIterator, class OutputIterator>
    static OutputIterator to_2D(InputIterator begin, InputIterator end, OutputIterator begin_target, const Matrix3& R);

    template <class InputIterator, class OutputIterator>
    static OutputIterator to_3D(InputIterator begin, InputIterator end, OutputIterator begin_target, const Matrix3& RT, Real offset);

    template <class InputIterator, class T>
    static void determine_seg_endpoints(InputIterator begin, InputIterator end, std::pair<T, T>& endpoints);

    template <class InputIterator>
    static PolyhedronPtr calc_convex_hull_3D(InputIterator first, InputIterator last);

    template <class InputIterator>
    static Vector3 calc_centroid_3D(InputIterator first, InputIterator last);
    
    template <class InputIterator>
    static PolyhedronPtr calc_hs_intersection(InputIterator start, InputIterator end, const VectorN& interior_point);

    template <class InputIterator>
    static Real find_hs_interior_point(InputIterator start, InputIterator end, Vector3& point);

    template <class InputIterator>
    static unsigned calc_dimensionality(InputIterator first, InputIterator last, Real tol = NEAR_ZERO);

    template <class InputIterator>
    static bool intersect_seg_convex_polygon(InputIterator begin, InputIterator end, const LineSeg2& seg, Real& te, Real& tl, Real tol = NEAR_ZERO);

    template <class ForwardIterator, class OutputIterator>
    static OutputIterator intersect_seg_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, OutputIterator outbegin);

    template <class BidirectionalIterator, class OutputIterator>
    static OutputIterator triangulate_polygon_2D(BidirectionalIterator begin, BidirectionalIterator end, OutputIterator begintarget, Real tol = NEAR_ZERO);

    /// Returns the point a*t + (1-t)*b on the line segment ab
    template <class T>
    static T point_on_seg(T a, T b, Real t) { return b + (a - b)*t; }

  private:
    static void update_box(const Vector2& lP, const Vector2& rP, const Vector2& bP, const Vector2& tP, const Vector2& U, const Vector2& V, Real& min_area_div4, Vector2& center, Vector2 axis[2], Real extent[2]);
    static void project_onto_axis(const Triangle& t, const Vector3& v, Real& rfmin, Real& rfmax);
    static PolygonLocationType in_tri(OrientationType ori1, OrientationType ori2, OrientationType ori3);
    static bool intersect(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d, Real tol = NEAR_ZERO);
    static bool intersect_prop(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d, Real tol = NEAR_ZERO);
    static void clip_convex_polygon_against_line(const Vector2& rkN, Real fC, unsigned& riQuantity, Vector2 akV[6]);

    template <class BidirectionalIterator>
    static bool in_cone(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, Real tol = NEAR_ZERO);

    template <class BidirectionalIterator>
    static bool diagonal(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, Real tol = NEAR_ZERO);

    template <class BidirectionalIterator>
    static bool diagonalie(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, Real tol = NEAR_ZERO);

    template <class OutputIterator>
    static unsigned advance(unsigned a, unsigned* aa, unsigned n, bool inside, const Vector2& v, OutputIterator& current);

    static bool same_side(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& b, Real tol = std::sqrt(std::numeric_limits<Real>::epsilon()));
    static pthread_mutex_t _qhull_mutex;
    static bool compute_intervals_isectline(const Triangle& t1, Real vv0, Real vv1, Real vv2, Real d0, Real d1, Real d2, Real d0d1, Real d0d2, Real& isect0, Real& isect1, Vector3& isectpoint0, Vector3& isectpoint1);
    static void isect2X(const Vector3& vtx0, const Vector3& vtx1, const Vector3& vtx2, Real vv0, Real vv1, Real vv2, Real d0, Real d1, Real d2, Real& isect0, Real& isect1, Vector3& isectpoint0, Vector3& isectpoint1);
    static SegSegIntersectType get_parallel_intersect_type(const LineSeg2& s1, const LineSeg2& s2, Vector2& isect, Vector2& isect2);
    static bool between(const Vector2& a, const Vector2& b, const Vector2& c, Real tol = NEAR_ZERO);
    static SegTriIntersectType seg_tri_cross(const LineSeg3& seg, const Triangle& t, Real tol = NEAR_ZERO);
    static SegTriIntersectType intersect_seg_tri_in_plane(const LineSeg3& seg, const Triangle& t, Vector3& isect, Vector3& isect2, Real tol = NEAR_ZERO);
    static PolygonLocationType in_tri(const Triangle& t, const Vector3& p, unsigned int skip_dim, Real tol = NEAR_ZERO);
    static void find_plane_coeff(const Triangle& t, Vector3& N, unsigned int& m, Real& dot);
    static Vector3 normal_vec(const Vector3& a, const Vector3& b, const Vector3& c);
    static OrientationType area_sign(const Vector2& a, const Vector2& b, const Vector2& c, Real tol = NEAR_ZERO);
    static LongReal area(const Vector2& a, const Vector2& b, const Vector2& c);
    static bool test_edge_edge(const Vector3& p, const Vector3& a, const Vector3& b, Real Ax, Real Ay, unsigned i0, unsigned i1);
    static bool test_edge_tri(const Vector3& a, const Vector3& b, const Triangle& t, unsigned i0, unsigned i1);
    static bool test_point_in_tri(const Vector3& p, const Triangle& t, unsigned i0, unsigned i1);
    static bool test_coplanar_tri_tri(const Vector3& N, const Triangle& t1, const Triangle& t2);
}; // end class

// include inline functions
#include "CompGeom.inl"

} // end namespace

#endif
