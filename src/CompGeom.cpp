/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <iterator>
#include <cmath>
#include <set>
#include <stack>
#include <fstream>
#include <Ravelin/Origin2d.h>
#include <Moby/Polyhedron.h>
#include <Moby/CompGeom.h>

/// Needed for qhull
#ifdef THREADSAFE
pthread_mutex_t Moby::CompGeom::_qhull_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

using namespace Ravelin;
using namespace Moby;

/// Helper function for calc_min_area_rect()
void CompGeom::update_box(const Point2d& lP, const Point2d& rP, const Point2d& bP, const Point2d& tP, const Vector2d& U, const Vector2d& V, double& min_area_div4, Point2d& center, Vector2d axis[2], double extent[2])
{
  Vector2d rlDiff = rP - lP;
  Vector2d tbDiff = tP - bP;
  double next[2] = { (double) 0.5 * (U.dot(rlDiff)), (double) 0.5 * (V.dot(tbDiff)) };
  double area_div4 = next[0] * next[1];
  if (area_div4 < min_area_div4)
  {
    min_area_div4 = area_div4;
    axis[0] = U;
    axis[1] = V;
    extent[0] = next[0];
    extent[1] = next[1];
    Vector2d lbDiff = lP - bP;
    center = lP + U*extent[0] + V*(extent[1] - V.dot(lbDiff));
  }
}

/// Gets the number of intersections for a seg/seg intersection code
unsigned CompGeom::get_num_intersects(CompGeom::SegSegIntersectType t)
{
  switch (t)
  {
    case eSegSegNoIntersect: return 0;
    case eSegSegIntersect: return 1;
    case eSegSegVertex: return 1;
    case eSegSegEdge: return 2;
    default:
      assert(false);
  }

  // make compiler happy
  return 0;
}

/// Gets the number of intersections for a line/line intersection code
unsigned CompGeom::get_num_intersects(CompGeom::LineLineIntersectType t)
{
  switch (t)
  {
    case eLineLineNoIntersect: return 0;
    case eLineLineIntersect: return 1;
    case eLineLineVertex: return 1;
    case eLineLineEdge: return 2;
    default:
      assert(false);
  }

  // make compiler happy
  return 0;
}

/// Gets the number of intersections for a seg/tri intersection code
unsigned CompGeom::get_num_intersects(CompGeom::SegTriIntersectType t)
{
  switch (t)
  {
    case eSegTriNoIntersect: return 0;
    case eSegTriInclFace: return 1;
    case eSegTriFace: return 1;
    case eSegTriInclVertex: return 1;
    case eSegTriVertex: return 1;
    case eSegTriInclEdge: return 1;
    case eSegTriEdge: return 1;
    case eSegTriInside: return 2;
    case eSegTriEdgeOverlap: return 2;
    case eSegTriPlanarIntersect: return 2;
    default:
      assert(false);
  }

  // make compiler happy
  return 0;
}


/// Utility method for triangulate_polygon_2D()
bool CompGeom::intersect_prop(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d, double tol)
{
  // eliminate improper cases
  if (collinear(a,b,c,tol) || collinear(a,b,d,tol) || 
      collinear(c,d,a,tol) || collinear(c,d,b,tol))
    return false;

  return (!area_sign(a,b,c,tol) == eLeft ^ !area_sign(a,b,d,tol) == eLeft) && (!area_sign(c,d,a,tol) == eLeft ^ !area_sign(c,d,b,tol) == eLeft);
}

/// Utility method for triangulate_polygon_2D()
bool CompGeom::intersect(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d, double tol)
{
  if (intersect_prop(a,b,c,d,tol))
    return true;
  else if (between(a,b,c,tol) || between(a,b,d,tol) || between(c,d,a,tol) || between(c,d,b,tol))
    return true;
  else
    return false;
}

/// Computes the shortest line segment between two line segments in 3D
/**
 * \param s1 the first segment
 * \param s2 the second segment
 * \param p1 the closest point on s1 to s2 (on return)
 * \param p2 the closest point on s2 to s1 (on return)
 * \return the distance between p1 and p2
 * \note algorithm taken from http://geometricalgorithms.com
 */
double CompGeom::calc_closest_points(const LineSeg3& s1, const LineSeg3& s2, Point3d& p1, Point3d& p2)
{
  double sc, sN, tc, tN;

  Vector3d u = s1.second - s1.first;
  Vector3d v = s2.second - s2.first;
  Vector3d w = s1.first - s2.first;
  double a = Vector3d::dot(u,u);
  double b = Vector3d::dot(u,v);
  double c = Vector3d::dot(v,v);
  double d = Vector3d::dot(u,w);
  double e = Vector3d::dot(v,w);
  double D = a*c - b*b;
  double sD = D, tD = D;

  // compute the line parameters of the two closest points
  if (D < NEAR_ZERO)
  {
    // the lines are almost parallel; force using endpoint on s1 to
    // prevent possible division by zero later
    sN = 0.0;
    sD = 1.0;
    tN = e;
    tD = c;
  }
  else
  {
    // get the closest points on the infinite lines
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0)
    {
      sN = 0.0;
      tN = e;
      tD = c;
    }
    else if (sN > sD)
    {
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < 0.0)
  {
    tN = 0.0;
    if (-d < 0.0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else
    {
      sN = -d;
      sD = a;
    }
  }
  else if (tN > tD)
  {
    tN = tD;
    if ((-d + b) < 0.0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else
    {
      sN = (-d + b);
      sD = a;
    }
  }

  // do the division to get sc and tc
  sc = (std::fabs(sN) < NEAR_ZERO) ? (double) 0.0 : sN / sD;
  tc = (std::fabs(tN) < NEAR_ZERO) ? (double) 0.0 : tN / tD;

  // compute p1 and p2
  p1 = s1.first + (s1.second - s1.first)*sc;
  p2 = s2.first + (s2.second - s2.first)*tc;

  // compute the distance
  return (p1 - p2).norm();
}

/// Calculates the Euclidean distance between a line segment and a point in 3D
/**
 * \param line_seg the two endpoints of the line segment
 * \param point the query point
 * \param t gives the parametric form line_seg.first + (line_seg.second - line_seg.first) * t 
 * \param closest point on the line segment to point 
 */
double CompGeom::calc_dist(const LineSeg3& line_seg, const Point3d& point, double& t, Point3d& closest)
{
  Vector3d v = line_seg.second - line_seg.first;
  Vector3d w = point - line_seg.first;
  double c1 = Vector3d::dot(w, v);
  if (c1 <= 0)
  {
    closest = line_seg.first;
    return w.norm();
  }
  double c2 = Vector3d::dot(v, v);
  if (c2 <= c1)
  {
    closest = line_seg.second;
    return (point - line_seg.second).norm();
  }
  t = c1 / c2;
  closest = line_seg.first + (line_seg.second - line_seg.first)*t;
  return (closest - point).norm();
}

/// Determines whether four points are coplanar
bool CompGeom::coplanar(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d, double tol)
{
  return volume_sign(a, b, c, d, tol) == eCoplanar;
}

/// Generates a point on a plane
Point3d CompGeom::generate_point_on_plane(const Vector3d& normal, double d)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the absolute values of the normal
  double abs_x = std::fabs(normal[X]);
  double abs_y = std::fabs(normal[Y]);
  double abs_z = std::fabs(normal[Z]);

  // solve using the largest component of the normal
  if (abs_x > abs_y)
  {
    if (abs_x > abs_z)
      return Point3d(d/normal[X], 0, 0);
    else
      return Point3d(0, 0, d/normal[Z]);
  }
  else
  {
    if (abs_y > abs_z)
      return Point3d(0, d/normal[Y], 0);
    else
      return Point3d(0, 0, d/normal[Z]);
  } 
}

/// Determines whether two triangles are coplanar
bool CompGeom::coplanar(const Triangle& t1, const Triangle& t2, double tol)
{
  if (volume_sign(t1.a, t1.b, t1.c, t2.a, tol) != eCoplanar)
    return false;
  if (volume_sign(t1.a, t1.b, t1.c, t2.b, tol) != eCoplanar)
    return false;
  if (volume_sign(t1.a, t1.b, t1.c, t2.c, tol) != eCoplanar)
    return false;
  return true;
}

/// Determines whether two line segments in 3D intersect
/**
 * \param s1 the points of the first line segment
 * \param s2 the points of the second line segment
 * \param isect the point of intersection (if any)
 * \param isect2 the second point of intersection (for edge intersection)
 * \return the intersection type (edge, vertex, true_intersect, no_intersect); 
 *         an edge intersection indicates that the segments collinearly overlap;
 *         vertex intersection indicates that an endpoint of
 *         one segment is on the other segment, but edge intersection does not
 *         hold; true_intersect indicates that the segments intersect at a
 *         single point.
 */
CompGeom::SegSegIntersectType CompGeom::intersect_segs(const LineSeg3& s1, const LineSeg3& s2, Point3d& isect, Point3d& isect2)
{
  const int X = 0, Y = 1, Z = 2;
  
  // get the individual points 
  const Point3d& a = s1.first;
  const Point3d& b = s1.second;
  const Point3d& c = s2.first;
  const Point3d& d = s2.second;

  // compute the two vectors
  Vector3d v1 = b - a;
  Vector3d v2 = d - c;
  
  // compute the squared norms for the two vectors
  double v1_sq_norm = Vector3d::dot(v1, v1);
  double v2_sq_norm = Vector3d::dot(v2, v2);

  // compute the dot-product of (v1, v2)
  double v1_v2 = Vector3d::dot(v1, v2);
  
  // if the segments are parallel, handle them specially
  if (std::fabs(v1_sq_norm * v2_sq_norm - v1_v2 * v1_v2) < NEAR_ZERO)
  {
    // if the segments are not collinear, can return now
    if (!collinear(a, b, c))
      return eSegSegNoIntersect;
    
    // get the dimension with the most variance
    unsigned bd;
    double v1x = v1[X] * v1[X];
    double v1y = v1[Y] * v1[Y];
    double v1z = v1[Z] * v1[Z];
    if (v1x > v1y)
      bd = (v1x > v1z) ? X : Z;
    else
      bd = (v1y > v1z) ? Y : Z;
    
    // get minimum and maximum points for v1
    const Point3d* v1_min, * v1_max;
    if (a[bd] > b[bd]) 
    {
      v1_min = &b;
      v1_max = &a;
    }
    else
    {
      v1_min = &a;
      v1_max = &b;
    }
    
    // get minimum and maximum points for v2
    const Point3d* v2_min, * v2_max;
    if (c[bd] > d[bd]) 
    {
      v2_min = &d;
      v2_max = &c;
    }
    else
    {
      v2_min = &c;
      v2_max = &d;
    }

    // see whether the two do not overlap
    if ((*v1_max)[bd] < (*v2_min)[bd] || (*v2_max)[bd] < (*v1_min)[bd])
      return eSegSegNoIntersect;

    // see whether the two meet at only a vertex
    if (std::fabs((*v1_max)[bd] - (*v2_min)[bd]) < NEAR_ZERO)
    {
      isect = *v1_max;

      return eSegSegVertex;
    }
    if (std::fabs((*v2_max)[bd] - (*v1_min)[bd]) < NEAR_ZERO)
    {
      isect = *v2_max;
      return eSegSegVertex;
    }    
  
    // the two overlap
    if ((*v1_max)[bd] > (*v2_max)[bd] && (*v1_min)[bd] < (*v2_max)[bd])
    {
      isect = ((*v1_min)[bd] > (*v2_min)[bd]) ? *v1_min : *v2_min;
      isect2 = *v2_max;
    }
    else
    {
      isect = ((*v2_min)[bd] > (*v1_min)[bd]) ? *v2_min : *v1_min;
      isect2 = *v1_max;
    }
    return eSegSegEdge;
  }

  // determine the point of intersection for the two lines
  Vector3d a_minus_c = a - c;
  Vector3d g1 = Vector3d::cross(a_minus_c, v2);
  Vector3d g2 = Vector3d::cross(v1, v2);
  double p = std::sqrt(Vector3d::dot(g1, g1) / Vector3d::dot(g2, g2));
  if (Vector3d::dot(g1, g2) > 0)
    p = -p;  
  isect = a + p * v1;
  
  // verify that the point is within the first segment
  Vector3d isect_m_b = isect - b;
  double t1 = std::sqrt(Vector3d::dot(isect_m_b, isect_m_b) / Vector3d::dot(-v1, -v1));
  if (Vector3d::dot(isect_m_b, -v1) < 0)
    t1 = -t1;
  if (Vector3d::dot(isect_m_b, -v1) < 0 || t1 > 1)
    return eSegSegNoIntersect;
  
  // verify that the point is within the second segment
  Vector3d isect_m_d = isect - d;
  double t2 = std::sqrt(Vector3d::dot(isect_m_d, isect_m_d) / Vector3d::dot(-v2, -v2));
  if (Vector3d::dot(isect_m_d, -v2) < 0)
    t2 = -t2;
  if (Vector3d::dot(isect_m_d, -v2) < 0 || t2 > 1)
    return eSegSegNoIntersect;

  // true intersection, return the type
  if (t1 < NEAR_ZERO || std::fabs(t1 - 1) < NEAR_ZERO)
    return eSegSegVertex;
  if (t2 < NEAR_ZERO || std::fabs(t2 - 1) < NEAR_ZERO)
    return eSegSegVertex;
  
  return eSegSegIntersect;  
}

/// Utility method for determining triangle / segment intersection
/**
 * This method determines the type of intersection for segment that crosses
 * triangles plane.
 * \note Adapted from O'Rourke, p. 237.
 * \throw NumericalException if segment lies in triangle's plane
 */
CompGeom::SegTriIntersectType CompGeom::seg_tri_cross(const LineSeg3& seg, const Triangle& t, double tol)
{
  VisibilityType v0 = volume_sign(seg.first, t.a, t.b, seg.second, tol);
  VisibilityType v1 = volume_sign(seg.first, t.b, t.c, seg.second, tol);
  VisibilityType v2 = volume_sign(seg.first, t.c, t.a, seg.second, tol);

  // look for the case of intersection with the face of the triangle
  if (((v0 == eInvisible) && (v1 == eInvisible) && (v2 == eInvisible)) ||
      ((v0 == eVisible) && (v1 == eVisible) && (v2 == eVisible)))
    return eSegTriInclFace;

  // look for the case of no intersection
  if (((v0 == eInvisible) || (v1 == eInvisible) || (v2 == eInvisible)) &&
      ((v0 == eVisible) || (v1 == eVisible) || (v2 == eVisible)))
    return eSegTriNoIntersect;
  else if ((v0 == eCoplanar) && (v1 == eCoplanar) && (v2 == eCoplanar))
    throw NumericalException("Segment and triangle unexpectedly coplanar"); 
  else if ((v0 == eCoplanar && v1 == eCoplanar) ||
           (v0 == eCoplanar && v2 == eCoplanar) ||
           (v1 == eCoplanar && v2 == eCoplanar))
    return eSegTriInclVertex;
  else if (v0 == eCoplanar || v1 == eCoplanar || v2 == eCoplanar)
    return eSegTriInclEdge;

  // should not make it here..
  throw std::runtime_error("Unexpected control path");
  return eSegTriNoIntersect;
}

/// Determines the intersection between a line segment and a plane in 3D
/**
 * \param nrm the plane normal 
 * \param d the plane offset such that nrm*x - d = 0
 * \param seg the line segment 
 * \return the parameter of the line, t, such that the point of intersection is
 *         determined by seg.first + (seg.second - seg.first)*t, or NaN if the 
 *         line and plane are parallel
 */
double CompGeom::intersect_seg_plane(const Vector3d& nrm, double d, const LineSeg3& seg)
{
  // generate a point on the plane
  Point3d V0 = generate_point_on_plane(nrm, d);

  // compute t
  return Vector3d::dot(nrm, V0 - seg.first)/Vector3d::dot(nrm, seg.second - seg.first);
}

/// Determines the intersection between a line segment and a plane in 3D
/**
 * \param tri a 3D triangle that is fully contained in the plane to be 
 *            intersected
 * \param seg a 3D line segment
 * \param isect the point of intersection (on return)
 * \param big_dim the dimension of the largest component (on return)
 * \param NEAR_ZERO the tolerance to test for zero
 * \return the intersection type: eSegPlaneInPlane (segment lies wholly within the
 *         plane), eSegPlaneFirst (first vertex of segment is on the plane, but not
 *         the second), eSegPlaneSecond (second vertex of segment is on the plane, but
 *         not the first), eSegPlaneToSide (the segment lies strictly to one side
 *         or the other of the plane), or eSegPlaneOtherIntersect (segment intersects
 *         the plane, and in-plane, first, and second do not hold).
 */
CompGeom::SegPlaneIntersectType CompGeom::intersect_seg_plane(const Triangle& tri, const LineSeg3& seg, Point3d& isect, unsigned int& big_dim, double tol)
{
  // find the plane coefficients
  double dot;
  Vector3d N;
  find_plane_coeff(tri, N, big_dim, dot);

  // check to see whether we are only testing a point; if so, see whether the
  // point is on the plane
  Vector3d rq = seg.second - seg.first;
  double norm_rq = rq.norm();
  if (norm_rq == (double) 0.0)
  {
    if (std::fabs(Vector3d::dot(seg.first, N) - dot) < NEAR_ZERO)
      return eSegPlaneInPlane;
    else
      return eSegPlaneToSide;
  }    

  // check whether segment is parallel to the plane
  long double num = dot - Vector3d::dot(seg.first, N);
  long double denom = Vector3d::dot(rq, N);
  if (denom < tol && denom > -tol)
  {
    // check whether both points are on plane
    if (num < tol && num > -tol)
      return eSegPlaneInPlane;
    else
      return eSegPlaneToSide;
  }
  long double t = num / denom;

  // compute point of intersection
  isect = seg.first + (t * (seg.second - seg.first));

  if ((t > tol) && (t < (double) 1.0 - tol))
    return eSegPlaneOtherIntersect;
  else if (num < tol && num > -tol)
    return eSegPlaneFirst;
  else if (std::fabs(num-denom) < tol)
    return eSegPlaneSecond;

  return eSegPlaneToSide;
}

/// Determines the feature for a line segment given a parameter
CompGeom::SegLocationType CompGeom::determine_seg_location(const LineSeg3& seg, double t)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Point3d& p = seg.first;
  const Point3d& q = seg.second;

  // get the largest value for the first point
  double lg_p = std::max(std::max(std::fabs(p[X]), std::fabs(p[Y])), std::fabs(p[Z]));

  // determine the tolerance for the first point
  double p_tol = (lg_p > 0.0) ? NEAR_ZERO : NEAR_ZERO/lg_p;

  // see whether on the origin 
  if (std::fabs(t) < p_tol)
    return eSegOrigin;

  // get the largest value for the second point
  double lg_q = std::max(std::max(std::fabs(q[X]), std::fabs(q[Y])), std::fabs(q[Z]));

  // determine the tolerance for the second point
  double q_tol = (lg_q > 0.0) ? NEAR_ZERO : NEAR_ZERO/lg_q;

  // see whether on the endpoint
  if (std::fabs(t - 1.0) < q_tol)
    return eSegEndpoint;

  // on neither endpoint; see whether it is in the interior of the segment 
  return (t > 0.0 && t < 1.0) ? eSegInterior : eSegOutside;
}

/// Determines the feature for a line segment given a parameter
CompGeom::SegLocationType CompGeom::determine_seg_location(const LineSeg2& seg, double t)
{
  const unsigned X = 0, Y = 1;
  const Point2d& p = seg.first;
  const Point2d& q = seg.second;

  // get the largest value for the first point
  double lg_p = std::max(std::fabs(p[X]), std::fabs(p[Y]));

  // determine the tolerance for the first point
  double p_tol = (lg_p > 0.0) ? NEAR_ZERO : NEAR_ZERO/lg_p;

  // see whether on the origin 
  if (std::fabs(t) < p_tol)
    return eSegOrigin;

  // get the largest value for the second point
  double lg_q = std::max(std::fabs(q[X]), std::fabs(q[Y]));

  // determine the tolerance for the second point
  double q_tol = (lg_q > 0.0) ? NEAR_ZERO : NEAR_ZERO/lg_q;

  // see whether on the endpoint
  if (std::fabs(t - 1.0) < q_tol)
    return eSegEndpoint;

  // on neither endpoint; see whether it is in the interior of the segment 
  return (t > 0.0 && t < 1.0) ? eSegInterior : eSegOutside;
}

/// Determines the intersection between a line segment and triangle in 2D
/**
 * \param seg a line segment in 2D
 * \param tri a triangle in 2D
 * \param isect contains the point of intersection, if any (on return)
 * \param isect2 contains a second point of intersection, if the intersection
 *               type is eSegTriInside, eSegTriEdgeOverlap, or 
 *               eSegTriPlanarIntersect 
 *               (on return)
 * \return eSegTriInside (if the segment lies wholly within the triangle face
 *         eSegTriVertex (if an endpoint of the segment coincides with a vertex
 *         of the triangle), eSegTriEdge (if an endpoint of the segment is in
 *         the relative interior of an edge of the triangle),  
 *         eSegTriEdgeOverlap (if the segment overlaps one edge of the
 *         triangle [or vice versa]), eSegTriPlanarIntersect (if the
 *         segment lies partially within the face of the triangle), or
 *         eSegTriNoIntersect (if there is no intersection).
 * \note assumes triangle is oriented ccw 
 */
CompGeom::SegTriIntersectType CompGeom::intersect_seg_tri(const LineSeg2& seg, const Point2d tri[3], Point2d& isect, Point2d& isect2, double tol)
{
  // get the two points
  const Point2d& p = seg.first;
  const Point2d& q = seg.second;

  // make segments out of the edges
  LineSeg2 e1(tri[0], tri[1]);
  LineSeg2 e2(tri[1], tri[2]);
  LineSeg2 e3(tri[2], tri[0]);

  // get the orientations w.r.t. the edges
  OrientationType pe1 = area_sign(tri[0], tri[1], p);
  OrientationType pe2 = area_sign(tri[1], tri[2], p);
  OrientationType pe3 = area_sign(tri[2], tri[0], p);
  OrientationType qe1 = area_sign(tri[0], tri[1], q);
  OrientationType qe2 = area_sign(tri[1], tri[2], q);
  OrientationType qe3 = area_sign(tri[2], tri[0], q);

  // see whether points are outside (most common)
  if ((pe1 == eRight && qe1 == eRight) ||
      (pe2 == eRight && qe2 == eRight) ||
      (pe3 == eRight && qe3 == eRight))
    return eSegTriNoIntersect;

  // check for edge overlap
  if (pe1 == eOn && qe1 == eOn)
  {
    // do parallel segment / segment intersection
    Vector2d dir = e1.second - e1.first;
    double t1 = determine_line_param(e1.first, dir, p);
    double t2 = determine_line_param(e1.first, dir, q);
    if (t2 < t1)
      std::swap(t1, t2); 
    t1 = std::max(t1, (double) 0.0);
    t2 = std::min(t2, (double) 1.0);
    if (rel_equal(t1, t2))
    {
      // vertex only
      isect = e1.first + (e1.second - e1.first)*t1;
      return eSegTriVertex;
    }
    else if (t2 > t1)
    {
      // edge overlap
      isect = e1.first + dir*t1;
      isect2 = e1.first + dir*t2;
      return eSegTriEdgeOverlap;
    }
    else
      // no overlap
      return eSegTriNoIntersect;  
  }
  if (pe2 == eOn && qe2 == eOn)
  {
    // do parallel segment / segment intersection
    Vector2d dir = e2.second - e2.first;
    double t1 = determine_line_param(e2.first, dir, p);
    double t2 = determine_line_param(e2.first, dir, q);
    if (t2 < t1)
      std::swap(t1, t2); 
    t1 = std::max(t1, (double) 0.0);
    t2 = std::min(t2, (double) 1.0);
    if (rel_equal(t1, t2))
    {
      // vertex only
      isect = e2.first + (e2.second - e2.first)*t1;
      return eSegTriVertex;
    }
    else if (t2 > t1)
    {
      // edge overlap
      isect = e2.first + dir*t1;
      isect2 = e2.first + dir*t2;
      return eSegTriEdgeOverlap;
    }
    else
      // no overlap
      return eSegTriNoIntersect;  
  }
  if (pe3 == eOn && qe3 == eOn)
  {
    Vector2d dir = e3.second - e3.first; 
    double t1 = determine_line_param(e3.first, dir, p);
    double t2 = determine_line_param(e3.first, dir, q);
    if (t2 < t1)
      std::swap(t1, t2); 
    t1 = std::max(t1, (double) 0.0);
    t2 = std::min(t2, (double) 1.0);
    if (rel_equal(t1, t2))
    {
      // vertex only
      isect = e3.first + (e3.second - e3.first)*t1;
      return eSegTriVertex;
    }
    else if (t2 > t1)
    {
      // edge overlap
      isect = e3.first + dir*t1;
      isect2 = e3.first + dir*t2;
      return eSegTriEdgeOverlap;
    }
    else
      // no overlap
      return eSegTriNoIntersect;  
  }

   // check for edge or vertex intersection with edge #1
  if (pe1 == eRight)
  {
    // if q is on edge 1, it may be on a vertex, an edge, or neither
    if (qe1 == eOn)
    {
      Vector2d dir = e1.second - e1.first;
      double t = determine_line_param(e1.first, dir, q);
      SegLocationType feat = determine_seg_location(e1, t);
      isect = e1.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // q *should* be to the left of e1; see whether q is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType qloc = in_tri(tri, q);
      if (qloc == ePolygonOnVertex)
      {
        isect2 = q;
        return eSegTriVertex;
      }
      else if (qloc == ePolygonOnEdge)
      {
        isect2 = q;
        return eSegTriEdge;
      }
      else if (qloc == ePolygonInside)
      {
        // intersect seg vs. e1 
        SegSegIntersectType t = intersect_segs(seg, e1, isect, isect2);
        isect2 = q;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            break;

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }
    }
  }
  else if (qe1 == eRight)
  {
    if (pe1 == eOn)
    {
      Vector2d dir = e1.second - e1.first;
      double t = determine_line_param(e1.first, dir, p);
      SegLocationType feat = determine_seg_location(e1, t);
      isect = e1.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // p *should* be to the left of e1; see whether p is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType ploc = in_tri(tri, p);
      if (ploc == ePolygonOnVertex)
      {
        isect2 = p;
        return eSegTriVertex;
      }
      else if (ploc == ePolygonOnEdge)
      {
        isect2 = p;
        return eSegTriEdge;
      }
      else if (ploc == ePolygonInside)
      {
        // intersect seg vs. e1
        SegSegIntersectType t = intersect_segs(seg, e1, isect, isect2);
        isect2 = p;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            break; 

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }  
    }
  }

  // check for edge or vertex intersection with edge #2
  if (pe2 == eRight)
  {
    if (qe2 == eOn)
    {
      Vector2d dir = e2.second - e2.first;
      double t = determine_line_param(e2.first, dir, q);
      SegLocationType feat = determine_seg_location(e2, t);
      isect = e2.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // q *should* be to the left of e2; see whether q is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType qloc = in_tri(tri, q);
      if (qloc == ePolygonOnVertex)
      {
        isect2 = q;
        return eSegTriVertex;
      }
      else if (qloc == ePolygonOnEdge)
      {
        isect2 = q;
        return eSegTriEdge;
      }
      else if (qloc == ePolygonInside)
      {
        // intersect seg vs. e2
        SegSegIntersectType t = intersect_segs(seg, e2, isect, isect2);
        isect2 = q;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            break; 

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }
    }
  }
  else if (qe2 == eRight)
  {
    if (pe2 == eOn)
    {
      Vector2d dir = e2.second - e2.first;
      double t = determine_line_param(e2.first, dir, p);
      SegLocationType feat = determine_seg_location(e2, t);
      isect = e2.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // p *should* be to the left of e2; see whether p is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType ploc = in_tri(tri, p);
      if (ploc == ePolygonOnVertex)
      {
        isect2 = p;
        return eSegTriVertex;
      }
      else if (ploc == ePolygonOnEdge)
      {
        isect2 = p;
        return eSegTriEdge;
      }
      else if (ploc == ePolygonInside)
      {
        // intersect seg vs. e2
        SegSegIntersectType t = intersect_segs(seg, e2, isect, isect2);
        isect2 = p;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            return eSegTriEdge;

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }      
    }
  }

  // check for edge or vertex intersection with edge #3
  if (pe3 == eRight)
  {
    if (qe3 == eOn)
    {
      Vector2d dir = e3.second - e3.first;
      double t = determine_line_param(e3.first, dir, q);
      SegLocationType feat = determine_seg_location(e3, t);
      isect = e3.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // q *should* be to the left of e3; see whether q is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType qloc = in_tri(tri, q);
      if (qloc == ePolygonOnVertex)
      {
        isect2 = q;
        return eSegTriVertex;
      }
      else if (qloc == ePolygonOnEdge)
      {
        isect2 = q;
        return eSegTriEdge;
      }
      else if (qloc == ePolygonInside)
      {
        // intersect seg vs. e3
        SegSegIntersectType t = intersect_segs(seg, e3, isect, isect2);
        isect2 = q;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            return eSegTriEdge;

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }
    }
  }
  else if (qe3 == eRight)
  {
    if (pe3 == eOn)
    {
      Vector2d dir = e3.second - e3.first;
      double t = determine_line_param(e3.first, dir, p);
      SegLocationType feat = determine_seg_location(e3, t);
      isect = e3.first + dir*t;
      if (feat == eSegOrigin || feat == eSegEndpoint)
        return eSegTriVertex;
      else if (feat == eSegInterior)
        return eSegTriEdge;
    }
    else
    {
      // p *should* be to the left of e3; see whether p is on a vertex, an
      // edge, or inside the triangle 
      PolygonLocationType ploc = in_tri(tri, p);
      if (ploc == ePolygonOnVertex)
      {
        isect2 = p;
        return eSegTriVertex;
      }
      else if (ploc == ePolygonOnEdge)
      {
        isect2 = p;
        return eSegTriEdge;
      }
      else if (ploc == ePolygonInside)
      {
        // intersect seg vs. e3
        SegSegIntersectType t = intersect_segs(seg, e3, isect, isect2);
        isect2 = p;
        std::swap(isect, isect2);
        switch (t)
        {
          case eSegSegIntersect:  
            return eSegTriPlanarIntersect;

          case eSegSegNoIntersect:
            return eSegTriEdge;

          // none of these should occur, so we'll return the points
          // of intersection and fudge the type 
          case eSegSegVertex:  
          case eSegSegEdge:
            return eSegTriPlanarIntersect;
        }
      }
    }
  }
 
  // if we're here, one of two cases has occurred: both points are inside
  // the triangle or the segment does not intersect the triangle; find out
  // which it is
  if (pe1 != eRight && pe2 != eRight && pe3 != eRight && 
      qe1 != eRight && qe2 != eRight && qe3 != eRight)
  { 
    isect = p;
    isect2 = q;
    return eSegTriInside;
  }
  else
  {
//    assert(in_tri(tri, p) == ePolygonOutside && in_tri(tri, q) == ePolygonOutside);
    return eSegTriNoIntersect;
  }

// NOTE: this code from www.geometrictools.com does not work!
/*
  double afDist[3];
  int aiSign[3], iPositive = 0, iNegative = 0, iZero = 0;

  // do the simple test first
  if (rel_equal(seg.first, seg.second, tol))
  {
    PolygonLocationType ptype = in_tri(tri, seg.first);
    switch (ptype)
    {
      case ePolygonOnVertex: 
        isect = seg.first;
        return eSegTriVertex;

      case ePolygonOnEdge:
        isect = seg.first;
        return eSegTriEdge;

      case ePolygonInside:
      case ePolygonOutside:
        return eSegTriNoIntersect;

      default:
        assert(false);
    }
    
    // should not get here
    assert(false);
    return eSegTriNoIntersect;
  }

  // convert segment to necessary representation
  Vector2d origin = (seg.first+seg.second)*0.5;
  Vector2d dir = Vector2::normalize(seg.second - seg.first);

  for (unsigned i=0; i< 3; i++)
  {
    Vector2d kDiff = tri[i] - origin;
    afDist[i] = kDiff.dot_perp(dir);
    if (afDist[i] > tol)
    {
      aiSign[i] = 1;
      iPositive++;
    }
    else if (afDist[i] < -tol)
    {
      aiSign[i] = -1;
      iNegative++;
    }
    else
    {
      afDist[i] = 0.0;
      aiSign[i] = 0;
      iZero++;
    }
  }

  // look for no intersection
  if (iPositive == 3 || iNegative == 3)
    return eSegTriNoIntersect;

  double afParam[2];
  
  // *** get interval
  // project triangle onto line
  double afProj[3];
  for (unsigned i=0; i< 3; i++)
  {
    Vector2d kDiff = tri[i] - origin;
    afProj[i] = dir.dot(kDiff);
  }

  // compute transverse intersections of triangle edges with line
  double fNumer, fDenom;
  unsigned iQuantity = 0;
  for (unsigned i0=2, i1 = 0; i1< 3; i0 = i1++)
  {
    if (aiSign[i0]*aiSign[i1] < 0)
    {
      assert(iQuantity < 2);
      fNumer = afDist[i0]*afProj[i1] - afDist[i1]*afProj[i0];
      fDenom = afDist[i0] - afDist[i1];
      afParam[iQuantity++] = fNumer/fDenom;
    }
  }

  // check for grazing contact
  if (iQuantity < 2)
    for (unsigned i0=1, i1=2, i2 = 0; i2 < 3; i0 = i1, i1 = i2, i2++)
      if (aiSign[i2] == 0)
      {
        assert(iQuantity < 2);
        afParam[iQuantity++] = afProj[i2];
      }

  // sort
  assert(iQuantity >= 1);
  if (iQuantity == 2)
  {
    if (afParam[0] > afParam[1])
      std::swap(afParam[0], afParam[1]);
  }
  else
    afParam[1] = afParam[0];

  // *** intersection
  const double afU0 = afParam[0];
  const double afU1 = afParam[1];
  const double afV0 = -0.5;         // minimum segment parameter
  const double afV1 = 0.5;          // maximum segment parameter
  double afOverlap0, afOverlap1;
  unsigned m_iQuantity = 0;
  if (afU1 < afV0 || afU0 > afV1)
    m_iQuantity = 0;
  else if (afU1 > afV0)
  {
    if (afU0 < afV1)
    {
      m_iQuantity = 2;
      afOverlap0 = (afU0 < afV0) ? afV0 : afU0;
      afOverlap1 = (afU1 > afV1) ? afV1 : afU1;
      if (afOverlap0 == afOverlap1)
        m_iQuantity = 1;
    }
    else // afU0 == afV1
    {
      m_iQuantity = 1;
      afOverlap0 = afU0;
    }
  }
  else // afU1 == afV0
  {
    m_iQuantity = 1;
    afOverlap0 = afU1;
  }

  // ** finally, see what we've got
  if (m_iQuantity == 2)
  {
    // get the points of intersection
    isect = origin + dir*afOverlap0;
    isect2 = origin + dir*afOverlap1;

    // get the orientations for the segments
    OrientationType p1 = area_sign(tri[0], tri[1], isect);
    OrientationType p2 = area_sign(tri[1], tri[2], isect);
    OrientationType p3 = area_sign(tri[2], tri[0], isect);
    OrientationType q1 = area_sign(tri[0], tri[1], isect2);
    OrientationType q2 = area_sign(tri[1], tri[2], isect2);
    OrientationType q3 = area_sign(tri[2], tri[0], isect2);

    // get the point locations
    PolygonLocationType t1 = in_tri(p1, p2, p3);
    PolygonLocationType t2 = in_tri(q1, q2, q3);

    switch (t1)
    {
      case ePolygonOnVertex:
        switch (t2)
        {
          case ePolygonOnVertex: 
            return (!rel_equal(isect, isect2)) ? eSegTriEdgeOverlap : eSegTriVertex;
           
          case ePolygonOnEdge: return eSegTriEdgeOverlap;
          case ePolygonInside: return eSegTriPlanarIntersect;

          // this should never happen (but likely will, with f.p. arithmetic)
          case ePolygonOutside:
            return eSegTriEdgeOverlap;  // just a guess

          default: assert(false);
        }  
        break;
        
      case ePolygonOnEdge:
        switch (t2)
        {
          case ePolygonInside: return eSegTriInside; 
          case ePolygonOnVertex: return eSegTriEdgeOverlap;
          case ePolygonOnEdge:
            // if they are on the same edge, it is edge overlap; otherwise,
            // it is in the face
            if ((p1 == eOn && q1 == eOn) || (p2 == eOn && q2 == eOn) ||
                (p3 == eOn && q3 == eOn))
              return eSegTriEdgeOverlap;
            else
              return eSegTriInside;

          // again, this should never happen...
          case ePolygonOutside: 
            return eSegTriEdgeOverlap;  // again, just a guess
        }

      case ePolygonInside:
        switch (t2)
        {
          case ePolygonInside:
          case ePolygonOnVertex:
          case ePolygonOnEdge:
            return eSegTriInside;

          // again, this should never happen
          case ePolygonOutside:
            return eSegTriInside;  // again, just a guess

          default: assert(false);
        }

      // this case should never happen -- all of these will be guesses
      case ePolygonOutside:
        switch (t2)
        {
          case ePolygonInside: return eSegTriInside;
          case ePolygonOnVertex: return eSegTriEdgeOverlap;
          case ePolygonOnEdge:   return eSegTriEdgeOverlap;
          case ePolygonOutside: return eSegTriEdgeOverlap;

          default: assert(false);
        }

        default: assert(false);
    } 
  }
  else if (m_iQuantity == 1)
  {
    // get the point of intersection
    isect = origin + dir*afOverlap0;

    switch (in_tri(tri, isect))
    {
      case ePolygonOnVertex: return eSegTriVertex;
      case ePolygonOnEdge: return eSegTriEdge;

      // shouldn't be here..  we'll just return edge
      default: return eSegTriEdge;
    } 
  }
  else
    return eSegTriNoIntersect;

  // should not get here
  assert(false);
  return eSegTriNoIntersect;
*/
}

/// Determines the intersection between a line segment and triangle in 3D
/**
 * Adapted from O'Rourke, p. 238.
 * \param seg a line segment in 3D
 * \param t a triangle in 3D
 * \param isect contains the point of intersection, if any (on return)
 * \param isect2 contains a second point of intersection, if the intersection
 *               type is seg_tri_inside, seg_tri_edge_overlap, or 
 *               seg_tri_planar_intersect 
 *               (on return)
 * \return eSegTriInside (if the segment lies wholly within the face of
 *         the triangle), eSegTriVertex (if an endpoint of the segment 
 *         coincides with a vertex
 *         of the triangle), eSegTriEdge (if an endpoint of the segment is in
 *         the relative interior of an edge of the triangle), eSegTriFace
 *         (if an endpoint of the segment is in the relative interior of the
 *         face of the triangle), eSegTriInclVertex (if the open segment
 *         includes a vertex of the triangle), eSegTriInclEdge (if the open
 *         segment includes a point in the relative interior of an edge of the
 *         triangle), eSegTriInclFace (if the open segment includes a 
 *         point in the relative interior of a face of the triangle), 
 *         eSegTriEdgeOverlap (if the segment overlaps one edge of the
 *         triangle [or vice versa]), eSegTriPlanarIntersect (if the
 *         segment lies partially within the face of the triangle), or
 *         eSegTriNoIntersect (if there is no intersection).
 */
CompGeom::SegTriIntersectType CompGeom::intersect_seg_tri(const LineSeg3& seg, const Triangle& t, Point3d& isect, Point3d& isect2, double tol)
{
  unsigned m;

  // intersect segment w/plane; also find point of intersection (if segment
  // not coplanar)
  SegPlaneIntersectType code = intersect_seg_plane(t, seg, isect, m, tol);

  if (code == eSegPlaneToSide)
    return eSegTriNoIntersect;
  
  if (code == eSegPlaneFirst)
  {
    PolygonLocationType loc = in_tri(t, seg.first, m, tol);
    switch (loc)
    {
      case ePolygonInside:   return eSegTriFace;
      case ePolygonOnVertex: return eSegTriVertex;
      case ePolygonOnEdge:   return eSegTriEdge;
      case ePolygonOutside:  return eSegTriNoIntersect;
    }
  }

  if (code == eSegPlaneSecond)
  {
    PolygonLocationType loc = in_tri(t, seg.second, m, tol);
    switch (loc)
    {
      case ePolygonInside:   return eSegTriFace;
      case ePolygonOnVertex: return eSegTriVertex;
      case ePolygonOnEdge:   return eSegTriEdge;
      case ePolygonOutside:  return eSegTriNoIntersect;
    }

    // should never get here
    throw std::runtime_error("Undefined execution path");
  }

  if (code == eSegPlaneInPlane)
    return intersect_seg_tri_in_plane(seg, t, isect, isect2, tol);
  else
  {
    assert(code == eSegPlaneOtherIntersect); 
    try
    {
      return seg_tri_cross(seg, t, tol);
    }
    catch (NumericalException e)
    {
      return intersect_seg_tri_in_plane(seg, t, isect, isect2, tol);
    }
  }
}

/// Intersects a segment and a triangle in 3D that lie in the same plane
CompGeom::SegTriIntersectType CompGeom::intersect_seg_tri_in_plane(const LineSeg3& seg, const Triangle& t, Point3d& isect, Point3d& isect2, double tol)
{
  // project to dimension with least variance
  const unsigned X = 0, Y = 1, Z = 2;
  Point2d tri_2D[3];
  LineSeg2 seg_2D;

  // determine which coordinate has the minimum variance
  std::pair<Point3d, Point3d> aabb = t.calc_AABB();
  const double varx = aabb.second[X] - aabb.first[X];
  const double vary = aabb.second[Y] - aabb.first[Y];
  const double varz = aabb.second[Z] - aabb.first[Z];

  // project to 2D
  if (varx < vary && varx < varz)
  {
    // minimum variance in the x direction; project to yz plane
    tri_2D[0] = Point2d(t.a[Y], t.a[Z]);
    tri_2D[1] = Point2d(t.b[Y], t.b[Z]);
    tri_2D[2] = Point2d(t.c[Y], t.c[Z]);
    seg_2D.first = Point2d(seg.first[Y], seg.first[Z]);
    seg_2D.second = Point2d(seg.second[Y], seg.second[Z]);
  }
  else if (vary < varx && vary < varz)
  {
    // minimum variance in the y direction; project to xz plane
    tri_2D[0] = Point2d(t.a[X], t.a[Z]);
    tri_2D[1] = Point2d(t.b[X], t.b[Z]);
    tri_2D[2] = Point2d(t.c[X], t.c[Z]);
    seg_2D.first = Point2d(seg.first[X], seg.first[Z]);
    seg_2D.second = Point2d(seg.second[X], seg.second[Z]);
  }
  else
  {
    // minimum variance in the z direction; project to xy plane
    tri_2D[0] = Point2d(t.a[X], t.a[Y]);
    tri_2D[1] = Point2d(t.b[X], t.b[Y]);
    tri_2D[2] = Point2d(t.c[X], t.c[Y]);
    seg_2D.first = Point2d(seg.first[X], seg.first[Y]);
    seg_2D.second = Point2d(seg.second[X], seg.second[Y]);
  }

  // get the code
  Point2d isect_2D, isect2_2D;
  SegTriIntersectType code = intersect_seg_tri(seg_2D, tri_2D, isect_2D, isect2_2D, tol);
  if (code != eSegTriNoIntersect)
  {
    // determine the points of intersection in 3D
    double len_d2d = (seg_2D.second - seg_2D.first).norm();
    double t1 = (isect_2D - seg_2D.first).norm()/len_d2d;
    double t2 = (isect2_2D - seg_2D.first).norm()/len_d2d;
    Vector3d d = seg.second - seg.first;
    isect = seg.first + d*t1;
    isect2 = seg.first + d*t2;
  }
  return code;
}

/// Intersects two non-coplanar triangles in 3D and returns the line of intersection
/**
 * \param t1 the first triangle
 * \param t2 the second triangle
 * \param p1 an endpoint of the line of intersection (on return)
 * \param p2 an endpoint of the line of intersection (on return)
 * \return <b>true</b> if the triangles intersect, <b>false</b> if they don't
 *         intersect or are coplanar
 */
bool CompGeom::intersect_noncoplanar_tris(const Triangle& t1, const Triangle& t2, Point3d& p1, Point3d& p2)
{
  // compute plane equation of t1 
  Vector3d e1 = t1.b - t1.a;
  Vector3d e2 = t1.c - t1.a;
  Vector3d n1 = Vector3d::cross(e1, e2);
  double d1 = -Vector3d::dot(n1, t1.a);

  // compute signed distances to plane of verts of second tri
  double du0 = Vector3d::dot(n1, t2.a) + d1;
  double du1 = Vector3d::dot(n1, t2.b) + d1;
  double du2 = Vector3d::dot(n1, t2.c) + d1;

  // coplanarity robustness check
  if (std::fabs(du0) < NEAR_ZERO) du0 = 0.0;
  if (std::fabs(du1) < NEAR_ZERO) du1 = 0.0;
  if (std::fabs(du2) < NEAR_ZERO) du2 = 0.0;
  double du0du1 = du0*du1;
  double du0du2 = du0*du2;
  if (du0du1 > 0.0 && du0du2 > 0.0)
    return false;

  // compute plane equation of t2 
  e1 = t2.b - t2.a;
  e2 = t2.c - t2.a;
  Vector3d n2 = Vector3d::cross(e1, e2);
  double d2 = -Vector3d::dot(n2, t2.a);

  // compute signed distances to plane of verts of first tri
  double dv0 = Vector3d::dot(n2, t1.a) + d2;
  double dv1 = Vector3d::dot(n2, t1.b) + d2;
  double dv2 = Vector3d::dot(n2, t1.c) + d2;

  // coplanarity robustness check
  if (std::fabs(dv0) < NEAR_ZERO) dv0 = 0.0;
  if (std::fabs(dv1) < NEAR_ZERO) dv1 = 0.0;
  if (std::fabs(dv2) < NEAR_ZERO) dv2 = 0.0;
  double dv0dv1 = dv0*dv1;
  double dv0dv2 = dv0*dv2;
  if (dv0dv1 > 0.0 && dv0dv2 > 0.0)
    return false;
  
  // compute direction of intersection line
  Vector3d D = Vector3d::cross(n1, n2);

  // compute index to the largest component of D
  double mx = std::fabs(D[0]);
  unsigned index = 0;
  double b = std::fabs(D[1]);
  double c = std::fabs(D[2]);
  if (b > mx)
  {
    mx = b;
    index = 1;
  }
  if (c > mx)
  {
    mx = c;
    index = 2;
  }

  // this is the simplified projection onto L
  double vp0 = t1.a[index];
  double vp1 = t1.b[index];
  double vp2 = t1.c[index];
  double up0 = t2.a[index];
  double up1 = t2.b[index];
  double up2 = t2.c[index];

  // compute interval for triangle 1
  Point3d isectpointA1, isectpointA2;
  Point2d isect1;
  if (!compute_intervals_isectline(t1, vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, isect1[0], isect1[1], isectpointA1, isectpointA2))
    // triangles are coplanar
    return false;

  // compute interval for triangle 2
  Point3d isectpointB1, isectpointB2;
  Point2d isect2;
  if (!compute_intervals_isectline(t2, up0, up1, up2, du0, du1, du2, du0du1, du0du2, isect2[0], isect2[1], isectpointB1, isectpointB2))
    // triangles are coplanar
    return false;

  // sort points
  bool swapped1 = false, swapped2 = false;
  if (isect1[0] > isect1[1])
  {
    std::swap(isect1[0], isect1[1]);
    swapped1 = true;
  }
  if (isect2[0] > isect2[1])
  {
    std::swap(isect2[0], isect2[1]);
    swapped2 = true;
  }
  
  // check whether triangles intersect
  if (isect1[1] < isect2[0] || isect2[1] < isect1[0])
    return false;

  // at this point, we know that the triangles intersect
  if (isect2[0] < isect1[0])
  {
    if (!swapped1)
      p1 = isectpointA1;
    else
      p1 = isectpointA2;

    if (isect2[1] < isect1[1])
    {
      if (!swapped2)
        p2 = isectpointB2;
      else
        p2 = isectpointB1;
    }
    else
    {
      if (!swapped1)
        p2 = isectpointA2;
      else
        p2 = isectpointA1;  
    }
  }
  else
  {
    if (!swapped2)
      p1 = isectpointB1;
    else
      p1 = isectpointB2;

    if (isect2[1] > isect1[1])
    {
      if (!swapped1)
        p2 = isectpointA2;
      else
        p2 = isectpointA1;
    }
    else
    {
      if (!swapped2)
        p2 = isectpointB2;
      else
        p2 = isectpointB1;
    }
  }

  return true;
}

/// Function for testing whether a point is in a triangle
bool CompGeom::point_in_tri(const Point3d& p, const Triangle& t, double tol)
{
  return (same_side(p, t.a, t.b, t.c, tol) &&
          same_side(p, t.b, t.a, t.c, tol) &&
          same_side(p, t.c, t.a, t.b, tol));
}

/// Utility function for testing whether a point is in a triangle
bool CompGeom::same_side(const Point3d& p1, const Point3d& p2, const Point3d& a, const Point3d& b, double tol)
{
  assert(tol >= 0.0);
  Vector3d cp1 = Vector3d::cross(b-a, p1-a);
  Vector3d cp2 = Vector3d::cross(b-a, p2-a);
  return (Vector3d::dot(cp1,cp2) >= tol);
}  

/**
 * Utility function for use by intersect_coplanar_tris()
 * Adapted from ODE.
 * \return <b>true</b> normally, <b>false</b> if triangles are coplanar
 */
bool CompGeom::compute_intervals_isectline(const Triangle& t, double vv0, double vv1, double vv2, double d0, double d1, double d2, double d0d1, double d0d2, double& isect0, double& isect1, Point3d& isectpoint0, Point3d& isectpoint1)
{
  if (d0d1 > 0.0)
    // here we know that d0d2 <= 0.0; that is d0, d1 are on the same side,
    // d2 or on the other or on other plane
    isect2X(t.c, t.a, t.b, vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpoint0, isectpoint1);
  else if (d0d2 > 0.0)
    // here we know that d0d1 <= 0.0 or that d0 != 0.0
    isect2X(t.b, t.a, t.c, vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpoint0, isectpoint1);
  else if (d1*d2 > 0.0 || d0 != 0.0)
    isect2X(t.a, t.b, t.c, vv0, vv1, vv2, d0, d1, d2, isect0, isect1, isectpoint0, isectpoint1);
  else if (d1 != 0.0)
    isect2X(t.b, t.a, t.c, vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpoint0, isectpoint1);
  else if (d2 != 0.0)
    isect2X(t.c, t.a, t.b, vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpoint0, isectpoint1);
  else
    // triangles are coplanar
    return false;

  return true;
}

/**
 * Utility function for use by intersect_coplanar_tris()
 * Adapted from ODE.
 */
void CompGeom::isect2X(const Point3d& vtx0, const Point3d& vtx1, const Point3d& vtx2, double vv0, double vv1, double vv2, double d0, double d1, double d2, double& isect0, double& isect1, Point3d& isectpoint0, Point3d& isectpoint1)
{
  double tmp = d0/(d0-d1);
  isect0 = vv0 + (vv1-vv0)*tmp;
  Vector3d diff = (vtx1 - vtx0) * tmp;
  isectpoint0 = diff + vtx0;
  tmp = d0/(d0-d2);
  isect1 = vv0 + (vv2-vv0)*tmp;
  diff = (vtx2 - vtx0) * tmp;
  isectpoint1 = diff + vtx0;
}

/// Computes the projection matrix to convert from 3D to 2D
Matrix3d CompGeom::calc_3D_to_2D_matrix(const Vector3d& normal)
{
  const unsigned X = 0, Y = 1, Z = 2;

  Vector3d v1, v2;
  Vector3d::determine_orthonormal_basis(normal, v1, v2);
  Matrix3d R;
  R.set_row(X, Origin3d(v1));
  R.set_row(Y, Origin3d(v2));
  R.set_row(Z, Origin3d(normal));

  return R;
}

/// Computes the offset in converting a 3D vector to 2D
/**
 * \param v a vector in 3D
 * \param R the 3D-to-2D projection matrix
 * \return the offset necessary to add to the Z-coordinate of a 3D vector 
 *         converted from 2D using the transpose of R
 */
double CompGeom::determine_3D_to_2D_offset(const Origin3d& o, const Matrix3d& R)
{
  const unsigned Z = 2;
  return (R * o)[Z];
}

/// Determines the location of a point on a triangle
/**
 * \pre it has been determined that p lies in t's plane
 */
CompGeom::PolygonLocationType CompGeom::in_tri(const Triangle& t, const Point3d& p, double tol)
{
  // find the plane coefficients
  unsigned big_dim;
  double dot;
  Vector3d N;
  find_plane_coeff(t, N, big_dim, dot);

  // determine the intersection type
  return in_tri(t, p, big_dim, NEAR_ZERO);

/*
  SegTriIntersectType itype = in_tri(t, p, big_dim, NEAR_ZERO);
  switch (itype)
  {
    case eSegTriNoIntersect:
      return ePolygonOutside;

    case eSegTriFace:
    case eSegTriInclFace:
    case eSegTriInside:
      return ePolygonInside;

    case eSegTriVertex:
    case eSegTriInclVertex:
      return ePolygonOnVertex;

    case eSegTriEdge:
    case eSegTriInclEdge:
    case eSegTriEdgeOverlap:
    case eSegTriPlanarIntersect:
      return ePolygonOnEdge;

    default:
      assert(false);
  }

  std::cerr << "CompGeom::in_tri() - should never have gotten here!" << std::endl;
  assert(false);
  return ePolygonOutside;
*/
}

/// Determines the location of a point on a triangle
/**
 * Adapted from O'Rourke, p. 233.
 * \param t a 3D triangle
 * \param p a 3D point
 * \param skip_dim the dimension to skip
 * \return on_vertex if p coincides with a vertex of t, on_edge if p is in the 
 *         relative interior of an edge of t, on_face if p is in the relative
 *         interior of a face of t, or outside_tri if p does not intersect t
 */
CompGeom::PolygonLocationType CompGeom::in_tri(const Triangle& t, const Point3d& p, unsigned int skip_dim, double tol)
{
  const unsigned THREE_D = 3, TRI_PTS = 3;
  Point2d pp;
  Triangle tri_2D;
  Point2d Tp[TRI_PTS];
  assert(tol >= (double) 0.0);

  // if the skipped dimension was not specified, find it
  if (skip_dim >= THREE_D)
  {
    Vector3d N = normal_vec(t.a, t.b, t.c);
    long double biggest = 0.0;
    for (unsigned i=0; i< THREE_D; i++)
    {
      long double t = std::fabs(N[i]);
      if (t > biggest)
      {
        biggest = t;
        skip_dim = i;
      }
    }
  }
  
  for (unsigned i=0, j=0; i< THREE_D; i++)
  {
    if (i == skip_dim)
      continue;
    pp[j] = p[i];
    for (unsigned k=0; k< THREE_D; k++)
      Tp[k][j] = t.get_vertex(k)[i];
    j++;
  }

  return in_tri(Tp, pp, tol);
}

/// Determines the intersection of a point on a 2D triangle
/**
 * Adapted from O'Rourke, p. 236.
 * \param t the 2D triangle 
 * \param p a 2D point
 * \param tol the tolerance with which to test for zero
 * \return ePolygonOnVertex if p coincides with a vertex of t, ePolygonOnEdge
 *         if p is in the relative interior of an edge of t, ePolygonInside if 
 *         p is in the relative interior of t, or ePolygonOutside otherwise
 */
CompGeom::PolygonLocationType CompGeom::in_tri(const Point2d t[3], const Point2d& p, double tol)
{
  assert(tol >= 0.0);

  const Point2d& ta = t[0];
  const Point2d& tb = t[1];
  const Point2d& tc = t[2];

  // get the signed areas
  OrientationType ori0 = area_sign(p, ta, tb, tol);
  OrientationType ori1 = area_sign(p, tb, tc, tol);
  OrientationType ori2 = area_sign(p, tc, ta, tol);

  // call the main method
  return in_tri(ori0, ori1, ori2);  
}

/// Determines the intersection of a point on a 2D triangle
/**
 * Adapted from O'Rourke, p. 236.
 * \param ta the first vertex of a 2D triangle
 * \param tb the second vertex of a 2D triangle
 * \param tc the third vertex of a 2D triangle
 * \param p a 2D point
 * \param NEAR_ZERO the tolerance with which to test for zero
 * \return ePolygonOnVertex if p coincides with a vertex of t, ePolygonOnEdge
 *         if p is in the relative interior of an edge of t, ePolygonInside if 
 *         p is in the relative interior of t, or ePolygonOutside otherwise
 */
CompGeom::PolygonLocationType CompGeom::in_tri(OrientationType ori0, OrientationType ori1, OrientationType ori2)
{
  if ((ori0 == eOn && ori1 == eLeft && ori2 == eLeft) ||
      (ori1 == eOn && ori0 == eLeft && ori2 == eLeft) ||
      (ori2 == eOn && ori0 == eLeft && ori1 == eLeft))
    return ePolygonOnEdge;

  if ((ori0 == eOn && ori1 == eRight && ori2 == eRight) ||
      (ori1 == eOn && ori0 == eRight && ori2 == eRight) ||
      (ori2 == eOn && ori0 == eRight && ori1 == eRight))
    return ePolygonOnEdge;

  if ((ori0 == eLeft && ori1 == eLeft && ori2 == eLeft) ||
      (ori0 == eRight && ori1 == eRight && ori2 == eRight))
    return ePolygonInside;

  if ((ori0 == eOn && ori1 == eOn) ||
      (ori0 == eOn && ori2 == eOn) ||
      (ori1 == eOn && ori2 == eOn))
    return ePolygonOnVertex;

  // some cases are clearly impossible-- we will handle those here

  // we're assuming that one of these areas is close to zero (right on the
  // zero side of our tolerance) -- it could be the case where two or even 
  // all three areas are off 
  if (ori0 == eOn && ori1 == eOn && ori2 == eOn)
    return ePolygonOnVertex;

  // must be outside
  return ePolygonOutside;
}

/// Gets the coefficients for a plane
/**
 * Utility method for intersect_seg_plane().  Adapted from O'Rourke, p. 229.
 * \param t the triangle that represents the plane
 * \param N contains the normal vector to the plane (on return)
 * \param m contains the biggest dimension (0-2) of the normal (on return)
 * \param dot contains the dot-product of the normal vector and the first
 *        triangle vertex (on return)
 */
void CompGeom::find_plane_coeff(const Triangle& t, Vector3d& N, unsigned int& m, double& dot)
{
  const unsigned int THREE_D = 3;
  N = t.calc_normal();
  dot = t.calc_offset(N);

  // find the largest component of N
  long double biggest = -1;
  for (unsigned i=0; i< THREE_D; i++)
  {
    long double t = std::fabs(N[i]);
    if (t > biggest)
    {
      biggest = t;
      m = i;
    }
  }
}

/// Determines the cross-product of (b-a) and (c-a)
Vector3d CompGeom::normal_vec(const Point3d& a, const Point3d& b, const Point3d& c)
{
  unsigned X = 0, Y = 1, Z = 2;
  
  Vector3d N;
  N[X] = (c[Z] - a[Z]) * (b[Y] - a[Y]) - (b[Z] - a[Z]) * (c[Y] - a[Y]);
  N[Y] = (b[Z] - a[Z]) * (c[X] - a[X]) - (b[X] - a[X]) * (c[Z] - a[Z]);
  N[Z] = (b[X] - a[X]) * (c[Y] - a[Y]) - (b[Y] - a[Y]) * (c[X] - a[X]);

  return N;
}

/// Determines whether three points are collinear
bool CompGeom::collinear(const Point2d& a, const Point2d& b, const Point2d& c, double tol)
{
  return area_sign(a, b, c, tol) == eOn;
}

/// Determines whether three points are collinear
bool CompGeom::collinear(const Point3d& a, const Point3d& b, const Point3d& c, double tol)
{  
  const unsigned X = 0, Y = 1, Z = 2;
  assert(tol >= 0.0);

   return (rel_equal((c[Z]-a[Z])*(b[Y]-a[Y]), (b[Z]-a[Z])*(c[Y]-a[Y]), tol) &&
          rel_equal((b[Z]-a[Z])*(c[X]-a[X]), (b[X]-a[X])*(c[Z]-a[Z]), tol) && 
          rel_equal((b[X]-a[X])*(c[Y]-a[Y]), (b[Y]-a[Y])*(c[X]-a[X]), tol));
}

/// Gets the volume of a tetrahedron composed of vertices a, b, c, d
long double CompGeom::volume(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  long double ax = a[X] - d[X];
  long double ay = a[Y] - d[Y];
  long double az = a[Z] - d[Z];
  long double bx = b[X] - d[X];
  long double by = b[Y] - d[Y];
  long double bz = b[Z] - d[Z];
  long double cx = c[X] - d[X];
  long double cy = c[Y] - d[Y];
  long double cz = c[Z] - d[Z];
  return ax * (by*cz-bz*cy) + ay * (bz*cx-bx*cz) + az * (bx*cy-by*cx);
}

/// Gets the volume of a tetrahedron composed of vertices a, b, c, d
long double CompGeom::volume(const Origin3d& a, const Origin3d& b, const Origin3d& c, const Origin3d& d)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  long double ax = a[X] - d[X];
  long double ay = a[Y] - d[Y];
  long double az = a[Z] - d[Z];
  long double bx = b[X] - d[X];
  long double by = b[Y] - d[Y];
  long double bz = b[Z] - d[Z];
  long double cx = c[X] - d[X];
  long double cy = c[Y] - d[Y];
  long double cz = c[Z] - d[Z];
  return ax * (by*cz-bz*cy) + ay * (bz*cx-bx*cz) + az * (bx*cy-by*cx);
}


/// Gets the sign of the volume of a tetrahedron composed of vertices a, b, c, d
/**
 * \returns -1 if d is "visible" from the plane of (a,b,c); 0 if the a,b,c, and
 *          d are coplanar, and +1 otherwise
 */
CompGeom::VisibilityType CompGeom::volume_sign(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d, double tol)
{
  assert(tol >= 0.0);

  long double v = volume(a, b, c, d);
  if (v < -tol)
    return eVisible;
  else if (v > tol)
    return eInvisible;
  else
    return eCoplanar;
}

/// Gets the sign of the volume of a tetrahedron composed of vertices a, b, c, d
/**
 * \returns -1 if d is "visible" from the plane of (a,b,c); 0 if the a,b,c, and
 *          d are coplanar, and +1 otherwise
 */
CompGeom::VisibilityType CompGeom::volume_sign(const Origin3d& a, const Origin3d& b, const Origin3d& c, const Origin3d& d, double tol)
{
  assert(tol >= 0.0);

  long double v = volume(a, b, c, d);
  if (v < -tol)
    return eVisible;
  else if (v > tol)
    return eInvisible;
  else
    return eCoplanar;
}

/// Gets the area of two 2D vectors with respect to an arbitrary center
long double CompGeom::area(const Point2d& a, const Point2d& b, const Point2d& c)
{
  const unsigned X = 0, Y = 1;
  long double a1 = (b[X] - a[X]) * (c[Y] - a[Y]);
  long double a2 = (c[X] - a[X]) * (b[Y] - a[Y]);
  return a1 - a2;
}

/// Gets the sign of the area of two vectors with respect to an arbitrary center
CompGeom::OrientationType CompGeom::area_sign(const Point2d& a, const Point2d& b, const Point2d& c, double tol)
{
  assert(tol >= 0.0);  

  const unsigned X = 0, Y = 1;
  long double a1 = (b[X] - a[X]) * (c[Y] - a[Y]);
  long double a2 = (c[X] - a[X]) * (b[Y] - a[Y]);

  // see whether a1 and a2 are relatively equal
//  if (rel_equal(a1, a2, tol))
//    return eOn;

  // get the area
  double x = a1 - a2;

if (x < tol && x > -tol)
  return eOn;

  // return the proper sign
  return (x > 0.0) ? eLeft : eRight;
}

/// Determines whether two lines/rays/segments in 2D intersect
/**
 * Taken from O'Rourke, p. 224.
 * \param o1 the origin of the first line segment/ray (or a point on the line)
 * \param dir1 the direction of the first line
 * \param min1 the minimum parameter of the first line segment/ray/line
 * \param max1 the maximum parameter of the first line segment/ray/line
 * \param o2 the origin of the first line segment/ray (or a point on the line)
 * \param dir2 the direction of the second line
 * \param min2 the minimum parameter of the second line segment/ray/line
 * \param max2 the maximum parameter of the second line segment/ray/line
 * \param s1 on return of eLineLineIntersect/eLineLineVertex/eLineLineEdge
 *        contains the parameter on the first line of the first point of
 *        intersection (o1 + dir1*s1)
 * \param s2 on return of eLineLineEdge contains the parameter on the first
 *        line of the second point of intersection (o1 + dir1*s2)
 * \return the intersection type (edge, vertex, true_intersect, no_intersect); 
 *         an edge intersection indicates that the lines/rays/segments 
 *         collinearly overlap; vertex intersection indicates that an endpoint 
 *         of one ray/segment is on the other line/ray/segment, but edge 
 *         intersection does not hold; true_intersect indicates that the 
 *         lines/rays/segments intersect at a single point.
 */
CompGeom::LineLineIntersectType CompGeom::intersect_lines(const Point2d& o1, const Vector2d& dir1, double min1, double max1, const Point2d& o2, const Vector2d& dir2, double min2, double max2, double& s, double& t)
{
  const int X = 0, Y = 1;
  LineLineIntersectType code = eLineLineNoIntersect;
  const double INF = std::numeric_limits<double>::max();
  
  // make sure that line parameters are correct (NOTE: if we wish, we can
  // just swap the parameters)
  assert(min1 <= max1);
  assert(min2 <= max2);

  // get two points on the line/ray/segment
  const Point2d& a = o1;
  const Vector2d b = o1+dir1;
  const Point2d& c = o2;
  const Vector2d d = o2+dir2;

  // determine whether the lines/rays/segments are parallel  
  long double denom = a[X] * (d[Y] - c[Y]) + b[X] * (c[Y] - d[Y]) + 
          d[X] * (b[Y] - a[Y]) + c[X] * (a[Y] - b[Y]);

  // if the denominator is zero, the segments are parallel; handle separately
  if (denom < NEAR_ZERO && denom > -NEAR_ZERO)
  {
    // compute parameters of second line/ray/seg w.r.t the first
    if (min2 == -INF)
      s = -INF;
    else
      s = (o1 - (o2 + dir2*min2)).norm_sq() / dir1.norm_sq();
    if (max2 == INF)
      t = INF;
    else
      t = (o1 - (o2 + dir2*max2)).norm_sq() / dir1.norm_sq();
    if (s > t)
      std::swap(s,t);

    // determine overlap of [s,t] and [min1,max1]
    s = std::max(s, min1);
    t = std::min(t, max1);
    if (t < s)
      return eLineLineNoIntersect;
    else
      return eLineLineEdge;  
  }

  // calculate the numerator
  long double num = a[X] * (d[Y] - c[Y]) + c[X] * (a[Y] - d[Y]) + d[X] * (c[Y] - a[Y]);

  // lines intersect; determine parameter of first line
  if (std::fabs(num) < NEAR_ZERO || std::fabs(num-denom) < NEAR_ZERO)
    code = eLineLineVertex;
  s = num / denom;

  // determine parameter of second line
  num = -(a[X] * (c[Y] - b[Y]) + b[X] * (a[Y] - c[Y]) + c[X] * (b[Y] - a[Y]));
  if ((num < NEAR_ZERO && num > -NEAR_ZERO) || (std::fabs(num -denom) < NEAR_ZERO)) 
    code = eLineLineVertex;
  t = num / denom;

  // check that they are within the ray/segment
  if ((min1 < s) && (s < max1) && (min2 < t) && (t < max2))
    code = eLineLineIntersect;
  else if ((min1 > s) || (s > max1) || (min2 > t) || (t > max2))
    code = eLineLineNoIntersect;

  return code;  
}

/// Determines whether two line segments in 2D intersect
/**
 * Taken from O'Rourke, p. 224.
 * \param s1 the points of the first line segment
 * \param s2 the points of the second line segment
 * \param isect the point of intersection (if any)
 * \param isect2 the second point of intersection (for edge intersection)
 * \return the intersection type (edge, vertex, true_intersect, no_intersect); 
 *         an edge intersection indicates that the segments collinearly overlap;
 *         vertex intersection indicates that an endpoint of
 *         one segment is on the other segment, but edge intersection does not
 *         hold; true_intersect indicates that the segments intersect at a
 *         single point.
 */
CompGeom::SegSegIntersectType CompGeom::intersect_segs(const LineSeg2& s1, const LineSeg2& s2, Point2d& isect, Point2d& isect2)
{
  const int X = 0, Y = 1;
  SegSegIntersectType code = eSegSegNoIntersect;
  
  // get the individual points 
  const Point2d& a = s1.first;
  const Point2d& b = s1.second;
  const Point2d& c = s2.first;
  const Point2d& d = s2.second;
  
  long double denom = a[X] * (d[Y] - c[Y]) + b[X] * (c[Y] - d[Y]) + 
          d[X] * (b[Y] - a[Y]) + c[X] * (a[Y] - b[Y]);

  // if the denominator is near, the segments are parallel; handle separately
  if (denom < NEAR_ZERO && denom > -NEAR_ZERO)
    return get_parallel_intersect_type(s1, s2, isect, isect2);

  // calculate the numerator
  long double num = a[X] * (d[Y] - c[Y]) + c[X] * (a[Y] - d[Y]) + d[X] * (c[Y] - a[Y]);

  if (std::fabs(num) < NEAR_ZERO || std::fabs(num-denom) < NEAR_ZERO)
    code = eSegSegVertex;
  long double s = num / denom;

  num = -(a[X] * (c[Y] - b[Y]) + b[X] * (a[Y] - c[Y]) + c[X] * (b[Y] - a[Y]));
  if ((num < NEAR_ZERO && num > -NEAR_ZERO) || (std::fabs(num -denom) < NEAR_ZERO)) 
    code = eSegSegVertex;
  
  long double t = num / denom;
  if ((0.0 < s) && (s < 1.0) && (0.0 < t) && (t < 1.0))
    code = eSegSegIntersect;
  else if ((0.0 > s) || (s > 1.0) || (0.0 > t) || (t > 1.0))
    code = eSegSegNoIntersect;

  isect[X] = a[X] + s * (b[X] - a[X]);
  isect[Y] = a[Y] + s * (b[Y] - a[Y]);
  
  return code;  
}

/// Function for determining the type of parallel intersection
/**
 * Taken from O'Rourke, p. 263.
 */
CompGeom::SegSegIntersectType CompGeom::get_parallel_intersect_type(const LineSeg2& s1, const LineSeg2& s2, Point2d& isect, Point2d& isect2)
{
  // get the individual points 
  const Point2d& a = s1.first;
  const Point2d& b = s1.second;
  const Point2d& c = s2.first;
  const Point2d& d = s2.second;
  
  // check for collinearity
  if (area_sign(a,b,c) == eRight)
    return eSegSegNoIntersect;

  if (between(a, b, c) && between(a, b, d))
  {
    isect = c;
    isect2 = d;
    return eSegSegEdge;
  }

  if (between(c, d, a) && between(c, d, b))
  {
    isect = a;
    isect2 = b;
    return eSegSegEdge;
  }

  if (between(a, b, c) && between(c, d, b))
  {
    isect = c;
    isect2 = b;
    return eSegSegEdge;
  }

  if (between(a, b, c) && between(c, d, a))
  {
    isect = c;
    isect2 = a;
    return eSegSegEdge;
  }

  if (between(a, b, d) && between(c, d, b))
  {
    isect = d;
    isect2 = b;
    return eSegSegEdge;
  }
  
  if (between(a, b, d) && between(c, d, a))
  {
    isect = d;
    isect2 = a;
    return eSegSegEdge;
  }
  
  return eSegSegNoIntersect;
}

/// Utility function used by get_parallel_intersect_type() and triangulate_polygon_2D()
bool CompGeom::between(const Point2d& a, const Point2d& b, const Point2d& c, double tol)
{
  const int X = 0, Y = 1;
  assert(tol >= 0.0);

  if (!collinear(a, b, c, tol))
    return false;
  
  if (std::fabs(a[X] - b[X]) < tol)
    return ((a[X] <= c[X]) && (c[X] <= b[X])) || ((a[X] >= c[X]) && (c[X] >= b[X]));
  else
    return ((a[Y] <= c[Y]) && (c[Y] <= b[Y])) || ((a[Y] >= c[Y]) && (c[Y] >= b[Y]));
}

/// Converts a point in 2D to a point in 3D
/**
 * \param p the point to transform
 * \param RT the matrix that projects 2D points to 3D
 * \param offset the last dimension of the 3D point is set to this value
 */
Origin3d CompGeom::to_3D(const Origin2d& p, const Matrix3d& RT, double offset)
{
  const unsigned X = 0, Y = 1;

  // transform the vector
  return RT * Origin3d(p[X], p[Y], offset);
}

/// Converts a 2D point to a 3D point 
/**
 * \param p the point to transform
 * \param RT the matrix that projects 2D points to 3D
 * \param offset the last dimension of the 3D point is set to this value
 */
Origin3d CompGeom::to_3D(const Origin2d& p, const Matrix3d& RT)
{
  const unsigned X = 0, Y = 1;

  // get elements of RT
  const double rxx = RT.xx();
  const double rxy = RT.xy();
  const double ryx = RT.yx();
  const double ryy = RT.yy();
  const double rzx = RT.zx();
  const double rzy = RT.zy();
  
  // get elements of p
  const double vx = p[X];
  const double vy = p[Y];

  // transform the point 
  return Origin3d(rxx*vx + rxy*vy, ryx*vx + ryy*vy, rzx*vx + rzy*vy);
}

/// Converts a 3D point to a 2D point 
/**
 * \param p the point to transform
 * \param R the matrix that projects 3D points to 2D
 * \note the last dimension is stripped
 */
Origin2d CompGeom::to_2D(const Point3d& p, const Matrix3d& R)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get elements of R
  const double rxx = R.xx();
  const double rxy = R.xy();
  const double rxz = R.xz();
  const double ryx = R.yx();
  const double ryy = R.yy();
  const double ryz = R.yz();
  
  // get elements of p
  const double vx = p[X];
  const double vy = p[Y];
  const double vz = p[Z];

  // transform the vector without the last row
  return Origin2d(rxx*vx + rxy*vy + rxz*vz, ryx*vx + ryy*vy + ryz*vz);
}

/**
 * Utility function for intersect_tri_tri()
 */
bool CompGeom::test_edge_edge(const Point3d& p, const Point3d& a, const Point3d& b, double Ax, double Ay, unsigned i0, unsigned i1)
{
  double Bx = a[i0] - b[i0];
  double By = a[i1] - b[i1];
  double Cx = p[i0] - a[i0];
  double Cy = p[i1] - a[i1];
  double f = Ay*Bx - Ax*By;
  double d = By*Cx - Bx*Cy;
  if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f))
  {
    double e = Ax*Cy-Ay*Cx;
    if (f > 0)
    {
      if (e >= 0 && e <= f)
        return true;
    }
    else if (e <= 0 && e >= f)
      return true;
  }

  return false;
}

/**
 * Utility function for intersect_tri_tri()
 */
bool CompGeom::test_edge_tri(const Point3d& a, const Point3d& b, const Triangle& t, unsigned i0, unsigned i1)
{
  double Ax = b[i0] - a[i0];
  double Ay = b[i1] - a[i1];
  if (test_edge_edge(a, t.a, t.b, Ax, Ay, i0, i1))
    return true;
  if (test_edge_edge(a, t.b, t.c, Ax, Ay, i0, i1))
    return true;
  if (test_edge_edge(a, t.c, t.a, Ax, Ay, i0, i1))
    return true;

  return false;
}

/**
 * Utility function for intersect_tri_tri()
 */
bool CompGeom::test_point_in_tri(const Point3d& p, const Triangle& t, unsigned i0, unsigned i1)
{
  double a = t.b[i1] - t.a[i1];
  double b = -(t.b[i0] - t.a[i0]);
  double c = -a*t.a[i0] - b*t.a[i1];
  double d0 = a*p[i0] + b*p[i1]+c;

  a=t.c[i1]-t.b[i1];
  b=-(t.c[i0]-t.b[i0]);
  c=-a*t.b[i0]-b*t.b[i1];
  double d1=a*p[i0]+b*p[i1]+c;

  a=t.a[i1]-t.c[i1];
  b=-(t.a[i0]-t.c[i0]);
  c=-a*t.c[i0]-b*t.c[i1];
  double d2=a*p[i0]+b*p[i1]+c;

  if(d0*d1 > 0.0 && d0*d2 > 0.0) 
    return true;

  return false;
}

// NOTE: this function is used by the disabled version for query_intersect_tri_tri()
/**
 * Utility function for query_intersect_tri_tri()
 */
bool CompGeom::test_coplanar_tri_tri(const Vector3d& N, const Triangle& t1, const Triangle& t2)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Vector3d A;
  unsigned i0,i1;

  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  A[X]=std::fabs(N[X]);
  A[Y]=std::fabs(N[Y]);
  A[Z]=std::fabs(N[Z]);
  if(A[X]>A[Y])
  {
    if(A[X]>A[Z])  
    {
      i0=Y;      /* A[0] is greatest */
      i1=Z;
    }
    else
    {
      i0=X;      /* A[2] is greatest */
      i1=Y;
    }
  }
  else   /* A[0]<=A[1] */
  {
    if(A[Z]>A[Y])
    {
      i0=X;      /* A[2] is greatest */
      i1=Y;                                           
    }
    else
    {
      i0=X;      /* A[1] is greatest */
      i1=Z;
    }
  }               
                
  /* test all edges of triangle 1 against the edges of triangle 2 */
  if (test_edge_tri(t1.a, t1.b, t2, i0, i1))
    return true;
  if (test_edge_tri(t1.b, t1.c, t2, i0, i1))
    return true;
  if (test_edge_tri(t1.c, t1.a, t2, i0, i1))
    return true;

  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  if (test_point_in_tri(t1.a, t2, i0, i1) || test_point_in_tri(t2.a, t1, i0, i1))
    return true;

  return false;
}

// NOTE: This version has been disabled b/c it sometimes yields incorrect
// results (original code from Moeller does as well).  Using slower version
// based on separating axes.
/**
 * Determines whether two triangles intersect.
 * \return <b>true</b> if the triangles intersect
 * \note adapted from [Moeller, 1997]
 */
/*
bool CompGeom::query_intersect_tri_tri(const Triangle& t1, const Triangle& t2)
{
  double a, b, c, d, e, f, x0, x1, y0, y1;
  double isect1[2], isect2[2];

  // compute plane equation of triangle(*t1.a,*t1.b,*t1.c)
  Vector3 E1 = t1.b - t1.a;
  Vector3 E2 = t1.c - t1.a;
  Vector3 N1 = Vector3::cross(E1, E2);
  double d1=-Vector3::dot(N1,t1.a);
  // plane equation 1: N1.X+d1=0

  // put *t2.a,*t2.b,*t2.c into plane equation 1 to compute signed distances to the plane
  double du0=Vector3::dot(N1,t2.a)+d1;
  double du1=Vector3::dot(N1,t2.b)+d1;
  double du2=Vector3::dot(N1,t2.c)+d1;

  // coplanarity robustness check 
  if (std::fabs(du0) < NEAR_ZERO) du0=0.0;
  if (std::fabs(du1) < NEAR_ZERO) du1=0.0;
  if (std::fabs(du2) < NEAR_ZERO) du2=0.0;
  double du0du1=du0*du1;
  double du0du2=du0*du2;

  // same sign on all of them + not equal 0 ? 
  if (du0du1 > 0.0 && du0du2 > 0.0) 
    return false;                    // no intersection occurs 

  // compute plane of triangle (*t2.a,*t2.b,*t2.c) 
  E1 = t2.b - t2.a;
  E2 = t2.c - t2.a;
  Vector3 N2 = Vector3::cross(E1, E2);
  double d2=-Vector3::dot(N2,t2.a);
  // plane equation 2: N2.X+d2=0 

  // put *t1.a,*t1.b,*t1.c into plane equation 2 
  double dv0 = Vector3::dot(N2,t1.a)+d2;
  double dv1 = Vector3::dot(N2,t1.b)+d2;
  double dv2 = Vector3::dot(N2,t1.c)+d2;

  if (std::fabs(dv0) < NEAR_ZERO) dv0=0.0;
  if (std::fabs(dv1) < NEAR_ZERO) dv1=0.0;
  if (std::fabs(dv2) < NEAR_ZERO) dv2=0.0;

  double dv0dv1=dv0*dv1;
  double dv0dv2=dv0*dv2;
  
  // same sign on all of them + not equal 0 ? 
  if (dv0dv1 > 0.0 && dv0dv2 > 0.0) 
    return false; 

  // compute direction of intersection line 
  Vector3 D = Vector3::cross(N1,N2);

  // compute and index to the largest component of D 
  double max=std::fabs(D[0]);
  unsigned index=0;
  double bb=std::fabs(D[1]);
  double cc=std::fabs(D[2]);
  if (bb>max) max=bb,index=1;
  if (cc>max) max=cc,index=2;

  // this is the simplified projection onto L
  double vp0=t1.a[index];
  double vp1=t1.b[index];
  double vp2=t1.c[index];

  double up0=t2.a[index];
  double up1=t2.b[index];
  double up2=t2.c[index];

  // compute interval for triangle 1 
  if (dv0dv1 > 0.0)
  {
    // here we know that D0D2<=0.0  
    // that is D0, D1 are on the same side, D2 on other or on the plane  
    a=vp0; 
    b=(vp0-vp2)*dv2; 
    c=(vp1-vp2)*dv2; 
    x0=dv2-dv0; 
    x1=dv2-dv1; 
  }
  else if (dv0dv2 > 0.0)
  { 
    // here we know that d0d1<=0.0  
    a=vp1; 
    b=(vp0-vp1)*dv1; 
    c=(vp2-vp1)*dv1; 
    x0 = dv1-dv0; 
    x1 = dv1-dv2; 
  } 
  else if (dv1*dv2 > 0.0 || dv0 !=0.0) 
  { 
    // here we know that d0d1<=0.0 or that D0!=0.0  
    a = vp0; 
    b = (vp1-vp0)*dv0; 
    c=(vp2-vp0)*dv0; 
    x0=dv0-dv1; 
    x1=dv0-dv2; 
  } 
  else if(dv1!=0.0) 
  { 
    a=vp1; 
    b=(vp0-vp1)*dv1; 
    c=(vp2-vp1)*dv1; 
    x0=dv1-dv0; 
    x1=dv1-dv2; 
  } 
  else if(dv2!=0.0) 
  { 
    a=vp2; 
    c=(vp0-vp2)*dv2; 
    c=(vp1-vp2)*dv2; 
    x0=dv2-dv0; 
    x1=dv2-dv1; 
  } 
  else 
  { 
    // triangles are coplanar  
    return test_coplanar_tri_tri(N1,t1,t2); 
  } 

  // compute interval for triangle 2 
  if (du0du1 > 0.0)
  {
    // here we know that D0D2<=0.0 
    // that is D0, D1 are on the same side, D2 on the other on the plane 
    d = up0; 
    e = (up0-up2)*du2; 
    f = (up1-up2)*du2; 
    y0 = du2-du0; 
    y1 = du2-du1;
  }
  else if (du0du2 > 0.0)
  {
    // here we know that d0d1<=0.0 
    d = up1; 
    e = (up0-up1)*du1; 
    f = (up2-up1)*du1; 
    y0 = du1-du0; 
    y1 = du1-du2;
  }
  else if (du1*du2 > 0.0 || du0 != 0.0)
  {
    // here we know that d0d1<=0.0 or that D0!=0.0 
    d = up0; 
    e = (up1-up0)*du0; 
    f = (up2-up0)*du0; 
    y0 = du0-du1; 
    y1 = du0-du2;
  }
  else if (du1 !=0.0)
  {
    d = up1; 
    e = (up0-up1)*du1; 
    f = (up2-up1)*du1; 
    y0 = du1-du0; 
    y1 = du1-du2;
  }
  else if(du2 !=0.0)
  {
    d=up2; 
    e=(up0-up2)*du2; 
    f=(up1-up2)*du2; 
    y0=du2-du0; 
    y1=du2-du1; 
  } 
  else 
  {
    // triangles are coplanar  
    return test_coplanar_tri_tri(N1,t1,t2); 
  }

  double xx=x0*x1;
  double yy=y0*y1;
  double xxyy=xx*yy;

  double tmp=a*xxyy;
  isect1[0]=tmp+b*x1*yy;
  isect1[1]=tmp+c*x0*yy;

  tmp=d*xxyy;
  isect2[0]=tmp+e*xx*y1;
  isect2[1]=tmp+f*xx*y0;

  if (isect1[0] > isect1[1])
    std::swap(isect1[0], isect1[1]);
  if (isect2[0] > isect2[1])
    std::swap(isect2[0], isect2[1]);

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return false;

  return true;
}
*/

/// Utility method for separating axis test version of query_intersect_tri_tri()
void CompGeom::project_onto_axis(const Triangle& rkTri, const Vector3d& rkAxis, double& rfMin, double& rfMax)
{
  double fdot0 = rkAxis.dot(rkTri.a);
  double fdot1 = rkAxis.dot(rkTri.b);
  double fdot2 = rkAxis.dot(rkTri.c);
  rfMin = fdot0;
  rfMax = rfMin;

  if (fdot1 < rfMin)
    rfMin = fdot1;
  else if (fdot1 > rfMax)
    rfMax = fdot1;

  if (fdot2 < rfMin)
    rfMin = fdot2;
  else if (fdot2 > rfMax)
    rfMax = fdot2;
}

/**
 * Determines whether two triangles intersect.
 * \return <b>true</b> if the triangles intersect
 * \note adapted from www.geometrictools.com 
 */
bool CompGeom::query_intersect_tri_tri(const Triangle& t0, const Triangle& t1)
{
  // get edge vectors for triangle0
  Vector3d akE0[3] =
  {
    t0.b - t0.a,
    t0.c - t0.b,
    t0.a - t0.c
  };

  // get normal vector of triangle0
  Vector3d kN0 = Vector3d::cross(akE0[0], akE0[1]);
  normalize_or_zero(kN0);

  // project triangle1 onto normal line of triangle0, test for separation
  double fN0dT0V0 = kN0.dot(t0.a);
  double fMin1, fMax1;
  project_onto_axis(t1,kN0,fMin1,fMax1);
  if (fN0dT0V0 < fMin1 || fN0dT0V0 > fMax1)
    return false;

  // get edge vectors for triangle1
  Vector3d akE1[3] =
  {
    t1.b - t1.a,
    t1.c - t1.b,
    t1.a - t1.c
  };

  // get normal vector of triangle1
  Vector3d kN1 = Vector3d::cross(akE1[0], akE1[1]);
  normalize_or_zero(kN1);

  Vector3d kDir;
  double fMin0, fMax0;
  int i0, i1;

  Vector3d kN0xN1 = Vector3d::cross(kN0, kN1);
  normalize_or_zero(kN0xN1);
  if (kN0xN1.dot(kN0xN1) >= NEAR_ZERO)
  {
    // triangles are not parallel

    // Project triangle0 onto normal line of triangle1, test for
    // separation.
    double fN1dT1V0 = kN1.dot(t1.a);
    project_onto_axis(t0,kN1,fMin0,fMax0);
    if (fN1dT1V0 < fMin0 || fN1dT1V0 > fMax0)
      return false;

    // directions E0[i0]xE1[i1]
    for (i1 = 0; i1 < 3; i1++)
    {
      for (i0 = 0; i0 < 3; i0++)
      {
        kDir = Vector3d::cross(akE0[i0], akE1[i1]);
        normalize_or_zero(kDir);
        project_onto_axis(t0,kDir,fMin0,fMax0);
        project_onto_axis(t1,kDir,fMin1,fMax1);
        if (fMax0 < fMin1 || fMax1 < fMin0)
          return false;
      }
    }
  }
  else  // triangles are parallel (and, in fact, coplanar)
  {
    // directions N0xE0[i0]
    for (i0 = 0; i0 < 3; i0++)
    {
      kDir = Vector3d::cross(kN0, akE0[i0]);
      normalize_or_zero(kDir);
      project_onto_axis(t0,kDir,fMin0,fMax0);
      project_onto_axis(t1,kDir,fMin1,fMax1);
      if (fMax0 < fMin1 || fMax1 < fMin0)
        return false;
    }

    // directions N1xE1[i1]
    for (i1 = 0; i1 < 3; i1++)
    {
      kDir = Vector3d::cross(kN1, akE1[i1]);
      project_onto_axis(t0,kDir,fMin0,fMax0);
      project_onto_axis(t1,kDir,fMin1,fMax1);
      if (fMax0 < fMin1 || fMax1 < fMin0)
        return false;
    }
  }

  return true;
}

/// Normalizes a vector or sets it to zero
void CompGeom::normalize_or_zero(Vector3d& v)
{
  double nrm = v.norm();
  if (nrm < NEAR_ZERO)
    v.set_zero();
  else
    v /= nrm;
}

/// Utility method for intersect_tris_2D()
/**
 * Method adapted from www.geometrictools.com
 */
void CompGeom::clip_convex_polygon_against_line(const Vector2d& rkN, double fC, unsigned& ri, Point2d akV[6])
{
  // input vertices assumed to be ccw

  // test on which side of line the vertices are
  int iPositive = 0, iNegative = 0, iPIndex = -1;
  double afTest[6];
  for (unsigned i=0; i< ri; i++)
  {
    afTest[i] = rkN.dot(akV[i]) - fC;
    if (afTest[i] > (double) 0.0)
    {
      iPositive++;
      if (iPIndex < 0)
        iPIndex = (unsigned) i;
    }
    else if (afTest[i] < (double) 0.0)
     iNegative++;
  }

  if (iPositive > 0)
  {
    if (iNegative > 0)
    {
      // line transversely intersects polygon
      Vector2d akCV[6];
      unsigned iCQuantity = 0;
      int iCur, iPrv;
      double fT;

      if (iPIndex > 0)
      {
        // first clip vertex on line
        iCur = iPIndex;
        iPrv = iCur-1;
        fT = afTest[iCur]/(afTest[iCur] - afTest[iPrv]);
        akCV[iCQuantity++] = akV[iCur]+fT*(akV[iPrv]-akV[iCur]);

        // vertices on positive side of line
        while (iCur < (int) ri && afTest[iCur] > (double) 0.0)
          akCV[iCQuantity++] = akV[iCur++];

        // last clip vertex on line
        if (iCur < (int) ri)
          iPrv = iCur - 1;
        else
        {
          iCur = 0;
          iPrv = ri - 1;
        }
        fT = afTest[iCur]/(afTest[iCur] - afTest[iPrv]);
        akCV[iCQuantity++] = akV[iCur]+fT*(akV[iPrv]-akV[iCur]);
      }
      else // iPIndex is 0
      {
        // vertices on positive side of line
        iCur = 0;
        while (iCur < (int) ri && afTest[iCur] > (double) 0.0)
          akCV[iCQuantity++] = akV[iCur++];

        // last clip vertex on line
        iPrv = iCur-1;
        fT = afTest[iCur]/(afTest[iCur] - afTest[iPrv]);
        akCV[iCQuantity++] = akV[iCur]+fT*(akV[iPrv]-akV[iCur]);

        // skip vertices on negative side
        while (iCur < (int) ri && afTest[iCur] <= (double) 0.0)
          iCur++;

        // first clip vertex on line
        if (iCur < (int) ri)
        {
          iPrv = iCur-1;
          fT = afTest[iCur]/(afTest[iCur] - afTest[iPrv]);
          akCV[iCQuantity++] = akV[iCur]+fT*(akV[iPrv]-akV[iCur]);

          // vertices on positive side of line
          while (iCur < (int) ri && afTest[iCur] > (double) 0.0)
            akCV[iCQuantity++] = akV[iCur++];
        }
        else
        {
          // iCur = 0;
          iPrv = ri - 1;
          fT = afTest[0]/(afTest[0] - afTest[iPrv]);
          akCV[iCQuantity++] = akV[0]+fT*(akV[iPrv]-akV[0]);
        }
      }

      ri = iCQuantity;
      for (unsigned i=0; i< iCQuantity; i++)
        akV[i] = akCV[i];
    }
    // else polygon fully on positive side of line, nothing to do
  }
  else // polygon does not intersect positive side of line, clip all
    ri = 0;
}

