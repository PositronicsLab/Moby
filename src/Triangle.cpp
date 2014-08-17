/****************************************************************************
 * Copyright 2XX5 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <Moby/Plane.h>
#include <Moby/InvalidIndexException.h>
#include <Moby/Constants.h>
#include <Moby/DegenerateTriangleException.h>
#include <Moby/Triangle.h>

using std::pair;
using namespace Ravelin;
using namespace Moby;

/// Constructs the triangle, making copies of the vertices
Triangle::Triangle(const Point3d& a, const Point3d& b, const Point3d& c)
{
  this->a = a;
  this->b = b;
  this->c = c;
}

/// Overloads the assignment operator
void Triangle::operator=(const Triangle& t)
{
  // copy the references
  this->a = t.a;
  this->b = t.b;
  this->c = t.c;
}

/// Destroys the triangle
Triangle::~Triangle()
{
}

/// Calculates the signed distance of a point from the triangle's plane
double Triangle::calc_signed_dist(const Point3d& p) const
{
  Vector3d normal = calc_normal();
  double offset = calc_offset(normal);
  Plane plane(normal, offset);
  return plane.calc_signed_distance(p);
}

/// Transforms the given triangle using the specified pose 
Triangle Triangle::transform(const Triangle& t, const Transform3d& m)
{
  Point3d a = m.transform_point(t.a);
  Point3d b = m.transform_point(t.b);
  Point3d c = m.transform_point(t.c);
  
  return Triangle(a, b, c);
}

/// Gets the specified vertex of the triangle
/**
 * \param i a zero-index in the range [0,2]
 */
const Point3d& Triangle::get_vertex(unsigned i) const
{
  switch (i)
  {
    case 0:
      return a;
    case 1:
      return b;
    case 2:
      return c;
    default:
      throw InvalidIndexException();
  }
 
  // make compiler happy
  return a; 
}

/// Determines what feature of the triangle (vertex, edge, face) a point is on
/**
 * \note uses relative tolerance
 */
Triangle::FeatureType Triangle::determine_feature(const Point3d& point) const
{
  // get the point in barycentric coordinates
  double s, t;
  determine_barycentric_coords(point, s, t);

  // get the feature
  return determine_feature(s, t);
}

/// Determines what feature of the triangle (vertex, edge, face) a point in barycentric coordinates is on
/**
 * \note uses relative tolerance
 */
Triangle::FeatureType Triangle::determine_feature(double s, double t) const
{
  const unsigned X = 0, Y= 1, Z = 2;

  // make sure we have a correct feature
  if (s < -NEAR_ZERO || t < -NEAR_ZERO || s+t > 1.0+NEAR_ZERO)
    return eNone;

  // get the largest value for the first vertex
  double lg_a = std::max(std::max(std::fabs(a[X]), std::fabs(a[Y])), std::fabs(a[Z]));

  // determine the tolerance for the first vertex
  double s_tol = (lg_a > 1.0) ? NEAR_ZERO : NEAR_ZERO/lg_a;

  // see whether on vertex a
  if (std::fabs(s - 1.0) < s_tol)
    return eVertexA;

  // get the largest value for the second vertex
  double lg_b = std::max(std::max(std::fabs(b[X]), std::fabs(b[Y])), std::fabs(b[Z]));

  // determine the tolerance for the second vertex
  double t_tol = (lg_b > 1.0) ? NEAR_ZERO : NEAR_ZERO/lg_b;

  // see whether on vertex b
  if (std::fabs(t - 1.0) < t_tol)
    return eVertexB;

  // get the largest value for the third vertex
  double lg_c = std::max(std::max(std::fabs(c[X]), std::fabs(c[Y])), std::fabs(c[Z]));

  // determine the tolerance for the third vertex
  double u_tol = (lg_c > 1.0) ? NEAR_ZERO : NEAR_ZERO/lg_c;

  // see whether on vertex c
  const double u = 1.0 - s - t;
  if (std::fabs(u - 1.0) < u_tol)
    return eVertexC;

  // see whether on edge bc, ac, or ab (that order)
  if (std::fabs(s) < s_tol)
    return eEdgeBC;
  else if (std::fabs(t) < t_tol)
    return eEdgeAC;
  else if (std::fabs(u) < u_tol)
    return eEdgeAB;

  // still here?  within the face 
  return eFace;
} 

/// Calculates the outward facing normal for this triangle (triangle is assumed to be oriented counter-clockwise)
Vector3d Triangle::calc_normal() const
{
  Vector3d cr = Vector3d::cross(b - a, c - b);
  double cr_nrm = cr.norm();
  if (cr_nrm > std::numeric_limits<double>::epsilon())
    return cr/cr_nrm; 
  else
    throw DegenerateTriangleException(); 
}

/// Computes the area of a triangle
double Triangle::calc_area() const
{
  return Vector3d::norm(Vector3d::cross(b - a, c - a))/2.0;
}

/// Computes the AABB for this triangle
std::pair<Point3d, Point3d> Triangle::calc_AABB() const
{
  const unsigned TRI_VERTS = 3, THREE_D = 3;

  // init box corners
  std::pair<Point3d, Point3d> corners;
  for (unsigned i=0; i< THREE_D; i++)
  {
    corners.first[i] = std::numeric_limits<double>::max();
    corners.second[i] = -std::numeric_limits<double>::max();
  }

  // compute AABB
  for (unsigned i=0; i< TRI_VERTS; i++)
  {
    const Point3d& v = get_vertex(i);
    for (unsigned j=0; j< THREE_D; j++)
    {
      corners.first[j] = std::min(corners.first[j], v[j]);
      corners.second[j] = std::max(corners.second[j], v[j]);
    }
  }

  return corners;
}

/// Calculates 2D barycentric coordinates of a triangle
/**
 * \note adapted from http://www.geometrictools.com
 */
void Triangle::determine_barycentric_coords(const Point2d tri[3], const Point2d& v, double& s, double& t)
{
/* OLD CODE: seemed to work fine, but Eberly's appears more robust numerically
  const double r1 = v[X] - tri[2][X];
  const double r2 = v[Y] - tri[2][Y];
  const double x1 = tri[0][X];
  const double x2 = tri[1][X];
  const double x3 = tri[2][X];
  const double y1 = tri[0][Y];
  const double y2 = tri[1][Y];
  const double y3 = tri[2][Y];
  const double denom = 1.0/(-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
  double rs1 = r1*denom;
  double rs2 = r2*denom;
  s = (y2-y3)*rs1 + (x3-x2)*rs2;
  t = (y3-y1)*rs1 + (x1-x3)*rs2;
*/

  // Compute the vectors relative to V2 of the triangle.
  Vector2d diff[3] =
  {
    tri[0] - tri[2],
    tri[1] - tri[2],
    v - tri[2] 
  };

  // If the vertices have large magnitude, the linear system of equations
  // for computing barycentric coordinates can be ill-conditioned.  To avoid
  // this, uniformly scale the triangle edges to be of order 1.  The scaling
  // of all differences does not change the barycentric coordinates.
  double fmax = (double)0.0;
  for (unsigned i = 0; i < 2; i++)
    for (unsigned j = 0; j < 2; j++)
    {
      double val = std::fabs(diff[i][j]);
      if (val > fmax)
        fmax = val;
     }

  // Scale down only large data.
  double inv_max = (double)0.0;
  if (fmax > (double)1.0)
  {
    inv_max = ((double)1.0)/fmax;
    for (unsigned i = 0; i < 3; i++)
      diff[i] *= inv_max;
  }

  double det = diff[0].dot_perp(diff[1]);
  if (std::fabs(det) > NEAR_ZERO)
  {
    double inv_det = ((double)1.0)/det;
    s = diff[2].dot_perp(diff[1])*inv_det;
    t = diff[0].dot_perp(diff[2])*inv_det;
  }
  else
  {
    // The triangle is a sliver.  Determine the longest edge and
    // compute barycentric coordinates with respect to that edge.
    Vector2d kE2 = tri[0] - tri[1];
    if (inv_max != (double) 0.0)
      kE2 *= inv_max;

    double max_len = kE2.norm();
    unsigned imaxidx = 2;
    double len = diff[1].norm();
    if (len > max_len)
    {
      imaxidx = 1;
      max_len = len;
    }
    len = diff[0].norm();
    if (len > max_len)
    {
      imaxidx = 0;
      max_len = len;
    }

    if (max_len > NEAR_ZERO)
    {
      double inv_len = ((double) 1.0)/max_len;
      if (imaxidx == 0)
      {
        // P-V2 = t(V0-V2)
        s = diff[2].dot(diff[0])*inv_len;
        t = (double) 0.0;
      }
      else if (imaxidx == 1)
      {
        // P-V2 = t(V1-V2)
        s = (double) 0.0;
        t = diff[2].dot(diff[1])*inv_len;
      }
      else if (imaxidx == 2)
      {
        // P-V1 = t(V0-V1)
        diff[2] = v - tri[1];
        if (inv_max != (double) 0.0)
           diff[2] *= inv_max;

        s = diff[2].dot(kE2)*inv_len;
        t = (double) 1.0 - s;
      }
      else
      {
        // The triangle is a nearly a point, just return equal weights.
        s = ((double) 1.0)/(double) 3.0;
        t = s;
      }
    }
  }
}

/// Determines the barycentric coordinates corresponding to a point for this triangle
void Triangle::determine_barycentric_coords(const Point3d& p, double& s, double& t) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  const Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
  double d00 = v0.dot(v0);
  double d01 = v0.dot(v1);
  double d11 = v1.dot(v1);
  double d20 = v2.dot(v0);
  double d21 = v2.dot(v1);
  double denom = d00 * d11 - d01 * d01;
  double v = (d11 * d20 - d01 * d21)/denom;
  double w = (d00 * d21 - d01 * d20)/denom;
  double u = 1.0 - v - w;
  s = u;
  t = v;
  assert((a*s + b*t + c*(1.0-s-t) - p).norm() < NEAR_ZERO);

/*
  // determine which coordinate has the minimum variance
  std::pair<Point3d, Point3d> aabb = calc_AABB();
  const double varx = aabb.second[X] - aabb.first[X];
  const double vary = aabb.second[Y] - aabb.first[Y];
  const double varz = aabb.second[Z] - aabb.first[Z];

  // project to 2D
  Point2d tri_2D[3];
  Point2d v_2D;

  if (varx < vary && varx < varz)
  {
    // minimum variance in the x direction; project to yz plane
    tri_2D[0] = Point2d(a[Y], a[Z]);
    tri_2D[1] = Point2d(b[Y], b[Z]);
    tri_2D[2] = Point2d(c[Y], c[Z]);
    v_2D = Point2d(v[Y], v[Z]);
  }
  else if (vary < varx && vary < varz)
  {
    // minimum variance in the y direction; project to xz plane
    tri_2D[0] = Point2d(a[X], a[Z]);
    tri_2D[1] = Point2d(b[X], b[Z]);
    tri_2D[2] = Point2d(c[X], c[Z]);
    v_2D = Point2d(v[X], v[Z]);
  }
  else
  {
    // minimum variance in the z direction; project to xy plane
    tri_2D[0] = Point2d(a[X], a[Y]);
    tri_2D[1] = Point2d(b[X], b[Y]);
    tri_2D[2] = Point2d(c[X], c[Y]);
    v_2D = Point2d(v[X], v[Y]);
  }
  
  // now determine barycentric coordinates for 2D triangle
  determine_barycentric_coords(tri_2D, v_2D, s, t);
/*
/*  
  // minimum variance is in x direction
  if (varx < vary && varx < varz)
  {
    const double r1 = v[Y] - _c[Y];
    const double r2 = v[Z] - _c[Z];
    const double y1 = _a[Y];
    const double y2 = _b[Y];
    const double y3 = _c[Y];
    const double z1 = _a[Z];
    const double z2 = _b[Z];
    const double z3 = _c[Z];
    const double denom = 1.0/(-y2*z1 + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3);
    double rs1 = r1*denom;
    double rs2 = r2*denom;
    s = (z2-z3)*rs1 + (y3-y2)*rs2;
    t = (z3-z1)*rs1 + (y1-y3)*rs2;
  }
  // minimum variance is in y direction
  else if (vary < varx && vary < varz)
  {
    const double r1 = v[X] - _c[X];
    const double r2 = v[Z] - _c[Z];
    const double x1 = _a[X];
    const double x2 = _b[X];
    const double x3 = _c[X];
    const double z1 = _a[Z];
    const double z2 = _b[Z];
    const double z3 = _c[Z];
    const double denom = 1.0/(-x2*z1 + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3);
    double rs1 = r1*denom;
    double rs2 = r2*denom;
    s = (z2-z3)*rs1 + (x3-x2)*rs2;
    t = (z3-z1)*rs1 + (x1-x3)*rs2;
  }
  // minimum variance is in z direction
  else
  {
    const double r1 = v[X] - _c[X];
    const double r2 = v[Y] - _c[Y];
    const double x1 = _a[X];
    const double x2 = _b[X];
    const double x3 = _c[X];
    const double y1 = _a[Y];
    const double y2 = _b[Y];
    const double y3 = _c[Y];
    const double denom = 1.0/(-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
    double rs1 = r1*denom;
    double rs2 = r2*denom;
    s = (y2-y3)*rs1 + (x3-x2)*rs2;
    t = (y3-y1)*rs1 + (x1-x3)*rs2;
  }
*/
}

/// Sends the triangle to the specified stream using VRML
void Triangle::to_vrml(std::ostream& o, Point3d diffuse_color, bool wireframe) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  o << "Shape" << std::endl;
  o << "{" << std::endl;
  o << "  appearance Appearance { material Material { diffuseColor " << diffuse_color[0] << " " << diffuse_color[1] << " " << diffuse_color[2] << " } }" << std::endl;
  o << "  geometry Indexed";
  if (!wireframe)
    o << "Face";
  else
    o << "Line";
  o << "Set" << std::endl;
  o << "  {" << std::endl;
  if (!wireframe)
    o << "  solid FALSE" << std::endl;
  o << "    coord Coordinate { point [ ";
  o << a[X] << " " << a[Y] << " " << a[Z] << ", ";
  o << b[X] << " " << b[Y] << " " << b[Z] << ", ";
  o << c[X] << " " << c[Y] << " " << c[Z] << " ]}" << std::endl;
  o << "    coordIndex [ ";
  if (!wireframe)
    o << " 0 1 2 -1 ] } }" << std::endl;
  else
    o << " 0 1 2 0 -1 ] } }" << std::endl;
}

/// Determines the distance between a line and a line segment
double Triangle::calc_sq_dist(const Point3d& origin, const Vector3d& dir, const LineSeg3& seg, Point3d& cp_line, Point3d& cp_seg, double& t_line, double& t_seg)
{
  // determine the origins of the segment
  Point3d seg_origin = (seg.first + seg.second)*0.5;

  // determine the directions of the segment
  Point3d seg_dir = seg.second - seg.first;
  double seg_dir_len = seg_dir.norm();
  seg_dir /= seg_dir_len;

  // determine the extents of the segment
  double seg_extent = seg_dir_len * 0.5;

  Vector3d kDiff = origin - seg_origin;
  double fA01 = -dir.dot(seg_dir);
  double fB0 = kDiff.dot(dir);
  double fC = kDiff.norm_sq();
  double fDet = std::fabs(1.0 - fA01*fA01);
  double fB1, fS0, fS1, fSqrDist, fExtDet;

  if (fDet >= NEAR_ZERO)
  {
    // The line and segment are not parallel.
    fB1 = -kDiff.dot(seg_dir);
    fS1 = fA01*fB0-fB1;
    fExtDet = seg_extent*fDet;

    if (fS1 >= -fExtDet)
    {
      if (fS1 <= fExtDet)
      {
        // Two interior points are closest, one on the line and one
        // on the segment.
        double fInvDet = (1.0)/fDet;
        fS0 = (fA01*fB1-fB0)*fInvDet;
        fS1 *= fInvDet;
        fSqrDist = fS0*(fS0+fA01*fS1+(2.0)*fB0)+fS1*(fA01*fS0+fS1+(2.0)*fB1)+fC;
      }
      else
      {
        // The end point e1 of the segment and an interior point of
        // the line are closest.
        fS1 = seg_extent;
        fS0 = -(fA01*fS1+fB0);
        fSqrDist = -fS0*fS0+fS1*(fS1+(2.0)*fB1)+fC;
      }
    }
    else
    {
      // The end point e0 of the segment and an interior point of the
      // line are closest.
      fS1 = -seg_extent;
      fS0 = -(fA01*fS1+fB0);
      fSqrDist = -fS0*fS0+fS1*(fS1+(2.0)*fB1)+fC;
    }
  }
  else
  {
    // The line and segment are parallel.  Choose the closest pair so that
    // one point is at segment origin.
    fS1 = 0.0;
    fS0 = -fB0;
    fSqrDist = fB0*fS0+fC;
  }

  cp_line = origin + fS0*dir;
  cp_seg = seg_origin + fS1*seg_dir;
  t_line = fS0;
  t_seg = fS1;
  return std::fabs(fSqrDist);
}

/// Determines the distance between a triangle and a line
/**
 * \param tri the query triangle
 * \param origin the origin of the line
 * \param dir the direction of the line 
 * \param cp_line the closest point on the line on return
 * \param cp_tri the closest point on the triangle on return
 * \return the Euclidean distance between this triangle and the point
 * \note code adapted from www.geometrictools.com 
 */
double Triangle::calc_sq_dist(const Triangle& tri, const Point3d& origin, const Vector3d& dir, Point3d& cp_line, Point3d& cp_tri, double& t_line)
{
  double t, t_seg;
  Point3d tmp1, tmp2;

  // Test if line intersects triangle.  If so, the squared distance is zero.
  Vector3d edge0 = tri.b - tri.a;
  Vector3d edge1 = tri.c - tri.a;
  Vector3d normal = tri.calc_normal();

  // get the direction and origin of the line segment
  double fNdD = normal.dot(dir);
  if (std::fabs(fNdD) > NEAR_ZERO)
  {
    // The line and triangle are not parallel, so the line intersects
    // the plane of the triangle.
    Vector3d kDiff = origin - tri.a;
    Vector3d kU, kV;
    Vector3d::determine_orthonormal_basis(dir,kU,kV);
    double fUdE0 = kU.dot(edge0);
    double fUdE1 = kU.dot(edge1);
    double fUdDiff = kU.dot(kDiff);
    double fVdE0 = kV.dot(edge0);
    double fVdE1 = kV.dot(edge1);
    double fVdDiff = kV.dot(kDiff);
    double fInvDet = (1.0)/(fUdE0*fVdE1 - fUdE1*fVdE0);

    // Barycentric coordinates for the point of intersection.
    double fB1 = (fVdE1*fUdDiff - fUdE1*fVdDiff)*fInvDet;
    double fB2 = (fUdE0*fVdDiff - fVdE0*fUdDiff)*fInvDet;
    double fB0 = 1.0 - fB1 - fB2;

    if (fB0 >= 0.0 && fB1 >= 0.0 && fB2 >= 0.0)
    {
      // Line parameter for the point of intersection.
      double fDdE0 = dir.dot(edge0);
      double fDdE1 = dir.dot(edge1);
      double fDdDiff = dir.dot(kDiff);
      t_line = fB1*fDdE0 + fB2*fDdE1 - fDdDiff;

      // The intersection point is inside or on the triangle.
      cp_line = origin + t*dir;
      cp_tri = tri.a + fB1*edge0 + fB2*edge1;
      return 0.0;
    }
  }

  // Either (1) the line is not parallel to the triangle and the point of
  // intersection of the line and the plane of the triangle is outside the
  // triangle or (2) the line and triangle are parallel.  Regardless, the
  // closest point on the triangle is on an edge of the triangle.  Compare
  // the line to all three edges of the triangle.
  double fSqrDist = std::numeric_limits<double>::max();
  for (unsigned i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
  {
    LineSeg3 edge(tri.get_vertex(i0), tri.get_vertex(i1));
    double fSqrDistTmp = calc_sq_dist(origin, dir, edge, tmp1, tmp2, t, t_seg);
    if (fSqrDistTmp < fSqrDist)
    {
      cp_line = tmp1;
      cp_tri = tmp2;
      fSqrDist = fSqrDistTmp;
      t_line = t;
    }
  }

  return fSqrDist;
}

/// Computes the closest distance between a triangle and a line segment
double Triangle::calc_sq_dist(const Triangle& tri, const LineSeg3& seg, Point3d& cp_tri, Point3d& cp_seg)
{
  // compute the segment origin, direction, and extent
  Point3d origin = (seg.first + seg.second) * 0.5;
  Vector3d dir = seg.second - seg.first;
  double dir_len = dir.norm();
  double extent = dir_len * 0.5;
  dir /= dir_len;

  // get the squared distance between a line containing the segment and
  // the triangle
  double t_line;
  double fSqrDist = calc_sq_dist(tri, origin, dir, cp_seg, cp_tri, t_line);

  if (t_line >= -extent)
  {
    if (t_line <= extent)
      return fSqrDist;
    else
      cp_seg = seg.second;
  }
  else
    cp_seg = seg.first;

  // still here?  need to compute the closest point to the triangle
  return calc_sq_dist(tri, cp_seg, cp_tri);
}

/// Determines the squared distance between two triangles
/**
 * \param point the query point
 * \param closest_point the closest point on the triangle is returned here
 * \return the squared Euclidean distance between this triangle and the point
 * \note code adapted from www.geometrictools.com 
 */
double Triangle::calc_sq_dist(const Triangle& t1, const Triangle& t2, Point3d& cp1, Point3d& cp2)
{
  Point3d tmp1, tmp2;

  // compare edges of t1 to the interior of t2
  double fSqrDist = std::numeric_limits<double>::max(), fSqrDistTmp;
  LineSeg3 seg;
  for (unsigned i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
  {
    // compute the distance between the triangle and the line segment
    seg.first = t1.get_vertex(i0);
    seg.second = t1.get_vertex(i1);
    fSqrDistTmp = calc_sq_dist(t2, seg, tmp2, tmp1);
    if (fSqrDistTmp < fSqrDist)
    {
      cp1 = tmp2;
      cp2 = tmp1;
      fSqrDist = fSqrDistTmp;
    }
  }

  // compare edges of t2 to the interior of t1
  for (unsigned i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
  {
    seg.first = t2.get_vertex(i0);
    seg.second = t2.get_vertex(i1);
    fSqrDistTmp = calc_sq_dist(t1, seg, tmp1, tmp2);
    if (fSqrDistTmp < fSqrDist)
    {
      cp1 = tmp1;
      cp2 = tmp2;
      fSqrDist = fSqrDistTmp;
    }
  }

  return fSqrDist;
}

/// Determines the distance between this triangle and a given point
/**
 * \param point the query point
 * \param s the Barycentric coordinate s s.t. a*s + b*t + c*(1-s-t) = p
 * \param t the Barycentric coordinate t 
 * \return the Euclidean distance between this triangle and the point
 */
double Triangle::calc_sq_dist(const Triangle& tri, const Point3d& point, double& s, double& t)
{
  // compute needed quantities
  Vector3d k_diff = tri.a - point;
  Vector3d E0 = tri.b - tri.a;
  Vector3d E1 = tri.c - tri.a;
  double A00 = E0.norm_sq();
  double A01 = Vector3d::dot(E0, E1);
  double A11 = E1.norm_sq(); 
  double B0 = Vector3d::dot(k_diff, E0);
  double B1 = Vector3d::dot(k_diff, E1);
  double C = k_diff.norm_sq();
  double det = std::fabs(A00*A11-A01*A01);
  s = A01*B1-A11*B0;
  t = A01*B0-A00*B1;
  double sqr_dist;

  if (s + t <= det)
  {
    if (s < 0.0)
    {
      if (t < 0.0)  
      {
        // region 4
        if (B0 < 0.0)
        {
          t = 0.0;
          if (-B0 >= A00)
          {
            s = 1.0;
            sqr_dist = A00+(2.0*B0)+C;
          }
          else
          {
            s = -B0/A00;
            sqr_dist = (B0*s)+C;
          }
        }
        else
        {
          s = 0.0;
          if (B1 >= 0.0)
          {
            t = 0.0;
            sqr_dist = C;
          }
          else if (-B1 >= A11)
          {
            t = 1.0;
            sqr_dist = A11+(2.0*B1)+C;
          }
          else
          {
            t = -B1/A11;
            sqr_dist = (B1*t)+C;
          }
        }  
      }
      else  
      {
        // region 3
        s = 0.0;
        if (B1 >= 0.0)
        {
          t = 0.0;
          sqr_dist = C;
        }
        else if (-B1 >= A11)
        {
          t = 1.0;
          sqr_dist = A11+(2.0*B1)+C;
        }
        else
        {
          t = -B1/A11;
          sqr_dist = (B1*t)+C;
        }
      }
    }
    else if (t < 0.0)  
    {
      // region 5
      t = 0.0;
      if (B0 >= 0.0)
      {
        s = 0.0;
        sqr_dist = C;
      }
      else if (-B0 >= A00)
      {
        s = 1.0;
        sqr_dist = A00+(2.0*B0)+C;
      }
      else
      {
        s = -B0/A00;
        sqr_dist = (B0*s)+C;
      }
    }
    else  
    {
      // region 0
      double inv_det = 1.0/det;
      s *= inv_det;
      t *= inv_det;
      sqr_dist = s*(A00*s+A01*t+(2.0*B0)) + t*(A01*s+A11*t+(2.0*B1)) + C;
    }
  }
  else
  {
    double tmp0, tmp1, numer, denom;

    if (s < 0.0)  
    {
      // region 2
      tmp0 = A01 + B0;
      tmp1 = A11 + B1;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = A00-2.0*A01+A11;
        if (numer >= denom)
        {
          s = 1.0;
          t = 0.0;
          sqr_dist = A00+(2.0*B0)+C;
        }
        else
        {
          s = numer/denom;
          t = 1.0 - s;
          sqr_dist = s*(A00*s+A01*t+2.0*B0) + t*(A01*s+A11*t+2.0*B1)+C;
        }
      }
      else
      {
        s = 0.0;
        if (tmp1 <= 0.0)
        {
          t = 1.0;
          sqr_dist = A11+2.0*B1+C;
        }
        else if (B1 >= 0.0)
        {
          t = 0.0;
          sqr_dist = C;
        }
        else
        {
          t = -B1/A11;
          sqr_dist = B1*t+C;
        }
      }
    }
    else if (t < 0.0)  
    {
      // region 6
      tmp0 = A01 + B1;
      tmp1 = A00 + B0;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = A00-2.0*A01+A11;
        if (numer >= denom)
        {
          t = 1.0;
          s = 0.0;
          sqr_dist = A11+2.0*B1+C;
        }
        else
        {
          t = numer/denom;
          s = 1.0 - t;
          sqr_dist = s*(A00*s+A01*t+2.0*B0) + t*(A01*s+A11*t+2.0*B1)+C;
        }
      }
      else
      {
        t = 0.0;
        if (tmp1 <= 0.0)
        {
          s = 1.0;
          sqr_dist = A00+2.0*B0+C;
        }
        else if (B0 >= 0.0)
        {
          s = 0.0;
          sqr_dist = C;
        }
        else
        {
          s = -B0/A00;
          sqr_dist = B0*s+C;
        }
      }
    }
    else
    {
      // region 1
      numer = A11 + B1 - A01 - B0;
      if (numer <= 0.0)
      {
        s = 0.0;
        t = 1.0;
        sqr_dist = A11+2.0*B1+C;
      }
      else
      {
        denom = A00-2.0*A01+A11;
        if (numer >= denom)
        {
          s = 1.0;
          t = 0.0;
          sqr_dist = A00+2.0*B0+C;
        }
        else
        {
          s = numer/denom;
          t = 1.0 - s;
          sqr_dist = s*(A00*s+A01*t+2.0*B0) + t*(A01*s+A11*t+2.0*B1)+C;
        }
      }
    }
  }

  // convert to standard notation 
  const double SPRIME = 1-s-t;
  const double TPRIME = s;
  s = SPRIME;
  t = TPRIME;
  return sqr_dist;
}

/// Determines the signed distance between this triangle and a given point
/**
 * \param point the query point
 * \param closest_point the closest point on the triangle is returned here
 * \return the Euclidean distance between this triangle and the point
 */
double Triangle::calc_signed_dist(const Point3d& point, Point3d& cp) const
{
  return Triangle::calc_signed_dist(*this, point, cp);
}

/// Determines the signed distance between this triangle and a given point
/**
 * \param point the query point
 * \param closest_point the closest point on the triangle is returned here
 * \return the Euclidean distance between this triangle and the point
 */
double Triangle::calc_signed_dist(const Triangle& tri, const Point3d& point, Point3d& cp) 
{
  double sq_dist = calc_sq_dist(tri, point, cp);
  if (tri.calc_signed_dist(point) >= 0.0)
    return std::sqrt(sq_dist);
  else
    return -std::sqrt(sq_dist);
}

/// Determines the distance between this triangle and a given point
/**
 * \param point the query point
 * \param closest_point the closest point on the triangle is returned here
 * \return the Euclidean distance between this triangle and the point
 */
double Triangle::calc_sq_dist(const Triangle& tri, const Point3d& point, Point3d& cp)
{
  // get the closest point in barycentric coordinates
  double s, t;
  double sqr_dist = calc_sq_dist(tri, point, s, t);

  // compute the closest point
  cp = tri.a*s + tri.b*t + tri.c*(1-s-t);

  // return the distance
  return sqr_dist;
}

/// Sends the triangle to the specified stream
std::ostream& Moby::operator<<(std::ostream& o, const Triangle& t)
{
  o << t.a << "; " << t.b << "; " << t.c;
  return o;
}

