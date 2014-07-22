/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/SSR.h>

using namespace Ravelin;
using std::pair;
using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using namespace Moby;

/// Initializes an empty SSR
SSR::SSR()
{
  center = Point3d::zero();
  R = IDENTITY_3x3;
  l = Vector2d((double) 0.0, (double) 0.0);
  radius = (double) 0.0;
}

/// Initializes a SSR from the given values
SSR::SSR(const Point3d& center, const Matrix3d& R, const Vector2d& l, double radius)
{
  this->center = center;
  this->R = R;
  this->l = l;
  this->radius = radius;
}

/// Initializes a SSR from an SSR (s) extruded along a direction (v)
SSR::SSR(const SSR& s, const Vector3d& v)
{
  assert(s.get_relative_pose() == v.pose);
  geom = s.geom;
  center = s.center;
  R = s.R;
  l = s.l;
  radius = s.radius + v.norm();
}

/// Transforms the SSR using the given transform
void SSR::transform(const Transform3d& T, BV* result) const
{
  // get the SSR
  SSR& s = *((SSR*) result);

  // copy this
  s = *this;

  // transform the center
  s.center = T.transform_point(center);

  // transform the orientation
  s.R = Matrix3d(T.q) * R;
}

/// Copies one SSR to another
SSR& SSR::operator=(const SSR& s)
{
  center = s.center;
  geom = s.geom;
  R = s.R;
  l = s.l;
  radius = s.radius;
  return *this;
}

/// Calculates the velocity-expanded bounding volume
BVPtr SSR::calc_swept_BV(CollisionGeometryPtr g, const SVelocityd& v) const
{
  throw std::runtime_error("SSR::calc_swept_BV() not implemented!");
  return BVPtr();
}

/// Calculates the distance of the sphere-swept rectangle from a line segment
double SSR::calc_dist(const SSR& s, const LineSeg3& seg)
{
  const unsigned Y = 1, Z = 2;

  // setup necessary variables to do distance from line seg to rectangle calc
  const Point3d& rect_cent = s.center;
  Vector3d axis1(s.R.get_column(Y), s.get_relative_pose());
  Vector3d axis2(s.R.get_column(Z), s.get_relative_pose());
  const Vector2d& lengths = s.l;
  Point3d cp_rect, cp_seg;

  // do the calculation
  double sq_dist = calc_sq_dist(seg, rect_cent, axis1, axis2, lengths, cp_seg, cp_rect);

  return std::sqrt(sq_dist) - s.radius;
}

/// Calculates the distance of the sphere-swept rectangle from a point
double SSR::calc_dist(const SSR& o, const Point3d& p)
{
  const unsigned Y = 1, Z = 2;

  // setup necessary variables to do distance from point to rectangle calc
  const Point3d& rect_cent = o.center;
  Vector3d axis1(o.R.get_column(Y), o.get_relative_pose());
  Vector3d axis2(o.R.get_column(Z), o.get_relative_pose());
  const Vector2d& lengths = o.l;
  Point3d cp_rect;

  // do the calculation
  double sq_dist = calc_sq_dist(p, rect_cent, axis1, axis2, lengths, cp_rect);

  return std::sqrt(sq_dist) - o.radius;
}

/// Computes the squared distance and closest point between a rectangle in 3D and a point
double SSR::calc_sq_dist(const Point3d& p, const Point3d& rect_center, const Vector3d& axis1, const Vector3d& axis2, const Vector2d& lengths, Point3d& cp_rect)
{
  assert(p.pose == rect_center.pose);
  assert(axis1.pose == rect_center.pose);
  assert(axis2.pose == rect_center.pose);

  // setup extents
  double extents[2] = { lengths[0]*0.5, lengths[1]*0.5 };

  // ** determine the distance to the rectangle **
  // following code adapted from www.geometrictools.com - thanks Dave Eberly!
  // setup the rectangle axes
  Vector3d diff = rect_center - p;
  double b0 = diff.dot(axis1);
  double b1 = diff.dot(axis2);
  double s0 = -b0, s1 = -b1;
  double sq_dist = diff.norm_sq();
  if (s0 < -extents[0])
  {
    s0 = -extents[0];
  }
  else if (s0 > extents[0])
  {
    s0 = extents[0];
  }
  sq_dist += s0*(s0 + ((double) 2.0)*b0);

  if (s1 < -extents[1])
  {
    s1 = -extents[1];
  }
  else if (s1 > extents[1])
  {
    s1 = extents[1];
  }
  sq_dist += s1*(s1 + ((double) 2.0)*b1);

  // Account for numerical round-off error.
  if (sq_dist < (double) 0.0)
  {
    sq_dist = (double) 0.0;
  }

  // closest point calculation is below
  cp_rect = rect_center + axis1*s0 + axis2*s1;

  return sq_dist;
}

/// Calculates the squared distance between a line and a line segment
double SSR::calc_sq_dist(const Point3d& origin, const Vector3d& dir, const LineSeg3& seg, Point3d& cp_line, Point3d& cp_seg, double& line_param)
{
  assert(origin.pose == dir.pose);
  assert(origin.pose == seg.first.pose);
  assert(origin.pose == seg.second.pose);

  // setup the line segment center and direction
  Point3d center = (seg.first + seg.second)*0.5;
  
  Vector3d seg_dir_unnormalized = seg.second - seg.first;
  double seg_len = seg_dir_unnormalized.norm();
  Vector3d seg_dir = seg_dir_unnormalized / seg_len;
  double seg_ext = seg_len * 0.5;

  // start calculating
  Vector3d diff = origin - center;

  double a01 = -dir.dot(seg_dir);
  double b0 = diff.dot(dir);
  double c = diff.norm_sq();
  double det = std::fabs((double) 1.0 - a01*a01);
  double b1, s0, s1, sq_dist, ext_det;

  if (det >= NEAR_ZERO)
  {
    // The line and segment are not parallel.
    b1 = -diff.dot(seg_dir);
    s1 = a01*b0 - b1;
    ext_det = seg_ext*det;

    if (s1 >= -ext_det)
    {
      if (s1 <= ext_det)
      {
        // Two interior points are closest, one on the line and one
        // on the segment.
        double inv_det = ((double) 1.0)/det;
        s0 = (a01*b1 - b0)*inv_det;
        s1 *= inv_det;
        sq_dist = s0*(s0 + a01*s1 + ((double) 2.0)*b0) +
                   s1*(a01*s0 + s1 + ((double) 2.0)*b1) + c;
      }
      else
      {
        // The endpoint e1 of the segment and an interior point of
        // the line are closest.
        s1 = seg_ext;
        s0 = -(a01*s1 + b0);
        sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
      }
    }
    else
    {
      // The end point e0 of the segment and an interior point of the
      // line are closest.
      s1 = -seg_ext;
      s0 = -(a01*s1 + b0);
      sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
    }
  }
  else
  {
    // The line and segment are parallel.  Choose the closest pair so that
    // one point is at segment center.
    s1 = (double) 0.0;
    s0 = -b0;
    sq_dist = b0*s0 + c;
  }

  // setup closest points
  cp_line = origin + dir*s0;
  cp_seg = center + seg_dir*s1;
  line_param = s0;
//  seg_param = s1;

  // Account for numerical round-off errors.
  if (sq_dist < (double) 0.0)
    sq_dist = (double) 0.0;
  return sq_dist;
}

/// Calculates the square distance between a line and a rectangle in 3D
double SSR::calc_sq_dist(const Point3d& origin, const Vector3d& dir, const Point3d& rect_center, const Vector3d& axis1, const Vector3d& axis2, const Vector2d& lengths, Point3d& cp_line, Point3d& cp_rect, double& line_param)
{
  const unsigned N_RECT_EDGES = 4;

  assert(origin.pose == dir.pose);
  assert(origin.pose == rect_center.pose);
  assert(axis1.pose == rect_center.pose);
  assert(axis2.pose == rect_center.pose);

  // setup the extents
  double extents[2] = { lengths[0]*0.5, lengths[1]*0.5 };

  // ** Routine for calculating distance between line and rectangle courtesy
  //    www.geometricalgorithms.com (thanks Dave Eberly!)
  // test whether line intersects rectangle; if so, the squared distance is zero
  Vector3d N = Vector3d::cross(axis1, axis2);
  double NdD = N.dot(dir);
  if (std::fabs(NdD) > NEAR_ZERO)
  {
    // the line and rectangle are not parallel, so the line intersects the
    // plane of the rectangle
    Vector3d diff = origin - rect_center;
    Vector3d U, V;
    Vector3d::determine_orthonormal_basis(dir, U, V);
    double UdD0 = U.dot(axis1);
    double UdD1 = U.dot(axis2);
    double UdPmC = U.dot(diff);
    double VdD0 = V.dot(axis1);
    double VdD1 = V.dot(axis2);
    double VdPmC = V.dot(diff);
    double inv_det = (double) 1.0/(UdD0*VdD1 - UdD1*VdD0);

    // rectangle coordinates for the point of intersection
    double s0 = (VdD1*UdPmC - UdD1*VdPmC)*inv_det;
    double s1 = (UdD0*VdPmC - VdD0*UdPmC)*inv_det;

    if (std::fabs(s0) <= extents[0] && std::fabs(s1) <= extents[1])
    {
      // line parameter for the point of intersection
      double DdD0 = dir.dot(axis1);
      double DdD1 = dir.dot(axis2);
      double DdDiff = dir.dot(diff);
      line_param = s0*DdD0 + s1*DdD1 - DdDiff;

      // the intersection point is inside or on the rectangle
      cp_line = origin + dir*line_param;
      cp_rect = rect_center + axis1*s0 + axis2*s1;

      return 0.0;
    }
  }

  // either (1) the line is not parallel to the rectangle and the point of
  // intersection of the line and the plane of the rectangle is outside the
  // rectangle or (2) the line and rectangle are parallel.  Regardless, the
  // closest point on the rectangle is on an edge of the rectangle.  Compare
  // the line to all four edges of the rectangle
  double min_sq_dist = std::numeric_limits<double>::max();

  // ** here we deviate from Eberly's code
  // form the four edges of the rectangle
  LineSeg3 rect_edges[N_RECT_EDGES];
  rect_edges[0].first = rect_center - axis1*extents[0] - axis2*extents[1]; 
  rect_edges[0].second = rect_center + axis1*extents[0] - axis2*extents[1];
  rect_edges[1].first = rect_edges[0].second;
  rect_edges[1].second = rect_center + axis1*extents[0] + axis2*extents[1];
  rect_edges[2].first = rect_edges[1].second;
  rect_edges[2].second = rect_center - axis1*extents[0] + axis2*extents[1];
  rect_edges[3].first = rect_edges[2].second;
  rect_edges[3].second = rect_edges[0].first;

  for (unsigned i=0; i< N_RECT_EDGES; i++)  
  {
    double lparam;
    Point3d closest_point_line, closest_point_seg;
    double sq_dist = calc_sq_dist(origin, dir, rect_edges[i], closest_point_line, closest_point_seg, lparam);
    if (sq_dist < min_sq_dist)
    {
      min_sq_dist = sq_dist;
      cp_line = closest_point_line;
      cp_rect = closest_point_seg;
      line_param = lparam;
    }
  }

  return min_sq_dist;
}

/// Calculates the squared distance between a line segment and a rectangle in 3D
double SSR::calc_sq_dist(const LineSeg3& seg, const Point3d& rect_center, const Vector3d& axis1, const Vector3d& axis2, const Vector2d& lengths, Point3d& cp_seg, Point3d& cp_rect)
{
  assert(rect_center.pose == seg.first.pose);
  assert(rect_center.pose == seg.second.pose);
  assert(rect_center.pose == axis1.pose);
  assert(rect_center.pose == axis2.pose);

  // calculate the line segment extents
  Vector3d seg_dir_unnormalized = seg.second - seg.first;
  double seg_len = seg_dir_unnormalized.norm();
  double seg_ext = seg_len * 0.5;

  // form a line from the segment
  Vector3d line_dir = seg_dir_unnormalized / seg_len;
  
  // determine the squared distance and closest point on the line to the rect
  double line_param;
  Point3d cp_line;
  double sq_dist_line = calc_sq_dist(seg.first, line_dir, rect_center, axis1, axis2, lengths, cp_line, cp_rect, line_param);
  if (line_param >= -seg_ext)
  {
    if (line_param <= seg_ext)
    {
      cp_seg = cp_line;
      return sq_dist_line;
    }
    else
    {
      // closest point to the line is at the segment endpoint; recompute
      // distance and closest point to segment endpoint
      cp_seg = seg.second;
      return calc_sq_dist(cp_seg, rect_center, axis1, axis2, lengths, cp_rect); 
    }
  }
  else
  {
    cp_seg = seg.first;
    return calc_sq_dist(cp_seg, rect_center, axis1, axis2, lengths, cp_rect);
  }
}

/// Calculates the squared distance between two rectangles in 3D
double SSR::calc_sq_dist(const Point3d& a_center, const Vector3d& aaxis1, const Vector3d& aaxis2, const Vector2d& alengths, const Point3d& b_center, const Vector3d& baxis1, const Vector3d& baxis2, const Vector2d& blengths, Point3d& cpa, Point3d& cpb)
{
  assert(a_center.pose == aaxis1.pose);
  assert(a_center.pose == aaxis2.pose);
  assert(a_center.pose == b_center.pose);
  assert(a_center.pose == baxis1.pose);
  assert(a_center.pose == baxis2.pose);

  const unsigned N_RECT_EDGES = 4;

  // setup the extents for the four rectangles
  double aextents[2] = { alengths[0]*0.5, alengths[1]*0.5 };
  double bextents[2] = { blengths[0]*0.5, blengths[1]*0.5 };

  // form the four edges of rectangle a
  LineSeg3 rect_a_edges[N_RECT_EDGES];
  rect_a_edges[0].first = a_center - aaxis1*aextents[0] - aaxis2*aextents[1]; 
  rect_a_edges[0].second = a_center + aaxis1*aextents[0] - aaxis2*aextents[1];
  rect_a_edges[1].first = rect_a_edges[0].second;
  rect_a_edges[1].second = a_center + aaxis1*aextents[0] + aaxis2*aextents[1];
  rect_a_edges[2].first = rect_a_edges[1].second;
  rect_a_edges[2].second = a_center - aaxis1*aextents[0] + aaxis2*aextents[1];
  rect_a_edges[3].first = rect_a_edges[2].second;
  rect_a_edges[3].second = rect_a_edges[0].first;

  // form the four edges of rectangle b 
  LineSeg3 rect_b_edges[N_RECT_EDGES];
  rect_b_edges[0].first = b_center - baxis1*bextents[0] - baxis2*bextents[1]; 
  rect_b_edges[0].second = b_center + baxis1*bextents[0] - baxis2*bextents[1];
  rect_b_edges[1].first = rect_b_edges[0].second;
  rect_b_edges[1].second = b_center + baxis1*bextents[0] + baxis2*bextents[1];
  rect_b_edges[2].first = rect_b_edges[1].second;
  rect_b_edges[2].second = b_center - baxis1*bextents[0] + baxis2*bextents[1];
  rect_b_edges[3].first = rect_b_edges[2].second;
  rect_b_edges[3].second = rect_b_edges[0].first;

  // now determine which has the smallest distance
  double min_sq_dist = std::numeric_limits<double>::max();
  Point3d closest_a, closest_b;  

  // first test edges of a against rectangle b
  for (unsigned i=0; i< N_RECT_EDGES; i++)
  {
    double sq_dist = calc_sq_dist(rect_a_edges[i], b_center, baxis1, baxis2, blengths, closest_a, closest_b);
    if (sq_dist < min_sq_dist)
    {
      min_sq_dist = sq_dist;
      cpa = closest_a;
      cpb = closest_b;
    }
  }

  // now test edges of b against rectangle a
  for (unsigned i=0; i< N_RECT_EDGES; i++)
  {
    double sq_dist = calc_sq_dist(rect_b_edges[i], a_center, aaxis1, aaxis2, alengths, closest_b, closest_a);
    if (sq_dist < min_sq_dist)
    {
      min_sq_dist = sq_dist;
      cpa = closest_a;
      cpb = closest_b;
    }
  }

  return min_sq_dist;
}

/// Calculates the distance between two sphere-swept rectangles and also calculates closest points
double SSR::calc_dist(const SSR& a, const SSR& b, Point3d& cpa, Point3d& cpb)
{
  const unsigned Y = 1, Z = 2;

  // verify the two are in the same frame
  shared_ptr<const Pose3d> P = a.get_relative_pose();
  assert(P == b.get_relative_pose());

  // calculate the squared distance and closest points between the two rects
  const Point3d& acenter = a.center;
  const Point3d& bcenter = b.center;
  Vector3d aaxis1(a.R.get_column(Y), P);
  Vector3d aaxis2(a.R.get_column(Z), P);
  Vector3d baxis1(b.R.get_column(Y), P);
  Vector3d baxis2(b.R.get_column(Z), P);
  const Vector2d& alengths = a.l;
  const Vector2d& blengths = b.l;
  double sq_dist = calc_sq_dist(acenter, aaxis1, aaxis2, alengths, 
                              bcenter, baxis1, baxis2, blengths, cpa, cpb);

  // determine direction from cpa to cpb
  Vector3d dir = cpb - cpa;
  double dir_norm = dir.norm();
  if (dir_norm > NEAR_ZERO)
  {
    // normalize the direction
    dir /= dir_norm;

    // move the point on cpa along dir by radius; do the same for cpb (neg)
    cpa += dir*a.radius;
    cpb -= dir*b.radius;
  }

  return std::sqrt(sq_dist) - (a.radius + b.radius);
}

/// Calculates the distance between two sphere-swept rectangles and also calculates closest points (uses a transformation from b to a)
/**
 * \param a the first sphere-swept volume
 * \param b the second sphere-swept volume
 * \param aTb the relative transform from b to a
 * \param cpa the closest point on SSR a (in a's frame)
 */
double SSR::calc_dist(const SSR& a, const SSR& b, const Transform3d& aTb, Point3d& cpa, Point3d& cpb)
{
  // create a new SSR (b in a's frame)
  SSR b_a;
  b_a.center = aTb.transform_point(b.center);
  Matrix3d R = aTb.q;
  b_a.R = R * b.R;
  b_a.l = b.l;
  b_a.radius = b.radius;

  return calc_dist(a, b_a, cpa, cpb); 
}

/// Determines whether two SSRs intersect
bool SSR::intersects(const SSR& a, const SSR& b)
{
  Point3d tmp;
  double dist = calc_dist(a, b, tmp, tmp);
  return (dist <= (double) 0.0);
}

/// Determines whether two SSRs intersect
bool SSR::intersects(const SSR& a, const SSR& b, const Transform3d& T)
{
  Point3d tmp;
  double dist = calc_dist(a, b, T, tmp, tmp);
  return (dist <= (double) 0.0);
}

/// Determines whether (and when) a line segment intersects with the SSR
bool SSR::intersects(const SSR& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q)
{
  const unsigned Y = 1, Z = 2;

  // verify everything is in the proper frame
  shared_ptr<const Pose3d> P = a.get_relative_pose();
  assert(seg.first.pose == seg.second.pose && seg.first.pose == P);

  // get the rectangle from the SSR
  Vector3d axis1(a.R.get_column(Y), P);
  Vector3d axis2(a.R.get_column(Z), P);
  const Point3d& rect_center = P;
  const Vector2d& lengths = a.l;

  // determine the distance between the line segment and the rectangle; also
  // get the closest points
  Point3d cp_seg, cp_rect;
  double sq_dist = calc_sq_dist(seg, rect_center, axis1, axis2, lengths, cp_seg, cp_rect);

  // check for no intersection
  if (std::sqrt(sq_dist) - a.radius > (double) 0.0)
    return false;

  // we have an intersection with the lozenge; determine the first point of
  // intersection
  Vector3d dir_unnormed = seg.first - cp_seg;
  double dir_unnormed_len = dir_unnormed.norm();
  if (dir_unnormed_len < NEAR_ZERO)
  {
    if (tmin > (double) 0.0)
      return false;
    else
    {
      tmin = (double) 0.0;
      q = seg.first;
      return true;
    }
  }
  Vector3d dir = dir_unnormed / dir_unnormed_len;
  q = cp_seg + dir*a.radius;

  // determine the parameter for q
  tmin = (q - seg.first).norm() / (seg.second - seg.first).norm();
  return true;
}

/// Determines whether a point is outside the SSR
bool SSR::outside(const SSR& a, const Point3d& point, double tol)
{
  double dist = calc_dist(a, point);
  return dist > tol;
}

/// Calculates the size (number of elements in) a SSR tree
unsigned SSR::calc_size() const
{
  unsigned sz = 1;
  BOOST_FOREACH(BVPtr child, children)
    sz += dynamic_pointer_cast<SSR>(child)->calc_size();

  return sz;
}

/// Gets the lower bounds of the SSR
Point3d SSR::get_lower_bounds() const
{
  const unsigned Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  const double SQRT3_3 = std::sqrt((double) 3.0)/3.0;

  // get the current pose
  shared_ptr<const Pose3d> P = get_relative_pose();

  // get the axes of the rectangle in 3D
  Vector3d axis1(R.get_column(Y), P);
  Vector3d axis2(R.get_column(Z), P);
  axis1.pose = axis2.pose = P;

  // setup the point at -INF
  Point3d ninf(-INF, -INF, -INF, P);
 
  // find the closest point on the rectangle to -inf
  Point3d cp;
  calc_sq_dist(ninf, center, axis1, axis2, this->l, cp);
  
  // move the closest point toward -inf by radius
  const Vector3d ones_norm(SQRT3_3, SQRT3_3, SQRT3_3, P);
  return cp - ones_norm*radius; 
}

/// Gets the upper bounds of the SSR
Point3d SSR::get_upper_bounds() const
{
  const unsigned Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  const double SQRT3_3 = std::sqrt((double) 3.0)/3.0;

   // get the current pose
  shared_ptr<const Pose3d> P = get_relative_pose();

 // get the axes of the rectangle in 3D
  Vector3d axis1(R.get_column(Y), P);
  Vector3d axis2(R.get_column(Z), P);
  axis1.pose = axis2.pose = P;

  // setup the point at INF
  Vector3d inf(INF, INF, INF, P);
 
  // find the closest point on the rectangle to inf
  Point3d cp;
  calc_sq_dist(inf, center, axis1, axis2, this->l, cp);
  
  // move the closest point toward inf by radius
  const Vector3d ones_norm(SQRT3_3, SQRT3_3, SQRT3_3, P);
  return cp + ones_norm*radius;
}

/// Gets the vertices of the rectangle
void SSR::get_rect_verts(Point3d rect_verts[4]) const
{
  const unsigned Y = 1, Z = 2;

  // setup the extents for the rectangle
  double extents[2] = { l[0]*0.5, l[1]*0.5 };

  // get the relative pose
  shared_ptr<const Pose3d> P = get_relative_pose();

  // get the axes for the rectangle
  Vector3d axis1(R.get_column(Y), P);
  Vector3d axis2(R.get_column(Z), P);

  // setup the vertices
  rect_verts[0] = center - axis1*extents[0] - axis2*extents[1]; 
  rect_verts[1] = center + axis1*extents[0] - axis2*extents[1]; 
  rect_verts[2] = center + axis1*extents[0] + axis2*extents[1]; 
  rect_verts[3] = center - axis1*extents[0] + axis2*extents[1]; 
}

/// Outputs the SSR to VRML (not yet implemented)
std::ostream& SSR::to_vrml(std::ostream& out, const Pose3d& T) const
{
  throw std::runtime_error("SSR::to_vrml() not implemented!");
  return out;
}

