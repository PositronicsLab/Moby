/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/SSL.h>

using namespace Ravelin;
using namespace Moby;
using std::pair;
using boost::dynamic_pointer_cast;

/// Initializes an empty SSL
SSL::SSL()
{
  p1 = p2 = Point3d::zero();
  radius = (double) 0.0;
}

/// Initializes a SSL from the given values
SSL::SSL(const Point3d& p0, const Point3d& p1, double radius)
{
  assert(p1.pose == p2.pose);
  this->p1 = p0;
  this->p2 = p1;
  this->radius = radius;
}

/// Initializes a SSL from an SSL (s) extruded along a direction (v)
SSL::SSL(const SSL& s, const Vector3d& v)
{
  assert(v.pose == s.get_relative_pose());
  geom = s.geom;
  p1 = s.p1; 
  p2 = s.p2;
  radius = s.radius + v.norm();
}

/// Transforms the SSL using the given transform
void SSL::transform(const Transform3d& T, BV* result) const
{
  // get the SSL
  SSL& s = *((SSL*) result);

  // copy this
  s = *this;

  // transform the two points on the line 
  s.p1 = T.transform_point(p1);
  s.p2 = T.transform_point(p2);
}

/// Copies one SSL to another
SSL& SSL::operator=(const SSL& s)
{
  p1 = s.p1;
  p2 = s.p2;
  radius = s.radius;
  geom = s.geom;
  return *this;
}

/// Calculates the velocity-expanded bounding volume
BVPtr SSL::calc_swept_BV(CollisionGeometryPtr g, const SVelocityd& v) const
{
  throw std::runtime_error("SSL::calc_swept_BV() not implemented!");
  return BVPtr();
}

/// Calculates the distance of the sphere-swept line from a point
double SSL::calc_dist(const SSL& o, const Point3d& p)
{
  Point3d cp;
  return std::sqrt(calc_sq_dist(LineSeg3(o.p1, o.p2), p, cp));
}

/// Calculates the squared distance between two line segments (and computes closest points)
/**
 * This code adapted from www.geometrictools.com (thanks Dave Eberly!)
 */
double SSL::calc_sq_dist(const LineSeg3& seg0, const LineSeg3& seg1, Point3d& cp0, Point3d& cp1)
{
  assert(seg0.first.pose == seg1.first.pose);

  const Point3d& p0 = seg0.first;
  const Point3d& p1 = seg0.second;
  const Point3d& q1 = seg1.first;
  const Point3d& q2 = seg1.second;

  // precompute some necessary things
  const Point3d center0 = (p0 + p1) * (double) 0.5;
  const Point3d center1 = (q1 + q2) * (double) 0.5; 
  const Vector3d dir0_unnorm = p1 - p0;
  const Vector3d dir1_unnorm = q2 - q1;
  const double len0 = dir0_unnorm.norm();
  const double len1 = dir1_unnorm.norm();
  assert(len0 > NEAR_ZERO && len1 > NEAR_ZERO);
  const Vector3d dir0 = dir0_unnorm / len0;
  const Vector3d dir1 = dir1_unnorm / len1;
  const double extent0 = len0 * (double) 0.5;
  const double extent1 = len1 * (double) 0.5;

  // now onto Eberly's procedure...
  const Vector3d diff = center0 - center1;
  const double a01 = -dir0.dot(dir1);
  const double b0 = diff.dot(dir0);
  const double b1 = -diff.dot(dir1);
  const double c = diff.norm_sq();
  const double det = std::fabs((double) 1.0 - a01*a01);
  double s0, s1, sq_dist;

  if (det >= NEAR_ZERO)
  {
    // segments not parallel
    s0 = a01*b1 - b0;
    s1 = a01*b0 - b1;
    double extdet0 = extent0*det;
    double extdet1 = extent1*det;

    if (s0 >= -extdet0)
    {
      if (s0 <= extdet0)
      {
        if (s1 >= -extdet1)
        {
          if (s1 <= extdet1)  // region 0 (interior)
          {
            // minimum at interior points of segments
            double invdet = ((double) 1.0) / det;
            s0 *= invdet;
            s1 *= invdet;
            sq_dist = s0*(s0 + a01*s1 + ((double) 2.0)*b0) + s1*(a01*s0 + s1 + ((double) 2.0) * b1) + c;
          }  
          else  // region 3 (side)
          {
            s1 = extent1;
            double tmps0 = -(a01*s1 + b0);
            if (tmps0 < -extent0)
            {
              s0 = -extent0;
              sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
            }
            else if (tmps0 <= extent0)
            {
              s0 = tmps0;
              sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
            }
            else
            {
              s0 = extent0;
              sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
            }
          }
        }
        else   // region 7 (side)
        {
          s1 = -extent1;
          double tmps0 = -(a01*s1 + b0);
          if (tmps0 < -extent0)
          {
            s0 = -extent0;
            sq_dist = s0*(s0 - ((double) 2.0) * tmps0) + s1*(s1 + ((double) 2.0) * b1) + c;
          }
          else if (tmps0 <= extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
          }
          else
          {
            s0 = extent0;
            sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
          }
        }
      }
      else
      {
        if (s1 >= -extdet1)
        {
          if (s1 <= extdet1)  // region 1 (side)
          {
            s0 = extent0;
            double tmps1 = -(a01*s0 + b1);
            if (tmps1 < -extent1)
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
            else if (tmps1 <= extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0)*b0) + c;
            }
            else
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
          }
          else // region 2 (corner)
          {
            s1 = extent1;
            double tmps0 = -(a01*s1 + b0);
            if (tmps0 < -extent0)
            {
              s0 = -extent0;
              sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
            }
            else if (tmps0 < extent0)
            {
              s0 = tmps0;
              sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
            }
            else
            {
              s0 = extent0;
              double tmps1 = -(a01*s0 + b1);
              if (tmps1 < -extent1)
              {
                s1 = -extent1;
                sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
              }
              else if (tmps1 <= extent1)
              {
                s1 = tmps1;
                sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0)*b0) + c;
              }
              else
              {
                s1 = extent1;
                sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
              }
            }
          }
        }
        else // region 8 (corner)
        {
          s1 = -extent1;
          double tmps0 = -(a01*s1 + b0);
          if (tmps0 < -extent0)
          {
            s0 = -extent0;
            sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
          }
          else if (tmps0 < extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
          }
          else
          {
            s0 = extent0;
            double tmps1 = -(a01*s0 + b1);
            if (tmps1 > extent1)
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
            else if (tmps1 >= -extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0) * b0) + c;
            }
            else
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
          }
        }
      }
    }
    else
    {
      if (s1 >= -extdet1)
      {
        if (s1 <= extdet1)  // region 5 (side)
        {
          s0 = -extent0;
          double tmps1 = -(a01*s0 + b1);
          if (tmps1 < -extent1)
          {
            s1 = -extent1;
            sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
          }
          else if (tmps1 <= extent1)
          {
            s1 = tmps1;
            sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0) * b0) + c;
          }
          else
          {
            s1 = extent1;
            sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
          }
        }
        else // region 4 (corner)
        {
          s1 = extent1;
          double tmps0 = -(a01*s1 + b0);
          if (tmps0 > extent0)
          {
            s0 = extent0;
            sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
          }
          else if (tmps0 >= -extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
          }
          else
          {
            s0 = -extent0;
            double tmps1 = -(a01*s0 + b1);
            if (tmps1 < -extent1)
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
            else if (tmps1 <= extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0)*b0) + c;
            }
            else
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
            }
          }
        }
      }
      else  // region 6 (corner)
      { 
        s1 = -extent1;
        double tmps0 = -(a01*s1 + b0);
        if (tmps0 > extent0)
        {
          s0 = extent0;
          sq_dist = s0*(s0 - ((double) 2.0)*tmps0) + s1*(s1 + ((double) 2.0)*b1) + c;
        }
        else if (tmps0 >= -extent0)
        {
          s0 = tmps0;
          sq_dist = -s0*s0 + s1*(s1 + ((double) 2.0)*b1) + c;
        }
        else
        {
          s0 = -extent0;
          double tmps1 = -(a01*s0 + b1);
          if (tmps1 < -extent1)
          {
            s1 = -extent1;
            sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
          }
          else if (tmps1 < extent1)
          {
            s1 = tmps1;
            sq_dist = -s1*s1 + s0*(s0 + ((double) 2.0)*b0) + c;
          }
          else
          {
            s1 = extent1;
            sq_dist = s1*(s1 - ((double) 2.0)*tmps1) + s0*(s0 + ((double) 2.0)*b0) + c;
          }
        }
      }
    }
  }
  else
  {
    // segemnts are parallel; average b0 term is designed to ensure symmetry
    // of the function.  That is, dist(seg1, seg2) and dist(seg2, seg1) should
    // produce the same number.
    const double e0pe1 = extent0 + extent1;
    const double sgn = (a01 > (double) 0.0 ? (double) -1.0 : (double) 1.0);
    const double b0Avr = ((double) 0.5)*(b0 - sgn*b1);
    double lambda = -b0Avr;
    if (lambda < -e0pe1)
      lambda = -e0pe1;
    else if (lambda > e0pe1)
      lambda = e0pe1;
    
    s1 = -sgn*lambda*extent1/e0pe1;
    s0 = lambda + sgn*s1;
    sq_dist = lambda*(lambda + ((double) 2.0)*b0Avr) + c;
  }

  cp0 = center0 + dir0*s0;
  cp1 = center1 + dir1*s1;

  // account for numerical round-off errors
  if (sq_dist < (double) 0.0)
    sq_dist = (double) 0.0;
  return sq_dist;
}

/// Calculates the squared distance between a line segment and point (and computes closest point on line segment)
/**
 * This code adapted from www.geometrictools.com (thanks Dave Eberly!)
 */
double SSL::calc_sq_dist(const LineSeg3& seg, const Point3d& q, Point3d& cp)
{
  const Point3d& p0 = seg.first;
  const Point3d& p1 = seg.second;
  const Vector3d dir_unnorm = p1 - p0;
  const double len = dir_unnorm.norm();
  assert(len > NEAR_ZERO);
  const Vector3d dir = dir_unnorm / len;
  const Point3d center = (p0 + p1) * (double) 0.5;
  const Vector3d diff = q - center;
  const double t = dir.dot(diff);
  const double extent = len * (double) 0.5;
  if (-extent < t)
  {
    if (t < extent)
      cp = center + dir*t;
    else
      cp = p1;
  }
  else
    cp = p0;

  return (cp - q).norm_sq();
}
 
/// Calculates the distance between two sphere-swept lines and also calculates closest points
double SSL::calc_dist(const SSL& a, const SSL& b, Point3d& cpa, Point3d& cpb)
{
  assert(a.get_relative_pose() == b.get_relative_pose());

  // get the two line segments and calculate the distance between them
  LineSeg3 s1(a.p1, a.p2);
  LineSeg3 s2(b.p1, b.p2);
  double dist = std::sqrt(calc_sq_dist(s1, s2, cpa, cpb)) - a.radius - b.radius;

  // determine closest points (on the SSL)
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

  // return the distance
  return dist;
}

/// Calculates the distance between two sphere-swept lines and also calculates closest points (uses a transformation from b to a)
/**
 * \param a the first sphere-swept volume
 * \param b the second sphere-swept volume
 * \param aTb the relative transform from b to a
 * \param cpa the closest point on SSL a (in a's frame)
 */
double SSL::calc_dist(const SSL& a, const SSL& b, const Transform3d& aTb, Point3d& cpa, Point3d& cpb)
{
  assert(aTb.target == a.get_relative_pose() && 
         aTb.source == b.get_relative_pose());

  // create a new SSL (b in a's frame)
  SSL b_a;
  b_a.p1 = aTb.transform_point(b.p1);
  b_a.p2 = aTb.transform_point(b.p2);
  b_a.radius = b.radius;

  return calc_dist(a, b_a, cpa, cpb); 
}

/// Determines whether two SSLs intersect
bool SSL::intersects(const SSL& a, const SSL& b)
{
  assert(a.get_relative_pose() == b.get_relative_pose());

  Point3d tmp;
  double dist = calc_dist(a, b, tmp, tmp);
  return (dist <= (double) 0.0);
}

/// Determines whether two SSLs intersect
bool SSL::intersects(const SSL& a, const SSL& b, const Transform3d& aTb)
{
  assert(aTb.target == a.get_relative_pose() && 
         aTb.source == b.get_relative_pose());

  Point3d tmp;
  double dist = calc_dist(a, b, aTb, tmp, tmp);
  return (dist <= (double) 0.0);
}

/// Determines whether (and when) a line segment intersects with the SSL
bool SSL::intersects(const SSL& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q)
{
  // determine the distance between seg and the underlying line segment; also
  // get the closest points
  Point3d cp_seg, cp_SSL;
  double sq_dist = calc_sq_dist(seg, LineSeg3(a.p1, a.p2), cp_seg, cp_SSL);

  // check for no intersection
  if (std::sqrt(sq_dist) - a.radius > (double) 0.0)
    return false;

  // we have an intersection with the capsule; determine the first point of
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

/// Determines whether a point is outside the SSL
bool SSL::outside(const SSL& a, const Point3d& point, double tol)
{
  double dist = calc_dist(a, point);
  return dist > tol;
}

/// Calculates the size (number of elements in) a SSL tree
unsigned SSL::calc_size() const
{
  unsigned sz = 1;
  BOOST_FOREACH(BVPtr child, children)
    sz += dynamic_pointer_cast<SSL>(child)->calc_size();

  return sz;
}

/// Gets the lower bounds of the SSL
Point3d SSL::get_lower_bounds() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // determine the lower left point
  Point3d ll(get_relative_pose());
  ll[X] = (p1[X] < p2[X]) ? p1[X] : p2[X];
  ll[Y] = (p1[Y] < p2[Y]) ? p1[Y] : p2[Y];
  ll[Z] = (p1[Z] < p2[Z]) ? p1[Z] : p2[Z];

  // move the closest point toward -inf by radius
  const Point3d ones((double) 1.0, (double) 1.0, (double) 1.0, get_relative_pose());
  return ll - ones*radius; 
}

/// Gets the upper bounds of the SSL
Point3d SSL::get_upper_bounds() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // determine the upper right point
  Point3d ur(get_relative_pose());
  ur[X] = (p1[X] > p2[X]) ? p1[X] : p2[X];
  ur[Y] = (p1[Y] > p2[Y]) ? p1[Y] : p2[Y];
  ur[Z] = (p1[Z] > p2[Z]) ? p1[Z] : p2[Z];

  // move the closest point toward inf by radius
  const Point3d ones((double) 1.0, (double) 1.0, (double) 1.0, get_relative_pose());
  return ur + ones*radius;
}

/// Outputs the SSL to VRML (not yet implemented)
std::ostream& SSL::to_vrml(std::ostream& out, const Pose3d& T) const
{
  throw std::runtime_error("SSL::to_vrml() not implemented!");
  return out;
}

