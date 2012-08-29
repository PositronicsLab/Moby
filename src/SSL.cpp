/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/SSL.h>

using namespace Moby;
using boost::dynamic_pointer_cast;

/// Initializes an empty SSL
SSL::SSL()
{
  p1 = p2 = ZEROS_3;
  radius = (Real) 0.0;
}

/// Initializes a SSL from the given values
SSL::SSL(const Vector3& p0, const Vector3& p1, Real radius)
{
  this->p1 = p0;
  this->p2 = p1;
  this->radius = radius;
}

/// Initializes a SSL from an SSL (s) extruded along a direction (v)
SSL::SSL(const SSL& s, const Vector3& v)
{
  p1 = s.p1; 
  p2 = s.p2;
  radius = s.radius + v.norm();
}

/// Copies one SSL to another
void SSL::operator=(const SSL& s)
{
  p1 = s.p1;
  p2 = s.p2;
  radius = s.radius;
}

/// Calculates the velocity-expanded bounding volume
BVPtr SSL::calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const
{
  throw std::runtime_error("SSL::calc_vel_exp_BV() not implemented!");
  return BVPtr();
}

/// Calculates the distance of the sphere-swept line from a point
Real SSL::calc_dist(const SSL& o, const Vector3& p)
{
  Vector3 cp;
  return std::sqrt(calc_sq_dist(LineSeg3(o.p1, o.p2), p, cp));
}

/// Calculates the squared distance between two line segments (and computes closest points)
/**
 * This code adapted from www.geometrictools.com (thanks Dave Eberly!)
 */
Real SSL::calc_sq_dist(const LineSeg3& seg0, const LineSeg3& seg1, Vector3& cp0, Vector3& cp1)
{
  const Vector3& p0 = seg0.first;
  const Vector3& p1 = seg0.second;
  const Vector3& q1 = seg1.first;
  const Vector3& q2 = seg1.second;

  // precompute some necessary things
  const Vector3 center0 = (p0 + p1) * (Real) 0.5;
  const Vector3 center1 = (q1 + q2) * (Real) 0.5; 
  const Vector3 dir0_unnorm = p1 - p0;
  const Vector3 dir1_unnorm = q2 - q1;
  const Real len0 = dir0_unnorm.norm();
  const Real len1 = dir1_unnorm.norm();
  assert(len0 > NEAR_ZERO && len1 > NEAR_ZERO);
  const Vector3 dir0 = dir0_unnorm / len0;
  const Vector3 dir1 = dir1_unnorm / len1;
  const Real extent0 = len0 * (Real) 0.5;
  const Real extent1 = len1 * (Real) 0.5;

  // now onto Eberly's procedure...
  const Vector3 diff = center0 - center1;
  const Real a01 = -dir0.dot(dir1);
  const Real b0 = diff.dot(dir0);
  const Real b1 = -diff.dot(dir1);
  const Real c = diff.norm_sq();
  const Real det = std::fabs((Real) 1.0 - a01*a01);
  Real s0, s1, sq_dist;

  if (det >= NEAR_ZERO)
  {
    // segments not parallel
    s0 = a01*b1 - b0;
    s1 = a01*b0 - b1;
    Real extdet0 = extent0*det;
    Real extdet1 = extent1*det;

    if (s0 >= -extdet0)
    {
      if (s0 <= extdet0)
      {
        if (s1 >= -extdet1)
        {
          if (s1 <= extdet1)  // region 0 (interior)
          {
            // minimum at interior points of segments
            Real invdet = ((Real) 1.0) / det;
            s0 *= invdet;
            s1 *= invdet;
            sq_dist = s0*(s0 + a01*s1 + ((Real) 2.0)*b0) + s1*(a01*s0 + s1 + ((Real) 2.0) * b1) + c;
          }  
          else  // region 3 (side)
          {
            s1 = extent1;
            Real tmps0 = -(a01*s1 + b0);
            if (tmps0 < -extent0)
            {
              s0 = -extent0;
              sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
            }
            else if (tmps0 <= extent0)
            {
              s0 = tmps0;
              sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
            }
            else
            {
              s0 = extent0;
              sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
            }
          }
        }
        else   // region 7 (side)
        {
          s1 = -extent1;
          Real tmps0 = -(a01*s1 + b0);
          if (tmps0 < -extent0)
          {
            s0 = -extent0;
            sq_dist = s0*(s0 - ((Real) 2.0) * tmps0) + s1*(s1 + ((Real) 2.0) * b1) + c;
          }
          else if (tmps0 <= extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
          }
          else
          {
            s0 = extent0;
            sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
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
            Real tmps1 = -(a01*s0 + b1);
            if (tmps1 < -extent1)
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
            else if (tmps1 <= extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
            else
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
          }
          else // region 2 (corner)
          {
            s1 = extent1;
            Real tmps0 = -(a01*s1 + b0);
            if (tmps0 < -extent0)
            {
              s0 = -extent0;
              sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
            }
            else if (tmps0 < extent0)
            {
              s0 = tmps0;
              sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
            }
            else
            {
              s0 = extent0;
              Real tmps1 = -(a01*s0 + b1);
              if (tmps1 < -extent1)
              {
                s1 = -extent1;
                sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
              }
              else if (tmps1 <= extent1)
              {
                s1 = tmps1;
                sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0)*b0) + c;
              }
              else
              {
                s1 = extent1;
                sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
              }
            }
          }
        }
        else // region 8 (corner)
        {
          s1 = -extent1;
          Real tmps0 = -(a01*s1 + b0);
          if (tmps0 < -extent0)
          {
            s0 = -extent0;
            sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
          }
          else if (tmps0 < extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
          }
          else
          {
            s0 = extent0;
            Real tmps1 = -(a01*s0 + b1);
            if (tmps1 > extent1)
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
            else if (tmps1 >= -extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0) * b0) + c;
            }
            else
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
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
          Real tmps1 = -(a01*s0 + b1);
          if (tmps1 < -extent1)
          {
            s1 = -extent1;
            sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
          }
          else if (tmps1 <= extent1)
          {
            s1 = tmps1;
            sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0) * b0) + c;
          }
          else
          {
            s1 = extent1;
            sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
          }
        }
        else // region 4 (corner)
        {
          s1 = extent1;
          Real tmps0 = -(a01*s1 + b0);
          if (tmps0 > extent0)
          {
            s0 = extent0;
            sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
          }
          else if (tmps0 >= -extent0)
          {
            s0 = tmps0;
            sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
          }
          else
          {
            s0 = -extent0;
            Real tmps1 = -(a01*s0 + b1);
            if (tmps1 < -extent1)
            {
              s1 = -extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
            else if (tmps1 <= extent1)
            {
              s1 = tmps1;
              sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
            else
            {
              s1 = extent1;
              sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
            }
          }
        }
      }
      else  // region 6 (corner)
      { 
        s1 = -extent1;
        Real tmps0 = -(a01*s1 + b0);
        if (tmps0 > extent0)
        {
          s0 = extent0;
          sq_dist = s0*(s0 - ((Real) 2.0)*tmps0) + s1*(s1 + ((Real) 2.0)*b1) + c;
        }
        else if (tmps0 >= -extent0)
        {
          s0 = tmps0;
          sq_dist = -s0*s0 + s1*(s1 + ((Real) 2.0)*b1) + c;
        }
        else
        {
          s0 = -extent0;
          Real tmps1 = -(a01*s0 + b1);
          if (tmps1 < -extent1)
          {
            s1 = -extent1;
            sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
          }
          else if (tmps1 < extent1)
          {
            s1 = tmps1;
            sq_dist = -s1*s1 + s0*(s0 + ((Real) 2.0)*b0) + c;
          }
          else
          {
            s1 = extent1;
            sq_dist = s1*(s1 - ((Real) 2.0)*tmps1) + s0*(s0 + ((Real) 2.0)*b0) + c;
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
    const Real e0pe1 = extent0 + extent1;
    const Real sgn = (a01 > (Real) 0.0 ? (Real) -1.0 : (Real) 1.0);
    const Real b0Avr = ((Real) 0.5)*(b0 - sgn*b1);
    Real lambda = -b0Avr;
    if (lambda < -e0pe1)
      lambda = -e0pe1;
    else if (lambda > e0pe1)
      lambda = e0pe1;
    
    s1 = -sgn*lambda*extent1/e0pe1;
    s0 = lambda + sgn*s1;
    sq_dist = lambda*(lambda + ((Real) 2.0)*b0Avr) + c;
  }

  cp0 = center0 + dir0*s0;
  cp1 = center1 + dir1*s1;

  // account for numerical round-off errors
  if (sq_dist < (Real) 0.0)
    sq_dist = (Real) 0.0;
  return sq_dist;
}

/// Calculates the squared distance between a line segment and point (and computes closest point on line segment)
/**
 * This code adapted from www.geometrictools.com (thanks Dave Eberly!)
 */
Real SSL::calc_sq_dist(const LineSeg3& seg, const Vector3& q, Vector3& cp)
{
  const Vector3& p0 = seg.first;
  const Vector3& p1 = seg.second;
  const Vector3 dir_unnorm = p1 - p0;
  const Real len = dir_unnorm.norm();
  assert(len > NEAR_ZERO);
  const Vector3 dir = dir_unnorm / len;
  const Vector3 center = (p0 + p1) * (Real) 0.5;
  const Vector3 diff = q - center;
  const Real t = dir.dot(diff);
  const Real extent = len * (Real) 0.5;
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
Real SSL::calc_dist(const SSL& a, const SSL& b, Vector3& cpa, Vector3& cpb)
{
  // get the two line segments and calculate the distance between them
  LineSeg3 s1(a.p1, a.p2);
  LineSeg3 s2(b.p1, b.p2);
  Real dist = std::sqrt(calc_sq_dist(s1, s2, cpa, cpb)) - a.radius - b.radius;

  // determine closest points (on the SSL)
  Vector3 dir = cpb - cpa;
  Real dir_norm = dir.norm();
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
Real SSL::calc_dist(const SSL& a, const SSL& b, const Matrix4& aTb, Vector3& cpa, Vector3& cpb)
{
  // create a new SSL (b in a's frame)
  SSL b_a;
  b_a.p1 = aTb.mult_point(b.p1);
  b_a.p2 = aTb.mult_point(b.p2);
  b_a.radius = b.radius;

  return calc_dist(a, b_a, cpa, cpb); 
}

/// Determines whether two SSLs intersect
bool SSL::intersects(const SSL& a, const SSL& b)
{
  Vector3 tmp;
  Real dist = calc_dist(a, b, tmp, tmp);
  return (dist <= (Real) 0.0);
}

/// Determines whether two SSLs intersect
bool SSL::intersects(const SSL& a, const SSL& b, const Matrix4& T)
{
  Vector3 tmp;
  Real dist = calc_dist(a, b, T, tmp, tmp);
  return (dist <= (Real) 0.0);
}

/// Determines whether (and when) a line segment intersects with the SSL
bool SSL::intersects(const SSL& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q)
{
  // determine the distance between seg and the underlying line segment; also
  // get the closest points
  Vector3 cp_seg, cp_SSL;
  Real sq_dist = calc_sq_dist(seg, LineSeg3(a.p1, a.p2), cp_seg, cp_SSL);

  // check for no intersection
  if (std::sqrt(sq_dist) - a.radius > (Real) 0.0)
    return false;

  // we have an intersection with the capsule; determine the first point of
  // intersection
  Vector3 dir_unnormed = seg.first - cp_seg;
  Real dir_unnormed_len = dir_unnormed.norm();
  if (dir_unnormed_len < NEAR_ZERO)
  {
    if (tmin > (Real) 0.0)
      return false;
    else
    {
      tmin = (Real) 0.0;
      q = seg.first;
      return true;
    }
  }
  Vector3 dir = dir_unnormed / dir_unnormed_len;
  q = cp_seg + dir*a.radius;

  // determine the parameter for q
  tmin = (q - seg.first).norm() / (seg.second - seg.first).norm();
  return true;
}

/// Determines whether a point is outside the SSL
bool SSL::outside(const SSL& a, const Vector3& point, Real tol)
{
  Real dist = calc_dist(a, point);
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
Vector3 SSL::get_lower_bounds(const Matrix4& T)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // transform the two points
  Vector3 p1T = T.mult_point(p1);
  Vector3 p2T = T.mult_point(p2);

  // determine the lower left point
  Vector3 ll;
  ll[X] = (p1T[X] < p2T[X]) ? p1T[X] : p2T[X];
  ll[Y] = (p1T[Y] < p2T[Y]) ? p1T[Y] : p2T[Y];
  ll[Z] = (p1T[Z] < p2T[Z]) ? p1T[Z] : p2T[Z];

  // move the closest point toward -inf by radius
  const Vector3 ones((Real) 1.0, (Real) 1.0, (Real) 1.0);
  return ll - ones*radius; 
}

/// Gets the upper bounds of the SSL
Vector3 SSL::get_upper_bounds(const Matrix4& T)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // transform the two points
  Vector3 p1T = T.mult_point(p1);
  Vector3 p2T = T.mult_point(p2);

  // determine the upper right point
  Vector3 ur;
  ur[X] = (p1T[X] > p2T[X]) ? p1T[X] : p2T[X];
  ur[Y] = (p1T[Y] > p2T[Y]) ? p1T[Y] : p2T[Y];
  ur[Z] = (p1T[Z] > p2T[Z]) ? p1T[Z] : p2T[Z];

  // move the closest point toward inf by radius
  const Vector3 ones((Real) 1.0, (Real) 1.0, (Real) 1.0);
  return ur + ones*radius;
}

/// Outputs the SSL to VRML (not yet implemented)
std::ostream& SSL::to_vrml(std::ostream& out, const Matrix4& T) const
{
  throw std::runtime_error("SSL::to_vrml() not implemented!");
  return out;
}

