/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <Moby/Log.h>
#include <Moby/NumericalException.h>
#include <Moby/ThickTriangle.h>

using namespace Moby;
using std::endl;

/// Constructs a thick triangle from a triangle
void ThickTriangle::construct_from_triangle(const Triangle& tri, Real epsilon)
{
  const unsigned TRI_VERTS = 3;

  // save the triangle
  this->tri = tri;

  // get the normal of the triangle
  _normal = tri.calc_normal();
  if (_normal.norm() < NEAR_ZERO)
    throw NumericalException("Triangle normal is zero length!");

  // get the direction of the three edges of the triangle
  for (unsigned i=0; i< TRI_VERTS; i++)
  {
    unsigned j = (i < TRI_VERTS-1) ? i+1 : 0;
    Vector3 e = Vector3::normalize(tri.get_vertex(j) - tri.get_vertex(i));
    Vector3 pn = Vector3::normalize(Vector3::cross(e, _normal));
    Real offset = Vector3::dot(pn, tri.get_vertex(j)) + epsilon;
    _planes.push_back(Plane(pn, offset));
  }

  // compute the offset of the triangle
  Real offset = tri.calc_offset(_normal);

  // add the biggest planes last
  _planes.push_back(Plane(_normal, offset+epsilon));
  _planes.push_back(Plane(-_normal, -offset+epsilon));
}

/// Determines the normal to this thick triangle
Vector3 ThickTriangle::determine_normal(const Vector3& p) const
{
  std::list<Plane>::const_reverse_iterator i = _planes.rbegin();

  // determine which plane is closest -- the negative plane
  Real d1 = std::fabs(i->calc_signed_distance(p));
  i++;
  // ... or the positive one
  Real d2 = std::fabs(i->calc_signed_distance(p));

  FILE_LOG(LOG_COLDET) << "ThickTriangle::determine_normal() +normal " << _normal << " d-: " << d1 << " d+: " << d2 << endl;

  return (d1 < d2) ? -_normal : _normal;
}

/// Determines whether a point is on or inside a thick triangle
bool ThickTriangle::point_inside(const Vector3& point) const
{
  BOOST_FOREACH(const Plane& plane, _planes)
    if (plane.calc_signed_distance(point) >= 0.0)
      return false;

  return true;
}

/// Determines the first point of intersection (if any) between the thick triangle and the line segment
/**
 * \param seg the line segment
 * \param tnear the parameter of the closest point of intersection (if any) on 
 *        success (isect = seg.first + (seg.second - seg.first)*tnear
 * \param isect the closest point of intersection (if any) on success
 * \return <b>true</b> if segment intersects thick triangle; 
 *         <b>false</b> otherwise
 * Algorithm taken from [Ericsson, 2005] 
 */
bool ThickTriangle::intersect_seg(const LineSeg3& seg, Real& tnear, Vector3& isect) const
{
  // determine the line segment parametrically
  const Vector3& p0 = seg.first;
  Vector3 dir = seg.second - seg.first;

  // init tnear and tfar
  tnear = 0; 
  Real tfar = 1; 

  FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() entered" << endl;

  // check all planes
  BOOST_FOREACH(const Plane& plane, _planes)
  {
    const Vector3& Pn = plane.get_normal();

    Real dist = -Pn.dot(p0) + plane.offset;
    Real denom = Pn.dot(dir);

    // early exit: seg parallel to plane
    if (denom == 0.0)
    {
      // if seg origin outside the plane's half space, exit
      if (dist < 0.0)
      {
        FILE_LOG(LOG_COLDET) << "  seg parallel to plane and seg origin outside of plane's halfspace" << endl;
        FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() exited" << endl;
        return false;
      }
    }
    else
    {
      // compute parameterized t value for intersection with current plane
      Real t = dist / denom;
  
      // categorize plane as front-facing or back facing
      if (denom < 0.0)
      {
        // update tfar
        if (t > tnear)
        {
          tnear = t;
          if (tfar < tnear)
          {
            FILE_LOG(LOG_COLDET) << "  tfar (" << tfar << ") < tnear (" << tnear << "): no intersection" << endl;
            FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() exited" << endl;
            return false;
          }
        }
      }
      else if (t < tfar)
      {
        tfar = t;
        if (tfar < tnear)
        {
          FILE_LOG(LOG_COLDET) << "  tfar (" << tfar << ") < tnear (" << tnear << "): no intersection" << endl;
          FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() exited" << endl;
          return false;
        }
      }
    }
  }

  // still here?  successful intersection
  isect = p0 + dir*tnear;

  FILE_LOG(LOG_COLDET) << "  point of intersection: " << isect << endl;
  FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() exited" << endl;

  return true;
}

