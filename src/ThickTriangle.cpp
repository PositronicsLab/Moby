/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <Moby/Log.h>
#include <Moby/NumericalException.h>
#include <Moby/ThickTriangle.h>

using namespace Ravelin;
using namespace Moby;
using std::endl;

/// Gets the pose from a thick triangle
boost::shared_ptr<const Ravelin::Pose3d> ThickTriangle::get_pose() const
{
  assert(!_planes.empty());
  boost::shared_ptr<const Ravelin::Pose3d> pose = _planes.front().get_pose();

  #ifndef NDEBUG
  for (std::list<Plane>::const_iterator i = _planes.begin(); i != _planes.end(); i++)
    assert(i->get_pose() == pose);
  #endif

  return pose;
}

/// Constructs a thick triangle from a triangle
void ThickTriangle::construct_from_triangle(const Triangle& tri, double epsilon)
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
    Vector3d e = Vector3d::normalize(tri.get_vertex(j) - tri.get_vertex(i));
    Vector3d pn = Vector3d::normalize(Vector3d::cross(e, _normal));
    double offset = Vector3d::dot(pn, tri.get_vertex(j)) + epsilon;
    _planes.push_back(Plane(pn, offset));
  }

  // compute the offset of the triangle
  double offset = tri.calc_offset(_normal);

  // add the biggest planes last
  _planes.push_back(Plane(_normal, offset+epsilon));
  _planes.push_back(Plane(-_normal, -offset+epsilon));
}

/// Determines the normal to this thick triangle
Vector3d ThickTriangle::determine_normal(const Point3d& point) const
{
  // convert the point to the thick triangle's space
  Point3d p = Pose3d::transform_point(get_pose(), point);

  std::list<Plane>::const_reverse_iterator i = _planes.rbegin();

  // determine which plane is closest -- the negative plane
  double d1 = std::fabs(i->calc_signed_distance(p));
  i++;
  // ... or the positive one
  double d2 = std::fabs(i->calc_signed_distance(p));

  FILE_LOG(LOG_COLDET) << "ThickTriangle::determine_normal() +normal " << _normal << " d-: " << d1 << " d+: " << d2 << endl;

  return (d1 < d2) ? -_normal : _normal;
}

/// Determines whether a point is on or inside a thick triangle
bool ThickTriangle::point_inside(const Point3d& point) const
{
  // convert the point to the thick triangle's space
  Point3d p = Pose3d::transform_point(get_pose(), point);

  BOOST_FOREACH(const Plane& plane, _planes)
    if (plane.calc_signed_distance(p) >= 0.0)
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
bool ThickTriangle::intersect_seg(const LineSeg3& seg, double& tnear, Point3d& isect) const
{
  // convert the segment to the thick triangle's space
  boost::shared_ptr<const Pose3d> pose = get_pose();
  Point3d s1 = Pose3d::transform_point(pose, seg.first);
  Point3d s2 = Pose3d::transform_point(pose, seg.second);

  // determine the line segment parametrically
  const Point3d& p0 = s1;
  Vector3d dir = s2 - s1;

  // init tnear and tfar
  tnear = 0; 
  double tfar = 1; 

  FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() entered" << endl;

  // check all planes
  BOOST_FOREACH(const Plane& plane, _planes)
  {
    const Vector3d& Pn = plane.get_normal();

    double dist = -Pn.dot(p0) + plane.offset;
    double denom = Pn.dot(dir);

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
      double t = dist / denom;
  
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
  isect = Pose3d::transform_point(seg.first.pose, p0 + dir*tnear);

  FILE_LOG(LOG_COLDET) << "  point of intersection: " << isect << endl;
  FILE_LOG(LOG_COLDET) << "ThickTriangle::intersect_seg() exited" << endl;

  return true;
}

