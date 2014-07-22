/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/OBB.h>
#include <Moby/AABB.h>

using namespace Ravelin;
using namespace Moby;
using std::pair;
using std::endl;

/// Copies an AABB
AABB& AABB::operator=(const AABB& a)
{
  geom = a.geom;
  minp = a.minp;
  maxp = a.maxp;
  return *this;
}

/// Determines whether an OBB and a line/ray/line segment intersect
/**
 * When intersecting, return intersection distance tmin and point q of
 * intersection.
 * \param a the OBB to check for intersection
 * \param seg the line segment to check for intersection
 * \param tmin on entry, contains the minimum value of the line parameter- for
 *        a line segment, this will generally be 0; when an intersection
 *        occurs, this will contain the distance of intersection from the
 *        tmin that was input on return
 * \param tmax the maximum value of the line parameter- for a line segment,
 *        this will generally be 1
 * \param q contains the point of intersection, if any, on return
 * \return <b>true</b> if the OBB and line intersect, <b>false</b> otherwise
 * \note code adapted from [Ericson, 2005], pp. 180-181
 */
bool AABB::intersects(const AABB& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q)
{
  // get the center of the AABB
  Point3d center = a.minp*0.5 + a.maxp*0.5;

  // convert the line segment to AABB space
  Point3d p = seg.first - center;
  Vector3d d = seg.second - center - p;

  // determine half-extents
  Vector3d l = center - a.minp;

  // compute zero tolerance for d
  const double ZERO_TOL = NEAR_ZERO * std::max((double) 1.0, std::max(l.norm(), d.norm()));

  FILE_LOG(LOG_BV) << "AABB::intersects() entered" << endl; 
  FILE_LOG(LOG_BV) << "  -- checking intersection between line segment " << seg.first << " / " << seg.second << " and AABB: " << endl << a;

  // for all three slabs
  for (unsigned i=0; i< 3; i++)
  {
    if (std::fabs(d[i]) < ZERO_TOL)
    {
      // line is parallel to slab; no hit if origin not within slab
      if (p[i] < -l[i] || p[i] > l[i])
      {
        FILE_LOG(LOG_BV) << "  -- seg parallel to slab " << i << " and origin not w/in slab = no intersection" << endl; 
        FILE_LOG(LOG_BV) << "    amin: " << -l << "  amax: " << l << "  p: " << p << endl;
        return false;
      }
    }
    else
    {
      // compute intersection value of line with near and far plane of slab
      double ood = 1.0 / d[i];
      double t1 = (-l[i] - p[i]) * ood;
      double t2 = (l[i] - p[i]) * ood;

      FILE_LOG(LOG_BV) << "  ood: " << ood << " t1: " << t1 << "  t2: " << t2 << endl; 

      // make t1 be the intersection with the near plane, t2 with far plane
      if (t1 > t2)
        std::swap(t1, t2);

      // compute the intersection of slab intersection intervals
      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      // exit with no collision as soon as slab intersection becomes empty
      if (tmin > tmax + ZERO_TOL)
      {
        FILE_LOG(LOG_BV) << "  tmin (" << tmin << ") > tmax (" << tmax << ") -- seg and OBB do not intersect" << endl;

        return false;
      }
    }
  }

  // ray intersects all three slabs; return point q and intersection t value
  q = p + d * tmin;

  FILE_LOG(LOG_BV) << "AABB::intersects() - seg and AABB intersect; first intersection: " << tmin << "(" << q << ")" << endl; 
  FILE_LOG(LOG_BV) << "AABB::intersects() exiting" << endl; 

  return true;
}

/// Transforms the OBB using the given transform
void AABB::transform(const Transform3d& T, BV* result) const
{
  // get the AABB
  AABB& aabb = *((AABB*) result);

  // copy this
  aabb = *this;

  // get the OBB version of this transformed
  OBB o;
  get_OBB().transform(T, &o);

  // now setup the AABB boundary points
  aabb.minp = o.get_lower_bounds();
  aabb.maxp = o.get_lower_bounds();
}

/// Determines whether a point is outside of a AABB
bool AABB::outside(const AABB& a, const Point3d& point, double tol)
{
  const unsigned THREE_D = 3;

  assert(a.get_relative_pose() == point.pose);

  for (unsigned i=0; i< THREE_D; i++)
    if (point[i] < a.minp[i] - tol || point[i] > a.maxp[i] + tol)
      return true;

  return false;
}

/// Determines whether two AABBs overlap
bool AABB::intersects(const AABB& a, const AABB& b)
{
  const unsigned THREE_D = 3;

  assert(a.get_relative_pose() == b.get_relative_pose());

  for (unsigned i=0; i< THREE_D; i++)
    if (a.minp[i] > b.maxp[i] || a.maxp[i] < b.minp[i])
      return false;
    else if (b.minp[i] > a.maxp[i] || b.maxp[i] < a.minp[i])
      return false;

  return true;
}

/// Determines whether two AABBs overlap
bool AABB::intersects(const AABB& a, const AABB& b, const Transform3d& aTb)
{
  // make OBBs out of a and b
  OBB oa = a.get_OBB();
  OBB ob = b.get_OBB();

  // check a against b
  return OBB::intersects(oa, ob, aTb);
}

/// Gets the AABB as an OBB
OBB AABB::get_OBB() const
{
  OBB o;
  o.R = IDENTITY_3x3;
  o.center = maxp*0.5 + minp*0.5;
  o.l = maxp*0.5 - minp*0.5;
  o.geom = geom;
  return o;
}

/// Sends this AABB to VRML
std::ostream& AABB::to_vrml(std::ostream& out, const Pose3d& T) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get translation and axis-angle for T
  const Origin3d& tx = T.x;
  AAngled rot = T.q;

  // pick a random color
  Point3d color((double) rand()/RAND_MAX, (double) rand()/RAND_MAX, (double) rand()/RAND_MAX);

  // setup the vertices of the box
  std::list<Point3d> vertices;
  vertices.push_back(Point3d(-1,-1,+1));
  vertices.push_back(Point3d(+1,-1,+1));
  vertices.push_back(Point3d(-1,+1,+1));
  vertices.push_back(Point3d(+1,+1,+1));
  vertices.push_back(Point3d(-1,+1,-1));
  vertices.push_back(Point3d(+1,+1,-1));
  vertices.push_back(Point3d(+1,-1,-1));
  vertices.push_back(Point3d(-1,-1,-1));

  // determine the scaling factor 
  Point3d l = maxp*0.5 - minp*0.5;

  out << "Transform {" << endl;
  out << "  translation " << tx[X] << " " << tx[Y] << " " << tx[Z] << endl;
  out << "  rotation " << rot.x << " " << rot.y << " " << rot.z << " " << rot.angle << endl;
  out << "  scale " << l[X] << " " << l[Y] << " " << l[Z] << endl;
  out << "  children Shape {" << endl;
  out << "    appearance Appearance { material Material { diffuseColor " << out << color[X] << " " << color[Y] << " " << color[Z] << " } }";
  out << "    geometry IndexedLineSet {" << endl;
  out << "      coord Coordinate { point [";
  BOOST_FOREACH(const Point3d& vertex, vertices)
    out << vertex[X] << " " << vertex[Y] << " " << vertex[Z] << ", ";
  out << "       ] }" << endl;
  out << "      coordIndex [ 0, 1, 3, 2, -1, 1, 6, 5, 3, -1, 0, 2, 4, 7, " << endl;
  out << "                   2, 3, 5, 4, -1, 0, 7, 6, 1, -1, 4, 5, 6, 7, -1] } }} " << endl;

  return out;
}

/// Constructs a velocity expanded OBB using the AABB
BVPtr AABB::calc_swept_BV(CollisionGeometryPtr g, const SVelocityd& v) const
{
  // construct an OBB and use it to determine the velocity-expanded bounding volume
  OBBPtr o(new OBB(get_OBB()));

  return o->calc_swept_BV(g, v);
}

/// Calculates the volume of this AABB
double AABB::calc_volume() const
{
  return (maxp[0]-minp[0])*(maxp[1]-minp[1])*(maxp[2]-minp[2]);
}

/// Gets the lower bounds for this AABB using an OBB
Point3d AABB::get_lower_bounds() const
{
  // construct an OBB and use it to determine the velocity-expanded bounding volume
  return get_OBB().get_lower_bounds();
}

/// Gets the upper bounds for this AABB using an OBB
Point3d AABB::get_upper_bounds() const
{
  // construct an OBB and use it to determine the velocity-expanded bounding volume
  return get_OBB().get_upper_bounds();
}

/// Gets the closest point on the AABB to a point
void AABB::get_closest_point(const AABB& a, const Point3d& p, Point3d& closest)
{
  const unsigned THREE_D = 3; 

  assert(a.get_relative_pose() == p.pose);

  // set closest to p initially
  closest = p;

  // clamp closest to box
  for (unsigned i=0; i< THREE_D; i++)
    if (closest[i] < a.minp[i])
      closest[i] = a.minp[i];
    else if (closest[i] > a.maxp[i])
      closest[i] = a.maxp[i];
}

/// Gets the farthest point on the AABB to a point
/**
 * Returns the squared distance to the farthest point
 */
double AABB::get_farthest_point(const AABB& a, const Point3d& p, Point3d& farthest)
{
  const unsigned X = 0, Y = 1, Z = 2, NVERTS = 8;

  assert(a.get_relative_pose() == p.pose);

  // get the eight vertices of the AABB
  Point3d vertices[NVERTS];
  vertices[0] = Point3d(a.minp[X], a.minp[Y], a.minp[Z]);
  vertices[1] = Point3d(a.minp[X], a.minp[Y], a.maxp[Z]);
  vertices[2] = Point3d(a.minp[X], a.maxp[Y], a.minp[Z]);
  vertices[3] = Point3d(a.minp[X], a.maxp[Y], a.maxp[Z]);
  vertices[4] = Point3d(a.maxp[X], a.minp[Y], a.minp[Z]);
  vertices[5] = Point3d(a.maxp[X], a.minp[Y], a.maxp[Z]);
  vertices[6] = Point3d(a.maxp[X], a.maxp[Y], a.minp[Z]);
  vertices[7] = Point3d(a.maxp[X], a.maxp[Y], a.maxp[Z]);

  // see which is farthest
  double farthest_dist = (vertices[0] - p).norm_sq();
  farthest = vertices[0];
  for (unsigned i=1; i< NVERTS; i++)
  {
    double dist = (vertices[i] - p).norm_sq();
    if (dist > farthest_dist)
    {
      farthest_dist = dist;
      farthest = vertices[i];
    }
  }

  return farthest_dist;
}

