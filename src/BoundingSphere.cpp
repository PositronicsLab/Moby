/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/OBB.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/SSL.h>
#include <Moby/BoundingSphere.h>

using std::pair;
using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::const_pointer_cast;
using std::endl;

BoundingSphere::BoundingSphere()
{
  center.set_zero();
  radius = 0.0;
}

BoundingSphere& BoundingSphere::operator=(const BoundingSphere& source)
{
  center = source.center;
  radius = source.radius;
  return *this;
}

/// Transforms the BoundingSphere using the given transform
void BoundingSphere::transform(const Transform3d& T, BV* result) const
{
  // get the BoundingSphere 
  BoundingSphere& s = *((BoundingSphere*) result);

  // copy this
  s = *this;

  // transform the center
  s.center = T.transform_point(center);
}

/// Sends the bounding sphere to VRML
std::ostream& BoundingSphere::to_vrml(std::ostream& out, const Pose3d& T) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // determine the new center for the bounding sphere
  Point3d c = center + T.x;

  // setup a random color
  Point3d color((double) rand()/RAND_MAX, (double) rand()/RAND_MAX,
                (double) rand()/RAND_MAX);

  // setup the VRML output
  out << "Transform {" << std::endl;
  out << "  translation " << c[X] << " " << c[Y] << " " << c[Z] << std::endl;
  out << "  children Shape {" << std::endl;
  out << "    appearance Appearance { material Material { " << std::endl;
  out << "      diffuseColor " << color[X] << " " << color[X] << " " << color[Z] << std::endl;
  out << "      transparency 0.9 } }" << std::endl;
  out << "    geometry Sphere { radius " << radius << " } } }" << std::endl;

  return out;
}

/// Calculates the velocity expanded bounding volume for the bounding sphere (calculates an OBB)
BVPtr BoundingSphere::calc_swept_BV(CollisionGeometryPtr g, const SVelocityd& v) const
{
  // get the corresponding body
  SingleBodyPtr b = g->get_single_body();

  // if the body does not move, just return the bounding sphere 
  if (!b->is_enabled())
  {
    FILE_LOG(LOG_BV) << "BoundingSphere::calc_swept_BV() entered" << endl;
    FILE_LOG(LOG_BV) << "  -- using original bounding sphere" << endl;
    FILE_LOG(LOG_BV) << "BoundingSphere::calc_swept_BV() exited" << endl;

    return const_pointer_cast<BoundingSphere>(get_this());
  }

  // get the velocity in the proper frame
  SVelocityd vx = Pose3d::transform(get_relative_pose(), v); 

  // otherwise, create a SSL 
  shared_ptr<SSL> ssl(new SSL);
  ssl->p1 = this->center;
  ssl->p2 = this->center + vx.get_linear();
  ssl->radius = this->radius;
  return ssl;
}

/// Determines whether two bounding spheres intersect
bool BoundingSphere::intersects(const BoundingSphere& s1, const BoundingSphere& s2)
{
  assert(s1.get_relative_pose() == s2.get_relative_pose());

  // get the squared distance between the two spheres centers
  double dist_sq = (s1.center - s2.center).norm_sq();

  // see whether the spheres are overlapping
  return dist_sq < (s1.radius + s2.radius)*(s1.radius + s2.radius);
}

/// Determines whether two bounding spheres intersect
bool BoundingSphere::intersects(const BoundingSphere& s1, const BoundingSphere& s2, const Transform3d& s1Ts2)
{
  assert(s1Ts2.target == s1.get_relative_pose() && 
         s1Ts2.source == s2.get_relative_pose());

  // determine transformed s2 center
  Point3d s2c = s1Ts2.x + s2.center;
  s2c.pose = s1.get_relative_pose();

  // get the squared distance between the two spheres centers
  double dist_sq = (s1.center - s2c).norm_sq();

  // see whether the spheres are overlapping
  return dist_sq < (s1.radius + s2.radius)*(s1.radius + s2.radius);
}

/// Calculates the signed distance between two bounding spheres
double BoundingSphere::calc_dist(const BoundingSphere& s1, const BoundingSphere& s2)
{
  assert(s1.get_relative_pose() == s2.get_relative_pose());

  // get the distance between the sphere centers
  double dist = (s1.center - s2.center).norm();

  return dist - s1.radius - s2.radius; 
}

/// Determines whether a point is outside of the bounding sphere
bool BoundingSphere::outside(const BoundingSphere& a, const Point3d& point, double tol)
{
  Point3d cx = point - a.center;
  return cx.norm_sq() > (double) a.radius*a.radius + std::sqrt(tol);
}

/// Determines whether there is an intersection between the primitive and a line segment
bool BoundingSphere::intersects(const BoundingSphere& s, const LineSeg3& seg, double& tmin, double tmax, Point3d& isect) 
{
  const unsigned X = 0, Y = 1, Z = 2;

  // account for sphere center in translation
  Point3d p = seg.first - s.center;
  Point3d q = seg.second - s.center;

  // look for:
  // (seg.first*t + seg.second*(1-t))^2 = r^2

  // use quadratic formula
  const double px = p[X];
  const double py = p[Y];
  const double pz = p[Z];
  const double qx = q[X];
  const double qy = q[Y];
  const double qz = q[Z];
  const double a = px*px + py*py + pz*pz - 2*px*qx + qx*qx - 2*py*qy + qy*qy -
                 2*pz*qz + qz*qz;
  const double b = 2*px*qx - 2*qx*qx + 2*py*qy - 2*qy*qy + 2*pz*qz - 2*qz*qz;
  const double c = qx*qx + qy*qy + qz*qz - s.radius*s.radius;

  // check for no solution
  double disc = b*b - 4*a*c;
  if (disc < 0.0)
    return false;

  // compute solutions
  disc = std::sqrt(disc);
  tmin = (-b + disc)/(2*a);
  tmax = (-b - disc)/(2*a);

  // look for minimum and maximum solutions
  if (tmax < tmin)
    std::swap(tmin, tmax);

  // compute isect
  isect = seg.first*tmin + seg.second*((double) 1.0 - tmin);

  return true;
}

