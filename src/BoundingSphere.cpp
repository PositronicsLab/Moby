/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/OBB.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/SSL.h>
#include <Moby/BoundingSphere.h>

using namespace Moby;
using boost::shared_ptr;
using boost::const_pointer_cast;
using std::endl;

BoundingSphere::BoundingSphere()
{
  center = ZEROS_3;
  radius = 0.0f;
}

BoundingSphere& BoundingSphere::operator=(const BoundingSphere& source)
{
  center = source.center;
  radius = source.radius;
  return *this;
}

/// Sends the bounding sphere to VRML
std::ostream& BoundingSphere::to_vrml(std::ostream& out, const Matrix4& T) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // determine the new center for the bounding sphere
  Vector3 c = center + T.get_translation();

  // setup a random color
  Vector3 color((Real) rand()/RAND_MAX, (Real) rand()/RAND_MAX,
                (Real) rand()/RAND_MAX);

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
BVPtr BoundingSphere::calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const
{
  // get the corresponding body
  SingleBodyPtr b = g->get_single_body();

  // if the body does not move, just return the OBB
  if (!b->is_enabled() || lv.norm()*dt < NEAR_ZERO)
  {
    FILE_LOG(LOG_BV) << "BoundingSphere::calc_vel_exp_BV() entered" << endl;
    FILE_LOG(LOG_BV) << "  -- using original bounding sphere" << endl;
    FILE_LOG(LOG_BV) << "BoundingSphere::calc_vel_exp_BV() exited" << endl;

    return const_pointer_cast<BoundingSphere>(get_this());
  }

  // otherwise, create a SSL 
  shared_ptr<SSL> ssl(new SSL);
  ssl->p1 = this->center;
  ssl->p2 = this->center + lv*dt;
  ssl->radius = this->radius;
  return ssl;
}

/// Determines whether two bounding spheres intersect
bool BoundingSphere::intersects(const BoundingSphere& s1, const BoundingSphere& s2)
{
  // get the squared distance between the two spheres centers
  Real dist_sq = (s1.center - s2.center).norm_sq();

  // see whether the spheres are overlapping
  return dist_sq < (s1.radius + s2.radius)*(s1.radius + s2.radius);
}

/// Determines whether two bounding spheres intersect
bool BoundingSphere::intersects(const BoundingSphere& s1, const BoundingSphere& s2, const Matrix4& s1Ts2)
{
  // determine transformed s2 center
  Vector3 s2c = s1Ts2.get_translation() + s2.center;

  // get the squared distance between the two spheres centers
  Real dist_sq = (s1.center - s2c).norm_sq();

  // see whether the spheres are overlapping
  return dist_sq < (s1.radius + s2.radius)*(s1.radius + s2.radius);
}

/// Calculates the signed distance between two bounding spheres
Real BoundingSphere::calc_dist(const BoundingSphere& s1, const BoundingSphere& s2)
{
  // get the distance between the sphere centers
  Real dist = (s1.center - s2.center).norm();

  return dist - s1.radius - s2.radius; 
}

/// Determines whether a point is outside of the bounding sphere
bool BoundingSphere::outside(const BoundingSphere& a, const Vector3& point, Real tol)
{
  Vector3 cx = point - a.center;
  return cx.norm_sq() > (Real) a.radius*a.radius + std::sqrt(tol);
}

/// Determines whether there is an intersection between the primitive and a line segment
bool BoundingSphere::intersects(const BoundingSphere& s, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& isect) 
{
  const unsigned X = 0, Y = 1, Z = 2;

  // account for sphere center in translation
  Vector3 p = seg.first - s.center;
  Vector3 q = seg.second - s.center;

  // look for:
  // (seg.first*t + seg.second*(1-t))^2 = r^2

  // use quadratic formula
  const Real px = p[X];
  const Real py = p[Y];
  const Real pz = p[Z];
  const Real qx = q[X];
  const Real qy = q[Y];
  const Real qz = q[Z];
  const Real a = px*px + py*py + pz*pz - 2*px*qx + qx*qx - 2*py*qy + qy*qy -
                 2*pz*qz + qz*qz;
  const Real b = 2*px*qx - 2*qx*qx + 2*py*qy - 2*qy*qy + 2*pz*qz - 2*qz*qz;
  const Real c = qx*qx + qy*qy + qz*qz - s.radius*s.radius;

  // check for no solution
  Real disc = b*b - 4*a*c;
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
  isect = seg.first*tmin + seg.second*((Real) 1.0 - tmin);

  return true;
}

