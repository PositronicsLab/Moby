/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Matrix3d.h>
#include <Moby/Triangle.h>
#include <Moby/Tetrahedron.h>

using namespace Ravelin;
using namespace Moby;

/// Calculates the signed distance from a point to the tetrahedron
double Tetrahedron::calc_signed_dist(const Point3d& p) const
{
  Point3d closest;

  // form necessary triangles
  Triangle abc(a, b, c);
  Triangle bdc(b, d, c);
  Triangle dac(d, a, c);
  Triangle dba(d, b, a);

  // calculate the squared distance
  double dist_abc = Triangle::calc_sq_dist(abc, p, closest);
  double dist_bdc = Triangle::calc_sq_dist(bdc, p, closest);
  double dist_dac = Triangle::calc_sq_dist(dac, p, closest);
  double dist_dba = Triangle::calc_sq_dist(dba, p, closest);

  // determine the minimum distance
  double min_dist = std::max((double) 0.0, std::min(dist_abc, std::min(dist_bdc, std::min(dist_dac, dist_dba))));
  min_dist = std::sqrt(min_dist);

  if (!outside(p))
    min_dist = -min_dist;

  return min_dist;
}

/// Determines whether a point is outside the tetrahedron
/**
 * \note assumes tetrahedron is oriented ccw
 */
bool Tetrahedron::outside(const Point3d& p, double tol) const
{
  // test triangle abc
  Triangle abc(a, b, c);
  Vector3d normal = abc.calc_normal();
  double offset = abc.calc_offset(normal);
  if (p.dot(normal) - offset > tol)
    return true;

  // test triangle bdc
  Triangle bdc(b, d, c);
  normal = bdc.calc_normal();
  offset = bdc.calc_offset(normal);
  if (p.dot(normal) - offset > tol)
    return true;

  // test triangle dac
  Triangle dac(d, a, c);
  normal = dac.calc_normal();
  offset = dac.calc_offset(normal);
  if (p.dot(normal) - offset > tol)
    return true;
  
  // test triangle dba
  Triangle dba(d, b, a);
  normal = dba.calc_normal();
  offset = dba.calc_offset(normal);
  return p.dot(normal) - offset > tol;
}

/// Determines the point corresponding to the barycentric coordinates
Point3d Tetrahedron::calc_point(double u, double v, double w) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  Origin3d bary(u, v, w);
  Matrix3d M;
  M.set_column(X, b - a);
  M.set_column(Y, c - a);
  M.set_column(Z, d - a);
  return (M * bary) + a;
}

/// Determines the barycentric coordinates of a point in space
void Tetrahedron::determine_barycentric_coords(const Point3d& px, double& u, double& v, double& w) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // transform point to tetrahedron space
  Point3d p = Pose3d::transform(a.pose, px);

  // Form system of linear equations
  // a(x)*u + b(x)*v + c(x)*w + d(x)*(1-u-v-w) = p(x)
  // a(y)*u + b(y)*v + c(y)*w + d(y)*(1-u-v-w) = p(y)
  // a(z)*u + b(z)*v + c(z)*w + d(z)*(1-u-v-w) = p(z)
  // This yields:
  // | ax-dx bx-dx cx-dx | * | u | = | px - dz |
  // | ay-dy by-dy cy-dy | * | v | = | py - dz |
  // | az-dz bz-dz cz-dz | * | w | = | pz - dz |
  Matrix3d M;
  M.set_column(X, a - d);
  M.set_column(Y, b - d);
  M.set_column(Z, c - d);
  M.invert();
  Origin3d bary = M * Origin3d(p - d);

/*
  // algorithm taken from Real Time Physics course notes (Mueller, et al.)
  Matrix3d M;
  M.set_column(X, b - a);
  M.set_column(Y, c - a);
  M.set_column(Z, d - a);
  M.invert();

  Point3d bary(M * Origin3d(p - a), a.pose);
*/

  // compute barycentric coordinates
  u = bary[X];
  v = bary[Y];
  w = bary[Z];

  // verify that barycentric coordinate is good and the points are nearby
  assert(!std::isnan(u));
  assert(!std::isnan(v));
  assert(!std::isnan(w));
  assert((calc_point(u, v, w) - p).norm() < NEAR_ZERO); 
}

/// Calculates the centroid of a tetrahedron
Point3d Tetrahedron::calc_centroid() const
{
  Point3d centroid = a + b + c + d;
  centroid *= 0.25;
  return centroid;
}

/// Calculates the volume of a tetrahedron
double Tetrahedron::calc_volume() const
{
  // create a triangle from vertices abc
  Triangle t(a,b,c);

  // get the normal to abc, and the offset
  Vector3d normal = t.calc_normal();
  double offset = t.calc_offset(normal);

  // project d onto the plane of abc to get the height
  double height = std::fabs(normal.dot(d) - offset);

  // return the volume
  return t.calc_area() * height / 3;
}

