/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/Matrix3d.h>
#include <Moby/Triangle.h>
#include <Moby/Tetrahedron.h>
#include <Moby/Log.h>

using namespace Ravelin;
using namespace Moby;

LinAlgd Tetrahedron::_LA;

/// Calculates the signed distance from a point to the tetrahedron
double Tetrahedron::calc_signed_dist(const Point3d& p, Point3d& closest) const
{
  Point3d closest_abc, closest_bdc, closest_dac, closest_dba;

  // form necessary triangles
  Triangle abc(a, b, c);
  Triangle bdc(b, d, c);
  Triangle dac(d, a, c);
  Triangle dba(d, b, a);

  // calculate the signed distances
  double dist_abc = abc.calc_signed_dist(p, closest_abc);
  double dist_bdc = bdc.calc_signed_dist(p, closest_bdc);
  double dist_dac = dac.calc_signed_dist(p, closest_dac);
  double dist_dba = dba.calc_signed_dist(p, closest_dba);

  // determine the minimum distance
  double min_dist = dist_abc;
  closest = closest_abc;
  if (std::fabs(dist_bdc) < std::fabs(min_dist))
  {
    if (min_dist > 0.0 || (min_dist < 0.0 && dist_bdc < 0.0))
    {
      min_dist = dist_bdc;
      closest = closest_bdc;
    }
  }
  if (std::fabs(dist_dac) < std::fabs(min_dist))
  {
    if (min_dist > 0.0 || (min_dist < 0.0 && dist_dac < 0.0))
    {
      min_dist = dist_dac;
      closest = closest_dac;
    }
  }
  if (std::fabs(dist_dba) < std::fabs(min_dist))
  {
    if (min_dist > 0.0 || (min_dist < 0.0 && dist_dba < 0.0))
    {
      min_dist = dist_dba;
      closest = closest_dba;
    }
  }

  FILE_LOG(LOG_COLDET) << "Tetrahedron::calc_signed_dist() entered " << std::endl;
  FILE_LOG(LOG_COLDET) << "  query point is " << p << std::endl;
  FILE_LOG(LOG_COLDET) << "  signed distance from triangle abc: " << dist_abc << std::endl;
  FILE_LOG(LOG_COLDET) << "  signed distance from triangle bdc: " << dist_bdc << std::endl;
  FILE_LOG(LOG_COLDET) << "  signed distance from triangle dac: " << dist_dac << std::endl;
  FILE_LOG(LOG_COLDET) << "  signed distance from triangle dba: " << dist_dba << std::endl;
  FILE_LOG(LOG_COLDET) << "  signed distance: " << min_dist << std::endl;
  FILE_LOG(LOG_COLDET) << "Tetrahedron::calc_signed_dist() exited" << std::endl;

  return min_dist;
}

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
  M.set_column(X, a - d);
  M.set_column(Y, b - d);
  M.set_column(Z, c - d);
  return (M * bary) + d;
}

/// Determines the barycentric coordinates of a point in space
void Tetrahedron::determine_barycentric_coords(const Point3d& px, double& u, double& v, double& w) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // transform point to tetrahedron space
  Point3d p = Pose3d::transform_point(a.pose, px);

  // Form system of linear equations
  // a(x)*u + b(x)*v + c(x)*w + d(x)*(1-u-v-w) = p(x)
  // a(y)*u + b(y)*v + c(y)*w + d(y)*(1-u-v-w) = p(y)
  // a(z)*u + b(z)*v + c(z)*w + d(z)*(1-u-v-w) = p(z)
  // This yields:
  // | ax-dx bx-dx cx-dx | * | u | = | px - dz |
  // | ay-dy by-dy cy-dy | * | v | = | py - dz |
  // | az-dz bz-dz cz-dz | * | w | = | pz - dz |

  // algorithm taken from Real Time Physics course notes (Mueller, et al.)
  Matrix3d M;
  M.set_column(X, a - d);
  M.set_column(Y, b - d);
  M.set_column(Z, c - d);
  Origin3d rhs(p - d);
  Origin3d bary = rhs;
  _LA.solve_fast(M, bary);
  if ((M * bary - rhs).norm() > NEAR_ZERO)
  {
    bary = rhs;
    M.set_column(X, a - d);
    M.set_column(Y, b - d);
    M.set_column(Z, c - d);
    Matrix3d U, V;
    Origin3d S;
    _LA.svd(M, U, S, V);
    _LA.solve_LS_fast(U, S, V, bary);
  }

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

/// Calculates the signed volume of a tetrahedron
double Tetrahedron::calc_volume() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup a matrix
  Matrix3d M;
  M.set_column(X, a - d);
  M.set_column(Y, b - d);
  M.set_column(Z, c - d);

  // compute the determinant
  double vol6 = M.det();

  return vol6/6.0;
}

