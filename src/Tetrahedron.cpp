/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Triangle.h>
#include <Moby/Matrix3.h>
#include <Moby/Tetrahedron.h>

using namespace Moby;

/// Calculates the signed distance from a point to the tetrahedron
Real Tetrahedron::calc_signed_dist(const Vector3& p) const
{
  Vector3 closest;

  // form necessary triangles
  Triangle abc(a, b, c);
  Triangle bdc(b, d, c);
  Triangle dac(d, a, c);
  Triangle dba(d, b, a);

  // calculate the squared distance
  Real dist_abc = Triangle::calc_sq_dist(abc, p, closest);
  Real dist_bdc = Triangle::calc_sq_dist(bdc, p, closest);
  Real dist_dac = Triangle::calc_sq_dist(dac, p, closest);
  Real dist_dba = Triangle::calc_sq_dist(dba, p, closest);

  // determine the minimum distance
  Real min_dist = std::max((Real) 0.0, std::min(dist_abc, std::min(dist_bdc, std::min(dist_dac, dist_dba))));
  min_dist = std::sqrt(min_dist);

  if (!outside(p))
    min_dist = -min_dist;

  return min_dist;
}

/// Determines whether a point is outside the tetrahedron
/**
 * \note assumes tetrahedron is oriented ccw
 */
bool Tetrahedron::outside(const Vector3& p, Real tol) const
{
  // test triangle abc
  Triangle abc(a, b, c);
  Vector3 normal = abc.calc_normal();
  Real offset = abc.calc_offset(normal);
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
Vector3 Tetrahedron::calc_point(Real u, Real v, Real w) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  Vector3 bary(u, v, w);
  Matrix3 M;
  M.set_column(X, b - a);
  M.set_column(Y, c - a);
  M.set_column(Z, d - a);
  return (M * bary) + a;
}

/// Determines the barycentric coordinates of a point in space
void Tetrahedron::determine_barycentric_coords(const Vector3& p, Real& u, Real& v, Real& w) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // algorithm taken from Real Time Physics course notes (Mueller, et al.)
  Matrix3 M;
  M.set_column(X, b - a);
  M.set_column(Y, c - a);
  M.set_column(Z, d - a);
  M.inverse();

  Vector3 bary = M * (p - a);
  u = bary[X];
  v = bary[Y];
  w = bary[Z];

  assert(!std::isnan(u));
  assert(!std::isnan(v));
  assert(!std::isnan(w));
}

/// Calculates the centroid of a tetrahedron
Vector3 Tetrahedron::calc_centroid() const
{
  Vector3 centroid = a + b + c + d;
  centroid *= 0.25;
  return centroid;
}

/// Calculates the volume of a tetrahedron
Real Tetrahedron::calc_volume() const
{
  // create a triangle from vertices abc
  Triangle t(a,b,c);

  // get the normal to abc, and the offset
  Vector3 normal = t.calc_normal();
  Real offset = t.calc_offset(normal);

  // project d onto the plane of abc to get the height
  Real height = std::fabs(normal.dot(d) - offset);

  // return the volume
  return t.calc_area() * height / 3;
}

