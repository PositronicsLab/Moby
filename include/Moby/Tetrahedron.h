/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _TETRAHEDRON_
#define _TETRAHEDRON_

#include <Ravelin/LinAlgd.h>
#include <Moby/Types.h>
#include <Moby/Constants.h>
#include <limits>

namespace Moby {

/// An indexed tetrahedron; should be oriented ccw
struct Tetrahedron
{
  Tetrahedron() {  }
  Tetrahedron(const Point3d& v1, const Point3d& v2, const Point3d& v3, const Point3d& v4) { a = v1; b = v2; c = v3; d = v4; }
  double calc_signed_dist(const Point3d& p) const;
  double calc_signed_dist(const Point3d& p, Point3d& cp) const;
  void determine_barycentric_coords(const Point3d& p, double& u, double& v, double& w) const;
  Point3d calc_point(double u, double v, double w) const;
  Point3d calc_centroid() const;
  double calc_volume() const;
  bool outside(const Point3d& p, double tol = NEAR_ZERO) const;

  /// First vertex of the tetrahedron
  Point3d a;

  /// Second vertex of the tetrahedron
  Point3d b;

  /// Third vertex of the tetrahedron
  Point3d c;

  /// Fourth vertex of the tetrahedron
  Point3d d;

  private:
    static Ravelin::LinAlgd _LA;
}; // end struct

} // end namespace

#endif

