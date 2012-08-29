/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _TETRAHEDRON_
#define _TETRAHEDRON_

#include <Moby/Constants.h>
#include <limits>

namespace Moby {

/// An indexed tetrahedron; should be oriented ccw
struct Tetrahedron
{
  Tetrahedron() {  }
  Tetrahedron(const Vector3& v1, const Vector3& v2, const Vector3& v3, const Vector3& v4) { a = v1; b = v2; c = v3; d = v4; }
  Real calc_signed_dist(const Vector3& p) const;
  void determine_barycentric_coords(const Vector3& p, Real& u, Real& v, Real& w) const;
  Vector3 calc_point(Real u, Real v, Real w) const;
  Vector3 calc_centroid() const;
  Real calc_volume() const;
  bool outside(const Vector3& p, Real tol = NEAR_ZERO) const;

  /// First vertex of the tetrahedron
  Vector3 a;

  /// Second vertex of the tetrahedron
  Vector3 b;

  /// Third vertex of the tetrahedron
  Vector3 c;

  /// Fourth vertex of the tetrahedron
  Vector3 d;
}; // end struct

} // end namespace

#endif

