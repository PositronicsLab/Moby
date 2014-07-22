/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _INDEXED_TETRA_
#define _INDEXED_TETRA_

#include <limits>

namespace Moby {

/// An indexed tetrahedron
struct IndexedTetra
{
  IndexedTetra() { a = b = c = d = std::numeric_limits<unsigned>::max(); }
  IndexedTetra(unsigned v1, unsigned v2, unsigned v3, unsigned v4) { a = v1; b = v2; c = v3; d = v4; }
  static void determine_barycentric_coords(const Point3d& v1, const Point3d& v2, const Point3d& v3, const Point3d& v4, const Point3d& p, double& u, double& v, double& w);
  static Point3d calc_point(const Point3d& v1, const Point3d& v2, const Point3d& v3, const Point3d& v4, double u, double v, double w);

  /// Index of the first vertex of the tetrahedron
  unsigned a;

  /// Index of the second vertex of the tetrahedron
  unsigned b;

  /// Index of the third vertex of the tetrahedron
  unsigned c;

  /// Index of the fourth vertex of the tetrahedron
  unsigned d;

}; // end struct

} // end namespace

#endif

