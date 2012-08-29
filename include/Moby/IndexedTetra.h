/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
  static void determine_barycentric_coords(const Vector3& v1, const Vector3& v2, const Vector3& v3, const Vector3& v4, const Vector3& p, Real& u, Real& v, Real& w);
  static Vector3 calc_point(const Vector3& v1, const Vector3& v2, const Vector3& v3, const Vector3& v4, Real u, Real v, Real w);

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

