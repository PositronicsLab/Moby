/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _INDEXED_TRI_
#define _INDEXED_TRI_

#include <limits>

namespace Moby {

/// An indexed triangle
struct IndexedTri
{
  IndexedTri() { a = b = c = std::numeric_limits<unsigned>::max(); }
  IndexedTri(unsigned v1, unsigned v2, unsigned v3) { a = v1; b = v2; c = v3; }

  /// Reverses the orientation of the triangle
  IndexedTri& reverse() { std::swap(b, c); return *this; }

  /// Index of the first vertex of the triangle
  unsigned a;

  /// Index of the second vertex of the triangle
  unsigned b;

  /// Index of the third vertex of the triangle
  unsigned c;

  bool operator<(const IndexedTri& t) const 
  {
    if (a < t.a)
      return true;
    else if (a == t.a)
    {
      if (b < t.b)
        return true;
      else if (b == t.b && c < t.c)
        return true;
    }
    return false;
  }
};

} // end namespace

#endif
