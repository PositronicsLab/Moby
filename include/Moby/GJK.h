/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _GJK_H
#define _GJK_H

#include <map>
#include <vector>
namespace Moby {

/// An implementation of the GJK algorithm 
class GJK
{
  private:
    struct SVertex
    {
      Ravelin::Origin3d v;      // vA - vB: the vertex in the simplex
      Ravelin::Origin3d vA, vB; // vertices from geometry A and B
    };

    class Simplex
    {
      public:
        enum SimplexType { ePoint, eSegment, eTriangle, eTetra };
        SimplexType get_simplex_type() const { return _type; }

        Simplex(const SVertex& p);

      private:
        SimplexType _type;
        SVertex _p1, _p2, _p3, _p4;
    }

  
}; // end class

} // end namespace

#endif
