/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_PAIRWISE_DIST_INFO_H_
#define _MOBY_PAIRWISE_DIST_INFO_H_

#include <Moby/Types.h>

namespace Moby {

struct PairwiseDistInfo
{
  CollisionGeometryPtr a;  // the first geometry
  CollisionGeometryPtr b;  // the second geometry
  double dist;             // the signed distance
  Point3d pa;              // the closest point on geometry A
  Point3d pb;              // the closest point on geometry B
}; // end struct

} // end namespace

#endif

