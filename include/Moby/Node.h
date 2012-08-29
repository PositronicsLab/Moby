/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _DEFORMABLE_BODY_NODE_H
#define _DEFORMABLE_BODY_NODE_H

#include <Moby/Constants.h>
#include <Moby/Vector3.h>

namespace Moby {

/// Defines a node in a deformable body
struct Node
{
  Node()
  {
    x = xd = xdd = f = ZEROS_3;
    mass = (Real) 0.0;
  }

  /// The global position of the node
  Vector3 x;

  /// The linear velocity of the node
  Vector3 xd;

  /// The linear acceleration of the node
  Vector3 xdd;

  /// The mass of the node (if any)
  Real mass;

  /// The force accumulator for the node
  Vector3 f;
}; // end class

} // end namespace

#endif

