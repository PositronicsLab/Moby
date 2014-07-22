/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _DEFORMABLE_BODY_NODE_H
#define _DEFORMABLE_BODY_NODE_H

#include <Moby/Constants.h>
#include <Ravelin/Vector3d.h>

namespace Moby {

/// Defines a node in a deformable body
struct Node
{
  Node()
  {
    x.set_zero();
    xd.set_zero();
    xdd.set_zero();
    f.set_zero();
    mass = (double) 0.0;
  }

  /// The global position of the node
  Point3d x;

  /// The linear velocity of the node
  Ravelin::Vector3d xd;

  /// The linear acceleration of the node
  Ravelin::Vector3d xdd;

  /// The mass of the node (if any)
  double mass;

  /// The force accumulator for the node
  Ravelin::Vector3d f;
}; // end class

} // end namespace

#endif

