/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _POLYHEDRAL_PRIMITIVE_H
#define _POLYHEDRAL_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

/// Defines a triangle-mesh-based primitive type used for inertial property calculation and geometry provisions
/**
 * The center-of-mass of derived types may be at the origin of the world,
 * or not.  Additionally, Primitive can take a transformation matrix in its
 * constructor, with which the primitive data (com, inertia matrix, and geometry)
 * can be transformed.
 */
class PolyhedralPrimitive : public Primitive 
{
  public:

    PolyhedralPrimitive() : Primitive() { }
    PolyhedralPrimitive(const Ravelin::Pose3d& T) : Primitive(T) { }

    // Gets the number of facets in this primitive
    virtual unsigned num_facets() const = 0;

    // Gets the facets in this primitive *in the global frame* such that a point x within the primitive obeys A*x <= b
    virtual void get_facets(boost::shared_ptr<const Ravelin::Pose3d> source_pose, Ravelin::MatrixNd& A, Ravelin::VectorNd& b) const = 0; 
}; // end class

} // end namespace

#endif
