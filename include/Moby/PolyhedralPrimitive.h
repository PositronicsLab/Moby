/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _POLYHEDRAL_PRIMITIVE_H
#define _POLYHEDRAL_PRIMITIVE_H

#include <Moby/Primitive.h>
#include <Moby/Polyhedron.h>

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
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const;

    /// Gets the polyhedron corresponding to this primitive (in its transformed state)
    const Polyhedron& get_polyhedron() const { return _poly; }

    // Gets the number of facets in this primitive
    virtual unsigned num_facets() const = 0;

/*
    // Gets the bounding radius of this primitive
    virtual double get_bounding_radius() const
    {
      // get the vertices
      std::vector<Point3d> verts;
      get_vertices(verts);
      if (verts.empty())
        return 0.0;

      // find which point is closest
      double max_dist = 0.0;
      for (unsigned i=0; i< verts.size(); i++)
        max_dist = std::max(max_dist, verts[i].norm()); 

      return max_dist;      
    }
*/

  template <class OutputIterator>
  static OutputIterator get_halfspaces(const Polyhedron& poly, boost::shared_ptr<const Ravelin::Pose3d> pose, const Ravelin::Transform3d& wTpose, OutputIterator output_begin);

  protected:
    double calc_signed_dist(boost::shared_ptr<const PolyhedralPrimitive> p, Point3d& pthis, Point3d& pp) const;
    Polyhedron _poly;
}; // end class

#include "PolyhedralPrimitive.inl"

} // end namespace

#endif
