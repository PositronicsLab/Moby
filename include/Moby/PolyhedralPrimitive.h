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

    /// Gets the polyhedron corresponding to this primitive (in its transformed state)
    Polyhedron get_polyhedron() const { return _poly; }

    /// Sets the pose for this polyhedral primitive, which transforms the underlying polyhedron 
    virtual void set_pose(const Ravelin::Pose3d& P)
    {
      boost::shared_ptr<Ravelin::Pose3d> Pp(new Ravelin::Pose3d(P));

      // get the transform from the global pose to the old pose and from
      // the new pose to the global frame

      // apply the transform from the global pose to the old pose 
      Ravelin::Transform3d wTPp = Ravelin::Pose3d::calc_relative_pose(Pp, GLOBAL);
      Ravelin::Transform3d FTw = Ravelin::Pose3d::calc_relative_pose(GLOBAL, _F);

      // get the vertices
      std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vertices = _poly.get_vertices();
      for (unsigned i=0; i< vertices.size(); i++)
      {
        Point3d p(vertices[i]->o, GLOBAL);
        p = FTw.transform_point(p);
        p = wTPp.transform_point(p);
        vertices[i]->o = Ravelin::Origin3d(p);
      }
    }

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
  protected:
    Polyhedron _poly;
}; // end class

} // end namespace

#endif
