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
    virtual double calc_dist_and_normal(const Point3d& p, std::vector<Ravelin::Vector3d>& normals) const;
    virtual osg::Node* create_visualization();
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void set_polyhedron(const Polyhedron& p);
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P);
    virtual bool is_convex() const;

    /// Gets the polyhedron corresponding to this primitive (in its transformed state)
    const Polyhedron& get_polyhedron() const { return _poly; }

    // Gets the number of facets in this primitive
    virtual unsigned num_facets() const { return _poly.get_faces().size();}

    // Gets the bounding radius of this primitive
    virtual double get_bounding_radius() const
    {
      // get the vertices
      const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& verts = _poly.get_vertices();
      if (verts.empty())
        return 0.0;

      // find which point is closest
      double max_dist = 0.0;
      for (unsigned i=0; i< verts.size(); i++)
        max_dist = std::max(max_dist, verts[i]->o.norm()); 

      return max_dist;      
    }

  template <class OutputIterator>
  static OutputIterator get_halfspaces(const Polyhedron& poly, boost::shared_ptr<const Ravelin::Pose3d> pose, const Ravelin::Transform3d& wTpose, OutputIterator output_begin);

  protected:
    void calc_mass_properties();
    double calc_signed_dist(boost::shared_ptr<const PolyhedralPrimitive> p, Point3d& pthis, Point3d& pp) const;
    Polyhedron _poly;
}; // end class

#include "PolyhedralPrimitive.inl"

} // end namespace

#endif
