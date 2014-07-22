/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _POLYHEDRON_H
#define _POLYHEDRON_H

#include <ostream>
#include <vector>
#include <list>
#include <map>
#include <Moby/IndexedTriArray.h>
#include <Moby/Constants.h>
#include <Moby/InvalidIndexException.h>

namespace Moby {

class Triangle;
  
/// Represents a three-dimensional polyhedron
/**
 * Though this class is used for representing polyhedron, it is possible to represent
 * arbitrary (i.e., non-closed) triangle meshes as well.
 */
class Polyhedron
{
  public:
    enum LocationType { eInside, eOutside, eOnVertex, eOnEdge, eOnFace };  

    Polyhedron() { _convexity_computed = false; }  
    Polyhedron(const IndexedTriArray& mesh);
    Polyhedron(const Polyhedron& p) { _convexity_computed = false; operator=(p); }
    static PolyhedronPtr minkowski(Polyhedron& p1, boost::shared_ptr<const Ravelin::Pose3d> T1, Polyhedron& p2, boost::shared_ptr<const Ravelin::Pose3d> T2, bool reflect_p2 = true);
    void operator=(const Polyhedron& p);
    const IndexedTriArray& get_mesh() const { return _mesh; }
    void transform(const Ravelin::Transform3d& T);
    const std::vector<Ravelin::Origin3d>& get_vertices() const { return _mesh.get_vertices(); }
    const std::vector<IndexedTri>& get_facets() const { return _mesh.get_facets(); }
    bool inside(const Ravelin::Origin3d& point, double tol = NEAR_ZERO);
    bool inside_or_on(const Ravelin::Origin3d& point, double tol = NEAR_ZERO);
    LocationType location(const Ravelin::Origin3d& point, double tol = NEAR_ZERO) const;
    static void to_vrml(std::ostream& out, const Polyhedron& p, Ravelin::Origin3d diffuse_color = Ravelin::Origin3d(1,1,1), bool wireframe = false);
    double calc_volume() const;
    bool consistent() const;
    bool degenerate() const;
    static IndexedTriArray construct_intersection(Polyhedron& p1, Polyhedron& p2);
    static IndexedTriArray construct_union(Polyhedron& p1, Polyhedron& p2);
    static IndexedTriArray construct_difference(Polyhedron& p1, Polyhedron& p2);

    template <class InputIterator1, class InputIterator2>
    Polyhedron(InputIterator1 verts_begin, InputIterator1 verts_end, InputIterator2 facets_begin, InputIterator2 facets_end);

    /// Gets the list of facets coincident to the i'th facet
    const std::list<unsigned>& get_incident_facets(unsigned i) const { return _mesh.get_incident_facets(i); } 
 
    /// Gets the signed distance and closest facet to a point
    double calc_signed_distance(const Ravelin::Origin3d& point, unsigned& closest_facet);

    /// Gets the signed distance from a point to the polyhedron
    double calc_signed_distance(const Ravelin::Origin3d& point) { unsigned discard; return calc_signed_distance(point, discard); }

    /// Gets the corners of the axis-aligned bounding box of this polyhedron
    std::pair<Ravelin::Origin3d, Ravelin::Origin3d> get_bounding_box_corners() const { return std::make_pair(_bb_min, _bb_max); }

    /// Determines whether this polyhedron convex (to w/in floating point tolerance)
    bool is_convex() { return convexity() < NEAR_ZERO; }

    /// Gets the convexity of this polyhedron 
    /**
     * Convexity values less than epsilon (where epsilon is some number near
     * zero) indicate that the polyhedron is convex; greater values indicate
     * that the polyhedron is non-convex.
     */ 
    double convexity() { if (!_convexity_computed) determine_convexity(); return _convexity; } 

  private:
    static void replace_edge(const std::vector<Ravelin::Origin3d>& v, std::vector<IndexedTri>& f, unsigned a, unsigned b, unsigned c, std::vector<unsigned>& del_list);
    static bool find_vertex(const std::vector<Ravelin::Origin3d>& vertices, const Ravelin::Origin3d& v);
    void calc_bounding_box();
    static void calc_subexpressions(double w0, double w1, double w2, double& f1, double& f2, double& f3, double& g0, double& g1, double& g2);
    void determine_convexity();  
    static Ravelin::Origin3d intersect_plane(const Ravelin::Vector3d& normal, double d, const Ravelin::Origin3d& p1, const Ravelin::Origin3d& p2);  
    static void slice(const Polyhedron& p1, const Polyhedron& p2, IndexedTriArray& mesh1, IndexedTriArray& mesh2);
    static void remove_outside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared);
    static void remove_inside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared);
    static bool bisects(const Triangle& a, const Triangle& b);
    static bool bisect(const Triangle& tbi, std::vector<Ravelin::Origin3d>& v, std::vector<IndexedTri>& f, unsigned i);
    static unsigned add_vertex(std::vector<Ravelin::Origin3d>& vertices, const Ravelin::Origin3d& v);

    Ravelin::Origin3d _bb_min, _bb_max;
    IndexedTriArray _mesh;
    double _convexity;
    bool _convexity_computed;
};

std::ostream& operator<<(std::ostream& out, const Polyhedron& p);

// include inline functions
#include "Polyhedron.inl"

} // end namespace

#endif
