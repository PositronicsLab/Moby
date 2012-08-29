/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POLYHEDRON_H
#define _POLYHEDRON_H

#include <ostream>
#include <vector>
#include <list>
#include <map>
#include <Moby/IndexedTriArray.h>
#include <Moby/Constants.h>
#include <Moby/Matrix4.h>
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
    static PolyhedronPtr minkowski(Polyhedron& p1, const Matrix4& T1, Polyhedron& p2, const Matrix4& T2, bool reflect_p2 = true);
    void operator=(const Polyhedron& p);
    const IndexedTriArray& get_mesh() const { return _mesh; }
    void transform(const Matrix4& T);
    void translate(const Vector3& p);
    const std::vector<Vector3>& get_vertices() const { return _mesh.get_vertices(); }
    const std::vector<IndexedTri>& get_facets() const { return _mesh.get_facets(); }
    bool inside(const Vector3& point, Real tol = NEAR_ZERO);
    bool inside_or_on(const Vector3& point, Real tol = NEAR_ZERO);
    LocationType location(const Vector3& point, Real tol = NEAR_ZERO) const;
    static void to_vrml(std::ostream& out, const Polyhedron& p, Vector3 diffuse_color = Vector3(1,1,1), bool wireframe = false);
    Real calc_volume() const;
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
    Real calc_signed_distance(const Vector3& point, unsigned& closest_facet);

    /// Gets the signed distance from a point to the polyhedron
    Real calc_signed_distance(const Vector3& point) { unsigned discard; return calc_signed_distance(point, discard); }

    /// Gets the corners of the axis-aligned bounding box of this polyhedron
    std::pair<Vector3, Vector3> get_bounding_box_corners() const { return std::make_pair(_bb_min, _bb_max); }

    /// Determines whether this polyhedron convex (to w/in floating point tolerance)
    bool is_convex() { return convexity() < NEAR_ZERO; }

    /// Gets the convexity of this polyhedron 
    /**
     * Convexity values less than epsilon (where epsilon is some number near
     * zero) indicate that the polyhedron is convex; greater values indicate
     * that the polyhedron is non-convex.
     */ 
    Real convexity() { if (!_convexity_computed) determine_convexity(); return _convexity; } 

  private:
    static void replace_edge(const std::vector<Vector3>& v, std::vector<IndexedTri>& f, unsigned a, unsigned b, unsigned c, std::vector<unsigned>& del_list);
    static bool find_vertex(const std::vector<Vector3>& vertices, const Vector3& v);
    void calc_bounding_box();
    static void calc_subexpressions(Real w0, Real w1, Real w2, Real& f1, Real& f2, Real& f3, Real& g0, Real& g1, Real& g2);
    void determine_convexity();  
    static Vector3 intersect_plane(const Vector3& normal, Real d, const Vector3& p1, const Vector3& p2);  
    static void slice(const Polyhedron& p1, const Polyhedron& p2, IndexedTriArray& mesh1, IndexedTriArray& mesh2);
    static void remove_outside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared);
    static void remove_inside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared);
    static bool bisects(const Triangle& a, const Triangle& b);
    static bool bisect(const Triangle& tbi, std::vector<Vector3>& v, std::vector<IndexedTri>& f, unsigned i);
    static unsigned add_vertex(std::vector<Vector3>& vertices, const Vector3& v);

    Vector3 _bb_min, _bb_max;
    IndexedTriArray _mesh;
    Real _convexity;
    bool _convexity_computed;
};

std::ostream& operator<<(std::ostream& out, const Polyhedron& p);

// include inline functions
#include "Polyhedron.inl"

} // end namespace

#endif
