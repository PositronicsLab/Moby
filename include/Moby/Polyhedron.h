/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _POLYHEDRON_H
#define _POLYHEDRON_H

#include <vector>
#include <Moby/Constants.h>
#include <Moby/Plane.h>
#include <Moby/Log.h>
#include <Moby/NumericalException.h>

// needed for qhull
extern "C"
{
  #include <qhull/qhull_a.h>
}


namespace Moby {

class PolyhedralPrimitive;

/// Represents a polyhedron using a winged-edge (type) data structure
/**
 * Assume the polygons below:
 *      A
 *    /  \
 *   B---C
 *   \  /
 *    D
 * with edge orientations given as (CA), (AB), (CB), (BD), (DC).
 * An original winged edge data structure use a clockwise traversal.
 * In my modified data structure, one looks at the edge orientation and the
 * face to determine the next edge in a traversal. For example, if one wants
 * to traverse the top polygon in clockwise from edge (CB), one would look at
 * the next edge on the right (nextR).
 */
class Polyhedron 
{
  public:
    class Edge;

    struct Feature
    {
      virtual ~Feature() {}
    };

    /// The vertex structure
    struct Vertex : public Feature
    {
      virtual ~Vertex() {}
      Ravelin::Origin3d o;                // the data
      std::list<boost::weak_ptr<Edge> > e; // edges coincident to this vertex 
      boost::shared_ptr<void> data;                // arbitrary user data
    };

    /// The face structure
    struct Face : public boost::enable_shared_from_this<Face>, public Feature
    {
      virtual ~Face() {}
      std::list<boost::weak_ptr<Edge> > e; // edges coincident to this face 
      void find_closest_feature(const Ravelin::Origin3d& p, std::map<boost::shared_ptr<Feature>, double>& distances, double& closest_dist, std::list<boost::shared_ptr<Polyhedron::Feature> >& closest_features);
      Plane get_plane() const;
      boost::shared_ptr<void> data;                // arbitrary user data
    };

    /// The edge structure
    struct Edge : public boost::enable_shared_from_this<Edge>, public Feature
    {
      virtual ~Edge() {}
      boost::shared_ptr<Vertex> v1, v2; // vertices of this edge 
      boost::shared_ptr<Face> faceR, faceL;   // faces coincident to this edge 
      boost::weak_ptr<Edge> prevR, nextR, prevL, nextL;   // faces coincident to edges 
      void find_closest_feature(const Ravelin::Origin3d& p, std::map<boost::shared_ptr<Feature>, double>& distances, double& closest_dist, std::list<boost::shared_ptr<Polyhedron::Feature> >& closest_features);
      boost::shared_ptr<void> data;                // arbitrary user data
    };

    // iterates over the edges in a face
    class EdgeFaceIterator
    {
      public:
        EdgeFaceIterator(boost::shared_ptr<Face> f);
        boost::shared_ptr<Edge> operator*();
        void advance_cw();
        void advance_ccw();
        bool has_cw();
        bool has_ccw();

      private:
        boost::shared_ptr<Face> f;
        boost::shared_ptr<Edge> e, term;
    };

    // iterates over the vertices in a face
    class VertexFaceIterator
    {
      public:
        VertexFaceIterator(boost::shared_ptr<Face> f, bool ccw);
        boost::shared_ptr<Vertex> operator*();
        void advance();
        bool has_next();

      private:
        boost::shared_ptr<Face> f;
        boost::shared_ptr<Edge> e, term;
        boost::shared_ptr<Vertex> v;
        bool ccw;        
    };

    enum LocationType { eInside, eOutside, eOnVertex, eOnEdge, eOnFace };  

    Polyhedron();
    Polyhedron(const Polyhedron& p) { _convexity_computed = false; operator=(p); }
    double find_closest_features(const Ravelin::Origin3d& p, std::list<boost::shared_ptr<Feature> >& closest_features, bool& inside) const;
    static Polyhedron calc_minkowski_diff(boost::shared_ptr<const PolyhedralPrimitive> pA, boost::shared_ptr<const PolyhedralPrimitive> pB, boost::shared_ptr<const Ravelin::Pose3d> poseA, boost::shared_ptr<const Ravelin::Pose3d> poseB);
    Polyhedron& operator=(const Polyhedron& p);
    std::vector<boost::shared_ptr<Vertex> >& get_vertices() { return _vertices; }
    const std::vector<boost::shared_ptr<Vertex> >& get_vertices() const { return _vertices; }
    const std::vector<boost::shared_ptr<Edge> >& get_edges() const { return _edges; }
    const std::vector<boost::shared_ptr<Face> >& get_faces() const { return _faces; }
    bool inside(const Ravelin::Origin3d& point, double tol = NEAR_ZERO);
    bool inside_or_on(const Ravelin::Origin3d& point, double tol = NEAR_ZERO);
    LocationType location(const Ravelin::Origin3d& point, double tol = NEAR_ZERO) const;
    static void to_vrml(std::ostream& out, const Polyhedron& p, Ravelin::Origin3d diffuse_color = Ravelin::Origin3d(1,1,1), bool wireframe = false);
    double calc_volume() const;
    bool degenerate() const;

    template <class ForwardIterator>
    static Polyhedron calc_convex_hull(ForwardIterator begin, ForwardIterator end);

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

    static double sqr(double x) { return x*x; }
    void calc_bounding_box();
    static void calc_subexpressions(double w0, double w1, double w2, double& f1, double& f2, double& f3, double& g0, double& g1, double& g2);
    void determine_convexity();  

    Ravelin::Origin3d _bb_min, _bb_max;
    double _convexity;
    bool _convexity_computed;

    std::vector<boost::shared_ptr<Vertex> > _vertices;
    std::vector<boost::shared_ptr<Face> > _faces;
    std::vector<boost::shared_ptr<Edge> > _edges;
};

std::ostream& operator<<(std::ostream& out, const Polyhedron& m);

#include "Polyhedron.inl"

} // end namespace

#endif
