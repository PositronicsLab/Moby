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
      std::list<boost::weak_ptr<Edge> > e; // edges coincident to this face (ccw ordering) 
      Plane get_plane() const;
      boost::shared_ptr<void> data;                // arbitrary user data
    };

    /// The edge structure
    struct Edge : public boost::enable_shared_from_this<Edge>, public Feature
    {
      virtual ~Edge() {}
      boost::shared_ptr<Vertex> v1, v2; // vertices of this edge 
      boost::shared_ptr<Face> face1, face2;   // faces coincident to this edge 
      boost::shared_ptr<void> data;                // arbitrary user data
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
        boost::shared_ptr<Vertex> v;
        std::list<boost::weak_ptr<Edge> >::const_reverse_iterator cw_iter;
        std::list<boost::weak_ptr<Edge> >::const_iterator ccw_iter;
        bool ccw;
    };

    enum LocationType { eInside, eOutside, eOnVertex, eOnEdge, eOnFace };  

    Polyhedron();
    Polyhedron(const Polyhedron& p) { _convexity_computed = false; operator=(p); }
    double vclip(boost::shared_ptr<const PolyhedralPrimitive> pA, boost::shared_ptr<const PolyhedralPrimitive> pB, boost::shared_ptr<const Ravelin::Pose3d> poseA, boost::shared_ptr<const Ravelin::Pose3d> poseB, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB);
    static Polyhedron calc_minkowski_diff(boost::shared_ptr<const PolyhedralPrimitive> pA, boost::shared_ptr<const PolyhedralPrimitive> pB, boost::shared_ptr<const Ravelin::Pose3d> poseA, boost::shared_ptr<const Ravelin::Pose3d> poseB); 
/*
    double find_closest_features(const Ravelin::Origin3d& p, std::list<boost::shared_ptr<Feature> >& closest_features, bool& inside) const;
*/
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

    /// Gets the Voronoi plane from two input features
    /**
     * The function takes the tyoes and the pointers of the two features and returns a plane
     * If a point is a positive distance away from the plane, then the point is closer to the first feature
     * If a point is a negative distance away from the plane, then the point is closer to the second feature
     */
   
    enum FeatureType { eVertex, eEdge, eFace };
    static double calc_dist(FeatureType fA, FeatureType fB, boost::shared_ptr<const Polyhedron::Feature> closestA, boost::shared_ptr<const Polyhedron::Feature> closestB, Ravelin::Transform3d& aTb);
  private:
    static boost::shared_ptr<Plane> voronoi_plane (FeatureType fA, FeatureType fB, boost::shared_ptr<const Ravelin::Pose3d> pose, boost::shared_ptr<const Polyhedron::Feature>& featureA, boost::shared_ptr<const Polyhedron::Feature>& featureB);
    static bool clip_edge(boost::shared_ptr<const Polyhedron::Edge> edge, Ravelin::Transform3d fTe, double& min_lambda, double& max_lambda, boost::shared_ptr<const Polyhedron::Feature >& min_N, boost::shared_ptr<const Polyhedron::Feature >& max_N, const std::list<std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > >& planes_neighbors);
    bool post_clip_deriv_check(FeatureType& fX, boost::shared_ptr<const Polyhedron::Feature >& X , boost::shared_ptr<const Polyhedron::Edge> edge, Ravelin::Transform3d& xTe, double& min_lambda, double& max_lambda, boost::shared_ptr<const Polyhedron::Feature >& min_N, boost::shared_ptr<const Polyhedron::Feature >& max_N);
    
    //void promote_featrues(FeatureType& fA, FeatureType& fB, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB, Ravelin::Transform3d& aTb)

    enum UpdateRule { eDone, eContinue, eInterpenetrating };

    UpdateRule handle_local_minimum(boost::shared_ptr<const Polyhedron::Vertex>& V, FeatureType& fF, boost::shared_ptr<const Polyhedron::Feature>& face,  const Polyhedron& face_poly, const Ravelin::Transform3d& vTf);

    UpdateRule update_vertex_vertex(FeatureType& fA, FeatureType& fB, Ravelin::Transform3d& aTb, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB);
    UpdateRule update_vertex_edge(FeatureType& fA, FeatureType& fB, Ravelin::Transform3d& aTb, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB);
    UpdateRule update_vertex_face(FeatureType& fA, FeatureType& fB, Ravelin::Transform3d& aTb, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB, const Polyhedron& face_poly);
    UpdateRule update_edge_edge(FeatureType& fA, FeatureType& fB, Ravelin::Transform3d& aTb, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB);
    UpdateRule update_edge_face(FeatureType& fA, FeatureType& fB, Ravelin::Transform3d& aTb, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB);
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
std::ostream& operator<<(std::ostream& out, const Polyhedron::Vertex& m);
std::ostream& operator<<(std::ostream& out, const Polyhedron::Edge& m);
std::ostream& operator<<(std::ostream& out, const Polyhedron::Face& m);

#include "Polyhedron.inl"

} // end namespace

#endif
