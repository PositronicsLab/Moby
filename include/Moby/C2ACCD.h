/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _C2A_CCD_H
#define _C2A_CCD_H

#include <list>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/sorted_pair>
#include <Moby/Log.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ThickTriangle.h>
#include <Moby/BV.h>
#include <Moby/Integrator.h>

namespace Moby {

class RigidBody;
class ArticulatedBody;
class CollisionGeometry;  
class DStruct;

/// Implements the ContactFinder abstract class to perform exact contact finding using vertices against triangles 
class C2ACCD : public CollisionDetection
{
  public:
    C2ACCD();
    virtual ~C2ACCD() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual bool is_collision(double epsilon = 0.0);
    virtual bool is_contact(double dt, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q0, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q1, std::vector<Event>& contacts);
    virtual void add_collision_geometry(CollisionGeometryPtr geom);
    virtual void remove_collision_geometry(CollisionGeometryPtr geom);
    virtual void remove_all_collision_geometries();
    virtual void add_rigid_body(RigidBodyPtr body); 
    virtual void remove_rigid_body(RigidBodyPtr body); 
    virtual void add_articulated_body(ArticulatedBodyPtr abody, bool disable_adjacent);
    virtual void remove_articulated_body(ArticulatedBodyPtr abody);

    template <class InputIterator>
    C2ACCD(InputIterator begin, InputIterator end);

    /// The distance tolerance for conservative advancement 
    double eps_tolerance;

    /// The smallest advancement step fraction of the step size
    double alpha_tolerance;

  private:

    class AThickTri : public ThickTriangle
    {
      public:
        AThickTri(const Triangle& tri, double tol) : ThickTriangle(tri, tol) {}
        boost::shared_ptr<const IndexedTriArray> mesh;  // the mesh that this triangle came from
        unsigned tri_idx;             // the index of this triangle
    };

    double calc_dist(boost::shared_ptr<SSR> a, boost::shared_ptr<SSR> b, const Ravelin::Transform3d& aTb, Point3d& cpa, Point3d& cpb) const;
    static bool query_intersect_seg_tri(const LineSeg3& seg, const Triangle& tri, double& t, Point3d& p);
    void determine_contacts(CollisionGeometryPtr a, CollisionGeometryPtr b, double toc, std::vector<Event>& contacts) const;
    bool check_collision(CollisionGeometryPtr a, CollisionGeometryPtr b, std::vector<std::pair<unsigned, unsigned> >& colliding_tris) const;
    double do_CA(double step_size, CollisionGeometryPtr a, CollisionGeometryPtr b, boost::shared_ptr<SSR> ssr_a, boost::shared_ptr<SSR> ssr_b, const Ravelin::Transform3d& aTb, double dt);
    double do_CAStep(double dist, const Ravelin::Vector3d& dab, CollisionGeometryPtr a, CollisionGeometryPtr b, boost::shared_ptr<SSR> ssr_a, boost::shared_ptr<SSR> ssr_b);
    double calc_mu(double dist, const Ravelin::Vector3d& n, CollisionGeometryPtr g, boost::shared_ptr<SSR> ssr, bool positive);
    void add_rigid_body_model(RigidBodyPtr body);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    void check_vertices(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr ob, const std::vector<const Point3d*>& a_verts, const Ravelin::Transform3d bTa_t0, const Ravelin::SVelocityd& a_vel, const Ravelin::SVelocityd& b_vel, double& earliest, std::vector<Event>& local_contacts) const;
    void check_geoms(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q0, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q1, std::vector<Event>& contacts); 
    void build_BV_tree(CollisionGeometryPtr geom);
    bool split(BVPtr source, BVPtr& tgt1, BVPtr& tgt2, const Ravelin::Vector3d& axis);
    void split_tris(const Point3d& point, const Ravelin::Vector3d& normal, const IndexedTriArray& orig_mesh, const std::list<unsigned>& ofacets, std::list<unsigned>& pfacets, std::list<unsigned>& nfacets);
    void get_vertices(BVPtr bv, std::vector<const Point3d*>& vertices) const;
    static bool coplanar_tris(const Triangle& ta, const Triangle& tb);
    static bool project_and_intersect(const Triangle& ta, const Triangle& tb, std::vector<Point3d>& contact_points);
    static bool project_and_intersect(const Triangle& t, const LineSeg3& s, std::vector<Point3d>& contact_points);
    static bool project_and_intersect(const Triangle& t, const Point3d& p);
    void determine_closest_features(const Triangle& ta, const Triangle& tb, Triangle::FeatureType& fa, Triangle::FeatureType& fb, std::vector<Point3d>& contact_points) const;
    void determine_closest_tris(CollisionGeometryPtr a, CollisionGeometryPtr b, const Ravelin::Transform3d& aTb, std::vector<std::pair<Triangle, Triangle> >& closest_tris) const;
    static DynamicBodyPtr get_super_body(CollisionGeometryPtr a);
    static unsigned find_body(const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q, DynamicBodyPtr body);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    template <class InputIterator, class OutputIterator>
    OutputIterator get_vertices(const IndexedTriArray& tris, InputIterator fselect_begin, InputIterator fselect_end, OutputIterator output);

    // lock for the contact set 
    pthread_mutex_t _contact_mutex;

    // mapping from CollisionGeometry pointers to root SSRs
    std::map<CollisionGeometryPtr, boost::shared_ptr<SSR> > _root_SSRs;

    // mapping from BVs to triangles contained within
    std::map<BVPtr, std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > > _meshes;
}; // end class

// include inline functions
#include "C2ACCD.inl"

} // end namespace

#endif

