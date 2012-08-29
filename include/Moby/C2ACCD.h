/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual bool is_collision(Real epsilon = 0.0);
    virtual bool is_contact(Real dt, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, std::vector<Event>& contacts);
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
    Real eps_tolerance;

    /// The smallest advancement step fraction of the step size
    Real alpha_tolerance;

  private:

    class AThickTri : public ThickTriangle
    {
      public:
        AThickTri(const Triangle& tri, Real tol) : ThickTriangle(tri, tol) {}
        boost::shared_ptr<const IndexedTriArray> mesh;  // the mesh that this triangle came from
        unsigned tri_idx;             // the index of this triangle
    };

    Real calc_dist(boost::shared_ptr<SSR> a, boost::shared_ptr<SSR> b, const Matrix4& aTb, Vector3& cpa, Vector3& cpb) const;
    static bool query_intersect_seg_tri(const LineSeg3& seg, const Triangle& tri, Real& t, Vector3& p);
    void determine_contacts(CollisionGeometryPtr a, CollisionGeometryPtr b, Real toc, std::vector<Event>& contacts) const;
    bool check_collision(CollisionGeometryPtr a, CollisionGeometryPtr b, std::vector<std::pair<unsigned, unsigned> >& colliding_tris) const;
    Real do_CA(Real step_size, CollisionGeometryPtr a, CollisionGeometryPtr b, boost::shared_ptr<SSR> ssr_a, boost::shared_ptr<SSR> ssr_b, const Matrix4& aTb, Real dt);
    Real do_CAStep(Real dist, const Vector3& dab, CollisionGeometryPtr a, CollisionGeometryPtr b, boost::shared_ptr<SSR> ssr_a, boost::shared_ptr<SSR> ssr_b);
    Real calc_mu(Real dist, const Vector3& n, CollisionGeometryPtr g, boost::shared_ptr<SSR> ssr, bool positive);
    void add_rigid_body_model(RigidBodyPtr body);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    void check_vertices(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr ob, const std::vector<const Vector3*>& a_verts, const Matrix4& bTa_t0, const std::pair<Vector3, Vector3>& a_vel, const std::pair<Vector3, Vector3>& b_vel, Real& earliest, std::vector<Event>& local_contacts) const;
    void check_geoms(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, std::vector<Event>& contacts); 
    void build_BV_tree(CollisionGeometryPtr geom);
    bool split(BVPtr source, BVPtr& tgt1, BVPtr& tgt2, const Vector3& axis, bool deformable);
    void split_tris(const Vector3& point, const Vector3& normal, const IndexedTriArray& orig_mesh, const std::list<unsigned>& ofacets, std::list<unsigned>& pfacets, std::list<unsigned>& nfacets);
    void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices) const;
    static bool coplanar_tris(const Triangle& ta, const Triangle& tb);
    static bool project_and_intersect(const Triangle& ta, const Triangle& tb, std::vector<Vector3>& contact_points);
    static bool project_and_intersect(const Triangle& t, const LineSeg3& s, std::vector<Vector3>& contact_points);
    static bool project_and_intersect(const Triangle& t, const Vector3& p);
    void determine_closest_features(const Triangle& ta, const Triangle& tb, Triangle::FeatureType& fa, Triangle::FeatureType& fb, std::vector<Vector3>& contact_points) const;
    void determine_closest_tris(CollisionGeometryPtr a, CollisionGeometryPtr b, const Matrix4& aTb, std::vector<std::pair<Triangle, Triangle> >& closest_tris) const;
    static DynamicBodyPtr get_super_body(CollisionGeometryPtr a);
    static unsigned find_body(const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q, DynamicBodyPtr body);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

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

