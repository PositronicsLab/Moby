/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MESH_DCD_H
#define _MESH_DCD_H

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
class MeshDCD : public CollisionDetection
{
  public:
    MeshDCD();
    virtual ~MeshDCD() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual bool is_collision(Real epsilon = 0.0);
    virtual bool is_contact(Real dt, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, std::vector<Event>& contact);
    virtual void add_collision_geometry(CollisionGeometryPtr geom);
    virtual void remove_collision_geometry(CollisionGeometryPtr geom);
    virtual void remove_all_collision_geometries();
    virtual void add_deformable_body(DeformableBodyPtr body);
    virtual void add_rigid_body(RigidBodyPtr body); 
    virtual void remove_rigid_body(RigidBodyPtr body); 
    virtual void remove_deformable_body(DeformableBodyPtr body);
    virtual void add_articulated_body(ArticulatedBodyPtr abody, bool disable_adjacent);
    virtual void remove_articulated_body(ArticulatedBodyPtr abody);

    template <class InputIterator>
    MeshDCD(InputIterator begin, InputIterator end);

    /// The epsilon tolerance [0,1]
    Real eps_tolerance;

    /// The intersection tolerance
    Real isect_tolerance;

    Real intersect_rects(const Vector3& normal, const Vector3 r1[4], const Vector3 r2[4], Vector3& isect1, Vector3& isect2);
    Real calc_first_isect(const Triangle& t, const LineSeg3& s1, const LineSeg3& s2, Vector3& p1, Vector3& p2); 
  private:
    Real intersect_rect(const Vector3& normal, const Vector3& axis1, const Vector3& axis2, const LineSeg3& rs1, const LineSeg3& rs2, const LineSeg3& s, Vector3& isect1, Vector3& isect2);
    static unsigned determine_cubic_roots(Real a, Real b, Real c, Real x[3]);
    void add_rigid_body_model(RigidBodyPtr body);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    static Event create_contact(Real toi, CollisionGeometryPtr a, CollisionGeometryPtr b, const Vector3& point, const Vector3& vpoint, const Triangle& t);
    void check_geoms(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, std::vector<Event>& contacts); 
    void check_geom(Real dt, CollisionGeometryPtr cg, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, std::vector<Event>& contacts); 
    void broad_phase(std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    void determine_contacts_rigid(CollisionGeometryPtr a, CollisionGeometryPtr b, Real t, Real dt, std::vector<Event>& contacts);
    void determine_contacts_deformable(CollisionGeometryPtr a, CollisionGeometryPtr b, Real t, Real dt, std::vector<Event>& contacts);
    void determine_contacts_rigid_deformable(CollisionGeometryPtr a, CollisionGeometryPtr b, Real t, Real dt, std::vector<Event>& contacts);
    void determine_contacts_deformable_rigid(CollisionGeometryPtr a, CollisionGeometryPtr b, Real t, Real dt, std::vector<Event>& contacts);
    Real calc_first_isect(const Triangle& t, const LineSeg3& seg, Vector3& p);
    static Real calc_first_isect(const Vector3& p, const Vector3& pdot, const Triangle& t, const Vector3& adot, const Vector3& bdot, const Vector3& cdot, Real dt);
    static Real calc_param(const LineSeg3& seg, const Vector3& p);
    bool is_collision(CollisionGeometryPtr a, CollisionGeometryPtr b);
    bool is_collision(CollisionGeometryPtr cg);
    static DynamicBodyPtr get_super_body(CollisionGeometryPtr a);
    static unsigned find_body(const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q, DynamicBodyPtr body);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    /// Indicates when bounds vectors need to be rebuilt
    bool _rebuild_bounds_vecs;

}; // end class

// include inline functions
#include "MeshDCD.inl"

} // end namespace

#endif

