/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CCD_H
#define _CCD_H

#include <list>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Ravelin/sorted_pair>
#include <Moby/Log.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/TorusPrimitive.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/CollisionDetection.h>
#include <Moby/BV.h>
#include <Moby/GJK.h>
#include <Moby/Polyhedron.h>

namespace Moby {

class RigidBody;
class ArticulatedBody;
class CollisionGeometry;

/// Implements the CollisionDetection abstract class to perform exact contact finding using abstract shapes
class CCD : public CollisionDetection
{
  public:
    CCD();
    virtual ~CCD() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi);
    virtual void find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      find_contacts(cgA, cgB, std::back_inserter(contacts), TOL);
    }
    static unsigned constrain_unsigned(int ii, int maxi){
     return (unsigned) std::min(std::max(ii,0),maxi);
    }

    template <class OutputIterator>
    OutputIterator find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL = NEAR_ZERO);

    /// Pairs of collision geometries that aren't checked for contact/collision
    /**
     * \note collisions between geometries for two disabled bodies and
     *       collisions between geometries for a single body are automatically
     *       not checked and do not need to be added to this set.
     */
    std::set<Ravelin::sorted_pair<CollisionGeometryPtr> > disabled_pairs;

  protected:
    virtual double calc_next_CA_Euler_step(const PairwiseDistInfo& pdi) { return calc_next_CA_Euler_step_generic(pdi); }

  private:
    // the 3 axes
    enum AxisType { eXAxis, eYAxis, eZAxis };

    double calc_next_CA_Euler_step_polyhedron_plane(boost::shared_ptr<PolyhedralPrimitive> p, const Ravelin::SVelocityd& rv, boost::shared_ptr<const Ravelin::Pose3d> P, const Ravelin::Vector3d& normal, double offset);
    double calc_next_CA_Euler_step_polyhedron_polyhedron(boost::shared_ptr<PolyhedralPrimitive> pA, boost::shared_ptr<PolyhedralPrimitive> pB, boost::shared_ptr<const Ravelin::Pose3d> poseA, boost::shared_ptr<const Ravelin::Pose3d> poseB, const Ravelin::SVelocityd& rvA, const Ravelin::SVelocityd& rvB, const Ravelin::Vector3d& n0, double offset);
    virtual double calc_next_CA_Euler_step_generic(const PairwiseDistInfo& pdi);
    virtual double calc_CA_Euler_step_generic(const PairwiseDistInfo& pdi);
    virtual double calc_CA_Euler_step_sphere(const PairwiseDistInfo& pdi);
    static double calc_max_dist(ArticulatedBodyPtr ab, RigidBodyPtr rb, const Ravelin::Vector3d& n, double rmax);
    static double calc_max_dist(RigidBodyPtr rb, const Ravelin::Vector3d& n, double rmax);
    static double calc_max_step(RigidBodyPtr rbA, RigidBodyPtr rbB, const Ravelin::Vector3d& n, double rmaxA, double rmaxB, double dist);

    // structure for doing broad phase collision detection
    struct BoundsStruct
    {
      bool end;                   // bounds is for start or end
      CollisionGeometryPtr geom;  // the geometry
      BVPtr bv;                   // the unexpanded bounding volume
      bool operator<(const BoundsStruct& bs) const { return (!end && bs.end); }
    };

    // gets the distance on farthest points
    std::map<CollisionGeometryPtr, double> _rmax;

    // see whether the bounds vectors need to be rebuilt
    bool _rebuild_bounds_vecs;

    // the bounding spheres
    std::map<CollisionGeometryPtr, BVPtr> _bounding_spheres;

    /// AABB bounds (x-axis)
    std::vector<std::pair<double, BoundsStruct> > _x_bounds;

    /// AABB bounds (y-axis)
    std::vector<std::pair<double, BoundsStruct> > _y_bounds;

    /// AABB bounds (z-axis)
    std::vector<std::pair<double, BoundsStruct> > _z_bounds;

    /// Swept BVs computed during last call to is_contact/update_contacts()
    std::map<CollisionGeometryPtr, BVPtr> _swept_BVs;

    /// Minimum observed distance between two bodies (to make conservative advancement faster in face of numerical error)
    std::map<Ravelin::sorted_pair<CollisionGeometryPtr>, double> _min_dist_observed;

    static BVPtr construct_bounding_sphere(CollisionGeometryPtr cg);
    void sort_AABBs(const std::vector<RigidBodyPtr>& rigid_bodies, double dt);
    void update_bounds_vector(std::vector<std::pair<double, BoundsStruct> >& bounds, AxisType axis, double dt, bool recreate_bvs);
    void build_bv_vector(const std::vector<RigidBodyPtr>& rigid_bodies, std::vector<std::pair<double, BoundsStruct> >& bounds);
    BVPtr get_swept_BV(CollisionGeometryPtr geom, BVPtr bv, double dt);

    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);

    template <class OutputIterator>
    OutputIterator find_contacts_polyhedron_polyhedron(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    template <class OutputIterator>
    OutputIterator find_contacts_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_plane_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_torus_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_sphere_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_cylinder_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_heightmap_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_sphere_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_sphere_heightmap(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_convex_heightmap(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class OutputIterator>
    OutputIterator find_contacts_box_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL);

    template <class RandomAccessIterator>
    void insertion_sort(RandomAccessIterator begin, RandomAccessIterator end);

    template <class OutputIterator>
    OutputIterator find_contacts_vertex_vertex(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> v1, boost::shared_ptr<Polyhedron::Vertex> v2, double signed_dist, OutputIterator output_begin);

    template <class OutputIterator>
    OutputIterator find_contacts_vertex_edge(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> v, boost::shared_ptr<Polyhedron::Edge> e, double signed_dist, OutputIterator output_begin);

    template <class OutputIterator>
    OutputIterator find_contacts_vertex_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> vA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin);

    template <class OutputIterator>
    OutputIterator find_contacts_edge_edge(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> e1, boost::shared_ptr<Polyhedron::Edge> e2, double signed_dist, OutputIterator output_begin);

    template <class OutputIterator>
    OutputIterator find_contacts_edge_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> eA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin);

    template <class OutputIterator>
    OutputIterator find_contacts_face_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Face> fA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin);

}; // end class

#include "CCD.inl"

} // end namespace

#endif

