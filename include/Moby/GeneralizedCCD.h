/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GENERALIZED_CCD_H
#define _GENERALIZED_CCD_H

#include <list>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/sorted_pair>
#include <Moby/Log.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ThickTriangle.h>
#include <Moby/BV.h>
#include <Moby/CompGeom.h>
#include <Moby/Integrator.h>

namespace Moby {

class RigidBody;
class ArticulatedBody;
class CollisionGeometry;  

/// Implements the CollisionDetection abstract class to perform exact contact finding using abstract shapes 
class GeneralizedCCD : public CollisionDetection
{
  public:
    GeneralizedCCD();
    virtual ~GeneralizedCCD() {}
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
    GeneralizedCCD(InputIterator begin, InputIterator end);

    /// The tolerance below which subdivision does not occur
    double eps_tolerance;

  private:

    // the 3 axes
    enum AxisType { eXAxis, eYAxis, eZAxis };

    // structure for doing broad phase collision detection
    struct BoundsStruct
    {
      bool end;                   // bounds is for start or end
      CollisionGeometryPtr geom;  // the geometry
      BVPtr bv;                   // the unexpanded bounding volume
      bool operator<(const BoundsStruct& bs) const { return (!end && bs.end); } 
    };

    // structure passed to determine_TOI
    struct DStruct
    {
      CollisionGeometryPtr ga; // geometry from which point is taken
      CollisionGeometryPtr gb; // geometry against which point is tested

      // gb's poses at time t0,tf 
      Ravelin::Pose3d Pb_t0, Pb_tf;

      // transforms from ga's frame to gb's frame at t0,tf
      Ravelin::Transform3d bTa_t0, bTa_tf;

      // vertex to test for intersection (in ga's frame) 
      Point3d u_a;     

      // relative orientations (orientation of ga in gb's frame) at times t0/tf 
      Ravelin::Quatd q0, qf;      

      // BV corresponding to part of gb
      BVPtr b_BV;              
    };

    /// Pair of poses (one at time t0, the other at time tf)
    struct PosePair
    {
      Ravelin::Pose3d t0;
      Ravelin::Pose3d tf;
    };

    /// Pair of transforms (one at time t0, the other at time tf)
    struct TransformPair
    { 
      Ravelin::Transform3d t0;
      Ravelin::Transform3d tf;
    };

    /// Structure used for computing deviation in a direction
    struct DeviationCalc
    {
      Ravelin::Quatd q1;            // orientation of body at t0
      Ravelin::Quatd q2;            // orientation of body at tf
      Ravelin::Vector3d u;          // vector from com to point (in local frame)
      Ravelin::Vector3d d;          // direction to compute deviation
    };

    /// Structure used for BV processing
    struct BVProcess
    {
      BVPtr bva;         // the first BV (velocity-expanded)
      BVPtr bvb;         // the second BV (velocity-expanded)
      unsigned nexp;     // the number of expansions
      BVPtr ax;          // the original first BV (not expanded)
      BVPtr bx;          // the original second BV (not expanded)
    };

    static bool bound_u(const Ravelin::Vector3d& u, const Ravelin::Quatd& q0, const Ravelin::Quatd& qf, Ravelin::Vector3d& normal1, Ravelin::Vector3d& normal2);
    static double calc_deviation(double t, void* params);
    static double calc_min_dev(const Ravelin::Vector3d& u, const Ravelin::Vector3d& d, const Ravelin::Quatd& q1, const Ravelin::Quatd& q2, double& t);
    static double calc_max_dev(const Ravelin::Vector3d& u, const Ravelin::Vector3d& d, const Ravelin::Quatd& q1, const Ravelin::Quatd& q2, double& t);
    static std::pair<double, double> calc_deviations(const Ravelin::Vector3d& u, const Ravelin::Vector3d& d, const Ravelin::Quatd& q1, const Ravelin::Quatd& q2, double ta, double tb);
    void populate_dstruct(DStruct* ds, CollisionGeometryPtr ga, const Ravelin::Pose3d& Pa_t0, const Ravelin::Pose3d& Pa_tf, CollisionGeometryPtr gb, const Ravelin::Pose3d& Pb_t0, const Ravelin::Pose3d& Pb_tf, BVPtr b_BV) const;
    void add_rigid_body_model(RigidBodyPtr body);
    double determine_TOI(double t0, double tf, const DStruct* ds, Point3d& pt, Ravelin::Vector3d& normal) const;
    double determine_TOI_fast(double t0, double tf, const DStruct* ds, Point3d& pt, Ravelin::Vector3d& normal) const;
    BVPtr get_swept_BV(CollisionGeometryPtr g, BVPtr bv, const PosePair& v);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    static Event create_contact(double toi, CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Ravelin::Vector3d& normal);
    void check_vertices(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr ob, const std::vector<const Point3d*>& a_verts, const PosePair& a_poses, const PosePair& b_poses, double& earliest, std::vector<Event>& local_contacts) const;
    void check_geoms(double dt, CollisionGeometryPtr a, const PosePair& a_poses, CollisionGeometryPtr b, const PosePair& b_poses, std::vector<Event>& contacts); 
    void broad_phase(const std::map<CollisionGeometryPtr, PosePair>& poses, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    void sort_AABBs(const std::map<CollisionGeometryPtr, PosePair>& poses);
    void update_bounds_vector(std::vector<std::pair<double, BoundsStruct> >& bounds, const std::map<CollisionGeometryPtr, PosePair>& poses, AxisType axis);
    void build_bv_vector(const std::map<CollisionGeometryPtr, PosePair>& poses, std::vector<std::pair<double, BoundsStruct> >& bounds);
    std::map<CollisionGeometryPtr, PosePair> get_poses(const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q0, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q1);
    static bool brent(double x_lower, double x_upper, double& x, double& fx, double (*f)(double, void*), double eps, void* params);

    template <class RandomAccessIterator>
    void insertion_sort(RandomAccessIterator begin, RandomAccessIterator end);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    /// Maps of transforms for a pair of geometries
    std::map<std::pair<CollisionGeometryPtr, CollisionGeometryPtr>, TransformPair> _transform_pairs;

    /// Swept BVs computed during last call to is_contact/update_contacts()
    std::map<CollisionGeometryPtr, std::map<BVPtr, BVPtr> > _swept_BVs;

    // lock for the contact map
    pthread_mutex_t _contact_mutex;

    // lock for the swept BVs
    pthread_mutex_t _swept_BVs_mutex;

    /// AABB bounds (x-axis)
    std::vector<std::pair<double, BoundsStruct> > _x_bounds;

    /// AABB bounds (y-axis)
    std::vector<std::pair<double, BoundsStruct> > _y_bounds;

    /// AABB bounds (z-axis)
    std::vector<std::pair<double, BoundsStruct> > _z_bounds;

    /// Indicates when bounds vectors need to be rebuilt
    bool _rebuild_bounds_vecs;

    /// Maximum depth of OBB expansions (default is inf)
    unsigned _max_dexp;
}; // end class

// include inline functions
#include "GeneralizedCCD.inl"

} // end namespace

#endif

