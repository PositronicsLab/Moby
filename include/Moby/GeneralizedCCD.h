/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
#include <Moby/Integrator.h>

namespace Moby {

class RigidBody;
class ArticulatedBody;
class CollisionGeometry;  
class DStruct;

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
      Ravelin::Vector3d xd;                 // the linear velocity of the rigid body
      Ravelin::Vector3d omega;              // the angular velocity of the rigid body
      BVPtr bv;                   // the unexpanded bounding volume
      bool operator<(const BoundsStruct& bs) const { return (!end && bs.end); } 
    };

    // structure passed to determine_TOI
    struct DStruct
    {
      CollisionGeometryPtr gb;    // body from which point is taken
      CollisionGeometryPtr gs;    // body against which point is tested
      Ravelin::Point3d p0;         // point in bs coordinates (before integration)
      Ravelin::Vector3d romega;     // relative angular velocity (3x1 vector, global frame)
      double thetad;        // norm of relative angular velocity
      double lvd;           // norm of relative linear velocity
      Ravelin::Vector3d k1;         // constant vector used to calculate point velocity
      Ravelin::Matrix3d k2;         // constant matrix used to calculate point velocity
      Ravelin::Twistd bs_v;      // velocity of body bs
      Ravelin::Vector3d u;          // vector from c.o.m. to point (in bb coordinates) 
      Ravelin::Quatd q0;            // initial relative orientation
      BVPtr s_BV;         // BV corresponding to part of gs
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
    static void populate_dstruct(DStruct* ds, CollisionGeometryPtr gb, CollisionGeometryPtr gs, const Ravelin::Twistd& b_vel, const Ravelin::Twistd& s_vel, BVPtr s_BV);
    void add_rigid_body_model(RigidBodyPtr body);
    double determine_TOI(double t0, double tf, const DStruct* ds, Ravelin::Point3d& pt, Ravelin::Vector3d& normal) const;
    BVPtr get_vel_exp_BV(CollisionGeometryPtr g, BVPtr bv, const Ravelin::Twistd& v);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    static Event create_contact(double toi, CollisionGeometryPtr a, CollisionGeometryPtr b, const Ravelin::Point3d& point, const Ravelin::Vector3d& normal);
    void check_vertices(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr ob, const std::vector<const Ravelin::Point3d*>& a_verts, const Ravelin::Transform3d& bTa_t0, const Ravelin::Twistd& a_vel, const Ravelin::Twistd& b_vel, double& earliest, std::vector<Event>& local_contacts) const;
    void check_geoms(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const Ravelin::Transform3d& aTb_t0, const Ravelin::Transform3d& bTa_t0, const Ravelin::Twistd& a_vel, const Ravelin::Twistd& b_vel, std::vector<Event>& contacts); 
    void broad_phase(const std::map<SingleBodyPtr, Ravelin::Twistd >& vel_map, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    void sort_AABBs(const std::map<SingleBodyPtr, Ravelin::Twistd >& vel_map);
    void update_bounds_vector(std::vector<std::pair<double, BoundsStruct> >& bounds, const std::map<SingleBodyPtr, Ravelin::Twistd >& vel_map, AxisType axis);
    void build_bv_vector(const std::map<SingleBodyPtr, Ravelin::Twistd >& vel_map, std::vector<std::pair<double, BoundsStruct> >& bounds);
    std::map<SingleBodyPtr, Ravelin::Twistd> get_velocities(const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q0, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q1, double dt) const;

    template <class RandomAccessIterator>
    void insertion_sort(RandomAccessIterator begin, RandomAccessIterator end);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    /// Velocity-expanded BVs computed during last call to is_contact/update_contacts()
    std::map<CollisionGeometryPtr, std::map<BVPtr, BVPtr> > _ve_BVs;

    // lock for the contact map
    pthread_mutex_t _contact_mutex;

    // lock for the velocity-expanded BVs
    pthread_mutex_t _ve_BVs_mutex;

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

