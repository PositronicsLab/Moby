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
    GeneralizedCCD(InputIterator begin, InputIterator end);

    /// The tolerance below which subdivision does not occur
    Real eps_tolerance;

  private:

    // the 3 axes
    enum AxisType { eXAxis, eYAxis, eZAxis };

    // structure for doing broad phase collision detection
    struct BoundsStruct
    {
      bool end;                   // bounds is for start or end
      CollisionGeometryPtr geom;  // the geometry
      Vector3 xd;                 // the linear velocity of the rigid body
      Vector3 omega;              // the angular velocity of the rigid body
      BVPtr bv;                   // the unexpanded bounding volume
      bool operator<(const BoundsStruct& bs) const { return (!end && bs.end); } 
    };

    // structure passed to determine_TOI
    struct DStruct
    {
      CollisionGeometryPtr gb;    // body from which point is taken
      CollisionGeometryPtr gs;    // body against which point is tested
      Vector3 p0;         // point in bs coordinates (before integration)
      Vector3 romega;     // relative angular velocity (3x1 vector, global frame)
      Real thetad;        // norm of relative angular velocity
      Real lvd;           // norm of relative linear velocity
      Vector3 k1;         // constant vector used to calculate point velocity
      Matrix3 k2;         // constant matrix used to calculate point velocity
      Vector3 bs_xd;      // linear velocity of body bs
      Vector3 bs_omega;   // angular velocity of body bs
      Vector3 u;          // vector from c.o.m. to point (in bb coordinates) 
      Quat q0;            // initial relative orientation
      BVPtr s_BV;         // BV corresponding to part of gs
    };

    /// Structure used for computing deviation in a direction
    struct DeviationCalc
    {
      Quat q1;            // orientation of body at t0
      Quat q2;            // orientation of body at tf
      Vector3 u;          // vector from com to point (in local frame)
      Vector3 d;          // direction to compute deviation
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

    static bool bound_u(const Vector3& u, const Quat& q0, const Quat& qf, Vector3& normal1, Vector3& normal2);
    static Real calc_deviation(Real t, void* params);
    static Real calc_min_dev(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real& t);
    static Real calc_max_dev(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real& t);
    static std::pair<Real, Real> calc_deviations(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real ta, Real tb);
    static void populate_dstruct(DStruct* ds, CollisionGeometryPtr gb, CollisionGeometryPtr gs, const Vector3& bs_lvel, const Vector3& bs_avel, const Vector3& bb_lvel, const Vector3& bb_avel, BVPtr s_BV);
    void add_rigid_body_model(RigidBodyPtr body);
    Real determine_TOI(Real t0, Real tf, const DStruct* ds, Vector3& pt, Vector3& normal) const;
    BVPtr get_vel_exp_BV(CollisionGeometryPtr g, BVPtr bv, const Vector3& lv, const Vector3& av);
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    static Event create_contact(Real toi, CollisionGeometryPtr a, CollisionGeometryPtr b, const Vector3& point, const Vector3& normal);
    void check_vertices(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr ob, const std::vector<const Vector3*>& a_verts, const Matrix4& bTa_t0, const std::pair<Vector3, Vector3>& a_vel, const std::pair<Vector3, Vector3>& b_vel, Real& earliest, std::vector<Event>& local_contacts) const;
    void check_geoms(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const Matrix4& aTb_t0, const Matrix4& bTa_t0, const std::pair<Vector3, Vector3>& a_vel, const std::pair<Vector3, Vector3>& b_vel, std::vector<Event>& contacts); 
    void broad_phase(const std::map<SingleBodyPtr, std::pair<Vector3, Vector3> >& vel_map, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check);
    void sort_AABBs(const std::map<SingleBodyPtr, std::pair<Vector3, Vector3> >& vel_map);
    void update_bounds_vector(std::vector<std::pair<Real, BoundsStruct> >& bounds, const std::map<SingleBodyPtr, std::pair<Vector3, Vector3> >& vel_map, AxisType axis);
    void build_bv_vector(const std::map<SingleBodyPtr, std::pair<Vector3, Vector3> >& vel_map, std::vector<std::pair<Real, BoundsStruct> >& bounds);
    std::map<SingleBodyPtr, std::pair<Vector3, Vector3> > get_velocities(const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q0, const std::vector<std::pair<DynamicBodyPtr, VectorN> >& q1, Real dt) const;

    template <class RandomAccessIterator>
    void insertion_sort(RandomAccessIterator begin, RandomAccessIterator end);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;

    /// Velocity-expanded BVs computed during last call to is_contact/update_contacts()
    std::map<CollisionGeometryPtr, std::map<BVPtr, BVPtr> > _ve_BVs;

    // lock for the contact map
    pthread_mutex_t _contact_mutex;

    // lock for the velocity-expanded BVs
    pthread_mutex_t _ve_BVs_mutex;

    /// AABB bounds (x-axis)
    std::vector<std::pair<Real, BoundsStruct> > _x_bounds;

    /// AABB bounds (y-axis)
    std::vector<std::pair<Real, BoundsStruct> > _y_bounds;

    /// AABB bounds (z-axis)
    std::vector<std::pair<Real, BoundsStruct> > _z_bounds;

    /// Indicates when bounds vectors need to be rebuilt
    bool _rebuild_bounds_vecs;

    /// Maximum depth of OBB expansions (default is inf)
    unsigned _max_dexp;
}; // end class

// include inline functions
#include "GeneralizedCCD.inl"

} // end namespace

#endif

