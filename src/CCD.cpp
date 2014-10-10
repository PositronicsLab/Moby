/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <stack>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <Moby/UnilateralConstraint.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/XMLTree.h>
#include <Moby/SSL.h>
#include <Moby/BoundingSphere.h>
#include <Moby/DummyBV.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/GaussianMixture.h>
#include <Moby/CSG.h>
#include <Moby/CCD.h>

using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using boost::shared_ptr;
using boost::tuple;
using boost::make_tuple;
using std::cerr;
using std::endl;
using std::set;
using std::queue;
using std::map;
using std::list;
using std::vector;
using std::priority_queue;
using std::pair;
using std::make_pair;
using namespace Ravelin;
using namespace Moby;

/// Constructs a collision detector with default tolerances
CCD::CCD()
{
}

/// Computes a conservative advancement step between two collision geometries
double CCD::calc_CA_step(const PairwiseDistInfo& pdi)
{
  double maxt = std::numeric_limits<double>::max();

  // get geometries, distance, and closest points
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;
  const Point3d& pA = pdi.pa;
  const Point3d& pB = pdi.pb;

  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());
  FILE_LOG(LOG_COLDET) << "rigid body A: " << rbA->id << "  rigid body B: " << rbB->id << std::endl;

  // if the distance is (essentially) zero, quit now
  if (pdi.dist <= 0.0)
  {
    FILE_LOG(LOG_COLDET) << "reported distance is: " << pdi.dist << std::endl;
    return 0.0;
  }

  // get the direction of the vector from body B to body A
  Vector3d d0 = Pose3d::transform_point(GLOBAL, pA) -
                Pose3d::transform_point(GLOBAL, pB);
  double d0_norm = d0.norm();
  FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
  FILE_LOG(LOG_COLDET) << "reported distance is: " << pdi.dist << std::endl;

  // get the direction of the vector (from body B to body A)
  Vector3d n0 = d0/d0_norm;
  Vector3d nA = Pose3d::transform_vector(rbA->get_pose(), n0);
  Vector3d nB = Pose3d::transform_vector(rbB->get_pose(), n0);

  // compute the distance that body A can move toward body B
  double dist_per_tA = calc_max_dist_per_t(rbA, -nA, _rmax[cgA]);

  // compute the distance that body B can move toward body A
  double dist_per_tB = calc_max_dist_per_t(rbB, nB, _rmax[cgB]);

  // compute the total distance
  double total_dist_per_t = dist_per_tA + dist_per_tB;
  if (total_dist_per_t < 0.0)
    total_dist_per_t = 0.0;

  FILE_LOG(LOG_COLDET) << "  distance: " << pdi.dist << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tA: " << dist_per_tA << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tB: " << dist_per_tB << std::endl;

  // compute the maximum safe step
  maxt = std::min(maxt, pdi.dist/total_dist_per_t);

  FILE_LOG(LOG_COLDET) << "  maxt: " << maxt << std::endl;

  // return the maximum safe step
  return maxt;
}

/// Computes a conservative advancement step between two collision geometries assuming that velocity is constant over the interval
double CCD::calc_CA_Euler_step(const PairwiseDistInfo& pdi)
{
  double maxt = std::numeric_limits<double>::max();

  // get geometries, distance, and closest points
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;
  const Point3d& pA = pdi.pa;
  const Point3d& pB = pdi.pb;

  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());
  FILE_LOG(LOG_COLDET) << "rigid body A: " << rbA->id << "  rigid body B: " << rbB->id << std::endl;

  // if the distance is (essentially) zero, quit now
  if (pdi.dist <= 0.0)
  {
    FILE_LOG(LOG_COLDET) << "reported distance is: " << pdi.dist << std::endl;
    return 0.0;
  }

  // get the direction of the vector from body B to body A
  Vector3d d0 = Pose3d::transform_point(GLOBAL, pA) -
                Pose3d::transform_point(GLOBAL, pB);
  double d0_norm = d0.norm();
  FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
  FILE_LOG(LOG_COLDET) << "reported distance is: " << pdi.dist << std::endl;

  // get the direction of the vector (from body B to body A)
  Vector3d n0 = d0/d0_norm;
  Vector3d nA = Pose3d::transform_vector(rbA->get_pose(), n0);
  Vector3d nB = Pose3d::transform_vector(rbB->get_pose(), n0);

  // compute the distance that body A can move toward body B
  double dist_per_tA = calc_max_velocity(rbA, -nA, _rmax[cgA]);

  // compute the distance that body B can move toward body A
  double dist_per_tB = calc_max_velocity(rbB, nB, _rmax[cgB]);

  // compute the total distance
  double total_dist_per_t = dist_per_tA + dist_per_tB;
  if (total_dist_per_t < 0.0)
    total_dist_per_t = 0.0;

  FILE_LOG(LOG_COLDET) << "  distance: " << pdi.dist << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tA: " << dist_per_tA << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tB: " << dist_per_tB << std::endl;

  // compute the maximum safe step
  maxt = std::min(maxt, pdi.dist/total_dist_per_t);

  FILE_LOG(LOG_COLDET) << "  maxt: " << maxt << std::endl;

  // return the maximum safe step
  return maxt;
}

/// Computes the maximum velocity along a particular direction (n)
double CCD::calc_max_velocity(RigidBodyPtr rb, const Vector3d& n, double rmax)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the velocities at t0
  const SVelocityd& v0 = Pose3d::transform(rb->get_pose(), rb->get_velocity());
  Vector3d xd0 = v0.get_linear();
  Vector3d w0 = v0.get_angular();
  return n.dot(xd0) + w0.norm()*rmax;
}

/// Solves the LP that maximizes <n, v + w x r>
double CCD::calc_max_dist_per_t(RigidBodyPtr rb, const Vector3d& n, double rlen)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the velocities at t0
  const SVelocityd& v0 = Pose3d::transform(rb->get_pose(), rb->get_velocity());
  Vector3d xd0 = v0.get_linear();
  Vector3d w0 = v0.get_angular();

  // get the bounds on velocity 
  const SVelocityd& v0p = rb->get_vel_upper_bounds();
  const SVelocityd& v0n = rb->get_vel_lower_bounds();
  assert(v0p.pose == v0.pose);
  assert(v0n.pose == v0.pose);

  // break velocity bounds into linear and angular components
  Vector3d xdp = v0p.get_linear();
  Vector3d xdn = v0n.get_linear();
  Vector3d omegap = v0p.get_angular();
  Vector3d omegan = v0n.get_angular();

  // (from Mirtich [1996], p. 131)
  double dist = n.dot(xd0) + w0.norm()*rlen;
  double wmaxx = std::max(std::fabs(omegap[X]), std::fabs(omegan[X]));
  double wmaxy = std::max(std::fabs(omegap[Y]), std::fabs(omegan[Y]));
  double wmaxz = std::max(std::fabs(omegap[Z]), std::fabs(omegan[Z]));
  dist += (n[X] < 0.0) ? xdn[X]*n[X] : xdp[X]*n[X];
  dist += (n[Y] < 0.0) ? xdn[Y]*n[Y] : xdp[Y]*n[Y];
  dist += (n[Z] < 0.0) ? xdn[Z]*n[Z] : xdp[Z]*n[Z];
  dist += rlen*(wmaxx*std::sqrt(1-n[X]*n[X]) + wmaxy*std::sqrt(1-n[Y]*n[Y]) +
                wmaxz*std::sqrt(1-n[Z]*n[Z]));

  // comparison using LP
  // LP is maximize n'(xd + w x r) = n'xd + (n x -r)'w
  #ifndef NDEBUG
  MatrixNd A(0,6);
  VectorNd b(0);
  static VectorNd c(6), l(6), u(6), x(6);
  c[0] = n[X];       c[1] = n[Y];       c[2] = n[Z];
  c[3] = n[X]*rlen; c[4] = n[Y]*rlen; c[5] = n[Z]*rlen;
  l[0] = xdn[X];     l[1] = xdn[Y];     l[2] = xdn[Z];
  u[0] = xdp[X];     u[1] = xdp[Y];     u[2] = xdp[Z];
  l[3] = omegan[X];  l[4] = omegan[Y];  l[5] = omegan[Z];
  u[3] = omegap[X];  u[4] = omegap[Y];  u[5] = omegap[Z];

  // solve the LP
  bool sol = lp_seidel(A, b, c, l, u, x);
  if (sol)
  {
    double dist2 = c.dot(x);
    if (std::fabs(dist2 - dist) > std::sqrt(NEAR_ZERO))
      FILE_LOG(LOG_COLDET) << "Difference detected: " << dist << " vs. " << dist2 << std::endl;
    dist = std::max(dist,dist2);
  }

  #endif

  return std::fabs(dist);
}

/// Implements Base::load_from_xml()
void CCD::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "CCD") == 0);
}

/// Implements Base::save_to_xml()
/**
 * \note neither the contact cache nor the pairs currently in collision are
 *       saved
 */
void CCD::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // (re)set the node name
  node->name = "CCD";
}

/****************************************************************************
 Methods for Drumwright-Shell algorithm begin
****************************************************************************/

/// Creates a contact constraint given the bare-minimum info
UnilateralConstraint CCD::create_contact(CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Vector3d& normal, double violation)
{
  UnilateralConstraint e;
  e.constraint_type = UnilateralConstraint::eContact;
  e.contact_point = point;
  e.contact_normal = normal;
  e.contact_geom1 = a;
  e.contact_geom2 = b;
  e.signed_violation = violation;

  // check for valid normal here
  assert(std::fabs(e.contact_normal.norm() - (double) 1.0) < NEAR_ZERO);

  // make the body first that comes first alphabetically
  if (LOGGING(LOG_COLDET))
  {
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();
    if (sb2->id < sb1->id)
    {
      std::swap(e.contact_geom1, e.contact_geom2);
      e.contact_normal = -e.contact_normal;
    }
  }

  // transform contact point and normal to global frame
  e.contact_point = Pose3d::transform_point(GLOBAL, e.contact_point);
  e.contact_normal = Pose3d::transform_vector(GLOBAL, e.contact_normal);

  return e;
}

/****************************************************************************
 Methods for static geometry intersection testing end
****************************************************************************/

/****************************************************************************
 Methods for broad phase begin
****************************************************************************/
void CCD::broad_phase(double dt, const vector<DynamicBodyPtr>& bodies, vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
{
  FILE_LOG(LOG_COLDET) << "CCD::broad_phase() entered" << std::endl;

  // clear the swept BVs
  _swept_BVs.clear();

  // get the set of rigid bodies
  vector<RigidBodyPtr> rbs;
  for (unsigned i=0; i< bodies.size(); i++)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
    if (ab)
      rbs.insert(rbs.end(), ab->get_links().begin(), ab->get_links().end());
    else
      rbs.push_back(dynamic_pointer_cast<RigidBody>(bodies[i]));
  }

  // look to see whether the bounds vector needs to be rebuilt
  unsigned count = 0;
  for (unsigned i=0; i< rbs.size(); i++)
    count += rbs[i]->geometries.size();
  if (count != _bounding_spheres.size())
  {
    // clear the map of bounding spheres
    _bounding_spheres.clear();

    // indicate the bounding vectors need to be rebuilt
    _rebuild_bounds_vecs = true;

    for (unsigned j=0; j< rbs.size(); j++)
    BOOST_FOREACH(CollisionGeometryPtr i, rbs[j]->geometries)
    {
      // get farthest distance on each geometry while we're at it
      _rmax[i] = i->get_farthest_point_distance();

      // get the primitive for the geometry
      PrimitivePtr p = i->get_geometry();

      // construct a bounding sphere for the geometry
      BVPtr bv = construct_bounding_sphere(i);

      // store the bounding sphere
      _bounding_spheres[i] = bv;
    }
  }
  else
  {
    // update any bounding spheres that need to be rebuild
    for (std::map<CollisionGeometryPtr, BVPtr>::iterator i = _bounding_spheres.begin(); i != _bounding_spheres.end(); i++)
    {
      RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(i->first->get_single_body());
      if (rb->is_enabled())
      {
        BVPtr old_bv = i->second;
        i->second = construct_bounding_sphere(i->first);
        for (unsigned j=0; j< _x_bounds.size(); j++)
        {
          if (_x_bounds[j].second.bv == old_bv) _x_bounds[j].second.bv = i->second;
          if (_y_bounds[j].second.bv == old_bv) _y_bounds[j].second.bv = i->second;
          if (_z_bounds[j].second.bv == old_bv) _z_bounds[j].second.bv = i->second;
        }
      }
    }
  }

  // clear the vector of pairs to check
  to_check.clear();

  // sort the AABBs
  sort_AABBs(rbs, dt);

  // store how many overlaps we have for pairs
  map<sorted_pair<CollisionGeometryPtr>, unsigned> overlaps;

  // set of active bounds
  set<CollisionGeometryPtr> active_bounds;

  // scan through the x-bounds
  for (unsigned i=0; i< _x_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_x_bounds[i].second.end)
    {
      assert(active_bounds.find(_x_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_x_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _x_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_x_bounds[i].second.geom);
    }
  }

  // scan through the y-bounds
  for (unsigned i=0; i< _y_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_y_bounds[i].second.end)
    {
      assert(active_bounds.find(_y_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_y_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _y_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_y_bounds[i].second.geom);
    }
  }

  // scan through the z-bounds
  for (unsigned i=0; i< _z_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_z_bounds[i].second.end)
    {
      assert(active_bounds.find(_z_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_z_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _z_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_z_bounds[i].second.geom);
    }
  }

  // now setup pairs to check
  for (map<sorted_pair<CollisionGeometryPtr>, unsigned>::const_iterator i = overlaps.begin(); i != overlaps.end(); i++)
  {
    FILE_LOG(LOG_COLDET) << i->second << " overlaps between " << i->first.first << " (" << i->first.first->get_single_body()->id << ") and " << i->first.second << " (" << i->first.second->get_single_body()->id << ")" << std::endl;

    // only consider the pair if they overlap in all three dimensions
    if (i->second < 3)
      continue;

    // if the pair is disabled, continue looping
    if (this->disabled_pairs.find(i->first) != this->disabled_pairs.end())
      continue;

    // get the rigid bodies corresponding to the geometries
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(i->first.first->get_single_body());
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(i->first.second->get_single_body());

    // don't check pairs from the same rigid body
    if (rb1 == rb2)
      continue;

    // if both rigid bodies are disabled, don't check
    if (!rb1->is_enabled() && !rb2->is_enabled())
      continue;

    // if we're here, we have a candidate for the narrow phase
    to_check.push_back(make_pair(i->first.first, i->first.second));
    FILE_LOG(LOG_COLDET) << "  ... checking pair" << std::endl;
  }

  FILE_LOG(LOG_COLDET) << "CCD::broad_phase() exited" << std::endl;
}

/// Gets the swept BV, creating it if necessary
BVPtr CCD::get_swept_BV(CollisionGeometryPtr cg, BVPtr bv, double dt)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // verify that the map for the geometry has already been setup
  map<CollisionGeometryPtr, BVPtr>::iterator vi;
  vi = _swept_BVs.find(cg);
  if (vi != _swept_BVs.end())
    return vi->second;

  // get the rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(cg->get_single_body());

  // get the current velocity and the velocity limits
  const SVelocityd& v = rb->get_velocity();
  Vector3d v_lo = rb->get_vel_lower_bounds().get_linear();
  Vector3d v_hi = rb->get_vel_upper_bounds().get_linear();

  // compute the swept BV
  BVPtr swept_bv = bv->calc_swept_BV(cg, v*dt);
  FILE_LOG(LOG_BV) << "new BV: " << swept_bv << std::endl;

  // attempt to get the swept BV as a SSL
  shared_ptr<SSL> ssl = dynamic_pointer_cast<SSL>(swept_bv);
  if (ssl)
  {
    // update the radius
    ssl->radius += dt*std::max(std::fabs(v_lo[X]), std::fabs(v_hi[X]));
    ssl->radius += dt*std::max(std::fabs(v_lo[Y]), std::fabs(v_hi[Y]));
    ssl->radius += dt*std::max(std::fabs(v_lo[Z]), std::fabs(v_hi[Z]));
  }

  // store the bounding volume
  _swept_BVs[cg] = swept_bv;

  return swept_bv;
}

void CCD::sort_AABBs(const vector<RigidBodyPtr>& rigid_bodies, double dt)
{
  // rebuild the vector of bounds
  if (_rebuild_bounds_vecs)
  {
    build_bv_vector(rigid_bodies, _x_bounds);
    _y_bounds = _x_bounds;
    _z_bounds = _x_bounds;
  }

  // update bounds vectors
  update_bounds_vector(_x_bounds, eXAxis, dt, true);
  update_bounds_vector(_y_bounds, eYAxis, dt, false);
  update_bounds_vector(_z_bounds, eZAxis, dt, false);

  // if geometry was added or removed, do standard sorts of vectors
  if (_rebuild_bounds_vecs)
  {
    std::sort(_x_bounds.begin(), _x_bounds.end());
    std::sort(_y_bounds.begin(), _y_bounds.end());
    std::sort(_z_bounds.begin(), _z_bounds.end());

    // now indicate that bounds vectors have been (re)built
    _rebuild_bounds_vecs = false;
  }
  else
  {
    // bounds are nearly sorted; do insertion sort
    insertion_sort(_x_bounds.begin(), _x_bounds.end());
    insertion_sort(_y_bounds.begin(), _y_bounds.end());
    insertion_sort(_z_bounds.begin(), _z_bounds.end());
  }
}

void CCD::update_bounds_vector(vector<pair<double, BoundsStruct> >& bounds, AxisType axis, double dt, bool recreate_bvs)
{
  const unsigned X = 0, Y = 1, Z = 2;

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() entered (axis=" << axis << ")" << std::endl;

  // iterate over bounds vector
  for (unsigned i=0; i< bounds.size(); i++)
  {
    // get the bounding volume and collision geometry
    BVPtr bv = bounds[i].second.bv;
    CollisionGeometryPtr geom = bounds[i].second.geom;

    // recreate the BV if necessary
    if (recreate_bvs)
      _swept_BVs.erase(geom);

    // get the swept bounding volume (should be defined in global frame)
    BVPtr swept_bv = get_swept_BV(geom, bv, dt);
    assert(swept_bv->get_relative_pose() == GLOBAL);

    // get the bound for the bounding volume
    Point3d bound = (bounds[i].second.end) ? swept_bv->get_upper_bounds() : swept_bv->get_lower_bounds();
    FILE_LOG(LOG_COLDET) << "  updating collision geometry: " << geom << "  rigid body: " << geom->get_single_body()->id << std::endl;

    // update the bounds for the given axis
    switch (axis)
    {
      case eXAxis:
        bounds[i].first = bound[X];
        break;

      case eYAxis:
        bounds[i].first = bound[Y];
        break;

      case eZAxis:
        bounds[i].first = bound[Z];
        break;

      default:
        assert(false);
    }

    if (bounds[i].second.end)
      FILE_LOG(LOG_COLDET) << "    upper bound: " << bounds[i].first << std::endl;
    if (!bounds[i].second.end)
      FILE_LOG(LOG_COLDET) << "    lower bound: " << bounds[i].first << std::endl;
  }

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() exited" << std::endl;
}

void CCD::build_bv_vector(const vector<RigidBodyPtr>& rigid_bodies, vector<pair<double, BoundsStruct> >& bounds)
{
  const double INF = std::numeric_limits<double>::max();

  // clear the vector
  bounds.clear();

  // iterate over all collision geometries
  for (unsigned j=0; j< rigid_bodies.size(); j++)
    BOOST_FOREACH(CollisionGeometryPtr i, rigid_bodies[j]->geometries)
    {
      // get the primitive for the geometry
      PrimitivePtr p = i->get_geometry();

      // get the bounding sphere
      BVPtr bv = _bounding_spheres.find(i)->second;

      // setup the bounds structure
      BoundsStruct bs;
      bs.end = false;
      bs.geom = i;
      bs.bv = bv;

      // add the lower bound
      bounds.push_back(make_pair(-INF, bs));

      // modify the bounds structure to indicate the end bound
      bs.end = true;
      bounds.push_back(make_pair(INF, bs));
    }
}

/// Constructs a bounding sphere for a given primitive type
BVPtr CCD::construct_bounding_sphere(CollisionGeometryPtr cg)
{
  // create the bounding sphere
  shared_ptr<BoundingSphere> sph(new BoundingSphere);

  // get the primitive type
  PrimitivePtr p = cg->get_geometry();

  // setup the relative pose
  shared_ptr<Pose3d> pose(new Pose3d(*p->get_pose()));
  assert(pose->rpose == GLOBAL);
  pose->rpose = cg->get_pose();
  sph->center.set_zero(pose);

  // now convert the center to the GLOBAL frame
  sph->center = Pose3d::transform_point(GLOBAL, sph->center);

  // look for sphere primitive (easiest case)
  shared_ptr<SpherePrimitive> sph_p = dynamic_pointer_cast<SpherePrimitive>(p);
  if (sph_p)
  {
    sph->radius = sph_p->get_radius();
    return sph;
  }

  // look for box primitive (also an easy case)
  shared_ptr<BoxPrimitive> box_p = dynamic_pointer_cast<BoxPrimitive>(p);
  if (box_p)
  {
    sph->radius = Origin3d(box_p->get_x_len()/2.0, box_p->get_y_len()/2.0, box_p->get_z_len()/2.0).norm();
    return sph;
  }

  // look for cone primitive
  shared_ptr<ConePrimitive> cone_p = dynamic_pointer_cast<ConePrimitive>(p);
  if (cone_p)
  {
    sph->radius = std::max(cone_p->get_height(), cone_p->get_radius());
    return sph;
  }

  // look for cylinder primitive
  shared_ptr<CylinderPrimitive> cyl_p = dynamic_pointer_cast<CylinderPrimitive>(p);
  if (cyl_p)
  {
    sph->radius = std::max(cyl_p->get_height(), cyl_p->get_radius());
    return sph;
  }

  // look for heightmap primitive
  shared_ptr<HeightmapPrimitive> hm_p = dynamic_pointer_cast<HeightmapPrimitive>(p);
  if (hm_p)
    return BVPtr(new DummyBV);

  // look for heightmap primitive
  shared_ptr<PlanePrimitive> plane_p = dynamic_pointer_cast<PlanePrimitive>(p);
  if (plane_p)
    return BVPtr(new DummyBV);

  // shouldn't still be here..
  assert(false);
  return BVPtr();
}

/****************************************************************************
 Methods for broad phase end
****************************************************************************/

/****************************************************************************
 Methods for solving LP start
****************************************************************************/

/// Solves a linear program using the method of Seidel
/**
 * This method exhibits complexity of O(d!n), where d is the dimension of the
 * variables and n is the number of constraints.
 * \param A the matrix for which Ax < b
 * \param b the vector for which Ax < b
 * \param c the optimization vector (maximizes c'x)
 * \param l the lower variable constraints on x
 * \param u the upper variable constraints on x
 * \param x the optimal solution (on successful return)
 * \return <b>true</b> if successful, <b>false</b> otherwise
 * \note using limits of +/- inf can result in overflow with this algorithm
 *       and is not recommended; use lower limits
 */
bool CCD::lp_seidel(const MatrixNd& A, const VectorNd& b, const VectorNd& c, const VectorNd& l, const VectorNd& u, VectorNd& x)
{
  // get number of rows and columns in A
  unsigned n = A.rows();
  unsigned d = A.columns();

  FILE_LOG(LOG_OPT) << "LP::lp() entered" << endl;
  FILE_LOG(LOG_OPT) << "A: " << endl << A;
  FILE_LOG(LOG_OPT) << "b: " << b << endl;
  FILE_LOG(LOG_OPT) << "c: " << c << endl;
  FILE_LOG(LOG_OPT) << "l: " << l << endl;
  FILE_LOG(LOG_OPT) << "u: " << u << endl;

  // base case d = 1
  if (d == 1)
  {
    FILE_LOG(LOG_OPT) << "base case, d = 1" << endl;

    double high = u[0]; 
    double low = l[0];

    for (unsigned i=0; i< n; i++)
    {
      if (A(i,0) > std::numeric_limits<double>::epsilon())
        high = std::min(high, b[i]/A(i,0));
      else if (A(i,0) < -std::numeric_limits<double>::epsilon())
        low = std::max(low, b[i]/A(i,0));
      else if (b[i] < -std::numeric_limits<double>::epsilon())
      {
        FILE_LOG(LOG_OPT) << "infeasible; b is negative and A is zero" << endl;

        return false; 
      }
    }

    // do a check for infeasibility
    if (high < low)
    {
      FILE_LOG(LOG_OPT) << "infeasible; high (" << high << ") < low (" << low << ")" << endl;

      return false; 
    }
  
    // set x
    x.resize(1);
    x[0] = (c[0] >= 0.0) ? high : low;

    FILE_LOG(LOG_OPT) << "optimal 1D x=" << x << endl;
    FILE_LOG(LOG_OPT) << "LP::lp() exited" << endl;

    // otherwise, good return
    return true; 
  }

  // all work variables
  #ifdef REENTRANT
  vector<VectorNd> aac;
  vector<VectorNd> aa;
  vector<unsigned> permut;
  vector<double> bbc;
  MatrixNd Aprime;
  VectorNd bb, aak, cprime, lprime, uprime, bprime, workv, f, g;
  #else 
  static vector<VectorNd> aac;
  static vector<VectorNd> aa;
  static vector<unsigned> permut;
  static vector<double> bbc;
  static MatrixNd Aprime;
  static VectorNd bb, aak, cprime, lprime, uprime, bprime, workv, f, g;
  #endif

  // pick a random shuffle for A and b
  permut.resize(n);
  for (unsigned i=0; i< n; i++)
    permut[i] = i;
  std::random_shuffle(permut.begin(), permut.end());

  // setup aa and bb
  aa.resize(n);
  bb.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    A.get_row(permut[i], aa[i]);
    bb[i] = b[permut[i]];
  }

  FILE_LOG(LOG_OPT) << "A (permuted): " << endl;
  for (unsigned i=0; i< n; i++)
    FILE_LOG(LOG_OPT) << aa[i] << endl;
  FILE_LOG(LOG_OPT) << "b (permuted): " << bb << endl;

  // setup optimum vector
  x.resize(d);
  for (unsigned i=0; i< d; i++)
    if (c[i] > 0.0)
        x[i] = u[i];
    else if (c[i] < 0.0) 
        x[i] = l[i];
    else
        x[i] = (std::fabs(l[i]) < std::fabs(u[i])) ? l[i] : u[i];

  FILE_LOG(LOG_OPT) << "initial x=" << x << endl;

  // process half-space constraints 
  for (unsigned i=0; i< n; i++)
  {
    FILE_LOG(LOG_OPT) << "-- processing halfspace constraint " << i << endl;

    // if x respects new halfspace constraint, nothing else to do..
    double val = bb[i] - VectorNd::dot(aa[i], x);
    if (val >= -EPS_DOUBLE)
    {    
      FILE_LOG(LOG_OPT) << "  ** constraint already satisfied!" << endl;
      continue;
    }

    FILE_LOG(LOG_OPT) << "  -- constraint not satisfied (" << val << "); solving recursively" << endl;

    // search for maximum value in the a vector
    unsigned k = std::numeric_limits<unsigned>::max();
    double maximal = -std::numeric_limits<double>::max();
    for (unsigned j=0; j< d; j++)
      if (std::fabs(aa[i][j]) > maximal && aa[i][j] != 0.0)
      {
        maximal = std::fabs(aa[i][j]);
        k = j;
      }

    FILE_LOG(LOG_OPT) << "  -- determined k: " << k << endl;

    // look for infeasibility
    if (k == std::numeric_limits<unsigned>::max())
    {
      FILE_LOG(LOG_OPT) << "  -- k is infeasible; problem is infeasible" << endl;

      return false; 
    }

    // setup useful vector and constant
    aak = aa[i];
    aak /= aa[i][k];
    double bak = bb[i]/aa[i][k];

    FILE_LOG(LOG_OPT) << "  -- vector a/a(k): " << aak << " " << bak << endl;

    // copy vectors aa and bb
    aac.resize(i);
    bbc.resize(i);
    for (unsigned j=0; j< i; j++)
    {
      aac[j] = aa[j];
      bbc[j] = bb[j];
    }

    // modify copy of vector aa
    for (unsigned j=0; j< i; j++)
    {
      workv = aak;
      workv *= aa[j][k];
      aac[j] -= workv;
      assert(std::fabs(aac[j][k]) < EPS_DOUBLE);
      aac[j] = remove_component(aac[j], k, workv);
    }

    // modify copy of vector bb 
    for (unsigned j=0; j< i; j++)
    {
      bbc[j] -= bak * aa[j][k];
      if (std::isinf(bbc[j]))
        bbc[j] = finitize(bbc[j]);
    }

    // modify copy of c
    assert(std::fabs((c[k] - aak[k]*c[k])) < EPS_DOUBLE);
    workv = aak;
    workv *= c[k];
    workv -= c;
    workv.negate(); 
    remove_component(workv, k, cprime);  

    // generate new lower and upper bounds for variables
    remove_component(l, k, lprime);
    remove_component(u, k, uprime);

    // setup two new constraints 
    f.set_zero(d);
    g.set_zero(d);
    f[k] = 1;
    g[k] = -1;
    assert(std::fabs((f[k] - aak[k])) < EPS_DOUBLE);
    assert(std::fabs((g[k] + aak[k])) < EPS_DOUBLE);
    f = remove_component(f -= aak, k, workv);
    g = remove_component(g += aak, k, workv);
    double bf = u[k] - bak;
    double bg = -l[k] + bak;
    if (std::isinf(bf))
      bf = finitize(bf);
    if (std::isinf(bg))
      bg = finitize(bg);
    aac.push_back(f);
    bbc.push_back(bf);
    aac.push_back(g);
    bbc.push_back(bg);

    // create the Aprime matrix from aac
    Aprime.resize(aac.size(), d-1);  
    for (unsigned j=0; j< aac.size(); j++)
      Aprime.set_row(j, aac[j]);

    // create the bprime vector from bbc
    bprime.resize(bbc.size());
    for (unsigned j=0; j< bbc.size(); j++)
      bprime[j] = bbc[j];

    FILE_LOG(LOG_OPT) << "  -- A': " << endl << Aprime;
    FILE_LOG(LOG_OPT) << "  -- b': " << bprime << endl;
    FILE_LOG(LOG_OPT) << "  -- c': " << cprime << endl;
    FILE_LOG(LOG_OPT) << "  -- u': " << uprime << endl;
    FILE_LOG(LOG_OPT) << "  -- l': " << lprime << endl;
    FILE_LOG(LOG_OPT) << "  -- f: " << f << " " << bf << endl;
    FILE_LOG(LOG_OPT) << "  -- g: " << g << " " << bg << endl;
    FILE_LOG(LOG_OPT) << " + solving recursive subproblem" << endl;

    // solve the (d-1)-dimensional problem and ``lift'' the solution
    if (!lp_seidel(Aprime,bprime,cprime,lprime,uprime,x))
      return false;

    FILE_LOG(LOG_OPT) << "  -- recursively determined x: " << x << endl;
    FILE_LOG(LOG_OPT) << "  -- k: " << k << endl;

    // insert a zero into the k'th dimension of x
    x = insert_component(x, k, workv);
    FILE_LOG(LOG_OPT) << "  -- x w/inserted component at k: " << x << endl;

    // solve for the proper k'th value of x
    x[k] = (bb[i] - VectorNd::dot(aa[i], x))/aa[i][k];
    FILE_LOG(LOG_OPT) << "  -- optimal x (to this point): " << x << endl;

    // verify that half-plane constraints still met
    for (unsigned j=0; j<= i; j++)
      assert(VectorNd::dot(aa[j], x) - bb[j] <= EPS_DOUBLE);

    // verify that lower constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] >= l[j] - EPS_DOUBLE);

    // verify that upper constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] <= u[j] + EPS_DOUBLE);
  }

  // verify that half-plane constraints still met
  if (LOGGING(LOG_OPT))
  {
    for (unsigned i=0; i< b.rows(); i++)
      FILE_LOG(LOG_OPT) << i << ": b - A*x = " << (b[i] - VectorNd::dot(A.row(i), x)) << endl;
  }
  for (unsigned j=0; j < n; j++)
    assert(VectorNd::dot(aa[j], x) - bb[j] <= EPS_DOUBLE);

  // verify that lower constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] >= l[i] - EPS_DOUBLE);

  // verify that upper constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] <= u[i] + EPS_DOUBLE);
  FILE_LOG(LOG_OPT) << "all halfspace constraints satisfied; optimum found!" << endl;
  FILE_LOG(LOG_OPT) << "optimum = " << x << endl;
  FILE_LOG(LOG_OPT) << "LP::lp_seidel() exited" << endl;

  return true; 
}

/// Inserts a component (value will be zero) at the k'th position in the given vector and returns a new vector
VectorNd& CCD::insert_component(const VectorNd& x, unsigned k, VectorNd& xn)
{
  xn.resize(x.size()+1);

  if (k == 0)
    xn.set_sub_vec(1, x);
  else if (k == x.size())
    xn.set_sub_vec(0, x);
  else
  {
    xn.set_sub_vec(0,x.segment(0,k));
    xn.set_sub_vec(k+1,x.segment(k,x.size()));
  }  
  xn[k] = 0;
  
  return xn;
}

/// Gets a subvector with the k'th component removed
/**
 * \note for use by lp()
 */
VectorNd& CCD::remove_component(const VectorNd& v, unsigned k, VectorNd& vp)
{
  const unsigned d = v.size();

  if (k == 0)
    return v.get_sub_vec(1,d,vp);
  else if (k == d-1)
    return v.get_sub_vec(0,d-1,vp);
  else
  {
    // setup v w/component k removed
    vp.resize(d-1);    
    vp.set_sub_vec(0, v.segment(0, k));
    vp.set_sub_vec(k, v.segment(k+1, d)); 
    return vp;
  }
}

/// Makes a (possibly) infinite value finite again
/**
 * Converts values of -inf to -DBL_MAX and inf to DBL_MAX
 * \note utility function for lp()
 */
double CCD::finitize(double x)
{
  if (x == std::numeric_limits<double>::infinity())
    return std::numeric_limits<double>::max();
  else if (x == -std::numeric_limits<double>::infinity())
    return -std::numeric_limits<double>::max();
  else
    return x;
}

/****************************************************************************
 Methods for solving LP end 
****************************************************************************/

