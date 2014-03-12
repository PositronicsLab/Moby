/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
#include <Moby/Event.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>  
#include <Moby/XMLTree.h>
#include <Moby/SSL.h>
#include <Moby/BoundingSphere.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
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

/// Finds the next event time between two rigid bodies
double CCD::find_next_contact_time(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB)
{
  Point3d pA, pB;

  FILE_LOG(LOG_COLDET) << "CCD::find_next_contact_time() entered" << std::endl;

  // compute distance and closest points
  double dist = CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB);
  FILE_LOG(LOG_COLDET) << " -- CCD: reported distance: " << dist << std::endl;

  // special case: distance is zero or less
  if (dist <= NEAR_ZERO)
    return std::numeric_limits<double>::max(); 

  // get the closest points in the global frame
  Point3d pA0 = Pose3d::transform_point(GLOBAL, pA);
  Point3d pB0 = Pose3d::transform_point(GLOBAL, pB);

  // get the direction of the vector from body B to body A
  Vector3d d0 =  pA0 - pB0;
  double d0_norm = d0.norm();
  FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
  FILE_LOG(LOG_COLDET) << "reported distance is: " << dist << std::endl; 
 
  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());

  // get the direction of the vector (from body B to body A)
  Vector3d n0 = d0/d0_norm;
  Vector3d nA = Pose3d::transform_vector(rbA->get_pose(), n0);
  Vector3d nB = Pose3d::transform_vector(rbB->get_pose(), n0);

  // compute the distance that body A can move toward body B
  double velA = calc_max_velocity(rbA, -nA, _rmax[cgA]);

  // compute the distance that body B can move toward body A
  double velB = calc_max_velocity(rbB, nB, _rmax[cgB]);

  // compute the total velocity 
  double total_vel = velA + velB;
  if (total_vel < 0.0)
    total_vel = 0.0;

  FILE_LOG(LOG_COLDET) << " -- CCD: normal: " << n0 << std::endl;
  FILE_LOG(LOG_COLDET) << " -- CCD: projected velocity from A: " << velA << std::endl;
  FILE_LOG(LOG_COLDET) << " -- CCD: projected velocity from B: " << velB << std::endl;

  // compute the maximum safe step
  return dist/total_vel;
}

/// Computes a conservative advancement step between two collision geometries 
double CCD::calc_CA_step(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB)
{
  double maxt = std::numeric_limits<double>::max();
  Point3d pA, pB;

  // compute distance and closest points
  double dist = CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB);

  // if the distance is zero, quit now
  if (dist < NEAR_ZERO)
    return 0.0;

  // get the direction of the vector from body B to body A
  Vector3d d0 = Pose3d::transform_point(GLOBAL, pA) - 
                Pose3d::transform_point(GLOBAL, pB);
  double d0_norm = d0.norm();
  FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
  FILE_LOG(LOG_COLDET) << "reported distance is: " << dist << std::endl; 

  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());

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

  FILE_LOG(LOG_COLDET) << "  distance: " << dist << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tA: " << dist_per_tA << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tB: " << dist_per_tB << std::endl;

  // compute the maximum safe step
  maxt = std::min(maxt, dist/total_dist_per_t);

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

  // get the bounds on acceleration
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

/// Creates a contact event given the bare-minimum info
Event CCD::create_contact(CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Vector3d& normal)
{
  Event e;
  e.event_type = Event::eContact;
  e.contact_point = point;
  e.contact_normal = normal;
  e.contact_geom1 = a;  
  e.contact_geom2 = b;  

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

/*
/// Determines whether there is a collision at the current position and orientation of the bodies
bool CCD::is_collision(double epsilon)
{
  // clear the set of colliding pairs and list of colliding triangles
  colliding_pairs.clear();
  colliding_tris.clear();

  // iterate over geometries 
  for (std::set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
  {
    // get the first geometry, its primitive, and its bounding volume 
    CollisionGeometryPtr g1 = *i;
    PrimitivePtr g1_primitive = g1->get_geometry();
    BVPtr bv1 = g1_primitive->get_BVH_root(g1);

    // get the pose for the geometry
    shared_ptr<const Pose3d> Pg1 = g1->get_pose();

    // loop through all other geometries
    std::set<CollisionGeometryPtr>::const_iterator j = i;
    j++;
    for (; j != _geoms.end(); j++)
    {
      // get the second geometry, its primitive, and its bounding volume 
      CollisionGeometryPtr g2 = *j;
      PrimitivePtr g2_primitive = g2->get_geometry();
      BVPtr bv2 = g2_primitive->get_BVH_root(g2);

      // see whether to check
      if (!is_checked(g1, g2))
        continue; 

      // get the pose for the second geometry 
      shared_ptr<const Pose3d> Pg2 = g2->get_pose();

      // compute the transform from g2 to g1 
      Transform3d g1Tg2 = Pose3d::calc_relative_pose(Pg2, Pg1); 

      // if intersects, add to colliding pairs
      if (intersect_BV_trees(bv1, bv2, g1Tg2, g1, g2))
        colliding_pairs.insert(make_sorted_pair(g1, g2));
    } 
  }

  return !colliding_pairs.empty();
}

/// Intersects two BV trees; returns <b>true</b> if one (or more) pair of the underlying triangles intersects
bool CCD::intersect_BV_trees(BVPtr a, BVPtr b, const Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b) 
{
  std::queue<tuple<BVPtr, BVPtr, bool> > q;

  // get address of the last colliding triangle pair on the queue
  CollidingTriPair* last = (colliding_tris.empty()) ? NULL : &colliding_tris.back();

  FILE_LOG(LOG_COLDET) << "CCD::intersect_BV_trees() entered" << endl;

  // intersect the BVs at the top level
  if (!BV::intersects(a, b, aTb))
  {
    FILE_LOG(LOG_COLDET) << "  no intersection at top-level BVs" << endl;
    FILE_LOG(LOG_COLDET) << "CCD::intersect_BV_trees() exited" << endl;

    return false;
  }

  // add a and b to the queue
  q.push(make_tuple(a, b, false));

  // drill down alternatingly until both trees are exhausted
  while (!q.empty())
  {
    // get the two nodes off of the front of the queue
    BVPtr bv1 = q.front().get<0>();
    BVPtr bv2 = q.front().get<1>();
    bool rev = q.front().get<2>();

    // no longer need the elm on the front of the queue
    q.pop();

    // check for bv1 and bv2 both leafs
    if (bv1->is_leaf() && bv2->is_leaf())
    {
      if (!rev)
        intersect_BV_leafs(bv1, bv2, aTb, geom_a, geom_b, std::back_inserter(colliding_tris));
      else
        intersect_BV_leafs(bv2, bv1, aTb, geom_a, geom_b, std::back_inserter(colliding_tris));

      // see whether we want to exit early
      if (mode == eFirstContact && !colliding_tris.empty() && last != &colliding_tris.back())
        return true;
    }

    // drill down through bv2, if possible
    if (bv2->is_leaf())
    {
      // check the children of o1
      BOOST_FOREACH(BVPtr child, bv1->children)
        if ((!rev && BV::intersects(child, bv2, aTb)) || (rev && BV::intersects(bv2, child, aTb)))
          q.push(make_tuple(child, bv2, rev));
    }
    else
      BOOST_FOREACH(BVPtr child, bv2->children)
      {
        if ((!rev && BV::intersects(bv1, child, aTb)) || (rev && BV::intersects(child, bv1, aTb)))
          q.push(make_tuple(child, bv1, !rev));
      }
  }

  // see whether we have an intersection
  if (!colliding_tris.empty() && last != &colliding_tris.back())
    return true;

  FILE_LOG(LOG_COLDET) << "  -- all intersection checks passed; no intersection" << endl;
  FILE_LOG(LOG_COLDET) << "CCD::intersect_BV_trees() exited" << endl;

  return false;
} 
*/

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
  update_bounds_vector(_x_bounds, eXAxis, dt);
  update_bounds_vector(_y_bounds, eYAxis, dt);
  update_bounds_vector(_z_bounds, eZAxis, dt);
  
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

void CCD::update_bounds_vector(vector<pair<double, BoundsStruct> >& bounds, AxisType axis, double dt)
{
  const unsigned X = 0, Y = 1, Z = 2;

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() entered (axis=" << axis << ")" << std::endl;

  // iterate over bounds vector
  for (unsigned i=0; i< bounds.size(); i++)
  {
    // get the bounding volume and collision geometry
    BVPtr bv = bounds[i].second.bv;
    CollisionGeometryPtr geom = bounds[i].second.geom;

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

  // shouldn't still be here..
  assert(false);
  return BVPtr();
}

/****************************************************************************
 Methods for broad phase end 
****************************************************************************/


