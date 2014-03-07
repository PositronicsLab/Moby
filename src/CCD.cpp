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
  _l.resize(6);
  _u.resize(6);
  _c.resize(6);
}

/// Computes a conservative advancement step between two dynamic bodies
double CCD::calc_CA_step(DynamicBodyPtr b1, DynamicBodyPtr b2)
{
  // drill down to pairs of rigid bodies
  ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(b1);
  ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(b2);
  if (ab1 && ab2)
  {
    double dt = std::numeric_limits<double>::max();
    BOOST_FOREACH(RigidBodyPtr rb1, ab1->get_links())
      BOOST_FOREACH(RigidBodyPtr rb2, ab2->get_links())
      {
        dt = std::min(dt, calc_CA_step(rb1, rb2));
        if (dt < NEAR_ZERO)
          return dt;
      }

    return dt;
  }
  else if (ab1)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(b2);
    BOOST_FOREACH(RigidBodyPtr rb1, ab1->get_links())
    {
      dt = std::min(dt, calc_CA_step(rb1, rb2));
      if (dt < NEAR_ZERO)
        return dt;
    }

    return dt;
  }
  else if (ab2)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(b1);
    BOOST_FOREACH(RigidBodyPtr rb2, ab2->get_links())
    {
      dt = std::min(dt, calc_CA_step(rb1, rb2));
      if (dt < NEAR_ZERO)
        return dt;
    }

    return dt;
  }
  else
  {
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(b1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(b2);
    return calc_CA_step(rb1, rb2);
  }
}

/// Computes a conservative advancement step between two dynamic bodies
double CCD::find_next_contact_time(DynamicBodyPtr b1, DynamicBodyPtr b2)
{
  // drill down to pairs of rigid bodies
  ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(b1);
  ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(b2);
  if (ab1 && ab2)
  {
    double dt = std::numeric_limits<double>::max();
    BOOST_FOREACH(RigidBodyPtr rb1, ab1->get_links())
      BOOST_FOREACH(RigidBodyPtr rb2, ab2->get_links())
        dt = std::min(dt, find_next_contact_time(rb1, rb2));

    return dt;
  }
  else if (ab1)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(b2);
    BOOST_FOREACH(RigidBodyPtr rb1, ab1->get_links())
      dt = std::min(dt, find_next_contact_time(rb1, rb2));

    return dt;
  }
  else if (ab2)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(b1);
    BOOST_FOREACH(RigidBodyPtr rb2, ab2->get_links())
      dt = std::min(dt, find_next_contact_time(rb1, rb2));

    return dt;
  }
  else
  {
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(b1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(b2);
    return find_next_contact_time(rb1, rb2);
  }
}

/// Finds the next event time between two rigid bodies
double CCD::find_next_contact_time(RigidBodyPtr rbA, RigidBodyPtr rbB)
{
  double maxt = std::numeric_limits<double>::max();
  Point3d pA, pB;

  // don't do CA step if rA == rbB or both bodies are disabled
  if (rbA == rbB || (!rbA->is_enabled() && !rbB->is_enabled()))
    return maxt;

  FILE_LOG(LOG_COLDET) << "CCD::find_next_contact_time() entered" << std::endl;

  // get the distance and closest points between the two bodies
  BOOST_FOREACH(CollisionGeometryPtr cgA, rbA->geometries)
    BOOST_FOREACH(CollisionGeometryPtr cgB, rbB->geometries)
    {
      // compute distance and closest points
      double dist = CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB);
      FILE_LOG(LOG_COLDET) << " -- CCD: reported distance: " << dist << std::endl;

      // special case: distance is zero or less
      if (dist <= NEAR_ZERO)
        continue; 

      // get the closest points in the global frame
      Point3d pA0 = Pose3d::transform_point(GLOBAL, pA);
      Point3d pB0 = Pose3d::transform_point(GLOBAL, pB);

      // get the direction of the vector from body B to body A
      Vector3d d0 =  pA0 - pB0;
      double d0_norm = d0.norm();
      FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
      FILE_LOG(LOG_COLDET) << "reported distance is: " << dist << std::endl; 
 
      // get the direction of the vector (from body B to body A)
      Vector3d n0 = d0/d0_norm;
      Vector3d nA = Pose3d::transform_vector(rbA->get_pose(), n0);
      Vector3d nB = Pose3d::transform_vector(rbB->get_pose(), n0);

      // get rA and rB
      Point3d rA = Pose3d::transform_point(rbA->get_pose(), pA);
      Point3d rB = Pose3d::transform_point(rbB->get_pose(), pB);

      // compute the distance that body A can move toward body B
      double velA = calc_max_velocity(rbA, -nA, rA);

      // compute the distance that body B can move toward body A
      double velB = calc_max_velocity(rbB, nB, rB);

      // compute the total velocity 
      double total_vel = velA + velB;
      if (total_vel < 0.0)
        total_vel = 0.0;

      FILE_LOG(LOG_COLDET) << " -- CCD: normal: " << n0 << std::endl;
      FILE_LOG(LOG_COLDET) << " -- CCD: projected velocity from A: " << velA << std::endl;
      FILE_LOG(LOG_COLDET) << " -- CCD: projected velocity from B: " << velB << std::endl;

      // compute the maximum safe step
      maxt = std::min(maxt, dist/total_vel);
    }

  // return the maximum safe step
  return maxt;
}

/// Computes a conservative advancement step between two rigid bodies
double CCD::calc_CA_step(RigidBodyPtr rbA, RigidBodyPtr rbB)
{
  double maxt = std::numeric_limits<double>::max();
  Point3d pA, pB;

  // don't do CA step if rA == rbB
  if (rbA == rbB || (!rbA->is_enabled() && !rbB->is_enabled()))
    return maxt;

  // get the distance and closest points between the two bodies
  BOOST_FOREACH(CollisionGeometryPtr cgA, rbA->geometries)
    BOOST_FOREACH(CollisionGeometryPtr cgB, rbB->geometries)
    {
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
 
      // get the direction of the vector (from body B to body A)
      Vector3d n0 = d0/d0_norm;
      Vector3d nA = Pose3d::transform_vector(rbA->get_pose(), n0);
      Vector3d nB = Pose3d::transform_vector(rbB->get_pose(), n0);

      // get rA and rB
      Point3d rA = Pose3d::transform_point(rbA->get_pose(), pA);
      Point3d rB = Pose3d::transform_point(rbB->get_pose(), pB);

      // compute the distance that body A can move toward body B
      double dist_per_tA = calc_max_dist_per_t(rbA, -nA, rA);

      // compute the distance that body B can move toward body A
      double dist_per_tB = calc_max_dist_per_t(rbB, nB, rB);

      // compute the total distance
      double total_dist_per_t = dist_per_tA + dist_per_tB;
      if (total_dist_per_t < 0.0)
        total_dist_per_t = 0.0;

      FILE_LOG(LOG_COLDET) << "  distance: " << dist << std::endl;
      FILE_LOG(LOG_COLDET) << "  dist per tA: " << dist_per_tA << std::endl;
      FILE_LOG(LOG_COLDET) << "  dist per tB: " << dist_per_tB << std::endl;

      // compute the maximum safe step
      maxt = std::min(maxt, dist/total_dist_per_t);
    }

  FILE_LOG(LOG_COLDET) << "  maxt: " << maxt << std::endl;

  // return the maximum safe step
  return maxt;
}

/// Computes the maximum velocity along a particular direction (n) 
double CCD::calc_max_velocity(RigidBodyPtr rb, const Vector3d& n, const Vector3d& r)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the velocities at t0
  const SVelocityd& v0 = Pose3d::transform(rb->get_pose(), rb->get_velocity());
  Vector3d xd0 = v0.get_linear();
  Vector3d w0 = v0.get_angular();
  return n.dot(xd0 + Vector3d::cross(w0, r)); 
}

/// Solves the LP that maximizes <n, v + w x r>
double CCD::calc_max_dist_per_t(RigidBodyPtr rb, const Vector3d& n, const Vector3d& r)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the velocities at t0
  const SVelocityd& v0 = Pose3d::transform(rb->get_pose(), rb->get_velocity());
  Vector3d xd0 = v0.get_linear();
  Vector3d w0 = v0.get_angular();

  // get the bounds on acceleration
  const SAcceld& a0p = rb->get_accel_upper_bounds(); 
  const SAcceld& a0n = rb->get_accel_lower_bounds(); 
  assert(a0p.pose == v0.pose);
  assert(a0n.pose == v0.pose);

  // break acceleration bounds in linear and angular components
  Vector3d xddp = a0p.get_linear();
  Vector3d xddn = a0n.get_linear();
  Vector3d alphap = a0p.get_angular();
  Vector3d alphan = a0n.get_angular();

  // setup lower and upper bounds
  _l[0] = xddn[X];   _u[0] = xddp[X];
  _l[1] = xddn[Y];   _u[1] = xddp[Y];
  _l[2] = xddn[Z];   _u[2] = xddp[Z];
  _l[3] = alphan[X]; _u[3] = alphap[X];
  _l[4] = alphan[Y]; _u[4] = alphap[Y];
  _l[5] = alphan[Z]; _u[5] = alphap[Z];

  // setup c
  _c[0] = n[X];
  _c[1] = n[Y];
  _c[2] = n[Z]; 
  _c[3] = n[Z]*r[Y] - n[Y]*r[Z];
  _c[4] = -n[Z]*r[X] + n[X]*r[Z];
  _c[5] = n[Y]*r[X] - n[X]*r[Y];

  // solve the LP
  return n.dot(xd0 + Vector3d::cross(w0, r)) + solve_lp(_c, _l, _u, _x); 
}

/// Determines the appropriate elements of the bounds vector
void CCD::to_binary(unsigned i, const VectorNd& l, const VectorNd& u, VectorNd& x)
{
  const unsigned N = u.size();
  unsigned n2 = 1 << (N-1);
  x.resize(N);

  // setup the binary number
  for (unsigned j=0; j< N; j++)
  {
    x[j] = (i / n2 > 0) ? u[j] : l[j];
    i = i % n2;
    n2 /= 2;
  }
}

/// Solves a linear program with only box constraints 
/**
 * \param c the optimization vector (maximizes c'r)
 * \param l the lower variable constraints on r
 * \param u the upper variable constraints on r
 * \param x the optimal solution on return
 * \return the optimal value
 */
double CCD::solve_lp(const VectorNd& c, const VectorNd& l, const VectorNd& u, VectorNd& x)
{
  const unsigned X = 0, Y = 1, Z = 2;

/*
  // setup the first possibility first
  double best = -std::numeric_limits<double>::max();

  // enumerate through all 2^6 = 64 possibilities
  for (unsigned i=0; i< 64; i++)
  {
    to_binary(i, l, u, _workv);
    double value = c.dot(_workv);
    if (value > best)
    {
      best = value;
      x = _workv;
    }
  }

  return best;
*/
  x = u;
  return c.dot(x);
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
void CCD::broad_phase(const map<CollisionGeometryPtr, PosePair>& poses, vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
{
  FILE_LOG(LOG_COLDET) << "CCD::broad_phase() entered" << std::endl;

  // clear the vector of pairs to check
  to_check.clear();

  // sort the AABBs
  sort_AABBs(bodies);

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

void CCD::sort_AABBs(const map<CollisionGeometryPtr, PosePair>& poses)
{
  // if a geometry was added or removed, rebuild the vector of bounds
  // and copy it to x, y, z dimensions
  if (_rebuild_bounds_vecs)
  {
    build_bv_vector(poses, _x_bounds);
    _y_bounds = _x_bounds;
    _z_bounds = _x_bounds;
  }

  // update bounds vectors
  update_bounds_vector(_x_bounds, poses, eXAxis);
  update_bounds_vector(_y_bounds, poses, eYAxis);
  update_bounds_vector(_z_bounds, poses, eZAxis);
  
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

void CCD::update_bounds_vector(vector<pair<double, BoundsStruct> >& bounds, RigidBodyPtr body, AxisType axis)
{
  const unsigned X = 0, Y = 1, Z = 2;

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() entered (axis=" << axis << ")" << std::endl;

  // iterate over bounds vector
  for (unsigned i=0; i< bounds.size(); i++)
  {
    // get the bounding volume and collision geometry
    BVPtr bv = bounds[i].second.bv;
    CollisionGeometryPtr geom = bounds[i].second.geom;

    // create an axis-aligned bounding box around the bounding volume 

    // get the current pose of the geometry in the global frame 

    // compute the swept BV as an axis-aligned bounding box

    // get the bound for the bounding volume
    Point3d bound = (bounds[i].second.end) ? swept.get_upper_bounds() : swept.get_lower_bounds();
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

void CCD::build_bv_vector(const map<CollisionGeometryPtr, PosePair>& poses, vector<pair<double, BoundsStruct> >& bounds)
{
  const double INF = std::numeric_limits<double>::max();

  // clear the vector
  bounds.clear();

  // iterate over all collision geometries
  for (set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
  {
    // get the primitive for the geometry
    PrimitivePtr p = (*i)->get_geometry();

    // get the top-level BV for the geometry
    BVPtr bv = p->get_BVH_root(*i);

    // get the poses for the geometry 
    assert(poses.find(*i) != poses.end());

    // setup the bounds structure
    BoundsStruct bs;
    bs.end = false;
    bs.geom = *i;
    bs.bv = bv;

    // add the lower bound
    bounds.push_back(make_pair(-INF, bs));

    // modify the bounds structure to indicate the end bound
    bs.end = true;
    bounds.push_back(make_pair(INF, bs));
  }
}

/****************************************************************************
 Methods for broad phase end 
****************************************************************************/


