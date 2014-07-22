/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <stack>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <Moby/CompGeom.h>
#include <Moby/CollisionDetection.h>
#include <Moby/Event.h>
#include <Moby/Constants.h>
#include <Moby/Polyhedron.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>  
#include <Moby/XMLTree.h>
#include <Moby/Integrator.h>
#include <Moby/SSR.h>
#include <Moby/C2ACCD.h>

// To delete
#include <Moby/SpherePrimitive.h>
#include <Moby/BoxPrimitive.h>

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
using std::stack;
using namespace Ravelin;
using namespace Moby;

/// For Alistair's priority queue modification
struct SSRPair
{
  SSRPair(shared_ptr<SSR> a, shared_ptr<SSR> b)
  {
    first = a;
    second = b;
    priority = (double) 0.0;
  }

  SSRPair(shared_ptr<SSR> a, shared_ptr<SSR> b, double priority)
  {
    first = a;
    second = b;
    this->priority = priority;
  }

  bool operator<(const SSRPair& s) const
  {
    return priority < s.priority;
  }

  shared_ptr<SSR> first, second;
  double priority;
};

/// Constructs a collision detector with default tolerances
/**
 * eps_tolerance is set to NEAR_ZERO and TOI tolerance is set to 1e-2.
 */
C2ACCD::C2ACCD()
{
  eps_tolerance = 1e-5;
  alpha_tolerance = 1e-2;
  pthread_mutex_init(&_contact_mutex, NULL);
  return_all_contacts = true;
}

void C2ACCD::add_collision_geometry(CollisionGeometryPtr cg)
{
  CollisionDetection::add_collision_geometry(cg);
  build_BV_tree(cg);
}

void C2ACCD::add_rigid_body(RigidBodyPtr rb)
{
  CollisionDetection::add_rigid_body(rb);
}

void C2ACCD::add_articulated_body(ArticulatedBodyPtr abody, bool disable_adjacent)
{
  CollisionDetection::add_articulated_body(abody, disable_adjacent);
}

void C2ACCD::remove_collision_geometry(CollisionGeometryPtr cg)
{
  CollisionDetection::remove_collision_geometry(cg);
}

void C2ACCD::remove_all_collision_geometries()
{
  CollisionDetection::remove_all_collision_geometries();
}

void C2ACCD::remove_rigid_body(RigidBodyPtr rb)
{
  CollisionDetection::remove_rigid_body(rb);
}

void C2ACCD::remove_articulated_body(ArticulatedBodyPtr abody)
{
  CollisionDetection::remove_articulated_body(abody);
}

/// Determines whether there is a contact in the given time interval
/**
 * \pre body states are at time tf
 */
bool C2ACCD::is_contact(double dt, const vector<pair<DynamicBodyPtr, VectorNd> >& q0, const vector<pair<DynamicBodyPtr, VectorNd> >& q1,vector<Event>& contacts)
{
  typedef pair<CollisionGeometryPtr, BVPtr> CG_BV;

  // clear the vector of contacts 
  contacts.clear();

  FILE_LOG(LOG_COLDET) << "C2ACCD::is_contact() entered" << endl;

  // setup the set of collision geometries encountered
  vector<CollisionGeometryPtr> encountered;

  // check the geometries
  for (set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
    for (set<CollisionGeometryPtr>::const_iterator j = i; j != _geoms.end(); j++)
    {
      if (i == j)
        continue;

      // get the two geometries
      CollisionGeometryPtr a = *i;
      CollisionGeometryPtr b = *j;

      // verify that the pair should be checked
      if (!is_checked(a, b))
        continue;

      // get the two rigid bodies
      RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
      RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

      // verify that the two bodies are rigid
      if (!rba || !rbb)
        throw std::runtime_error("One or more bodies is not rigid; C2ACCD only works with rigid bodies");

      // test the geometries for contact
      check_geoms(dt, a, b, q0, q1, contacts);
    } 

  FILE_LOG(LOG_COLDET) << "contacts:" << endl;
  if (contacts.empty())
    FILE_LOG(LOG_COLDET) << " -- no contacts in narrow phase" << endl;
  if (LOGGING(LOG_COLDET))
    for (unsigned i=0; i< contacts.size(); i++)
      FILE_LOG(LOG_COLDET) << contacts[i] << std::endl;

  FILE_LOG(LOG_COLDET) << "C2ACCD::is_contact() exited" << endl << endl;

  // sort the contacts based on time
  std::sort(contacts.begin(), contacts.end());

  // indicate whether impact has occurred
  return !contacts.empty();
}

/// Gets the "super" body for a collision geometry
DynamicBodyPtr C2ACCD::get_super_body(CollisionGeometryPtr geom)
{
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(geom->get_single_body());
  assert(rb);
  ArticulatedBodyPtr ab = rb->get_articulated_body();
  if (ab)
    return ab;
  else
    return rb;
}

/// Finds the index of the body / state pair for the given body
unsigned C2ACCD::find_body(const vector<pair<DynamicBodyPtr, VectorNd> >& q, DynamicBodyPtr body)
{
  for (unsigned i=0; i< q.size(); i++)
    if (q[i].first == body)
      return i;

  return std::numeric_limits<unsigned>::max();
}

/// Does a collision check for a pair of geometries 
/**
 * \param dt the time interval
 * \param a the first geometry
 * \param b the second geometry
 * \param contacts on return
 */
void C2ACCD::check_geoms(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const vector<pair<DynamicBodyPtr, VectorNd> >& q0, const vector<pair<DynamicBodyPtr, VectorNd> >& q1, vector<Event>& contacts)
{
  VectorNd q, qtmp;

  FILE_LOG(LOG_COLDET) << "C2ACCD::check_geoms() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  checking geometry " << a->id << " for body " << a->get_single_body()->id << std::endl;
  FILE_LOG(LOG_COLDET) << "  against geometry " << b->id << " for body " << b->get_single_body()->id << std::endl;

  // get the SSR's for a and b
  shared_ptr<SSR> ssr_a = _root_SSRs[a];
  shared_ptr<SSR> ssr_b = _root_SSRs[b];

  // get bodies for a and b
  DynamicBodyPtr ba = get_super_body(a);
  DynamicBodyPtr bb = get_super_body(b);

  // get the states for bodies a and b
  const unsigned idx_a = find_body(q0, ba);
  const unsigned idx_b = find_body(q0, bb);
  assert(idx_a == find_body(q1, ba));
  assert(idx_b == find_body(q1, bb));
  assert(idx_a != std::numeric_limits<unsigned>::max());
  assert(idx_b != std::numeric_limits<unsigned>::max());
  const VectorNd& qa0 = q0[idx_a].second;
  const VectorNd& qa1 = q1[idx_a].second;
  const VectorNd& qb0 = q0[idx_b].second;
  const VectorNd& qb1 = q1[idx_b].second;

  // set bodies to states qa0 and qb0
  ba->set_generalized_coordinates(DynamicBody::eEuler, qa0);
  bb->set_generalized_coordinates(DynamicBody::eEuler, qb0);

  // step to TOC
  double TOC = (double) 0.0;
  double h;
  do
  {
    // get the transforms from a to b
    shared_ptr<const Pose3d> Ta = a->get_pose();
    shared_ptr<const Pose3d> Tb = b->get_pose();
    Transform3d aTb = Pose3d::calc_relative_pose(Tb, Ta);

    // determine a safe step we can take
    double dcur = std::numeric_limits<double>::max();
    h = do_CA(dt, a, b, ssr_a, ssr_b, aTb, (double) 1.0 - TOC);
    FILE_LOG(LOG_COLDET) << " -- determined h: " << h << "  TOC: " << TOC << endl;

    // update the TOC
    TOC += h;

    // advance the bodies' states to time TOC
    (q = qa0) *= ((double) 1.0 - TOC);
    (qtmp = qa1) *= TOC;
    q += qtmp;
    ba->set_generalized_coordinates(DynamicBody::eEuler, q);
    (q = qb0) *= ((double) 1.0 - TOC);
    (qtmp = qb1) *= TOC;
    q += qtmp;
    bb->set_generalized_coordinates(DynamicBody::eEuler, q);
  }
  while (h > alpha_tolerance && TOC < (double) 1.0);

  // if there is an impact, determine contacts
  if (TOC < (double) 1.0)
    determine_contacts(a, b, TOC, contacts);

  FILE_LOG(LOG_COLDET) << "TOC: " << TOC << endl;
  FILE_LOG(LOG_COLDET) << "C2ACCD::check_geoms() exited" << endl;
}

/// Attempts to do line segment / triangle intersection and returns the first point of contact (p) and associated parameter (t)
bool C2ACCD::query_intersect_seg_tri(const LineSeg3& seg, const Triangle& tri, double& t, Point3d& p)
{
  Point3d p1, p2;
  CompGeom::SegTriIntersectType itype = CompGeom::intersect_seg_tri(seg, tri, p1, p2);
  if (itype == CompGeom::eSegTriNoIntersect || 
      itype == CompGeom::eSegTriInclEdge ||
      itype == CompGeom::eSegTriInclVertex ||
      itype == CompGeom::eSegTriVertex ||
      itype == CompGeom::eSegTriEdge ||
      itype == CompGeom::eSegTriEdgeOverlap)
    return false;

  // determine parameters of both
  double seglen = (seg.second - seg.first).norm();
  double t1 = (p1 - seg.first).norm() / seglen;
  double t2 = (p2 - seg.first).norm() / seglen;
  if (t1 < t2)
  {
    t = t1;
    p = p1;
  }
  else
  {
    t = t2;
    p = p2;
  }

  return true;
}

/// Determines the contacts between two geometries (closest features approach)
/**
 * The bodies should be right at the point of contact.
 */
void C2ACCD::determine_contacts(CollisionGeometryPtr a, CollisionGeometryPtr b, double toc, vector<Event>& contacts) const
{
  FILE_LOG(LOG_COLDET) << "C2ACCD::determine_contacts() entered" << endl;

  // get the rigid bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

  FILE_LOG(LOG_COLDET) << "  body 1: " << rba->id << endl;
  FILE_LOG(LOG_COLDET) << "  body 2: " << rbb->id << endl;

  // setup vectors of triangles in collision
  vector<pair<unsigned, unsigned> > colliding_tris;

  // determine aTb
  shared_ptr<const Pose3d> Ta = a->get_pose();
  shared_ptr<const Pose3d> Tb = b->get_pose();
  Transform3d aTb = Pose3d::calc_relative_pose(Tb, Ta);

  // setup transformation from Ta to global frame
  Transform3d wTa = Pose3d::calc_relative_pose(Ta, GLOBAL);

  // determine closest triangles
  vector<pair<Triangle, Triangle> > closest_tris;
  determine_closest_tris(a, b, aTb, closest_tris);

  // determine features from closest triangles
  for (unsigned i=0; i< closest_tris.size(); i++)
  {
    // transform the two triangles to global frame & determine closest features
    Triangle tA = Triangle::transform(closest_tris[i].first, wTa);
    Triangle tB = Triangle::transform(closest_tris[i].second, wTa);
    Triangle::FeatureType fA, fB; 
    vector<Point3d> contact_points;
    determine_closest_features(tA, tB, fA, fB, contact_points);

    FILE_LOG(LOG_COLDET) << "  -- examining triangles " << tA << " and " << tB << endl;
    FILE_LOG(LOG_COLDET) << "  closest features: " << fA << " and " << fB << endl;

    // determine the normal
    const Triangle* tri = NULL;
    if (fB == Triangle::eFace)
    {
      tri = &tB; 
      FILE_LOG(LOG_COLDET) << "    normal comes from triangle B" << endl;
    }
    else
    {
      tri = &tA;
      FILE_LOG(LOG_COLDET) << "    normal comes from triangle A" << endl;
    }

    // process all contact points
    for (unsigned j=0; j< contact_points.size(); j++)
    {
      // get the normal
      Vector3d normal = tri->calc_normal();

      // setup the contact 
      Event e;
      e.t = toc;
      e.event_type = Event::eContact;
      e.contact_geom1 = a;
      e.contact_geom2 = b;
      e.contact_normal = normal;
      e.contact_point = contact_points[j];

      // get the signed distance
      double sdist = tri->calc_signed_dist(e.contact_point);

      // get the relative velocity at the contact point
      Vector3d pva = rba->calc_point_vel(e.contact_point);
      Vector3d pvb = rbb->calc_point_vel(e.contact_point);
      Vector3d rvel = Pose3d::transform_vector(normal.pose, pva) - Pose3d::transform_vector(normal.pose, pvb);

      FILE_LOG(LOG_COLDET) << " -- normal: " << normal << std::endl;
      FILE_LOG(LOG_COLDET) << " -- contact point: " << e.contact_point << " sdist: " << sdist << endl;
      FILE_LOG(LOG_COLDET) << " -- relative velocity: " << rvel << " r'n: " << normal.dot(rvel) << endl;

      // if the point is effectively on the plane, skip it
      if (std::fabs(sdist) < NEAR_ZERO)
        continue;

      // get the normal dotted with the point and the relative velocity
/*
      if (sdist < (double) 0.0 && normal.dot(rvel) > (double) 0.0)
      {
        FILE_LOG(LOG_COLDET) << "reversing normal!" << std::endl;
        c.normal = -c.normal;
      }
*/
      if (sdist < (double) 0.0)
      {
        if (tri == &tB)
          e.contact_normal = -e.contact_normal;
      }
      else if (tri == &tA)
        e.contact_normal = -e.contact_normal;

      // insert the contact into the map
      contacts.push_back(e);
    }
/*
    // determine the normal
    Vector3 normal .set_zero();
    if (fB == Triangle::eFace)
    {
      Vector3 tA_center = (tA.a + tA.b + tA.c)*0.333333;
      normal = tB.calc_normal();
      double sdist = tB.calc_signed_dist(tA_center);
      FILE_LOG(LOG_COLDET) << "    normal " << normal << " comes from triangle B" << endl;
      FILE_LOG(LOG_COLDET) << "      signed distance from center of A: " << sdist << endl;
      if (sdist < (double) 0.0)
      {
        normal = -normal;
        FILE_LOG(LOG_COLDET) << "      ** reversing normal!" << endl;
      }
    }
    else
    {
      Vector3 tB_center = (tB.a + tB.b + tB.c)*0.333333;
      normal = tA.calc_normal();
      double sdist = tA.calc_signed_dist(tB_center);
      FILE_LOG(LOG_COLDET) << "    normal " << normal << " comes from triangle A" << endl;
      FILE_LOG(LOG_COLDET) << "      signed distance from center of B: " << sdist << endl;
      if (tA.calc_signed_dist(tB_center) > (double) 0.0)
      {
        normal = -normal;
        FILE_LOG(LOG_COLDET) << "      ** reversing normal!" << endl;
      }
    }

    // create contacts
    for (unsigned j=0; j< contact_points.size(); j++)
    {
      Contact c;
      c.geom1 = a;
      c.geom2 = b;
      c.normal = normal;
      c.contact_geometry = Contact::ePoint;
      c.point = contact_points[j];
      contacts.insert(make_pair(toc, c));
    }
*/
  }
}

void C2ACCD::determine_closest_tris(CollisionGeometryPtr a, CollisionGeometryPtr b, const Transform3d& aTb, vector<pair<Triangle, Triangle> >& closest_tris) const
{
  Point3d dummy;

  // get the root SSRs
  shared_ptr<SSR> ssr_a = _root_SSRs.find(a)->second;
  shared_ptr<SSR> ssr_b = _root_SSRs.find(b)->second;

  // add them to a queue
  queue<pair<shared_ptr<SSR>, shared_ptr<SSR> > > Q;
  Q.push(make_pair(ssr_a, ssr_b));

  // setup the distance
  double max_dist = std::numeric_limits<double>::max();

  // process until the queue is empty
  while (!Q.empty())
  {
    // pop the pair off the front of the queue
    pair<shared_ptr<SSR>, shared_ptr<SSR> > ssrs = Q.front();
    Q.pop();

    // get the distance between the two SSRs
    double dist = SSR::calc_dist(*ssrs.first, *ssrs.second, aTb, dummy, dummy);
    if (dist > max_dist)
      continue;

    // if both SSRs are leafs, get the distance between the primitives
    if (ssrs.first->is_leaf() && ssrs.second->is_leaf())
    {
      // get triangles in ssr_a and ssr_b
      assert(_meshes.find(ssrs.first) != _meshes.end());
      assert(_meshes.find(ssrs.second) != _meshes.end());
      const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_a = _meshes.find(ssrs.first)->second;
      const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_b = _meshes.find(ssrs.second)->second;
      assert(!mesh_a.second.empty());
      assert(!mesh_b.second.empty());

      // find closest of all pairs of triangles in ssrs
      BOOST_FOREACH(unsigned ta_idx, mesh_a.second)
      {
        Triangle ta = mesh_a.first->get_triangle(ta_idx, aTb.target);
        BOOST_FOREACH(unsigned tb_idx, mesh_b.second)
        {
          Triangle utb = mesh_b.first->get_triangle(tb_idx, aTb.source);
          Triangle tb = Triangle::transform(utb, aTb);
          double dist = std::sqrt(Triangle::calc_sq_dist(ta, tb, dummy, dummy));
          if (dist + NEAR_ZERO < max_dist)
          {
            // distance is appreciably lower...  clear the vector of closest
            // triangles
            closest_tris.clear();
            closest_tris.push_back(make_pair(ta, tb));
            max_dist = dist;
          }
          else if (dist - eps_tolerance < max_dist)
          {
            // distance is as low/lower, but not so much to clear the vector
            closest_tris.push_back(make_pair(ta, tb));
            if (dist < max_dist)
              max_dist = dist;
          }
        }
      }
    }
    else if (ssrs.first->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, ssrs.second->children)
        Q.push(make_pair(ssrs.first, dynamic_pointer_cast<SSR>(child)));
    }
    else if (ssrs.second->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, ssrs.first->children)
        Q.push(make_pair(dynamic_pointer_cast<SSR>(child), ssrs.second));
    }
    else
    {
      // neither is a leaf; recurse into the larger volume
      if (ssrs.first->calc_volume() > ssrs.second->calc_volume())
        BOOST_FOREACH(BVPtr child, ssrs.first->children)
          Q.push(make_pair(dynamic_pointer_cast<SSR>(child), ssrs.second));
      else
        BOOST_FOREACH(BVPtr child, ssrs.second->children)
          Q.push(make_pair(ssrs.first, dynamic_pointer_cast<SSR>(child)));
    }
  }
}

/// Determines the contacts between two geometries (intersecting version)
/*
void C2ACCD::determine_contacts(CollisionGeometryPtr a, CollisionGeometryPtr b, double dt, double toc, vector<Event>& contacts) const
{
  const unsigned N_TRI_VERTS = 3;
  double beta = 0.001;

  FILE_LOG(LOG_COLDET) << "C2ACCD::determine_contacts() entered" << endl;

  // get the rigid bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

  // setup vectors of triangles in collision
  vector<pair<unsigned, unsigned> > colliding_tris;

  // save the bodies' states
  VectorNd qa, qadot, qb, qbdot;
  rba->get_position_state(qa);
  rbb->get_position_state(qb);

  do
  {
    // step the bodies positions forward by toc + alpha*dt
    rba->get_position_state_deriv(qadot);
    rbb->get_position_state_deriv(qbdot);
    qadot *= (toc + dt*beta);
    qbdot *= (toc + dt*beta);
    qadot += qa;
    qbdot += qb;
    rba->set_position_state(qadot);
    rbb->set_position_state(qbdot);
    beta *= 2.0;
    if (beta > 1.0)
    {
      FILE_LOG(LOG_COLDET) << " *** beta too large!  exiting..." << endl;
      FILE_LOG(LOG_COLDET) << "C2ACCD::determine_contacts() exited" << endl;
      return;
    }
  }
  while (!check_collision(a, b, colliding_tris));

  // now, restore the rigid body states and step to right before collision
  rba->set_position_state(qa);
  rbb->get_position_state(qb);
  rba->get_position_state_deriv(qadot);
  rbb->get_position_state_deriv(qbdot);
  qadot *= toc;
  qbdot *= toc;
  qadot += qa;
  qbdot += qb;
  rba->set_position_state(qadot);
  rbb->set_position_state(qbdot);

  // get the root SSRs
  shared_ptr<SSR> ssr_a = _root_SSRs.find(a)->second;
  shared_ptr<SSR> ssr_b = _root_SSRs.find(b)->second;

  // get the meshes
  assert(_meshes.find(ssr_a) != _meshes.end());
  assert(_meshes.find(ssr_b) != _meshes.end());
  shared_ptr<const IndexedTriArray> mesh_a = _meshes.find(ssr_a)->second.first;
  shared_ptr<const IndexedTriArray> mesh_b = _meshes.find(ssr_b)->second.first;

  // get the transforms for rigid bodies a and b
  const Pose3d& Ta = rba->get_pose();
  const Pose3d& Tb = rbb->get_pose();

  // setup a multimap for storing contact parameters and data
  vector<Event> contact_store;

  // for each pair of triangles, determine contacts
  for (unsigned i=0; i< colliding_tris.size(); i++)
  {
    // get the two indices
    unsigned ta_idx = colliding_tris[i].first;
    unsigned tb_idx = colliding_tris[i].second;

    // get the triangles in the global frame
    Triangle ta = Triangle::transform(mesh_a->get_triangle(ta_idx), Ta);
    Triangle tb = Triangle::transform(mesh_b->get_triangle(tb_idx), Tb);

    // check vertices of ta against tb
    for (unsigned j=0; j< N_TRI_VERTS; j++)
    {
      // get the vertex
      const Vector3& v = ta.get_vertex(j);

      // determine the relative point velocity
      Vector3 v_dot = rba->calc_point_vel(v) - rbb->calc_point_vel(v);

      // calculate the first intersection, if any
      double tparam;
      Vector3 isect;
      FILE_LOG(LOG_COLDET) << " -- checking intersection between triangle " << tb << " and line seg " << v << " / " << (v+v_dot*dt) << endl;
      if (query_intersect_seg_tri(LineSeg3(v, v+v_dot*dt), tb, tparam, isect))
      {
        Contact c;
        c.geom1 = a;
        c.geom2 = b;
        c.contact_geometry = Contact::ePoint;
        c.point = isect;
        c.normal = tb.calc_normal();
        if (tb.calc_signed_dist(v) < (double) 0.0)
          c.normal = -c.normal;
        contact_store.insert(make_pair(tparam, c));
      }
    }

    // check vertices of tb against ta
    for (unsigned j=0; j< N_TRI_VERTS; j++)
    {
      // get the vertex
      const Vector3& v = tb.get_vertex(j);

      // determine the relative point velocity
      Vector3 v_dot = rbb->calc_point_vel(v) - rba->calc_point_vel(v);

      // calculate the first intersection, if any
      double tparam;
      Vector3 isect;
      FILE_LOG(LOG_COLDET) << " -- checking intersection between triangle " << ta << " and line seg " << v << " / " << (v+v_dot*dt) << endl;
      if (query_intersect_seg_tri(LineSeg3(v, v+v_dot*dt), ta, tparam, isect))
      {
        Contact c;
        c.geom1 = a;
        c.geom2 = b;
        c.contact_geometry = Contact::ePoint;
        c.point = isect;
        c.normal = -ta.calc_normal();
        if (ta.calc_signed_dist(v) < (double) 0.0)
          c.normal = -c.normal;
        contact_store.insert(make_pair(tparam, c));
      }
    }
  }

  // restore the bodies' states 
  rba->set_position_state(qa);
  rbb->get_position_state(qb);

  // move contacts into the contact store
  if (contact_store.empty())
  {
    FILE_LOG(LOG_COLDET) << "  no intersections found!" << endl;
    FILE_LOG(LOG_COLDET) << "C2ACCD::determine_contacts() exited" << endl;
    return;
  }
  double first_t = contact_store.begin()->first;
  for (vector<Event>::const_iterator i = contact_store.begin(); i != contact_store.end(); i++)
  {
    if (i->first > first_t + NEAR_ZERO)
      break;
    contacts.insert(make_pair(toc, i->second));
  }

  FILE_LOG(LOG_COLDET) << "C2ACCD::determine_contacts() exited" << endl;
}
*/

/// Checks two geometries for collision at their current transforms and stores pairs of colliding triangles
bool C2ACCD::check_collision(CollisionGeometryPtr a, CollisionGeometryPtr b, vector<pair<unsigned, unsigned> >& colliding_tris) const
{
  // clear vector of colliding triangles
  colliding_tris.clear();

  // get the root SSRs
  shared_ptr<SSR> ssr_a = _root_SSRs.find(a)->second;
  shared_ptr<SSR> ssr_b = _root_SSRs.find(b)->second;

  // compute transformation from b's frame to a's frame
  shared_ptr<const Pose3d> Ta = a->get_pose();
  shared_ptr<const Pose3d> Tb = b->get_pose();
  Transform3d aTb = Pose3d::calc_relative_pose(Tb, Ta);

  // setup a queue for checking
  queue<pair<shared_ptr<SSR>, shared_ptr<SSR> > > Q;
  Q.push(make_pair(ssr_a, ssr_b));

  // process until the queue is empty
  while (!Q.empty())
  {
    // get the pair off of the front of the queue
    pair<shared_ptr<SSR>, shared_ptr<SSR> > test = Q.front();
    Q.pop();

    // check whether they intersect
    if (!BV::intersects(test.first, test.second, aTb))
      continue;

    // they do intersect...  if they are both leaf nodes, check individual
    // triangles for intersection
    if (test.first->is_leaf() && test.second->is_leaf())
    {
      // get the triangles enclosed by both SSRs
      const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_a = _meshes.find(test.first)->second;
      const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_b = _meshes.find(test.second)->second;

      // do pairwise triangle intersections
      BOOST_FOREACH(unsigned ta_idx, mesh_a.second)
      {
        Triangle ta = mesh_a.first->get_triangle(ta_idx, aTb.target);
        BOOST_FOREACH(unsigned tb_idx, mesh_b.second)
        {
          Triangle utb = mesh_b.first->get_triangle(tb_idx, aTb.source);
          Triangle tb = Triangle::transform(utb, aTb);
          if (CompGeom::query_intersect_tri_tri(ta, tb))
            colliding_tris.push_back(make_pair(ta_idx, tb_idx));
        }
      }
    }
    else if (test.first->is_leaf())
    {
      // test against children of the SSR from b
      BOOST_FOREACH(BVPtr child, test.second->children)
        Q.push(make_pair(test.first, dynamic_pointer_cast<SSR>(child)));
    }
    else if (test.second->is_leaf())
    {
      // test against children of the SSR from a
      BOOST_FOREACH(BVPtr child, test.first->children)
        Q.push(make_pair(dynamic_pointer_cast<SSR>(child), test.second));
    }
    else
    {
      // neither is a leaf; test against children of SSR with larger volume
      double vol_a = test.first->calc_volume();
      double vol_b = test.second->calc_volume();
      if (vol_a < vol_b)
      {
        BOOST_FOREACH(BVPtr child, test.second->children)
          Q.push(make_pair(test.first, dynamic_pointer_cast<SSR>(child)));
      }
      else
      {
        BOOST_FOREACH(BVPtr child, test.first->children)
          Q.push(make_pair(dynamic_pointer_cast<SSR>(child), test.second));
      }
    }
  }

  return !colliding_tris.empty();
}

/// Determines the motion bound on a SSR
double C2ACCD::calc_mu(double dist, const Vector3d& n, CollisionGeometryPtr g, boost::shared_ptr<SSR> ssr, bool positive)
{
  const unsigned N_RECT_VERTS = 4;

  // get the rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(g->get_single_body());

  // get the velocity of the rigid body in the global frame
  SVelocityd v = Pose3d::transform(GLOBAL, rb->get_velocity());

  // calculate linear contribution to mu
  double mu = n.dot(v.get_linear());
  if (!positive)
    mu = -mu;
//  mu = std::fabs(mu);
// NOTE: uncomment the line 'mu = std::fabs(mu)' and comment out the line 
// 'mu = -mu' to use the old C2A motion bound

  FILE_LOG(LOG_COLDET) << " linear contribution to mu: " << mu << std::endl;

  // get the angular velocity of rb and the angular speed of rb
  const Vector3d& omega = v.get_linear();
  double aspeed = omega.norm();

  // two cases: ssr is an internal node or a leaf node
  if (ssr->is_leaf())
  {
    FILE_LOG(LOG_COLDET) << "  ** SSR is a leaf node" << endl;

    // ** calculate directed motion bound on triangle vertices
    // get all triangle vertices (local frame)
    std::vector<Point3d> tri_verts;
    assert(_meshes.find(ssr) != _meshes.end());
    pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_tris = _meshes.find(ssr)->second;
    BOOST_FOREACH(unsigned tri_idx, mesh_tris.second)
    {
      Triangle t = mesh_tris.first->get_triangle(tri_idx, g->get_pose());
      tri_verts.push_back(t.a);
      tri_verts.push_back(t.b);    
      tri_verts.push_back(t.c);    
      FILE_LOG(LOG_COLDET) << "    triangle vertex (local): " << t.a << endl; 
      FILE_LOG(LOG_COLDET) << "    triangle vertex (local): " << t.b << endl; 
      FILE_LOG(LOG_COLDET) << "    triangle vertex (local): " << t.c << endl; 
    }

    // determine the maximum projection
    double mx_proj = (double) 0.0;
    for (unsigned i=0; i< tri_verts.size(); i++)
    {
      double proj = Vector3d::cross(Vector3d(tri_verts[i]), omega).norm();
      if (proj > mx_proj)
        mx_proj = proj;
    } 

    // normalize by angular speed
    mx_proj /= aspeed;

    FILE_LOG(LOG_COLDET) << "  -- maximum projection: " << mx_proj << endl;

    // update mu
    if (mx_proj > (double) 0.0)
      mu += Vector3d::cross(omega, n).norm() * mx_proj;
  }
  else
  {
    FILE_LOG(LOG_COLDET) << "  ** SSR is an internal node" << endl;

    // ** calculate direction motion bound on SSR rectangle vertices
    // get the rectangle vertices (local frame)
    Point3d rect_verts[N_RECT_VERTS];
    ssr->get_rect_verts(rect_verts);

    FILE_LOG(LOG_COLDET) << "    rectangle vertex 1: " << rect_verts[0] << endl;
    FILE_LOG(LOG_COLDET) << "    rectangle vertex 2: " << rect_verts[1] << endl;
    FILE_LOG(LOG_COLDET) << "    rectangle vertex 3: " << rect_verts[2] << endl;
    FILE_LOG(LOG_COLDET) << "    rectangle vertex 4: " << rect_verts[3] << endl;

    // determine the maximum projection
    double mx_proj = (double) 0.0;
    for (unsigned i=0; i< N_RECT_VERTS; i++)
      mx_proj = std::max(mx_proj, Vector3d::cross(Vector3d(rect_verts[i]), omega).norm());

    // normalize with aspeed
    mx_proj /= aspeed;

    FILE_LOG(LOG_COLDET) << "  -- maximum projection: " << mx_proj << endl;

    // update mu
    if (mx_proj > (double) 0.0)
      mu += Vector3d::cross(omega, n).norm()*(mx_proj + ssr->radius);
  }

  FILE_LOG(LOG_COLDET) << "  -- determined mu: " << mu << endl;

  return mu;
}

/// CA step (Mirtich's approach)
double C2ACCD::do_CAStep(double dist, const Vector3d& dab, CollisionGeometryPtr a, CollisionGeometryPtr b, shared_ptr<SSR> ssr_a, shared_ptr<SSR> ssr_b)
{
  FILE_LOG(LOG_COLDET) << "C2ACCD::do_CAStep() entered" << endl;

  if (dist <= (double) 0.0)
  {
    FILE_LOG(LOG_COLDET) << "  distance is zero; exiting..." << endl;
    return 0.0;
  } 

  // determine n
  double dab_len = dab.norm();
  Vector3d n = dab/dab_len;

  // get the individual bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

  // determine the motion bound on each 
  double mu_a = calc_mu(dist, n, a, ssr_a, true);
  double mu_b = calc_mu(dist, n, b, ssr_b, false);
  double mu = mu_a + mu_b;

  // if mu is negative, there is no max step
  if (mu <= (double) 0.0)
    return std::numeric_limits<double>::max();

  // solve for the maximum step 
  double mdt = std::max(dist - eps_tolerance, (double) 0.0)/mu; 

  FILE_LOG(LOG_COLDET) << "C2ACCD::do_CAStep() determined max dt: " << mdt << endl;
  return mdt;
}

/// CA algorithm (no "control")
double C2ACCD::do_CA(double step_size, CollisionGeometryPtr a, CollisionGeometryPtr b, shared_ptr<SSR> ssr_a, shared_ptr<SSR> ssr_b, const Transform3d& aTb, double dt)
{
  const double TTOL = alpha_tolerance / step_size;
  FILE_LOG(LOG_COLDET) << "C2ACCD::do_CA() entered" << endl;

  // place the two SSR's onto a priority queue
  priority_queue<SSRPair> Q;
  Q.push(SSRPair(ssr_a, ssr_b));

  while (!Q.empty())
  {
    // get the pair off of the front of the queue
    SSRPair ssrs = Q.top();
    Q.pop();

    // if the SSR's are both leafs, we may update dt
    if (ssrs.first->is_leaf() && ssrs.second->is_leaf())
    {
      // compute the distance between the triangles in the SSRs and get closest
      // points in global frame
      Point3d cpa, cpb;
      double dist = calc_dist(ssrs.first, ssrs.second, aTb, cpa, cpb);

      // indicate closest points are in a's frame
      shared_ptr<const Pose3d> Ta = a->get_pose();
      cpa.pose = Ta;
      cpb.pose = Ta;

      // transform closest points to global frame
      cpa = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpa);
      cpb = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpb);
      FILE_LOG(LOG_COLDET) << " -- SSR leafs detected, distance: " << dist << " closest points: " << cpa << " and " << cpb << endl;

      // compute dt for the primitives 
      if ((cpa - cpb).norm() > NEAR_ZERO)
      {
        double dtstar = do_CAStep(dist, cpb - cpa, a, b, ssrs.first, ssrs.second);
        if (dtstar < dt)
          dt = dtstar;
      }
      else
        return 0.0; 

      continue;
    }
    else if (ssrs.first->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, ssrs.second->children)
      {
        // cast the child
        shared_ptr<SSR> bchild = dynamic_pointer_cast<SSR>(child);

        // compute the distance between the SSRs and get closest points in 
        // global frame
        Point3d cpa, cpb;
        double dist = SSR::calc_dist(*ssrs.first, *bchild, aTb, cpa, cpb);

        // indicate closest points are in a's frame
        shared_ptr<const Pose3d> Ta = a->get_pose();
        cpa.pose = Ta;
        cpb.pose = Ta;

        // transform closest points to global frame
        cpa = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpa);
        cpb = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpb);

        // compute dt for the pair; if dt* > dt, ignore
        double dtstar = (double) 0.0;
        if ((cpa - cpb).norm() > NEAR_ZERO)
        {
          dtstar = do_CAStep(dist, cpa - cpb, b, a, bchild, ssrs.first);
          if (dtstar > dt + TTOL)
            continue;
        }

        // smaller bound..  add it for processing
        Q.push(SSRPair(ssrs.first, bchild, (double) 1.0/dtstar));
      }
    }
    else if (ssrs.second->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, ssrs.first->children)
      {
        // cast the child
        shared_ptr<SSR> achild = dynamic_pointer_cast<SSR>(child);

        // compute the distance between the SSRs and get closest points in 
        // global frame
        Point3d cpa, cpb;
        double dist = SSR::calc_dist(*achild, *ssrs.second, aTb, cpa, cpb);

        // indicate closest points are in a's frame
        shared_ptr<const Pose3d> Ta = a->get_pose();
        cpa.pose = Ta;
        cpb.pose = Ta;

        // transform closest points to global frame
        cpa = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpa);
        cpb = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpb);

        // compute dt for the pair; if dt* > dt, ignore
        double dtstar = (double) 0.0;
        if ((cpa - cpb).norm() > NEAR_ZERO)
        {
          dtstar = do_CAStep(dist, cpb - cpa, a, b, achild, ssrs.second);
          if (dtstar > dt + TTOL)
            continue;
        }

        // smaller bound..  add it for processing
        Q.push(SSRPair(achild, ssrs.second, (double) 1.0/dtstar));
      }
    }
    else
    {
      // both are internal nodes..  drill down one with the largest volume 
      if (ssrs.first->calc_volume() > ssrs.second->calc_volume())
      {
        BOOST_FOREACH(BVPtr child, ssrs.first->children)
        {
          // cast the child
          shared_ptr<SSR> achild = dynamic_pointer_cast<SSR>(child);

          // compute the distance between the SSRs and get closest points in 
          // global frame
          Point3d cpa, cpb;
          double dist = SSR::calc_dist(*achild, *ssrs.second, aTb, cpa, cpb);

          // indicate closest points are in a's frame
          shared_ptr<const Pose3d> Ta = a->get_pose();
          cpa.pose = Ta;
          cpb.pose = Ta;

          // transform closest points to global frame
          cpa = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpa);
          cpb = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpb);

          // compute dt for the pair; if dt* > dt, ignore
          double dtstar = (double) 0.0;
          if ((cpa - cpb).norm() > NEAR_ZERO)
          {
            dtstar = do_CAStep(dist, cpb - cpa, a, b, achild, ssrs.second);
            if (dtstar > dt + TTOL)
              continue;
          }

          // smaller bound..  add it for processing
          Q.push(SSRPair(achild, ssrs.second, (double) 1.0/dtstar));
        }
      }
      else
      {
        BOOST_FOREACH(BVPtr child, ssrs.second->children)
        {
          // cast the child
          shared_ptr<SSR> bchild = dynamic_pointer_cast<SSR>(child);

          // compute the distance between the SSRs and get closest points in 
          // global frame
          Point3d cpa, cpb;
          double dist = SSR::calc_dist(*ssrs.first, *bchild, aTb, cpa, cpb);

          // indicate closest points are in a's frame
          shared_ptr<const Pose3d> Ta = a->get_pose();
          cpa.pose = Ta;
          cpb.pose = Ta;

          // transform closest points to global frame
          cpa = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpa);
          cpb = Pose3d::transform_point(shared_ptr<const Pose3d>(), cpb);

          // compute dt for the pair; if dt* > dt, ignore
          double dtstar = (double) 0.0;
          if ((cpa - cpb).norm() > NEAR_ZERO)
          {
            dtstar = do_CAStep(dist, cpa - cpb, b, a, bchild, ssrs.first);
            if (dtstar > dt + TTOL)
              continue;
          }

          // smaller bound..  add it for processing
          Q.push(SSRPair(ssrs.first, bchild, (double) 1.0/dtstar));
        }
      } 
    }
  } 

  FILE_LOG(LOG_COLDET) << "C2ACCD::do_CA() exiting" << endl;
  return dt;
}

/*
/// C2A algorithm
double C2ACCD::do_C2A(CollisionGeometryPtr a, CollisionGeometryPtr b, shared_ptr<SSR> ssr_a, shared_ptr<SSR> ssr_b, const Pose3d& aTb, double& dcur, double dtcur, double w)
{
  const double DTHRESH = eps_tolerance;

  FILE_LOG(LOG_COLDET) << "C2ACCD::do_C2A() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  SSR a: " << ssr_a << "  SSR b: " << ssr_b << endl;
  FILE_LOG(LOG_COLDET) << "  dcur: " << dcur << "  dtcur: " << dtcur << "  w: " << w << endl;

  // increase w near the end of the conservative advancement
  if (dcur < DTHRESH)
    w = (double) 1.0;

  if (ssr_a->is_leaf() && ssr_b->is_leaf())
  {
    FILE_LOG(LOG_COLDET) << " *** SSR's are both leaf nodes!" << endl;
    FILE_LOG(LOG_COLDET) << " -- finding closest triangle pair..." << endl;

    // calculate the distance and closest points between ssr_a and ssr_b
    Vector3 cpa, cpb, cpa_cand, cpb_cand;
    double d = std::numeric_limits<double>::max();

    // get triangles in ssr_a and ssr_b
    assert(_meshes.find(ssr_a) != _meshes.end());
    assert(_meshes.find(ssr_b) != _meshes.end());
    pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_a = _meshes[ssr_a];
    pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_b = _meshes[ssr_b];
    assert(!mesh_a.second.empty());
    assert(!mesh_b.second.empty());

    // find closest of all pairs of triangles in ssr_a and ssr_b
    BOOST_FOREACH(unsigned ta_idx, mesh_a.second)
    {
      Triangle ta = mesh_a.first->get_triangle(ta_idx);
      BOOST_FOREACH(unsigned tb_idx, mesh_b.second)
      {
        Triangle tb = Triangle::transform(mesh_b.first->get_triangle(tb_idx), aTb);
        double dist = Triangle::calc_sq_dist(ta, tb, cpa_cand, cpb_cand);
        if (dist < d)
        {
          d = dist;
          cpa = cpa_cand;
          cpb = cpb_cand;
        }
      }
    }

    // we want actual distance, not square distance..
    d = std::sqrt(d);

    if (d < dcur)
      dcur = d;
    
    // determine cpa and cpb in the global frame, and determine dab
    cpa = a->get_pose().mult_point(cpa);
    cpb = a->get_pose().mult_point(cpb);
    Vector3 dab = cpb - cpa;

    FILE_LOG(LOG_COLDET) << " -- closest triangle pair distance: " << d << endl;
    FILE_LOG(LOG_COLDET) << " -- closest points (global frame): " << cpa << " and " << cpb << endl;

    // calculate the conservative advancement step
    double dt = do_CAStep(d, dab, a, b, ssr_a, ssr_b);

    if (dt < dtcur)
      dtcur = dt;

    FILE_LOG(LOG_COLDET) << "  -- dtcur = " << dtcur << endl;
    FILE_LOG(LOG_COLDET) << "C2ACCD::do_C2A() exiting" << endl;

    return dtcur;
  }

  vector<tuple<double, Vector3, Vector3, shared_ptr<SSR> > > distances;
  if (!ssr_a->is_leaf())
  {
    // compute the distances
    BOOST_FOREACH(BVPtr bv, ssr_a->children)
    {
      Vector3 cpa, cpb;
      shared_ptr<SSR> achild = dynamic_pointer_cast<SSR>(bv);
      double d = SSR::calc_dist(*achild, *ssr_b, aTb, cpa, cpb);
      distances.push_back(make_tuple(d, cpa, cpb, achild));
    }

    // find the smallest distance
    vector<tuple<double, Vector3, Vector3, shared_ptr<SSR> > >::const_iterator closest = std::min_element(distances.begin(), distances.end());
    if (closest->get<0>() > w*dcur)
      return do_C2A(a, b, closest->get<3>(), ssr_b, aTb, dcur, dtcur, w);
    else
    {
      Vector3 cpa = a->get_pose().mult_point(closest->get<1>());
      Vector3 cpb = a->get_pose().mult_point(closest->get<2>());
      Vector3 dab = cpb - cpa;

      // calculate conservative advancement
      double dt = do_CAStep(closest->get<0>(), dab, a, b, closest->get<3>(), ssr_b);
      if (dt < dtcur)
        dtcur = dt;
    }
  }
  else
  {
    // compute the distances
    BOOST_FOREACH(BVPtr bv, ssr_b->children)
    {
      Vector3 cpa, cpb;
      shared_ptr<SSR> bchild = dynamic_pointer_cast<SSR>(bv);
      double d = SSR::calc_dist(*ssr_a, *bchild, aTb, cpa, cpb);
      distances.push_back(make_tuple(d, cpa, cpb, bchild));
    }

    // find the smallest distance
    vector<tuple<double, Vector3, Vector3, shared_ptr<SSR> > >::const_iterator closest = std::min_element(distances.begin(), distances.end());
    if (closest->get<0>() > w*dcur)
      return do_C2A(a, b, ssr_a, closest->get<3>(), aTb, dcur, dtcur, w);
    else
    {
      Vector3 cpa = a->get_pose().mult_point(closest->get<1>());
      Vector3 cpb = a->get_pose().mult_point(closest->get<2>());
      Vector3 dab = cpb - cpa;

      // calculate conservative advancement
      double dt = do_CAStep(closest->get<0>(), dab, a, b, ssr_a, ssr_b);
      if (dt < dtcur)
        dtcur = dt;
    }
  }

  return dtcur;


  // declare nodes and distances
  shared_ptr<SSR> A, B, C, D;

  if (!ssr_a->is_leaf())
  {
    // setup nodes
    A = dynamic_pointer_cast<SSR>(ssr_a->children.front());
    B = dynamic_pointer_cast<SSR>(ssr_a->children.back());
    C = D = ssr_b;
  }
  else
  {
    // setup nodes
    A = B = ssr_a;
    C = dynamic_pointer_cast<SSR>(ssr_b->children.front());
    D = dynamic_pointer_cast<SSR>(ssr_b->children.back());
  }

  // calculate distances
  Vector3 cpa1, cpa2, cpb1, cpb2;
  double d1 = SSR::calc_dist(*A, *C, aTb, cpa1, cpb1);
  double d2 = SSR::calc_dist(*B, *D, aTb, cpa2, cpb2);
  FILE_LOG(LOG_COLDET) << " -- A: " << A << " B: " << B << " C: " << C << " D: " << D << endl;
  FILE_LOG(LOG_COLDET) << " -- d1: " << d1 << "  d2: " << d2 << endl;

  if (d2 < d1)
  {
    if (d2 > w*dcur)
      return do_C2A(a, b, B, D, aTb, dcur, dtcur, w);
    else
    {
      // determine cpa2 and cpb2 in global frame and determine dab
      cpa2 = a->get_pose().mult_point(cpa2);
      cpb2 = b->get_pose().mult_point(cpb2);
      Vector3 dab = cpb2 - cpa2;

      // calculate conservative advancement
      double dt = do_CAStep(d2, dab, a, b, B);

      if (dt < dtcur)
        dtcur = dt;
    }
    if (d1 > w*dcur)
      return do_C2A(a, b, A, C, aTb, dcur, dtcur, w);
    else
    {
      // get cpa1 and cpb1 in global frame and determine dab
      cpa1 = a->get_pose().mult_point(cpa1);
      cpb1 = a->get_pose().mult_point(cpb1);
      Vector3 dab = cpb1 - cpa1;

      // calculate conservative advancement
      double dt = do_CAStep(d1, dab, a, b, A);

      if (dt < dtcur)
        dtcur = dt;
    }
  }
  else
  {
    if (d1 > w*dcur)
      return do_C2A(a, b, A, C, aTb, dcur, dtcur, w);
    else
    {
      // get cpa1 and cpb1 in global frame and determine dab
      cpa1 = a->get_pose().mult_point(cpa1);
      cpb1 = a->get_pose().mult_point(cpb1);
      Vector3 dab = cpb1 - cpa1;

      // calculate conservative advancement
      double dt = do_CAStep(d1, dab, a, b, A);

      if (dt < dtcur)
        dtcur = dt;
    }
    if (d2 > w*dcur)
      return do_C2A(a, b, B, D, aTb, dcur, dtcur, w);
    else
    {
      // determine cpa2 and cpb2 in global frame and determine dab
      cpa2 = a->get_pose().mult_point(cpa2);
      cpb2 = a->get_pose().mult_point(cpb2);
      Vector3 dab = cpb2 - cpa2;

      // calculate conservative advancement
      double dt = do_CAStep(d2, dab, a, b, B);

      if (dt < dtcur)
        dtcur = dt;
    }
  }

  FILE_LOG(LOG_COLDET) << "  -- dtcur = " << dtcur << endl;
  FILE_LOG(LOG_COLDET) << "C2ACCD::do_C2A() exiting" << endl;

  return dtcur;
}
*/

/// Calculates the smallest distance between triangles in two SSR leaf nodes
double C2ACCD::calc_dist(shared_ptr<SSR> a, shared_ptr<SSR> b, const Transform3d& aTb, Point3d& cpa, Point3d& cpb) const
{
  FILE_LOG(LOG_COLDET) << "C2ACCD::calc_dist() entered" << endl;

  // get the meshes for the two SSRs
  assert(_meshes.find(a) != _meshes.end());
  assert(_meshes.find(b) != _meshes.end());
  const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_a = _meshes.find(a)->second;
  const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_b = _meshes.find(b)->second;

  // find closest points on triangles
  Point3d cpa_cand, cpb_cand;
  double min_dist = std::numeric_limits<double>::max();
  BOOST_FOREACH(unsigned ta_idx, mesh_a.second)
  {
    Triangle ta = mesh_a.first->get_triangle(ta_idx, aTb.target);
    FILE_LOG(LOG_COLDET) << "   ta: " << ta << endl; 
    BOOST_FOREACH(unsigned tb_idx, mesh_b.second)
    {
      Triangle utb = mesh_b.first->get_triangle(tb_idx, aTb.source);
      Triangle tb = Triangle::transform(utb, aTb);
      FILE_LOG(LOG_COLDET) << "    tb: " << tb << endl;
      double dist = Triangle::calc_sq_dist(ta, tb, cpa_cand, cpb_cand);
      FILE_LOG(LOG_COLDET) << "    sq dist: " << dist << endl;
      if (dist < min_dist)
      {
        min_dist = dist;
        cpa = cpa_cand;
        cpb = cpb_cand;
      }
    }
  }

  return std::sqrt(min_dist);
}

/// Implements Base::load_from_xml()
void C2ACCD::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "C2ACCD") == 0);

  // call parent method
  CollisionDetection::load_from_xml(node, id_map);

  // get the epsilon tolerance, if specified
  XMLAttrib* eps_tol_attr = node->get_attrib("eps-tolerance");
  if (eps_tol_attr)
    this->eps_tolerance = eps_tol_attr->get_real_value();

  // get the alpha tolerance, if specified
  XMLAttrib* alpha_tol_attr = node->get_attrib("alpha-tolerance");
  if (alpha_tol_attr)
    this->alpha_tolerance = alpha_tol_attr->get_real_value();
}

/// Implements Base::save_to_xml()
/**
 * \note neither the contact cache nor the pairs currently in collision are 
 *       saved
 */
void C2ACCD::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call parent save_to_xml() method first
  CollisionDetection::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "C2ACCD";

  // save the eps tolerance
  node->attribs.insert(XMLAttrib("eps-tolerance", eps_tolerance));

  // save the alpha tolerance
  node->attribs.insert(XMLAttrib("alpha-tolerance", alpha_tolerance));
}

/****************************************************************************
 Methods for C2A algorithm begin 
****************************************************************************/

/****************************************************************************
 Methods for C2A algorithm end 
****************************************************************************/

/****************************************************************************
 Methods for static geometry intersection testing begin 
****************************************************************************/

/// Determines whether there is a collision at the current position and orientation of the bodies
/**
 * \note the epsilon parameter is ignored
 */
bool C2ACCD::is_collision(double epsilon)
{
  // clear the set of colliding pairs and list of colliding triangles
  colliding_pairs.clear();
  colliding_tris.clear();

  // iterate over geometries 
  for (std::set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
  {
    // get the first geometry and its primitive 
    CollisionGeometryPtr g1 = *i;

    // get the transform and the inverse transform for this geometry
    shared_ptr<const Pose3d> T1 = g1->get_pose();
    Transform3d g1Tw = Pose3d::calc_relative_pose(GLOBAL, T1);

    // loop through all other geometries
    std::set<CollisionGeometryPtr>::const_iterator j = i;
    j++;
    for (; j != _geoms.end(); j++)
    {
      // get the second geometry and its primitive 
      CollisionGeometryPtr g2 = *j;

      // see whether to check
      if (!is_checked(g1, g2))
        continue; 

      // get the two BV trees
      BVPtr bv1 = _root_SSRs[g1];
      BVPtr bv2 = _root_SSRs[g2];

      // get the transform for g2 and its inverse
      shared_ptr<const Pose3d> T2 = g2->get_pose();
      Transform3d g1Tg2 = Pose3d::calc_relative_pose(T2, T1);

      // if intersects, add to colliding pairs
      if (intersect_BV_trees(bv1, bv2, g1Tg2, g1, g2))
        colliding_pairs.insert(make_sorted_pair(g1, g2));
    } 
  }

  return !colliding_pairs.empty();
}

/// Intersects two BV trees; returns <b>true</b> if one (or more) pair of the underlying triangles intersects
bool C2ACCD::intersect_BV_trees(BVPtr a, BVPtr b, const Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b) 
{
  std::queue<tuple<BVPtr, BVPtr, bool> > q;

  // get address of the last colliding triangle pair on the queue
  CollidingTriPair* last = (colliding_tris.empty()) ? NULL : &colliding_tris.back();

  FILE_LOG(LOG_COLDET) << "C2ACCD::intersect_BV_trees() entered" << endl;

  // intersect the BVs at the top level
  if (!BV::intersects(a, b, aTb))
  {
    FILE_LOG(LOG_COLDET) << "  no intersection at top-level BVs" << endl;
    FILE_LOG(LOG_COLDET) << "C2ACCD::intersect_BV_trees() exited" << endl;

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
  FILE_LOG(LOG_COLDET) << "C2ACCD::intersect_BV_trees() exited" << endl;

  return false;
} 

/*
/// Calculates the distance between two geometries as well as the closest points
double C2ACCD::calc_distance(CollisionGeometryPtr a, CollisionGeometryPtr b, Vector3& cpa, Vector3& cpb)
{
  // junk variables (variables necessary for calling a function, but whose
  // value we don't use)
  Vector3 junkA, junkB;

  // get the two bounding volumes
  shared_ptr<SSR> ssra = _root_SSRs[a];
  shared_ptr<SSR> ssrb = _root_SSRs[b];

  // determine the relative transform from b to a
  shared_ptr<const Pose3d> Ta = a->get_pose();
  shared_ptr<const Pose3d> Tb = b->get_pose();
  Transform3d aTb = Pose3d::calc_relative_pose(Tb, Ta);

  // setup the minimum distance
  double min_dist = std::numeric_limits<double>::max(); 

  // now, add the pair to the queue for processing
  queue<pair<shared_ptr<SSR>, shared_ptr<SSR> > Q;
  Q.push(make_pair(ssra, ssrb));

  // process until the queue is empty
  while (!Q.empty())
  {
    // get the distance between the bounding volumes at the front of the queue
    pair<shared_ptr<SSR>, shared_ptr<SSR> > front = Q.front();
    Q.pop();
    double dist = SSR::calc_dist(*front.first, *front.second, aTb, junkA, junkB);

    // if the distance is greater than the minimum encountered distance, don't
    // process this pair further
    if (dist > min_dist)
      continue;

    // determine which branch to go down
    if (front.first->is_leaf())
    {
      // SSR's are both leaves, find the minimum distance between triangles
      if (front.second->is_leaf())
      {
        // get the sets of triangles from both BVs
        pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_a = _meshes[front.first];
        pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_b = _meshes[front.second];

        // process all pairs of triangles
        BOOST_FOREACH(unsigned t1_idx, mesh_a.second)
        {
          Triangle t1 = mesh_a.first->get_triangle(t1_idx);
          BOOST_FOREACH(unsigned t2_idx, mesh_b.second)
          {
            Triangle t2 = mesh_b.first->get_triangle(t2_idx);
            double dist = Triangle::calc_sq_dist(t1, t2, junkA, junkB);
            if (dist < min_dist)
            {
              min_dist = dist;
              cpa = junkA;
              cpb = junkB;
            }
          }
        }
      }
      // drill down b
      else
      {
        BOOST_FOREACH(BVPtr child, front.second->children)
          Q.push(make_pair(front.first, child));
      }
    }
    else
    {
      // if SSR for b is a leaf, drill down a
      if (front.second->is_leaf())
      {
        BOOST_FOREACH(BVPtr child, front.first->children)
         Q.push(make_pair(child, front.second));
      }
      // neither is a leaf; drill down BV with greater volume
      else
      {
        double vol_a = front.first->calc_volume();
        double vol_b = front.second->calc_volume();
        if (vol_a > vol_b)
          BOOST_FOREACH(BVPtr child, front.first->children)
            Q.push(make_pair(child, front.second));
        else
          BOOST_FOREACH(BVPtr child, front.second->children)
            Q.push(make_pair(front.first, child));
      }
    }
  }
}
*/

/****************************************************************************
 Methods for static geometry intersection testing end 
****************************************************************************/


/****************************************************************************
 Methods for building bounding box trees begin 
****************************************************************************/

/// Builds an bounding box tree (OBB or AABB) from an indexed triangle mesh using a top-down approach 
/*
 * \return the root of the bounding box tree
 */
void C2ACCD::build_BV_tree(CollisionGeometryPtr geom)
{
  const unsigned THREE_D = 3;
  BVPtr child1, child2;

  FILE_LOG(LOG_BV) << "C2ACCD::build_BB_tree() entered" << endl;

  // get the intersection tolerance for the primitive
  double intersection_tolerance = geom->get_geometry()->get_intersection_tolerance();

  // get the mesh
  shared_ptr<const IndexedTriArray> mesh = geom->get_geometry()->get_mesh();

  // get the vertices from the mesh
  const vector<Origin3d>& verts = mesh->get_vertices();

  // make vertices using points
  vector<Point3d> vertices(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = Point3d(verts[i], geom->get_pose());

  // build an BV around all vertices 
  BVPtr root = BVPtr(new SSR(vertices.begin(), vertices.end()));

  // set root to point to all facet indices
  list<unsigned> tris_idx;
  for (unsigned i=0; i< mesh->get_facets().size(); i++)
    tris_idx.push_back(i);

  // setup mapping from BV to mesh
  _meshes[root] = make_pair(mesh, tris_idx);

  FILE_LOG(LOG_BV) << "  -- created root: " << root << endl;

  // add the root to a splitting queue
  queue<BVPtr> Q;
  Q.push(root);

  // split until we can't split further
  while (!Q.empty())
  {
    // get the bounding volume off of the top of the queue
    BVPtr bv = Q.front();
    Q.pop();

    FILE_LOG(LOG_BV) << "  -- splitting BV: " << bv << endl;

    // split the bounding volume across each of the three axes
    for (unsigned i=0; i< 3; i++)
    {
      Vector3d axis;

      // get the i'th column of R if an RSS
      shared_ptr<SSR> ssr = dynamic_pointer_cast<SSR>(bv);
      assert(ssr);
      ssr->R.get_column(i, axis);

      // split the bounding volume across the axis
      if (split(bv, child1, child2, axis))
        break;
    }

    // make sure that this BV was divisible
    if (!child1)
      continue;

    // setup child pointers
    bv->children.push_back(child1);
    bv->children.push_back(child2);

    // get lists of triangles for children
    assert(_meshes.find(child1) != _meshes.end());
    assert(_meshes.find(child2) != _meshes.end());
    const std::list<unsigned>& c1tris = _meshes.find(child1)->second.second;
    const std::list<unsigned>& c2tris = _meshes.find(child2)->second.second;

    // add children to the queue for processing if they have more than one tri
    if (c1tris.size() > 1)
      Q.push(child1);
    if (c2tris.size() > 1)
      Q.push(child2);

    FILE_LOG(LOG_BV) << "  -- built children " << child1 << ", " << child2 << " from " << bv << endl;
    FILE_LOG(LOG_BV) << "  -- ID: " << child1 << endl;
    FILE_LOG(LOG_BV) << "  -- ID: " << child2 << endl;
  }

  // now, collapse the tree
  Q.push(root);
  while (!Q.empty())
  {
    // for any children with a greater volume than the bv in question,
    // remove the grandchildren and add them as children
    BVPtr bv = Q.front();
    double vol = bv->calc_volume();
    bool erased_one = false;
    for (list<BVPtr>::iterator i = bv->children.begin(); i != bv->children.end(); )
    {
      // get the volume of this child
      double voli = (*i)->calc_volume();
      if (!(*i)->is_leaf() && voli > vol + NEAR_ZERO)
      {
        erased_one = true;
        BOOST_FOREACH(BVPtr gchild, (*i)->children)
          bv->children.push_back(gchild);
        i = (*i)->children.erase(i);
      }
      else
        i++;
    }

    if (!erased_one)
    {
      Q.pop();
      BOOST_FOREACH(BVPtr child, bv->children)
        if (!child->is_leaf())
          Q.push(child);
    }
  }

  Q.push(root);
  while (!Q.empty())
  {
    // get the element off of the front of the queue
    BVPtr bv = Q.front();
    Q.pop();

    // fatten the bounding volume
    // cast it to an SSR
    shared_ptr<SSR> ssr = dynamic_pointer_cast<SSR>(bv);
    assert(ssr);

    // fatten the SSR
    ssr->radius += intersection_tolerance;

    // add all children to the queue
    if (!bv->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, bv->children)
        Q.push(child);
    }

    // wipe out userdata
    bv->userdata = shared_ptr<void>();
  }

  // save the root
  _root_SSRs[geom] = dynamic_pointer_cast<SSR>(root);

  // output how many triangles are in each bounding volume
  if (LOGGING(LOG_BV))
  {
    stack<pair<BVPtr, unsigned> > S;
    S.push(make_pair(root, 0));
    while (!S.empty())
    {
      // get the node off of the top of the stack
      BVPtr node = S.top().first;
      unsigned depth = S.top().second;
      S.pop();
      
      // get the triangles in this BV
      const std::list<unsigned>& tris = _meshes.find(node)->second.second;
      std::ostringstream out;
      for (unsigned i=0; i< depth; i++)
        out << " ";
      out << "triangles covered by BV " << node << ": " << tris.size() << endl;
      FILE_LOG(LOG_BV) << out.str();

      // put all children onto the stack
      BOOST_FOREACH(BVPtr child, node->children)
        S.push(make_pair(child, depth+1));
    }
  }

  FILE_LOG(LOG_BV) << "C2ACCD::build_BV_tree() exited" << endl;
}

/// Splits a collection of triangles along a splitting plane into 2 new meshes 
void C2ACCD::split_tris(const Point3d& point, const Vector3d& normal, const IndexedTriArray& orig_mesh, const list<unsigned>& ofacets, list<unsigned>& pfacets, list<unsigned>& nfacets) 
{
  // get original vertices and facets
  const vector<Origin3d>& vertices = orig_mesh.get_vertices();
  const vector<IndexedTri>& facets = orig_mesh.get_facets();

  // determine the splitting plane: ax + by + cz = d
  double offset = Vector3d::dot(point, normal);

  // determine the side of the splitting plane of the triangles
  Plane plane(normal, offset);
  BOOST_FOREACH(unsigned i, ofacets)
  {
    // get the vertices from the facet in the same plane as the normal
    Point3d pa(vertices[facets[i].a], normal.pose);    
    Point3d pb(vertices[facets[i].b], normal.pose);    
    Point3d pc(vertices[facets[i].c], normal.pose);    

    // get the three signed distances
    double sa = plane.calc_signed_distance(pa);
    double sb = plane.calc_signed_distance(pb);
    double sc = plane.calc_signed_distance(pc);
    double min_s = std::min(sa, std::min(sb, sc));
    double max_s = std::max(sa, std::max(sb, sc));    

    // see whether we can cleanly put the triangle into one side
    if (min_s > 0)
      pfacets.push_back(i);
    else if (max_s < 0)
      nfacets.push_back(i);
    else
    {
      // triangle is split down the middle; get its centroid
      Triangle tri(pa, pb, pc);
      Point3d tri_centroid = tri.calc_centroid();
      double scent = plane.calc_signed_distance(tri_centroid);
      if (scent > 0)
        pfacets.push_back(i);
      else
        nfacets.push_back(i);
    }
  }
}

/// Splits a bounding box  along a given axis into two new bounding boxes; returns true if split successful
bool C2ACCD::split(shared_ptr<BV> source, shared_ptr<BV>& tgt1, shared_ptr<BV>& tgt2, const Vector3d& axis) 
{
  // setup two lists of triangles
  list<unsigned> ptris, ntris;

  // clear both targets
  tgt1 = shared_ptr<BV>();
  tgt2 = shared_ptr<BV>();

  // get the mesh and the list of triangles
  assert(_meshes.find(source) != _meshes.end());
  const pair<shared_ptr<const IndexedTriArray>, list<unsigned> >& mesh_data = _meshes.find(source)->second;
  shared_ptr<const IndexedTriArray> mesh = mesh_data.first;
  const list<unsigned>& tris = mesh_data.second;

  // make sure that not trying to split a single triangle
  assert(tris.size() > 1); 

  // determine the centroid of this set of triangles
  list<Triangle> t_tris;
  BOOST_FOREACH(unsigned idx, tris)
    t_tris.push_back(mesh->get_triangle(idx, source->get_relative_pose())); 
  Point3d centroid = CompGeom::calc_centroid_3D(t_tris.begin(), t_tris.end());

  // get the side of the splitting plane of the triangles
  split_tris(centroid, axis, *mesh, tris, ptris, ntris);
  if (ptris.empty() || ntris.empty())
    return false;

  // get vertices from both meshes
  vector<Origin3d> porigins, norigins;
  get_vertices(*mesh, ptris.begin(), ptris.end(), std::back_inserter(porigins));
  get_vertices(*mesh, ntris.begin(), ntris.end(), std::back_inserter(norigins));

  // setup the vertices
  vector<Point3d> pverts(porigins.size()), nverts(norigins.size());
  for (unsigned i=0; i< porigins.size(); i++)
    pverts[i] = Point3d(porigins[i], centroid.pose);
  for (unsigned i=0; i< norigins.size(); i++)
    nverts[i] = Point3d(norigins[i], centroid.pose);

  // create two new BVs 
  tgt1 = shared_ptr<SSR>(new SSR(pverts.begin(), pverts.end()));
  tgt2 = shared_ptr<SSR>(new SSR(nverts.begin(), nverts.end()));

  // setup mesh data for the BVs
  _meshes[tgt1] = make_pair(mesh, ptris);
  _meshes[tgt2] = make_pair(mesh, ntris);

  return true;
}


/****************************************************************************
 Methods for building bounding box trees end 
****************************************************************************/

/****************************************************************************
 Methods for closest triangle features determination begin 
****************************************************************************/

/// Determines whether two triangles are coplanar
bool C2ACCD::coplanar_tris(const Triangle& ta, const Triangle& tb)
{
   // if the triangles are coplanar, we have face / face
    return (CompGeom::coplanar(ta.a, ta.b, ta.c, tb.a) &&
            CompGeom::coplanar(ta.a, ta.b, ta.c, tb.b) &&
            CompGeom::coplanar(ta.a, ta.b, ta.c, tb.c));
}

/// Projects two triangles to 2D and attempts to intersect them
bool C2ACCD::project_and_intersect(const Triangle& ta, const Triangle& tb, vector<Point3d>& contact_points)
{
  const unsigned N_TRI_VERTS = 3;

  // setup useful poses
  shared_ptr<const Pose2d> GLOBAL_2D;
  shared_ptr<const Pose3d> P = ta.a.pose;

  // verify that the poses are all identical
  assert(P== tb.a.pose);
  assert(ta.b.pose == tb.b.pose);
  assert(ta.c.pose == tb.c.pose);
  assert(ta.a.pose == ta.b.pose);
  assert(P == ta.c.pose);

  // project the triangles to 2D
  Matrix3d R = CompGeom::calc_3D_to_2D_matrix(ta.calc_normal());
  double offset1 = CompGeom::determine_3D_to_2D_offset(Origin3d(ta.a), R);
  double offset2 = CompGeom::determine_3D_to_2D_offset(Origin3d(tb.a), R);
  Point2d t1[N_TRI_VERTS], t2[N_TRI_VERTS];
  for (unsigned i=0; i< N_TRI_VERTS; i++)
  {
    t1[i] = Point2d(CompGeom::to_2D(ta.get_vertex(i), R), GLOBAL_2D);
    t2[i] = Point2d(CompGeom::to_2D(tb.get_vertex(i), R), GLOBAL_2D);
  }

  // make t1 and t2 ccw
  if (!CompGeom::ccw(t1, t1+N_TRI_VERTS))
    std::swap(t1[1], t1[2]);
  if (!CompGeom::ccw(t2, t2+N_TRI_VERTS))
    std::swap(t2[1], t2[2]);

  // now attempt to intersect
  Point2d output[6];
  Point2d* end = CompGeom::intersect_tris(t1, t2, output);

  // determine contact points
  double offset = (offset1+offset2)*0.5;
  Matrix3d RT = Matrix3d::transpose(R);
  for (Point2d* begin = output; begin != end; begin++)
  {
    Origin2d o(*begin);
    contact_points.push_back(Point3d(CompGeom::to_3D(o, RT, offset), P));
  }

  return end != output;
}

/// Projects a triangle and a line segment to 2D and attempts to intersect them
bool C2ACCD::project_and_intersect(const Triangle& t, const LineSeg3& s, vector<Point3d>& contact_points)
{
  const unsigned N_TRI_VERTS = 3;
  shared_ptr<const Pose2d> GLOBAL_2D;

  // verify that poses are equal
  shared_ptr<const Pose3d> P = t.a.pose;
  assert(P == t.b.pose);
  assert(P == t.c.pose);
  assert(P == s.first.pose);
  assert(P == s.second.pose);

  // project the triangle to 2D
  Matrix3d R = CompGeom::calc_3D_to_2D_matrix(t.calc_normal());
  double offset1 = CompGeom::determine_3D_to_2D_offset(Origin3d(t.a), R);
  double offset2 = CompGeom::determine_3D_to_2D_offset(Origin3d(s.first), R);
  Point2d t0[N_TRI_VERTS];
  for (unsigned i=0; i< N_TRI_VERTS; i++)
    t0[i] = Point2d(CompGeom::to_2D(t.get_vertex(i), R), GLOBAL_2D);

  // project s to 2D
  LineSeg2 s0(Point2d(CompGeom::to_2D(s.first, R), GLOBAL_2D),
              Point2d(CompGeom::to_2D(s.second, R), GLOBAL_2D));

  // if t0 is not ccw, make it so
  if (!CompGeom::ccw(t0, t0+N_TRI_VERTS))
    std::swap(t0[1], t0[2]);

  // now attempt to intersect
  Point2d output[2];
  Point2d* end = CompGeom::intersect_seg_tri(s0, t0, output);

  // determine contact points
  double offset = (offset1+offset2)*0.5;
  Matrix3d RT = Matrix3d::transpose(R);
  for (Point2d* begin = output; begin != end; begin++)
  {
    Origin2d o(*begin);
    contact_points.push_back(Point3d(CompGeom::to_3D(o, RT, offset), P));
  }

  return end != output;
}

/// Determines whether a point lies inside the 2D projection of the triangle
bool C2ACCD::project_and_intersect(const Triangle& t, const Point3d& p)
{
  const unsigned N_TRI_VERTS = 3;
  shared_ptr<const Pose2d> GLOBAL_2D;

  // project the triangle to 2D
  Matrix3d R = CompGeom::calc_3D_to_2D_matrix(t.calc_normal());
  Point2d t0[N_TRI_VERTS];
  for (unsigned i=0; i< N_TRI_VERTS; i++)
    t0[i] = Point2d(CompGeom::to_2D(t.get_vertex(i), R), GLOBAL_2D);

  // project p to 2D
  Point2d p0(CompGeom::to_2D(p, R), GLOBAL_2D);

  // if t0 is not ccw, make it so
  if (!CompGeom::ccw(t0, t0+N_TRI_VERTS))
    std::swap(t0[1], t0[2]);

  // now see whether p0 inside t0
  return CompGeom::in_tri(t0, p0) != CompGeom::ePolygonOutside;
}

/// Determines the closest features between two triangles
void C2ACCD::determine_closest_features(const Triangle& ta, const Triangle& tb, Triangle::FeatureType& fa, Triangle::FeatureType& fb, vector<Point3d>& contact_points) const
{
  const double TOL = NEAR_ZERO;

  FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() entered" << endl;
  FILE_LOG(LOG_COLDET) << "    examining triangle: " << ta << endl;
  FILE_LOG(LOG_COLDET) << "          and triangle: " << tb << endl;

  // FACE / FACE check 
  double dot = std::fabs(ta.calc_normal().dot(tb.calc_normal()));
  if (std::fabs(dot - (double) 1.0) < TOL)
  {
    // project ta and tb to 2D for further checking
    if (project_and_intersect(ta, tb, contact_points))
    {
      fa = Triangle::eFace;
      fb = Triangle::eFace;
      FILE_LOG(LOG_COLDET) << "    face / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }

  // determine distances of vertices to faces
  double dist_ta_tba = ta.calc_signed_dist(tb.a);
  double dist_ta_tbb = ta.calc_signed_dist(tb.b);
  double dist_ta_tbc = ta.calc_signed_dist(tb.c);
  double dist_taa_tb = tb.calc_signed_dist(ta.a);
  double dist_tab_tb = tb.calc_signed_dist(ta.b);
  double dist_tac_tb = tb.calc_signed_dist(ta.c);

  FILE_LOG(LOG_COLDET) << " distance from a on triangle b to triangle a: " << dist_ta_tba << endl;
  FILE_LOG(LOG_COLDET) << " distance from b on triangle b to triangle a: " << dist_ta_tbb << endl;
  FILE_LOG(LOG_COLDET) << " distance from c on triangle b to triangle a: " << dist_ta_tbc << endl;
  FILE_LOG(LOG_COLDET) << " distance from a on triangle a to triangle b: " << dist_taa_tb << endl;
  FILE_LOG(LOG_COLDET) << " distance from b on triangle a to triangle b: " << dist_tab_tb << endl;
  FILE_LOG(LOG_COLDET) << " distance from c on triangle a to triangle b: " << dist_tac_tb << endl;

  // FACE / EDGE check
  if (std::fabs(dist_ta_tba - dist_ta_tbb) < NEAR_ZERO && std::fabs(dist_ta_tbc) > std::fabs(dist_ta_tba))
  {
    // project ta and edge (a,b) of tb to 2D for further checking
    if (project_and_intersect(ta, LineSeg3(tb.a, tb.b), contact_points))
    {
      fa = Triangle::eFace;
      fb = Triangle::eEdgeAB;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }
  if (std::fabs(dist_ta_tba - dist_ta_tbc) < NEAR_ZERO && std::fabs(dist_ta_tbb) > std::fabs(dist_ta_tba))
  {
    // project ta and edge (a,c) of tb to 2D for further checking
    if (project_and_intersect(ta, LineSeg3(tb.a, tb.c), contact_points))
    {
      fa = Triangle::eFace;
      fb = Triangle::eEdgeAC;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }
  if (std::fabs(dist_ta_tbb - dist_ta_tbc) < NEAR_ZERO && std::fabs(dist_ta_tba) > std::fabs(dist_ta_tbb))
  {
    // project ta and edge (b,c) of tb to 2D for further checking
    if (project_and_intersect(ta, LineSeg3(tb.b, tb.c), contact_points))
    {
      fa = Triangle::eFace;
      fb = Triangle::eEdgeBC;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }
  if (std::fabs(dist_taa_tb - dist_tab_tb) < NEAR_ZERO && std::fabs(dist_tac_tb) > std::fabs(dist_taa_tb))
  {
    // project tb and edge (a,b) of ta to 2D for further checking
    if (project_and_intersect(tb, LineSeg3(ta.a, ta.b), contact_points))
    {
      fa = Triangle::eEdgeAB;
      fb = Triangle::eFace;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }
  if (std::fabs(dist_taa_tb - dist_tac_tb) < NEAR_ZERO && std::fabs(dist_tab_tb) > std::fabs(dist_taa_tb))
  {
    // project tb and edge (a,c) of ta to 2D for further checking
    if (project_and_intersect(tb, LineSeg3(ta.a, ta.c), contact_points))
    {
      fa = Triangle::eEdgeAC;
      fb = Triangle::eFace;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }
  if (std::fabs(dist_tab_tb - dist_tac_tb) < NEAR_ZERO && std::fabs(dist_taa_tb) > std::fabs(dist_tab_tb))
  {
    // project tb and edge (b,c) of ta to 2D for further checking
    if (project_and_intersect(tb, LineSeg3(ta.b, ta.c), contact_points))
    {
      fa = Triangle::eEdgeBC;
      fb = Triangle::eFace;
      FILE_LOG(LOG_COLDET) << "    edge / face contact determined" << endl;
      for (unsigned i=0; i< contact_points.size(); i++)
        FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points[i] << endl;
      FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;
      return;
    }
  }

  // setup distances and faces
  vector<tuple<double, Point3d, const Triangle*, bool> > df;
  df.push_back(make_tuple(std::fabs(dist_ta_tba), tb.a, &ta, true));
  df.push_back(make_tuple(std::fabs(dist_ta_tbb), tb.b, &ta, true));
  df.push_back(make_tuple(std::fabs(dist_ta_tbc), tb.c, &ta, true));
  df.push_back(make_tuple(std::fabs(dist_taa_tb), ta.a, &tb, false));
  df.push_back(make_tuple(std::fabs(dist_tab_tb), ta.b, &tb, false));
  df.push_back(make_tuple(std::fabs(dist_tac_tb), ta.c, &tb, false));
  std::sort(df.begin(), df.end());

  // reverse so we can pop off the back
  std::reverse(df.begin(), df.end());

  // process in order of closest to the triangle
  while (!df.empty())
  {
    // get the tuple
    tuple<double, Point3d, const Triangle*, bool> tup = df.back();
    df.pop_back();
    
    // get the point and the triangle
    const Point3d& p = tup.get<1>();
    const Triangle& tri = *tup.get<2>();
    bool tri_from_A = tup.get<3>();
    if (project_and_intersect(tri, p))
    {
      // get the closest point between the triangle and p
      Point3d closest;
      Triangle::calc_sq_dist(tri, p, closest);

      fa = Triangle::eFace;
      fb = Triangle::eVertexA;
      if (!tri_from_A)
        std::swap(fa, fb);
      contact_points.push_back((p+closest)*(double) 0.5);
      FILE_LOG(LOG_COLDET) << "    vertex / face contact determined" << endl;
      FILE_LOG(LOG_COLDET) << "      contact point: " << contact_points.back() << endl;
      return;
    }
  }

  FILE_LOG(LOG_COLDET) << "  ** only degenerate contact types remaining" << endl;
  FILE_LOG(LOG_COLDET) << " -- C2ACCD::determine_closest_features() exited" << endl;

  // FACE / FACE, FACE / EDGE, FACE / VERTEX possibilities eliminated
  // now must check EDGE / EDGE, EDGE / VERTEX, VERTEX / VERTEX
/*
  // setup minimum distance
  double min_dist = std::numeric_limits<double>::max();

  // EDGE / EDGE: check each edge of ta against each edge of tb
  for (unsigned i=0; i< N_TRI_EDGES; i++)
    for (unsigned j=0; j< N_TRI_EDGES; j++)
    {
      if (parallel(ta_edges[i], tb_edges[j]) && 
          project_and_intersect(ta_edges[i], tb_edges[j]))
      {
        // get distance between the two edges
        Vector3 dummy1, dummy2;
        double sq_dist = CompGeom::calc_closest_points(ta_edges[i], tb_edges[j], dummy1, dummy2);
        if (sq_dist > min_dist)
          continue;

        // store new minimum distance
        min_dist = sq_dist;

        // assign features
        switch (i)
        {
          case 0:  fa = Triangle::eEdgeAB;  break;
          case 1:  fa = Triangle::eEdgeBC;  break;
          case 2:  fa = Triangle::eEdgeAC;  break;
        }
        switch (j)
        {
          case 0:  fb = Triangle::eEdgeAB;  break;
          case 1:  fb = Triangle::eEdgeBC;  break;
          case 2:  fb = Triangle::eEdgeAC;  break;
        }      
      }
    }

  // EDGE / VERTEX:

  // FACE / VERTEX: 
  // check each vertex of ta against tb
  Vector3 dummy;
  double sq_dist_tb_taa = CompGeom::calc_sq_dist(tb, ta.a, dummy);
  double sq_dist_tb_tab = CompGeom::calc_sq_dist(tb, ta.b, dummy);
  double sq_dist_tb_tac = CompGeom::calc_sq_dist(tb, ta.c, dummy);
  if (sq_dist_tb_taa < min_dist)
  {
    fa = Triangle::eVertexA;
    fb = Triangle::eFace;
    min_dist = sq_dist_tb_taa;
  }
  if (sq_dist_tb_tab < min_dist)
  {
    fa = Triangle::eVertexB;
    fb = Triangle::eFace;
    min_dist = sq_dist_tb_tab;
  }
  if (sq_dist_tb_tac < min_dist)
  {
    fa = Triangle::eVertexC;
    fb = Triangle::eFace;
    min_dist = sq_dist_tb_tac;
  }
  // check each vertex of tb against ta
  double sq_dist_ta_tba = CompGeom::calc_sq_dist(ta, tb.a, dummy);
  double sq_dist_ta_tbb = CompGeom::calc_sq_dist(ta, tb.b, dummy);
  double sq_dist_ta_tbc = CompGeom::calc_sq_dist(ta, tb.c, dummy);
  if (sq_dist_ta_tba < min_dist)
  {
    fa = Triangle::eFace;
    fb = Triangle::eVertexA;
    min_dist = sq_dist_ta_tba;
  }
  if (sq_dist_ta_tbb < min_dist)
  {
    fa = Triangle::eFace;
    fb = Triangle::eVertexB;
    min_dist = sq_dist_ta_tbb;
  }
  if (sq_dist_ta_tbc < min_dist)
  {
    fa = Triangle::eFace;
    fb = Triangle::eVertexC;
    min_dist = sq_dist_ta_tbc;
  }
*/
}

/****************************************************************************
 Methods for closest triangle features determination end 
****************************************************************************/
