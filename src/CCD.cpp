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

// TODO: remove this as integrator is Euler 8/11/15
/*
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
*/

/// Computes a conservative advancement step between two collision geometries assuming that velocity is constant over the interval
double CCD::calc_CA_Euler_step(const PairwiseDistInfo& pdi)
{
  // get primitives 
  PrimitivePtr pA = pdi.a->get_geometry(); 
  PrimitivePtr pB = pdi.b->get_geometry();

  // look for case of sphere
  if (dynamic_pointer_cast<SpherePrimitive>(pA) || 
      dynamic_pointer_cast<SpherePrimitive>(pB))
    return calc_CA_Euler_step_sphere(pdi);  

  // no special cases apply: call generic
  return calc_CA_Euler_step_generic(pdi);
}

/// Computes the conservative advancement time for a sphere
double CCD::calc_CA_Euler_step_sphere(const PairwiseDistInfo& pdi)
{
  // reset the minimum observed distance, if possible
  if (pdi.dist >= 0.0)
    _min_dist_observed[make_sorted_pair(pdi.a, pdi.b)] = 0.0;

  // if the distance is greater than zero, use standard conservative
  // advancement
  if (pdi.dist > NEAR_ZERO)
    return calc_CA_Euler_step_generic(pdi);

  // if the relative velocity at the point of contact is zero, return infinity
  std::vector<UnilateralConstraint> contacts;
  find_contacts(pdi.a, pdi.b, std::back_inserter(contacts), NEAR_ZERO);
  if ((contacts.size() == 1 && 
      std::fabs(contacts.front().calc_constraint_vel()) < NEAR_ZERO*10))
  {
    FILE_LOG(LOG_SIMULATOR) << "-- sphere/primitive contact with relative velocity of " << contacts.front().calc_constraint_vel() << "; reporting infinite conservative advancement time" << std::endl;
    return std::numeric_limits<double>::max();
  }
  else
  {
    if (contacts.size() >= 1)
      FILE_LOG(LOG_SIMULATOR) << "-- sphere/primitive contact with relative velocity of " << contacts.front().calc_constraint_vel() << std::endl;
  }
 
  // otherwise, use standard conservative advancement 
  return calc_CA_Euler_step_generic(pdi);
}

/// Generic method for conservative step calculation
double CCD::calc_CA_Euler_step_generic(const PairwiseDistInfo& pdi)
{
  double maxt = std::numeric_limits<double>::max();

  // reset the minimum observed distance, if possible
  if (pdi.dist >= 0.0)
    _min_dist_observed[make_sorted_pair(pdi.a, pdi.b)] = 0.0;

  // get geometries, distance, and closest points
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;
  const Point3d& pA = pdi.pa;
  const Point3d& pB = pdi.pb;

  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());
  FILE_LOG(LOG_COLDET) << "rigid body A: " << rbA->id << "  rigid body B: " << rbB->id << std::endl;

  // if the distance is (essentially) zero, do process for bodies in contact 
  if (pdi.dist <= 0.0)
    return calc_next_CA_Euler_step_generic(pdi);

  // get the direction of the vector from body B to body A
  Vector3d d0 = Pose3d::transform_point(GLOBAL, pA) -
                Pose3d::transform_point(GLOBAL, pB);
  double d0_norm = d0.norm();
  FILE_LOG(LOG_COLDET) << "closest point on " << rbA->id << ": " << Pose3d::transform_point(GLOBAL, pA) << std::endl;
  FILE_LOG(LOG_COLDET) << "closest point on " << rbB->id << ": " << Pose3d::transform_point(GLOBAL, pB) << std::endl;
  FILE_LOG(LOG_COLDET) << "distance between closest points is: " << d0_norm << std::endl;
  FILE_LOG(LOG_COLDET) << "reported distance is: " << pdi.dist << std::endl;

  // get the direction of the vector (from body B to body A)
  Vector3d n0 = d0/d0_norm;

  // if bodies are interpenetrating, reverse n0
  double dist = pdi.dist;
  if (pdi.dist < 0.0)
  {
    double& min_dist = _min_dist_observed[make_sorted_pair(pdi.a, pdi.b)];
    min_dist = std::min(min_dist, pdi.dist);
    dist = NEAR_ZERO + (pdi.dist - min_dist);
  }

  // compute the distance that body A can move toward body B
  double dist_per_tA = calc_max_dist(rbA, -n0, _rmax[cgA]);

  // compute the distance that body B can move toward body A
  double dist_per_tB = calc_max_dist(rbB, n0, _rmax[cgB]); 

  // compute the total distance
  double total_dist_per_t = dist_per_tA + dist_per_tB;
  if (total_dist_per_t < 0.0)
    total_dist_per_t = 0.0;

  FILE_LOG(LOG_COLDET) << "  distance: " << pdi.dist << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tA: " << dist_per_tA << std::endl;
  FILE_LOG(LOG_COLDET) << "  dist per tB: " << dist_per_tB << std::endl;

  // compute the maximum safe step
  maxt = std::min(maxt, dist/total_dist_per_t);

  FILE_LOG(LOG_COLDET) << "  maxt: " << maxt << std::endl;

  // return the maximum safe step
  return maxt;
}

/// Generic method for *next* conservative step calculation (only for bodies in contact)
double CCD::calc_next_CA_Euler_step_generic(const PairwiseDistInfo& pdi)
{
  const double INF = std::numeric_limits<double>::max();
  const double MIN_STEP = 1e-5;
  const double ZERO_VEL = NEAR_ZERO;

  // get the contacts
  vector<UnilateralConstraint> contacts;
  find_contacts(pdi.a, pdi.b, contacts);

  // ensure that at least one contact was found
  if (contacts.empty())
  {
//    throw std::runtime_error("No contacts found and at least one expected");
      std::cerr << "No contacts found and at least one expected" << std::endl;
      return INF;
  }

  // get the contact offset <n, x> = d
  const UnilateralConstraint& c = contacts.front();
  double d = c.contact_normal.dot(c.contact_point);

  // get the contact points
  vector<Vector3d> contact_points(contacts.size());
  for (unsigned i=0; i< contacts.size(); i++)
    contact_points[i] = contacts[i].contact_point;

  if (LOGGING(LOG_COLDET))
  {
    for (unsigned i=0; i< contact_points.size(); i++)
      FILE_LOG(LOG_COLDET) << "CCD::::calc_next_CA_Euler_step_generic() contact point: " << contact_points[i] << std::endl;
  }

  // if the constraint velocity is non-zero at a contact point, contact 
  // manifold is changing- return minimum step size (for now) 
  for (unsigned i=0; i< contacts.size(); i++)
  {
    const double CVEL = contacts[i].calc_contact_vel(contacts[i].contact_normal);
    FILE_LOG(LOG_COLDET) << "CCD::calc_next_CA_Euler_step_generic() - contact point velocity: " << CVEL << std::endl;
/*
    if (CVEL > ZERO_VEL)
      return MIN_STEP;
    else if (CVEL < -ZERO_VEL)
      return 0.0;
*/
    if (CVEL < -ZERO_VEL)
      return 0.0;
  }

  // if there are three or more contact points, contact points define a plane, 
  // and constraint velocity at all points of contact is zero, return infinity
  if (contacts.size() >= 3)
  {
    // ensure that all contact points have the same offset and normal  
    for (unsigned i=1; i< contacts.size(); i++)
    {
      // compute the dot product between the original normal and this one
      double dot = contacts[i].contact_normal.dot(c.contact_normal);
      if (false && std::fabs(dot - 1.0) > NEAR_ZERO)
      {
        std::cerr << "CCD::calc_next_CA_Euler_step_generic() - contact normal unexpectedly mis-aligned (returning minimum step)" << std::endl;
        return MIN_STEP;
      }

      // compute the offset
      double dprime = c.contact_normal.dot(contacts[i].contact_point);
      if (false && std::fabs(dprime - d) > NEAR_ZERO)
      {
        std::cerr << "CCD::calc_next_CA_Euler_step_generic() - unexpected contact offset (returning minimum step)" << std::endl;
        return MIN_STEP;
      }
    }

    FILE_LOG(LOG_COLDET) << "CCD::calc_next_CA_Euler_step_generic() - checking whether contact points are a 2-simplex" << std::endl;

    // just need to find one non-collinear set of points
    bool twosimplex = false;
    for (unsigned i=2; i< contact_points.size(); i++)
      if (!CompGeom::collinear(contact_points[0], contact_points[1],
                               contact_points[2]))
      {
        twosimplex = true;
        break;
      }

    // if it's not a two simplex, go through polygonal
    // process; otherwise, continue below
    if (twosimplex)
    {
      FILE_LOG(LOG_COLDET) << "CCD::calc_next_CA_Euler_step_generic() - contact points lie on a proper 2-simplex" << std::endl;

      return INF;
    }
  }

  // get the two geometries
  shared_ptr<CollisionGeometry> cgA = c.contact_geom1;
  shared_ptr<CollisionGeometry> cgB = c.contact_geom2; 

  // contacts do not form a plane; first get two rigid bodies 
  shared_ptr<RigidBodyd> rbA = dynamic_pointer_cast<RigidBodyd>(cgA->get_single_body());
  shared_ptr<RigidBodyd> rbB = dynamic_pointer_cast<RigidBodyd>(cgB->get_single_body());

  // now get relative velocities at centers-of-mass
  SVelocityd vA = Pose3d::transform(rbA->get_mixed_pose(), rbA->get_velocity());
  SVelocityd vB = Pose3d::transform(rbB->get_mixed_pose(), rbB->get_velocity());

  // contacts do not form a plane; many pairwise cases follow
  shared_ptr<PolyhedralPrimitive> pA = dynamic_pointer_cast<PolyhedralPrimitive>(cgA->get_geometry());
  if (pA)
  {
    // look for another polyhedron
    shared_ptr<PolyhedralPrimitive> pB = dynamic_pointer_cast<PolyhedralPrimitive>(cgB->get_geometry());
    if (pB)
    {
      // compute the relative velocity at each polyhedron's pose
      shared_ptr<const Pose3d> poseA = pA->get_pose(cgA);
      shared_ptr<const Pose3d> poseB = pB->get_pose(cgB);
      Transform3d aTb = Pose3d::calc_relative_pose(poseB, poseA);

      // get the velocity in the A pose
      SVelocityd rvA = Pose3d::transform(poseA, rbA->get_velocity()) -
                       Pose3d::transform(poseA, rbB->get_velocity());

      // now transform it to the B pose
      SVelocityd rvB = aTb.inverse_transform(rvA);

      return calc_next_CA_Euler_step_polyhedron_polyhedron(pA, pB, poseA, poseB, rvA, rvB, c.contact_normal, d);
    }

    // look for a plane
    shared_ptr<PlanePrimitive> planeB = dynamic_pointer_cast<PlanePrimitive>(cgB->get_geometry());
    if (planeB)
    {
      // compute the relative velocity at the polyhedron's pose
      shared_ptr<const Pose3d> poseA = pA->get_pose(cgA);
      SVelocityd rv = Pose3d::transform(poseA, rbA->get_velocity()) -
                      Pose3d::transform(poseA, rbB->get_velocity());

      return calc_next_CA_Euler_step_polyhedron_plane(pA, rv, poseA, c.contact_normal, d);
    }
  }

  // look for plane case
  shared_ptr<PlanePrimitive> planeA = dynamic_pointer_cast<PlanePrimitive>(cgA->get_geometry());
  if (planeA)
  {
    // look for a polyhedron
    shared_ptr<PolyhedralPrimitive> pB = dynamic_pointer_cast<PolyhedralPrimitive>(cgB->get_geometry());
    if (pB)
    {
      // compute the relative velocity at the polyhedron's pose
      shared_ptr<const Pose3d> poseB = pB->get_pose(cgB);
      SVelocityd rv = Pose3d::transform(poseB, rbA->get_velocity()) -
                      Pose3d::transform(poseB, rbB->get_velocity());

      return calc_next_CA_Euler_step_polyhedron_plane(pB, -rv, poseB, -c.contact_normal, -d);
    }
  }

  // still here? we don't know what to do- return infinity and count on
  // constraint stabilization to patch things up
  return INF;  
}

/// Computes the next step between a polyhedron and a plane
/**
 * \param nA the contact normal in frame A, given that polyhedron is defined in frame A
 * \param offset the offset of the contact plane
 */
double CCD::calc_next_CA_Euler_step_polyhedron_plane(shared_ptr<PolyhedralPrimitive> p, const SVelocityd& rv, shared_ptr<const Pose3d> P, const Vector3d& normal, double offset0)
{
  const double INF = std::numeric_limits<double>::max();
  double max_step = INF;

  // for illustration purposes, assume that polyhedron is body A, plane is
  // body B, contact normal is [0 1 0], and A's linear velocity is [0 -1 0]

  // get the polyhedron and its vertices
  const Polyhedron& poly = p->get_polyhedron();
  const std::vector<shared_ptr<Polyhedron::Vertex> >& v = poly.get_vertices();

  // get the normal in P's frame
  Vector3d nP = Pose3d::transform_vector(P, normal);

  // get the new offset
  Vector3d p0 = normal * offset0;
  const double offset = nP.dot(Pose3d::transform_point(P, p0));

  // get the Euclidean norm of the angular velocity
  double av_norm = rv.get_angular().norm();

  // compute the dot product of the relative linear velocity and the normal
  double lv_dot_n = -nP.dot(rv.get_linear());

  FILE_LOG(LOG_COLDET) << "CCD::calc_next_CA_Euler_step_polyhedron_plane(): contact normal: " << normal << " offset: " << offset0 << std::endl;  
  FILE_LOG(LOG_COLDET) << "  lv'n: " << lv_dot_n << std::endl;  

  // loop over all vertices of the polyhedron (assumed to be in P frame)
  for (unsigned i=0; i< v.size(); i++)
  {
    // get the distance of the vertex from the origin
    double r = v[i]->o.norm();

    // setup the vertex
    Vector3d vertex(v[i]->o, P);

    // compute the distance from the vertex to the contact plane <n, x> = d
    double dist = nP.dot(vertex) - offset;
    FILE_LOG(LOG_COLDET) << "vertex: " << Pose3d::transform_point(GLOBAL, vertex) << " distance: " << dist << " speed: " << std::max(0.0, lv_dot_n + av_norm*r) << "  max step: " << dist/std::max(0.0, lv_dot_n + av_norm*r) << std::endl;  

    // if the distance is effectively zero, ignore the vertex
    if (dist < NEAR_ZERO)
      continue;

    // compute the speed of the vertex
    double speed = std::max(0.0, lv_dot_n + av_norm*r);

    // divide the distance by the maximum speed of that vertex   
    max_step = std::min(max_step, dist/speed);
  }

  return max_step;
}

/// Computes the next step between two polyhedra
/**
 * \param nA the contact normal in frame A, given that polyhedron is defined in frame A
 * \param offset the offset of the contact plane
 */
double CCD::calc_next_CA_Euler_step_polyhedron_polyhedron(shared_ptr<PolyhedralPrimitive> pA, shared_ptr<PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB, const SVelocityd& rvA, const SVelocityd& rvB, const Vector3d& n0, double offset0)
{
  const double INF = std::numeric_limits<double>::max();
  double max_step = INF;

  // get the normal in A and B's frames
  Vector3d nA = Pose3d::transform_vector(poseA, n0);
  Vector3d nB = Pose3d::transform_vector(poseB, -n0);

  // compute offsets w.r.t. nA/nB frame
  Vector3d p0 = n0 * offset0;
  const double offsetA = nA.dot(Pose3d::transform_point(poseA, p0));
  const double offsetB = nB.dot(Pose3d::transform_point(poseB, p0));

  // get the polyhedra and vertices
  const Polyhedron& polyA = pA->get_polyhedron();
  const Polyhedron& polyB = pB->get_polyhedron();
  const std::vector<shared_ptr<Polyhedron::Vertex> >& vA = polyA.get_vertices();
  const std::vector<shared_ptr<Polyhedron::Vertex> >& vB = polyB.get_vertices();

  // get the Euclidean norms of the angular velocities
  double avA_norm = rvA.get_angular().norm();
  double avB_norm = rvB.get_angular().norm();

  // compute the dot products of the relative linear velocities and the normals
  double lvA_dot_n = -nA.dot(rvA.get_linear());
  double lvB_dot_n = nB.dot(rvB.get_linear());

  // loop over all vertices of polyhedron A 
  for (unsigned i=0; i< vA.size(); i++)
  {
    // get the distance of the vertex from the origin
    double r = vA[i]->o.norm();

    // get the vertex in the pose
    Vector3d vertex(vA[i]->o, poseA);

    // compute the distance from the vertex to the contact plane <n, x> = d
    double dist = nA.dot(vertex) - offsetA;

    // if the distance is effectively zero, ignore the vertex
    if (dist < NEAR_ZERO)
      continue;

    // compute the speed of the vertex
    double speed = std::max(0.0, lvA_dot_n + avA_norm*r);

    // divide the distance by the maximum speed of that vertex   
    max_step = std::min(max_step, dist/speed);
  }

  // loop over all vertices of polyhedron B 
  for (unsigned i=0; i< vB.size(); i++)
  {
    // get the distance of the vertex from the origin
    double r = vB[i]->o.norm();

    // get the vertex in the pose
    Vector3d vertex(vB[i]->o, poseB);

    // compute the distance from the vertex to the contact plane <n, x> = d
    double dist = nB.dot(vertex) - offsetB;

    // if the distance is effectively zero, ignore the vertex
    if (dist < NEAR_ZERO)
      continue;

    // compute the speed of the vertex
    double speed = std::max(0.0, lvB_dot_n + avB_norm*r);

    // divide the distance by the maximum speed of that vertex   
    max_step = std::min(max_step, dist/speed);
  }

  return max_step;
}

/// Computes the maximum velocity along a particular direction (n)
double CCD::calc_max_dist(ArticulatedBodyPtr ab, RigidBodyPtr rb, const Vector3d& n, double rmax)
{
  // get the base link
  RigidBodyPtr base = dynamic_pointer_cast<RigidBody>(ab->get_base_link());

  // get the base link's velocity
  const SVelocityd& base_v0 = Pose3d::transform(GLOBAL, base->get_velocity());
  Vector3d base_xd0 = base_v0.get_linear();

  // setup the initial movement
  double mvmt = n.dot(base_xd0);

  // get the inner joint
  JointPtr inner = rb->get_inner_joint_explicit();

  // add the movement in for the rigid body
  mvmt += 2.0 * rmax * inner->qd.norm();

  // get the joint pose
  Pose3d joint_pose = *inner->get_pose();
  joint_pose.update_relative_pose(GLOBAL);

  // keep looping until we arrive at the base link
  while (true)
  {
    rb = inner->get_inboard_link();
    if (rb == base)
      break;
    JointPtr next_inner = rb->get_inner_joint_explicit();
    Pose3d next_joint_pose = *next_inner->get_pose();
    next_joint_pose.update_relative_pose(GLOBAL);
    mvmt += next_inner->qd.norm() * (next_joint_pose.x - joint_pose.x).norm(); 
    joint_pose = next_joint_pose;
    inner = next_inner;
  }

  return mvmt; 
}

/// Computes the maximum velocity along a particular direction (n)
double CCD::calc_max_dist(RigidBodyPtr rb, const Vector3d& n, double rmax)
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (!rb->is_enabled())
    return 0.0;

  // if the body is part of an articulated body, do that calculation instead
  if (rb->get_articulated_body() && rb->get_articulated_body()->get_base_link() != rb)
    return calc_max_dist(dynamic_pointer_cast<ArticulatedBody>(rb->get_articulated_body()), rb, n, rmax);

  // get the velocities at t0
  const SVelocityd& v0 = Pose3d::transform(GLOBAL, rb->get_velocity());
  Vector3d xd0 = v0.get_linear();
  Vector3d w0 = v0.get_angular();
  FILE_LOG(LOG_COLDET) << "CCD::calc_max_velocity() called on " << rb->id << std::endl;
  FILE_LOG(LOG_COLDET) << "  n = " << n << std::endl;
  FILE_LOG(LOG_COLDET) << "  xd0 = " << xd0 << std::endl;
  FILE_LOG(LOG_COLDET) << "  w0 = " << w0 << std::endl;
  FILE_LOG(LOG_COLDET) << "  <n, xd0> = " << n.dot(xd0) << std::endl;
  FILE_LOG(LOG_COLDET) << "  ||w0 x n|| * r = " << Vector3d::cross(w0, n).norm()*rmax << std::endl;
  return n.dot(xd0) + Vector3d::cross(w0, n).norm()*rmax;
}

// NOTE: migration of dynamics computations to Ravelin break this function; slated for removal on 8/10/15
/*
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
*/

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
 Methods for broad phase begin
****************************************************************************/
void CCD::broad_phase(double dt, const vector<ControlledBodyPtr>& bodies, vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
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
    {
      BOOST_FOREACH(shared_ptr<RigidBodyd> rb, ab->get_links())
        rbs.push_back(dynamic_pointer_cast<RigidBody>(rb));
    }
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
    FILE_LOG(LOG_COLDET) << i->second << " overlaps between " << i->first.first << " (" << i->first.first->get_single_body()->body_id << ") and " << i->first.second << " (" << i->first.second->get_single_body()->body_id << ")" << std::endl;

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

  // compute the swept BV
  BVPtr swept_bv = bv->calc_swept_BV(cg, v*dt);
  FILE_LOG(LOG_BV) << "new BV: " << swept_bv << std::endl;

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
    FILE_LOG(LOG_COLDET) << "  updating collision geometry: " << geom << "  rigid body: " << geom->get_single_body()->body_id << std::endl;

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

  // look for torus primitive (also an easy case)
  shared_ptr<TorusPrimitive> torus_p = dynamic_pointer_cast<TorusPrimitive>(p);
  if (torus_p)
  {
    sph->radius = torus_p->get_major_radius() + torus_p->get_minor_radius();
    return sph;
  }

  // look for box primitive (also an easy case)
  shared_ptr<BoxPrimitive> box_p = dynamic_pointer_cast<BoxPrimitive>(p);
  if (box_p)
  {
    sph->radius = Origin3d(box_p->get_x_len()/2.0, box_p->get_y_len()/2.0, box_p->get_z_len()/2.0).norm();
    return sph;
  }

  // look for generic polyhedral primitive
  shared_ptr<PolyhedralPrimitive> pp = dynamic_pointer_cast<PolyhedralPrimitive>(p);
  if (pp)
  {
    sph->radius = pp->get_bounding_radius();
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


