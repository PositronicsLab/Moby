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
#include <Ravelin/PrismaticJointd.h>
#include <Ravelin/SphericalJointd.h>
#include <Ravelin/UniversalJointd.h>
#include <Ravelin/RevoluteJointd.h>
#include <Ravelin/FixedJointd.h>
#include <Moby/Constraint.h>
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
#include <Moby/ConstraintSimulator.h>
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

/// Gets a counter-clockwise oriented triangle from a Polyhedron::Face
Triangle CCD::get_triangle(shared_ptr<Polyhedron::Face> f, const Transform3d& wTf)
{
  shared_ptr<Polyhedron::Vertex> v;

  // setup a counter clockwise iterator
  Polyhedron::VertexFaceIterator vfi(f, true);

  // get first vertex 
  assert(vfi.has_next());
  v = *vfi;
  Point3d a = wTf.transform_point(Point3d(v->o, wTf.source)); 

  // get second vertex
  vfi.advance();
  v = *vfi;
  Point3d b = wTf.transform_point(Point3d(v->o, wTf.source)); 

  // get third vertex
  assert(vfi.has_next());
  vfi.advance();
  v = *vfi;
  Point3d c = wTf.transform_point(Point3d(v->o, wTf.source));   

  // verify that the triangle has only three vertices
  assert(!vfi.has_next());

  // setup the triangle
  return Triangle(a, b, c);
}

/// Conversative advancement approach for when two bodies are touching
double CCD::calc_next_CA_Euler_step(const PairwiseDistInfo& pdi, double epsilon)
{
  // get the primitive types
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;
  PrimitivePtr pA = cgA->get_geometry();
  PrimitivePtr pB = cgB->get_geometry();

  // look for polyhedron vs. half-space or vs. polyhedron
  if (boost::dynamic_pointer_cast<PlanePrimitive>(pA))
  {
    // shouldn't be a half-space but you never know...
    if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return 0.0;
    else if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pB))
    {
      PairwiseDistInfo pdi2;
      pdi2.a = pdi.b;
      pdi2.b = pdi.a;
      pdi2.dist = pdi.dist;
      pdi2.pa = pdi.pb;
      pdi2.pb = pdi.pa;
      return calc_next_CA_Euler_step_polyhedron_plane(pdi2, epsilon);
    }
  }
  else if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pA))
  {
    // look for polyhedron vs. half-space
    if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return calc_next_CA_Euler_step_polyhedron_plane(pdi, epsilon);
    else if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pB))
      return calc_next_CA_Euler_step_polyhedron_polyhedron(pdi, epsilon);
  }

  // still here? Bodies are too close- return 0
  return 0.0;
}

/// Computes the Euler step for a polyhedron vs. a plane
double CCD::calc_next_CA_Euler_step_polyhedron_plane(const PairwiseDistInfo& pdi, double epsilon)
{
  const unsigned Y = 1;
  vector<unsigned> contact_points, high_points;
  vector<double> dists;

  // get the primitive types
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;

  // get the polyhedron primitive
  shared_ptr<PolyhedralPrimitive> ppoly = dynamic_pointer_cast<PolyhedralPrimitive>(cgA->get_geometry());

  // get the plane primitive 
  shared_ptr<PlanePrimitive> pplane = dynamic_pointer_cast<PlanePrimitive>(cgB->get_geometry());

  // get the two poses
  boost::shared_ptr<const Ravelin::Pose3d> poly_pose = ppoly->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> plane_pose = pplane->get_pose(cgB);

  // get the transforms from polyhedron pose to plane pose and from plane
  // pose to global frame
  Transform3d planeTpoly = Pose3d::calc_relative_pose(poly_pose, plane_pose); 
  Transform3d wTplane = Pose3d::calc_relative_pose(plane_pose, GLOBAL); 

  // get the polyhedron
  const Polyhedron& poly = ppoly->get_polyhedron();

  // get all vertices from the polyhedron
  const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& v = poly.get_vertices();

  // setup minimum distance 
  double min_dist = std::numeric_limits<double>::max(); 

  // get the normal to the plane in the global frame
  Vector3d y_vec(0.0, 1.0, 0.0, plane_pose);
  Vector3d n0 = wTplane.transform_vector(y_vec);

  // prepare for separating bodies
  bool separating = false;
  double max_step = std::numeric_limits<double>::max();

  // process all vertices from the polyhedron
  for (unsigned i=0; i< v.size(); i++)
  {
    // transform the vertex to the half-space pose
    Point3d p = planeTpoly.transform_point(Point3d(v[i]->o, poly_pose));

    // get the distance from the plane
    double dist = p[Y]; 

    // if the distance is less than the minimum, save it
    if (dist < min_dist)
      min_dist = dist;
    
    // if distance is within tolerance, it's a contact point- measure the
    // current projected velocity
    if (dist < epsilon)
    {
      // get the contact in the global frame
      Point3d p0 = wTplane.transform_point(p);

      // create a unilateral constraint
      Constraint c = create_contact(cgA, cgB, p0, n0, dist);

      // output the velocity
      FILE_LOG(LOG_CONSTRAINT) << "Contact point found: " << p0 << " velocity: " << c.calc_constraint_vel(0) << std::endl;

      // if velocity at a point of contact is greater than zero, return 0.0
      // (let simulator take minimum size) 
      double cv = c.calc_constraint_vel(0);
      if (cv > NEAR_ZERO)
      {
        separating = true;
        max_step = std::min(max_step, epsilon - dist)/cv; 
        continue;
      }

      // save the index of the contact point
      contact_points.push_back(i);
    }
    // else the vertex is away from the plane
    else
    {
      // save the distance to the plane, and the index of the point
      dists.push_back(dist);
      high_points.push_back(i);
    }
  }

  // look for separating condition
  if (separating)
    return max_step;

  // get the rigid body corresponding to the polyhedron
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  
  // at least one vertex must be in contact, or why are we here? 
  assert(min_dist < epsilon);

  // get the angular velocity vector in the half-space pose
  Vector3d omega = Pose3d::transform(plane_pose, rb->get_velocity()).get_angular();  

  // setup the divisor
  double divisor = 0.0;

  // if rb is part of an articulated body
  shared_ptr<ArticulatedBodyd> ab = rb->get_articulated_body();
  if (ab)
  {
    double omega = 0.0;

    // compute the angular velocity of every joint up to the base
    shared_ptr<Jointd> inner = rb->get_inner_joint_explicit();
    while (inner)
    {
      omega += inner->qd.norm();
      shared_ptr<RigidBodyd> inboard = inner->get_inboard_link();
      if (inboard->is_base() && ab->is_floating_base())
        omega += inboard->get_velocity().get_angular().norm();
      inner = inboard->get_inner_joint_explicit();
    }

    // update the divisor
    divisor += omega;
  }
  else
  {
    // else, get the angular velocity vector in the global frame
    Vector3d omega = Pose3d::transform(GLOBAL, rb->get_velocity()).get_angular();  
  
    // update the divisor with the cross product of the norm and the angular
    // velocity
    divisor += Vector3d::cross(omega, n0).norm();
  }
 
  // setup a minimum time
  double min_t = std::numeric_limits<double>::max();

  // loop through all vertices above the half-space
  for (unsigned i=0; i< high_points.size(); i++)
  {
    // setup a maximum distance
    double max_dist = 0.0;

    // loop through all contact points
    for (unsigned j=0; j< contact_points.size(); j++)
    {
      // get the distance between the two points
      double dist = (v[high_points[i]]->o - v[contact_points[j]]->o).norm();
      if (dist > max_dist)
        max_dist = dist; 
    }

    FILE_LOG(LOG_CONSTRAINT) << "Distance from adjacent vertex " << v[high_points[i]]->o << " to contact point: " << max_dist << std::endl;

    // compute the conservative time
    double t = dists[i]/max_dist;
    FILE_LOG(LOG_CONSTRAINT) << "t for this point: " << (t / divisor) << std::endl;
    if (t < min_t)
      min_t = t;
  }

  // divide min_t by the cross term
  min_t /= divisor; 

  return min_t;
}

/// Computes the next Euler step for a polyhedron vs. a polyhedron 
double CCD::calc_next_CA_Euler_step_polyhedron_polyhedron(const PairwiseDistInfo& pdi, double epsilon)
{
  vector<unsigned> contact_pointsA, contact_pointsB;
  vector<unsigned> high_pointsA, high_pointsB;
  vector<double> distsA, distsB;
  const unsigned Y = 1;

  // get the primitive types
  CollisionGeometryPtr cgA = pdi.a; 
  CollisionGeometryPtr cgB = pdi.b;

  // get the polyhedron primitives
  shared_ptr<PolyhedralPrimitive> ppA = dynamic_pointer_cast<PolyhedralPrimitive>(cgA->get_geometry());
  shared_ptr<PolyhedralPrimitive> ppB = dynamic_pointer_cast<PolyhedralPrimitive>(cgB->get_geometry());

  // get the two poses
  boost::shared_ptr<const Ravelin::Pose3d> poseA = ppA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> poseB = ppB->get_pose(cgB);

  // get the transforms from two polyhedra to global frame 
  // pose to global frame
  Transform3d wTA = Pose3d::calc_relative_pose(poseA, GLOBAL); 
  Transform3d wTB = Pose3d::calc_relative_pose(poseB, GLOBAL); 

  // get the polyhedra
  const Polyhedron& polyA = ppA->get_polyhedron();
  const Polyhedron& polyB = ppB->get_polyhedron();

  // get all vertices from the polyhedra
  const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vA = polyA.get_vertices();
  const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vB = polyB.get_vertices();

  // get last contacts for these geometries
  bool normal_set = false;
  Vector3d n0;  // plane defined by <n0, x> - d = 0
  double d;
  shared_ptr<ConstraintSimulator> sim(_simulator);
  vector<Constraint>& c = sim->impact_constraint_handler.last_constraints;

  // loop may have to go multiple times
  double eps_plus = epsilon; 
  for (unsigned j=0; j< 10; j++, eps_plus *= 10.0)
  {
    // loop through 
    for (unsigned i=0; i< c.size(); i++)
    {
      // only look for contact constraints
      if (c[i].constraint_type != Constraint::eContact)
        continue;

      // look for contact between the two geometries
      if (c[i].contact_geom1 == cgA && c[i].contact_geom2 == cgB)
      { 
        // if the normal has already been set, verify it's the same
        if (normal_set)
          assert(std::fabs(n0.dot(c[i].contact_normal) - 1.0) < NEAR_ZERO);
        // otherwise, save the normal and compute the point on the plane
        else
        {
          n0 = c[i].contact_normal;
          d = n0.dot(c[i].contact_point);
          normal_set = true;
        }
      }
      else if (c[i].contact_geom1 == cgB && c[i].contact_geom2 == cgA)
      {
        // if the normal has already been set, verify it's the same
        if (normal_set)
          assert(std::fabs(n0.dot(-c[i].contact_normal) - 1.0) < NEAR_ZERO);
        // otherwise, save the reverse of the normal
        else
        {
          n0 = -c[i].contact_normal;
          d = n0.dot(c[i].contact_point);
          normal_set = true;
        }
      }
    }

    // if the normal has not been set, compute constraints
    if (!normal_set)
    {
      // compute constraints
      find_contacts(cgA, cgB, std::back_inserter(c), eps_plus);
    
      // redo the calculations
      continue; 
    }
    else
      break;
  }

  // verify that the normal has been set
  if (!normal_set)
  {
    FILE_LOG(LOG_CONSTRAINT) << "CCD::calc_next_CA_Euler_step_polyhedron_polyhedron()- warning! interpenetration indicated but no contact points found!" << std::endl;
    return std::numeric_limits<double>::max();
  }

  // prepare for separating bodies
  bool separating = false;
  double max_step = std::numeric_limits<double>::max();

  // process all vertices from the first polyhedron
  for (unsigned i=0; i< vA.size(); i++)
  {
    // get the vertex in the global frame 
    Point3d p0 = wTA.transform_point(Point3d(vA[i]->o, poseA));

    // get the distance from the plane
    double dist = p0.dot(n0) - d; 

    // if distance is within tolerance, it's a contact point- measure the
    // current projected velocity
    if (dist < epsilon)
    {
      // create a unilateral constraint
      Constraint c = create_contact(cgA, cgB, p0, n0, dist);

      // if velocity at a point of contact is greater than zero, return 0.0
      // (let simulator take minimum size) 
      double cv = c.calc_constraint_vel(0);
      FILE_LOG(LOG_CONSTRAINT) << "Contact point found: " << p0 << " velocity: " << cv << std::endl;
      if (cv > NEAR_ZERO)
      {
        separating = true;
        max_step = std::min(max_step, epsilon - dist)/cv; 
        continue;
      }

      // save the index of the contact point
      contact_pointsA.push_back(i);
    }
    // else the vertex is away from the plane
    else
    {
      // save the distance to the plane, and the index of the point
      distsA.push_back(dist);
      high_pointsA.push_back(i);
    }
  }

  // process all vertices from the second polyhedron
  for (unsigned i=0; i< vB.size(); i++)
  {
    // get the vertex in the global frame 
    Point3d p0 = wTB.transform_point(Point3d(vB[i]->o, poseB));

    // get the distance from the plane
    double dist = p0.dot(n0) - d; 

    // if distance is within tolerance, it's a contact point- measure the
    // current projected velocity
    if (dist < epsilon)
    {
      // create a unilateral constraint
      Constraint c = create_contact(cgA, cgB, p0, n0, dist);

      // if velocity at a point of contact is greater than zero, return 0.0
      // (let simulator take minimum size) 
      double cv = c.calc_constraint_vel(0);
      FILE_LOG(LOG_CONSTRAINT) << "Contact point found: " << p0 << " velocity: " << cv << std::endl;
      if (cv > NEAR_ZERO)
      {
        separating = true;
        max_step = std::min(max_step, epsilon - dist)/cv; 
        continue;
      }

      // save the index of the contact point
      contact_pointsB.push_back(i);
    }
    // else the vertex is away from the plane
    else
    {
      // save the distance to the plane, and the index of the point
      distsB.push_back(dist);
      high_pointsB.push_back(i);
    }
  }

  // look for separating condition
  if (separating)
    return max_step;

  // get the rigid bodies corresponding to the two polyhedra
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());
  
  // setup the divisor
  double divisor = 0.0;

  // if rbA is part of an articulated body
  shared_ptr<ArticulatedBodyd> abA = rbA->get_articulated_body();
  if (abA)
  {
    double omega = 0.0;

    // compute the angular velocity of every joint up to the base
    shared_ptr<Jointd> inner = rbA->get_inner_joint_explicit();
    while (inner)
    {
      omega += inner->qd.norm();
      shared_ptr<RigidBodyd> inboard = inner->get_inboard_link();
      if (inboard->is_base() && abA->is_floating_base())
        omega += inboard->get_velocity().get_angular().norm();
      else
        inner = inboard->get_inner_joint_explicit();
    }

    // update the divisor
    divisor += omega;
  }
  else
  {
    // else, get the angular velocity vector in the global frame
    Vector3d omegaA = Pose3d::transform(GLOBAL, rbA->get_velocity()).get_angular();  
  
    // update the divisor with the cross product of the norm and the angular
    // velocity
    divisor += Vector3d::cross(omegaA, n0).norm();
  }

  // if rbB is part of an articulated body
  shared_ptr<ArticulatedBodyd> abB = rbB->get_articulated_body();
  if (abA)
  {
    double omega = 0.0;

    // compute the angular velocity of every joint up to the base
    shared_ptr<Jointd> inner = rbB->get_inner_joint_explicit();
    while (inner)
    {
      omega += inner->qd.norm();
      shared_ptr<RigidBodyd> inboard = inner->get_inboard_link();
      if (inboard->is_base() && abB->is_floating_base())
        omega += inboard->get_velocity().get_angular().norm();
      else
        inner = inboard->get_inner_joint_explicit();
    }

    // uddate the divisor
    divisor += omega;
  }
  else
  {
    // else, get the angular velocity vector in the global frame
    Vector3d omegaB = Pose3d::transform(GLOBAL, rbB->get_velocity()).get_angular();  
  
    // update the divisor with the cross product of the norm and the angular
    // velocity
    divisor += Vector3d::cross(omegaB, n0).norm();
  }

  // setup a minimum time
  double min_t = std::numeric_limits<double>::max();

  // loop through all vertices from A above the half-space
  for (unsigned i=0; i< high_pointsA.size(); i++)
  {
    // setup a maximum distance
    double max_dist = 0.0;

    // loop through all contact points
    for (unsigned j=0; j< contact_pointsA.size(); j++)
    {
      // get the distance between the two points
      double dist = (vA[high_pointsA[i]]->o - vA[contact_pointsA[j]]->o).norm();
      FILE_LOG(LOG_CONSTRAINT) << "Distance from adjacent vertex " << vA[high_pointsA[i]]->o << " to contact point: " << dist << std::endl;
      if (dist > max_dist)
        max_dist = dist; 
    }

    // compute the conservative time
    double t = distsA[i]/max_dist;
    FILE_LOG(LOG_CONSTRAINT) << "t for this point: " << (t / divisor) << std::endl;
    if (t < min_t)
      min_t = t;
  }

  // loop through all vertices from B above the half-space
  for (unsigned i=0; i< high_pointsB.size(); i++)
  {
    // setup a maximum distance
    double max_dist = 0.0;

    // loop through all contact points
    for (unsigned j=0; j< contact_pointsB.size(); j++)
    {
      // get the distance between the two points
      double dist = (vB[high_pointsB[i]]->o - vB[contact_pointsB[j]]->o).norm();
      if (dist > max_dist)
        max_dist = dist; 

      FILE_LOG(LOG_CONSTRAINT) << "Distance from adjacent vertex " << vB[high_pointsB[i]]->o << " to contact point: " << dist << std::endl;
    }

    // compute the conservative time
    double t = distsB[i]/max_dist;
    FILE_LOG(LOG_CONSTRAINT) << "t for this point: " << (t / divisor) << std::endl;
    if (t < min_t)
      min_t = t;
  }

  // compute the real min_t 
  min_t /= divisor; 

  return min_t;
}

/// Method for conservative step calculation
/**
 * \param pdi the distance information for a pair of geometries
 * \param epsilon the (small) rigid boundary to maintain between geometries
 */
double CCD::calc_CA_Euler_step(const PairwiseDistInfo& pdi, double epsilon)
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

  // get the compliant layer depth
  double delta = cgA->compliant_layer_depth + cgB->compliant_layer_depth;

  // get the two underlying bodies
  RigidBodyPtr rbA = dynamic_pointer_cast<RigidBody>(cgA->get_single_body());
  RigidBodyPtr rbB = dynamic_pointer_cast<RigidBody>(cgB->get_single_body());
  FILE_LOG(LOG_COLDET) << "rigid body A: " << rbA->id << "  rigid body B: " << rbB->id << std::endl;

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

  // get the current distance 
  double dist = pdi.dist;
  if (dist < epsilon)
    return calc_next_CA_Euler_step(pdi, epsilon);

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
  if (dist > delta+NEAR_ZERO)
    maxt = std::min(maxt, (dist-delta)/total_dist_per_t);
  else
    maxt = std::min(maxt, dist/total_dist_per_t);

  FILE_LOG(LOG_COLDET) << "  maxt: " << maxt << std::endl;

  // return the maximum safe step
  return maxt;
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
  static map<ArticulatedBodyPtr, vector<double> > link_length_map;

  // get the base link
  RigidBodyPtr base = dynamic_pointer_cast<RigidBody>(ab->get_base_link());

  // get the explicit joints in the articulated body
  const vector<shared_ptr<Jointd> >& joints = ab->get_explicit_joints(); 

  // see whether the articulated body has been processed into the link length
  // map already
  vector<double>& link_lengths = link_length_map[ab];
  if (link_lengths.size() != joints.size())
  {
    // resize link lengths
    link_lengths.resize(joints.size());

    // store the multibody's current configuration
    VectorNd q, qzero;
    ab->get_generalized_coordinates_euler(q);

    // set the multibody to a zero configuration
    qzero.set_zero(q.size());
    ab->set_generalized_coordinates_euler(qzero);

    // process all joints 
    for (unsigned i=0; i< joints.size(); i++)
    {
      // get the joint index
      unsigned j = joints[i]->get_index();

      // get the joint pose in the global frame
      Pose3d joint_pose = *joints[i]->get_pose(); 
      joint_pose.update_relative_pose(GLOBAL);

      // get the inboard link; if the inboard link is the base, compute the
      // distance from the joint to the base c.o.m. 
      shared_ptr<RigidBodyd> inboard = joints[i]->get_inboard_link();
      if (inboard == ab->get_base_link())
      {
        // update the length from the floating-base
        Pose3d base_x = *base->get_pose();
        base_x.update_relative_pose(GLOBAL);
        link_lengths[j] = (joint_pose.x - base_x.x).norm(); 
        continue;
      }

      // get the preceding joint
      shared_ptr<Jointd> prec_joint = inboard->get_inner_joint_explicit();

      // get the poses for both joints
      Pose3d prec_joint_pose = *prec_joint->get_pose(); 
      prec_joint_pose.update_relative_pose(GLOBAL);

      // get the distance between the joint positions 
      link_lengths[j] = (joint_pose.x - prec_joint_pose.x).norm(); 
    }

    // restore the multibody's previous configuration
    ab->set_generalized_coordinates_euler(q);
  }

  // get the base link's velocity
  const SVelocityd& base_v = Pose3d::transform(base->get_gc_pose(), base->get_velocity());
  Vector3d base_xd = base_v.get_linear();
  Vector3d base_omega = base_v.get_angular();

  // limit the base omega
  double base_omega_nrm = base_omega.norm();
  if (base_omega_nrm > 2.0*M_PI)
    base_omega *= 2.0*M_PI/base_omega_nrm;

  // transform n
  Vector3d n0 = Pose3d::transform_vector(base->get_gc_pose(), n);

  // setup the initial movement
  double mvmt = std::fabs(n0.dot(base_xd));

  // get the inner joint
  JointPtr inner = rb->get_inner_joint_explicit();

  // get the length so far
  double len = 2.0 * rmax;

  // compute the movement for the rigid body
  if (dynamic_pointer_cast<PrismaticJointd>(inner))
    mvmt += len + std::fabs(inner->q[0]) + std::fabs(inner->qd[0]);
  else
  {
    assert(dynamic_pointer_cast<SphericalJointd>(inner) ||
           dynamic_pointer_cast<UniversalJointd>(inner) ||
           dynamic_pointer_cast<RevoluteJointd>(inner) ||
           dynamic_pointer_cast<FixedJointd>(inner));
    mvmt += len * inner->qd.norm();
  }

  // keep looping until we arrive at the base link
  while (true)
  {
    // get the inboard link
    rb = inner->get_inboard_link();

    // handle the base specially
    if (rb == base)
    {
      // add in base angular velocity
      mvmt += link_lengths[inner->get_index()] * base_omega.norm();
      break;
    }

    // update the poses 
    JointPtr next_inner = rb->get_inner_joint_explicit();

    // update the movement
    if (dynamic_pointer_cast<PrismaticJointd>(next_inner))
    {
      len += std::fabs(next_inner->q[0]) + std::fabs(next_inner->qd[0]);
      mvmt += std::fabs(next_inner->qd[0]);
    }
    else
    {
      len += link_lengths[inner->get_index()]; 
      assert(dynamic_pointer_cast<SphericalJointd>(next_inner) ||
             dynamic_pointer_cast<UniversalJointd>(next_inner) ||
             dynamic_pointer_cast<RevoluteJointd>(next_inner) ||
             dynamic_pointer_cast<FixedJointd>(next_inner));
      mvmt += len * next_inner->qd.norm();
    }

    // update the joint 
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

  // limit the angular velocity 
  double w0_nrm = w0.norm();
  if (w0_nrm > 2.0*M_PI)
    w0 *= 2.0*M_PI/w0_nrm;

  return std::fabs(n.dot(xd0)) + Vector3d::cross(w0, n).norm()*rmax;
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

  // do not verify that the node name is correct; class may be subclassed
  // assert(strcasecmp(node->name.c_str(), "CCD") == 0);

  // call parent
  CollisionDetection::load_from_xml(node, id_map);
}

/// Implements Base::save_to_xml()
/**
 * \note neither the contact cache nor the pairs currently in collision are
 *       saved
 */
void CCD::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // do not (re)set the node name - there can be derived classes
  // node->name = "CCD";

  // call the parent method 
  CollisionDetection::save_to_xml(node, shared_objects);
}

// Gets an arbitrary point on a feature
Origin3d CCD::get_arbitrary_point(shared_ptr<const Polyhedron::Feature> feat)
{
  if (dynamic_pointer_cast<const Polyhedron::Vertex>(feat)) {
    shared_ptr<const Polyhedron::Vertex> v = boost::static_pointer_cast<const Polyhedron::Vertex>(feat);
    return v->o;
  }
  else if (dynamic_pointer_cast<const Polyhedron::Edge>(feat)) {
    shared_ptr<const Polyhedron::Edge> e = boost::static_pointer_cast<const Polyhedron::Edge>(feat);
    return e->v1->o;
  }
  else {
    shared_ptr<const Polyhedron::Face> fconst = boost::static_pointer_cast<const Polyhedron::Face>(feat);
    shared_ptr<Polyhedron::Face> f = boost::const_pointer_cast<Polyhedron::Face>(fconst);
    Polyhedron::VertexFaceIterator vfi(f, true);
    return (*vfi)->o;
  }
}

void CCD::find_closest_points(boost::shared_ptr<const Polyhedron::Feature> fA, boost::shared_ptr<const Polyhedron::Feature> fB, const Ravelin::Transform3d& wTa, const Ravelin::Transform3d& wTb, Point3d& pA, Point3d& pB){
  if (dynamic_pointer_cast<const Polyhedron::Vertex>(fA)) {
    
    // fA is an vertex
    shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(fA);

    // Setup the point from A.
    Ravelin::Vector3d vAa(vA->o, wTa.source);
    pA = wTa.transform_point(vAa);

    if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {
   
      // fB is a vertex
      shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(fB);

      // Setup the point from B
      Ravelin::Vector3d vBb(vB->o, wTb.source);
    
      //Transforming point into world frame
      pB = wTb.transform_point(vBb);
    }
    else if (dynamic_pointer_cast<const Polyhedron::Edge>(fB)) {

      // fB is an edge
      shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(fB);

      Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), wTb.source);
      Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), wTb.source);
    
      // transform the edge into the world frame
      Ravelin::Vector3d vB1w = wTb.transform_point(vB1b);
      Ravelin::Vector3d vB2w = wTb.transform_point(vB2b);

      // Create line segment
      LineSeg3 line(vB1w,vB2w);
     
      // Compute distance and closest point
      Point3d p;
      double t;
      double dist = CompGeom::calc_dist(line,pA,t,p);
      pB = p;
    }
    else {

      // fB is a face
      boost::shared_ptr<const Pose3d> GLOBAL3D;

      // cast features to non-constant
      boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(fA);
      boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(fB);
      boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);     

      // Transform point from A into B's frame.
      Ravelin::Vector3d vAb = wTb.inverse_transform_point(wTa.transform_point(vAa));
      vAb.pose = GLOBAL3D;                 // hack around plane being in B's frame

      // Find the minimum
      Plane planeB = faceB->get_plane();   // plane will be in B's frame
      double dist = planeB.calc_signed_distance(vAb);

      // project the point onto the plane
      Ravelin::Vector3d vAb_on_planeB = vAb - planeB.get_normal()*dist;
      vAb_on_planeB.pose = wTb.source;
      // this is correct because when v-clip ends and a face vertex case
      // the vertex will always be in the voronoi region, and therefore,
      // the vertex projection is always on face B.
      vAb_on_planeB.pose = wTb.source;
      pB = wTb.transform_point(vAb_on_planeB);
    }
  }
  else if (dynamic_pointer_cast<const Polyhedron::Edge>(fA)) {

    // fA is an edge
    if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {

      // already implemented, just need to flip it.
      find_closest_points(fB, fA, wTb, wTa, pB, pA);
    }
    else if (dynamic_pointer_cast<const Polyhedron::Edge>(fB)) {
      // Features are two edges.

      // Cast pointers
      boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(fA);
      boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(fB);
    
      //create vectors
      Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), wTa.source);
      Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), wTa.source);

      Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), wTb.source);
      Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), wTb.source);
    
      //transform features to global frame
      Ravelin::Vector3d vA1w = wTa.transform_point(vA1a);
      Ravelin::Vector3d vA2w = wTa.transform_point(vA2a);
      Ravelin::Vector3d vB1w = wTb.transform_point(vB1b);
      Ravelin::Vector3d vB2w = wTb.transform_point(vB2b);

      //create line segment
      LineSeg3 lineA(vA1w,vA2w);
      LineSeg3 lineB(vB1w,vB2w);

      // compute distance and closest point
      Point3d p1;
      Point3d p2;
      double dist = CompGeom::calc_closest_points(lineA,lineB,p1,p2);
      pA = p1;
      pB = p2;
    }
    else {

      // should not be reached b/c V-Clip does not return face/edge.
      assert(false);
    }
  }
  else {
     if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {

      // already implemented, just need to flip it.
      find_closest_points(fB, fA, wTb, wTa, pB, pA);
    }
    else{
      // should not be reached
      assert(false);
    }

  }
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
      _rmax[i] = i->get_farthest_point_distance() + i->get_geometry()->get_pose()->x.norm();

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


