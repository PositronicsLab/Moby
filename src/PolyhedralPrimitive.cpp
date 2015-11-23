/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/SpherePrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/PolyhedralPrimitive.h>

using std::vector;
using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Computes the signed distance between two polyhedra
double PolyhedralPrimitive::calc_signed_dist(shared_ptr<const PolyhedralPrimitive> p, Point3d& pthis, Point3d& pp) const
{
  const double INF = std::numeric_limits<double>::max();
  shared_ptr<TessellatedPolyhedron> tpoly;

  // verify that p is convex
  if (!is_convex() && !p->is_convex())
    throw std::runtime_error("Polyhedron is not convex!");

  // if the primitive is polyhedral and convex, can use vclip 
  shared_ptr<const PolyhedralPrimitive> bthis = dynamic_pointer_cast<const PolyhedralPrimitive>(shared_from_this());
  shared_ptr<const Pose3d> poseA = pthis.pose;
  shared_ptr<const Pose3d> poseB = pp.pose;
  shared_ptr<const Polyhedron::Feature> closestA, closestB;

  // attempt to use vclip
  double dist = Polyhedron::vclip(bthis, p, poseA, poseB, closestA, closestB); 
  if (dist >= 0.0)
    return dist;

  // compute transforms
  Transform3d wTa = Pose3d::calc_relative_pose(poseA, GLOBAL);
  Transform3d wTb = Pose3d::calc_relative_pose(poseB, GLOBAL);

  // compute half-spaces
  vector<std::pair<Vector3d, double> > hs;
  PolyhedralPrimitive::get_halfspaces(bthis->get_polyhedron(), poseA, wTa, std::back_inserter(hs));
  PolyhedralPrimitive::get_halfspaces(p->get_polyhedron(), poseB, wTb, std::back_inserter(hs));

  // find a point interior to both
  Ravelin::Origin3d ip;
  double dist2 = -CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);
  assert(dist2 < NEAR_ZERO);

  // attempt to calculate the half-space intersection
  try
  {
    tpoly = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
  }
  catch (NumericalException e)
  {
    // if we're here, then the volume of intersection is an area, find
    // closest points and quit
    Vector3d ip0(ip, GLOBAL);
    pthis = wTa.inverse_transform_point(ip0);
    pp = wTb.inverse_transform_point(ip0);

    return 0.0;
  } 

 // setup minimum distance and the contact plane offset
  double min_dist = INF;
  double offset = INF;

  // get the vertices from the polyhedron
  const std::vector<Ravelin::Origin3d>& vertices = tpoly->get_vertices();

  // determine the signed distance 
  const std::vector<IndexedTri>& facets = tpoly->get_facets();
  for (unsigned i=0; i< facets.size(); i++)
  {
    // setup a triangle
    Triangle tri(Point3d(vertices[facets[i].a], GLOBAL), 
                 Point3d(vertices[facets[i].b], GLOBAL),
                 Point3d(vertices[facets[i].c], GLOBAL));

    // TODO: ensure that the plane containing this triangle is from 
    // polyhedron B by checking that signed distances from all vertices of B 
    // are negative
    // NOTE: this is currently an O(n) operation, but it could be turned into
    //       an O(lg N) one

    // get the reverse of the normal
    Ravelin::Vector3d ncand = -tri.calc_normal();

    // TODO: get the extremal point in the direction of the inverse normal
    // NOTE: this is currently an O(n) operation, but it could be turned into
    //       an O(lg N) one
    Point3d p;
    assert(false);

    // compute the distance of the extremal point from the face
    double dist = ncand.dot(p);
    assert(dist > -NEAR_ZERO);
    if (dist < min_dist)
    {
      // if the absolute distance is less than the minimum, we've found our
      // normal
      min_dist = dist;
      offset = -ncand.dot(Point3d(vertices[facets[i].a], GLOBAL));
    }
  }

  // return the negative minimum distance; don't worry about closest points
  // during interpenetration
  return -min_dist;
}

/// Computes the signed distance between this primitive and another
double PolyhedralPrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // now try polyhedron/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
  {
    throw std::runtime_error("Sphere/polyhedron distance not yet implemented");
    return 0.0;
  } 

  // now try polyhedron/plane
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return planep->calc_signed_dist(bthis, pp, pthis);
  }

  // now try polyhedron/heightmap
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
  {
    throw std::runtime_error("Heightmap/polyhedron distance not yet implemented");
    return 0.0;
  }

  // now try convex polyhedron/convex polyhedron
  shared_ptr<const PolyhedralPrimitive> polyp = dynamic_pointer_cast<const PolyhedralPrimitive>(p);
  if (polyp)
    return calc_signed_dist(polyp, pthis, pp);

  throw std::runtime_error("Unanticipated signed distance types!");
  return 0.0;
}
