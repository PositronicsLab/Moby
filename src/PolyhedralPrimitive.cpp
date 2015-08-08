/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/SpherePrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/PolyhedralPrimitive.h>

using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

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
  if (is_convex() && polyp && p->is_convex())
  {
    shared_ptr<const PolyhedralPrimitive> bthis = dynamic_pointer_cast<const PolyhedralPrimitive>(shared_from_this());
    shared_ptr<const Pose3d> poseA = pthis.pose;
    shared_ptr<const Pose3d> poseB = pp.pose;
    shared_ptr<const Polyhedron::Feature> closestA, closestB;

    // try v-clip
    double dist = Polyhedron::vclip(bthis, polyp, poseA, poseB, closestA, closestB); 

    // if dist is negative, get precise distance and closest points using
    // calc_minkowski(.)
    if (dist < 0.0)
    {
      // compute the Minkowski difference
      Polyhedron mdiff = Polyhedron::calc_minkowski_diff(bthis, polyp, poseA, poseB);

      // TODO: Bjoern, this is where you find the closest feature to the origin,
      // return the negation of its distance to the origin, and determine
      // closest points on the two polyhedral primitives

      // 3. find closest feature to the origin
      shared_ptr<const Pose3d> GLOBAL3D;
      Ravelin::Origin3d o(0,0,0);
      Ravelin::Vector3d origin_vector(o,GLOBAL3D);
      
      std::vector<boost::shared_ptr<Polyhedron::Vertex> > vertices = mdiff.get_vertices();
      double min_dist = std::numeric_limits<double>::max();
      boost::shared_ptr<std::pair<int, int> > min_pair;
      for(int i = 0; i < vertices.size(); i++)
      {
        Ravelin::Vector3d vertex_vector (vertices[i]->o,GLOBAL3D);
        double dist = (origin_vector - vertex_vector).norm();
        if(dist < min_dist)
        {
          min_dist = dist;
          min_pair= boost::static_pointer_cast<std::pair<int, int> >(vertices[i]->data);
        }
      }

      std::vector<boost::shared_ptr<Polyhedron::Vertex> > vvthis = bthis->get_polyhedron().get_vertices();
      std::vector<boost::shared_ptr<Polyhedron::Vertex> > vvp = bthis->get_polyhedron().get_vertices();
      boost::shared_ptr<Polyhedron::Vertex> vthis = vvthis[min_pair->first];
      boost::shared_ptr<Polyhedron::Vertex> vp = vvp[min_pair->second];
      pthis = vthis->o;
      pp = vp->o;
      return -min_dist;
  }
  else
    return dist;
 } 
}
