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
    }
  }
}

