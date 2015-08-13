/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <limits>
#include <stack>
#include <list>
#include <set>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/CollisionDetection.h>

using namespace Ravelin;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::vector;
using std::make_pair;
using std::list;
using std::pair;
using namespace Moby;

/// Default broad phase function (checks everything)
void CollisionDetection::broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
{
  // clear vector of bodies to be checked
  to_check.clear();

  // get the set of rigid bodies
  vector<RigidBodyPtr> rbs;
  for (unsigned i=0; i< bodies.size(); i++)
  {
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(bodies[i]);
    if (ab)
    {
      BOOST_FOREACH(shared_ptr<RigidBodyd> link, ab->get_links()) 
        rbs.push_back(dynamic_pointer_cast<RigidBody>(link));
    }
    else
      rbs.push_back(dynamic_pointer_cast<RigidBody>(bodies[i]));
  }

  // process each pair of bodies
  for (unsigned i=0; i< rbs.size(); i++)
    BOOST_FOREACH(CollisionGeometryPtr gi, rbs[i]->geometries)
      for (unsigned j=i+1; j< rbs.size(); j++)
        if (rbs[i]->is_enabled() || rbs[j]->is_enabled())
          BOOST_FOREACH(CollisionGeometryPtr gj, rbs[j]->geometries)
            to_check.push_back(std::make_pair(gi, gj));
}

/// Creates a contact constraint given the bare-minimum info
UnilateralConstraint CollisionDetection::create_contact(CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Vector3d& normal, double violation)
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
    shared_ptr<SingleBodyd> sb1 = e.contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = e.contact_geom2->get_single_body();
    if (sb2->body_id < sb1->body_id)
    {
      std::swap(e.contact_geom1, e.contact_geom2);
      e.contact_normal = -e.contact_normal;
    }
  }

  // transform contact point and normal to global frame
  e.contact_point = Pose3d::transform_point(GLOBAL, e.contact_point);
  e.contact_normal = Pose3d::transform_vector(GLOBAL, e.contact_normal);

  // compute contact tangents
  e.determine_contact_tangents();

  // compute normal and tangent time derivatives
  e.compute_contact_dots();

  return e;
}

