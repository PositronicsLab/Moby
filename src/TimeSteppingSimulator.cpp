/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <unistd.h>
#include <boost/tuple/tuple.hpp>
#include <Moby/XMLTree.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Dissipation.h>
#include <Moby/ControlledBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/SustainedUnilateralConstraintSolveFailException.h>
#include <Moby/InvalidStateException.h>
#include <Moby/InvalidVelocityException.h>
#include <Moby/TimeSteppingSimulator.h>

#ifdef USE_OSG
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/PositionAttitudeTransform>
#include <osg/Quat>
#endif // USE_OSG

using std::endl;
using std::set;
using std::list;
using std::vector;
using std::map;
using std::make_pair;
using std::multimap;
using std::pair;
using boost::tuple;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Default constructor
TimeSteppingSimulator::TimeSteppingSimulator()
{
  min_step_size = 1e-7;
}

/// Steps the simulator forward by the given step size
double TimeSteppingSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();

  // determine the set of collision geometries
  determine_geometries();

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif

  // clear stored derivatives
  _current_dx.resize(0);

  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << " by " << step_size << std::endl;
  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q, qd;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(cb);
      db->get_generalized_coordinates_euler(q);
      db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " Euler coordinates (before): " << q << std::endl;
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " spatial velocity (before): " << qd << std::endl;
    }
  }

  // do broad phase collision detection (must be done before any Euler steps)
  broad_phase(step_size);

  // compute pairwise distances at the current configuration
  calc_pairwise_distances();

  // do the Euler step
  step_si_Euler(step_size);

  // call the callback
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  // do constraint stabilization
  shared_ptr<ConstraintSimulator> simulator = dynamic_pointer_cast<ConstraintSimulator>(shared_from_this());
  FILE_LOG(LOG_SIMULATOR) << "stabilization started" << std::endl;
  cstab.stabilize(simulator);
  FILE_LOG(LOG_SIMULATOR) << "stabilization done" << std::endl;

  // write out constraint violation
  #ifndef NDEBUG
  std::ofstream cvio("cvio.dat", std::ostream::app);
  double vio = 0.0;
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
    vio = std::min(vio, _pairwise_distances[i].dist);
  cvio << vio << std::endl;
  cvio.close();
  #endif

  return step_size;
}

void TimeSteppingSimulator::update_visualization()
{
  // Update transforms first.
  Simulator::update_visualization();

  // Only do the following part if built with OSG support.
  #ifdef USE_OSG

  // Get the set of all geometries.
  vector<CollisionGeometryPtr> cgs;
  const auto& bodies = get_dynamic_bodies();
  for (auto body : bodies)
  {
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(body);
    if (rb)
      cgs.insert(cgs.end(), rb->geometries.begin(), rb->geometries.end());
    else
    {
      ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);
      for (unsigned i = 0; i < ab->get_links().size(); ++i)
      {
        RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(ab->get_links()[i]);
        cgs.insert(cgs.end(), rb->geometries.begin(), rb->geometries.end());
      }
    }
  }

  // Process all of the polyhedral geometries with non-zero compliant layer
  // depths.
  for (auto geom : cgs)
  {
    // Fast check: skip bodies without a compliant layer.
    if (geom->compliant_layer_depth == 0)
      continue;

    // Try to cast to a polyhedral body.
    shared_ptr<PolyhedralPrimitive> pp = dynamic_pointer_cast<PolyhedralPrimitive>(geom->get_geometry());
    if (!pp)
      continue;

    // Get the polyhedron.
    const Polyhedron& polyhedron = pp->get_polyhedron();

    // Process each vertex of the polyhedron.
    for (unsigned i = 0; i< polyhedron.get_vertices().size(); ++i)
    {
      shared_ptr<Polyhedron::Vertex> v = polyhedron.get_vertices()[i];
      Origin3d extrusion_direction = v->sum_coincident_normals();
      Origin3d new_vert = v->o + extrusion_direction * geom->compliant_layer_depth;
      osg::Vec3& osg_vert = *pp->get_visualization_vertex(v);
      osg_vert = osg::Vec3(new_vert.x(), new_vert.y(), new_vert.z());
    }
    pp->dirty_vertex_visualization_data();
  }

  // Process each pairwise distance information structure.
  for (const auto& pdi : _pairwise_distances)
  {
    // TODO(edrumwri): Do this properly.
    // Determine the contact plane.
    Plane contact_plane(0, 1, 0, 0);

    // Get the first geometry, see whether it is polyhedral and has a compliant
    // layer. If so, adjust it.
    if (pdi.a->compliant_layer_depth > 0)
    {
      shared_ptr<PolyhedralPrimitive> pp = dynamic_pointer_cast<PolyhedralPrimitive>(pdi.a->get_geometry());
      if (pp)
        adjust_compliant_polyhedral_geometry(pdi.a, contact_plane);
    }

    if (pdi.b->compliant_layer_depth > 0)
    {
      shared_ptr<PolyhedralPrimitive> pp = dynamic_pointer_cast<PolyhedralPrimitive>(pdi.b->get_geometry());
      if (pp)
        adjust_compliant_polyhedral_geometry(pdi.b, -contact_plane);
    }
  }

  #endif // USE_OSG

  /*
  // Get the set of contacting bodies.
  for (const auto& pdi : _pairwise_distances)
  {
    // Check whether the bodies are intersecting.
    if (pdi.dist < 0)
    {

    }
  }

  // TODO: Cache closest features to speed v-clip.

  // For all polyhedral bodies with a compliant layer, extrude the faces by
  // the true depth - ε. Intersecting two extruded such polyhedral bodies, or
  // one polyhedral body with a compliant layer against another body without
  // a compliant layer, should yield no intersection at that configuration.
  // Intersecting a polyhedral body with a non-polyhedral body- both of which
  // use a compliant layer- means that we can't update the visualization based
  // on their interaction.

  // For all polyhedral bodies with a compliant layer, extrude the faces by
  // the true depth.
  for (const auto body : polyhedral_bodies)
  {

  }
   */
}

// Adjusts a polyhedral geometry with a compliant layer with respect to the
// given contact plane.
void TimeSteppingSimulator::adjust_compliant_polyhedral_geometry(CollisionGeometryPtr cg, const Plane& contact_plane)
{
  // Get the polyhedral geometry.
  shared_ptr<PolyhedralPrimitive> pp = static_pointer_cast<PolyhedralPrimitive>(cg->get_geometry());

  // Get all vertices below the contact plane.
  const vector<shared_ptr<Polyhedron::Vertex>>& vertices = pp->get_polyhedron().get_vertices();
  Transform3d wTpp = Pose3d::calc_relative_pose(pp->get_pose(cg), GLOBAL);
  vector<shared_ptr<Polyhedron::Vertex>> below_verts;
  for (auto v : vertices)
  {
    Point3d v_pp(v->o + v->sum_coincident_normals() * cg->compliant_layer_depth,
                 pp->get_pose(cg));
    Point3d v_w = wTpp.transform_point(v_pp);
    double signed_dist = contact_plane.calc_signed_distance(v_w);
    if (signed_dist < -NEAR_ZERO)
      below_verts.push_back(v);
  }

  // Get the set of normals (faces) common to all of these vertices.
  const vector<shared_ptr<Polyhedron::Face>>& all_faces = pp->get_polyhedron().get_faces();
  set<shared_ptr<Polyhedron::Face>> common_faces(all_faces.begin(), all_faces.end());
  for (auto v : below_verts)
  {
    // Get all faces from this vertex.
    set<shared_ptr<Polyhedron::Face>> v_faces;
    for (auto weak_edge : v->e)
    {
      shared_ptr<Polyhedron::Edge> edge(weak_edge);
      v_faces.insert(edge->face1);
      v_faces.insert(edge->face2);
    }

    // Do the set intersection.
    vector<shared_ptr<Polyhedron::Face>> output;
    std::set_intersection(common_faces.begin(), common_faces.end(),
                          v_faces.begin(), v_faces.end(), std::back_inserter(output));
    common_faces.clear();
    common_faces.insert(output.begin(), output.end());
  }

  // See which normal requires the least compression to move all points above
  // the contact plane. Assume that we want to use a normal of direction n to
  // push a point below the plane *onto* the contact plane. In other words,
  // point p currently lies at distance q below the contact plane <y, d>:
  // y'p - d = q
  // Direction n from some face would push p further below the contact plane.
  // Instead, projecting p along -n should allow us to to put the point right
  // on the contact plane:
  // αn'p = -q
  // Putting α on the left hand side:
  // α = (-y'p + d)/(-n'p)
  // α's interpretation is how far a piston has to be compressed to travel one
  // unit of distance along the separation vector. If the piston is parallel to
  // the separation vector, α will reach its optimal (minimal) value of 1.
  // As the vectors become less parallel- including up to the case of
  // orthogonality, for which α will become infinity- α will become larger.
  double min_max_alpha = std::numeric_limits<double>::infinity();
  shared_ptr<Polyhedron::Face> best_face;
  for (auto f : common_faces)
  {
    double max_alpha = 1.0;
    for (auto v : below_verts)
    {
      // Get the vertex in the global frame.
      const Point3d v_pp(v->o + v->sum_coincident_normals() * cg->compliant_layer_depth,
                         pp->get_pose(cg));
      const Point3d p = wTpp.transform_point(v_pp);

      // Get the signed distance of the point below the plane.
      const double q = contact_plane.calc_signed_distance(p);

      // TODO(edrumwri): Why is p not a vector for the operation below?
      // Compute alpha.
      const double alpha = -q / f->get_plane().get_normal().dot(p);
      assert(alpha >= 1.0 - std::numeric_limits<double>::epsilon());

      // Update the maximum alpha.
      max_alpha = std::max(alpha, max_alpha);
    }

    // Update the minimum alpha.
    if (max_alpha < min_max_alpha)
    {
      min_max_alpha = max_alpha;
      best_face = f;
    }
  }

  // Get the piston direction.
  const Vector3d n = best_face->get_plane().get_normal();

  // Compress the vertices.
  for (auto v : below_verts)
  {
    // Get the vertex in the global frame.
    const Point3d v_pp(v->o + v->sum_coincident_normals() * cg->compliant_layer_depth,
                       pp->get_pose(cg));
    Point3d p = wTpp.transform_point(v_pp);

    // Get the signed distance of the point below the plane.
    const double q = contact_plane.calc_signed_distance(p);

    // Recompute alpha.
    const double alpha = -q / n.dot(p);

    // Compute the new vertex position (in the global frame).
    p += n * alpha * q;

    // Transform the vertex location back to the primitive pose.
    Point3d vnew_pp = wTpp.inverse_transform_point(p);

    // Update the visualization vertex.
    osg::Vec3& vv = *pp->get_visualization_vertex(v);
    vv = osg::Vec3(vnew_pp.x(), vnew_pp.y(), vnew_pp.z());
  }
}

/*
void TimeSteppingSimulator::update_polyhedral_visualization(CollisionGeometryPtr cg)
{
  // Verify that this geometry is polyhedral.
  PrimitivePtr geometry = cg->get_geometry();
  auto poly_prim = dynamic_pointer_cast<PolyhedralPrimitive>(geometry);
  assert(poly_prim);

  // If the geometry does not have a compliant layer, quit now.
  if (cg->compliant_layer_depth == 0.0)
    return;

  // Test this geometry against every other.
  for (const auto& pdi : _pairwise_distances)
  {
    // See whether one of the geometries is cg.
    if (pdi.a != cg && pdi.b != cg)
      continue;

    // Get the other geometry.
    CollisionGeometryPtr cg_other = pdi.a;
    if (cg_other == cg)
      cg_other = pdi.b;

    // If the geometries are not intersecting, continue.
    if (pdi.dist > 0)
      continue;

    // Get the contact plane.

    //

    // If the geometries are intersecting, the extension quantity is (a) the
    // signed distance between the rigid cores minus the compliant layer (if
    // only this body has a compliant layer) or (b) half the signed distance
    // between the rigid cores minus the compliant layer depths (if both
    // geometries have a compliant layer).
    if (pdi.dist < 0)
    {
      // TODO: What if one of the geometries can't extend to the full extension
      // length?
      const double core_dist = pdi.dist + cg->compliant_layer_depth + cg_other->compliant_layer_depth;
      const double extension = (cg_other->compliant_layer_depth > 0) ? core_dist/2 : core_dist;
      min_extension = std::min(extension, min_extension);
    }
  }

  // Compress the pistons slightly.
  min_extension -= std::sqrt(std::numeric_limits<double>::epsilon());
  assert(min_extension >= 0);

  // Verify that cg does not intersect with any other (geometry? polyhedral
  // geometry?) at this extension.

  // Loop over all faces of cg's primitive, extruding each face to its full
  // extent. If the extruded face 
}
*/

/// Does a full integration cycle (but not necessarily a full step)
double TimeSteppingSimulator::do_mini_step(double dt)
{
  VectorNd q, qd, qdd;
  std::vector<VectorNd> qsave;

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::do_mini_step() - current time is " << current_time << std::endl;

  // init qsave to proper size and save generalized coordinates for all bodies
  qsave.resize(_bodies.size());
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->get_generalized_coordinates_euler(qsave[i]);
  }

  // ****
  // prepare to compute the first conservative step

  // do broad phase collision detection
  broad_phase(dt);

  // compute pairwise distances
  calc_pairwise_distances();

  // compute the first conservative step, unless the minimum step is equal
  // to dt.
  double h = dt;
  if (dt == min_step_size)
  {
    h = calc_next_CA_Euler_step();
    if (h < std::min(dt, min_step_size))
      h = std::min(dt, min_step_size);
    else if (h > dt)
      h = dt;
  }

  // integrate the bodies' positions by h
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->set_generalized_coordinates_euler(qsave[i]);
    db->get_generalized_velocity(DynamicBodyd::eEuler, q);
    q *= h;
    q += qsave[i];
    db->set_generalized_coordinates_euler(q);
  }

  FILE_LOG(LOG_SIMULATOR) << "Position integration ended w/h = " << h << std::endl;

  // only integrate velocity if h > 0
  if (h > 0.0)
  { 
    // prepare to calculate forward dynamics
    precalc_fwd_dyn();

    // compute forward dynamics
    calc_fwd_dyn(h);

    // integrate the bodies' velocities forward by h
    for (unsigned i=0; i< _bodies.size(); i++)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
      db->get_generalized_acceleration(qdd);
      qdd *= h;
      db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
      FILE_LOG(LOG_DYNAMICS) << "old velocity: " << qd << std::endl; 
      qd += qdd;
      db->set_generalized_velocity(DynamicBodyd::eSpatial, qd);
      FILE_LOG(LOG_DYNAMICS) << "new velocity: " << qd << std::endl; 
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "Integrated velocity by " << h << std::endl;

  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints();

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(h);

  // update the time
  current_time += h;

  // do a mini-step callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn((ConstraintSimulator*) this);

  return h;
}

/// Checks to see whether all constraints are met
bool TimeSteppingSimulator::constraints_met(const std::vector<PairwiseDistInfo>& current_pairwise_distances)
{
  // see whether constraints are met to specified tolerance
  for (unsigned i=0; i< _pairwise_distances.size(); i++)
  {
    const PairwiseDistInfo& pdi = _pairwise_distances[i];
    if (current_pairwise_distances[i].dist < 0.0 &&
        pdi.dist < current_pairwise_distances[i].dist - NEAR_ZERO)
    {
      // check whether one of the bodies is compliant
      RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
      RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
      if (rba->compliance == RigidBody::eCompliant || 
          rbb->compliance == RigidBody::eCompliant)
        continue;

      FILE_LOG(LOG_SIMULATOR) << "signed distance between " << rba->id << " and " << rbb->id << "(" << pdi.dist << ") below tolerance: " << (pdi.dist - current_pairwise_distances[i].dist) << std::endl;

      return false;
    }
  }

  return true;
}

/// Gets the current set of contact geometries
std::set<sorted_pair<CollisionGeometryPtr> > TimeSteppingSimulator::get_current_contact_geoms() const
{
  std::set<sorted_pair<CollisionGeometryPtr> > contact_geoms;

  for (unsigned i=0; i< _rigid_constraints.size(); i++)
    contact_geoms.insert(make_sorted_pair(_rigid_constraints[i].contact_geom1, _rigid_constraints[i].contact_geom2));

  return contact_geoms; 
}

/// Finds the next event time assuming constant velocity
/**
 * This method returns the next possible time of contact, discarding current
 * contacts from consideration. 
 * \note proper operation of this function is critical. If the function
 *       improperly designates an event as not occuring at the current time,
 *       calc_next_CA_Euler_step(.) will return a small value and prevent large
 *       integration steps from being taken. If the function improperly
 *       designates an event as occuring at the current time, constraint
 *       violation could occur.
 */
double TimeSteppingSimulator::calc_next_CA_Euler_step() const
{
  const double INF = std::numeric_limits<double>::max();
  double next_event_time = INF;

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::calc_next_CA_Euler_step entered" << std::endl; 

  // process each articulated body, looking for next joint events
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    // see whether the i'th body is articulated
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(_bodies[i]);
    if (!ab)
      continue;

    // get limit events in [t, t+dt] (if any)
    const vector<shared_ptr<Jointd> >& joints = ab->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
    {
      JointPtr joint = dynamic_pointer_cast<Joint>(joints[i]);
      for (unsigned j=0; j< joint->num_dof(); j++)
      {
        if (joint->q[j] < joint->hilimit[j] && joint->qd[j] > 0.0)
        {
          double t = (joint->hilimit[j] - joint->q[j])/joint->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
        if (joint->q[j] > joint->lolimit[j] && joint->qd[j] < 0.0)
        {
          double t = (joint->lolimit[j] - joint->q[j])/joint->qd[j];
          next_event_time = std::min(next_event_time, t);
        }
      }
    }
  }

  // if the distance between any pair of bodies is sufficiently small
  // get next possible event time
  BOOST_FOREACH(PairwiseDistInfo pdi, _pairwise_distances)
  {
    // only process if neither of the bodies is compliant
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(pdi.a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(pdi.b->get_single_body());
    if (rba->compliance == RigidBody::eCompliant || 
        rbb->compliance == RigidBody::eCompliant)
      continue; 

      // compute an upper bound on the event time
      double event_time = _coldet->calc_CA_Euler_step(pdi, contact_dist_thresh);

      FILE_LOG(LOG_SIMULATOR) << "Next contact time between " << pdi.a->get_single_body()->body_id << " and " << pdi.b->get_single_body()->body_id << ": " << event_time << std::endl;

      // not a current event, find when it could become active
      next_event_time = std::min(next_event_time, event_time);
    }

  FILE_LOG(LOG_SIMULATOR) << "TimeSteppingSimulator::calc_next_CA_Euler_step exited" << std::endl; 

  return next_event_time;
}

/// Does a semi-implicit step
/*
void TimeSteppingSimulator::step_si_Euler(double dt)
{
  FILE_LOG(LOG_SIMULATOR) << "-- doing semi-implicit Euler step" << std::endl;
  const double INF = std::numeric_limits<double>::max();

  VectorNd q, qd, qdd;
  std::vector<VectorNd> qsave;

  // init qsave to proper size
  qsave.resize(_bodies.size());

  // save generalized coordinates for all bodies
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->get_generalized_coordinates_euler(qsave[i]);
  }

  // do broad phase collision detection
  broad_phase(dt);

  // compute pairwise distances
  calc_pairwise_distances();

  // integrate the bodies' positions by dt 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->set_generalized_coordinates_euler(qsave[i]);
    db->get_generalized_velocity(DynamicBodyd::eEuler, q);
    q *= dt;
    q += qsave[i];
    db->set_generalized_coordinates_euler(q);
  }

  // prepare to calculate forward dynamics
  precalc_fwd_dyn();

  // apply compliant unilateral constraint forces
  calc_compliant_unilateral_constraint_forces();

  // compute forward dynamics
  calc_fwd_dyn(dt);

  // integrate the bodies' velocities forward by dt 
  for (unsigned i=0; i< _bodies.size(); i++)
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(_bodies[i]);
    db->get_generalized_acceleration(qdd);
    qdd *= dt;
    db->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
    FILE_LOG(LOG_DYNAMICS) << "old velocity: " << qd << std::endl; 
    qd += qdd;
    db->set_generalized_velocity(DynamicBodyd::eSpatial, qd);
    FILE_LOG(LOG_DYNAMICS) << "new velocity: " << qd << std::endl; 
  }

  // dissipate some energy
  if (_dissipator)
  {
    vector<shared_ptr<DynamicBodyd> > bodies;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
      bodies.push_back(dynamic_pointer_cast<DynamicBodyd>(cb));
    _dissipator->apply(bodies);
  }

  // recompute pairwise distances
  calc_pairwise_distances();

  // find unilateral constraints
  find_unilateral_constraints(contact_dist_thresh);

  // handle any impacts
  calc_impacting_unilateral_constraint_forces(dt);

  // update the time
  current_time += dt;

  // do a mini-step callback
  if (post_mini_step_callback_fn)
    post_mini_step_callback_fn((ConstraintSimulator*) this);

  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(cb);
      db->get_generalized_coordinates_euler(q);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " position (after integration): " << q << std::endl;
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "-- semi-implicit Euler step completed" << std::endl;
}
*/

/// Does a semi-implicit step (version with conservative advancement)
void TimeSteppingSimulator::step_si_Euler(double dt)
{
  FILE_LOG(LOG_SIMULATOR) << "-- doing semi-implicit Euler step" << std::endl;
  const double INF = std::numeric_limits<double>::max();

  // do a number of mini-steps until integrated forward fully
  double h = 0.0;
  while (h < dt)
    h += do_mini_step(dt-h);

  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q;
    BOOST_FOREACH(ControlledBodyPtr cb, _bodies)
    {
      shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(cb);
      db->get_generalized_coordinates_euler(q);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->body_id << " position (after integration): " << q << std::endl;
    }
  }

  FILE_LOG(LOG_SIMULATOR) << "-- semi-implicit Euler step completed" << std::endl;
}

/// Implements Base::load_from_xml()
void TimeSteppingSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // do not verify node name b/c this may be a parent class
  // assert(strcasecmp(node->name.c_str(), "TimeSteppingSimulator") == 0);

  // first, load all data specified to the ConstraintSimulator object
  ConstraintSimulator::load_from_xml(node, id_map);

  // read the minimum step size
  XMLAttrib* min_step_attrib = node->get_attrib("min-step-size");
  if (min_step_attrib)
    min_step_size = min_step_attrib->get_real_value();
}

/// Implements Base::save_to_xml()
void TimeSteppingSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call ConstraintSimulator's save method 
  ConstraintSimulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "TimeSteppingSimulator";

  // save the minimum step size
  node->attribs.insert(XMLAttrib("min-step-size", min_step_size));
}


