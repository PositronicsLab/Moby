#include <list>
#include <map>
#include <algorithm>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/Joint.h>
#include <Moby/SignedDistDot.h>

using boost::dynamic_pointer_cast;
using std::vector;
using std::map;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Computes the signed distance between two single bodies
double SignedDistDot::calc_signed_dist(shared_ptr<SingleBodyd> sb1, shared_ptr<SingleBodyd> sb2)
{
  Point3d dummy1, dummy2;

  // set the signed distance to be the minimum
  double min_dist = std::numeric_limits<double>::max();

  // get the two rigid bodies
  shared_ptr<RigidBody> rb1 = dynamic_pointer_cast<RigidBody>(sb1);
  shared_ptr<RigidBody> rb2 = dynamic_pointer_cast<RigidBody>(sb2);

  BOOST_FOREACH(CollisionGeometryPtr cg1, rb1->geometries)
    BOOST_FOREACH(CollisionGeometryPtr cg2, rb2->geometries)
      min_dist = std::min(min_dist, CollisionGeometry::calc_signed_dist(cg1, cg2, dummy1, dummy2));

  return min_dist;
}

/// Computes the Jacobian for the time derivative of the signed distance functions vs. impulses applied at contact points
/**
 * Assume the signed distance function is Phi(q(t)), so 
 * Phi(q(t_0 + \Delta t))_{t_0} \approx Phi(q(t_0)) + d/dt Phi(q(t_0)) * dt ==>
 * d/dt Phi(q(t_0)) \approx (Phi(q(t_0 + \Delta t)) - Phi(q(t_0))/dt 
 */
void SignedDistDot::compute_signed_dist_dot_Jacobians(UnilateralConstraintProblemData& q, MatrixNd& Cdot_iM_CnT, MatrixNd& Cdot_iM_CsT, MatrixNd& Cdot_iM_CtT, MatrixNd& Cdot_iM_LT, VectorNd& Cdot_v)
{
  const double DT = NEAR_ZERO;
  vector<shared_ptr<DynamicBodyd> > tmp_supers1, tmp_supers2, isect;
  VectorNd gc, gv;

  // get all pairs of bodies involved in contact
  vector<shared_ptr<DynamicBodyd> > ubodies;
  for (unsigned i=0; i< q.signed_distances.size(); i++)
  {
    // get the two single bodies
    shared_ptr<SingleBodyd> s1 = q.signed_distances[i].a->get_single_body();
    shared_ptr<SingleBodyd> s2 = q.signed_distances[i].b->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> sb1 = ImpactConstraintHandler::get_super_body(s1); 
    shared_ptr<DynamicBodyd> sb2 = ImpactConstraintHandler::get_super_body(s2); 

    // add the bodies to ubodies
    ubodies.push_back(sb1);
    ubodies.push_back(sb2);
  } 

  // get all unique bodies involved in contact
  std::sort(ubodies.begin(), ubodies.end());
  ubodies.erase(std::unique(ubodies.begin(), ubodies.end()), ubodies.end());

  // save all configurations for all bodies involved in contact
  map<shared_ptr<DynamicBodyd>, VectorNd> gc_map;
  for (unsigned i=0; i< ubodies.size(); i++)
    ubodies[i]->get_generalized_coordinates_euler(gc_map[ubodies[i]]);

  // save all velocities for all bodies involved in contact
  map<shared_ptr<DynamicBodyd>, VectorNd> gv_map;
  for (unsigned i=0; i< ubodies.size(); i++)
    ubodies[i]->get_generalized_velocity(DynamicBodyd::eSpatial, gv_map[ubodies[i]]);

  // resize Cdot(v)
  Cdot_v.resize(q.signed_distances.size());

  // for each pair of bodies
  for (unsigned k=0; k< q.signed_distances.size(); k++)
  {
    // get the two single bodies
    shared_ptr<SingleBodyd> s1 = q.signed_distances[k].a->get_single_body();
    shared_ptr<SingleBodyd> s2 = q.signed_distances[k].b->get_single_body();

    // get the signed distance between the two bodies
    double phi = q.signed_distances[k].dist;

    // integrates bodies' positions forward
    shared_ptr<DynamicBodyd> sup1 = ImpactConstraintHandler::get_super_body(s1);
    shared_ptr<DynamicBodyd> sup2 = ImpactConstraintHandler::get_super_body(s2);
    tmp_supers1.clear();
    tmp_supers1.push_back(sup1);
    if (sup1 != sup2)
      tmp_supers1.push_back(sup2);
    integrate_positions(tmp_supers1, DT);

    // compute the signed distance function
    double phi_new = calc_signed_dist(s1, s2);

    // set the appropriate entry of Cdot(v) 
    Cdot_v[k] = (phi_new - phi)/DT;  

    // restore coordinates and velocities
    restore_coords_and_velocities(tmp_supers1, gc_map, gv_map);
  }  

  // resize the Jacobians
  Cdot_iM_CnT.resize(q.signed_distances.size(), q.N_CONTACTS);
  Cdot_iM_CsT.resize(q.signed_distances.size(), q.N_CONTACTS);
  Cdot_iM_CtT.resize(q.signed_distances.size(), q.N_CONTACTS);
  Cdot_iM_LT.resize(q.signed_distances.size(), q.N_LIMITS);

  // prepare iterators for contacts
  ColumnIteratord Cn_iter = Cdot_iM_CnT.column_iterator_begin();
  ColumnIteratord Cs_iter = Cdot_iM_CsT.column_iterator_begin();
  ColumnIteratord Ct_iter = Cdot_iM_CtT.column_iterator_begin();
  ColumnIteratord L_iter =  Cdot_iM_LT.column_iterator_begin();

  // for each pair of bodies
  for (unsigned k=0; k< q.signed_distances.size(); k++)
  {
    // get the two single bodies
    shared_ptr<SingleBodyd> s1 = q.signed_distances[k].a->get_single_body();
    shared_ptr<SingleBodyd> s2 = q.signed_distances[k].b->get_single_body();

    // get the two bodies involved
    shared_ptr<DynamicBodyd> sup1 = ImpactConstraintHandler::get_super_body(s1);
    shared_ptr<DynamicBodyd> sup2 = ImpactConstraintHandler::get_super_body(s2);
    tmp_supers1.clear();
    tmp_supers1.push_back(sup1);
    tmp_supers1.push_back(sup2);

    // sort the vector so we can do the intersection
    std::sort(tmp_supers1.begin(), tmp_supers1.end());

    // get the signed distance between the two bodies
    double phi = q.signed_distances[k].dist;

    // for each contact constraint 
    for (unsigned i=0; i< q.contact_constraints.size(); i++)
    {
      // zero the Cn, Cs, and Ct iterators
      *Cn_iter = 0.0;
      *Cs_iter = 0.0;
      *Ct_iter = 0.0;

      // see whether constraint will have any effect on this pair of bodies
      isect.clear();
      tmp_supers2.clear();
      q.contact_constraints[i]->get_super_bodies(std::back_inserter(tmp_supers2));

      // sort the vector so we can do the intersection
      std::sort(tmp_supers2.begin(), tmp_supers2.end());

      // do the intersection
      std::set_intersection(tmp_supers1.begin(), tmp_supers1.end(),
                            tmp_supers2.begin(), tmp_supers2.end(),
                            std::back_inserter(isect));
      if (isect.empty())
        continue;                      

      // apply a test impulse in the normal direction
      apply_impulse(*q.contact_constraints[i], 
                    q.contact_constraints[i]->contact_normal);

      // integrates bodies' positions forward
      integrate_positions(isect, DT);

      // compute the signed distance function
      double phi_new = calc_signed_dist(s1, s2);

      // set the appropriate entry of the Jacobian
      *Cn_iter = (phi_new - phi)/DT - q.Cdot_v[k];  Cn_iter++;

      // restore coordinates and velocities
      restore_coords_and_velocities(isect, gc_map, gv_map);

      // apply a test impulse in the first tangent direction
      apply_impulse(*q.contact_constraints[i], 
                    q.contact_constraints[i]->contact_tan1);

      // integrates bodies' positions forward
      integrate_positions(isect, DT);

      // compute the signed distance function
      phi_new = calc_signed_dist(s1, s2);

      // set the appropriate entry of the Jacobian
      *Cs_iter = (phi_new - phi)/DT - q.Cdot_v[k];  Cs_iter++;

      // restore coordinates and velocities
      restore_coords_and_velocities(isect, gc_map, gv_map);

      // apply a test impulse in the second tangent direction
      apply_impulse(*q.contact_constraints[i], 
                    q.contact_constraints[i]->contact_tan2);

      // integrates bodies' positions forward
      integrate_positions(isect, DT);

      // compute the signed distance function
      phi_new = calc_signed_dist(s1, s2);

      // set the appropriate entry of the Jacobian
      *Ct_iter = (phi_new - phi)/DT - q.Cdot_v[k];  Ct_iter++;

      // restore coordinates and velocities
      restore_coords_and_velocities(isect, gc_map, gv_map);
    }

    // for each limit constraint 
    for (unsigned i=0; i< q.limit_constraints.size(); i++)
    {
      // zero the LT iterator
      *L_iter = 0.0;

      // see whether constraint will have any effect on this pair of bodies
      isect.clear();
      tmp_supers2.clear();
      q.limit_constraints[i]->get_super_bodies(std::back_inserter(tmp_supers2));
      std::set_intersection(tmp_supers1.begin(), tmp_supers1.end(),
                            tmp_supers2.begin(), tmp_supers2.end(),
                            std::back_inserter(isect));
      if (isect.empty())
        continue;                      

      // apply a test impulse in the limit direction
      apply_impulse(*q.limit_constraints[i]);

      // integrates bodies' positions forward
      integrate_positions(isect, DT);

      // compute the signed distance function
      double phi_new = calc_signed_dist(s1, s2);

      // set the appropriate entry of the Jacobian
      *L_iter = (phi_new - phi)/DT - q.Cdot_v[k];  L_iter++;

      // restore coordinates and velocities
      restore_coords_and_velocities(isect, gc_map, gv_map);
    }
  }
}

/// Applies a test impulse to a limit contact
void SignedDistDot::apply_impulse(const UnilateralConstraint& contact_constraint, const Vector3d& dir)
{
  VectorNd gf1, gf2;

  // get the single bodies
  shared_ptr<SingleBodyd> s1 = contact_constraint.contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> s2 = contact_constraint.contact_geom2->get_single_body();

  // get the super bodies
  shared_ptr<DynamicBodyd> sb1 = ImpactConstraintHandler::get_super_body(s1); 
  shared_ptr<DynamicBodyd> sb2 = ImpactConstraintHandler::get_super_body(s2); 

  // setup the contact force
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = Origin3d(Pose3d::transform_point(GLOBAL, contact_constraint.contact_point)); 
  SForced f(P);
  f.set_force(Pose3d::transform_vector(P, dir));

  // convert the force to generalized forces
  sb1->convert_to_generalized_force(s1, f, gf1);
  sb2->convert_to_generalized_force(s2, -f, gf2);

  // apply the impulse
  sb1->apply_generalized_impulse(gf1);
  sb2->apply_generalized_impulse(gf2);
}

/// Applies a test impulse to a limit contact
void SignedDistDot::apply_impulse(const UnilateralConstraint& limit_constraint)
{
  // get the limit joint
  JointPtr j = limit_constraint.limit_joint;

  // get the articulated body
  shared_ptr<ArticulatedBody> ab = j->get_articulated_body();

  // get the direction of the impulse application
  double dir = (limit_constraint.limit_upper) ? -1.0 : 1.0;

  // setup a generalized force
  VectorNd gf(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));
  gf.set_zero();
  gf[j->get_index() + limit_constraint.limit_dof] = dir;

  // apply the impulse
  ab->apply_generalized_impulse(gf);
}

/// Restore coordinates and velocities
void SignedDistDot::restore_coords_and_velocities(const vector<shared_ptr<DynamicBodyd> >& isect, map<shared_ptr<DynamicBodyd>, VectorNd>& gc_map, map<shared_ptr<DynamicBodyd>, VectorNd>& gv_map)
{
  // restore the configurations and velocities of the involved bodies
  for (unsigned j=0; j< isect.size(); j++)
  {
    isect[j]->set_generalized_coordinates_euler(gc_map[isect[j]]);
    isect[j]->set_generalized_velocity(DynamicBodyd::eSpatial, gv_map[isect[j]]);
  }
}

// integrate the bodies in the intersection of the sets forward
void SignedDistDot::integrate_positions(const vector<shared_ptr<DynamicBodyd> >& isect, double dt)
{
  VectorNd gv, gc;

  for (unsigned j=0; j< isect.size(); j++)
  {
    isect[j]->get_generalized_coordinates_euler(gc);
    isect[j]->get_generalized_velocity(DynamicBodyd::eEuler, gv);
    gv *= dt;
    gc += gv;
    isect[j]->set_generalized_coordinates_euler(gc);
  }
}


