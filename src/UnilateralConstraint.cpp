/***************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <fstream>

#ifdef USE_OSG
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/MatrixTransform>
#include <osg/Material>
#endif

#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/RigidBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/Log.h>
#include <Moby/UnilateralConstraint.h>

using namespace Ravelin;
using namespace Moby;
using std::pair;
using std::list;
using std::vector;
using std::map;
using std::multimap;
using std::set;
using std::endl;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

// static declarations
MatrixNd UnilateralConstraint::J1, UnilateralConstraint::J2, UnilateralConstraint::workM1, UnilateralConstraint::workM2;
MatrixNd UnilateralConstraint::JJ, UnilateralConstraint::J, UnilateralConstraint::Jx, UnilateralConstraint::Jy, UnilateralConstraint::dJ1, UnilateralConstraint::dJ2;
VectorNd UnilateralConstraint::v, UnilateralConstraint::workv, UnilateralConstraint::workv2; 

/// Creates an empty constraint 
UnilateralConstraint::UnilateralConstraint()
{
  _contact_frame = shared_ptr<Pose3d>(new Pose3d);
  tol = NEAR_ZERO;              // default collision tolerance
  compliance = eRigid;
  constraint_type = eNone;
  signed_violation = 0.0;
  limit_dof = std::numeric_limits<unsigned>::max();
  limit_epsilon = (double) 0.0;
  limit_upper = false;
  limit_impulse = (double) 0.0;
  contact_normal.set_zero(GLOBAL);
  contact_impulse.set_zero(GLOBAL);
  contact_point.set_zero(GLOBAL);
  contact_mu_coulomb = (double) 0.0;
  contact_mu_viscous = (double) 0.0;
  contact_penalty_Kp = (double) 0.0;
  contact_penalty_Kv = (double) 0.0;
  contact_epsilon = (double) 0.0;
  contact_NK = 4;
  compliance = eCompliant;
}

UnilateralConstraint& UnilateralConstraint::operator=(const UnilateralConstraint& e)
{
  tol = e.tol;
  signed_violation = e.signed_violation;
  constraint_type = e.constraint_type;
  compliance = e.compliance;
  limit_epsilon = e.limit_epsilon;
  limit_dof = e.limit_dof;
  limit_upper = e.limit_upper;
  limit_impulse = e.limit_impulse;
  limit_joint = e.limit_joint;
  contact_normal = e.contact_normal;
  contact_geom1 = e.contact_geom1;
  contact_geom2 = e.contact_geom2;
  contact_point = e.contact_point;
  contact_impulse = e.contact_impulse;
  contact_mu_coulomb = e.contact_mu_coulomb;
  contact_mu_viscous = e.contact_mu_viscous;
  contact_penalty_Kp = e.contact_penalty_Kp;
  contact_penalty_Kv = e.contact_penalty_Kv;
  contact_epsilon = e.contact_epsilon;
  contact_NK = e.contact_NK;
  contact_tan1 = e.contact_tan1;
  contact_tan2 = e.contact_tan2;

  return *this;
}

/// Computes the constraint data
void UnilateralConstraint::compute_constraint_data(MatrixNd& M, VectorNd& q) const
{
  if (constraint_type == eContact)
  {
    // setup useful indices
    const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

    // get the two single bodies
    shared_ptr<SingleBodyd> sb1 = contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> su1 = dynamic_pointer_cast<DynamicBodyd>(sb1->get_super_body());
    shared_ptr<DynamicBodyd> su2 = dynamic_pointer_cast<DynamicBodyd>(sb2->get_super_body());

    // verify the contact point, normal, and tangents are in the global frame
    assert(contact_point.pose == GLOBAL);
    assert(contact_normal.pose == GLOBAL);
    assert(contact_tan1.pose == GLOBAL);
    assert(contact_tan2.pose == GLOBAL);

    // setup the contact frame
    _contact_frame->q.set_identity();
    _contact_frame->x = contact_point;

    // get the numbers of generalized coordinates for the two super bodies
    const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBodyd::eSpatial);
    const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // resize the Jacobians 
    J1.set_zero(THREE_D, NGC1);
    J2.set_zero(THREE_D, NGC2);

    // get the directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

    // setup a matrix of contact directions
    // R' transforms contact orientation to global orientation
    Matrix3d R;
    R.set_column(N, normal);
    R.set_column(S, tan1);
    R.set_column(T, tan2);

    // compute the Jacobians for the two bodies; Jacobian transforms velocities
    // in mixed frame to velocities in contact frame
    su1->calc_jacobian(su1->get_gc_pose(), _contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin1, J1);
    su2->calc_jacobian(su2->get_gc_pose(), _contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin2, J2);

    // compute the constraint inertia matrix for the first body
    su1->transpose_solve_generalized_inertia(J1, workM1);
    J1.mult(workM1, M);

    // compute the constraint inertia matrix for the second body
    su2->transpose_solve_generalized_inertia(J2, workM1);
    J2.mult(workM1, workM2);
    M += workM2;

    // compute the constraint velocity
    su1->get_generalized_velocity(DynamicBodyd::eSpatial, v);
    J1.mult(v, q);

    // free v1 and allocate v2 and workv
    su2->get_generalized_velocity(DynamicBodyd::eSpatial, v);
    q += J2.mult(v, workv);
  }
  else if (constraint_type == eLimit)
  {
    // get the super body
    ArticulatedBodyPtr ab = limit_joint->get_articulated_body();
    RCArticulatedBodyPtr su = dynamic_pointer_cast<RCArticulatedBody>(ab);

    // case 1: reduced-coordinate articulated body
    if (su)
    {
      // determine the joint limit index
      unsigned idx = limit_joint->get_coord_index() + limit_dof;

      // setup a vector to solve
      v.set_zero(su->num_generalized_coordinates(DynamicBodyd::eSpatial));
      v[idx] = 1.0;

      // solve
      su->solve_generalized_inertia(v, workv); 
      M.resize(1,1);
      M(0,0) = workv[idx];
    }
    else
    {
      // TODO: handle absolute coordinate articulated bodies here
      // note: to do this constraint handler also needs to setup constraint Jac
      //       as an equality constraint

      // setup joint velocity Jacobian here (Dx)

      // we need to compute:
      // | M  Jx' | x | delta xd | = | j |
      // | Jx 0   |   | lambda   | = | 0 |
      // such that:
      // Dx*xd^+ >= 0

      // 
    }

    // get the joint velocity
    q.resize(1);
    q[0] = limit_joint->qd[limit_dof];

    // if we're at an upper limit, negate q
    if (limit_upper)
      q.negate(); 
  }
} 

/// Determines whether two constraints are linked
bool UnilateralConstraint::is_linked(const UnilateralConstraint& e1, const UnilateralConstraint& e2)
{
  if (e1.constraint_type == eContact)
  {
    // get the two single bodies
    shared_ptr<SingleBodyd> e1sb1 = e1.contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> e1sb2 = e1.contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> e1s1 = dynamic_pointer_cast<DynamicBodyd>(e1sb1->get_super_body());
    shared_ptr<DynamicBodyd> e1s2 = dynamic_pointer_cast<DynamicBodyd>(e1sb2->get_super_body());

    // examine against other constraint type
    if (e2.constraint_type == eContact)
    {
      // get the two single bodies
      shared_ptr<SingleBodyd> e2sb1 = e2.contact_geom1->get_single_body();
      shared_ptr<SingleBodyd> e2sb2 = e2.contact_geom2->get_single_body();

      // get the two super bodies
      shared_ptr<DynamicBodyd> e2s1 = dynamic_pointer_cast<DynamicBodyd>(e2sb1->get_super_body());
      shared_ptr<DynamicBodyd> e2s2 = dynamic_pointer_cast<DynamicBodyd>(e2sb2->get_super_body());

      // see whether there are any bodies in common
      return e1s1 == e2s1 || e1s1 == e2s2 || e1s2 == e2s1 || e1s2 == e2s2;
    }
    else if (e2.constraint_type == eLimit)
    {
      ArticulatedBodyPtr ab = e2.limit_joint->get_articulated_body();
      return e1s1 == ab || e1s2 == ab; 
    }
    else
    {
      assert(false);

      // even though we shouldn't be here, we'll return true (it's conservative)
      return true;
    }
  }
  else if (e1.constraint_type == eLimit)
  {
    if (e2.constraint_type == eContact)
      return is_linked(e2, e1);
    else if (e2.constraint_type == eLimit)
    {
      ArticulatedBodyPtr ab1 = e1.limit_joint->get_articulated_body();
      ArticulatedBodyPtr ab2 = e2.limit_joint->get_articulated_body();
      return ab1 == ab2;
    }
    else
    {
      assert(false);

      // even though we shouldn't be here, we'll return true (it's conservative)
      return true;
    }
  }
  else
  {
    assert(false);

    // even though we shouldn't be here, we'll return true (it's conservative)
    return true;
  }
}

/// Updates the constraint data
void UnilateralConstraint::compute_cross_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // verify that the constraints are linked
  if (!is_linked(*this, e))
    return;
    
  switch (constraint_type)
  {
    case eContact:
      switch (e.constraint_type)
      {
        case eContact: compute_cross_contact_contact_constraint_data(e, M); break;
        case eLimit:   compute_cross_contact_limit_constraint_data(e, M); break;
        case eNone:    M.resize(0,0); break;
      }
      break;

    case eLimit:
      switch (e.constraint_type)
      {
        case eContact: compute_cross_limit_contact_constraint_data(e, M); break;
        case eLimit:   compute_cross_limit_limit_constraint_data(e, M); break;
        case eNone:    M.resize(0,0); break;
      }
      break;

    case eNone:
      M.resize(0,0);
      break;
  }
}

/// Updates contact/contact cross constraint data
/**
 * From two contact points, we can have up to three separate super bodies. 
 */
void UnilateralConstraint::compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // get the unique super bodies
  shared_ptr<DynamicBodyd> bodies[4];
  shared_ptr<DynamicBodyd>* end = get_super_bodies(bodies);
  end = e.get_super_bodies(end);
  std::sort(bodies, end);
  end = std::unique(bodies, end);

  // determine how many unique super bodies we have
  const unsigned NSUPER = end - bodies;

  // clear M
  M.set_zero(3,3);

  // if we have exactly two super bodies, process them individually
  if (NSUPER == 1)
    compute_cross_contact_contact_constraint_data(e, M, bodies[0]);
  if (NSUPER == 2)
  {
    compute_cross_contact_contact_constraint_data(e, M, bodies[0]);
    compute_cross_contact_contact_constraint_data(e, M, bodies[1]);
  }
  else if (NSUPER == 3)
  {
    // find the one common super body
    shared_ptr<DynamicBodyd> bodies1[2], bodies2[2], isect[1];
    shared_ptr<DynamicBodyd>* end1 = get_super_bodies(bodies1);
    shared_ptr<DynamicBodyd>* end2 = e.get_super_bodies(bodies2);
    std::sort(bodies1, end1);
    std::sort(bodies2, end2);
    shared_ptr<DynamicBodyd>* isect_end = std::set_intersection(bodies1, end1, bodies2, end2, isect);
    assert(isect_end - isect == 1);
    compute_cross_contact_contact_constraint_data(e, M, isect[0]);
  }
  else if (NSUPER == 4)
    assert(false);
}

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, MatrixNd& M, shared_ptr<DynamicBodyd> su) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the first two single bodies
  shared_ptr<SingleBodyd> sba1 = contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sba2 = contact_geom2->get_single_body();

  // get the first two super bodies
  shared_ptr<DynamicBodyd> sua1 = dynamic_pointer_cast<DynamicBodyd>(sba1->get_super_body());
  shared_ptr<DynamicBodyd> sua2 = dynamic_pointer_cast<DynamicBodyd>(sba2->get_super_body());

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // resize Jacobian 
  J.resize(THREE_D, NGC);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = contact_point;

  // get the directions in the constraint frame
  Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
  Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
  Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

  // setup a matrix of contact directions
  Matrix3d R;
  R.set_column(N, normal);
  R.set_column(S, tan1);
  R.set_column(T, tan2);

  // compute the Jacobians, checking to see whether necessary
  if (sua1 == su)
  {
    su->calc_jacobian(su->get_gc_pose(), _contact_frame, sba1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, J);
    compute_cross_contact_contact_constraint_data(e, M, su, J);
  }
  if (sua2 == su)
  {
    su->calc_jacobian(su->get_gc_pose(), _contact_frame, sba2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, J);
    compute_cross_contact_contact_constraint_data(e, M, su, J);
  }
} 

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, MatrixNd& M, shared_ptr<DynamicBodyd> su, const MatrixNd& J) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the second two single bodies
  shared_ptr<SingleBodyd> sbb1 = e.contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sbb2 = e.contact_geom2->get_single_body();

  // get the second two super bodies
  shared_ptr<DynamicBodyd> sub1 = dynamic_pointer_cast<DynamicBodyd>(sbb1->get_super_body());
  shared_ptr<DynamicBodyd> sub2 = dynamic_pointer_cast<DynamicBodyd>(sbb2->get_super_body());

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // resize Jacobian 
  Jx.resize(THREE_D, NGC);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = e.contact_point;

  // get the directions in the constraint frame
  Vector3d normal = Pose3d::transform_vector(_contact_frame, e.contact_normal);
  Vector3d tan1 = Pose3d::transform_vector(_contact_frame, e.contact_tan1);
  Vector3d tan2 = Pose3d::transform_vector(_contact_frame, e.contact_tan2);

  // setup a matrix of contact directions
  Matrix3d R;
  R.set_column(N, normal);
  R.set_column(S, tan1);
  R.set_column(T, tan2);

  // compute the Jacobians, checking to see whether necessary
  if (sub1 == su)
  {
    // first compute the Jacobian
    su->calc_jacobian(su->get_gc_pose(), _contact_frame, sbb1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, Jx);

    // now update M 
    su->transpose_solve_generalized_inertia(Jx, workM1);
    M += J.mult(workM1, workM2);
  }
  if (sub2 == su)
  {
    su->calc_jacobian(su->get_gc_pose(), _contact_frame, sbb2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, Jx);

    // now update M
    su->transpose_solve_generalized_inertia(Jx, workM1);
    M += J.mult(workM1, workM2);
  }
}

/// Updates contact/limit cross constraint data
void UnilateralConstraint::compute_cross_contact_limit_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the articulated body of the constraint
  ArticulatedBodyPtr ab = e.limit_joint->get_articulated_body();
  RCArticulatedBodyPtr su = dynamic_pointer_cast<RCArticulatedBody>(ab);
  assert(su);

  // get the index of the limit joint
  unsigned idx = e.limit_joint->get_coord_index() + e.limit_dof;

  // get the two single bodies
  shared_ptr<SingleBodyd> sb1 = contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sb2 = contact_geom2->get_single_body();

  // get the two super bodies
  shared_ptr<DynamicBodyd> su1 = dynamic_pointer_cast<DynamicBodyd>(sb1->get_super_body());
  shared_ptr<DynamicBodyd> su2 = dynamic_pointer_cast<DynamicBodyd>(sb2->get_super_body());

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = contact_point;
  _contact_frame->rpose = GLOBAL;

  // get the directions in the constraint frame
  Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
  Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
  Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

  // setup a matrix of contact directions
  Matrix3d R;
  R.set_column(N, normal);
  R.set_column(S, tan1);
  R.set_column(T, tan2);

  // get the numbers of generalized coordinates for the two super bodies
  const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBodyd::eSpatial);
  const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // see whether limit is equal to su1
  if (su == su1)
  {
    // resize Jacobian
    J1.resize(THREE_D, NGC1);

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(su1->get_gc_pose(), _contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, J1);

    // compute the constraint inertia matrix for the first body
    su1->transpose_solve_generalized_inertia(J1, workM1);

    // get the appropriate row of workM
    if (M.columns() == 3)
      M.row(0) = workM1.row(idx);
    else
      M.column(0) = workM1.row(idx);

    // determine whether to negate the row
    if (e.limit_upper)
      M.negate();
  }
  else
    // setup M
    M.set_zero(1, 3);

  // handle case of articulated body equal to contact one super body 
  if (ab == su2)
  {
    // resize Jacobian
    J1.resize(THREE_D, NGC2);

    // compute the Jacobians for the two bodies
    su2->calc_jacobian(su2->get_gc_pose(), _contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, J1);

    // compute the constraint inertia matrix for the first body
    su2->transpose_solve_generalized_inertia(J1, workM1);

    // get the appropriate row of workM
    if (M.columns() == 3)
      M.row(0) += workM1.row(idx); 
    else
      M.column(0) = workM1.row(idx);
    // determine whether to negate the row
    if (e.limit_upper)
      M.negate();
  }
} 

/// Updates limit/contact cross constraint data
void UnilateralConstraint::compute_cross_limit_contact_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // compute the cross constraint data
  e.compute_cross_contact_limit_constraint_data(*this, workM2);

  // transpose the matrix
  MatrixNd::transpose(workM2, M);
} 

/// Updates limit/limit cross constraint data
void UnilateralConstraint::compute_cross_limit_limit_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // get the super body
  ArticulatedBodyPtr ab = limit_joint->get_articulated_body();
  RCArticulatedBodyPtr su = dynamic_pointer_cast<RCArticulatedBody>(ab);
  assert(su);

  // determine the joint limit indices
  unsigned idx1 = limit_joint->get_coord_index() + limit_dof;
  unsigned idx2 = e.limit_joint->get_coord_index() + e.limit_dof;

  // case 1: reduced-coordinate articulated body
  if (su)
  {
    // setup a vector to solve
    workv.set_zero(su->num_generalized_coordinates(DynamicBodyd::eSpatial));
    workv[idx1] = 1.0;

    // solve
    su->solve_generalized_inertia(workv, workv2); 

    // determine whether to negate
    double value = workv2[idx2];
    if ((limit_upper && !e.limit_upper) ||
        (!limit_upper && e.limit_upper))
      value = -value;

    // setup M
    M.resize(1,1);
    M.data()[0] = value;
  }
  else
  {
      // TODO: handle absolute coordinate articulated bodies here
      // note: to do this constraint handler also needs to setup constraint Jac
      //       as an equality constraint

      // setup joint velocity Jacobian here (Dx)

      // we need to compute:
      // | M  Jx' | x | delta xd | = | j |
      // | Jx 0   |   | lambda   | = | 0 |
      // such that:
      // Dx*xd^+ >= 0

      // 
  }
} 

/// Sets the contact parameters for this constraint
void UnilateralConstraint::set_contact_parameters(const ContactParameters& cparams)
{
  contact_mu_coulomb = cparams.mu_coulomb;
  contact_mu_viscous = cparams.mu_viscous;
  contact_penalty_Kp = cparams.penalty_Kp;
  contact_penalty_Kv = cparams.penalty_Kv;
  contact_epsilon = cparams.epsilon;
  contact_NK = cparams.NK;

  // redetermine contact tangents
  determine_contact_tangents();
  
  assert(contact_NK >= 4);
}

double calc_constraint_vel2(const UnilateralConstraint& e)
{
  assert (e.constraint_type == UnilateralConstraint::eContact);
  shared_ptr<SingleBodyd> sba = e.contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sbb = e.contact_geom2->get_single_body();

  // get the vels 
  const SVelocityd& va = sba->get_velocity(); 
  const SVelocityd& vb = sbb->get_velocity(); 

  // get the bodies as rigid bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(sba);
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(sbb);

  // transform velocity to mixed frames
  shared_ptr<const Pose3d> Pa = rba->get_mixed_pose(); 
  shared_ptr<const Pose3d> Pb = rbb->get_mixed_pose(); 
  SVelocityd ta = Pose3d::transform(Pa, va);
  SVelocityd tb = Pose3d::transform(Pb, vb);

  // transform normal to mixed frame
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = e.contact_point;
  P->rpose = GLOBAL;
  Vector3d normal = Pose3d::transform_vector(P, e.contact_normal);

  // get the linear velocities and project against the normal
  Vector3d ra(e.contact_point - Pa->x);
  Vector3d rb(e.contact_point - Pb->x);
  Vector3d xda = ta.get_linear(); 
  Vector3d xdb = tb.get_linear(); 
  Vector3d wa = ta.get_angular();
  Vector3d wb = tb.get_angular();
  ra.pose = GLOBAL;
  rb.pose = GLOBAL;
  xda.pose = GLOBAL;
  xdb.pose = GLOBAL;
  wa.pose = GLOBAL;
  wb.pose = GLOBAL;
  Vector3d v(xda - xdb + Vector3d::cross(wa, ra) - Vector3d::cross(wb, rb));
  v.pose = normal.pose;
  return v.dot(normal);
}

/// Computes the velocity of this constraint
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
double UnilateralConstraint::calc_constraint_vel() const
{
  if (constraint_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    shared_ptr<SingleBodyd> sba = contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 

    // setup the constraint frame
    _contact_frame->x = contact_point;
    _contact_frame->q.set_identity();
    _contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(_contact_frame);
    SVelocityd ta = Pose3d::transform(const_contact_frame, va); 
    SVelocityd tb = Pose3d::transform(const_contact_frame, vb); 

    // get the contact normal in the correct pose
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);

    FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::calc_constraint_vel() entered" << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "normal (constraint frame): " << normal << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "tangent 1 (constraint frame): " << Pose3d::transform_vector(_contact_frame, contact_tan1) << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "tangent 2 (constraint frame): " << Pose3d::transform_vector(_contact_frame, contact_tan2) << std::endl;
/*
    FILE_LOG(LOG_CONSTRAINT) << "spatial velocity (mixed frame) for body A: " << Pose3d::transform(dynamic_pointer_cast<RigidBody>(sba)->get_mixed_pose(), ta) << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "spatial velocity (constraint frame) for body A: " << ta << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "spatial velocity (mixed frame) for body B: " << Pose3d::transform(dynamic_pointer_cast<RigidBody>(sbb)->get_mixed_pose(), tb) << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "spatial velocity (constraint frame) for body B: " << tb << std::endl;
*/
    FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::calc_constraint_vel() exited" << std::endl;

    // get the linear velocities and project against the normal
    //assert(std::fabs(normal.dot(ta.get_linear() - tb.get_linear())) < NEAR_ZERO || (std::fabs(normal.dot(ta.get_linear() - tb.get_linear()) - calc_constraint_vel2(*this)))/std::fabs(normal.dot(ta.get_linear() - tb.get_linear())) < NEAR_ZERO);
    return normal.dot(ta.get_linear() - tb.get_linear());
  }
  else if (constraint_type == eLimit)
  {
    double qd = limit_joint->qd[limit_dof];
    return (limit_upper) ? -qd : qd;
  }
  else
  {
    assert(false);
    return 0.0;
  }
}  

/// Sends the constraint to the specified stream
std::ostream& Moby::operator<<(std::ostream& o, const UnilateralConstraint& e)
{
  switch (e.constraint_type)
  {
    case UnilateralConstraint::eNone:
      o << "(constraint type: none)" << std::endl;
      return o;

    case UnilateralConstraint::eLimit:
      o << "(constraint type: joint limit)" << std::endl;
      break;

    case UnilateralConstraint::eContact:
      o << "(constraint type: contact)" << std::endl;
      break;
  }
	 
  if (e.compliance == UnilateralConstraint::eRigid) 
    o << "compliance: rigid" << std::endl;
  else
    o << "compliance: compliant, kp = " << e.contact_penalty_Kp << ", kv = " << e.contact_penalty_Kv << std::endl;

  if (e.constraint_type == UnilateralConstraint::eLimit)
  {
    o << "limit joint ID: " << e.limit_joint->id << std::endl;
    o << "limit joint coordinate index: " << e.limit_joint->get_coord_index() << std::endl;
    o << "limit joint DOF: " << e.limit_dof << std::endl;
    o << "upper limit? " << e.limit_upper << std::endl;
    o << "limit velocity: " << e.calc_constraint_vel() << std::endl;
  }
  else if (e.constraint_type == UnilateralConstraint::eContact)
  {
    if (e.contact_geom1)
    {
      shared_ptr<SingleBodyd> sb1(e.contact_geom1->get_single_body());
      if (sb1)
      {
	o << "body1: " << sb1->body_id << std::endl;
      }
      else
	o << "body1: (undefined)" << std::endl;
    }
    else
      o << "geom1: (undefined)" << std::endl;
  
    if (e.contact_geom2)
    {
      shared_ptr<SingleBodyd> sb2(e.contact_geom2->get_single_body());
      if (sb2)
      {
	o << "body2: " << sb2->body_id << std::endl;
      }    
      else
	o << "body2: (undefined)" << std::endl;
    }
    else
      o << "geom2: (undefined)" << std::endl;

    o << "contact point / normal pose: " << ((e.contact_point.pose) ? Pose3d(*e.contact_point.pose).update_relative_pose(GLOBAL) : GLOBAL) << std::endl;
    o << "contact point: " << e.contact_point << std::endl;
    shared_ptr<SingleBodyd> sba = e.contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sbb = e.contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 

    // setup the constraint frame
    shared_ptr<Pose3d> constraint_frame(new Pose3d);
    constraint_frame->x = e.contact_point;
    constraint_frame->q.set_identity();
    constraint_frame->rpose = GLOBAL;
    shared_ptr<const Pose3d> const_constraint_frame = boost::const_pointer_cast<const Pose3d>(constraint_frame);

    // compute the velocities at the contact point
    SVelocityd ta = Pose3d::transform(const_constraint_frame, va); 
    SVelocityd tb = Pose3d::transform(const_constraint_frame, vb); 

    // get the contact normal in the correct pose
    Vector3d normal = Pose3d::transform_vector(constraint_frame, e.contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(constraint_frame, e.contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(constraint_frame, e.contact_tan2);

    // get the linear velocities and project against the normal
    Vector3d rvlin = ta.get_linear() - tb.get_linear();
    //assert(std::fabs(normal.dot(rvlin)) < NEAR_ZERO || std::fabs(normal.dot(rvlin) - calc_constraint_vel2(e))/std::fabs(normal.dot(rvlin)) < NEAR_ZERO);
    o << "relative normal velocity: " << normal.dot(rvlin) << std::endl;
    o << "relative tangent 1 velocity: " << tan1.dot(rvlin) << std::endl;
    o << "relative tangent 2 velocity: " << tan2.dot(rvlin) << std::endl;
    o << "calc_constraint_vel() reports: " << std::endl;
    e.calc_constraint_vel();
  }

  return o;
}

#ifdef USE_OSG
/// Copies this matrix to an OpenSceneGraph Matrixd object
static void to_osg_matrix(const Pose3d& src, osg::Matrixd& tgt)
{
  // get the rotation matrix
  Matrix3d M = src.q;

  // setup the rotation components of tgt
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = M(i,j);

  // setup the translation components of tgt
  for (unsigned i=X; i<= Z; i++)
    tgt(W,i) = src.x[i];

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (double) 0.0;
  tgt(W,W) = (double) 1.0;
}
#endif

/// Makes a contact visualizable
osg::Node* UnilateralConstraint::to_visualization_data() const
{
  #ifdef USE_OSG
  const float CONE_HEIGHT = .2f;
  const float CONE_RADIUS = .2f;
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the transformation matrix for the cone
  Vector3d x_axis, z_axis;
  Vector3d::determine_orthonormal_basis(contact_normal, x_axis, z_axis);
  Matrix3d R;
  R.set_column(X, x_axis);
  R.set_column(Y, contact_normal);
  R.set_column(Z, -z_axis);
  Vector3d x = contact_point + contact_normal;
  Pose3d T;
  T.q = R;
  T.x = Origin3d(x);

  // setup the transform node for the cone
  osg::Matrixd m;
  to_osg_matrix(T, m);
  osg::MatrixTransform* transform = new osg::MatrixTransform;
  transform->setMatrix(m);

  // create the new color
  osg::Material* mat = new osg::Material;
  const float RED = (float) rand() / RAND_MAX;
  const float GREEN = (float) rand() / RAND_MAX;
  const float BLUE = (float) rand() / RAND_MAX;
  mat->setColorMode(osg::Material::DIFFUSE);
  mat->setDiffuse(osg::Material::FRONT, osg::Vec4(RED, GREEN, BLUE, 1.0f));
  transform->getOrCreateStateSet()->setAttribute(mat);

  // create the line
  osg::Geometry* linegeom = new osg::Geometry;
  osg::Vec3Array* varray = new osg::Vec3Array;
  linegeom->setVertexArray(varray);  
  varray->push_back(osg::Vec3((float) contact_point[X], (float) contact_point[Y], (float) contact_point[Z]));
  varray->push_back(osg::Vec3((float) contact_point[X] + (float) contact_normal[X], (float) contact_point[Y] + (float) contact_normal[Y], (float) contact_point[Z] + (float) contact_normal[Z]));
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(linegeom);

  // create the cone
  osg::Cone* cone = new osg::Cone;
  cone->setRadius(CONE_RADIUS);
  cone->setHeight(CONE_HEIGHT);
  geode->addDrawable(new osg::ShapeDrawable(cone));

  // add the geode
  transform->addChild(geode);

  return transform;
  #else
  return NULL;
  #endif
}

/// Given a vector of constraints, determines all of the sets of connected constraints
/**
 * A set of connected constraints is the set of all constraints such that, for a
 * given constraint A in the set, there exists another constraint B for which A
 * and B share at least one rigid body.  
 * \param constraints the list of constraints
 * \param groups the islands of connected constraints on return
 */
void UnilateralConstraint::determine_connected_constraints(const vector<UnilateralConstraint>& constraints, const vector<JointPtr>& implicit_joints, list<list<UnilateralConstraint*> >& groups, list<list<shared_ptr<DynamicBodyd> > >& remaining_islands)
{
  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_connected_contacts() entered" << std::endl;

  // clear the groups
  groups.clear();
  remaining_islands.clear();

  // copy the list of constraints -- only ones with geometry
  list<UnilateralConstraint*> constraints_copy;
  BOOST_FOREACH(const UnilateralConstraint& e, constraints)
    if (e.constraint_type != UnilateralConstraint::eNone)
      constraints_copy.push_back((UnilateralConstraint*) &e);
  
  // The way that we'll determine the constraint islands is to treat each rigid
  // body present in the constraints as a node in a graph; nodes will be connected
  // to other nodes if (a) they are both present in constraint or (b) they are
  // part of the same articulated body.  Nodes will not be created for disabled
  // bodies.
  set<shared_ptr<SingleBodyd> > nodes;
  multimap<shared_ptr<SingleBodyd>, shared_ptr<SingleBodyd> > edges;
  typedef multimap<shared_ptr<SingleBodyd>, shared_ptr<SingleBodyd> >::const_iterator EdgeIter;

  // get all single bodies present in the unilateral constraints
  for (list<UnilateralConstraint*>::const_iterator i = constraints_copy.begin(); i != constraints_copy.end(); i++)
  {
    if ((*i)->constraint_type == UnilateralConstraint::eContact)
    {
      shared_ptr<SingleBodyd> sb1((*i)->contact_geom1->get_single_body());
      shared_ptr<SingleBodyd> sb2((*i)->contact_geom2->get_single_body());
      if (sb1->is_enabled())
        nodes.insert(sb1);
      if (sb2->is_enabled())
        nodes.insert(sb2);
      if (sb1->is_enabled() && sb2->is_enabled())
      {
        edges.insert(std::make_pair(sb1, sb2));
        edges.insert(std::make_pair(sb2, sb1));
      }
    }
    else if ((*i)->constraint_type == UnilateralConstraint::eLimit)
    {
      shared_ptr<RigidBodyd> inboard = (*i)->limit_joint->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = (*i)->limit_joint->get_outboard_link();
      nodes.insert(inboard);
      nodes.insert(outboard);
    }
    else 
      assert(false);
  }

  // get all single bodies present in the implicit bilateral constraints
  // add nodes and create edges between them
  for (unsigned i=0; i< implicit_joints.size(); i++)
  {
    JointPtr j = implicit_joints[i];
    shared_ptr<RigidBodyd> inboard = j->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = j->get_outboard_link();
    if (inboard->is_enabled() && outboard->is_enabled())
    {
      nodes.insert(inboard);
      nodes.insert(outboard);
      edges.insert(std::make_pair(inboard, outboard));
      edges.insert(std::make_pair(outboard, inboard));
    }
  } 

  FILE_LOG(LOG_CONSTRAINT) << " -- single bodies in constraints:" << std::endl;
  if (LOGGING(LOG_CONSTRAINT))
    for (set<shared_ptr<SingleBodyd> >::const_iterator i = nodes.begin(); i != nodes.end(); i++)
      FILE_LOG(LOG_CONSTRAINT) << "    " << (*i)->body_id << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << std::endl;

  // add connections between articulated rigid bodies -- NOTE: don't process
  // articulated bodies twice!
  set<shared_ptr<ArticulatedBodyd> > ab_processed;
  BOOST_FOREACH(shared_ptr<SingleBodyd> sb, nodes)
  {
    // if the body is not part of an articulated body, skip it
    shared_ptr<ArticulatedBodyd> abody = sb->get_articulated_body();
    if (!abody)
      continue;

    // see whether it has already been processed
    if (ab_processed.find(abody) != ab_processed.end())
      continue;

    // indicate that the articulated body will now have been processed
    ab_processed.insert(abody);

     // process all implicit joints in the articulated body
    const vector<shared_ptr<Jointd> >& implicit_joints = abody->get_implicit_joints();
    for (unsigned i=0; i< implicit_joints.size(); i++)
    {
      shared_ptr<Jointd> j = implicit_joints[i];
      shared_ptr<RigidBodyd> inboard = j->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = j->get_outboard_link();
      if (inboard->is_enabled() && outboard->is_enabled())
      {
        nodes.insert(inboard);
        nodes.insert(outboard);
        edges.insert(std::make_pair(inboard, outboard));
        edges.insert(std::make_pair(outboard, inboard));
      }
    }   

    // get all links in the articulated body
    const vector<shared_ptr<RigidBodyd> >& links = abody->get_links();

    // add edges between all pairs for which there are links
    vector<shared_ptr<RigidBodyd> >::const_iterator j, k;
    for (j = links.begin(); j != links.end(); j++)
    {
      // no sense iterating over all other links if link pointed to by j is
      // not a node
      if (nodes.find(*j) == nodes.end())
        continue;

      // iterate over all other nodes
      k = j;
      for (k++; k != links.end(); k++)
        if (nodes.find(*k) != nodes.end())
        {
          edges.insert(std::make_pair(*j, *k));
          edges.insert(std::make_pair(*k, *j));
        }
    }      
  }

  // Now, we'll remove nodes from the set until there are no more nodes.
  // For each removed node, we'll get add all constraints that contain the single 
  // body to the group; all neighboring nodes will then be processed.
  while (!nodes.empty())
  {
    // get the node from the front
    shared_ptr<SingleBodyd> node = *nodes.begin();

    // add a list to the contact groups
    groups.push_back(list<UnilateralConstraint*>());
    FILE_LOG(LOG_CONSTRAINT) << " -- constraints in group: " << std::endl;

    // create a node queue, with this node added
    std::queue<shared_ptr<SingleBodyd> > node_q;
    node_q.push(node);

    // setup a set of processed nodes
    std::set<shared_ptr<SingleBodyd> > processed_nodes;

    // loop until the queue is empty
    while (!node_q.empty())
    {
      // get the node off of the front of the node queue
      node = node_q.front();
      node_q.pop();

      // indicate that the node has now been processed
      processed_nodes.insert(node);

      // add all neighbors of the node that have not been processed already 
      // to the node queue
      std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
      for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
        if (processed_nodes.find(i->second) == processed_nodes.end())
          node_q.push(i->second);

      // loop through all remaining constraints
      for (list<UnilateralConstraint*>::iterator i = constraints_copy.begin(); i != constraints_copy.end(); )
      {
        if ((*i)->constraint_type == UnilateralConstraint::eContact)
        {
          shared_ptr<SingleBodyd> sb1((*i)->contact_geom1->get_single_body());
          shared_ptr<SingleBodyd> sb2((*i)->contact_geom2->get_single_body());

          // see whether one of the bodies is equal to the node
          if (sb1 == node || sb2 == node)
          {
            groups.back().push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->constraint_type == UnilateralConstraint::eLimit)
        {
          RigidBodyPtr inboard = (*i)->limit_joint->get_inboard_link();
          RigidBodyPtr outboard = (*i)->limit_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            groups.back().push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else
          assert(false);
      }

      // if no unilateral constraints have been added, add to remaining islands
      if (groups.back().empty())
      {
        // don't need an empty group of unilateral constraints
        groups.pop_back();

        // create a new island
        remaining_islands.push_back(list<shared_ptr<DynamicBodyd> >());

        // create a secondary node q and secondary processing set
        std::queue<shared_ptr<SingleBodyd> > node_q2;
        std::set<shared_ptr<SingleBodyd> > processed_nodes2;        

        // add the node to the queue
        node_q2.push(node);

        // find all connected bodies
        while (!node_q2.empty())
        {
          // get the node off of the front of the node queue
          node = node_q2.front();
          node_q2.pop();

          // indicate that the node has now been processed
          processed_nodes2.insert(node);

          // add the super body of this node to the island
          remaining_islands.back().push_back(node->get_super_body());

          // add all neighbors of the node that have not been processed already 
          // to the node queue
          std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
          for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
            if (processed_nodes2.find(i->second) == processed_nodes2.end())
              node_q2.push(i->second);
        }

        // finally, make the island of super bodies unique
        list<shared_ptr<DynamicBodyd> >& island = remaining_islands.back();
        island.sort();
        island.erase(std::unique(island.begin(), island.end()), island.end());
      } 
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_connected_constraints() exited" << std::endl;
}

/// Removes groups of contacts that contain no active contacts 
void UnilateralConstraint::remove_inactive_groups(list<list<UnilateralConstraint*> >& groups)
{
  typedef list<list<UnilateralConstraint*> >::iterator ListIter;

  for (ListIter i = groups.begin(); i != groups.end(); )
  {
    // look for impact in list i
    bool active_detected = false;
    BOOST_FOREACH(UnilateralConstraint* e, *i)
    {
      if (e->determine_constraint_class() == UnilateralConstraint::eNegative)
      {
        active_detected = true;
        break;
      }
    }

    // if no active constraint in the list, remove the list
    if (!active_detected)
    {
      ListIter j = i;
      j++;
      groups.erase(i);
      i = j;
    }
    else
      i++;
  }
}

/// Writes an constraint to the specified filename in VRML format for visualization
/**
 * \todo add a cone onto the arrows
 */
void UnilateralConstraint::write_vrml(const std::string& fname, double sphere_radius, double normal_length) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  std::ofstream out;
  
  // open the file for writing
  out.open(fname.c_str());
  if (out.fail())
    throw std::runtime_error("Unable to open file for writing in UnilateralConstraint::write_vrml()");

  // write the VRML header
  out << "#VRML V2.0 utf8" << std::endl << std::endl;

  // *************************************************
  // first, write the contact point 
  // *************************************************

  // determine a random color that will be used for contact and normal
  double c_x = (double) rand() / RAND_MAX;
  double c_y = (double) rand() / RAND_MAX;
  double c_z = (double) rand() / RAND_MAX;

  // write the transform for the contact point
  out << "Transform {" << std::endl;
  out << "  translation "; 
  out << contact_point[X] << " " << contact_point[Y] << " " << contact_point[Z] << std::endl;
  out << "  children " << endl;

  // write the shape node, using default appearance
  out << "  Shape {" << std::endl;
  out << "    appearance Appearance { material Material {" << std::endl;
  out << "      transparency 0" << std::endl;
  out << "      shininess 0.2" << std::endl;
  out << "      ambientIntensity 0.2" << std::endl;
  out << "      emissiveColor 0 0 0" << std::endl;
  out << "      specularColor 0 0 0" << std::endl;
  out << "      diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "      }}" << std::endl;

  // write the geometry (a sphere)
  out << "  geometry Sphere {" << std::endl; 
  out << "    radius " << sphere_radius << " }}} # end sphere, shape, transform " << std::endl;

  // *************************************************
  // now, write the normal
  // *************************************************

  // determine the normal edge
  Vector3d normal_start = contact_point;
  Vector3d normal_stop = normal_start + contact_normal*normal_length;

  // write the shape node, using default appearance
  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material {" << std::endl;
  out << "    transparency 0" << std::endl;
  out << "    shininess 0.2" << std::endl;
  out << "    ambientIntensity 0.2" << std::endl;
  out << "    emissiveColor 0 0 0" << std::endl;
  out << "    specularColor 0 0 0" << std::endl;
  out << "    diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "    }}" << std::endl;

  // write the geometry
  out << "  geometry IndexedLineSet {" << std::endl; 
  out << "    coord Coordinate { point [ ";
  out << normal_start[X] << " " << normal_start[Y] << " " << normal_start[Z] << ", ";
  out << normal_stop[X] << " " << normal_stop[Y] << " " << normal_stop[Z] << " ] } " << std::endl;
  out << "    coordIndex [ 0, 1, -1 ] }}" << std::endl;

  // **********************************************
  // determine the axis-angle rotation for the cone
  // **********************************************

  // first compose an arbitrary vector d
  Vector3d d(1,1,1);
  if (std::fabs(contact_normal[X]) > std::fabs(contact_normal[Y]))
  {
    if (std::fabs(contact_normal[X]) > std::fabs(contact_normal[Z]))
      d[X] = 0;
    else
      d[Z] = 0;
  }
  else
  {
    if (std::fabs(contact_normal[Y]) > std::fabs(contact_normal[Z]))
      d[Y] = 0;
    else
      d[Z] = 0;
  }
    
  // compute the cross product of the normal and the vector
  Vector3d x = Vector3d::normalize(Vector3d::cross(contact_normal, d));
  Vector3d y;
  y = contact_normal;
  Vector3d z = Vector3d::normalize(Vector3d::cross(x, contact_normal));

  // compute theta and the axis of rotation
  double theta = std::acos((x[X] + y[Y] + z[Z] - 1)/2);
  Vector3d axis(z[Y] - y[Z], x[Z] - z[X], y[X] - x[Y]);
  axis *= -(1.0/(2 * std::sin(theta)));
    
  // finally, write the cone to show the normal's direction
  out << "Transform {" << std::endl;
  out << "  rotation ";
   out  << axis[X] <<" "<< axis[1] <<" "<< axis[Z] <<" "<< theta << std::endl;
  out << "  translation ";
   out << normal_stop[X] <<" "<< normal_stop[Y] <<" "<< normal_stop[Z];
  out << std::endl;
  out << "  children [" << std::endl;
  out << "    Shape {" << std::endl;
  out << "      appearance Appearance { material Material {" << std::endl;
  out << "        transparency 0" << std::endl;
  out << "        shininess 0.2" << std::endl;
  out << "        ambientIntensity 0.2" << std::endl;
  out << "        emissiveColor 0 0 0" << std::endl;
  out << "        specularColor 0 0 0" << std::endl;
  out << "        diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "        }}" << std::endl;
  out << "      geometry Cone {" << std::endl;
  out << "        bottomRadius " << sphere_radius << std::endl;
  out << "        height " << (normal_length * .1) << std::endl;
  out << "      } } ] }" << std::endl;
  out.close();
}

// computes the velocity of a contact constraint in a particular direction 
double UnilateralConstraint::calc_contact_vel(const Vector3d& v) const
{
  // verify that this is a contact
  assert(constraint_type == eContact);

  shared_ptr<SingleBodyd> sba = contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sbb = contact_geom2->get_single_body();
  assert(sba && sbb);

  // get the vels 
  const SVelocityd& va = sba->get_velocity();
  const SVelocityd& vb = sbb->get_velocity();

  // setup the constraint frame
  _contact_frame->x = contact_point;
  _contact_frame->q.set_identity();
  _contact_frame->rpose = GLOBAL;

  // compute the velocities at the contact point
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(_contact_frame);
  SVelocityd ta = Pose3d::transform(const_contact_frame, va);
  SVelocityd tb = Pose3d::transform(const_contact_frame, vb);

  // transform the vector
  Vector3d vx = Pose3d::transform_vector(_contact_frame, v);

  return vx.dot(ta.get_linear() - tb.get_linear());
}

/// Determines the set of contact tangents
void UnilateralConstraint::determine_contact_tangents()
{
  // get the two bodies of the contact
  assert(constraint_type == UnilateralConstraint::eContact);
  assert(contact_geom1 && contact_geom2);
  shared_ptr<SingleBodyd> sba = contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sbb = contact_geom2->get_single_body();
  assert(sba && sbb);

  // verify the contact point, normal, and tangents are in the global frame
  assert(contact_point.pose == GLOBAL);
  assert(contact_normal.pose == GLOBAL);
  assert(contact_tan1.pose == GLOBAL);
  assert(contact_tan2.pose == GLOBAL);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = contact_point;

  // get the velocities at the point of contat
  const SVelocityd& va = sba->get_velocity(); 
  const SVelocityd& vb = sbb->get_velocity();
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(_contact_frame);
  SVelocityd ta = Pose3d::transform(const_contact_frame, va);
  SVelocityd tb = Pose3d::transform(const_contact_frame, vb);
  Vector3d rvel = ta.get_linear() - tb.get_linear();

  // get the normal in the same frame
  Vector3d normal_cp = Pose3d::transform_vector(_contact_frame, contact_normal);

  // now remove the normal components from this relative velocity
  double dot = normal_cp.dot(rvel);
  rvel -= (normal_cp * dot);

  // see whether we can use this vector as a contact tangent and set the
  // friction type 
  double tan_norm = rvel.norm();
  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_contact_tangents() - tangent velocity magnitude: " << tan_norm << std::endl;

  // determine an orthonormal basis using the two contact tangents
  Vector3d::determine_orthonormal_basis(contact_normal, contact_tan1, contact_tan2);
  assert(!std::isnan(contact_tan1.norm()));
  assert(!std::isnan(contact_tan2.norm()));
}

/// Determines the type of constraint 
UnilateralConstraint::UnilateralConstraintClass UnilateralConstraint::determine_constraint_class() const
{
  // get the constraint velocity
  double vel = calc_constraint_vel();

  FILE_LOG(LOG_SIMULATOR) << "-- constraint type: " << constraint_type << " velocity: " << vel << std::endl;

  if (vel > tol)
    return ePositive;
  else if (vel < -tol)
    return eNegative;
  else
    return eZero;
}

/// Computes the constraint tolerance
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
double UnilateralConstraint::calc_constraint_tol() const
{
  if (constraint_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    shared_ptr<SingleBodyd> sba = contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 

    // setup the constraint frame
    _contact_frame->x = contact_point;
    _contact_frame->q.set_identity();
    _contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
    shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(_contact_frame);
    SVelocityd ta = Pose3d::transform(const_contact_frame, va); 
    SVelocityd tb = Pose3d::transform(const_contact_frame, vb); 

    // compute the difference in linear velocities
    return std::max((ta.get_linear() - tb.get_linear()).norm(), (double) 1.0);
  }
  else if (constraint_type == eLimit)
  {
    double qd = limit_joint->qd[limit_dof];
    return std::max((double) 1.0, std::fabs(qd));
  }
  else
  {
    assert(false);
    return 0.0;
  }
}

/// Gets the super bodies for the constraint
unsigned UnilateralConstraint::get_super_bodies(shared_ptr<DynamicBodyd>& db1, shared_ptr<DynamicBodyd>& db2) const
{
  // look for empty constraint
  if (constraint_type == UnilateralConstraint::eNone)
    return 0;

  // look for limit constraint
  if (constraint_type == UnilateralConstraint::eLimit)
  {
    RigidBodyPtr outboard = limit_joint->get_outboard_link();
    db1 = dynamic_pointer_cast<DynamicBodyd>(outboard->get_articulated_body());
    return 1;
  }
  else if (constraint_type == UnilateralConstraint::eContact)
  {
    shared_ptr<SingleBodyd> sb1 = contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(sb1->get_articulated_body());
    ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(sb2->get_articulated_body());
    if (ab1)
      db1 = ab1;
    else
    {
      if (sb1->is_enabled())
        db1 = dynamic_pointer_cast<DynamicBodyd>(sb1);
    }
    if (ab2)
      db2 = ab2;
    else
    {
      if (sb2->is_enabled())
        db2 = dynamic_pointer_cast<DynamicBodyd>(sb2);
    }
    return 2;
  }
  else
  {
    assert(false);
    return 0;
  }
}

