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
#include <Moby/SingleBody.h>
#include <Moby/Spatial.h>
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
  stick_tol = NEAR_ZERO;
  compliance = eRigid;
  constraint_type = eNone;
  signed_violation = 0.0;
  limit_dof = std::numeric_limits<unsigned>::max();
  limit_epsilon = (double) 0.0;
  limit_upper = false;
  limit_impulse = (double) 0.0;
  contact_normal.set_zero();
  contact_normal_dot.set_zero();
  contact_tan1_dot.set_zero();
  contact_tan2_dot.set_zero();
  contact_impulse.set_zero();
  contact_point.set_zero();
  contact_mu_coulomb = (double) 0.0;
  contact_mu_viscous = (double) 0.0;
  contact_penalty_Kp = (double) 0.0;
  contact_penalty_Kv = (double) 0.0;
  contact_epsilon = (double) 0.0;
  contact_NK = 4;
  _ftype = eUndetermined;
  deriv_type = eVel;
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
  contact_normal_dot = e.contact_normal_dot;
  contact_tan1_dot = e.contact_tan1;
  contact_tan2_dot = e.contact_tan2;
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
  stick_tol = e.stick_tol;
  _ftype = e._ftype;
  deriv_type = e.deriv_type;

  return *this;
}

/// Computes the constraint data
void UnilateralConstraint::compute_constraint_data(MatrixNd& M, VectorNd& q) const
{
  if (deriv_type == eVel)
    compute_vconstraint_data(M, q);
  else
    compute_aconstraint_data(M, q);
}

/// Computes the cross constraint data
void UnilateralConstraint::compute_cross_constraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  assert(deriv_type == e.deriv_type);

  if (deriv_type == eVel)
    compute_cross_vconstraint_data(e, M);
  else
    compute_cross_aconstraint_data(e, M);
}

/// Computes the acceleration constraint data
void UnilateralConstraint::compute_aconstraint_data(MatrixNd& M, VectorNd& q) const
{
  assert(constraint_type == eContact);

  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the two single bodies
  SingleBodyPtr sb1 = contact_geom1->get_single_body();
  SingleBodyPtr sb2 = contact_geom2->get_single_body();

  // get the two super bodies
  DynamicBodyPtr su1 = sb1->get_super_body();
  DynamicBodyPtr su2 = sb2->get_super_body();

  // get the numbers of generalized coordinates for the two super bodies
  const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBody::eSpatial);

  // verify the contact point, normal, and tangents are in the global frame
  assert(contact_point.pose == GLOBAL);
  assert(contact_normal.pose == GLOBAL);
  assert(contact_tan1.pose == GLOBAL);
  assert(contact_tan2.pose == GLOBAL);

  // verify that the friction type has been set
  assert(_ftype != eUndetermined);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = contact_point;

  // case 1: sticking friction
  if (_ftype == eSticking)
  {
    // get the directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

    // setup a matrix of contact directions
    Matrix3d R;
    R.set_column(N, normal);
    R.set_column(S, tan1);
    R.set_column(T, tan2);

    // resize the Jacobians 
    J1.resize(THREE_D, NGC1);
    J2.resize(THREE_D, NGC2);

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin1, J1);
    su2->calc_jacobian(_contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin2, J2);

    // compute the contact inertia matrix for the first body
    su1->transpose_solve_generalized_inertia(J1, workM1);
    J1.mult(workM1, M);

    // compute the contact inertia matrix for the second body
    su2->transpose_solve_generalized_inertia(J2, workM1);
    J2.mult(workM1, workM2);
    M += workM2;

    // compute the directional accelerations
    su1->get_generalized_acceleration(v);
    J1.mult(v, q);
    su2->get_generalized_acceleration(v);
    q += J2.mult(v, workv);

    // update the contact vector data
    compute_dotv_data(q);
  }
  else
  {
    // get the normal and sliding directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);

    // resize the Jacobians 
    J1.resize(1,NGC1);   // normal Jacobian for body 1
    J2.resize(1,NGC2);   // normal Jacobian for body 2
    dJ1.resize(1,NGC1);  // sliding Jacobian for body 1
    dJ2.resize(1,NGC2);  // sliding Jacobian for body 2

    // get shared vectors
    SharedVectorNd J1n = J1.row(0);  // normal direction Jacobian for body 1 
    SharedVectorNd J1s = dJ1.row(0); // sliding direction Jacobian for body 1
    SharedVectorNd J2n = J2.row(0);  // normal direction Jacobian for body 2
    SharedVectorNd J2s = dJ2.row(0); // sliding direction Jacobian for body 2

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    Jlin1.transpose_mult(normal, J1n);
    Jlin1.transpose_mult(tan1, J1s);
    su2->calc_jacobian(_contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    Jlin2.transpose_mult(-normal, J2n);
    Jlin2.transpose_mult(-tan1, J2s);

    // setup the first solution vector (N - u_s*Q)
    J1s *= -contact_mu_coulomb; 
    J1s += J1n; 

    // compute N*inv(M)*(N - u_s*Q)' for the first body
    su1->solve_generalized_inertia(J1s, workv);
    M.resize(1,1);
    M(0,0) = J1n.dot(workv);

    // setup the second solution vector (N - u_s*Q)
    J2s *= -contact_mu_coulomb; 
    J2s += J2n;

    // compute N*inv(M)*(N - u_s*Q)' for the second body
    su2->solve_generalized_inertia(J2s, workv);
    M(0,0) += J2n.dot(workv);

    // compute the normal acceleration
    su1->get_generalized_acceleration(v);
    q[0] = J1n.dot(v);
    su2->get_generalized_acceleration(v);
    q[0] += J2n.dot(v);

    // update the contact vector data
    compute_dotv_data(q);
  }
} 

/// Computes the contact vector data (\dot{N}v and Na)
void UnilateralConstraint::compute_dotv_data(VectorNd& q) const
{
  assert(constraint_type == eContact);

  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the two single bodies
  SingleBodyPtr sb1 = contact_geom1->get_single_body();
  SingleBodyPtr sb2 = contact_geom2->get_single_body();

  // get the two super bodies
  DynamicBodyPtr su1 = sb1->get_super_body();
  DynamicBodyPtr su2 = sb2->get_super_body();

  // get the numbers of generalized coordinates for the two super bodies
  const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBody::eSpatial);

  // verify the derivative of the direction vectors are in the global frame
  assert(contact_normal_dot.pose == GLOBAL);
  assert(contact_tan1_dot.pose == GLOBAL);
  assert(contact_tan2_dot.pose == GLOBAL);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = contact_point;

  // case 1: sticking friction
  if (_ftype == eSticking)
  {
    // get the directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

    // get the directional derivatives in the constraint frame 
    Vector3d dnormal = Pose3d::transform_vector(_contact_frame, contact_normal_dot);
    Vector3d dtan1 = Pose3d::transform_vector(_contact_frame, contact_tan1_dot);
    Vector3d dtan2 = Pose3d::transform_vector(_contact_frame, contact_tan2_dot);

    // setup a matrices of contact directions and directional derivatives
    Matrix3d R, dR;
    R.set_column(N, normal);
    R.set_column(S, tan1);
    R.set_column(T, tan2);
    dR.set_column(N, dnormal);
    dR.set_column(S, dtan1);
    dR.set_column(T, dtan2);

    // resize the Jacobians 
    J1.resize(THREE_D, NGC1);
    J2.resize(THREE_D, NGC2);
    dJ1.resize(THREE_D, NGC1);
    dJ2.resize(THREE_D, NGC2);

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    dR.transpose_mult(Jlin1, J1);
    su2->calc_jacobian(_contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    (-dR).transpose_mult(Jlin2, J2);

    // compute the time-derivatives of the Jacobians for the two bodies
    su1->calc_jacobian_dot(_contact_frame, sb1, JJ);
    SharedConstMatrixNd dJlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(dJlin1, dJ1); 
    su2->calc_jacobian_dot(_contact_frame, sb2, JJ);
    SharedConstMatrixNd dJlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(dJlin2, dJ2);

    // update J1 and J2
    if (dJ1.columns() > 0) J1 += dJ1;
    if (dJ2.columns() > 0) J2 += dJ2;

    // scale J
    J1 *= 2.0;
    J2 *= 2.0;

    // update v using 2*\dot{J}*[n t1 t2]
    su1->get_generalized_velocity(DynamicBody::eSpatial, v);
    FILE_LOG(LOG_CONSTRAINT) << "Body 1 generalized velocity: " << v << std::endl;
    q += J1.mult(v, workv);
    su2->get_generalized_velocity(DynamicBody::eSpatial, v);
    FILE_LOG(LOG_CONSTRAINT) << "Body 2 generalized velocity: " << v << std::endl;
    q += J2.mult(v, workv);

    FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::compute_dotv_data() exited" << std::endl;
  }
  else
  {
    // get the normal and its time derivative in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d dnormal = Pose3d::transform_vector(_contact_frame, contact_normal_dot);

    // resize the Jacobians 
    J1.resize(1, NGC1);
    J2.resize(1, NGC2);
    dJ1.resize(1, NGC1);
    dJ2.resize(1, NGC2);

    // get the rows of the Jacobians for output
    SharedVectorNd J1n = J1.row(N); 
    SharedVectorNd J2n = J2.row(N); 
    SharedVectorNd dJ1n = dJ1.row(N); 
    SharedVectorNd dJ2n = dJ2.row(N); 

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    Jlin1.transpose_mult(dnormal, J1n);
    su2->calc_jacobian(_contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    Jlin2.transpose_mult(-dnormal, J2n);

    // compute the time-derivatives of the Jacobians for the two bodies
    su1->calc_jacobian_dot(_contact_frame, sb1, JJ);
    SharedConstMatrixNd dJlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    dJlin1.transpose_mult(normal, dJ1n);
    su2->calc_jacobian_dot(_contact_frame, sb2, JJ);
    SharedConstMatrixNd dJlin2 = JJ.block(0, THREE_D, 0, JJ.columns());
    dJlin2.transpose_mult(-normal, dJ2n);

    // update J1 and J2
    J1 += dJ1;
    J2 += dJ2;

    // scale J
    J1 *= 2.0;
    J2 *= 2.0;

    // update v using 2*\dot{J}*[n t1 t2]
    su1->get_generalized_velocity(DynamicBody::eSpatial, v);
    FILE_LOG(LOG_CONSTRAINT) << "Body 1 generalized velocity: " << v << std::endl;
    q += J1.mult(v, workv);
    su2->get_generalized_velocity(DynamicBody::eSpatial, v);
    FILE_LOG(LOG_CONSTRAINT) << "Body 2 generalized velocity: " << v << std::endl;
    q += J2.mult(v, workv);

    FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::compute_dotv_data() exited" << std::endl;
  }
}

/// Computes the constraint data
void UnilateralConstraint::compute_vconstraint_data(MatrixNd& M, VectorNd& q) const
{
  if (constraint_type == eContact)
  {
    // setup useful indices
    const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

    // get the two single bodies
    SingleBodyPtr sb1 = contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr su1 = sb1->get_super_body();
    DynamicBodyPtr su2 = sb2->get_super_body();

    // verify the contact point, normal, and tangents are in the global frame
    assert(contact_point.pose == GLOBAL);
    assert(contact_normal.pose == GLOBAL);
    assert(contact_tan1.pose == GLOBAL);
    assert(contact_tan2.pose == GLOBAL);

    // setup the contact frame
    _contact_frame->q.set_identity();
    _contact_frame->x = contact_point;

    // get the numbers of generalized coordinates for the two super bodies
    const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBody::eSpatial);
    const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBody::eSpatial);

    // resize the Jacobians 
    J1.set_zero(THREE_D, NGC1);
    J2.set_zero(THREE_D, NGC2);

    // get the directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(_contact_frame, contact_tan2);

    // setup a matrix of contact directions
    Matrix3d R;
    R.set_column(N, normal);
    R.set_column(S, tan1);
    R.set_column(T, tan2);

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin1 = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin1, J1);
    su2->calc_jacobian(_contact_frame, sb2, JJ);
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
    su1->get_generalized_velocity(DynamicBody::eSpatial, v);
    J1.mult(v, q);

    // free v1 and allocate v2 and workv
    su2->get_generalized_velocity(DynamicBody::eSpatial, v);
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
      v.set_zero(su->num_generalized_coordinates(DynamicBody::eSpatial));
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
    SingleBodyPtr e1sb1 = e1.contact_geom1->get_single_body();
    SingleBodyPtr e1sb2 = e1.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr e1s1 = e1sb1->get_super_body();
    DynamicBodyPtr e1s2 = e1sb2->get_super_body();

    // examine against other constraint type
    if (e2.constraint_type == eContact)
    {
      // get the two single bodies
      SingleBodyPtr e2sb1 = e2.contact_geom1->get_single_body();
      SingleBodyPtr e2sb2 = e2.contact_geom2->get_single_body();

      // get the two super bodies
      DynamicBodyPtr e2s1 = e2sb1->get_super_body();
      DynamicBodyPtr e2s2 = e2sb2->get_super_body();

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
void UnilateralConstraint::compute_cross_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // verify that the constraints are linked
  if (!is_linked(*this, e))
    return;
    
  switch (constraint_type)
  {
    case eContact:
      switch (e.constraint_type)
      {
        case eContact: compute_cross_contact_contact_vconstraint_data(e, M); break;
        case eLimit:   compute_cross_contact_limit_vconstraint_data(e, M); break;
        case eNone:    M.resize(0,0); break;
      }
      break;

    case eLimit:
      switch (e.constraint_type)
      {
        case eContact: compute_cross_limit_contact_vconstraint_data(e, M); break;
        case eLimit:   compute_cross_limit_limit_vconstraint_data(e, M); break;
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
void UnilateralConstraint::compute_cross_contact_contact_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // get the unique super bodies
  DynamicBodyPtr bodies[4];
  DynamicBodyPtr* end = get_super_bodies(bodies);
  end = e.get_super_bodies(end);
  std::sort(bodies, end);
  end = std::unique(bodies, end);

  // determine how many unique super bodies we have
  const unsigned NSUPER = end - bodies;

  // clear M
  M.set_zero(3,3);

  // if we have exactly two super bodies, process them individually
  if (NSUPER == 1)
    compute_cross_contact_contact_vconstraint_data(e, M, bodies[0]);
  if (NSUPER == 2)
  {
    compute_cross_contact_contact_vconstraint_data(e, M, bodies[0]);
    compute_cross_contact_contact_vconstraint_data(e, M, bodies[1]);
  }
  else if (NSUPER == 3)
  {
    // find the one common super body
    DynamicBodyPtr bodies1[2], bodies2[2], isect[1];
    DynamicBodyPtr* end1 = get_super_bodies(bodies1);
    DynamicBodyPtr* end2 = e.get_super_bodies(bodies2);
    std::sort(bodies1, end1);
    std::sort(bodies2, end2);
    DynamicBodyPtr* isect_end = std::set_intersection(bodies1, end1, bodies2, end2, isect);
    assert(isect_end - isect == 1);
    compute_cross_contact_contact_vconstraint_data(e, M, isect[0]);
  }
  else if (NSUPER == 4)
    assert(false);
}

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M, DynamicBodyPtr su) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the first two single bodies
  SingleBodyPtr sba1 = contact_geom1->get_single_body();
  SingleBodyPtr sba2 = contact_geom2->get_single_body();

  // get the first two super bodies
  DynamicBodyPtr sua1 = sba1->get_super_body();
  DynamicBodyPtr sua2 = sba2->get_super_body();

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBody::eSpatial);

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
    su->calc_jacobian(_contact_frame, sba1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, J);
    compute_cross_contact_contact_vconstraint_data(e, M, su, J);
  }
  if (sua2 == su)
  {
    su->calc_jacobian(_contact_frame, sba2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, J);
    compute_cross_contact_contact_vconstraint_data(e, M, su, J);
  }
} 

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M, DynamicBodyPtr su, const MatrixNd& J) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the second two single bodies
  SingleBodyPtr sbb1 = e.contact_geom1->get_single_body();
  SingleBodyPtr sbb2 = e.contact_geom2->get_single_body();

  // get the second two super bodies
  DynamicBodyPtr sub1 = sbb1->get_super_body();
  DynamicBodyPtr sub2 = sbb2->get_super_body();

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBody::eSpatial);

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
    su->calc_jacobian(_contact_frame, sbb1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, Jx);

    // now update M 
    su->transpose_solve_generalized_inertia(Jx, workM1);
    M += J.mult(workM1, workM2);
  }
  if (sub2 == su)
  {
    su->calc_jacobian(_contact_frame, sbb2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, Jx);

    // now update M
    su->transpose_solve_generalized_inertia(Jx, workM1);
    M += J.mult(workM1, workM2);
  }
}

/// Updates contact/limit cross constraint data
void UnilateralConstraint::compute_cross_contact_limit_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M) const
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
  SingleBodyPtr sb1 = contact_geom1->get_single_body();
  SingleBodyPtr sb2 = contact_geom2->get_single_body();

  // get the two super bodies
  DynamicBodyPtr su1 = sb1->get_super_body();
  DynamicBodyPtr su2 = sb2->get_super_body();

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
  const unsigned NGC1 = su1->num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned NGC2 = su2->num_generalized_coordinates(DynamicBody::eSpatial);

  // see whether limit is equal to su1
  if (su == su1)
  {
    // resize Jacobian
    J1.resize(THREE_D, NGC1);

    // compute the Jacobians for the two bodies
    su1->calc_jacobian(_contact_frame, sb1, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    R.transpose_mult(Jlin, J1);

    // compute the constraint inertia matrix for the first body
    su1->transpose_solve_generalized_inertia(J1, workM1);

    // get the appropriate row of workM
    M = workM1.row(idx);

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
    su2->calc_jacobian(_contact_frame, sb2, JJ);
    SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
    (-R).transpose_mult(Jlin, J1);

    // compute the constraint inertia matrix for the first body
    su2->transpose_solve_generalized_inertia(J1, workM1);

    // get the appropriate row of workM
    M += workM1.row(idx); 

    // determine whether to negate the row
    if (e.limit_upper)
      M.negate();
  }
} 

/// Updates limit/contact cross constraint data
void UnilateralConstraint::compute_cross_limit_contact_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M) const
{
  // compute the cross constraint data
  e.compute_cross_contact_limit_vconstraint_data(*this, workM2);

  // transpose the matrix
  MatrixNd::transpose(workM2, M);
} 

/// Updates limit/limit cross constraint data
void UnilateralConstraint::compute_cross_limit_limit_vconstraint_data(const UnilateralConstraint& e, MatrixNd& M) const
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
    workv.set_zero(su->num_generalized_coordinates(DynamicBody::eSpatial));
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

/// Updates the contact data
void UnilateralConstraint::compute_cross_aconstraint_data(const UnilateralConstraint& c, MatrixNd& M) const
{
  // verify that the contacts are linked
  if (!is_linked(*this, c))
    return;
   
  if (constraint_type == eContact && c.constraint_type == eContact)
    compute_cross_contact_contact_aconstraint_data(c, M);
  else
    M.resize(0,0); 
}

/// Updates contact/contact cross contact data
/**
 * From two contact points, we can have up to three separate super bodies. 
 */
void UnilateralConstraint::compute_cross_contact_contact_aconstraint_data(const UnilateralConstraint& c, MatrixNd& M) const
{
  // get the unique super bodies
  DynamicBodyPtr bodies[4];
  DynamicBodyPtr* end = get_super_bodies(bodies);
  end = c.get_super_bodies(end);
  std::sort(bodies, end);
  end = std::unique(bodies, end);

  // determine how many unique super bodies we have
  const unsigned NSUPER = end - bodies;

  // determine how many rows and columns
  const unsigned ROWS = (_ftype == eSlipping) ? 1 : 3;
  const unsigned COLS = (c._ftype == eSlipping) ? 1 : 3;

  // clear M
  M.set_zero(ROWS,COLS);

  // if we have exactly two super bodies, process them individually
  if (NSUPER == 1)
    compute_cross_contact_contact_aconstraint_data(c, M, bodies[0]);
  if (NSUPER == 2)
  {
    compute_cross_contact_contact_aconstraint_data(c, M, bodies[0]);
    compute_cross_contact_contact_aconstraint_data(c, M, bodies[1]);
  }
  else if (NSUPER == 3)
  {
    // find the one common super body
    DynamicBodyPtr bodies1[2], bodies2[2], isect[1];
    DynamicBodyPtr* end1 = get_super_bodies(bodies1);
    DynamicBodyPtr* end2 = c.get_super_bodies(bodies2);
    std::sort(bodies1, end1);
    std::sort(bodies2, end2);
    DynamicBodyPtr* isect_end = std::set_intersection(bodies1, end1, bodies2, end2, isect);
    assert(isect_end - isect == 1);
    compute_cross_contact_contact_aconstraint_data(c, M, isect[0]);
  }
  else if (NSUPER == 4)
    assert(false);
}

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_aconstraint_data(const UnilateralConstraint& c, MatrixNd& M, DynamicBodyPtr su) const
{
  static MatrixNd J;

  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the first two single bodies
  SingleBodyPtr sba1 = contact_geom1->get_single_body();
  SingleBodyPtr sba2 = contact_geom2->get_single_body();

  // get the first two super bodies
  DynamicBodyPtr sua1 = sba1->get_super_body();
  DynamicBodyPtr sua2 = sba2->get_super_body();

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBody::eSpatial);

  // verify that the Coulomb friction type has been determined
  assert(_ftype != eUndetermined);

  // handle the two types of friction separately
  if (_ftype == eSticking)
  {
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
      su->calc_jacobian(_contact_frame, sba1, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      R.transpose_mult(Jlin, J);
      compute_cross_contact_contact_aconstraint_data(c, M, su, J);
    }
    if (sua2 == su)
    {
      su->calc_jacobian(_contact_frame, sba2, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      (-R).transpose_mult(Jlin, J);
      compute_cross_contact_contact_aconstraint_data(c, M, su, J);
    }
  }
  else // sliding contact
  {
    // resize Jacobian 
    J.resize(1, NGC);

    // get the row of the Jacobian for output
    SharedVectorNd Jn = J.row(N); 

    // setup the contact frame
    _contact_frame->q.set_identity();
    _contact_frame->x = contact_point;

    // get the normal in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);

    // compute the Jacobians, checking to see whether necessary
    if (sua1 == su)
    {
      su->calc_jacobian(_contact_frame, sba1, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      Jlin.transpose_mult(normal, Jn);
      compute_cross_contact_contact_aconstraint_data(c, M, su, J);
    }
    if (sua2 == su)
    {
      su->calc_jacobian(_contact_frame, sba2, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      Jlin.transpose_mult(-normal, Jn);
      compute_cross_contact_contact_aconstraint_data(c, M, su, J);
    }
  }
} 

/// Computes cross contact data for one super body
void UnilateralConstraint::compute_cross_contact_contact_aconstraint_data(const UnilateralConstraint& c, MatrixNd& M, DynamicBodyPtr su, const MatrixNd& J) const
{
  // setup useful indices
  const unsigned N = 0, S = 1, T = 2, THREE_D = 3;

  // get the second two single bodies
  SingleBodyPtr sbb1 = c.contact_geom1->get_single_body();
  SingleBodyPtr sbb2 = c.contact_geom2->get_single_body();

  // get the second two super bodies
  DynamicBodyPtr sub1 = sbb1->get_super_body();
  DynamicBodyPtr sub2 = sbb2->get_super_body();

  // get the gc pose for the super body
  shared_ptr<const Pose3d> P = su->get_gc_pose();

  // get the number of generalized coordinates for the super body
  const unsigned NGC = su->num_generalized_coordinates(DynamicBody::eSpatial);

  // verify that the friction type is given
  assert(_ftype != eUndetermined);

  // setup the contact frame
  _contact_frame->q.set_identity();
  _contact_frame->x = c.contact_point;

  if (c._ftype == eSticking)
  {
    // resize Jacobian 
    Jx.resize(THREE_D, NGC);

    // get the directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, c.contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, c.contact_tan1);
    Vector3d tan2 = Pose3d::transform_vector(_contact_frame, c.contact_tan2);

    // setup a matrix of contact directions
    Matrix3d R;
    R.set_column(N, normal);
    R.set_column(S, tan1);
    R.set_column(T, tan2);

    // compute the Jacobians, checking to see whether necessary
    if (sub1 == su)
    {
      // first compute the Jacobian
      su->calc_jacobian(_contact_frame, sbb1, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      R.transpose_mult(Jlin, Jx);

      // now update M 
      su->transpose_solve_generalized_inertia(Jx, workM1);
      M += J.mult(workM1, workM2);
    }
    if (sub2 == su)
    {
      // first compute the Jacobian
      su->calc_jacobian(_contact_frame, sbb2, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      (-R).transpose_mult(Jlin, Jx);

      // now update M
      su->transpose_solve_generalized_inertia(Jx, workM1);
      M += J.mult(workM1, workM2);
    }
  }
  else  // sliding contact
  {
    // resize Jacobians 
    Jx.resize(1, NGC);
    Jy.resize(1, NGC);

    // setup the shared vectors
    SharedVectorNd Jxn = Jx.row(N);
    SharedVectorNd Jyn = Jy.row(N);

    // get the normal and sliding directions in the constraint frame
    Vector3d normal = Pose3d::transform_vector(_contact_frame, c.contact_normal);
    Vector3d tan1 = Pose3d::transform_vector(_contact_frame, c.contact_tan1);

    // compute the Jacobians, checking to see whether necessary
    if (sub1 == su)
    {
      // first compute the Jacobian
      su->calc_jacobian(_contact_frame, sbb1, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      Jlin.transpose_mult(normal, Jxn);
      Jlin.transpose_mult(tan1, Jyn);

      // setup the first solution vector (N - u_s*Q)
      Jy *= -contact_mu_coulomb; 
      Jy += Jx; 

      // now update M 
      su->transpose_solve_generalized_inertia(Jy, workM1);
      M += J.mult(workM1, workM2);
    }
    if (sub2 == su)
    {
      su->calc_jacobian(_contact_frame, sbb2, JJ);
      SharedConstMatrixNd Jlin = JJ.block(0, THREE_D, 0, JJ.columns());
      Jlin.transpose_mult(-normal, Jxn);
      Jlin.transpose_mult(-tan1, Jyn);

      // setup the first solution vector (N - u_s*Q)
      Jy *= -contact_mu_coulomb; 
      Jy += Jx; 

      // now update M
      su->transpose_solve_generalized_inertia(Jx, workM1);
      M += J.mult(workM1, workM2);
    }
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
  assert(contact_NK >= 4);
}

double calc_constraint_accel2(const UnilateralConstraint& e)
{
  assert (e.constraint_type == UnilateralConstraint::eContact);
  SingleBodyPtr sba = e.contact_geom1->get_single_body();
  SingleBodyPtr sbb = e.contact_geom2->get_single_body();

  // get the vels and accelerations
  const SVelocityd& va = sba->get_velocity(); 
  const SVelocityd& vb = sbb->get_velocity(); 
  const SAcceld& aa = sba->get_accel(); 
  const SAcceld& ab = sbb->get_accel(); 

  // get the bodies as rigid bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(sba);
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(sbb);

  // transform velocity and acceleration to mixed frames
  shared_ptr<const Pose3d> Pa = rba->get_mixed_pose(); 
  shared_ptr<const Pose3d> Pb = rbb->get_mixed_pose(); 
  SVelocityd tva = Pose3d::transform(Pa, va);
  SVelocityd tvb = Pose3d::transform(Pb, vb);
  SAcceld taa = Pose3d::transform(Pa, aa);
  SAcceld tab = Pose3d::transform(Pb, ab);

  // transform normal and derivative to mixed frame
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = e.contact_point;
  P->rpose = GLOBAL;
  Vector3d normal = Pose3d::transform_vector(P, e.contact_normal);
  Vector3d normal_dot = Pose3d::transform_vector(P, e.contact_normal_dot);

  // setup terms
  Vector3d ra(e.contact_point - Pa->x);
  Vector3d rb(e.contact_point - Pb->x);
  Vector3d xda = tva.get_linear(); 
  Vector3d xdb = tvb.get_linear(); 
  Vector3d xdda = taa.get_linear(); 
  Vector3d xddb = tab.get_linear(); 
  Vector3d wa = tva.get_angular();
  Vector3d wb = tvb.get_angular();
  Vector3d ala = taa.get_angular();
  Vector3d alb = tab.get_angular();
  ra.pose = normal.pose;
  rb.pose = normal.pose;
  xda.pose = normal.pose;
  xdb.pose = normal.pose;
  wa.pose = normal.pose;
  wb.pose = normal.pose;
  xdda.pose = normal.pose;
  xddb.pose = normal.pose;
  ala.pose = normal.pose;
  alb.pose = normal.pose;
  Vector3d v1(xdda - xddb + Vector3d::cross(ala, ra) - Vector3d::cross(alb, rb) + Vector3d::cross(wa, Vector3d::cross(wa, ra)) - Vector3d::cross(wb, Vector3d::cross(wb, rb)));
  Vector3d v2(xda - xdb + Vector3d::cross(wa, ra) - Vector3d::cross(wb, rb));
  v1.pose = normal.pose;
  v2.pose = normal.pose;

  // get the linear velocities and project against the normal
  return normal.dot(v1) + 2.0*normal_dot.dot(v2);
}


/// Computes the acceleration of this contact
/**
 * Positive acceleration indicates acceleration away, negative acceleration
 * indicates acceleration that will lead to impact/interpenetration.
 */
double UnilateralConstraint::calc_constraint_accel() const
{
  if (constraint_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    SingleBodyPtr sba = contact_geom1->get_single_body();
    SingleBodyPtr sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the velocities and accelerations 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 
    const SAcceld& aa = sba->get_accel(); 
    const SAcceld& ab = sbb->get_accel(); 

    // setup the constraint frame
    _contact_frame->x = contact_point;
    _contact_frame->q.set_identity();
    _contact_frame->rpose = GLOBAL;

    // compute the velocities and accelerations at the contact point
    SVelocityd tva = Pose3d::transform(_contact_frame, va); 
    SVelocityd tvb = Pose3d::transform(_contact_frame, vb); 
    SAcceld taa = Pose3d::transform(_contact_frame, aa); 
    SAcceld tab = Pose3d::transform(_contact_frame, ab); 

    // get the contact normal and derivative in the correct pose
    Vector3d normal = Pose3d::transform_vector(_contact_frame, contact_normal);
    Vector3d normal_dot = Pose3d::transform_vector(_contact_frame, contact_normal_dot);

    // compute: d<n, dx/dt + w x r>/dt =  
    // compute: <n,  d^2x/dt^2 + dw/dt x r> + <dn/dt, dx/dt + w x r>  
    double ddot = normal.dot(taa.get_linear() - tab.get_linear());
    ddot += 2.0*normal_dot.dot(tva.get_linear() - tvb.get_linear());

    #ifndef NDEBUG
    static bool displayed_once = false;
    if (!displayed_once && !CompGeom::rel_equal(ddot, calc_constraint_accel2(*this), 1e-4))
    {
      displayed_once = true;
      std::cerr << "UnilateralConstraint::calc_constraint_accel() warning: accelerations do not match to desired tolerance" << std::endl;
      std::cerr << " -- computed acceleration: " << ddot << std::endl;
      std::cerr << " -- checked acceleration: " << calc_constraint_accel2(*this) << std::endl;
    }
    #endif
    return ddot;
  }
  else if (constraint_type == eLimit)
  {
    double qdd = limit_joint->qdd[limit_dof];
    return (limit_upper) ? -qdd : qdd;
  }
  else
  {
    assert(false);
    return 0.0;
  }
}  

double calc_constraint_vel2(const UnilateralConstraint& e)
{
  assert (e.constraint_type == UnilateralConstraint::eContact);
  SingleBodyPtr sba = e.contact_geom1->get_single_body();
  SingleBodyPtr sbb = e.contact_geom2->get_single_body();

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
    SingleBodyPtr sba = contact_geom1->get_single_body();
    SingleBodyPtr sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 

    // setup the constraint frame
    _contact_frame->x = contact_point;
    _contact_frame->q.set_identity();
    _contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
    SVelocityd ta = Pose3d::transform(_contact_frame, va); 
    SVelocityd tb = Pose3d::transform(_contact_frame, vb); 

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
      SingleBodyPtr sb1(e.contact_geom1->get_single_body());
      if (sb1)
      {
        o << "body1: " << sb1->id << std::endl;
      }
      else
        o << "body1: (undefined)" << std::endl;
    }
    else
      o << "geom1: (undefined)" << std::endl;
  
    if (e.contact_geom2)
    {
      SingleBodyPtr sb2(e.contact_geom2->get_single_body());
      if (sb2)
      {
        o << "body2: " << sb2->id << std::endl;
      }    
      else
        o << "body2: (undefined)" << std::endl;
    }
    else
      o << "geom2: (undefined)" << std::endl;

    o << "contact point / normal pose: " << ((e.contact_point.pose) ? Pose3d(*e.contact_point.pose).update_relative_pose(GLOBAL) : GLOBAL) << std::endl;
    o << "contact point: " << e.contact_point << " frame: " << std::endl;
    o << "normal: " << e.contact_normal << " frame: " << std::endl;
    if (e.deriv_type == UnilateralConstraint::eVel)
    {
      SingleBodyPtr sba = e.contact_geom1->get_single_body();
      SingleBodyPtr sbb = e.contact_geom2->get_single_body();
      assert(sba && sbb);

      // get the vels 
      const SVelocityd& va = sba->get_velocity(); 
      const SVelocityd& vb = sbb->get_velocity(); 

      // setup the constraint frame
      shared_ptr<Pose3d> constraint_frame(new Pose3d);
      constraint_frame->x = e.contact_point;
      constraint_frame->q.set_identity();
      constraint_frame->rpose = GLOBAL;

      // compute the velocities at the contact point
      SVelocityd ta = Pose3d::transform(constraint_frame, va); 
      SVelocityd tb = Pose3d::transform(constraint_frame, vb); 

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
    else
      o << "relative normal acceleration: " << e.calc_constraint_accel() << " (1)   " << calc_constraint_accel2(e) << " (2)" << std::endl;
//      o << "relative tangent acceleration: " << (e.tan1.dot(ralin) +  << " (1)   " << calc_constraint_accel2(e) << " (2)" << std::endl;
//      o << "relative tangent acceleration: " << e.calc_constraint_accel() << " (1)   " << calc_constraint_accel2(e) << " (2)" << std::endl;
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
void UnilateralConstraint::determine_connected_constraints(const vector<UnilateralConstraint>& constraints, list<list<UnilateralConstraint*> >& groups)
{
  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_connected_contacts() entered" << std::endl;

  // clear the groups
  groups.clear();

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
  set<SingleBodyPtr> nodes;
  multimap<SingleBodyPtr, SingleBodyPtr> edges;
  typedef multimap<SingleBodyPtr, SingleBodyPtr>::const_iterator EdgeIter;

  // get all single bodies present in the constraints
  for (list<UnilateralConstraint*>::const_iterator i = constraints_copy.begin(); i != constraints_copy.end(); i++)
  {
    if ((*i)->constraint_type == UnilateralConstraint::eContact)
    {
      SingleBodyPtr sb1((*i)->contact_geom1->get_single_body());
      SingleBodyPtr sb2((*i)->contact_geom2->get_single_body());
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
      RigidBodyPtr inboard = (*i)->limit_joint->get_inboard_link();
      RigidBodyPtr outboard = (*i)->limit_joint->get_outboard_link();
      nodes.insert(inboard);
      nodes.insert(outboard);
    }
    else 
      assert(false);
  }

  FILE_LOG(LOG_CONSTRAINT) << " -- single bodies in constraints:" << std::endl;
  if (LOGGING(LOG_CONSTRAINT))
    for (set<SingleBodyPtr>::const_iterator i = nodes.begin(); i != nodes.end(); i++)
      FILE_LOG(LOG_CONSTRAINT) << "    " << (*i)->id << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << std::endl;

  // add connections between articulated rigid bodies -- NOTE: don't process
  // articulated bodies twice!
  set<ArticulatedBodyPtr> ab_processed;
  BOOST_FOREACH(SingleBodyPtr sb, nodes)
  {
    // if the body is not part of an articulated body, skip it
    ArticulatedBodyPtr abody = sb->get_articulated_body();
    if (!abody)
      continue;

    // see whether it has already been processed
    if (ab_processed.find(abody) != ab_processed.end())
      continue;

    // indicate that the articulated body will now have been processed
    ab_processed.insert(abody);

    // get all links in the articulated body
    const vector<RigidBodyPtr>& links = abody->get_links();

    // add edges between all pairs for which there are links
    vector<RigidBodyPtr>::const_iterator j, k;
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
    SingleBodyPtr node = *nodes.begin();

    // add a list to the contact groups
    groups.push_back(list<UnilateralConstraint*>());
    FILE_LOG(LOG_CONSTRAINT) << " -- constraints in group: " << std::endl;

    // create a node queue, with this node added
    std::queue<SingleBodyPtr> node_q;
    node_q.push(node);

    // loop until the queue is empty
    while (!node_q.empty())
    {
      // get the node off of the front of the node queue
      node = node_q.front();
      node_q.pop();

      // erase the node from the set of nodes
      nodes.erase(node);

      // add all neighbors of the node that have not been processed already 
      // to the node queue
      std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
      for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
        if (nodes.find(i->second) != nodes.end())
          node_q.push(i->second);

      // loop through all remaining constraints
      for (list<UnilateralConstraint*>::iterator i = constraints_copy.begin(); i != constraints_copy.end(); )
      {
        if ((*i)->constraint_type == UnilateralConstraint::eContact)
        {
          SingleBodyPtr sb1((*i)->contact_geom1->get_single_body());
          SingleBodyPtr sb2((*i)->contact_geom2->get_single_body());

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
    }
  }

  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_connected_constraints() exited" << std::endl;
}

/*
/// Computes normal and contact Jacobians for a body
void UnilateralConstraint::compute_contact_jacobian(const UnilateralConstraint& e, MatrixN& Jc, MatrixN& iM_JcT, MatrixN& iM_DcT, unsigned ci, const map<DynamicBodyPtr, unsigned>& gc_indices)
{
  map<DynamicBodyPtr, unsigned>::const_iterator miter;
  SAFESTATIC FastThreadable<VectorN> tmpv, tmpv2, workv, workv2;

  // get the two bodies
  SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
  SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

  // get the super bodies
  DynamicBodyPtr ab1 = sb1->get_articulated_body();
  DynamicBodyPtr ab2 = sb2->get_articulated_body();
  DynamicBodyPtr super1 = (ab1) ? ab1 : sb1;
  DynamicBodyPtr super2 = (ab2) ? ab2 : sb2;

  // process the first body
  miter = gc_indices.find(super1);
  if (miter != gc_indices.end())
  {
    const unsigned index = miter->second;

   // convert the normal force to generalized forces
    super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_normal, ZEROS_3, tmpv());
    Jc.set_sub_mat(ci, index, tmpv(), true);

    // convert the tangent forces to generalized forces
    super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_tan1, ZEROS_3, workv());
    super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_tan2, ZEROS_3, workv2());

    // compute iM_JcT and iM_DcT components
    super1->solve_generalized_inertia(DynamicBody::eSpatial, tmpv(), tmpv2());
    iM_JcT.set_sub_mat(index, ci, tmpv2());
    super1->solve_generalized_inertia(DynamicBody::eSpatial, workv(), tmpv2());
    iM_DcT.set_sub_mat(index, ci*2, tmpv2());
    super1->solve_generalized_inertia(DynamicBody::eSpatial, workv2(), tmpv2());
    iM_DcT.set_sub_mat(index, ci*2+1, tmpv2());
  }

  // process the second body
  miter = gc_indices.find(super2);
  if (miter != gc_indices.end())
  {
    const unsigned index = miter->second;

    // convert the normal force to generalized forces
    super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_normal, ZEROS_3, tmpv());
    Jc.set_sub_mat(ci, index, tmpv(), true);

    // compute iM_JcT components
    super2->solve_generalized_inertia(DynamicBody::eSpatial, tmpv(), tmpv2());
    iM_JcT.set_sub_mat(index, ci, tmpv2());

    // convert the tangent forces to generalized forces
    super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_tan1, ZEROS_3, workv());
    super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_tan2, ZEROS_3, workv2());

    // compute iM_JcT and iM_DcT components
    super2->solve_generalized_inertia(DynamicBody::eSpatial, tmpv(), tmpv2());
    iM_JcT.set_sub_mat(index, ci, tmpv2());
    super2->solve_generalized_inertia(DynamicBody::eSpatial, workv(), tmpv2());
    iM_DcT.set_sub_mat(index, ci*2, tmpv2());
    super2->solve_generalized_inertia(DynamicBody::eSpatial, workv2(), tmpv2());
    iM_DcT.set_sub_mat(index, ci*2+1, tmpv2());
  }
}

/// Computes normal and contact Jacobians for a body
void UnilateralConstraint::compute_contact_jacobians(const UnilateralConstraint& e, VectorN& Nc, VectorN& Dcs, VectorN& Dct)
{
  SAFESTATIC FastThreadable<VectorN> Nc1, Nc2, Dcs1, Dcs2, Dct1, Dct2;

  // get the two bodies
  SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
  SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

  // make sure that the two bodies are ordered
  if (sb2 < sb1)
    std::swap(sb1, sb2);

  // get the super bodies
  DynamicBodyPtr ab1 = sb1->get_articulated_body();
  DynamicBodyPtr ab2 = sb2->get_articulated_body();
  DynamicBodyPtr super1 = (ab1) ? ab1 : sb1;
  DynamicBodyPtr super2 = (ab2) ? ab2 : sb2;

  // get the total number of GC's
  const unsigned GC1 = super1->num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned GC2 = super2->num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned NGC = (super1 != super2) ? GC1 + GC2 : GC1;

  // zero the Jacobian vectors
  Nc.set_zero(NGC);
  Dcs.set_zero(NGC);
  Dct.set_zero(NGC);

  // process the first body
  // compute the 'r' vector
  Vector3d r1 = e.contact_point - sb1->get_position();

  // convert the normal force to generalized forces
  super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_normal, ZEROS_3, Nc1());

  // convert first tangent direction to generalized forces
  super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_tan1, ZEROS_3, Dcs1());

  // convert second tangent direction to generalized forces
  super1->convert_to_generalized_force(DynamicBody::eSpatial, sb1, e.contact_point, e.contact_tan2, ZEROS_3, Dct1());

  // convert the normal force to generalized forces
  super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_normal, ZEROS_3, Nc2());

  // convert first tangent direction to generalized forces
  super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_tan1, ZEROS_3, Dcs2());

  // convert second tangent direction to generalized forces
  super2->convert_to_generalized_force(DynamicBody::eSpatial, sb2, e.contact_point, -e.contact_tan2, ZEROS_3, Dct2());

  // now, set the proper elements in the Jacobian
  if (super1 == super2)
  {
    Nc1() += Nc2();
    Dcs1() += Dcs2();
    Dct1() += Dct2();
    Nc.copy_from(Nc1());
    Dcs.copy_from(Dcs1());
    Dct.copy_from(Dct1());
  }
  else
  {
    Nc.set_sub_vec(0, Nc1());
    Dcs.set_sub_vec(0, Dcs1());
    Dct.set_sub_vec(0, Dct1());
    Nc.set_sub_vec(GC1, Nc2());
    Dcs.set_sub_vec(GC1, Dcs2());
    Dct.set_sub_vec(GC1, Dct2());
  }
}
*/

/// Uses the convex hull of the contact manifold to reject contact points
void UnilateralConstraint::determine_convex_set(list<UnilateralConstraint*>& group)
{
return;
}

void UnilateralConstraint::process_convex_set_group(list<UnilateralConstraint*>& group)
{
  vector<Point3d*> hull;

  // get all points
  vector<Point3d*> points;
  BOOST_FOREACH(UnilateralConstraint* e, group)
  {
    assert(e->constraint_type == UnilateralConstraint::eContact);
    points.push_back(&e->contact_point);
  }

  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_convex_set() entered" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << " -- initial number of contact points: " << points.size() << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << " coefficients of friction: ";
  BOOST_FOREACH(const UnilateralConstraint* e, group)
    FILE_LOG(LOG_CONSTRAINT) << e->contact_mu_coulomb << " ";
  FILE_LOG(LOG_CONSTRAINT) << std::endl;

  // determine whether points are collinear
  const Point3d& pA = *points.front(); 
  const Point3d& pZ = *points.back();
  bool collinear = true;
  for (unsigned i=1; i< points.size()-1; i++)
    if (!CompGeom::collinear(pA, pZ, *points[i]))
    {
      collinear = false;
      break;
    }

  // easiest case: collinear
  if (collinear)
  {
    FILE_LOG(LOG_CONSTRAINT) << " -- contact points are all collinear" << std::endl;

    // just get endpoints
    pair<Point3d*, Point3d*> ep;
    CompGeom::determine_seg_endpoints(points.begin(), points.end(), ep);

    // iterate through, looking for the contact points
    for (list<UnilateralConstraint*>::iterator i = group.begin(); i != group.end(); )
    {
      if (&(*i)->contact_point == ep.first || &(*i)->contact_point == ep.second)
        i++;
      else
        i = group.erase(i);
    }
    assert(!group.empty());

    FILE_LOG(LOG_CONSTRAINT) << " -- remaining contact points after removal: " << std::endl;
    if (LOGGING(LOG_CONSTRAINT))
    {
      BOOST_FOREACH(const UnilateralConstraint* e, group)
        FILE_LOG(LOG_CONSTRAINT) << *e << std::endl;
    }

    return;
  }
  // determine whether the contact manifold is 2D or 3D
  else if (is_contact_manifold_2D(group))
  { 
    FILE_LOG(LOG_CONSTRAINT) << " -- contact points appear to be on a 2D contact manifold" << std::endl;

    try
    {
      // attempt to fit a plane to the points
      Vector3d normal;
      double offset;
      CompGeom::fit_plane(points.begin(), points.end(), normal, offset);

      // compute the 2D convex hull
      CompGeom::calc_convex_hull(points.begin(), points.end(), normal, std::back_inserter(hull));
      if (hull.empty())
        throw NumericalException();
    }
    catch (NumericalException e)
    {
      FILE_LOG(LOG_CONSTRAINT) << " -- unable to compute 2D convex hull; falling back to computing line endpoints" << std::endl;

      // compute the segment endpoints
      pair<Point3d*, Point3d*> ep;
      CompGeom::determine_seg_endpoints(points.begin(), points.end(), ep);

      // iterate through, looking for the contact points
      for (list<UnilateralConstraint*>::iterator i = group.begin(); i != group.end(); )
      {
        if (&(*i)->contact_point == ep.first || &(*i)->contact_point == ep.second)
          i++;
        else
          i = group.erase(i);
      }

      FILE_LOG(LOG_CONSTRAINT) << " -- remaining contact points after removal: " << std::endl;
      if (LOGGING(LOG_CONSTRAINT))
      {
        BOOST_FOREACH(const UnilateralConstraint* e, group)
          FILE_LOG(LOG_CONSTRAINT) << *e << std::endl;
      }

      return;
    }
  }
  else
  {
    try
    {
      FILE_LOG(LOG_CONSTRAINT) << " -- contact points appear to be on a 3D contact manifold" << std::endl;

      // compute the 3D convex hull
      CompGeom::calc_convex_hull(points.begin(), points.end(), std::back_inserter(hull));
      if (hull.empty())
        throw NumericalException();
    }
    catch (NumericalException e)
    {
      try
      {
        FILE_LOG(LOG_CONSTRAINT) << " -- 3D convex hull failed; trying 2D convex hull" << std::endl;

        // attempt to fit a plane to the points
        Vector3d normal;
        double offset;
        CompGeom::fit_plane(points.begin(), points.end(), normal, offset);

        // compute the 2D convex hull
        CompGeom::calc_convex_hull(points.begin(), points.end(), normal, std::back_inserter(hull));
        if (hull.empty())
          throw NumericalException();

/*
        // hull was successful, fix normals as necessary
        BOOST_FOREACH(UnilateralConstraint* e, group)
        {
          double dot = e->contact_normal.dot(normal);
          if (std::fabs(dot - 1.0) > NEAR_ZERO && std::fabs(dot + 1.0) > NEAR_ZERO)
          {
            // need to fix the normal
          }
        }
*/
      }
      catch (NumericalException e)
      {
        // compute the segment endpoints
        pair<Point3d*, Point3d*> ep;
        CompGeom::determine_seg_endpoints(points.begin(), points.end(), ep);

        // iterate through, looking for the contact points
        for (list<UnilateralConstraint*>::iterator i = group.begin(); i != group.end(); )
        {
          if (&(*i)->contact_point == ep.first || &(*i)->contact_point == ep.second)
            i++;
          else
            i = group.erase(i);
        }

        FILE_LOG(LOG_CONSTRAINT) << " -- unable to compute 2D convex hull; falling back to computing line endpoints" << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << " -- remaining contact points after removal: " << std::endl;
        if (LOGGING(LOG_CONSTRAINT))
        {
          BOOST_FOREACH(const UnilateralConstraint* e, group)
            FILE_LOG(LOG_CONSTRAINT) << *e << std::endl;
        }

        return;
      }      
    }
  }

  // if we're here, convex hull was successful. now sort the hull
  std::sort(hull.begin(), hull.end());

  // remove points
  for (list<UnilateralConstraint*>::iterator i = group.begin(); i != group.end(); )
  {
    if (std::binary_search(hull.begin(), hull.end(), &((*i)->contact_point)))
      i++;
    else
      i = group.erase(i);
  }

  FILE_LOG(LOG_CONSTRAINT) << " -- remaining contact points after removal using convex hull: " << group.size() << std::endl;
}

/// Determines whether all constraints in a set are 2D or 3D
bool UnilateralConstraint::is_contact_manifold_2D(const list<UnilateralConstraint*>& constraints)
{
  // get the first contact as a plane
  assert(constraints.front()->constraint_type == UnilateralConstraint::eContact);
  Plane plane(constraints.front()->contact_normal, constraints.front()->contact_point);

  // iterate over the remaining contacts
  for (list<UnilateralConstraint*>::const_iterator i = ++(constraints.begin()); i != constraints.end(); i++)
  {
    assert((*i)->constraint_type == UnilateralConstraint::eContact);
    if (!plane.on_plane((*i)->contact_point))
      return false;
  }

  return true;
}

/**
 * Complexity of computing a minimal set:
 * N = # of contacts, NGC = # of generalized coordinates
 * NGC << N
 *
 * Cost of computing J*inv(M)*J', J*v for one contacts: NGC^3
 *                                    for R contacts: NGC^3 + 2*NGC^2*R
 * Cost of Modified Gauss elimination for M contacts (M < NGC), 
        M x NGC matrix: M^2*NGC
 *
 * Overall cost: 2*NGC^2*R (for R > NGC, where many redundant contact points
                            present) + NGC^3
 * therefore generalized coordinates are the limiting factor...
 */
/// Computes a minimal set of contact constraints
void UnilateralConstraint::determine_minimal_set(list<UnilateralConstraint*>& group)
{
  // if there are very few constraints, quit now
  if (group.size() <= 4)
    return;

  FILE_LOG(LOG_CONSTRAINT) << "UnilateralConstraint::determine_minimal_set() entered" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << " -- initial number of constraints: " << group.size() << std::endl;

  // setup a mapping from pairs of single bodies to groups of constraints
  map<sorted_pair<SingleBodyPtr>, list<UnilateralConstraint*> > contact_groups;

  // move all contact constraints into separate groups
  for (list<UnilateralConstraint*>::iterator i = group.begin(); i != group.end(); )
  {
    if ((*i)->constraint_type == UnilateralConstraint::eContact)
    {
      // get the two bodies
      SingleBodyPtr sb1 = (*i)->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = (*i)->contact_geom2->get_single_body();

      // move the contact to the group
      contact_groups[make_sorted_pair(sb1, sb2)].push_back(*i);
      i = group.erase(i);
    }
    else
      i++;
  }

  // process each group independently, then recombine
  for (map<sorted_pair<SingleBodyPtr>, list<UnilateralConstraint*> >::iterator i = contact_groups.begin(); i != contact_groups.end(); i++)
  {
    determine_convex_set(i->second);
    group.insert(group.end(), i->second.begin(), i->second.end()); 
  }
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

  SingleBodyPtr sba = contact_geom1->get_single_body();
  SingleBodyPtr sbb = contact_geom2->get_single_body();
  assert(sba && sbb);

  // get the vels 
  const SVelocityd& va = sba->get_velocity();
  const SVelocityd& vb = sbb->get_velocity();

  // setup the constraint frame
  _contact_frame->x = contact_point;
  _contact_frame->q.set_identity();
  _contact_frame->rpose = GLOBAL;

  // compute the velocities at the contact point
  SVelocityd ta = Pose3d::transform(_contact_frame, va);
  SVelocityd tb = Pose3d::transform(_contact_frame, vb);

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
  SingleBodyPtr sba = contact_geom1->get_single_body();
  SingleBodyPtr sbb = contact_geom2->get_single_body();
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
  SVelocityd ta = Pose3d::transform(_contact_frame, va);
  SVelocityd tb = Pose3d::transform(_contact_frame, vb);
  Vector3d rvel = ta.get_linear() - tb.get_linear();

  // get the normal in the same frame
  Vector3d normal_cp = Pose3d::transform_vector(_contact_frame, contact_normal);

  // now remove the normal components from this relative velocity
  double dot = normal_cp.dot(rvel);
  rvel -= (normal_cp * dot);

  // see whether we can use this vector as a contact tangent and set the
  // friction type 
  double tan_norm = rvel.norm();

  if (tan_norm < stick_tol)
  {
    _ftype = eSticking;

    // determine an orthonormal basis using the two contact tangents
    Vector3d::determine_orthonormal_basis(contact_normal, contact_tan1, contact_tan2);
  }
  else
  {
    _ftype = eSlipping;

    contact_tan1 = rvel / tan_norm;
    contact_tan1.pose = GLOBAL;
    contact_tan2 = Vector3d::cross(contact_normal, contact_tan1);
    contact_tan2.normalize();
  }
}

/// Determines the type of constraint 
UnilateralConstraint::UnilateralConstraintClass UnilateralConstraint::determine_constraint_class() const
{
  if (deriv_type == eVel)
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
  else
  {
    // get the constraint acceleration
    double acc = calc_constraint_accel();

    FILE_LOG(LOG_SIMULATOR) << "-- constraint type: " << constraint_type << " acceleration: " << acc << std::endl;

    if (acc > tol)
      return ePositive;
    else if (acc < -tol)
      return eNegative;
    else
      return eZero;
  }
}

/// Computes the constraint tolerance
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
double UnilateralConstraint::calc_vconstraint_tol() const
{
  if (constraint_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    SingleBodyPtr sba = contact_geom1->get_single_body();
    SingleBodyPtr sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 

    // setup the constraint frame
    _contact_frame->x = contact_point;
    _contact_frame->q.set_identity();
    _contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
    SVelocityd ta = Pose3d::transform(_contact_frame, va); 
    SVelocityd tb = Pose3d::transform(_contact_frame, vb); 

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

/// Computes the constraint tolerance
/**
 * Positive velocity indicates separation, negative acceleration indicates
 * contact must be treated, zero acceleration indicates rest, positive
 * acceleration indicates contact is separating.
 */
double UnilateralConstraint::calc_aconstraint_tol() const
{
  if (constraint_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    SingleBodyPtr sba = contact_geom1->get_single_body();
    SingleBodyPtr sbb = contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the velocities and accelerations 
    const SVelocityd& va = sba->get_velocity(); 
    const SVelocityd& vb = sbb->get_velocity(); 
    const SAcceld& aa = sba->get_accel(); 
    const SAcceld& ab = sbb->get_accel(); 

    // compute the velocities and accelerations at the contact point
    SVelocityd tva = Pose3d::transform(contact_point.pose, va); 
    SVelocityd tvb = Pose3d::transform(contact_point.pose, vb); 
    SAcceld taa = Pose3d::transform(contact_point.pose, aa); 
    SAcceld tab = Pose3d::transform(contact_point.pose, ab); 

    // get the relative velocity and acceleration norms
    double rv_norm = (tva.get_linear() - tvb.get_linear()).norm();
    double ra_norm = (taa.get_linear() - tab.get_linear()).norm();

    // compute the tolerance
    return std::max(std::max(rv_norm, ra_norm*contact_normal_dot.norm()*2.0), 1.0);
  }
  else if (constraint_type == eLimit)
  {
    double qdd = limit_joint->qdd[limit_dof];
    return std::max((double) 1.0, std::fabs(qdd));
  }
  else
  {
    assert(false);
    return 0.0;
  }
}

/// Gets the super bodies for the constraint
unsigned UnilateralConstraint::get_super_bodies(DynamicBodyPtr& db1, DynamicBodyPtr& db2) const
{
  // look for empty constraint
  if (constraint_type == UnilateralConstraint::eNone)
    return 0;

  // look for limit constraint
  if (constraint_type == UnilateralConstraint::eLimit)
  {
    RigidBodyPtr outboard = limit_joint->get_outboard_link();
    db1 = outboard->get_articulated_body();
    return 1;
  }
  else if (constraint_type == UnilateralConstraint::eContact)
  {
    SingleBodyPtr sb1 = contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
    if (ab1)
      db1 = ab1;
    else
    {
      if (sb1->is_enabled())
        db1 = sb1;
    }
    if (ab2)
      db2 = ab2;
    else
    {
      if (sb2->is_enabled())
        db2 = sb2;
    }
    return 2;
  }
  else
  {
    assert(false);
    return 0;
  }
}

