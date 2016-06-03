/***************************************************************************
 * Copyright 2016 Evan Drumwright
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
#include <Moby/Constraint.h>

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

/// Above this friction coefficient, there is no slip
double Constraint::NO_SLIP_COEFF = 100.0;

/// Creates an empty constraint 
Constraint::Constraint()
{
  tol = NEAR_ZERO;              // default collision tolerance
  constraint_type = eNone;
  signed_violation = 0.0;
  contact_stiffness = 0.0;
  contact_damping = 0.0;
  limit_stiffness = 0.0;
  limit_damping = 0.0;
  limit_dof = std::numeric_limits<unsigned>::max();
  limit_epsilon = (double) 0.0;
  limit_upper = false;
  limit_impulse = (double) 0.0;
  contact_normal.set_zero(GLOBAL);
  contact_impulse.set_zero(GLOBAL);
  contact_point.set_zero(GLOBAL);
  contact_mu_coulomb = (double) 0.0;
  contact_mu_viscous = (double) 0.0;
  contact_epsilon = (double) 0.0;
  contact_NK = 4;
}

Constraint& Constraint::operator=(const Constraint& e)
{
  tol = e.tol;
  signed_violation = e.signed_violation;
  constraint_type = e.constraint_type;
  contact_stiffness = e.contact_stiffness;
  contact_damping = e.contact_damping;
  limit_stiffness = e.limit_stiffness;
  limit_damping = e.limit_damping;
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
  contact_epsilon = e.contact_epsilon;
  contact_NK = e.contact_NK;
  contact_tan1 = e.contact_tan1;
  contact_tan2 = e.contact_tan2;
  implicit_joint = e.implicit_joint;
  inv_dyn_joint = e.inv_dyn_joint;
  qdot_des = e.qdot_des;
  inv_dyn_impulse = e.inv_dyn_impulse;
  inv_dyn_fmax = e.inv_dyn_fmax;
  spring_damper_joint = e.spring_damper_joint;
  spring_damper_impulse = e.spring_damper_impulse;

  return *this;
}

/// Evaluates a spring 
double Constraint::eval_spring() const
{
  // TODO: implement me
}

/// Calculates the maximum force a spring/damper constraint can apply
double Constraint::calc_spring_damper_fmax() const
{
  // TODO: implement me
}

/// Applies a test impulse to a constraint
void Constraint::apply_test_impulse(unsigned var_index)
{
  // setup the impulse
  VectorNd x(num_variables());
  x.set_zero();
  x[var_index] = 1.0;

  // determine the bodies that the impulses are applied to
  map<shared_ptr<DynamicBodyd>, VectorNd> gj;
  apply_impulses(x, gj);

  // apply the impulses
  for (map<shared_ptr<DynamicBodyd>, VectorNd>::iterator i = gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/// Determines whether slack is allowed in the constraint
bool Constraint::is_constraint_slackable(unsigned i) const
{
  if (constraint_type == eInverseDynamics && inv_dyn_fmax.size() == inv_dyn_joint->num_dof())
  {
    const double INF = std::numeric_limits<double>::max();
    bool all_inf = true;
    for (unsigned i=0; i< inv_dyn_fmax.size(); i++)
      if (inv_dyn_fmax[i] < INF)
      {
        bool all_inf = false;
        break;
      }
    return !all_inf;
  }
  else if (constraint_type == eSpringDamper)
  {
    return true;
  }
  else
    return false;
} 

/// Gets the number of optimization / mathematical programming variables
unsigned Constraint::num_variables() const
{
  switch (constraint_type)
  {
    case eNone:
      return 0;

    case eContact:
      #ifdef USE_AP_MODEL
      return contact_NK+2;
      #else
      return (contact_mu_coulomb > 0.0) ? 5 : 1; 
      #endif

    case eLimit:
      return 1;

    case eInverseDynamics:
      return qdot_des.size();

    case eSpringDamper:
      return 1;

    case eImplicitJoint:
      return implicit_joint->num_constraint_eqns();
  }

  assert(false); 
}

/// Gets the limits on the specified constraint variable
double Constraint::get_lower_variable_limit(unsigned var_index) const
{
  const double INF = std::numeric_limits<double>::max();

  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      if (contact_mu_coulomb >= NO_SLIP_COEFF && var_index > 0)
        return -INF;
      return 0.0; 

    case eLimit:
      return 0.0;

    case eInverseDynamics:
      if (inv_dyn_fmax.size() > 0)
        return -inv_dyn_fmax[var_index]; 
      else
        return -INF;

    case eSpringDamper:
      return -calc_spring_damper_fmax();

    case eImplicitJoint:
      return -INF; 
  }

  assert(false); 
}

/// Gets the limits on the specified constraint variable
double Constraint::get_upper_variable_limit(unsigned var_index) const
{
  const double INF = std::numeric_limits<double>::max();

  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      return INF; 

    case eLimit:
      return INF;

    case eInverseDynamics:
      if (inv_dyn_fmax.size() > 0)
        return inv_dyn_fmax[var_index]; 
      else
        return INF;

    case eSpringDamper:
      return calc_spring_damper_fmax();

    case eImplicitJoint:
      return INF; 
  }

  assert(false); 
}

/// Gets the number of impulsive variables
unsigned Constraint::num_impulsive_variables() const
{
  // always return 'true' unless this is an LCP
  #ifdef USE_AP_MODEL
  if (constraint_type == eContact)
    return num_variables()-1;
  #else
  return num_variables();
  #endif
}

/// Gets the type of the specified constraint equation
Constraint::ConstraintEquationType Constraint::get_constraint_equation_type(unsigned constraint_eqn_index) const
{
  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      #ifdef USE_AP_MODEL
      return eComplementarity;
      #else
      return eInequality; 
      #endif

    case eLimit:
      #ifdef USE_AP_MODEL
      return eComplementarity;
      #else
      return eInequality;
      #endif

    case eInverseDynamics:
    case eSpringDamper:
    case eImplicitJoint:
      return eEquality; 
  }
}

/// Gets the damping applied to this constraint (if any) 
double Constraint::get_damping(unsigned constraint_eqn_index) const
{
  switch (constraint_type)
  {
    case eContact:
    {
      // do not damp if past the compliant layer
      if (constraint_eqn_index > 0)
        return 0.0;
      double clayer_depth = contact_geom1->compliant_layer_depth +
                            contact_geom2->compliant_layer_depth;
      return (clayer_depth + signed_violation > 0.0) ? contact_damping : 0.0; 
    }

    case eLimit:
    {
      if (constraint_eqn_index > 0)
        return 0.0;
      double clayer_depth = limit_joint->compliant_layer_depth; 
      return (clayer_depth + signed_violation > 0.0) ? limit_damping : 0.0; 
    }
   
    default:
      return 0.0;
  }
}

/// Applies restitution to a vector of impulsive variables
void Constraint::apply_restitution(VectorNd& x)
{
  switch (constraint_type)
  {
    case eContact:
      x[0] *= contact_epsilon;
      x.segment(1, x.size()).set_zero();
      break;

    case eLimit:
      x[0] *= limit_epsilon;
      x.segment(1, x.size()).set_zero();
      break;

    default:
      x.set_zero();
      break;
  }
}

/// Applies impulses using constraint variables, resulting in generalized impulse vectors to be applied to dynamic bodies 
void Constraint::apply_impulses(const VectorNd& x, map<shared_ptr<DynamicBodyd>, VectorNd>& gj)
{
  map<shared_ptr<DynamicBodyd>, VectorNd>::iterator gj_iter;
  VectorNd workv;

  // get the constraint type
  if (constraint_type == eContact)
  {
    // setup the contact frame
    shared_ptr<Pose3d> P(new Pose3d);
    P->q.set_identity();
    P->x = contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = contact_normal * x[0];
    #ifdef USE_AP_MODEL
    double scal = M_PI * 2.0 / contact_NK;
    for (unsigned i=0; i< contact_NK; i++)
      j += contact_tan 
    #else
    if (contact_mu_coulomb > 0.0)
    {
      j += contact_tan1 * x[1] - contact_tan1 * x[2];
      j += contact_tan2 * x[3] - contact_tan2 * x[4]; 
    }
    #endif

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);
    contact_impulse = Pose3d::transform(GLOBAL, jx);

    // setup the contact impulse 
    SForced w(contact_impulse);

    // get the two single bodies of the contact
    shared_ptr<SingleBodyd> sb1 = contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> sb2 = contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> b1 = dynamic_pointer_cast<DynamicBodyd>(sb1->get_super_body());
    shared_ptr<DynamicBodyd> b2 = dynamic_pointer_cast<DynamicBodyd>(sb2->get_super_body());

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, workv);
      gj_iter->second += workv;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, workv);
      gj_iter->second += workv;
    }
  }
  else if (constraint_type == eLimit)
  {
    // TODO: fix this so that limits need not be part of articulated bodies

    // get the articulated body and the iterator for the articulated body
    shared_ptr<RCArticulatedBody> ab = dynamic_pointer_cast<RCArticulatedBody>(limit_joint->get_articulated_body());
    assert(ab);

    // get the index of the joint
    unsigned idx = limit_joint->get_coord_index() + limit_dof;

    // setup the body if necessary
    gj_iter = gj.find(ab);
    if (gj_iter == gj.end())
    {
      gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));
      gj_iter = gj.find(ab);
    } 

    // set the limit impulse
    limit_impulse = (limit_upper) ? -x[0] : x[0];

    // set the limit force
    gj_iter->second[idx] += limit_impulse;
  }
  else if (constraint_type == eImplicitJoint)
  {
    // get the inboard and outboard links
    shared_ptr<RigidBodyd> inboard = implicit_joint->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = implicit_joint->get_outboard_link();

    // get the two super bodies
    shared_ptr<DynamicBodyd> ib = dynamic_pointer_cast<DynamicBodyd>(inboard->get_super_body());
    shared_ptr<DynamicBodyd> ob = dynamic_pointer_cast<DynamicBodyd>(outboard->get_super_body());

    // setup constraint impulse 
    implicit_joint_impulse = x;

    // each Jacobian is of size joint equations x number of spatial coordinates
    // apply the impulse to the inboard body
    if (inboard->is_enabled())
    {
      // get the super body
      shared_ptr<DynamicBodyd> sb = dynamic_pointer_cast<DynamicBodyd>(inboard->get_super_body());

      // get the Jacobian
      MatrixNd J;
      implicit_joint->calc_constraint_jacobian(true, J);

      // compute the generalized force
      VectorNd gf;
      J.transpose_mult(implicit_joint_impulse, gf);

      // setup the body in the map if necessary
      gj_iter = gj.find(sb);
      if (gj_iter == gj.end())
      {
        gj[sb].set_zero(sb->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(sb);
      } 

      // update the force
      gj_iter->second += gf;
    }

    // apply the impulse to the inboard body
    if (outboard->is_enabled())
    {
      // get the super body
      shared_ptr<DynamicBodyd> sb = dynamic_pointer_cast<DynamicBodyd>(outboard->get_super_body());

      // get the Jacobian
      MatrixNd J;
      implicit_joint->calc_constraint_jacobian(false, J);

      // compute the generalized force
      VectorNd gf;
      J.transpose_mult(implicit_joint_impulse, gf);

      // setup the body in the map if necessary
      gj_iter = gj.find(sb);
      if (gj_iter == gj.end())
      {
        gj[sb].set_zero(sb->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(sb);
      } 

      // update the force
      gj_iter->second += gf;
    } 
  }
  else if (constraint_type == eSpringDamper)
  {
    // get the super body from the inverse dynamics joint (if any)
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(spring_damper_joint->get_articulated_body());
    if (ab)
    {
      // setup a generalized force vector
      VectorNd gf;
      gf.set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));

      // get the spring and damper impulse
      spring_damper_impulse = x;

      // set the impulse appropriately
      const unsigned IDX_START = spring_damper_joint->get_coord_index();
      gf.segment(IDX_START, IDX_START + spring_damper_joint->num_dof()) = x;

      // setup the body in the map if necessary
      gj_iter = gj.find(ab);
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(ab);
      } 

      // update the force
      gj_iter->second += gf;
    }
    else
    {
      // spring damper forces joints not part of articulated bodies not yet
      // supported
      assert(false);
    }
  }
  else if (constraint_type == eInverseDynamics)
  {
    // get the super body from the inverse dynamics joint (if any)
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(inv_dyn_joint->get_articulated_body());
    if (ab)
    {
      // setup a generalized force vector
      VectorNd gf;
      gf.set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));

      // setup the inverse dynamics impulse
      inv_dyn_impulse = x;

      // set the impulse appropriately
      const unsigned IDX_START = inv_dyn_joint->get_coord_index();
      gf.segment(IDX_START, IDX_START + inv_dyn_joint->num_dof()) = x;

      // setup the body in the map if necessary
      gj_iter = gj.find(ab);
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBodyd::eSpatial));
        gj_iter = gj.find(ab);
      } 

      // update the force
      gj_iter->second += gf;
    }
    else
    {
      // inverse dynamics for joints not part of articulated bodies not yet
      // supported
      assert(false);
    }
  }
}

/// Gets the type of constraint coefficient
unsigned Constraint::get_velocity_constraint_index(unsigned constraint_eqn_index) const
{
  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      #ifdef USE_AP_MODEL
      assert(false); 
      #else
      return constraint_eqn_index; 
      #endif

    case eLimit:
    case eInverseDynamics:
    case eSpringDamper:
    case eImplicitJoint:
      return constraint_eqn_index;
  } 
}

/// Gets the type of constraint coefficient
Constraint::ConstraintCoeffType Constraint::get_constraint_coeff_type(unsigned constraint_eqn_index) const
{
  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      #ifdef USE_AP_MODEL
      return eVelocityConstraint;
      #else
      if (contact_mu_coulomb >= NO_SLIP_COEFF)
        return eVelocityConstraint;
      else
        return (constraint_eqn_index == 0) ? eVelocityConstraint : eImpulsiveConstraint; 
      #endif

    case eLimit:
    case eInverseDynamics:
    case eSpringDamper:
    case eImplicitJoint:
      return eVelocityConstraint;
  } 
}

/// Gets the coefficients of a constraint of the form a'*x >= b or a'*x = b
void Constraint::get_impulse_constraint_coeffs(unsigned constraint_eqn_index, SharedVectorNd& a)
{
  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      #ifdef USE_AP_MODEL
      assert(false); 
      #else
      // case 1: no slip 
      assert(contact_mu_coulomb < NO_SLIP_COEFF);
      assert(constraint_eqn_index > 0);
      // case 2: Coulomb friction (linearize the friction cone)
      if (contact_NK == 4)
      {
        a[0] = contact_mu_coulomb;
        a[1] = a[2] = a[3] = a[4] = -1.0;
      }
      else
      {
        double theta = (double) (constraint_eqn_index-1)/(contact_NK/2-1) * M_PI_2;
        a[0] = contact_mu_coulomb;
        a[1] = -std::cos(theta);
        a[2] = -std::sin(theta);
        a[3] = -std::cos(theta);
        a[4] = -std::sin(theta);
      }
      break; 
      #endif

    case eLimit:
    case eInverseDynamics:
    case eSpringDamper:
    case eImplicitJoint:
      assert(false); 
  }
}

/// Gets the *stabilized* projected velocity of the given constraint equation 
double Constraint::calc_projected_stab_vel(unsigned constraint_eqn_index, double inv_dt) const
{
  const unsigned SPATIAL_DIM = 6;

  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      if (constraint_eqn_index == 0)
      {
        double compliant_layer_depth = contact_geom1->compliant_layer_depth +
                                       contact_geom2->compliant_layer_depth;
        double gamma = contact_stiffness * -signed_violation * inv_dt;
        return calc_projected_vel(constraint_eqn_index) - gamma;  // N*v^+ >= gamma 
      }
      else
        return calc_projected_vel(constraint_eqn_index);

    case eLimit:
      assert(constraint_eqn_index == 0);
      {
        double gamma = limit_stiffness * -signed_violation * inv_dt;
        return -calc_projected_vel(constraint_eqn_index) + gamma;  // L*v^+ >= gamma 
      }

    case eInverseDynamics:
      assert(false); // NOTE: this must be tested
      return calc_projected_vel(constraint_eqn_index); // P*v^+ = \dot{q}_des

    case eSpringDamper:
      assert(false); // NOTE: this must be tested
      return eval_spring()*inv_dt; // K*v^+ = \varphi/dt - K*v 

    case eImplicitJoint:
    {
      double C[SPATIAL_DIM];
      implicit_joint->evaluate_constraints(C);
      double gamma = inv_dt*-C[constraint_eqn_index];
      return calc_projected_vel(constraint_eqn_index) - gamma;  // J*v^+ = gamma
    }
  }
}


/// Gets the right-hand side of a constraint of the form a'*x >= b or a'*x = b
double Constraint::get_constraint_rhs(unsigned constraint_eqn_index, double inv_dt)
{
  const unsigned SPATIAL_DIM = 6;

  switch (constraint_type)
  {
    case eNone:
      assert(false); 

    case eContact:
      if (constraint_eqn_index == 0)
      {
        double compliant_layer_depth = contact_geom1->compliant_layer_depth +
                                       contact_geom2->compliant_layer_depth;
        double gamma = contact_stiffness * -signed_violation * inv_dt;
        return -calc_projected_vel(constraint_eqn_index) + gamma;  // N*v^+ >= gamma 
      }
      #ifdef USE_AP_MODEL
      FIX ME
      #else
      // case 1: infinite friction
      if (contact_mu_coulomb >= NO_SLIP_COEFF)
      {
        switch (constraint_eqn_index)
        {
          case 1: return -calc_projected_vel(1);    // S*v^+ = 0
          case 2: return -calc_projected_vel(2);    // T*v^+ = 0
          default:
            assert(false);
        }
      }
      else
      {
        return 0.0;  // these are all linearized friction cone coeffs
      }
      #endif

    case eLimit:
      assert(constraint_eqn_index == 0);
      {
        double gamma = limit_stiffness * -signed_violation * inv_dt;
        return -calc_projected_vel(constraint_eqn_index) + gamma;  // L*v^+ >= gamma 
      }

    case eInverseDynamics:
      return qdot_des[constraint_eqn_index] - calc_projected_vel(constraint_eqn_index); // P*v^+ = \dot{q}_des

    case eSpringDamper:
      return eval_spring()*inv_dt; // K*v^+ = \varphi/dt - K*v 

    case eImplicitJoint:
    {
      double C[SPATIAL_DIM];
      implicit_joint->evaluate_constraints(C);
      double gamma = inv_dt*-C[constraint_eqn_index];
      return gamma-calc_projected_vel(constraint_eqn_index);  // J*v^+ = gamma
    }
  }
}

/// Gets the number of constraint equations
unsigned Constraint::num_constraint_equations() const
{
  switch (constraint_type)
  {
    case eNone:
      return 0;

    case eContact:
      #ifdef USE_AP_MODEL
      return contact_NK + 2;
      #else
      if (contact_mu_coulomb == 0.0)
        return 1;
      else if (contact_mu_coulomb < NO_SLIP_COEFF)
        return contact_NK/4 + 1;
      else
        return 1; 
      #endif

    case eLimit:
      return 1;

    case eInverseDynamics:
      return qdot_des.size();

    case eSpringDamper:
      return 1;

    case eImplicitJoint:
      return implicit_joint->num_constraint_eqns();
  }

  assert(false); 
}

/// Sets the contact parameters for this constraint
void Constraint::set_contact_parameters(const ContactParameters& cparams)
{
  contact_mu_coulomb = cparams.mu_coulomb;
  contact_mu_viscous = cparams.mu_viscous;
  contact_epsilon = cparams.epsilon;
  contact_stiffness = cparams.stiffness;
  contact_damping = cparams.damping;
  contact_NK = cparams.NK;

  // redetermine contact tangents
  determine_contact_tangents();
  
  assert(contact_NK >= 4);
}

/// Computes the velocity of this constraint projected along one of the variable dimensions
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
double Constraint::calc_projected_vel(unsigned var_index) const
{
  shared_ptr<Pose3d> contact_frame(new Pose3d);
  MatrixNd J;
  VectorNd qd, workv;

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
    contact_frame->x = contact_point;
    contact_frame->q.set_identity();
    contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(contact_frame);
    SVelocityd ta = Pose3d::transform(const_contact_frame, va); 
    SVelocityd tb = Pose3d::transform(const_contact_frame, vb); 

    // get the contact direction in the correct pose
    Vector3d dir;
    switch (var_index)
    {
      case 0: 
        dir = Pose3d::transform_vector(contact_frame, contact_normal);
        break;

      case 1:
        dir = Pose3d::transform_vector(contact_frame, contact_tan1);
        break;

      case 2:
        dir = -Pose3d::transform_vector(contact_frame, contact_tan1);
        break;

      case 3:
        dir = Pose3d::transform_vector(contact_frame, contact_tan2);
        break;

      case 4:
        dir = -Pose3d::transform_vector(contact_frame, contact_tan2);
        break;


      default:
        assert(false);
    }

    // get the linear velocities and project against the normal
    return dir.dot(ta.get_linear() - tb.get_linear());
  }
  else if (constraint_type == eLimit)
  {
    assert(var_index == 0);
    double qd = limit_joint->qd[limit_dof];
    return (limit_upper) ? -qd : qd;
  }
  else if (constraint_type == eImplicitJoint)
  {
    // get the inboard and outboard links
    shared_ptr<RigidBodyd> inboard = implicit_joint->get_inboard_link();
    shared_ptr<RigidBodyd> outboard = implicit_joint->get_outboard_link();

    // get the two super bodies
    shared_ptr<DynamicBodyd> ib = dynamic_pointer_cast<DynamicBodyd>(inboard->get_super_body());
    shared_ptr<DynamicBodyd> ob = dynamic_pointer_cast<DynamicBodyd>(outboard->get_super_body());

    // initialize the constraint velocity to zero
    double cvel = 0.0;

    // each Jacobian is of size joint equations x number of spatial coordinates
    // apply the impulse to the inboard body
    if (inboard->is_enabled())
    {
      // get the super body
      shared_ptr<DynamicBodyd> sb = dynamic_pointer_cast<DynamicBodyd>(inboard->get_super_body());

      // get the Jacobian
      implicit_joint->calc_constraint_jacobian(true, J);

      // compute the constraint velocity 
      sb->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
      J.mult(qd, workv);
      cvel += workv[var_index]; 
    }

    // apply the impulse to the inboard body
    if (outboard->is_enabled())
    {
      // get the super body
      shared_ptr<DynamicBodyd> sb = dynamic_pointer_cast<DynamicBodyd>(outboard->get_super_body());

      // get the Jacobian
      implicit_joint->calc_constraint_jacobian(false, J);

      // compute the constraint velocity 
      sb->get_generalized_velocity(DynamicBodyd::eSpatial, qd);
      J.mult(qd, workv);
      cvel += workv[var_index]; 
    } 

    return cvel;
  }
  else if (constraint_type == eInverseDynamics)
  {
    return inv_dyn_joint->qd[var_index];
  }
  else if (constraint_type == eSpringDamper)
  {
    return spring_damper_joint->qd[var_index];
  }
  else
  {
    assert(false);
    return 0.0;
  }
}  


/// Computes the velocity of this constraint
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
double Constraint::calc_constraint_vel(unsigned constraint_eqn_index) const
{
  shared_ptr<Pose3d> contact_frame(new Pose3d);

  if (constraint_type == eContact)
    assert(constraint_eqn_index == 0);
  return calc_projected_vel(constraint_eqn_index);
}  

/// Sends the constraint to the specified stream
std::ostream& Moby::operator<<(std::ostream& o, const Constraint& e)
{
  switch (e.constraint_type)
  {
    case Constraint::eNone:
      o << "(constraint type: none)" << std::endl;
      return o;

    case Constraint::eLimit:
      o << "(constraint type: joint limit)" << std::endl;
      break;

    case Constraint::eContact:
      o << "(constraint type: contact)" << std::endl;
      break;

    case Constraint::eInverseDynamics:
      o << "(constraint type: inverse dynamics)" << std::endl;
      break;

    case Constraint::eImplicitJoint:
      o << "(constraint type: implicit joint)" << std::endl;
      break;

    case Constraint::eSpringDamper:
      o << "(constraint type: spring/damper)" << std::endl;
      break;
  }
	 
  if (e.constraint_type == Constraint::eLimit)
  {
    o << "limit joint ID: " << e.limit_joint->id << std::endl;
    o << "limit joint coordinate index: " << e.limit_joint->get_coord_index() << std::endl;
    o << "limit joint DOF: " << e.limit_dof << std::endl;
    o << "upper limit? " << e.limit_upper << std::endl;
    o << "limit velocity: " << e.calc_constraint_vel(0) << std::endl;
  }
  else if (e.constraint_type == Constraint::eContact)
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

    Pose3d P;
    P.rpose = e.contact_point.pose;
    P.update_relative_pose(GLOBAL);
    o << "contact point / normal pose: " << P << std::endl;
    o << "contact point: " << e.contact_point << std::endl;
    o << "contact normal: " << e.contact_normal << std::endl;
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
  }
  else if (e.constraint_type == Constraint::eSpringDamper)
  {
    o << "spring/damper joint ID: " << e.spring_damper_joint->id << std::endl;
    o << "joint velocity: " << e.calc_constraint_vel(0) << std::endl;
  }
  else if (e.constraint_type == Constraint::eImplicitJoint)
  {
    o << "implicit joint ID: " << e.spring_damper_joint->id << std::endl;
    o << "constraint velocities:";
    for (unsigned i=0; i< e.num_constraint_equations(); i++)
      o << " " << e.calc_constraint_vel(i);
    o << std::endl;
  }
  else if (e.constraint_type == Constraint::eInverseDynamics)
  {
    o << "inverse dynamics joint ID: " << e.inv_dyn_joint->id << std::endl;
    o << "joint velocity: " << e.calc_constraint_vel(0) << std::endl;
  }
  else
    assert(false);

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
osg::Node* Constraint::to_visualization_data() const
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
void Constraint::determine_connected_constraints(const vector<Constraint>& constraints, list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > >& groups, list<vector<shared_ptr<DynamicBodyd> > >& remaining_islands)
{
  FILE_LOG(LOG_CONSTRAINT) << "Constraint::determine_connected_contacts() entered" << std::endl;

  // clear the groups
  groups.clear();
  remaining_islands.clear();

  // copy the list of constraints -- only ones with geometry
  list<Constraint*> constraints_copy;
  BOOST_FOREACH(const Constraint& e, constraints)
    if (e.constraint_type != Constraint::eNone)
      constraints_copy.push_back((Constraint*) &e);
  
  // The way that we'll determine the constraint islands is to treat each rigid
  // body present in the constraints as a node in a graph; nodes will be connected
  // to other nodes if (a) they are both present in constraint or (b) they are
  // part of the same articulated body.  Nodes will not be created for disabled
  // bodies.
  set<shared_ptr<SingleBodyd> > nodes;
  multimap<shared_ptr<SingleBodyd>, shared_ptr<SingleBodyd> > edges;
  typedef multimap<shared_ptr<SingleBodyd>, shared_ptr<SingleBodyd> >::const_iterator EdgeIter;

  // get all single bodies present in the unilateral constraints
  for (list<Constraint*>::const_iterator i = constraints_copy.begin(); i != constraints_copy.end(); i++)
  {
    if ((*i)->constraint_type == Constraint::eContact)
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
    else if ((*i)->constraint_type == Constraint::eLimit)
    {
      shared_ptr<RigidBodyd> inboard = (*i)->limit_joint->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = (*i)->limit_joint->get_outboard_link();
      if (inboard->is_enabled())
        nodes.insert(inboard);
      if (outboard->is_enabled())
        nodes.insert(outboard);
      if (inboard->is_enabled() && outboard->is_enabled())
      {
        edges.insert(std::make_pair(inboard, outboard));
        edges.insert(std::make_pair(outboard, inboard));
      }
    }
    else if ((*i)->constraint_type == Constraint::eInverseDynamics) 
    {
      shared_ptr<RigidBodyd> inboard = (*i)->inv_dyn_joint->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = (*i)->inv_dyn_joint->get_outboard_link();
      if (inboard->is_enabled())
        nodes.insert(inboard);
      if (outboard->is_enabled())
        nodes.insert(outboard);
      if (inboard->is_enabled() && outboard->is_enabled())
      {
        edges.insert(std::make_pair(inboard, outboard));
        edges.insert(std::make_pair(outboard, inboard));
      }
    }
    else if ((*i)->constraint_type == Constraint::eSpringDamper) 
    {
      shared_ptr<RigidBodyd> inboard = (*i)->spring_damper_joint->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = (*i)->spring_damper_joint->get_outboard_link();
      if (inboard->is_enabled())
        nodes.insert(inboard);
      if (outboard->is_enabled())
        nodes.insert(outboard);
      if (inboard->is_enabled() && outboard->is_enabled())
      {
        edges.insert(std::make_pair(inboard, outboard));
        edges.insert(std::make_pair(outboard, inboard));
      }
    }
    else if ((*i)->constraint_type == Constraint::eImplicitJoint)
    {
      shared_ptr<RigidBodyd> inboard = (*i)->implicit_joint->get_inboard_link();
      shared_ptr<RigidBodyd> outboard = (*i)->implicit_joint->get_outboard_link();
      if (inboard->is_enabled())
        nodes.insert(inboard);
      if (outboard->is_enabled())
        nodes.insert(outboard);
      if (inboard->is_enabled() && outboard->is_enabled())
      {
        edges.insert(std::make_pair(inboard, outboard));
        edges.insert(std::make_pair(outboard, inboard));
      }
    }
    else 
      assert(false);
  }

  FILE_LOG(LOG_CONSTRAINT) << " -- single bodies in constraints:" << std::endl;
  if (LOGGING(LOG_CONSTRAINT))
    for (set<shared_ptr<SingleBodyd> >::const_iterator i = nodes.begin(); i != nodes.end(); i++)
      FILE_LOG(LOG_CONSTRAINT) << "    " << (*i)->body_id << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << std::endl;

  // Now, we'll remove nodes from the set until there are no more nodes.
  // For each removed node, we'll get add all constraints that contain the single 
  // body to the group; all neighboring nodes will then be processed.
  while (!nodes.empty())
  {
    // get the node from the front
    shared_ptr<SingleBodyd> node = *nodes.begin();

    // add lists to the contact groups
    groups.push_back(std::pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > >());
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
      nodes.erase(node);
      node_q.pop();

      // add the node to the groups
      groups.back().second.push_back(node);

      // indicate that the node has now been processed
      processed_nodes.insert(node);

      // add all neighbors of the node that have not been processed already 
      // to the node queue
      std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
      for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
        if (processed_nodes.find(i->second) == processed_nodes.end())
          node_q.push(i->second);

      // loop through all remaining constraints
      for (list<Constraint*>::iterator i = constraints_copy.begin(); i != constraints_copy.end(); )
      {
        if ((*i)->constraint_type == Constraint::eContact)
        {
          shared_ptr<SingleBodyd> sb1((*i)->contact_geom1->get_single_body());
          shared_ptr<SingleBodyd> sb2((*i)->contact_geom2->get_single_body());

          // see whether one of the bodies is equal to the node
          if (sb1 == node || sb2 == node)
          {
            assert(!groups.empty());
            groups.back().first.push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->constraint_type == Constraint::eLimit)
        {
          RigidBodyPtr inboard = (*i)->limit_joint->get_inboard_link();
          RigidBodyPtr outboard = (*i)->limit_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            assert(!groups.empty());
            groups.back().first.push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->constraint_type == Constraint::eInverseDynamics) 
        {
          shared_ptr<RigidBodyd> inboard = (*i)->inv_dyn_joint->get_inboard_link();
          shared_ptr<RigidBodyd> outboard = (*i)->inv_dyn_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            assert(!groups.empty());
            groups.back().first.push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->constraint_type == Constraint::eSpringDamper) 
        {
          shared_ptr<RigidBodyd> inboard = (*i)->spring_damper_joint->get_inboard_link();
          shared_ptr<RigidBodyd> outboard = (*i)->spring_damper_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            assert(!groups.empty());
            groups.back().first.push_back(*i);
            i = constraints_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->constraint_type == Constraint::eImplicitJoint)
        {
          shared_ptr<RigidBodyd> inboard = (*i)->implicit_joint->get_inboard_link();
          shared_ptr<RigidBodyd> outboard = (*i)->implicit_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            assert(!groups.empty());
            groups.back().first.push_back(*i);
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

    // if no unilateral constraints have been added, add to remaining islands
    if (groups.back().first.empty())
    {
      // don't need an empty group of unilateral constraints
      groups.pop_back();

      // create a new island
      remaining_islands.push_back(vector<shared_ptr<DynamicBodyd> >());

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
      vector<shared_ptr<DynamicBodyd> >& island = remaining_islands.back();
      std::sort(island.begin(), island.end());
      island.erase(std::unique(island.begin(), island.end()), island.end());
    }
  }

  assert(remaining_islands.empty());

  FILE_LOG(LOG_CONSTRAINT) << "Constraint::determine_connected_constraints() exited" << std::endl;
}

/// Removes groups of contacts that contain no active contacts 
void Constraint::remove_inactive_groups(list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > >& groups)
{
  typedef list<pair<list<Constraint*>, list<shared_ptr<SingleBodyd> > > >::iterator ListIter;

  for (ListIter i = groups.begin(); i != groups.end(); )
  {
    // look for impact in list i
    bool active_detected = false;
    BOOST_FOREACH(Constraint* e, i->first)
    {
      if (e->determine_constraint_velocity_class() == Constraint::eNegative)
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
void Constraint::write_vrml(const std::string& fname, double sphere_radius, double normal_length) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  std::ofstream out;
  
  // open the file for writing
  out.open(fname.c_str());
  if (out.fail())
    throw std::runtime_error("Unable to open file for writing in Constraint::write_vrml()");

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
double Constraint::calc_contact_vel(const Vector3d& v) const
{
  shared_ptr<Pose3d> contact_frame(new Pose3d);

  // verify that this is a contact
  assert(constraint_type == eContact);

  shared_ptr<SingleBodyd> sba = contact_geom1->get_single_body();
  shared_ptr<SingleBodyd> sbb = contact_geom2->get_single_body();
  assert(sba && sbb);

  // get the vels 
  const SVelocityd& va = sba->get_velocity();
  const SVelocityd& vb = sbb->get_velocity();

  // setup the constraint frame
  contact_frame->x = contact_point;
  contact_frame->q.set_identity();
  contact_frame->rpose = GLOBAL;

  // compute the velocities at the contact point
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(contact_frame);
  SVelocityd ta = Pose3d::transform(const_contact_frame, va);
  SVelocityd tb = Pose3d::transform(const_contact_frame, vb);

  // transform the vector
  Vector3d vx = Pose3d::transform_vector(contact_frame, v);

  return vx.dot(ta.get_linear() - tb.get_linear());
}

/// Determines the set of contact tangents
void Constraint::determine_contact_tangents()
{
  shared_ptr<Pose3d> contact_frame(new Pose3d);

  // get the two bodies of the contact
  assert(constraint_type == Constraint::eContact);
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
  contact_frame->q.set_identity();
  contact_frame->x = contact_point;

  // get the velocities at the point of contat
  const SVelocityd& va = sba->get_velocity(); 
  const SVelocityd& vb = sbb->get_velocity();
  shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(contact_frame);
  SVelocityd ta = Pose3d::transform(const_contact_frame, va);
  SVelocityd tb = Pose3d::transform(const_contact_frame, vb);
  Vector3d rvel = ta.get_linear() - tb.get_linear();

  // get the normal in the same frame
  Vector3d normal_cp = Pose3d::transform_vector(contact_frame, contact_normal);

  // now remove the normal components from this relative velocity
  double dot = normal_cp.dot(rvel);
  rvel -= (normal_cp * dot);

  // see whether we can use this vector as a contact tangent and set the
  // friction type 
  double tan_norm = rvel.norm();
  FILE_LOG(LOG_CONSTRAINT) << "Constraint::determine_contact_tangents() - tangent velocity magnitude: " << tan_norm << std::endl;

  // determine an orthonormal basis using the two contact tangents
  Vector3d::determine_orthonormal_basis(contact_normal, contact_tan1, contact_tan2);
  assert(!std::isnan(contact_tan1.norm()));
  assert(!std::isnan(contact_tan2.norm()));
}

/// Determines the type of constraint 
Constraint::ConstraintVelocityClass Constraint::determine_constraint_velocity_class() const
{
  // get the constraint velocity
  double vel = calc_constraint_vel(0);

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
double Constraint::calc_constraint_tol() const
{
  shared_ptr<Pose3d> contact_frame(new Pose3d);

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
    contact_frame->x = contact_point;
    contact_frame->q.set_identity();
    contact_frame->rpose = GLOBAL;

    // compute the velocities at the contact point
    shared_ptr<const Pose3d> const_contact_frame = boost::const_pointer_cast<const Pose3d>(contact_frame);
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

/// Determines whether two constraints are linked
bool Constraint::is_linked(const Constraint& e1, const Constraint& e2)
{
  // get the super bodies from the two constraints and see whether there is
  // an intersection
  std::vector<shared_ptr<DynamicBodyd> > supers1, supers2;
  e1.get_super_bodies(std::back_inserter(supers1));
  e2.get_super_bodies(std::back_inserter(supers2));

  // see whether there is an intersection
  std::list<shared_ptr<DynamicBodyd> > isect;
  std::set_intersection(supers1.begin(), supers1.end(), supers2.begin(), supers2.end(), std::back_inserter(isect));
  return !isect.empty();   
}


