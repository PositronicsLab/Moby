/****************************************************************************
 * Copyright 2014 Samuel Zapolsky 
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/permute.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/PenaltyConstraintHandler.h>
#ifdef HAVE_IPOPT
#include <Moby/NQP_IPOPT.h>
#include <Moby/LCP_IPOPT.h>
#endif

using namespace Ravelin;
using namespace Moby;
using std::list;
using boost::shared_ptr;
using std::vector;
using std::map;
using std::endl;
using std::cerr;
using std::pair;
using std::min_element;
using boost::dynamic_pointer_cast;

/// Sets up the default parameters for the Penalty constraint handler
PenaltyConstraintHandler::PenaltyConstraintHandler(){}

// Processes Penaltys
void PenaltyConstraintHandler::process_constraints(const vector<UnilateralConstraint>& constraints) const
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "PenaltyConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // apply the method to all contacts
  if (!constraints.empty())
    apply_model(constraints);
  else
    FILE_LOG(LOG_CONSTRAINT) << " (no constraints?!)" << endl;

  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "PenaltyConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of constraints
/**
 * \param constraints a set of constraints
 */
void PenaltyConstraintHandler::apply_model(const vector<UnilateralConstraint>& constraints) const
{
  const double INF = std::numeric_limits<double>::max();
  std::map<sorted_pair<RigidBodyPtr>, unsigned> deepest;

  // ******************************************************************
  // Apply penalty model to only deepest constraint between two bodies 
  // ******************************************************************
  for (unsigned i=0; i<constraints.size();i++)
  {
    const UnilateralConstraint& e = constraints[i];
    if(e.compliance != UnilateralConstraint::eCompliant)
      continue;

    // only handle contacts for now
    if (e.constraint_type != UnilateralConstraint::eContact)
      continue;

    // get the two rigid bodies out
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(e.contact_geom1->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(e.contact_geom2->get_single_body());

    // if a body is disabled, clear its pointer
    if (!rba->is_enabled())
      rba = RigidBodyPtr();
    if (!rbb->is_enabled())
      rbb = RigidBodyPtr();

    // see whether the bodies have already been processed
    map<sorted_pair<RigidBodyPtr>, unsigned>::const_iterator deepest_iter = deepest.find(make_sorted_pair(rba, rbb));
    if (deepest_iter == deepest.end())
      deepest[make_sorted_pair(rba, rbb)] = i;
    else
    {
      if (constraints[i].signed_violation < constraints[deepest_iter->second].signed_violation)
      deepest[make_sorted_pair(rba, rbb)] = i;
    }  
  }

  // process all deepest pairs
  for (map<sorted_pair<RigidBodyPtr>, unsigned>::const_iterator i = deepest.begin(); i != deepest.end(); i++)
  {
    // process it
    const UnilateralConstraint& e = constraints[i->second];

    // setup the constraint frame
    shared_ptr<Pose3d> contact_frame(new Pose3d);
    contact_frame->x = e.contact_point;
    contact_frame->q.set_identity();
    contact_frame->rpose = GLOBAL;

    // get the contact normal in the correct pose
    Vector3d normal = Pose3d::transform_vector(contact_frame, e.contact_normal);
    Vector3d tan1;// = Pose3d::transform_vector(contact_frame, e.contact_tan1);
    Vector3d tan2;// = Pose3d::transform_vector(contact_frame, e.contact_tan2);
    Ravelin::Vector3d::determine_orthonormal_basis(normal,tan1,tan2);

    FILE_LOG(LOG_CONSTRAINT) << "relative normal: " << normal << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "relative tangent 1: " << tan1 << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "relative tangent 2: " << tan2<< std::endl;

    Vector3d penalty_force(0,0,0,contact_frame);

    // get contacting bodies
    SingleBodyPtr sba = e.contact_geom1->get_single_body();
    SingleBodyPtr sbb = e.contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels
    const SVelocityd& va = sba->get_velocity();
    const SVelocityd& vb = sbb->get_velocity();

    // compute the velocities at the contact point
    SVelocityd ta = Pose3d::transform(contact_frame, va);
    SVelocityd tb = Pose3d::transform(contact_frame, vb);

    // get the linear velocities and project against the normal
    Vector3d rvlin = ta.get_linear() - tb.get_linear();
//    assert(std::fabs(normal.dot(rvlin)) < NEAR_ZERO || std::fabs(normal.dot(rvlin) - calc_event_vel2(e))/std::fabs(normal.dot(rvlin)) < NEAR_ZERO);
      FILE_LOG(LOG_CONSTRAINT) << "relative normal velocity: " << normal.dot(rvlin) << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << "relative tangent 1 velocity: " << tan1.dot(rvlin) << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << "relative tangent 2 velocity: " << tan2.dot(rvlin) << std::endl;

    /// Position -----------------------------------------------------------------


    Point3d cpa = Pose3d::transform_point(sba->get_pose(),e.contact_point),
            cpb = Pose3d::transform_point(sbb->get_pose(),e.contact_point);
    Point3d cp = e.contact_point;

    FILE_LOG(LOG_CONSTRAINT) << "cpa: " << cpa << endl;
    FILE_LOG(LOG_CONSTRAINT) << "cpb: " << cpb << endl;

    double depth = e.signed_violation;

    FILE_LOG(LOG_CONSTRAINT) << "Contact Point: " << e.contact_point << endl;
    FILE_LOG(LOG_CONSTRAINT) << "I.P. Depth: " << depth << endl;
    FILE_LOG(LOG_CONSTRAINT) << "Contact Velocity [x y z]: " << rvlin << endl;
    FILE_LOG(LOG_CONSTRAINT) << "Contact Velocity [N S T]: " << Vector3d(normal.dot(rvlin), tan1.dot(rvlin), tan2.dot(rvlin)) << endl;
    FILE_LOG(LOG_CONSTRAINT) << "contact_depth_penalty (Kp): " << e.contact_penalty_Kp << endl;
    FILE_LOG(LOG_CONSTRAINT) << "contact_velocity_penalty (Kv): " << e.contact_penalty_Kv << endl;
    FILE_LOG(LOG_CONSTRAINT) << "contact_mu_viscous: " << e.contact_mu_viscous << endl;

    if(depth < 0){
      // Depth
      penalty_force -= depth * normal * e.contact_penalty_Kp;

      // Velocity
      penalty_force += -rvlin.dot(normal) * normal * e.contact_penalty_Kv;

      // Friction
      penalty_force += -tan1.dot(rvlin) * tan1 * e.contact_mu_viscous;
      penalty_force += -tan2.dot(rvlin) * tan2 * e.contact_mu_viscous;

      FILE_LOG(LOG_CONSTRAINT) << "Penalty Force: " << penalty_force << endl;

      Ravelin::VectorNd gf;
      if(sba->is_enabled()){
        sba->convert_to_generalized_force(sba,Ravelin::SForced(penalty_force,Vector3d(0,0,0,contact_frame),contact_frame),gf);
        sba->add_generalized_force(gf);
     }

      if(sbb->is_enabled()){
        sbb->convert_to_generalized_force(sbb,Ravelin::SForced(-penalty_force,Vector3d(0,0,0,contact_frame),contact_frame),gf);
        sbb->add_generalized_force(gf);
      }
    }
  }
}


