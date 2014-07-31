/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public
 * License (found in COPYING).
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
#include <Moby/Event.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/PenaltyEventHandler.h>
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

/// Sets up the default parameters for the Penalty event handler
PenaltyEventHandler::PenaltyEventHandler(){}

// Processes Penaltys
void PenaltyEventHandler::process_events(const vector<Event>& events) const
{
  FILE_LOG(LOG_EVENT) << "*************************************************************";
  FILE_LOG(LOG_EVENT) << endl;
  FILE_LOG(LOG_EVENT) << "PenaltyEventHandler::process_events() entered";
  FILE_LOG(LOG_EVENT) << endl;
  FILE_LOG(LOG_EVENT) << "*************************************************************";
  FILE_LOG(LOG_EVENT) << endl;

  // apply the method to all contacts
  if (!events.empty())
    apply_model(events);
  else
    FILE_LOG(LOG_EVENT) << " (no events?!)" << endl;

  FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
  FILE_LOG(LOG_EVENT) << "PenaltyEventHandler::process_events() exited" << endl;
  FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
}

/// Applies the model to a set of events
/**
 * \param events a set of events
 */
void PenaltyEventHandler::apply_model(const vector<Event>& events) const
{
  const double INF = std::numeric_limits<double>::max();

  // **********************************************************
  // Apply penalty model to each event
  // **********************************************************
  for (int i=0;i<events.size();i++)
  {
    const Event& e = events[i];
    if(!e.contact_compliant)
      continue;
    // setup the event frame
    shared_ptr<Pose3d> event_frame(new Pose3d);
    event_frame->x = e.contact_point;
    event_frame->q.set_identity();
    event_frame->rpose = GLOBAL;

    // get the contact normal in the correct pose
    Vector3d normal = Pose3d::transform_vector(event_frame, e.contact_normal);
    Vector3d tan1;// = Pose3d::transform_vector(event_frame, e.contact_tan1);
    Vector3d tan2;// = Pose3d::transform_vector(event_frame, e.contact_tan2);
    Ravelin::Vector3d::determine_orthonormal_basis(normal,tan1,tan2);

    FILE_LOG(LOG_EVENT) << "relative normal: " << normal << std::endl;
    FILE_LOG(LOG_EVENT) << "relative tangent 1: " << tan1 << std::endl;
    FILE_LOG(LOG_EVENT) << "relative tangent 2: " << tan2<< std::endl;


    Vector3d penalty_force(0,0,0,event_frame);

    // get contacting bodies
    SingleBodyPtr sba = e.contact_geom1->get_single_body();
    SingleBodyPtr sbb = e.contact_geom2->get_single_body();
    assert(sba && sbb);

    // get the vels
    const SVelocityd& va = sba->get_velocity();
    const SVelocityd& vb = sbb->get_velocity();

    // compute the velocities at the contact point
    SVelocityd ta = Pose3d::transform(event_frame, va);
    SVelocityd tb = Pose3d::transform(event_frame, vb);

    // get the linear velocities and project against the normal
    Vector3d rvlin = ta.get_linear() - tb.get_linear();
//    assert(std::fabs(normal.dot(rvlin)) < NEAR_ZERO || std::fabs(normal.dot(rvlin) - calc_event_vel2(e))/std::fabs(normal.dot(rvlin)) < NEAR_ZERO);
    FILE_LOG(LOG_EVENT) << "relative normal velocity: " << normal.dot(rvlin) << std::endl;
    FILE_LOG(LOG_EVENT) << "relative tangent 1 velocity: " << tan1.dot(rvlin) << std::endl;
    FILE_LOG(LOG_EVENT) << "relative tangent 2 velocity: " << tan2.dot(rvlin) << std::endl;

    /// Position -----------------------------------------------------------------


    Point3d cpa = Pose3d::transform_point(sba->get_pose(),e.contact_point),
            cpb = Pose3d::transform_point(sbb->get_pose(),e.contact_point);
    Point3d cp = e.contact_point;

    FILE_LOG(LOG_EVENT) << "cpa: " << cpa << endl;
    FILE_LOG(LOG_EVENT) << "cpb: " << cpb << endl;

    double depth = CollisionGeometry::calc_signed_dist(e.contact_geom1,
                                                       e.contact_geom2,
                                                       cpa,
                                                       cpb);

    FILE_LOG(LOG_EVENT) << "Contact Point: " << e.contact_point << endl;
    FILE_LOG(LOG_EVENT) << "I.P. Depth: " << depth << endl;
    FILE_LOG(LOG_EVENT) << "Contact Velocity [x y z]: " << rvlin << endl;
    FILE_LOG(LOG_EVENT) << "Contact Velocity [N S T]: " << Vector3d(normal.dot(rvlin), tan1.dot(rvlin), tan2.dot(rvlin)) << endl;
    FILE_LOG(LOG_EVENT) << "contact_depth_penalty (Kp): " << e.contact_penalty_Kp << endl;
    FILE_LOG(LOG_EVENT) << "contact_velocity_penalty (Kv): " << e.contact_penalty_Kv << endl;
    FILE_LOG(LOG_EVENT) << "contact_mu_viscous: " << e.contact_mu_viscous << endl;

    if(depth < 0){
      // Depth
      penalty_force += depth * normal * e.contact_penalty_Kp;

      // Velocity
      penalty_force += -rvlin.dot(normal) * normal * e.contact_penalty_Kv;

      // Friction
      penalty_force += -tan1.dot(rvlin) * tan1 * e.contact_mu_viscous;
      penalty_force += -tan2.dot(rvlin) * tan2 * e.contact_mu_viscous;

      FILE_LOG(LOG_EVENT) << "Penalty Force: " << penalty_force << endl;

      Ravelin::VectorNd gf;
      if(sba->is_enabled()){
        sba->convert_to_generalized_force(sba,Ravelin::SForced(penalty_force,Vector3d(0,0,0,event_frame),event_frame),gf);
        sba->add_generalized_force(gf);
      }

      if(sbb->is_enabled()){
        sbb->convert_to_generalized_force(sbb,Ravelin::SForced(-penalty_force,Vector3d(0,0,0,event_frame),event_frame),gf);
        sbb->add_generalized_force(gf);
      }

    }
  }
}


