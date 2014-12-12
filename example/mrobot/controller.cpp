/*****************************************************************************
 * Controller for mobile robot 
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
boost::shared_ptr<Moby::EventDrivenSimulator> sim;

using namespace Ravelin;
using namespace Moby;

// setup the controller callback
void controller(DynamicBodyPtr body, double t, void*)
{
  const unsigned LEFT = 0, RIGHT = 1;

  // get the robot
  RCArticulatedBodyPtr robot = boost::dynamic_pointer_cast<RCArticulatedBody>(body);

  // get the state and velocity
  VectorNd q, dq;
  robot->get_generalized_coordinates(DynamicBody::eEuler, q);
  robot->get_generalized_velocity(DynamicBody::eSpatial, dq);

  // setup the PD controller
  const double KP = , KV = ;

  // set q_des, dq_des, ddq_des;

  // compute inverse dynamics torques

  // setup the feedback torques
  VectorNd fleft(1), fright(1);
  fleft[0] = KP*(q_des[LEFT] - q[LEFT]) + KV*(dq_des[LEFT] - dq[LEFT]);
  fright[0] = KP*(q_des[RIGHT] - q[RIGHT]) + KV*(dq_des[RIGHT] - dq[RIGHT]);

  // apply the torques
  JointPtr left = robot->get_joints()[0];
  JointPtr right = robot->get_joints()[1];
  assert(left->id == "left_wheel_hinge");
  assert(right->id == "right_wheel_hinge");
  left->add_force(fleft);
  right->add_force(fright);
}

void post_event_callback_fn(const std::vector<Moby::UnilateralConstraint>& e,
                            boost::shared_ptr<void> empty)
{
  std::cout << ">> start post_event_callback_fn(.)" << std::endl;

  // PROCESS CONTACTS
  for(unsigned i=0;i<e.size();i++){
    if (e[i].constraint_type == Moby::UnilateralConstraint::eContact)
    {
      Moby::SingleBodyPtr sb1 = e[i].contact_geom1->get_single_body();
      Moby::SingleBodyPtr sb2 = e[i].contact_geom2->get_single_body();

      std::cout << "contact: " << sb1->id << " and " << sb2->id << std::endl;
      std::cout << "i = " << e[i].contact_impulse.get_linear() << std::endl;
      std::cout << "p = " << e[i].contact_point << std::endl;
      std::cout << "n = " << e[i].contact_normal << std::endl;
//      std::cout << "s = " << e[i].contact_tan1 << std::endl;
//      std::cout << "t = " << e[i].contact_tan2 << std::endl;
//      std::cout << "muC = " << e[i].contact_mu_coulomb << std::endl;
//      std::cout << "muV = " << e[i].contact_mu_viscous << std::endl;
    }
  }
  std::cout << "<< end post_event_callback_fn(.)" << std::endl;
}
void controller_callback(Moby::DynamicBodyPtr dbp, double t, void*)
{
/*
  std::cout << ">> start controller_callback(.)" << std::endl;
  Moby::RCArticulatedBodyPtr
      part = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(dbp);
  Ravelin::VectorNd x,xd;
  static double last_t;
  double h = t-last_t;
  last_t = t;
  part->get_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->get_generalized_velocity( Moby::DynamicBody::eEuler,xd);

  const std::vector<Moby::RigidBodyPtr>& links = part->get_links();
  std::cout << "Time = " << t << std::endl;

  for(int i=0;i<links.size();i++){
    boost::shared_ptr<const Ravelin::Pose3d> Ipose = links[i]->get_inertial_pose();
    boost::shared_ptr<const Ravelin::Pose3d> Lpose = links[i]->get_pose();

    std::cout << links[i]->id << std::endl;
    std::cout << "Ipose x = " << Ravelin::Pose3d::calc_relative_pose(Ipose,Moby::GLOBAL).x << std::endl;
    std::cout << "Lpose x = " << Ravelin::Pose3d::calc_relative_pose(Lpose,Moby::GLOBAL).x << std::endl;
    if(i != 0){
      boost::shared_ptr<const Ravelin::Pose3d> Jpose = links[i]->get_inner_joint_explicit()->get_pose();
      std::cout << "Jpose x = " << Ravelin::Pose3d::calc_relative_pose(Jpose,Moby::GLOBAL).x << std::endl;
    }
  }
  std::cout << "x = " << x << std::endl;
  std::cout << "v = " << xd << std::endl;



//  std::cout << "x =\n\t"
//            << RPY[0] << "\n\t"
//            << RPY[1] << "\n\t"
//            << RPY[2] << "\n\t"
//            << x[0] - THETA_SW_OFFSET << "\n\nxd =\n\t"
//            << dRPY[0] << "\n\t"
//            << dRPY[1] << "\n\t"
//            << dRPY[2] << "\n\t"
//            << xd[0] << "\n\t";
  std::cout << "<< end controller_callback(.)" << std::endl;
*/
}

// ============================================================================
// ================================ CALLBACKS =================================

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  Moby::RCArticulatedBodyPtr robot;

  // get a reference to the EventDrivenSimulator instance and the robot
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<Moby::EventDrivenSimulator>(i->second);

    // find the robot reference
    if (!robot)
      robot = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(i->second);
  }

  // make sure the robot was found
  assert(robot);

  // set the controller
  robot->controller = &controller;
}
} // end extern C

