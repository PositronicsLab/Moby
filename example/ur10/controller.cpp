/*****************************************************************************
 * "Controller" for grasping 
 ****************************************************************************/
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>
#include <stdlib.h>

// define below to use inverse dynamics
#define USE_INV_DYN

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;
using boost::dynamic_pointer_cast;

Moby::RCArticulatedBodyPtr robot;
boost::shared_ptr<TimeSteppingSimulator> sim;
boost::shared_ptr<GravityForce> grav;
std::map<std::string, double> qd_des;

void check_constraint_num(std::vector<Constraint>& constraints, boost::shared_ptr<void> data)
{
  std::vector<Constraint>& c_constraints = sim->get_rigid_constraints();

  // get the mapping from joint names to gc indices
  const std::vector<shared_ptr<Jointd> >& joints = robot->get_joints();

  // setup the mapping from joint ids to joints 
  std::map<std::string, shared_ptr<Jointd> > mapping;
  for (unsigned i=0; i< joints.size(); i++)
    mapping[joints[i]->joint_id] = joints[i];

  // setup inverse dynamics constraints
  for (std::map<std::string, double>::const_iterator i = qd_des.begin(); i != qd_des.end(); i++)
  {
    Constraint c;
    c.constraint_type = Constraint::eInverseDynamics;
    c.inv_dyn_joint = dynamic_pointer_cast<Joint>(mapping[i->first]);
    assert(c.inv_dyn_joint->num_dof() == 1);
    c.qdot_des.resize(1); 
    c.qdot_des[0] = i->second;
    constraints.push_back(c);
  }
}

VectorNd& controller(shared_ptr<ControlledBody> body, VectorNd& u, double t, void*)
{
  // get the generalized coordinates and velocity
  VectorNd q, qd;
  robot->get_generalized_coordinates_euler(q);
  robot->get_generalized_velocity(DynamicBodyd::eEuler, qd);

  // get the mapping from joint names to gc indices
  const std::vector<shared_ptr<Jointd> >& joints = robot->get_joints();

  // setup the mapping from joint ids to coordinate indices
  std::map<std::string, unsigned> mapping;
  for (unsigned i=0; i< joints.size(); i++)
    mapping[joints[i]->joint_id] = joints[i]->get_coord_index();

  // determine the desired position, velocity, and acceleration for the controller 
  const double PERIOD = 5.0;
  const double AMP = 0.5;
  const double SMALL_AMP = AMP*0.1;
  double sh_pan_q_des = std::sin(t)*AMP*PERIOD;
  double sh_pan_qd_des = std::cos(t)*AMP*PERIOD;
  qd_des["shoulder_pan_joint"] = sh_pan_qd_des; 
  double sh_lift_q_des = std::sin(t*2.0)*SMALL_AMP*PERIOD*2.0;
  double sh_lift_qd_des = std::cos(t*2.0)*SMALL_AMP*PERIOD*2.0;
  qd_des["shoulder_lift_joint"] = sh_lift_qd_des; 
  double elbow_q_des = std::sin(t*2.0/3.0)*AMP*PERIOD*2.0/3.0;
  double elbow_qd_des = std::cos(t*2.0/3.0)*AMP*PERIOD*2.0/3.0;
  qd_des["elbow_joint"] = elbow_qd_des;
  double wrist1_q_des = std::sin(t*1.0/7.0)*AMP*PERIOD*1.0/7.0;
  double wrist1_qd_des = std::cos(t*1.0/7.0)*AMP*PERIOD*1.0/7.0;
  qd_des["wrist_1_joint"] = wrist1_qd_des;
  double wrist2_q_des = std::sin(t*2.0/11.0)*AMP*PERIOD*2.0/11.0;
  double wrist2_qd_des = std::cos(t*2.0/11.0)*AMP*PERIOD*2.0/11.0;
  qd_des["wrist_2_joint"] = wrist2_qd_des;
  double wrist3_q_des = std::sin(t*3.0/13.0)*AMP*PERIOD*3.0/13.0;
  double wrist3_qd_des = std::cos(t*3.0/13.0)*AMP*PERIOD*3.0/13.0;
  qd_des["wrist_3_joint"] = wrist3_qd_des;

  // zero forces
  u.set_zero(robot->num_generalized_coordinates(DynamicBodyd::eSpatial));

  // set finger actuator forces
  u[mapping["l_finger_actuator"]] = 100.0;
  u[mapping["r_finger_actuator"]] = -100.0;

  #ifdef USE_INV_DYN
  return u;
  #endif

  // compute the errors
  double sh_pan_q_err = (sh_pan_q_des - q[mapping["shoulder_pan_joint"]]);
  double sh_pan_qd_err = (sh_pan_qd_des - qd[mapping["shoulder_pan_joint"]]);
  double sh_lift_q_err = (sh_lift_q_des - q[mapping["shoulder_lift_joint"]]);
  double sh_lift_qd_err = (sh_lift_qd_des - qd[mapping["shoulder_lift_joint"]]);
  double elbow_q_err = (elbow_q_des - q[mapping["elbow_joint"]]);
  double elbow_qd_err = (elbow_qd_des - qd[mapping["elbow_joint"]]);
  double wrist1_q_err = (wrist1_q_des - q[mapping["wrist_1_joint"]]);
  double wrist1_qd_err = (wrist1_qd_des - qd[mapping["wrist_1_joint"]]);
  double wrist2_q_err = (wrist2_q_des - q[mapping["wrist_2_joint"]]); 
  double wrist2_qd_err = (wrist2_qd_des - qd[mapping["wrist_2_joint"]]); 
  double wrist3_q_err = (wrist3_q_des - q[mapping["wrist_3_joint"]]); 
  double wrist3_qd_err = (wrist3_qd_des - qd[mapping["wrist_3_joint"]]); 

  // setup gains
  const double SH_KP = 300.0, SH_KV = 120.0;
  const double EL_KP = 60.0, EL_KV = 24.0;
  const double WR_KP = 15.0, WR_KV = 6.0;
 
  // compute the actuator forces
  double sh_pan_f = SH_KP*sh_pan_q_err + SH_KV*sh_pan_qd_err;
  double sh_lift_f = SH_KP*sh_lift_q_err + SH_KV*sh_lift_qd_err;
  double elbow_f = EL_KP*elbow_q_err + EL_KV*elbow_qd_err;
  double wrist1_f = WR_KP*wrist1_q_err + WR_KV*wrist1_qd_err;
  double wrist2_f = WR_KP*wrist2_q_err + WR_KV*wrist2_qd_err;
  double wrist3_f = WR_KP*wrist3_q_err + WR_KV*wrist3_qd_err;

  // set the actuator forces for the arm
  u.set_zero(robot->num_generalized_coordinates(DynamicBodyd::eSpatial));
  u[mapping["shoulder_pan_joint"]] = sh_pan_f;
  u[mapping["shoulder_lift_joint"]] = sh_lift_f;
  u[mapping["elbow_joint"]] = elbow_f;
  u[mapping["wrist_1_joint"]] = wrist1_f;
  u[mapping["wrist_2_joint"]] = wrist2_f;
  u[mapping["wrist_3_joint"]] = wrist3_f;

  return u; 
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  const unsigned Z = 2;

  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);
    
    if (i->first == "ur10_schunk_hybrid")
      robot = boost::dynamic_pointer_cast<RCArticulatedBody>(i->second);
  }
  assert(sim);
  #ifdef USE_INV_DYN
  sim->constraint_callback_fn = &check_constraint_num;
  #endif
  assert(robot);
  robot->controller = &controller; 

  // make the base fixed
  robot->set_floating_base(false);

  // sets the starting velocity for the robot joints
  const std::vector<shared_ptr<Jointd> >& joints = robot->get_joints();
  const double PERIOD = 5.0;
  const double AMP = 0.5;
  const double SMALL_AMP = AMP*0.1;
  std::map<std::string, double> qd_init;
  qd_init["shoulder_pan_joint"] = std::cos(0)*AMP*PERIOD;
  qd_init["shoulder_lift_joint"] = std::cos(0)*SMALL_AMP*PERIOD*2.0;
  qd_init["elbow_joint"] = std::cos(0)*AMP*PERIOD*2.0/3.0;
  qd_init["wrist_1_joint"] = std::cos(0)*AMP*PERIOD*1.0/7.0;
  qd_init["wrist_2_joint"] = std::cos(0)*AMP*PERIOD*2.0/11.0;
  qd_init["wrist_3_joint"] = std::cos(0)*AMP*PERIOD*3.0/13.0;
  VectorNd qd;
  robot->get_generalized_velocity(DynamicBodyd::eEuler, qd);
  for (unsigned i=0; i< joints.size(); i++)
    qd[joints[i]->get_coord_index()] = qd_init[joints[i]->joint_id];
  robot->set_generalized_velocity(DynamicBodyd::eEuler, qd);

  // get the poses for all links
  const std::vector<shared_ptr<RigidBodyd> >& links = robot->get_links();
  for (unsigned i=0; i< links.size(); i++)
  {
    Pose3d T = *links[i]->get_pose();
    T.update_relative_pose(GLOBAL);
    std::cout << "pose of " << links[i]->body_id << ": " << T << std::endl;
  }
}
} // end extern C

