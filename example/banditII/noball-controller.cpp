#include <cmath>
#include <sys/types.h>
#include <sys/times.h>
#include <fstream>
#include <Moby/XMLReader.h>
#include <Moby/Simulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/DeformableBody.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/GravityForce.h>
#include <Moby/LinAlg.h>
#include <Moby/Constants.h>
#include <Moby/RNEAlgorithm.h>
#include <Moby/EventDrivenSimulator.h>

#include <osgDB/WriteFile>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// plugin must be "extern C"

extern "C" {

int ITER;

/// Controls the robot
void control_PID(RCArticulatedBodyPtr robot, Real time)
{
  Real ACCEL = 1000.0;

  // setup desired joint positions
  std::map<std::string, Real> q_des;
  q_des["left-shoulder1-joint"] = 0.0;
  q_des["left-shoulder2-joint"] = 0.0;
  q_des["left-bicep-joint"] = 0.0; 
  q_des["left-elbow-joint"] = 0.0;
  q_des["left-forearm-joint"] = 0.0;
  q_des["left-hand-joint"] = 0.0;
  q_des["left-claw-left-joint"] = 0.1;
  q_des["left-claw-right-joint"] = -0.1;

  // setup desired joint velocities
  std::map<std::string, Real> qd_des;
  qd_des["left-shoulder1-joint"] = 0.0;
  qd_des["left-shoulder2-joint"] = 0.0;
  qd_des["left-bicep-joint"] = 0.0; 
  qd_des["left-elbow-joint"] = 0.0;
  qd_des["left-forearm-joint"] = 0.0;
  qd_des["left-hand-joint"] = 0.0;
  qd_des["left-claw-left-joint"] = 0.0;
  qd_des["left-claw-right-joint"] = 0.0;

  // setup gains
  std::map<std::string, std::pair<Real, Real> > gains;
  gains["left-shoulder1-joint"] = std::make_pair(300,120);
  gains["left-shoulder2-joint"] = std::make_pair(300,120);
  gains["left-bicep-joint"] = std::make_pair(100,40);
  gains["left-elbow-joint"] = std::make_pair(60,24);
  gains["left-forearm-joint"] = std::make_pair(25,10);
  gains["left-hand-joint"] = std::make_pair(15,6);
  gains["left-claw-left-joint"] = std::make_pair(15,6);
  gains["left-claw-right-joint"] = std::make_pair(15,6);

  // compute inverse dynamics
  std::map<RigidBodyPtr, RCArticulatedBodyInvDynData> inv_dyn_data;
  for (unsigned i=1; i< robot->get_links().size(); i++)
  {
    // setup the inverse dynamics data
    RCArticulatedBodyInvDynData id_data;
    id_data.fext = robot->get_links()[i]->sum_forces(); 
    id_data.text = robot->get_links()[i]->sum_torques(); 
    JointPtr joint(robot->get_links()[i]->get_inner_joint_implicit());
    id_data.qdd = VectorN(1);

    // update the acceleration of the joint
    if (joint->id == "left-claw-right-joint")
      id_data.qdd[0] = -ACCEL;
    else if (joint->id == "left-claw-left-joint")
      id_data.qdd[0] = ACCEL;
    else
      id_data.qdd[0] = 0;
    inv_dyn_data[robot->get_links()[i]] = id_data;
  }

  // compute inverse dynamics
  RNEAlgorithm rne;
  std::map<JointPtr, VectorN> actuator_forces = rne.calc_inv_dyn(robot, inv_dyn_data);

  // clear and set motor torques
  for (std::map<JointPtr, VectorN>::const_iterator i = actuator_forces.begin(); i != actuator_forces.end(); i++)
  {
    // reset motor torque
    i->first->reset_force();

    // add computed torque
    i->first->add_force(i->second);

    // get the two gains
    const Real KP = gains[i->first->id].first / 10.0;
    const Real KV = gains[i->first->id].second / 10.0;
    
    // for outputting desired position and velocity
    std::string fname1 = i->first->id + ".pos";
    std::string fname2 = i->first->id + ".vel";
    std::ofstream out(fname1.c_str(), std::ofstream::app);
    out << q_des[i->first->id] << " " << i->first->q[0] << std::endl;
    out.close();
    out.open(fname2.c_str(), std::ostream::app);
    out << qd_des[i->first->id] << " " << i->first->qd[0] << std::endl;
    out.close();
    
    // add feedback torque to joints
    Real perr = q_des[i->first->id] - i->first->q[0];
    Real derr = qd_des[i->first->id] - i->first->qd[0];
    VectorN fb_torque(1);
    fb_torque[0] = perr*KP + derr*KV;
    i->first->add_force(fb_torque);
    
  }
}

boost::shared_ptr<EventDrivenSimulator> sim;

static int step = 0;

// callback for EventDrivenSimulator
void event_callback_fn( std::vector<Event>& events, boost::shared_ptr<void> p ) {

  std::vector<Event>::iterator it;

  int i = 0;
  Real time;
  int iter = ITER;

  std::stringstream ss_state_trans;
  ss_state_trans << "state_trans_" << step << "_" << iter << ".stl";
  std::string trans_state_file_name = ss_state_trans.str();

  osg::Node* node_trans = sim->get_transient_vdata();
  osgDB::writeNodeFile(*node_trans, trans_state_file_name);

  std::stringstream ss_state_persist;
  ss_state_persist << "state_persist_" << step << "_" << iter << ".stl";
  std::string persist_state_file_name = ss_state_persist.str();

  osg::Node* node_persist = sim->get_persistent_vdata();
  osgDB::writeNodeFile(*node_persist, persist_state_file_name);

  step++;  
}

/// The main control loop
void controller(DynamicBodyPtr robot, Real time, void* data)
{
  control_PID(dynamic_pointer_cast<RCArticulatedBody>(robot), time);
}

void init(void* separator, const std::map<std::string, BasePtr>& read_map, Real time)
{
  if (read_map.find("0xa024ed4") == read_map.end())
    throw std::runtime_error("noball-controller.cpp:init()- unable to find simulator!");
  // get a reference to the EventDrivenSimulator instance
  sim = dynamic_pointer_cast<EventDrivenSimulator>( read_map.find( "0xa024ed4" )->second );
  // register the event callback function
  sim->event_callback_fn = &event_callback_fn;

  // find the robot
  if (read_map.find("bandit-arm") == read_map.end())
    throw std::runtime_error("noball-controller.cpp:init()- unable to find bandit arm!");
  DynamicBodyPtr robot = dynamic_pointer_cast<DynamicBody>(read_map.find("bandit-arm")->second);
  if (!robot)
    throw std::runtime_error("noball-controller.cpp:init()- unable to cast bandit-arm to type DynamicBody");

  // setup the controller
  robot->controller = controller;
}

} // end extern C

