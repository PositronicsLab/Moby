#include <cmath>
#include <sys/types.h>
#include <sys/times.h>
#include <fstream>
#include <Moby/XMLReader.h>
#include <Moby/Simulator.h>
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/DeformableBody.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/GravityForce.h>
#include <Moby/LinAlg.h>
#include <Moby/Constants.h>
#include <Moby/RNEAlgorithm.h>
#ifdef USE_OSG
#include <osgDB/WriteFile>
#endif
//#include <osgDB/ReaderWriterSTL>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/*****************************************************************************
 * Controller for getting the robot arm to grasp a ball.
 * (to be used with bandit-arm-kineout.xml)
 ****************************************************************************/

/// plugin must be "extern C"

extern "C" {

int ITER;

// The position of the ball in the grippers' coordinate systems
bool first = true;
Vector3 _ball_grip_left;
Vector3 _ball_grip_right;

// The ball body
RigidBodyPtr left_gripper, right_gripper;
DeformableBodyPtr ball;

/// Controls the robot
void control_PID(RCArticulatedBodyPtr robot, Real time)
{
  Real ACCEL = 1000.0;

  // determine joint positions for two dof of shoulder
  const Real SPEED = 100;
  const Real MAG = 0.25;
  Real lsh1 = 1.30 + std::sin(time*SPEED)*MAG;
  Real lsh2 = std::sin(2*time*SPEED)*MAG;
  Real lsh1_qd = std::cos(time*SPEED)*MAG*SPEED;
  Real lsh2_qd = std::cos(2*time*SPEED)*SPEED*MAG*2;
  Real lsh1_qdd = -std::sin(time*SPEED)*SPEED*SPEED*MAG;
  Real lsh2_qdd = -std::sin(2*SPEED*time)*SPEED*SPEED*MAG*4;

  // setup desired joint positions
  std::map<std::string, Real> q_des;
  //q_des["left-shoulder1-joint"] = 1.28;
  q_des["left-shoulder1-joint"] = 0.0;
  q_des["left-shoulder2-joint"] = 0.0;
  q_des["left-bicep-joint"] = 0.0; 
  q_des["left-elbow-joint"] = 0.0;
  q_des["left-forearm-joint"] = 0.0;
  q_des["left-hand-joint"] = 0.0;
  q_des["left-claw-left-joint"] = 0.1;
  q_des["left-claw-right-joint"] = -0.1;

  //q_des["left-claw-left-joint"] = -.198998;
  //q_des["left-claw-right-joint"] = .199998;
  //q_des["left-claw-left-joint"] = -0.25;
  //q_des["left-claw-right-joint"] = 0.25;
  //q_des["left-claw-left-joint"] = -0.37;
  //q_des["left-claw-right-joint"] = 0.37;

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
    out << time << " " << q_des[i->first->id] << " " << i->first->q[0] << std::endl;
    out.close();
    out.open(fname2.c_str(), std::ostream::app);
    out << time << " " << qd_des[i->first->id] << " " << i->first->qd[0] << std::endl;
    out.close();

    // add feedback torque to joints
    Real perr = q_des[i->first->id] - i->first->q[0];
    Real derr = qd_des[i->first->id] - i->first->qd[0];
    VectorN fb_torque(1);
    fb_torque[0] = perr*KP + derr*KV;
    i->first->add_force(fb_torque);

  }
}

/// The main control loop
void controller(DynamicBodyPtr robot, Real time, void* data)
{
  // determine coordinates of ball in gripper coordinate frames
  if (first)
  {
    _ball_grip_left = Matrix4::inverse_transform(left_gripper->get_transform()) .mult_point(ball->get_position());  
    _ball_grip_right = Matrix4::inverse_transform(right_gripper->get_transform()).mult_point(ball->get_position());  
    first = false;
  }
  else
  {
    // output the combined error from the starting position w.r.t. both grippers
    std::ofstream out("error.ball", std::ios::app);
    Vector3 ball_grip_left = Matrix4::inverse_transform(left_gripper->get_transform()).mult_point(ball->get_position());  
    Vector3 ball_grip_right = Matrix4::inverse_transform(right_gripper->get_transform()).mult_point(ball->get_position());  
    Real err = std::sqrt((ball_grip_left - _ball_grip_left).norm_sq() + (ball_grip_right - _ball_grip_right).norm_sq());
    out << time << " " << err << std::endl;
    out.close();
  }

  control_PID(dynamic_pointer_cast<RCArticulatedBody>(robot), time);
}

boost::shared_ptr<EventDrivenSimulator> sim;

static int step = 0;

// callback for EventDrivenSimulator
void event_callback_fn( std::vector<Event>& events, boost::shared_ptr<void> p ) {

  std::vector<Event>::iterator it;

  int i = 0;
  Real time;
  int iter = ITER;

  std::ofstream collog;
  std::stringstream ss_collog;
  ss_collog << "collog_" << step << "_" << iter;
  std::string collog_name = ss_collog.str();

  collog.open( collog_name.c_str() );
  if( collog.fail() )
    throw std::runtime_error("Unable to open file for writing collision log in bandit-controller.cpp::event_callback_fn");

  for ( it = events.begin() ; it < events.end(); it++ ) {
    Event event = *it;

    if( event.event_type != Event::eContact ) continue;

    time = event.t;

    std::stringstream ss_event;
    ss_event << "event_" << step << "_" << iter << "_" << i++ << ".wrl";
    std::string event_file_name = ss_event.str();


    SingleBodyPtr geom1 = event.contact_geom1->get_single_body();
    SingleBodyPtr geom2 = event.contact_geom2->get_single_body();

    event.write_vrml( event_file_name, 0.01, 0.1 );
    //collog << "event:" << event_file_name << " geom1:" << event.contact_geom1->id << " geom2:" << event.contact_geom2->id << "\n";

    //collog << "event:" << event_file_name << "\tgeom1:" << geom1->id << "\tgeom2:" << geom2->id;
    if (event.constraint_joint) {
    	VectorN q = event.constraint_joint->q;

	collog << "\tjoint_angles(" << q.size() << "):"; 
	for( unsigned int i = 0; i < q.size(); i++ )
	    collog <<  q[i] << " ";
    }
    collog << "\n";
  }
  collog.close();
/*
  std::stringstream ss_state;
  ss_state << "state_" << step << "_" << iter << ".stl";
  std::string state_file_name = ss_state.str();

  osg::Node* node = sim->get_persistent_vdata();
  osgDB::writeNodeFile(*node, state_file_name);
*/
  std::stringstream ss_state_trans;
  ss_state_trans << "state_trans_" << step << "_" << iter << ".stl";
  std::string trans_state_file_name = ss_state_trans.str();

  #ifdef USG_OSG
  osg::Node* node_trans = sim->get_transient_vdata();
  osgDB::writeNodeFile(*node_trans, trans_state_file_name);

  std::stringstream ss_state_persist;
  ss_state_persist << "state_persist_" << step << "_" << iter << ".stl";
  std::string persist_state_file_name = ss_state_persist.str();

  osg::Node* node_persist = sim->get_persistent_vdata();
  osgDB::writeNodeFile(*node_persist, persist_state_file_name);
  #endif

  step++;  
}

void init(void* separator, const std::map<std::string, BasePtr>& read_map, Real time)
{
  if (read_map.find("0x978939c") == read_map.end())
    throw std::runtime_error("bandit-controller.cpp:init()- unable to find simulator!");
  // get a reference to the EventDrivenSimulator instance
  sim = dynamic_pointer_cast<EventDrivenSimulator>( read_map.find( "0x978939c" )->second );
  // register the event callback function
  sim->event_callback_fn = &event_callback_fn;

  // find the robot
  if (read_map.find("bandit-arm") == read_map.end())
    throw std::runtime_error("bandit-controller.cpp:init()- unable to find bandit arm!");
  DynamicBodyPtr robot = dynamic_pointer_cast<DynamicBody>(read_map.find("bandit-arm")->second);
  if (!robot)
    throw std::runtime_error("bandit-controller.cpp:init()- unable to cast bandit-arm to type DynamicBody");

  // setup the controller
  robot->controller = controller;

  // find the ball
  if (read_map.find("ball") == read_map.end())
    throw std::runtime_error("bandit-controller.cpp:init()- unable to find ball!");
  ball = dynamic_pointer_cast<DeformableBody>(read_map.find("ball")->second);
  if (!ball)
    throw std::runtime_error("bandit-controller.cpp:init()- unable to cast ball to type DeformableBody");

  // find the grippers
  if (read_map.find("claw-left-left") == read_map.end())
    throw std::runtime_error("bandit-controller.cpp:init()- unable to find claw-left-left!");
  left_gripper = dynamic_pointer_cast<RigidBody>(read_map.find("claw-left-left")->second);
  if (!left_gripper)
    throw std::runtime_error("bandit-controller.cpp:init()- unable to cast claw-left-right to type RigidBody");
  if (read_map.find("claw-left-right") == read_map.end())
    throw std::runtime_error("bandit-controller.cpp:init()- unable to find claw-left-right!");
  right_gripper = dynamic_pointer_cast<RigidBody>(read_map.find("claw-left-right")->second);
  if (!right_gripper)
    throw std::runtime_error("bandit-controller.cpp:init()- unable to cast claw-left-right to type RigidBody");
}

} // end extern C

