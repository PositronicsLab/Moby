/*****************************************************************************
 * Controller for mobile robot
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>

#undef USE_INV_DYN
#ifdef USE_INV_DYN
#include <Pacer/controller.h>
#include <Pacer/robot.h>
#endif

#ifdef USE_INV_DYN
boost::shared_ptr<Pacer::Robot> pacer_robot;
#endif

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

// global vars
boost::shared_ptr<Moby::EventDrivenSimulator> sim;
RigidBodyPtr left_wheel_link, right_wheel_link, chassis_link;

// set the desired wheel speeds
const double UL = .10;
const double UR = 0.05;

// get the step size
double STEP_SIZE = 1e-5 * 5.0;
// calculates inverse dynamics torques under the no slip model
void calc_inverse_dynamics(RCArticulatedBodyPtr robot, const VectorNd& qdd, VectorNd& tau)
{
  #ifdef USE_INV_DYN
  const unsigned LEFT = 0, RIGHT = 1;

  // get the state and velocity
  VectorNd robot_q, robot_dq;
  robot->get_generalized_coordinates(DynamicBody::eEuler, robot_q);
  robot->get_generalized_velocity(DynamicBody::eSpatial, robot_dq);

  // setup the base pose
  shared_ptr<Pose3d> base_pose(new Pose3d);
  base_pose->x = Origin3d(robot_q[2], robot_q[3], robot_q[4]);
  base_pose->q = Quatd(robot_q[5], robot_q[6], robot_q[7], robot_q[8]);

  // setup the base velocity
  SVelocityd base_xd(robot_dq[5], robot_dq[6], robot_dq[7], robot_dq[2], robot_dq[3], robot_dq[4]);

  // setup q and qd
  std::map<std::string, double> q, qd;
  q["0left_wheel_hinge"] = robot_q[LEFT];
  q["0right_wheel_hinge"] = robot_q[RIGHT];
  qd["0left_wheel_hinge"] = robot_dq[LEFT];
  qd["0right_wheel_hinge"] = robot_dq[RIGHT];

  // get v, M, N, S/T, f
  boost::shared_ptr<const Pacer::RobotData> data = Pacer::Robot::gen_vars_from_model(q, qd, base_pose, base_xd, pacer_robot);
  const VectorNd& v = data->generalized_qd;
  const MatrixNd& N = data->N;
  const MatrixNd& M = data->M;
  const MatrixNd& ST = data->D;
  const VectorNd& f = data->generalized_fext;

  // define DT
  const double DT = STEP_SIZE;

  // call the method
  VectorNd cf;
  if (!Pacer::Controller::inverse_dynamics_no_slip(v, qdd, M, N, ST, f, DT, tau, cf))
    throw std::runtime_error("Could not solve inverse dynamics!");
  #else
  tau.set_zero(2);
  #endif
}

// setup the controller callback
void controller(DynamicBodyPtr body, double t, void*)
{
  const unsigned LEFT = 0, RIGHT = 1;
  static bool first_time = true;
  static Origin3d x0;

  // get the robot
  RCArticulatedBodyPtr robot = boost::dynamic_pointer_cast<RCArticulatedBody>(body);

  // get the state and velocity
  VectorNd q, dq;
  robot->get_generalized_coordinates(DynamicBody::eEuler, q);
  robot->get_generalized_velocity(DynamicBody::eSpatial, dq);

  // see whether this is the first time this function is run
  if (first_time)
  {
    first_time = false;

    // init x(t_0), y(t_0), theta(t_0)
    x0[0] = q[2];
    x0[1] = q[3];
    Quatd quat(q[5], q[6], q[7], q[8]);
    Matrix3d R = quat;
    x0[2] = std::atan2(R(2,0), R(0,0));
  }

  // setup the PD controller
  const double KV = .1;

  // set dq_des, ddq_des;
  double dq_des[2];
  VectorNd ddq_des(2);
  dq_des[LEFT] = UL;
  dq_des[RIGHT] = UR;
  ddq_des[LEFT] = ddq_des[RIGHT] = 0.0;

  const double DT = STEP_SIZE;
  ddq_des[LEFT] = (dq_des[LEFT] - dq[LEFT])/DT; 
  ddq_des[RIGHT] = (dq_des[RIGHT] - dq[RIGHT])/DT;

  // compute inverse dynamics torques
  VectorNd tau(2);
  calc_inverse_dynamics(robot, ddq_des, tau);

//std::cout << "L: " << dq[0] << " R: " << dq[1] << std::endl;
  // setup the joint torques
  VectorNd fleft(1), fright(1);
  fleft[0] = tau[LEFT] + KV*(dq_des[LEFT] - dq[LEFT]);
  fright[0] = tau[RIGHT] + KV*(dq_des[RIGHT] - dq[RIGHT]);

  // collect state data
  std::ofstream out("state.data", std::ostream::app);
  out << t;
  for (unsigned i=0; i< q.size(); i++)
    out << " " << q[i];
  for (unsigned i=0; i< dq.size(); i++)
    out << " " << dq[i];
  out << std::endl;
  out.close();

  // apply the torques
  JointPtr left = robot->get_joints()[0];
  JointPtr right = robot->get_joints()[1];
  assert(left->id == "left_wheel_hinge");
  assert(right->id == "right_wheel_hinge");
  left->add_force(fleft);
  right->add_force(fright);
}

void contact_callback_fn(std::vector<Moby::UnilateralConstraint>& e,
                            boost::shared_ptr<void> empty)
{
  const unsigned LEFT = 0, RIGHT = 1, CHASSIS = 2;

  #ifdef USE_INV_DYN
  // clear all existing contact data
  std::vector<Pacer::EndEffector>& eefs = pacer_robot->get_end_effectors();
  for (unsigned i=0; i< eefs.size(); i++)
  {
    eefs[i].normal.clear();
    eefs[i].tan1.clear();
    eefs[i].tan2.clear();
    eefs[i].point.clear();
    eefs[i].mu_coulomb.clear();
  }

  // process contacts
  for(unsigned i=0;i<e.size();i++){
    if (e[i].constraint_type == Moby::UnilateralConstraint::eContact)
    {
      Moby::SingleBodyPtr sb1 = e[i].contact_geom1->get_single_body();
      Moby::SingleBodyPtr sb2 = e[i].contact_geom2->get_single_body();
      if (sb1 == left_wheel_link || sb2 == left_wheel_link)
      {
        eefs[LEFT].point.push_back(e[i].contact_point);
        eefs[LEFT].normal.push_back(e[i].contact_normal);
        eefs[LEFT].tan1.push_back(e[i].contact_tan1);
        eefs[LEFT].tan2.push_back(e[i].contact_tan2);
        eefs[LEFT].mu_coulomb.push_back(100.0);
        if (sb2 == left_wheel_link)
          eefs[LEFT].normal.back() = eefs[LEFT].normal.back();
      }
      else if (sb1 == right_wheel_link || sb2 == right_wheel_link)
      {
        eefs[RIGHT].point.push_back(e[i].contact_point);
        eefs[RIGHT].normal.push_back(e[i].contact_normal);
        eefs[RIGHT].tan1.push_back(e[i].contact_tan1);
        eefs[RIGHT].tan2.push_back(e[i].contact_tan2);
        eefs[RIGHT].mu_coulomb.push_back(100.0);
        if (sb2 == right_wheel_link)
          eefs[RIGHT].normal.back() = eefs[RIGHT].normal.back();
      }
      else if (sb1 == chassis_link || sb2 == chassis_link)
      {
        eefs[CHASSIS].point.push_back(e[i].contact_point);
        eefs[CHASSIS].normal.push_back(e[i].contact_normal);
        eefs[CHASSIS].tan1.push_back(e[i].contact_tan1);
        eefs[CHASSIS].tan2.push_back(e[i].contact_tan2);
        eefs[CHASSIS].mu_coulomb.push_back(100.0);
        if (sb2 == chassis_link)
          eefs[CHASSIS].normal.back() = eefs[CHASSIS].normal.back();
      }
    }
  }
  #endif
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

  // setup the pre-event callback
  sim->constraint_callback_fn = &contact_callback_fn;

  // make sure the robot was found
  assert(robot);

  // set the controller
  robot->controller = &controller;

  // get the chassis, left, and right wheels
  const std::string CHASSIS_LINK_ID = "chassis";
  const std::string LEFT_WHEEL_LINK_ID = "left_wheel_link";
  const std::string RIGHT_WHEEL_LINK_ID = "right_wheel_link";
  chassis_link = robot->find_link(CHASSIS_LINK_ID);
  left_wheel_link = robot->find_link(LEFT_WHEEL_LINK_ID);
  right_wheel_link = robot->find_link(RIGHT_WHEEL_LINK_ID);

  // init Pacer
  #ifdef USE_INV_DYN
  const std::string model_file = "mrobot.sdf";
  const std::string pacer_file = "mrobot.pacer";
  pacer_robot = boost::shared_ptr<Pacer::Robot>(new Pacer::Robot(model_file, pacer_file));
  #endif

  // determine the velocity of the robot's base using initial conditions
  const double WHEEL_RAD = 0.11;
  const double AXLE_LEN = 0.34;
  const double THETA0 = 0.0;
  double xd0 = WHEEL_RAD * 0.5 * (UL + UR) * std::cos(THETA0);
  double yd0 = WHEEL_RAD * 0.5 * (UL + UR) * std::sin(THETA0);
  double thetad0 = WHEEL_RAD/AXLE_LEN * (UR - UL);

  // set the velocity of the robot's base
  const unsigned X = 0, Y = 1, THETA = 2;
  SVelocityd base_xd(GLOBAL);
  Vector3d lv(GLOBAL), av(GLOBAL);
  lv.set_zero();
  av.set_zero();
  lv[X] = xd0;
  lv[Y] = yd0;
  av[THETA] = thetad0;
  base_xd.set_angular(av);
  base_xd.set_linear(lv);
  robot->get_base_link()->set_velocity(base_xd);

  // set the velocities at the robot's wheels
  JointPtr left = robot->get_joints()[0];
  JointPtr right = robot->get_joints()[1];
  assert(left->id == "left_wheel_hinge");
  assert(right->id == "right_wheel_hinge");
  left->qd[0] = UL;
  right->qd[0] = UR;

  // update the robot's link velocities
  robot->update_link_velocities();

  // clear state data
  std::ofstream out("state.data");
  out.close();
}
} // end extern C

