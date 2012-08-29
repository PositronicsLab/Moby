#include <cmath>
#include <fstream>
#include <Moby/XMLReader.h>
#include <Moby/Simulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Joint.h>
#include <Moby/GravityForce.h>
#include <Moby/LinAlg.h>
#include <Moby/Constants.h>
#include <Moby/RNEAlgorithm.h>
#include <Moby/Vector2.h>

#ifdef USE_INVENTOR
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/So.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoMaterial.h>
#endif

// reads and controls a simulated two-wheel mobile robot
// uses computed torque control to move the robot forward at a constant
// velocity

using std::vector;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Moby;

extern unsigned ITER;
extern Real STEP_SIZE;

extern "C" {

// next two lines used only for forward/angular motion tests
Real ZZ = 0, ZZD = 0, ZZDD = .1;
Real ZTHETA = 0, ZTHETAD = 0, ZTHETADD = 0.1;
Real TOTAL_ERROR = 0;


#ifdef USE_INVENTOR
SoCylinder* circle;
SoSeparator* rpath, * dpath;
#endif

/// Computes the angle (around the y axis) of the robot
Real get_y_axis_angle(const Matrix3& R)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real NEAR_ZERO = std::sqrt(std::numeric_limits<Real>::epsilon());

  Vector3 YAXIS(0,1,0);

  // following code adapted from AAngle.cpp
  // ************************************************************
  Real acos_angle = (R(X,X) + R(Y,Y) + R(Z,Z) - 1)/2; 
  if (acos_angle > 1.0)
    acos_angle = 1.0;
  else if (acos_angle < -1.0)
    acos_angle = -1.0;

  // compute angle of rotation
  Real angle = std::acos(acos_angle);

  // if the angle of rotation is zero, don't need to do any more
  if (angle == 0.0)
    return 0.0;

  // if the angle of rotation is PI, solve for axis
  Real x,y,z;
  if (angle == M_PI)
  {
    x = std::sqrt((R(X,X)+1)/2);
    y = std::sqrt((R(Y,Y)+1)/2);
    z = std::sqrt((R(Z,Z)+1)/2);
  }
  else
  {
    Real constant = 1.0/(2.0*std::sin(angle));
    x = constant * (R(Z,Y) - R(Y,Z));
    y = constant * (R(X,Z) - R(Z,X));
    z = constant * (R(Y,X) - R(X,Y));
  }

  // make sure that x, y, and z are not NaN
  assert(!std::isnan(x));
  assert(!std::isnan(y));
  assert(!std::isnan(z));

  // get the length of determined axis [x y z]; if length is zero, angle is
  // zero...
  Real len = std::sqrt(x*x + y*y + z*z);
  if (len == 0.0)
    return 0.0;

  // normalize vector [x y z] (generally not necessary, but safe...)
  if (std::fabs(len - 1.0) > NEAR_ZERO)
  {
    Real leni = 1.0/len;
    x *= leni;
    y *= leni;
    z *= leni;
  }

  // form the determined axis and reverse the angle, if necessary
  Vector3 det_axis(x,y,z);
  return (Vector3::dot(det_axis, YAXIS) < 0) ? -angle : angle;
  // ************************************************************
  // adapted code ends...
}


/// Gets the minimum angle necessary to get to b from a
Real angle_diff(Real a, Real b)
{
  Real d = b - a;
  if (d < -M_PI) return d+2*M_PI;
  if (d > M_PI) return d-2*M_PI;
  return d;
}

/// Traces a figure-eight parametrically (period is for pi=[2pi*j,2pi*(j+1)] for j >= 0)
void fig8(Real t, Real& x, Real& xd, Real& y, Real& yd)
{
  x = -4*std::cos(t);
  y = 8*std::sin(t)*std::cos(t);
  xd = 4*std::sin(t);
  yd = 8*(std::cos(t)*std::cos(t) - std::sin(t)*std::sin(t));
}

/// Controls the base using errors and differentials
/**
 * \param x_des the desired x position (global frame)
 * \param xd_des the desired x velocity (global frame)
 * \param z_des the desired z position (global frame)
 * \param zd_des the desired z velocity (global frame)
 * \param theta_des the desired angle about the global y-axis 
 * \param thetad_des the desired angular velocity about the global y-axis 
 */
void controlBaseDiff(RCArticulatedBodyPtr robot, Real x_des, Real xd_des, Real z_des, Real zd_des, Real theta_des, Real thetad_des)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real KV = 50.0;
  const Real VSCALE = 0.001; 
  Real WHEEL_RADIUS, AXLE_LEN;

  // constants for pioneer
  WHEEL_RADIUS = 0.0825;
  AXLE_LEN = 0.3206;

  // get the base link
  RigidBodyPtr base = robot->get_links().front();
RigidBodyPtr lwheel = robot->get_links()[1];
RigidBodyPtr rwheel = robot->get_links()[2];

  // get the current base position, linear velocity, and linear acceleration
  const Vector3& base_x = base->get_position();
  const Vector3& base_xd = base->get_lvel();

  // get the current base orientation and angular velocity (around y-axis)
  Matrix3 R;
  base->get_transform().get_rotation(&R);
  Real base_theta = get_y_axis_angle(R);
  Real base_thetad = base->get_avel()[Y];

  // get left and right joints 
  JointPtr left = robot->find_joint("left-wheel-joint");
  JointPtr right = robot->find_joint("right-wheel-joint");

  // get the current joint positions for the joints
  Real left_q = left->q[0];
  Real right_q = right->q[0];

  // determine the direction that local z-axis is pointing
  Vector3 localZaxis;
  R.get_column(Z,localZaxis);
  Real vx = localZaxis[X];
  Real vz = localZaxis[Z];

  // determine how fast base rotates as a function of only one wheel moving
  const Real BASE_ROTATION_RATE = WHEEL_RADIUS / AXLE_LEN; 

  // setup the Jacobian for the base
  // column 1 is left wheel, column 2 is right
  MatrixN J(3,2);
  J(X,0) = vx*WHEEL_RADIUS/2.0;    J(X,1) = vx*WHEEL_RADIUS/2.0;
  J(Y,0) = vz*WHEEL_RADIUS/2.0;    J(Y,1) = vz*WHEEL_RADIUS/2.0;
  J(Z,0) = BASE_ROTATION_RATE;    J(Z,1) = -BASE_ROTATION_RATE;

  // log the desired base positions, velocities, and accelerations
  std::ofstream out("base.out", std::ios::app);
  out << x_des << " " << z_des << " " << theta_des << std::endl;
  out.close();

  // log the current base position and orientation
  out.open("base.true.out", std::ios::app);
  out << base_x[X] << " " << base_x[Z] << " " << base_theta << std::endl;
  out.close();
  
  // get the linear errors (x/z plane)
  Vector3 base_x_err = Vector3(x_des, base_x[Y], z_des) - base_x;
  Vector3 base_xd_err = Vector3(xd_des, base_x[Y], zd_des) - base_xd;

  // get the angular errors (about the y-axis)
  Real base_theta_err = angle_diff(base_theta, theta_des);
  Real base_thetad_err = thetad_des - base_thetad;

  // setup the operational-space differential vector; first dimension is
  // x-axis differential, second dimension is z-axis differential; third
  // dimension is y-axis angular differential
  VectorN xdiff(3);
  xdiff[0] = base_x_err[X] + VSCALE*base_xd_err[X];
  xdiff[1] = base_x_err[Z] + VSCALE*base_xd_err[Z];
  xdiff[2] = base_theta_err + VSCALE*base_thetad_err;

  // add to the total positional error
  TOTAL_ERROR += (base_x_err[X]*base_x_err[X]) + (base_x_err[Z]*base_x_err[Z]);
  
  // compute the SVD pseudo-inverse of the Jacobian
  MatrixN Jcross = Moby::LinAlg::pseudo_inverse(J);

  // compute the desired joint velocities
  VectorN qdiff = Jcross * xdiff;

  // ***********************************************
  // command the robot 
  // ***********************************************

  // do feedback control    
  VectorN force(1);  

  // left wheel joint
  JointPtr lw_joint = robot->find_joint("left-wheel-joint");
  Real lqd_t = lw_joint->qd[0];  
  force[0] = KV*(qdiff[0] - lqd_t);
  lw_joint->add_force(force);
  out.open("left",std::ios::app);
  out << force[0] << " " << lqd_t << std::endl;
  out.close();  

  // right wheel joint
  JointPtr rw_joint = robot->find_joint("right-wheel-joint");
  Real rqd_t = rw_joint->qd[0];  
  force[0] = KV*(qdiff[1] - rqd_t);
  rw_joint->add_force(force);
  out.open("right",std::ios::app);
  out << force[0] << " " << rqd_t << std::endl;
  out.close();  

  // write the base control 
  std::cout << "  desired base: " << x_des << " " << z_des << " " << theta_des << std::endl;
  std::cout << "  current base: " << base_x[X] << " " << base_x[Z] << " " << base_theta << std::endl;
  std::cout << "  desired base velocity: " << xd_des << " " << zd_des << " " << thetad_des << std::endl;
  std::cout << "  current base velocity: " << base_xd[X] << " " << base_xd[Z] << " " << base_thetad << std::endl;
  std::cout << "  forces (L/R): " << KV*(qdiff[0] - lqd_t) << " " << force[0] << std::endl;
  std::cout << "  wheel velocities (L/R): " << lqd_t << " " << rqd_t << std::endl;
/*
  std::cout << "  desired wheel velocities: " << qdiff << std::endl;
  std::cout << "  base velocity (predicted): " << (J*Vector2(lqd_t, rqd_t)) << std::endl;
  std::cout << "  jacobian: " << std::endl << J;
  std::cout << "  x-differential for IK: " << xdiff << std::endl;
  std::cout << "  mean positional error: " << TOTAL_ERROR/ITER << std::endl;
*/
}

/// Determines the value of theta for the robot
Real calc_theta(Real xdes, Real zdes, Real xcurrent, Real zcurrent)
{
  Vector2 dir(xdes-xcurrent, zdes-zcurrent);
  Real theta = 0.0;
  if (dir.norm() > 0.0)
  {
    dir.normalize();

    // NOTE: we have to negate the 'y' argument b/c positive y corresponds to
    // negative z the way things have been setup
    theta = std::atan2(-dir[1], dir[0]) - M_PI_2;
    if (theta < -M_PI)
      theta += 2 * M_PI;
  }

  return theta;
}

/// Controls the base by translating base commands into joint commands
/**
 * This method is retained for pedagogical reasons, though it is no longer
 * used.  It demonstrates how one would determine, analytically, the proper
 * joint accelerations necessary to achieve a desired base acceleration.
 */
// NOTE: this method has not been tested with new contact models
/*
std::map<std::string, Vector3> controlBaseAccel(RCArticulatedBodyPtr robot, Real x_des, Real xd_des, Real xdd_des, Real z_des, Real zd_des, Real zdd_des, Real theta_des, Real thetad_des, Real thetadd_des)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Real PERR_GAIN = 2.5, DERR_GAIN = 0, DDERR_GAIN = 0;
  const Real WHEEL_RADIUS = 0.11 * 1.0/8.0;
  const Real AXLE_LEN = 0.116*2;

  // compute forward dynamics
  robot->calc_fwd_dyn();

  // get the base link
  RigidBodyPtr base = robot->get_links().front();

  // get the current base position, linear velocity, and linear acceleration
  const Vector3& base_x_vec = base->get_position();
  const Vector3& base_xd_vec = base->get_lvel();
  const Vector3& base_xdd_vec = base->get_laccel();

  // get the current base orientation and angular velocity (around y-axis)
  Matrix3 R;
  base->get_transform().get_rotation(&R);
  AAngle aa(&R);
  if (aa._y < 0)
    aa._angle = -aa._angle;
  Real base_theta = aa._angle;
  Real base_thetad = base->get_avel()[Y];
  Real base_thetadd = base->get_aaccel()[Y];

  // form "true" base position, velocity, accel
  Vector3 base_x(base_x_vec[X], base_x_vec[Z], base_theta);
  Vector3 base_xd(base_xd_vec[X], base_xd_vec[Z], base_thetad);
  Vector3 base_xdd(base_xdd_vec[X], base_xdd_vec[Z], base_thetadd);

  // get left and right joints 
  JointPtr left = robot->find_joint("left-wheel-joint");
  JointPtr right = robot->find_joint("right-wheel-joint");

  // get the current joint positions, velocities, and accelerations for the 
  // joints
  Real left_q = left->q[0];
  Real right_q = right->q[0];
  Vector2 qd(left->qd[0], right->qd[0]);
  Vector2 qdd(left->qdd[0], right->qdd[0]);

  // log the desired base positions, velocities, and accelerations
  std::ofstream out("base-x.out", std::ios::app);
  out << x_des << " " << xd_des << " " << xdd_des << std::endl;
  out.close();
  out.open("base-z.out", std::ios::app);
  out << z_des << " " << zd_des << " " << zdd_des << std::endl;
  out.close();
  out.open("base-theta.out", std::ios::app);
  out << theta_des << " " << thetad_des << " " << thetadd_des << std::endl;
  out.close();

  // log the current base position and orientation
  out.open("base-x.true.out", std::ios::app);
  out << base_x[0] << " " << base_xd[0] << " " << base_xdd[0] << std::endl;
  out.close();
  out.open("base-z.true.out", std::ios::app);
  out << base_x[1] << " " << base_xd[1] << " " << base_xdd[1] << std::endl;
  out.close();
  out.open("base-theta.true.out", std::ios::app);
  out << base_x[2] << " " << base_xd[2] << " " << base_xdd[2] << std::endl;
  out.close();
  
  // determine the direction that local z-axis is pointing
  Vector3 localZaxis;
  R.get_column(Z,localZaxis);
  Real vx = localZaxis[X];
  Real vz = localZaxis[Z];

  // determine how fast base rotates as a function of only one wheel moving
  const Real BASE_ROTATION_RATE = WHEEL_RADIUS / AXLE_LEN; 

  // setup the Jacobian for the base
  // column 1 is left wheel, column 2 is right
  MatrixN J(3,2);
  J(X,0) = vx*WHEEL_RADIUS/2.0;    J(X,1) = vx*WHEEL_RADIUS/2.0;
  J(Y,0) = vz*WHEEL_RADIUS/2.0;    J(Y,1) = vz*WHEEL_RADIUS/2.0;
  J(Z,0) = BASE_ROTATION_RATE;    J(Z,1) = -BASE_ROTATION_RATE;

  // compute the SVD pseudo-inverse of the Jacobian
  MatrixN Jcross = Moby::LinAlg::pseudo_inverse(J);

  // setup vectors of desired base accelerations, velocities, and positions
  Vector3 base_xdd_des(xdd_des, zdd_des, thetadd_des);
  Vector3 base_xd_des(xd_des, zd_des, thetad_des);
  Vector3 base_x_des(x_des, z_des, theta_des);

  // determine errors 
  Vector3 base_x_err = base_x_des - base_x;
  Vector3 base_xd_err = base_xd_des - base_xd;
  Vector3 base_xdd_err = base_xdd_des - base_xdd;

  // determine desired joint positions -- since positions suffer from drift,
  // we'll utilize the positional error
  Vector2 q(left_q, right_q);
  VectorN q_diff = Jcross * base_x_err;
//  VectorN q_des = q + q_diff;
  VectorN q_des = q;

  // determine desired joint accelerations and velocities
  VectorN qdd_des = Jcross * (base_xdd_des + DDERR_GAIN*base_xdd_err);
//  VectorN qdd_des = Jcross * (base_xdd_des + PERR_GAIN*base_x_err + DERR_GAIN*base_xd_err);
  VectorN qd_des = Jcross * base_xd_des;
//  VectorN qd_des = Jcross * (base_xd_des + PERR_GAIN*base_x_err + DERR_GAIN*base_xd_err);
//  VectorN qd_des = qd;

  // log base control data
  out.open("base-control.out", std::ios::app);
  out << "Jacobian: " << std::endl << J;
  out << "xdd (des): " << base_xdd_des << std::endl;
  out << "xdd (current): " << base_xdd << std::endl;
  out << "xdd (err): " << (base_xdd_des - base_xdd) << std::endl;
  out << "xd (des): " << base_xd_des << std::endl;
  out << "xd (current): " << (base_xd_des - base_xd) << std::endl;
  out << "xd (error): " << (base_xd_des - base_xd) << std::endl;
  out << "x (des): " << base_x_des << std::endl;
  out << "x (current): " << base_x << std::endl;
  out << "x (error): " << base_x_err << std::endl;
  out << "qdd (des): " << qdd_des << std::endl;
  out << "qdd (current): " << qdd << std::endl;
  out << "qd (des): " << qd_des << std::endl;
  out << "qd (current): " << qd << std::endl;
  out << "q (des): " << q_des << std::endl;
  out << "q (current): " << q << std::endl;
  out << std::endl;
  out.close();

  // make a map of commands
  std::map<std::string, Vector3> commands;
//  commands["left-wheel-joint"] = Vector3(q_des[0], qd[0], qdd_des[0]);
//  commands["right-wheel-joint"] = Vector3(q_des[1], qd[1], qdd_des[1]);
  commands["left-wheel-joint"] = Vector3(q_des[0], qd[0], 1);
  commands["right-wheel-joint"] = Vector3(q_des[1], qd[1], 1);
  return commands;
}

/// Does computed torque control
void control_PID(RCArticulatedBodyPtr robot, const std::map<std::string, Vector3>& commands)
{
  // ********************************************************************
  // compute inverse dynamics
  // ********************************************************************
  std::map<RigidBodyPtr, RCArticulatedBodyInvDynData> inv_dyn_data;

  // get the commands
  Vector2 q_des, qd_des, qdd_des;
  q_des[0] = commands.find("left-wheel-joint")->second[0];
  qd_des[0] = commands.find("left-wheel-joint")->second[1];
  qdd_des[0] = commands.find("left-wheel-joint")->second[2];
  q_des[1] = commands.find("right-wheel-joint")->second[0];
  qd_des[1] = commands.find("right-wheel-joint")->second[1];
  qdd_des[1] = commands.find("right-wheel-joint")->second[2];

  // write the desired joint positions, velocities, and acceleration
  std::ofstream out;
  out.open("left-wheel-joint.out", std::ios::app);
  assert(!out.fail());
  out << q_des[0] << " " << qd_des[0] << " " << qdd_des[0] << std::endl;
  out.close();
  out.open("right-wheel-joint.out", std::ios::app);
  assert(!out.fail());
  out << q_des[1] << " " << qd_des[1] << " " << qdd_des[1] << std::endl;
  out.close();

  // get all links
  RigidBodyPtr base = robot->get_links().front();
  RigidBodyPtr left_wheel = robot->find_link("wheel-left");
  RigidBodyPtr right_wheel = robot->find_link("wheel-right");

  // setup inverse dynamics data for the left wheel
  RCArticulatedBodyInvDynData lw_data;
  lw_data._fext = left_wheel->sum_forces(); 
  lw_data._text = left_wheel->sum_torques(); 
  lw_data._qdd = VectorN(1);
  lw_data._qdd[0] = qdd_des[0];
  inv_dyn_data[left_wheel] = lw_data;

  // setup inverse dynamics data for the right wheel
  RCArticulatedBodyInvDynData rw_data;
  rw_data._fext = right_wheel->sum_forces(); 
  rw_data._text = right_wheel->sum_torques(); 
  rw_data._qdd = VectorN(1);
  rw_data._qdd[0] = qdd_des[0];
  inv_dyn_data[right_wheel] = rw_data;

  // setup inverse dynamics data for the base
  RCArticulatedBodyInvDynData base_data;
  base_data._fext = base->sum_forces();
  base_data._text = base->sum_torques();
  inv_dyn_data[base] = base_data;

  // compute inverse dynamics
  RNEAlgorithm rne;
  std::map<JointPtr, VectorN> actuator_forces = rne.calc_inv_dyn(robot, inv_dyn_data);

  // clear and set motor torques
  for (std::map<JointPtr, VectorN>::const_iterator i = actuator_forces.begin(); i != actuator_forces.end(); i++)
  {
    i->first->reset_force();
    i->first->add_force(i->second);
  }

  // do PID control    
  Real KP = 0, KV = 100;
  VectorN force(1);  

  // left wheel joint
  JointPtr lw_joint = robot->find_joint("left-wheel-joint");
  Real lq_t = lw_joint->q[0];
  Real lqd_t = lw_joint->qd[0];  
  Real lP = KP*(q_des[0] - lq_t);
  Real lV = KV*(qd_des[0] - lqd_t);
  force[0] = lP + lV;
  out.open("left-wheel-joint.torques.out", std::ofstream::app);
  out << lw_joint->get_force()[0] << " " << lP << " " << lV << " 0 " << (lw_joint->get_force()[0] + force[0]) << std::endl;
  out.close();
force[0] = 1;
  lw_joint->add_force(force);
  out.open("left", std::ios::app);
  out << force[0] << std::endl;
  out.close();

  // dump all of the true q's to a file
  out.open("right-wheel-joint.true.out", std::ofstream::app);
  out << lq_t << " " << lqd_t << " " << lw_joint->qdd[0] << std::endl;
  out.close();

  // right wheel joint
  JointPtr rw_joint = robot->find_joint("right-wheel-joint");
  Real rq_t = rw_joint->q[0];
  Real rqd_t = rw_joint->qd[0];  
  Real rP = KP*(q_des[1] - rq_t);
  Real rV = KV*(qd_des[1] - rqd_t);
  force[0] = rP + rV;
  out.open("right-wheel-joint.torques.out", std::ofstream::app);
  out << rw_joint->get_force()[0] << " " << rP << " " << rV << " 0 " << (rw_joint->get_force()[0] + force[0]) << std::endl;
  out.close();
force[0] = 1;
  rw_joint->add_force(force);
  out.open("right", std::ios::app);
  out << force[0] << std::endl;
  out.close();

  // dump all of the true q's to a file
  out.open("left-wheel-joint.true.out", std::ofstream::app);
  out << rq_t << " " << rqd_t << " " << rw_joint->qdd[0] << std::endl;
  out.close();
}
*/

/// Controller
void controller(DynamicBodyPtr pioneer, Real time, void* data)
{
  const unsigned X = 0, Z = 2;
  const Real DT = 1e-6;
  const Real PATH_FREQ = 1;
  static unsigned last_draw = 0;
  const Real NEAR_ZERO = std::sqrt(std::numeric_limits<Real>::epsilon());

  // get the pioneer as an articulated body
  RCArticulatedBodyPtr robot = dynamic_pointer_cast<RCArticulatedBody>(pioneer);

  // get the current base position and orientation for the robot
  RigidBodyPtr base = robot->get_links().front();
  Real current_x = base->get_position()[X];
  Real current_z = base->get_position()[Z];


  // uncomment following block for forward motion only
/*
  Real x = 0, xd = 0, theta = 0, thetad = 0;
  Real z = ZZ;
  Real zd = ZZD;  

  // control base
  controlBaseDiff(robot, x, xd, z, zd, theta, thetad);
*/

  // uncomment following block for angular motion only
/*
  Real x = 0, xd = 0, z = 0, zd = 0;
  Real theta = ZTHETA;
  Real thetad = ZTHETAD;

  // control base
  controlBaseDiff(robot, x, xd, z, zd, theta, thetad);
*/

  // uncomment following block for sinusoidal motion
///*
  // determine next position, orientation for the robot
  // NOTE: we use 3/2 * pi to start from (0,0)
  // Set the time scale smaller to get more accurate control
  const Real TIME_SCALE = 1.0/10;
  const Real THETA_ERR_SCALE = 0.01;
  Real t = time*TIME_SCALE + 3.0/2.0*M_PI;
  Real x, x0, z, z0, xd, zd, xd0, zd0, theta, thetad;

  // get the desired position of the robot now and shortly after now
  fig8(t, x, xd, z, zd);
  fig8(t+DT, x0, xd0, z0, zd0);

  // determine the error in the robot's orientation
  Real theta_err = calc_theta(x, z, current_x, current_z);

  // determine theta
  Vector2 dir(x0 - x, z0 - z);
  if (dir.norm() > NEAR_ZERO)
  {
    // get the unit vector pointing to the direction of travel
    dir.normalize();
    theta = std::atan2(-dir[1], dir[0]) - M_PI_2;
    if (theta < -M_PI)
      theta += 2 * M_PI;

    // get the next point in the trajectory
    Real x1, z1, xd1, zd1;
    fig8(t+2.0*DT, x1, xd1, z1, zd1);

    // get the unit vector for the next direction of travel
    Vector2 dir2(x1 - x0, z1 - z0);
    if (dir2.norm() > NEAR_ZERO)
    {
      dir2.normalize();
      Real thetah = std::atan2(-dir2[1], dir2[0]) - M_PI_2;
      if (thetah < -M_PI)
        thetah += 2 * M_PI;
      thetad = angle_diff(theta, thetah)/DT;
    }
    else
      thetad = 0.0;
  }
  else
  {
    Matrix3 R;
    base->get_transform().get_rotation(&R);
    Real base_theta = get_y_axis_angle(R);
    theta = base_theta;
    thetad = 0.0;
  }

  // add in the error to theta
  theta += theta_err * THETA_ERR_SCALE;

  // control base
  controlBaseDiff(robot, x, (xd0-xd)/DT, z, (zd0-zd)/DT, theta, thetad);
//*/


  // visualize desired path
  if (time > last_draw*PATH_FREQ)
  {
    last_draw++;

    #ifdef USE_INVENTOR
    SoSeparator* dpath_sep = new SoSeparator;
    SoTransform* dpath_trans = new SoTransform;
    dpath_trans->translation = SbVec3f(x, .05, z);
    dpath_sep->addChild(dpath_trans);
    dpath_sep->addChild(circle);
    dpath->addChild(dpath_sep);

    // visualize robot path
    SoSeparator* rpath_sep = new SoSeparator;
    SoTransform* rpath_trans = new SoTransform;
    rpath_trans->translation = SbVec3f(current_x, .05, current_z);
    rpath_sep->addChild(rpath_trans);
    rpath_sep->addChild(circle);
    rpath->addChild(rpath_sep);
    #endif
  }

  // update theta/z desireds
  ZZD += STEP_SIZE*ZZDD*0.1;
  ZZ += STEP_SIZE*ZZD*0.1;
  ZTHETAD += STEP_SIZE*ZTHETADD*0.1;
  ZTHETA += STEP_SIZE*ZTHETAD*0.1;
  if (ZTHETA > M_PI)
    ZTHETA = -2*M_PI + ZTHETA;
}

void init(void* main_separator, const std::map<std::string, BasePtr>& read_map, Real time)
{
  // find the robot
  if (read_map.find("pioneer2") == read_map.end())
    throw std::runtime_error("controller.cpp:init() - unable to find pioneer2!");
  DynamicBodyPtr robot = dynamic_pointer_cast<DynamicBody>(read_map.find("pioneer2")->second);
  if (!robot)
    throw std::runtime_error("controller.cpp:init() - unable to cast pioneer2 to type DynamicBody");

  // setup the controller
  robot->controller = controller;

  #ifdef USE_INVENTOR
  // draw disks 
  const Real RADIUS = 0.025;
  const Real HEIGHT = 0.003125;

  // get the main separator
  SoSeparator* main_sep = (SoSeparator*) main_separator;

  // create a circle
  circle = new SoCylinder;
  circle->radius = RADIUS;
  circle->height = HEIGHT;

  // init separators
  rpath = new SoSeparator;
  dpath = new SoSeparator;

  // add a complexity node to cut down on rendering time
  SoComplexity* complexity = new SoComplexity;
  rpath->addChild(complexity);
  dpath->addChild(complexity);

  // create materials
  SoMaterial* rmat = new SoMaterial();
  rmat->diffuseColor = SbColor(0,0,1);
  SoMaterial* dmat = new SoMaterial();
  dmat->diffuseColor = SbColor(0,1,0);
  rpath->addChild(rmat);
  dpath->addChild(dmat);
  
  // add the separators to the main separator  
  main_sep->addChild(rpath);
  main_sep->addChild(dpath);
  #endif
}

} // end extern "C"

