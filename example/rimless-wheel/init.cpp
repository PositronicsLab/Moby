/*****************************************************************************
 * Initializer for rimless wheel 
 ****************************************************************************/
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>
#include <stdlib.h>
#include "params.h"

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr ground;
Moby::RigidBodyPtr wheel;
boost::shared_ptr<TimeSteppingSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup simulator callback
void post_step_callback(Simulator* sim)
{
  const unsigned Z = 2;
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = wheel->calc_kinetic_energy();
  Transform3d gTw = Pose3d::calc_relative_pose(wheel->get_pose(), GLOBAL);
  double PE = wheel->get_inertia().m*gTw.x[Z]*-grav->gravity[Z];
  out << sim->current_time << " " << KE << " " << PE << " " << (KE+PE) << std::endl;
  out.close();

  // see whether there is significant undesired rotation
  AAngled aa = Pose3d::calc_relative_pose(wheel->get_pose(), GLOBAL).q;
  out.open("angular.dat", std::ostream::app);
  out << sim->current_time << " " << std::fabs(aa.x) << " " << std::fabs(aa.z) << " " << " " << std::fabs(aa.angle) << std::endl; 
  out.close(); 

  // see whether we are in a ballistic flight phase
  TimeSteppingSimulator* esim = (TimeSteppingSimulator*) sim;
  boost::shared_ptr<CollisionDetection> coldet = esim->get_collision_detection(); 
  CollisionGeometryPtr cgw = wheel->geometries.front();
  CollisionGeometryPtr cgg = ground->geometries.front();
  Point3d cpw, cpg;
  double dist = coldet->calc_signed_dist(cgw, cgg, cpw, cpg);
  if (dist > 1e-4)
    std::cerr << "-- in a ballistic flight phase at time " << sim->current_time << std::endl;

  // fast exit conditions
  if (FIND_MAP)
  {
    if (KE < NEAR_ZERO)
    {
      std::cerr << "kinetic energy too small!" << std::endl;
      out.open("system.state", std::ostream::app);
      out << "0.0" << std::endl;
      exit(0);
    }
    if (sim->current_time > 100.0)
    {
      std::cerr << "simulation ran too long!" << std::endl;
      out.open("system.state", std::ostream::app);
      out << "DnF" << std::endl;
      exit(0);
    }
  }

  // if we're finding the return map, see whether the next step has been
  // encountered 
  if (FIND_MAP)
  {
    // look to see whether the last processed contact was a new spoke
    std::ifstream in("IPC.token");
    if (in.fail())
      return;
    in.close();

    // it was, output the state of the system
    std::ofstream out("system.state", std::ostream::app);
    out << " " << wheel->get_velocity().get_angular()[1] << std::endl;
    out.close();
    exit(0);
  }
}

// setup simulator mini-callback
void post_ministep_callback(ConstraintSimulator* sim)
{
  const unsigned Y = 1;

  // see whether there is significant undesired rotation
  AAngled aa = Pose3d::calc_relative_pose(wheel->get_pose(), GLOBAL).q;
  if ((std::fabs(aa.x) > 1e-4 || std::fabs(aa.z) > 1e-4) && std::fabs(aa.angle) > 1e-6)
  {
    std::cerr << "excessive angular rotation out of the plane: " << std::max(std::fabs(aa.x), std::fabs(aa.z)) << std::endl;
  }

  // output the velocity
  std::ofstream out("velocity.dat", std::ostream::app);
  out << sim->current_time << " " << wheel->get_velocity().get_angular()[Y] << std::endl;
  out.close();

}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  const unsigned Z = 2;

  // overwrite the energy and velocity files
  std::ofstream out("energy.dat");
  out.close();
  out.open("velocity.dat");
  out.close();
  out.open("cvio.dat");
  out.close();

  // If use robot is active also init dynamixel controllers
  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);
    if (i->first == "WHEEL")
      wheel = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "GROUND")
      ground = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;
  sim->post_mini_step_callback_fn = &post_ministep_callback;

  // initialize the system to the fixed point
//  double theta = M_PI/N_SPOKES;
  double theta = 0.0;
  const double ALPHA = 0.2;
  const double LAMBDA_SQ = 2.0/3;
  const double MU = 1.0 + LAMBDA_SQ * (std::cos(2*M_PI/N_SPOKES) - 1.0);
  assert(MU > 0.0 && MU < 1.0);
//  double theta_dot = std::sqrt((4.0*MU*MU*LAMBDA_SQ*std::sin(M_PI/N_SPOKES)*std::sin(ALPHA))/(1.0 - MU*MU));

  // get theta dot
  char* theta_dot_str = getenv("RIMLESS_WHEEL_THETAD");
  if (!theta_dot_str) 
  {
    std::cerr << "RIMLESS_WHEEL_THETAD not defined!" << std::endl;
    exit(-1);
  }
  double theta_dot = std::atof(theta_dot_str);
  out.open("system.state", std::ostream::app);
  out << theta_dot << " ";
  out.close();

  // get the distance per revolution
  const double DIST_PER_REV = 2*M_PI*R;
  const double REV_PER_SEC = theta_dot / (M_PI*2.0);
  const double DIST_PER_SEC = DIST_PER_REV * REV_PER_SEC;

  // set the rotation about y
  Quatd q_wheel = Matrix3d::rot_Y(theta);
  Vector3d lvel(DIST_PER_SEC, 0.0, 0.0, wheel->get_velocity().pose); 
  Vector3d avel(0, theta_dot, 0, wheel->get_velocity().pose);
  Pose3d P;
  P.x.set_zero();
  P.x[Z] = 0.866025403784439; 
  P.q = q_wheel;
  wheel->set_pose(P);
  SVelocityd v;
  v.set_angular(avel);
  v.set_linear(lvel);
  wheel->set_velocity(v);
/*
  // Set initial conditions from ruina paper
  Ravelin::VectorNd x,xd;
  part->get_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  x[1] = 0; // x
  x[2] = 0; // y
  x[3] = 0.1236; // z

//     9.866765986740000e-002
//    -9.248610676160000e-003
//    -1.601658349552200e-001
//     3.435833890385830e+000
//    -1.322096551035500e-001
//    -1.990961987794000e-002
//     4.712423746697700e-001
//    -3.925591686648300e-001

  double  PHI = 9.866765986740000e-002,  // yaw
          THE = -1.601658349552200e-001, // pitch + alpha
          PSI = -9.248610676160000e-003; // roll
  Ravelin::Matrix3d Rz = Ravelin::Matrix3d::rot_Z(PHI),
                    Ry = Ravelin::Matrix3d::rot_Y(THE),
                    Rx = Ravelin::Matrix3d::rot_X(PSI),
                    Rzxy = Ry.mult(Rx).mult(Rz);
  Quatd q = Quatd(Rzxy);

  x[4] = q[0];
  x[5] = q[1];
  x[6] = q[2];
  x[7] = q[3];
  x[0] =  3.435833890385830e+000 + THETA_SW_OFFSET; // Theta_sw

  // Convert Time derivative of RPY to angular velocity (w)
  double  dPHI = -1.322096551035500e-001, // d yaw / dt
          dTHE =  4.712423746697700e-001, // d pitch / dt
          dPSI = -1.990961987794000e-002; // d roll / dt
  //  a_dot is time derivative of RPY
  // Numerically derive skew symmetric
  // angular velocity tensor matrix W
  double h = 1e-6;
  Matrix3d
      Rz2 = Ravelin::Matrix3d::rot_Z(PHI + h*dPHI),
      Rx2 = Ravelin::Matrix3d::rot_X(PSI + h*dPSI),
      Ry2 = Ravelin::Matrix3d::rot_Y(THE + h*dTHE);
   Matrix3d
      Rzxy2 = Ry2.mult(Rx2).mult(Rz2);
   Matrix3d
      W = ((Rzxy2-Rzxy)/h).mult_transpose(Rzxy);

  // w is angular velocity
  Vector3d w((W(2,1)-W(1,2))/2,(W(0,2)-W(2,0))/2,(W(1,0)-W(0,1))/2);
  xd.set_zero();
  xd.set_sub_vec(4,w);
  xd[0] = -3.925591686648300e-001; // Theta_sw

  part->set_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->set_generalized_velocity( Moby::DynamicBody::eSpatial,xd);
*/
}
} // end extern C
