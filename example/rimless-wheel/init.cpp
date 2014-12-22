/*****************************************************************************
 * Initializer for rimless wheel 
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>

using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr wheel;
boost::shared_ptr<EventDrivenSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup controller callback
void post_step_callback(Simulator* sim)
{
  const unsigned Z = 2;
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = wheel->calc_kinetic_energy();
  Transform3d gTw = Pose3d::calc_relative_pose(wheel->get_inertial_pose(), GLOBAL);
  double PE = wheel->get_inertia().m*gTw.x[Z]*-grav->gravity[Z];
  out << KE << " " << PE << " " << (KE+PE) << std::endl;
  out.close();
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{

  // overwrite the energy file
  std::ofstream out("energy.dat");
  out.close();

  // If use robot is active also init dynamixel controllers
  // get a reference to the EventDrivenSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<EventDrivenSimulator>(i->second);
    if (i->first == "WHEEL")
      wheel = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;

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
