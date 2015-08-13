/*****************************************************************************
 * "Controller" for sliding/rolling sphere example 
 ****************************************************************************/
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Moby/ContactParameters.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>
#include <stdlib.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr sphere;
Moby::RigidBodyPtr ground;
boost::shared_ptr<TimeSteppingSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup simulator callback
void post_step_callback(Simulator* s)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the sphere radius
  const double R = 1.0;

  // get the bottom of the sphere
  Transform3d wTs = Pose3d::calc_relative_pose(sphere->get_pose(), GLOBAL);

  shared_ptr<Pose3d> Pbot(new Pose3d);  
  Pbot->rpose = GLOBAL;
  Pbot->x = wTs.x;
  Pbot->x[Z] -= R;

  // get the velocity of the sphere at the contact point
  SVelocityd v = sphere->get_velocity();
  Transform3d botTv = Pose3d::calc_relative_pose(v.pose, Pbot);
  SVelocityd xd = botTv.transform(v);
  Vector3d linear = xd.get_linear();

/*
  SVelocityd v = sphere->get_velocity();
  Origin3d xd(v.get_linear());
  Origin3d omega(v.get_angular());
  Origin3d s(1.0, 0.0, 0.0);
  Origin3d t(0.0, 1.0, 0.0);
  Origin3d crosss = Origin3d::cross(-wTs.x, s);
  Origin3d crosst = Origin3d::cross(-wTs.x, t);
*/

  // output the sliding velocity at the contact 
  std::ofstream out("contactv.dat", std::ostream::app);
  out << sim->current_time << " " << linear[X] << " " << linear[Y] << " " << linear[Z] << std::endl;
//  out << sim->current_time << " " << (s.dot(xd) + crosss.dot(omega)) << " " << (t.dot(xd) + crosst.dot(omega)) << std::endl; 
//  out << sim->current_time << " " << v[3] << " " << v[4] << " " << v[5] << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
  out.close();

  out.open("velocity.dat", std::ostream::app);
  out << sim->current_time << " " << v[3] << " " << v[4] << " " << v[5] << " " << v[0] << " " << v[1] << " " << v[2] << std::endl; 
  out.close();

  out.open("ke.dat", std::ostream::app);
  out << sim->current_time << " " << sphere->calc_kinetic_energy() << std::endl;
  out.close();
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  const unsigned Z = 2;

  // wipe out contactv
  std::ofstream out("contactv.dat");
  out.close();
  out.open("ke.dat");
  out.close();

  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);
    if (i->first == "sphere")
      sphere = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "ground")
      ground = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
