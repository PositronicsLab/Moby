/*****************************************************************************
 * "Controller" for sliding/rolling sphere example 
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>
#include <stdlib.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr sphere;
boost::shared_ptr<EventDrivenSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup simulator callback
void post_step_callback(Simulator* sim)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the sphere radius
  const double R = 1.0;

  // get the bottom of the sphere
  Transform3d wTs = Pose3d::calc_relative_pose(sphere->get_inertial_pose(), GLOBAL);
  shared_ptr<Pose3d> Pbot(new Pose3d);  
  Pbot->rpose = GLOBAL;
  Pbot->x = wTs.x;
  Pbot->x[Z] -= R;

  // get the velocity of the sphere at the contact point
  SVelocityd v = sphere->get_velocity();
  SVelocityd xd = Pose3d::calc_relative_pose(v.pose, Pbot).transform(v);
  Vector3d linear = xd.get_linear();

  // output the sliding velocity at the contact 
  std::ofstream out("contactv.dat", std::ostream::app);
  out << sim->current_time << " " << linear[X] << " " << linear[Y] << std::endl; 
  out.close();
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  const unsigned Z = 2;

  // get a reference to the EventDrivenSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<EventDrivenSimulator>(i->second);
    if (i->first == "sphere")
      sphere = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
