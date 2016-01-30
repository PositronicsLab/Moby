/*****************************************************************************
 * "Controller" for sliding to sticking box example 
 ****************************************************************************/
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <fstream>
#include <stdlib.h>

using boost::shared_ptr;
using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr pendulum_link, pendulum_rod_link;
Moby::RCArticulatedBodyPtr pendulum;
boost::shared_ptr<Ravelin::Jointd> first_joint;
boost::shared_ptr<TimeSteppingSimulator> sim;
boost::shared_ptr<GravityForce> grav;

/// plugin must be "extern C"
extern "C" {

void controller(shared_ptr<ControlledBody> cb, double, void*)
{
  RCArticulatedBodyPtr rcab = dynamic_pointer_cast<RCArticulatedBody>(cb);
  for (unsigned i=0; i< rcab->get_links().size(); i++)
    std::cout << "link " << rcab->get_links()[i]->body_id << " pose: " << Pose3d::calc_relative_pose(rcab->get_links()[i]->get_pose(), GLOBAL) << std::endl;
  std::cout << "---" << std::endl;
}

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
    if (i->first == "pendulum")
      pendulum = boost::dynamic_pointer_cast<RCArticulatedBody>(i->second);
    if (i->first == "pendulum_link")
      pendulum_link = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "pendulum_rod_link")
      pendulum_rod_link = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
    if (i->first == "q")
      first_joint = boost::dynamic_pointer_cast<Jointd>(i->second);
  }

  pendulum->controller = &controller;
}
} // end extern C
