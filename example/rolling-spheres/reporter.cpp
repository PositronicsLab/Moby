/*****************************************************************************
 * Reporting function for rolling sphere 
 ****************************************************************************/

// to run moby-driver -r -p=libreporter.so -oi -s=0.01 rolling-spheres-polyhedral.xml
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/GravityForce.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

boost::shared_ptr<TimeSteppingSimulator> sim;
RigidBodyPtr sphere;

// setup simulator callback
void post_step_callback(Simulator* s)
{
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = sphere->calc_kinetic_energy(sphere->get_pose());
  out << KE << std::endl;
  out.close();
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, BasePtr>& read_map, double time)
{
  // overwrite files
  std::ofstream out("energy.dat");
  out.close();

  for (std::map<std::string, BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);

    // find the robot reference
    if (!sphere)
    {
      sphere = boost::dynamic_pointer_cast<RigidBody>(i->second);
      if (sphere && !sphere->is_enabled())
        sphere.reset();
    }
  }

  // setup the callback
  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
