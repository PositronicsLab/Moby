/*****************************************************************************
 * Controller plugin template 
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

// pointers to objects found within the XML file
boost::shared_ptr<TimeSteppingSimulator> sim;

// post-simulation step callback
void post_step_callback(Simulator* s)
{
}

// control function
void controller(ControlledBodyPtr body, double time, void* controller_data)
{
  
}

/// plugin must be "extern C"
extern "C" {

// initialization function is called by Moby; setup callbacks here 
void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);

    // example of looking for a body and setting up a callback function for it
//    if (i->first == "ground")
//    {
//      ControlledBodyPtr ground = boost::dynamic_pointer_cast<RigidBody>(i->second);
//	assert(ground);
//      ground->controller = &controller;
//    }
  }

  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
