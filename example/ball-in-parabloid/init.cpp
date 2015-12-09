/*****************************************************************************
 * Initializer for ball/parabloid example 
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
using namespace Ravelin;
using namespace Moby;

Moby::RigidBodyPtr parabloid;
Moby::RigidBodyPtr ball;
boost::shared_ptr<TimeSteppingSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup simulator callback
void post_step_callback(Simulator* sim)
{
  const unsigned Z = 2;

  // output the energy of the ball
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = ball->calc_kinetic_energy();
  Transform3d gTw = Pose3d::calc_relative_pose(ball->get_pose(), GLOBAL);
  double PE = ball->get_inertia().m*gTw.x[Z]*-grav->gravity[Z];
  out << KE << " " << PE << " " << (KE+PE) << std::endl;
  out.close();

  // output the position of the ball
  out.open("position.dat", std::ostream::app);
  out << gTw.x << std::endl;
  out.close();

  // see whether we are in a ballistic flight phase
  TimeSteppingSimulator* esim = (TimeSteppingSimulator*) sim;
  boost::shared_ptr<CollisionDetection> coldet = esim->get_collision_detection(); 
  CollisionGeometryPtr cgw = ball->geometries.front();
  CollisionGeometryPtr cgg = parabloid->geometries.front();
  Point3d cpw, cpg;
  double dist = coldet->calc_signed_dist(cgw, cgg, cpw, cpg);
  if (dist > 1e-4)
    std::cerr << "-- in a ballistic flight phase (dist=" << dist << ") at time " << sim->current_time << std::endl;
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  const unsigned Z = 2;

  // overwrite the energy and position files
  std::ofstream out("energy.dat");
  out.close();
  out.open("position.dat");
  out.close();

  // If use robot is active also init dynamixel controllers
  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);
    if (i->first == "ball")
      ball = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "parabloid")
      parabloid = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
