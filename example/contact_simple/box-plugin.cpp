/*****************************************************************************
 * "Controller" for sliding to sticking box example 
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

Moby::RigidBodyPtr box;
boost::shared_ptr<EventDrivenSimulator> sim;
boost::shared_ptr<GravityForce> grav;

// setup simulator callback
void post_step_callback(Simulator* sim)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the box height 
  const double H = 1.0;

  // get the bottoms of the box 
  Transform3d wTs = Pose3d::calc_relative_pose(box->get_pose(), GLOBAL);
  Vector3d p1 = wTs.transform_point(Vector3d(-.5, -.5, -.5, box->get_pose()));
  Vector3d p2 = wTs.transform_point(Vector3d(-.5, -.5, .5, box->get_pose()));
  Vector3d p3 = wTs.transform_point(Vector3d(.5, -.5, -.5, box->get_pose()));
  Vector3d p4 = wTs.transform_point(Vector3d(.5, -.5, .5, box->get_pose()));

  // get the bottoms of the box in the global frame
  shared_ptr<Pose3d> P1(new Pose3d), P2(new Pose3d), P3(new Pose3d), P4(new Pose3d);  
  P1->x = Origin3d(p1);
  P2->x = Origin3d(p2);
  P3->x = Origin3d(p3);
  P4->x = Origin3d(p4);

  // get the velocity of the box at the contact points
  SVelocityd v = box->get_velocity();
  Vector3d xd1 = Pose3d::calc_relative_pose(v.pose, P1).transform(v).get_linear();
  Vector3d xd2 = Pose3d::calc_relative_pose(v.pose, P2).transform(v).get_linear();
  Vector3d xd3 = Pose3d::calc_relative_pose(v.pose, P3).transform(v).get_linear();
  Vector3d xd4 = Pose3d::calc_relative_pose(v.pose, P4).transform(v).get_linear();

/*
  SVelocityd v = box->get_velocity();
  Origin3d xd(v.get_linear());
  Origin3d omega(v.get_angular());
  Origin3d s(1.0, 0.0, 0.0);
  Origin3d t(0.0, 1.0, 0.0);
  Origin3d crosss = Origin3d::cross(-wTs.x, s);
  Origin3d crosst = Origin3d::cross(-wTs.x, t);
*/

  // output the sliding velocity at the contact 
  std::ofstream out("contactv.dat", std::ostream::app);
  out << sim->current_time << " " << xd1[X] << " " << xd1[Y] << " " << xd2[X] << " " << xd2[Y] << " " << xd3[X] << " " << xd3[Y] << " " << xd4[X] << " " << xd4[Y] << std::endl;
//  out << sim->current_time << " " << (s.dot(xd) + crosss.dot(omega)) << " " << (t.dot(xd) + crosst.dot(omega)) << std::endl; 
//  out << sim->current_time << " " << v[3] << " " << v[4] << " " << v[5] << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
  out.close();

  out.open("ke.dat", std::ostream::app);
  out << sim->current_time << " " << box->calc_kinetic_energy() << std::endl;
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

  // get a reference to the EventDrivenSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<EventDrivenSimulator>(i->second);
    if (i->first == "box")
      box = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);
  }

  sim->post_step_callback_fn = &post_step_callback;
}
} // end extern C
