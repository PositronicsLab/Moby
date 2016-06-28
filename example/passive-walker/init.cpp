/*****************************************************************************
 * Initializer for passive dynamic walker 
 ****************************************************************************/

// to run moby-driver -r -p=libpassive-walker-init.so -oi -s=0.001 walker.xml
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
RCArticulatedBodyPtr walker;
shared_ptr<RigidBodyd> ground;
CollisionGeometryPtr lfoot_geom, rfoot_geom;
boost::shared_ptr<GravityForce> grav;

enum TouchdownEvent { eLeftFoot, eRightFoot };

TouchdownEvent touchdown_stop; 

const double ALPHA = 0.0702;

// setup simulator callback
void post_step_callback(Simulator* s)
{
/*
  const unsigned Z = 2;
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = walker->calc_kinetic_energy();
  shared_ptr<RigidBodyd> base = walker->get_links().front();
  shared_ptr<RigidBodyd> l1 = walker->get_links().back();
  Transform3d gTb = Pose3d::calc_relative_pose(base->get_pose(), GLOBAL);
  Transform3d gTl1 = Pose3d::calc_relative_pose(l1->get_pose(), GLOBAL);
  double PEb = base->get_inertia().m*gTb.x[Z]*-grav->gravity[Z];
  double PEl1 = l1->get_inertia().m*gTl1.x[Z]*-grav->gravity[Z];
  out << KE << " " << (PEb+PEl1) << " " << (KE+PEb+PEl1) << std::endl;
  out.close();

  // get the signed distance between the walker feet and the ground
  CollisionGeometryPtr ground_geom = dynamic_pointer_cast<RigidBody>(ground)->geometries.front();
  Point3d dummy1, dummy2;
  shared_ptr<CollisionDetection> coldet = sim->get_collision_detection();
  double dL = coldet->calc_signed_dist(lfoot_geom, ground_geom, dummy1, dummy2);
  double dR = coldet->calc_signed_dist(rfoot_geom, ground_geom, dummy1, dummy2);

  // if the signed distance is negative, project the walker upward by the
  // minimum distance
  if (std::min(dL, dR) < 0.0)
  {
    Pose3d P(*base->get_pose());
    P.x[2] -= std::min(dL, dR);
    base->set_pose(P);
    walker->update_link_poses();
  }
*/
}

// mirrors the walker's configuration
void mirror()
{
  VectorNd q, qd;

  // get the generalized coordinates and velocity
  walker->get_generalized_coordinates_euler(q);
  walker->get_generalized_velocity(DynamicBodyd::eSpatial, qd);

  // base upright, foot forward by 45 degrees should be mirrored to
  // base pitched backward 45 degrees, foot backward by -45 degrees
  const double x = q[0];
  const double xd = qd[0];

  // get the left foot pose 
  shared_ptr<const Pose3d> lf_pose = walker->get_base_link()->get_pose(); 

  // get the right foot pose 
  Pose3d rf_pose = *walker->get_links().back()->get_pose();

  // transform the velocity 
  shared_ptr<const Pose3d> base_vel_pose = walker->get_base_link()->get_velocity().pose;
  SVelocityd v = walker->get_links().back()->get_velocity();
  v.pose = base_vel_pose;
  walker->get_base_link()->set_velocity(v);
 
  // set the left foot pose to the right foot pose
  rf_pose.update_relative_pose(lf_pose->rpose);
  walker->get_base_link()->set_pose(rf_pose);

/*
  // update the angular velocity
  SVelocityd base_xd = base_link->get_velocity();
  base_xd.set_angular(Vector3d(R * Origin3d(base_xd.get_angular()), base_xd.pose));
  base_link->set_velocity(base_xd);
*/
  // get the joint
  shared_ptr<Jointd> joint = walker->get_joints().front();

  // update the joint
  joint->q = -x;
  joint->qd = -xd;

  // update the link transforms
  walker->update_link_poses();
  walker->update_link_velocities();
}

// called when the simulator completes a step
void post_event_callback_fn(const std::vector<Constraint>& e,
                            boost::shared_ptr<void> empty)
{
  VectorNd q, qd;
  std::cout << "> > start post_event_callback_fn(.)" << std::endl;

  // PROCESS CONTACTS
  for(unsigned i=0;i<e.size();i++){
    if (e[i].constraint_type == Constraint::eContact)
    {
      shared_ptr<SingleBodyd> sb1 = e[i].contact_geom1->get_single_body();
      shared_ptr<SingleBodyd> sb2 = e[i].contact_geom2->get_single_body();

      if (e[i].contact_geom1 == rfoot_geom || e[i].contact_geom2 == rfoot_geom)
        std::cout << "right foot contacting; signed distance: " << e[i].signed_distance << std::endl;
      if (e[i].contact_geom1 == lfoot_geom || e[i].contact_geom2 == lfoot_geom)
        std::cout << "left foot contacting; signed distance: " << e[i].signed_distance << std::endl;

      // look for contacting foot matching the stop
      if ((touchdown_stop == eRightFoot && 
          (e[i].contact_geom1 == rfoot_geom || e[i].contact_geom2 == rfoot_geom)) || (touchdown_stop == eLeftFoot &&
          (e[i].contact_geom1 == lfoot_geom || e[i].contact_geom2 == lfoot_geom)))
      {
/*
        // mirror the coordinates
        mirror();
        // get the generalized coordinates and velocity
        walker->get_generalized_coordinates_euler(q);
        walker->get_generalized_velocity(DynamicBodyd::eSpatial, qd);

        // write out the generalized coordinates and velocity
        std::ofstream output("output.dat");
        output << q[0];
        for (unsigned i=1; i< q.size(); i++) 
          output << " " << q[i]; 
        for (unsigned i=0; i< qd.size(); i++) 
          output << " " << qd[i]; 
        output << std::endl;
        output.close();

        // exit with no error message
        exit(0);
*/
      }
    }
  }
  std::cout << "<< end post_event_callback_fn(.)" << std::endl;
}

VectorNd& controller_callback(ControlledBodyPtr dbp, VectorNd& gf, double t, void*)
{
  shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(dbp);
  unsigned ngc = db->num_generalized_coordinates(DynamicBodyd::eSpatial);
  gf.set_zero(ngc);
  return gf;
/*
  std::cout << "> > start controller_callback(.)" << std::endl;
  RCArticulatedBodyPtr
      walker = boost::dynamic_pointer_cast<RCArticulatedBody>(dbp);
  Ravelin::VectorNd x,xd;
  static double last_t;
  double h = t-last_t;
  last_t = t;
  walker->get_generalized_coordinates_euler(x);
  walker->get_generalized_velocity( DynamicBodyd::eEuler,xd);

  const std::vector<shared_ptr<RigidBodyd> >& links = walker->get_links();
  std::cout << "Time = " << t << std::endl;

  for(int i=0;i<links.size();i++){
    boost::shared_ptr<const Ravelin::Pose3d> Ipose = links[i]->get_pose();
    boost::shared_ptr<const Ravelin::Pose3d> Lpose = links[i]->get_pose();

    std::cout << links[i]->body_id << std::endl;
    std::cout << "Ipose x = " << Ravelin::Pose3d::calc_relative_pose(Ipose,GLOBAL).x << std::endl;
    std::cout << "Lpose x = " << Ravelin::Pose3d::calc_relative_pose(Lpose,GLOBAL).x << std::endl;
    if(i != 0){
      boost::shared_ptr<const Ravelin::Pose3d> Jpose = links[i]->get_inner_joint_explicit()->get_pose();
      std::cout << "Jpose x = " << Ravelin::Pose3d::calc_relative_pose(Jpose,GLOBAL).x << std::endl;
    }
  }
  std::cout << "x = " << x << std::endl;
  std::cout << "v = " << xd << std::endl;



//  std::cout << "x =\n\t"
//            << RPY[0] << "\n\t"
//            << RPY[1] << "\n\t"
//            << RPY[2] << "\n\t"
//            << x[0] - THETA_SW_OFFSET << "\n\nxd =\n\t"
//            << dRPY[0] << "\n\t"
//            << dRPY[1] << "\n\t"
//            << dRPY[2] << "\n\t"
//            << xd[0] << "\n\t";
  std::cout << "<< end controller_callback(.)" << std::endl;
*/
}

// ============================================================================
// ================================ CALLBACKS =================================

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, BasePtr>& read_map, double time)
{
  // overwrite files
  std::ofstream out("energy.dat");
  out.close();

  // If use robot is active also init dynamixel controllers
  // get a reference to the TimeSteppingSimulator instance
  for (std::map<std::string, BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(i->second);

    // find the robot reference
    if (!walker)
      walker = boost::dynamic_pointer_cast<RCArticulatedBody>(i->second);

    // find the gravity vector
    if (!grav)
      grav = boost::dynamic_pointer_cast<GravityForce>(i->second);

    // find the ground
    if (!ground && i->first == "GROUND")
      ground = boost::dynamic_pointer_cast<RigidBody>(i->second);
  }

  // setup the callback
  sim->post_step_callback_fn = &post_step_callback;

  // read in the initial configuration
  std::ifstream in("input.dat");
  if (in.fail())
    throw std::runtime_error("No input file located!");
  VectorNd q, qd;
  q.resize(walker->num_generalized_coordinates(DynamicBodyd::eEuler));
  qd.resize(walker->num_generalized_coordinates(DynamicBodyd::eSpatial));
  for (unsigned i=2; i< q.size(); i++)
    in >> q[i]; assert(!in.eof());
  for (unsigned i=0; i< qd.size(); i++)
    in >> qd[i]; assert(!in.eof());

  // set the generalized coordinates and velocity
  walker->set_generalized_coordinates_euler(q);
  walker->set_generalized_velocity(DynamicBodyd::eSpatial, qd);

// mirror the configuration
//mirror();

  // indicate that we stop on left foot touchdown
  touchdown_stop = eRightFoot;

  // setup callbacks
  walker->controller                  = &controller_callback;
//  sim->constraint_post_callback_fn  = &post_event_callback_fn;

  // get the collision geometries for the feet
  lfoot_geom = dynamic_pointer_cast<RigidBody>(walker->find_link("LLEG"))->geometries.front();
  rfoot_geom = dynamic_pointer_cast<RigidBody>(walker->find_link("RLEG"))->geometries.front();
}
} // end extern C
