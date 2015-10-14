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

const double ALPHA = 0.0702;

Ravelin::Origin3d quat2rpy(const Ravelin::Quatd& q_){

  Ravelin::Origin3d rpy;
  Ravelin::VectorNd q(4);
  q[0] = q_.w;
  q[1] = q_.x;
  q[2] = q_.y;
  q[3] = q_.z;

  // Singularity Condition 2(q1q3 + q0q2) == +/- 1
  assert(fabs(2*(q[1]*q[3] + q[0]*q[2])) < (1.0 - NEAR_ZERO) || fabs(2*(q[1]*q[3] + q[0]*q[2])) > (1.0 + NEAR_ZERO));

  // (3-2-1) z-y-x Tait-Bryan rotation
  rpy[0] = atan2((q[2]*q[3] + q[0]*q[1]),0.5-(q[1]*q[1] + q[2]*q[2]));
  rpy[1] = asin(2.0*(-q[1]*q[3] + q[0]*q[2]));
  rpy[2] = atan2((q[1]*q[2] + q[0]*q[3]),0.5-(q[2]*q[2] + q[3]*q[3]));
  return rpy;
}

// setup simulator callback
void post_step_callback(Simulator* s)
{
  const unsigned Z = 2;
  std::ofstream out("energy.dat", std::ostream::app);
  double KE = walker->calc_kinetic_energy();
  shared_ptr<RigidBodyd> base = walker->get_links().front();
  shared_ptr<RigidBodyd> l1 = walker->get_links().back();
  Transform3d gTb = Pose3d::calc_relative_pose(base->get_inertial_pose(), GLOBAL);
  Transform3d gTl1 = Pose3d::calc_relative_pose(l1->get_inertial_pose(), GLOBAL);
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
}

Ravelin::Quatd rpy2quat(const Ravelin::Origin3d& rpy){

  const double PHI = rpy[0] * (double) 0.5;
  const double THE = rpy[1] * (double) 0.5;
  const double PSI = rpy[2] * (double) 0.5;

  // precompute trig fns
  const double CPHI = std::cos(PHI);
  const double SPHI = std::sin(PHI);
  const double CPSI = std::cos(PSI);
  const double SPSI = std::sin(PSI);
  const double CTHE = std::cos(THE);
  const double STHE = std::sin(THE);

  // construct Quaternion
  Ravelin::Quatd q;
  q.w = (CPHI * CTHE * CPSI + SPHI * STHE * SPSI);
  q.x = (SPHI * CTHE * CPSI - CPHI * STHE * SPSI);
  q.y = (CPHI * STHE * CPSI + SPHI * CTHE * SPSI);
  q.z = (CPHI * CTHE * SPSI - SPHI * STHE * CPSI);
  return q;
}

// called when the simulator completes a step

void post_event_callback_fn(const std::vector<UnilateralConstraint>& e,
                            boost::shared_ptr<void> empty)
{
  std::cout << "> > start post_event_callback_fn(.)" << std::endl;

  // PROCESS CONTACTS
  for(unsigned i=0;i<e.size();i++){
    if (e[i].constraint_type == UnilateralConstraint::eContact)
    {
      shared_ptr<SingleBodyd> sb1 = e[i].contact_geom1->get_single_body();
      shared_ptr<SingleBodyd> sb2 = e[i].contact_geom2->get_single_body();

      std::cout << "contact: " << sb1->body_id << " and " << sb2->body_id << std::endl;
      std::cout << "i = " << e[i].contact_impulse.get_linear() << std::endl;
      std::cout << "p = " << e[i].contact_point << std::endl;
      std::cout << "n = " << e[i].contact_normal << std::endl;
//      std::cout << "s = " << e[i].contact_tan1 << std::endl;
//      std::cout << "t = " << e[i].contact_tan2 << std::endl;
//      std::cout << "muC = " << e[i].contact_mu_coulomb << std::endl;
//      std::cout << "muV = " << e[i].contact_mu_viscous << std::endl;
    }
  }
  std::cout << "<< end post_event_callback_fn(.)" << std::endl;
}

void controller_callback(ControlledBodyPtr dbp, double t, void*)
{
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
    boost::shared_ptr<const Ravelin::Pose3d> Ipose = links[i]->get_inertial_pose();
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


//  walker->controller                  = &controller_callback;
//  sim->constraint_post_callback_fn  = &post_event_callback_fn;

  // Set initial conditions from ruina paper
  Ravelin::VectorNd x,xd;
  walker->get_generalized_coordinates_euler(x);
  walker->get_generalized_velocity( DynamicBodyd::eSpatial,xd);

  x[1] = 0; // x
  x[2] = 0; // y
// x[2] = 0.6969; // y
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

/*
  Ravelin::Matrix3d Rz;
  Rz(0,0) = std::cos(PHI);  Rz(0,1) = std::sin(PHI);  Rz(0,2) = 0.0;
  Rz(1,0) = -std::sin(PHI); Rz(1,1) = std::cos(PHI);  Rz(1,2) = 0.0;
  Rz(2,0) = 0.0;            Rz(2,1) = 0.0;            Rz(2,2) = 1.0;

  Ravelin::Matrix3d Rx;
  Rx(0,0) = 1.0;    Rx(0,1) = 0.0;            Rx(0,2) = 0.0;
  Rx(1,0) = 0.0;    Rx(1,1) = std::cos(PSI);  Rx(1,2) = std::sin(PSI);
  Rx(2,0) = 0.0;    Rx(2,1) = -std::sin(PSI); Rx(2,2) = std::cos(PSI);

  Ravelin::Matrix3d Ry;
  Ry(0,0) = std::cos(THE);  Ry(0,1) = 0.0;   Ry(0,2) = -std::sin(THE);
  Ry(1,0) = 0.0;            Ry(1,1) = 1.0;   Ry(1,2) = 0.0;
  Ry(2,0) = std::sin(THE);  Ry(1,2) = 0.0;   Ry(2,2) = std::cos(THE);
  Ravelin::Matrix3d Rzxy = Ry * Rx * Rz;
*/

  // Note: each rotation matrix defined in Michael Coleman's Matlab code
  // is a rotation about the negative axis
  Ravelin::Matrix3d Rz = Ravelin::Matrix3d::rot_Z(-PHI),
                    Ry = Ravelin::Matrix3d::rot_Y(-THE),
                    Rx = Ravelin::Matrix3d::rot_X(-PSI),
                    Rzxy = Ry.mult(Rx).mult(Rz);
  Quatd q = Quatd(Rzxy);
   
  // output the rotation
  Matrix3d R = q;
  std::cout << "initial rotation: " << std::endl << R;
  const double THETA_OFFSET = M_PI;
  //printf ("q0:%5.5f q1:%5.5f q2:%5.5f q3:%5.5f\n",q[0],q[1],q[2],q[3]);

  //x[4] = q[0];
  //x[5] = q[1];
  //x[6] = q[2];
  //x[7] = q[3];
  x[4] = 0;
  x[5] = 0;
  x[6] = 0.0871557427476582;
  x[7] = 0.9961946980917456;
  x[0] =  M_PI + THETA_OFFSET; // Theta_sw
  //x[4]=0;
  //x[5]=0;
  //x[6]=0;
  //x[7]=1;
  //x[0]=0;    

  // Convert Time derivative of RPY to angular velocity (w)
  double  dPHI = -1.322096551035500e-001, // d yaw / dt
          dPSI = -1.990961987794000e-002, // d roll / dt
          dTHE =  4.712423746697700e-001; // d pitch / dt
  //  a_dot is time derivative of RPY
  // Numerically derive skew symmetric
  // angular velocity tensor matrix W
  // can do this analytically also
  double h = 1e-6;
  Matrix3d
      Rz2 = Ravelin::Matrix3d::rot_Z(-PHI - h*dPHI),
      Rx2 = Ravelin::Matrix3d::rot_X(-PSI - h*dPSI),
      Ry2 = Ravelin::Matrix3d::rot_Y(-THE - h*dTHE);
   Matrix3d
      Rzxy2 = Ry2.mult(Rx2).mult(Rz2);
   Matrix3d
      W = ((Rzxy2-Rzxy)/h).mult_transpose(Rzxy);

  // w is angular velocity
  //Vector3d w((W(2,1)-W(1,2))/2,(W(0,2)-W(2,0))/2,(W(1,0)-W(0,1))/2);
  printf ("w1:%5.5f w2:%5.5f w3:%5.5f\n",(W(2,1)-W(1,2))/2,(W(0,2)-W(2,0))/2,(W(1,0)-W(0,1))/2);
  Vector3d w(0,0,0);
  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  //w[0]=0;
  //w[1]=0;
  //w[2]=0;
  xd.set_zero();
  xd.set_sub_vec(4,w);
  xd[0] = 0; // -Theta_sw
  //xd[0]=0;
  // get the collision geometries for the feet
  lfoot_geom = dynamic_pointer_cast<RigidBody>(walker->find_link("LLEG"))->geometries.front();
  rfoot_geom = dynamic_pointer_cast<RigidBody>(walker->find_link("RLEG"))->geometries.front();

  walker->set_generalized_coordinates_euler(x);
  walker->set_generalized_velocity( DynamicBodyd::eSpatial,xd);

  walker->set_floating_base(false);
//  walker->get_base_link()->set_enabled(false);

}
} // end extern C
