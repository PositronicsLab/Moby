/*****************************************************************************
 * Initializer for passive dynamic walker 
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
boost::shared_ptr<Moby::EventDrivenSimulator> sim;

const double ALPHA = 0.0702,
             THETA_SW_OFFSET = -M_PI;

using namespace Ravelin;

Ravelin::Origin3d quat2rpy(const Ravelin::Quatd& q_){

  Ravelin::Origin3d rpy;
  Ravelin::VectorNd q(4);
  q[0] = q_.w;
  q[1] = q_.x;
  q[2] = q_.y;
  q[3] = q_.z;

  // Singularity Condition 2(q1q3 + q0q2) == +/- 1
  assert(fabs(2*(q[1]*q[3] + q[0]*q[2])) < (1.0 - Moby::NEAR_ZERO) || fabs(2*(q[1]*q[3] + q[0]*q[2])) > (1.0 + Moby::NEAR_ZERO));

  // (3-2-1) z-y-x Tait-Bryan rotation
  rpy[0] = atan2((q[2]*q[3] + q[0]*q[1]),0.5-(q[1]*q[1] + q[2]*q[2]));
  rpy[1] = asin(2.0*(-q[1]*q[3] + q[0]*q[2]));
  rpy[2] = atan2((q[1]*q[2] + q[0]*q[3]),0.5-(q[2]*q[2] + q[3]*q[3]));
  return rpy;
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

void post_event_callback_fn(const std::vector<Moby::UnilateralConstraint>& e,
                            boost::shared_ptr<void> empty)
{
  std::cout << ">> start post_event_callback_fn(.)" << std::endl;

  // PROCESS CONTACTS
  for(unsigned i=0;i<e.size();i++){
    if (e[i].constraint_type == Moby::UnilateralConstraint::eContact)
    {
      Moby::SingleBodyPtr sb1 = e[i].contact_geom1->get_single_body();
      Moby::SingleBodyPtr sb2 = e[i].contact_geom2->get_single_body();

      std::cout << "contact: " << sb1->id << " and " << sb2->id << std::endl;
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
void controller_callback(Moby::DynamicBodyPtr dbp, double t, void*)
{
  std::cout << ">> start controller_callback(.)" << std::endl;
  Moby::RCArticulatedBodyPtr
      part = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(dbp);
  Ravelin::VectorNd x,xd;
  static double last_t;
  double h = t-last_t;
  last_t = t;
  part->get_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->get_generalized_velocity( Moby::DynamicBody::eEuler,xd);

  const std::vector<Moby::RigidBodyPtr>& links = part->get_links();
  std::cout << "Time = " << t << std::endl;

  for(int i=0;i<links.size();i++){
    boost::shared_ptr<const Ravelin::Pose3d> Ipose = links[i]->get_inertial_pose();
    boost::shared_ptr<const Ravelin::Pose3d> Lpose = links[i]->get_pose();

    std::cout << links[i]->id << std::endl;
    std::cout << "Ipose x = " << Ravelin::Pose3d::calc_relative_pose(Ipose,Moby::GLOBAL).x << std::endl;
    std::cout << "Lpose x = " << Ravelin::Pose3d::calc_relative_pose(Lpose,Moby::GLOBAL).x << std::endl;
    if(i != 0){
      boost::shared_ptr<const Ravelin::Pose3d> Jpose = links[i]->get_inner_joint_explicit()->get_pose();
      std::cout << "Jpose x = " << Ravelin::Pose3d::calc_relative_pose(Jpose,Moby::GLOBAL).x << std::endl;
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

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  Moby::RCArticulatedBodyPtr part;

  // If use robot is active also init dynamixel controllers
  // get a reference to the EventDrivenSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<Moby::EventDrivenSimulator>(i->second);

    // find the robot reference
    if (!part)
      part = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(i->second);
  }

//  part->controller                  = &controller_callback;
//  sim->constraint_post_callback_fn  = &post_event_callback_fn;

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

  part->set_floating_base(false);
//  part->get_base_link()->set_enabled(false);

}
} // end extern C
