/*****************************************************************************
 * Controller for LINKS robot
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
boost::shared_ptr<Moby::EventDrivenSimulator> sim;

using namespace Ravelin;

Ravelin::Vector3d quat2rpy(const Ravelin::Quatd& q_){

  Ravelin::Vector3d rpy;
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

Ravelin::Quatd rpy2quat(const Ravelin::Vector3d& rpy){

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
  part->get_generalized_coordinates( Moby::DynamicBody::eSpatial,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  std::cout << "Time = " << t << std::endl;
  std::cout << "x = " << x << std::endl;
  std::cout << "v = " << xd << std::endl;
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

  part->controller                  = &controller_callback;
  sim->constraint_post_callback_fn  = &post_event_callback_fn;

  // Set initial conditions from ruina paper

  Ravelin::VectorNd x,xd;
  part->get_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  x[1] =  0; // x
  x[2] = 0; // y
  x[3] = 1.1232; // z

  double  PHI = 0.09866765986740,  // yaw
          THE = -0.16016583495522, // pitch
          PSI = -0.00924861067616; // roll
  Vector3d a(PSI,THE,PHI);
  Quatd q = rpy2quat(a);

  x[4] = q[0];
  x[5] = q[1];
  x[6] = q[2];
  x[7] = q[3];
  x[0] =  3.43583389038583; // Theta_sw

  // Convert Time derivative of RPY to angular velocity (w)
  double  dPHI = -0.13220965510356, // d yaw / dt
          dTHE =  0.47124237466979, // d pitch / dt
          dPSI = -0.01990961987794; // d roll / dt
  //  a_dot is time derivative of RPY
  Vector3d a_dot(dPSI,dTHE,dPHI);

  // Numerically derive skew symmetric
  // angular velocity tensor matrix W
  double h = 1e-4;
  Matrix3d
      R1(rpy2quat(a)),
      R2(rpy2quat(a+h*a_dot)),
      W = ((R2-R1)/h).mult_transpose(R1);

  // w is angular velocity
  Vector3d w(W(2,1),W(0,2),W(1,0));

  xd[4] = w[0]; // phi = roll
  xd[5] = w[1]; // Theta_st = pitch
  xd[6] = w[2]; // psi = yaw
  xd[0] = -0.39255916866486; // Theta_sw

  part->set_generalized_coordinates( Moby::DynamicBody::eEuler,x);
  part->set_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

}
} // end extern C
