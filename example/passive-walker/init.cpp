/*****************************************************************************
 * Controller for LINKS robot
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>

void controller_callback(Moby::DynamicBodyPtr dbp, double t, void*)
{
  Moby::RCArticulatedBodyPtr
      part = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(dbp);
  Ravelin::VectorNd x,xd;
  part->get_generalized_coordinates( Moby::DynamicBody::eSpatial,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  std::cout << "x = " << x << std::endl;
  std::cout << "v = " << xd << std::endl;
}

// ============================================================================
// ================================ CALLBACKS =================================

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{

  boost::shared_ptr<Moby::EventDrivenSimulator> sim;
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

  part->controller                     = &controller_callback;

  // Set initial conditions from ruina paper

  Ravelin::VectorNd x,xd;
  part->get_generalized_coordinates( Moby::DynamicBody::eSpatial,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  x[1] =  0; // x
  x[2] = 0; // y
  x[3] = 1.1249; // z

  // Joint to RLEG
  x[6] =  0.09866765986740; // Phi
  x[4] = -0.00924861067616; // Psi
  x[5] = -0.16016583495522; // Theta_st
  x[0] =  3.43583389038583; // Theta_sw

  x[6] = -0.13220965510356; // Phi
  x[4] = -0.01990961987794; // Psi
  x[5] =  0.47124237466979; // Theta_st
  x[0] = -0.39255916866486; // Theta_sw

  part->set_generalized_coordinates( Moby::DynamicBody::eSpatial,x);
  part->set_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

}
} // end extern C
