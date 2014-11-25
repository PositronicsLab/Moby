/*****************************************************************************
 * Controller for LINKS robot
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
boost::shared_ptr<Moby::EventDrivenSimulator> sim;

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
      std::cout << "s = " << e[i].contact_tan1 << std::endl;
      std::cout << "t = " << e[i].contact_tan2 << std::endl;
      std::cout << "muC = " << e[i].contact_mu_coulomb << std::endl;
      std::cout << "muV = " << e[i].contact_mu_viscous << std::endl;

//      visualize_ray(  e[i].contact_point,
//                      e[i].contact_point + e[i].contact_impulse.get_linear()*10.0,
//                      Ravelin::Vector3d(1,0.5,0),
//                      0.1,
//                      sim
//                    );
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

  std::cout << "x = " << x << std::endl;
  std::cout << "v = " << xd << std::endl;
  std::cout << "<< end controller_callback(.)" << std::endl;
}

// ============================================================================
// ================================ CALLBACKS =================================

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
  part->get_generalized_coordinates( Moby::DynamicBody::eSpatial,x);
  part->get_generalized_velocity( Moby::DynamicBody::eSpatial,xd);

  x[1] =  0; // x
  x[2] = 0; // y
  x[3] = 1.1249; // z
//  x[3] = 1.4; // z

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
