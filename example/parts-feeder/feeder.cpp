/*****************************************************************************
 * Controller for LINKS robot
 ****************************************************************************/
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>

 boost::shared_ptr<Moby::EventDrivenSimulator> sim;
 Moby::RCArticulatedBodyPtr slide;

void controller_callback(Moby::ControlledBodyPtr dbp, double t, void*)
{
  double p =  0.001*cos(t*(2.0*M_PI)*500.0);
  double v =  -0.001*sin(t*1000.0) *(2.0*M_PI)*500.0;
  double Kp = 1e3
      ,Kv = 1e1
              ;

  for(int i=0;i<slide->get_joints().size();i++){
    if(slide->get_joints()[i]->num_dof() == 0) continue;
    static Ravelin::VectorNd U(1);
    U[0] = Kp*(p - slide->get_joints()[i]->q[0]) + Kv*(v - slide->get_joints()[i]->qd[0]);
    std::cout << "x = [ " << slide->get_joints()[i]->q[0] << " , " << slide->get_joints()[i]->qd[0] << " ];" << std::endl;
    std::cout << "x* = [ " << p << " , " << v << " ];" << std::endl;
    std::cout << "U = [ " << U[0] << " ];" << std::endl;
    slide->get_joints()[i]->add_force(U);
  }
}

// ============================================================================
// ================================ CALLBACKS =================================

/// plugin must be "extern C"
extern "C" {

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double time)
{
  // If use robot is active also init dynamixel controllers
  // get a reference to the EventDrivenSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<Moby::EventDrivenSimulator>(i->second);

    // find the robot reference
    if (!slide)
    {
      slide = boost::dynamic_pointer_cast<Moby::RCArticulatedBody>(i->second);
    }
  }

  slide->controller                     = &controller_callback;

}
} // end extern C
