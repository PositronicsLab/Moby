/*****************************************************************************
 * "Initializer" for ellipse 
 ****************************************************************************/
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <Moby/ConstraintSimulator.h>
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

Moby::RigidBodyPtr ellipse, ground, wall_pos, wall_neg;
boost::shared_ptr<ConstraintSimulator> sim;

// normal sample generator 
boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normal_gen(boost::mt19937(time(0)), boost::normal_distribution<>());

// exponential sample generator
boost::variate_generator<boost::mt19937, boost::exponential_distribution<> > exp_gen(boost::mt19937(time(0)), boost::exponential_distribution<>());

/// plugin must be "extern C"
extern "C" {

double randu() { return (double) rand() / RAND_MAX; }

void perturb_initial_conditions()
{
  const double LV_PERTURB = 0.1;
  const double AV_PERTURB = 0.05;
  const double ELLIPSE_RADIUS_PERTURB = 0.1;
  VectorNd gc, gv;
  Point3d dummy1, dummy2;

  // perturb the initial position
  ellipse->get_generalized_coordinates_euler(gc);
  gc[0] += normal_gen();
  gc[1] += normal_gen();

  // perturb the initial rotation
  Quatd quat = Matrix3d::rot_Z(randu() * M_PI * 2.0);
  gc[3] = quat.x;  
  gc[4] = quat.y;  
  gc[5] = quat.z;  
  gc[6] = quat.w;  

  // perturb the ellipse radius
  std::ostringstream x_rad, y_rad;
//  x_rad << (randu() * ELLIPSE_RADIUS_PERTURB * 2.0 - ELLIPSE_RADIUS_PERTURB);
//  y_rad << (randu() * ELLIPSE_RADIUS_PERTURB * 2.0 - ELLIPSE_RADIUS_PERTURB);
x_rad << "1.0";
y_rad << "2.0";
  setenv("ELLIPSE_X_RAD", x_rad.str().c_str(), 1);
  setenv("ELLIPSE_Y_RAD", y_rad.str().c_str(), 1);

  // ensure that the initial configuration is not in contact with the ground
  VectorNd newgc = gc;
  ellipse->set_generalized_coordinates_euler(newgc);

  // get the collision detector
  shared_ptr<CollisionDetection> coldet = sim->get_collision_detection();

  // get the collision geometries for the ellipse and the ground
  shared_ptr<CollisionGeometry> ground_cg = ground->geometries.front();
  shared_ptr<CollisionGeometry> ellipse_cg = ellipse->geometries.front();

  // check for contact between ground and ellipse
  if (coldet->calc_signed_dist(ground_cg, ellipse_cg, dummy1, dummy2) <= 0.0)
  {
    // reset the generalized coords
    ellipse->set_generalized_coordinates_euler(gc);

    // try again
    perturb_initial_conditions();
    return;
  }

  // perturb the initial velocity 
  ellipse->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  for (unsigned i=0; i< 3; i++)
    gv[i+3] += normal_gen()*LV_PERTURB; 
    gv[5] += normal_gen()*AV_PERTURB; 
  ellipse->set_generalized_velocity(DynamicBodyd::eSpatial, gv);
}

// perturbs contact parameters between the ground and the ellipse
void perturb_contact_parameters()
{
  // perturb the contact parameters between the ellipse and the ground
  ContactParameters& cp_ground = *sim->contact_params[make_sorted_pair(ellipse, ground)];
  cp_ground.epsilon = randu();
  cp_ground.mu_coulomb = randu();

  // perturb the contact parameters between the ellipse and the positive wall 
  ContactParameters& cp_pos = *sim->contact_params[make_sorted_pair(ellipse, wall_pos)];
  cp_pos.epsilon = randu();
  cp_pos.mu_coulomb = randu()*0.01;

  // perturb the contact parameters between the ellipse and the negative wall 
  ContactParameters& cp_neg = *sim->contact_params[make_sorted_pair(ellipse, wall_neg)];
  cp_neg.epsilon = randu();
  cp_neg.mu_coulomb = randu()*0.01;
}

void init(void* separator, const std::map<std::string, Moby::BasePtr>& read_map, double t)
{
  // randomize the seed
  srand(time(0));

  // get a reference to the ConstraintSimulator instance
  for (std::map<std::string, Moby::BasePtr>::const_iterator i = read_map.begin();
       i !=read_map.end(); i++)
  {
    // Find the simulator reference
    if (!sim)
      sim = boost::dynamic_pointer_cast<ConstraintSimulator>(i->second);
    if (i->first == "ellipse")
      ellipse = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "ground")
      ground = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "wall_pos")
      wall_pos = boost::dynamic_pointer_cast<RigidBody>(i->second);
    if (i->first == "wall_neg")
      wall_neg = boost::dynamic_pointer_cast<RigidBody>(i->second);
  }

  // perturb the initial conditions for the ellipse
  perturb_initial_conditions();

  // perturb the contact parameters for the ellipse
  perturb_contact_parameters();
}
} // end extern C
