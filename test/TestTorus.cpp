#include <Moby/XMLReader.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include "gtest/gtest.h"

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::map;
using std::vector;
using namespace Ravelin;
using namespace Moby;

static double sqr(double x) { return x*x; }

static void to_rpy(const Matrix3d& R, double& r, double& p, double& y)
{
  const unsigned X = 0, Y = 1, Z = 2;

  r = std::atan2(R(Y,X), R(X,X));
  p = std::atan2(-R(Z,X), std::sqrt(sqr(R(Z,Y)) + sqr(R(Z,Z))));
  y = std::atan2(R(Z,Y), R(Z,Z));
}

// outputs data for the torus
static void output(double t, RigidBodyPtr torus)
{
  const Origin3d YAXIS(0.0, 1.0, 0.0);
  const Origin3d ZAXIS(0.0, 0.0, 1.0);

  // get the orientation of the torus
  Pose3d P(*torus->get_pose());
  P.update_relative_pose(GLOBAL);
  Matrix3d R = P.q;

  // get the current angular velocity in the global frame
  Vector3d omega = Pose3d::transform(GLOBAL, torus->get_velocity()).get_angular();

  // determine time derivative of yaw
  double ydot = omega[2];

  // get roll, pitch, and yaw
  double r, p, y;
  to_rpy(R, r, p, y);

  // open the file for writing
  std::ofstream out("torus.dat", std::ostream::app);
  out << r << " " << p << " " << y << " " << ydot << std::endl;
  out.close();
}

void find(const map<std::string, BasePtr>& read_map, shared_ptr<RigidBody>& rb, const std::string& id)
{
  for (map<std::string, BasePtr>::const_iterator i = read_map.begin(); i != read_map.end(); i++)
  {
     rb = dynamic_pointer_cast<RigidBody>(i->second);
     if (rb && rb->id == id)
       return;
  }
}

void find(const map<std::string, BasePtr>& read_map, shared_ptr<TimeSteppingSimulator>& sim)
{
  for (map<std::string, BasePtr>::const_iterator i = read_map.begin(); i != read_map.end(); i++)
  {
     sim = dynamic_pointer_cast<TimeSteppingSimulator>(i->second);
     if (sim)
       return;
  }
}

TEST(Torus, KineticEnergy)
{
  const double TOL = 1e-6;
  const double DT = 1e-4;
  const unsigned N_STEPS = (unsigned) (3.0/DT);
  const std::string FNAME("torus.xml");
  const std::string TORUS_ID("torus");
  shared_ptr<TimeSteppingSimulator> sim;
  shared_ptr<RigidBody> torus;

  // log contact 
  Moby::Log<Moby::OutputToFile>::reporting_level = (LOG_SIMULATOR | LOG_CONSTRAINT);
  Moby::OutputToFile::stream.open("logging.out");

  // load in the torus file
  map<std::string, BasePtr> READ_MAP = XMLReader::read(FNAME);
  
  // look for the simulator
  find(READ_MAP, sim);

  // set the minimum step size on the simulator
  sim->min_step_size = std::min(DT, sim->min_step_size);

  // look for the torus 
  find(READ_MAP, torus, TORUS_ID);

  // modify the initial conditions

  // compute the kinetic energy in the global frame
  double KE1 = torus->calc_kinetic_energy(torus->get_mixed_pose());

  // delete the output file
  std::ofstream out("torus.dat");
  out.close();

  // integrate for a few seconds
  for (unsigned i=0; i< N_STEPS; i++)
  {
    output(sim->current_time, torus);
    sim->step(DT);
  }

  // compute the kinetic energy in the global frame again 
  double KE2 = torus->calc_kinetic_energy(torus->get_mixed_pose());

  Moby::OutputToFile::stream.close();

  // compare the two kinetic energies
  EXPECT_NEAR(KE1, KE2, TOL);
}

