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
  const double DT = 1e-3;
  const unsigned N_STEPS = (unsigned) (3.0/DT);
  const std::string FNAME("sphere.xml");
  const std::string SPHERE_ID("sphere");
  shared_ptr<TimeSteppingSimulator> sim;
  shared_ptr<RigidBody> sphere;

  // log contact 
  Moby::Log<Moby::OutputToFile>::reporting_level = (LOG_SIMULATOR | LOG_CONSTRAINT);
  Moby::OutputToFile::stream.open("logging.out");

  // load in the sphere file
  map<std::string, BasePtr> READ_MAP = XMLReader::read(FNAME);
  
  // look for the simulator
  find(READ_MAP, sim);

  // look for the sphere 
  find(READ_MAP, sphere, SPHERE_ID);

  // modify the initial conditions

  // compute the kinetic energy in the global frame
  double KE1 = sphere->calc_kinetic_energy(sphere->get_mixed_pose());

  // integrate for a few seconds
  for (unsigned i=0; i< N_STEPS; i++)
    sim->step(DT);

  // compute the kinetic energy in the global frame again 
  double KE2 = sphere->calc_kinetic_energy(sphere->get_mixed_pose());

  Moby::OutputToFile::stream.close();

  // compare the two kinetic energies
  EXPECT_NEAR(KE1, KE2, TOL);
}

