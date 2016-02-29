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

TEST(Boxes, ConstraintViolation)
{
  const double TOL = 1e-6;
  const double DT = 1e-2;
  const std::string FNAME("box.xml");
  const std::string BOX_ID("box");
  shared_ptr<TimeSteppingSimulator> sim;
  shared_ptr<RigidBody> box;
  const unsigned NTIMES = 100;

  // log contact 
  Moby::Log<Moby::OutputToFile>::reporting_level = (LOG_SIMULATOR | LOG_CONSTRAINT);
  Moby::OutputToFile::stream.open("logging.out");

  // log energy
  std::ofstream energy_out("energy.dat");

  // do this multiple times
  for (unsigned j=0; j< NTIMES; j++)
  {
    // load in the box file
    map<std::string, BasePtr> READ_MAP = XMLReader::read(FNAME);
  
    // look for the simulator
    find(READ_MAP, sim);

    // look for the box 
    find(READ_MAP, box, BOX_ID);

    // set tolerance for constraint stabilization and minimum step size
//    sim->min_step_size = TOL/10000.0;
//    sim->min_step_size = 1e-5;
    sim->min_step_size = 1e-1;
    sim->cstab.eps = -NEAR_ZERO;

    // modify the initial orientation of the box
    Quatd q;
    q.x = (double) rand() / RAND_MAX * 2.0 - 1.0;
    q.y = (double) rand() / RAND_MAX * 2.0 - 1.0;
    q.z = (double) rand() / RAND_MAX * 2.0 - 1.0;
    q.w = (double) rand() / RAND_MAX * 2.0 - 1.0;
    q.normalize();
    Pose3d P = *box->get_pose();
    P.q = q;
    box->set_pose(P);

    // modify the velocity of the box
    Vector3d xd(box->get_inertial_pose());
    Vector3d w(box->get_inertial_pose());
    SVelocityd v(box->get_inertial_pose());
    for (unsigned i=0; i< 3; i++)
    {
      xd[i] = (double) rand() / RAND_MAX * 2.0 - 1.0;
      w[i] = (double) rand() / RAND_MAX * 2.0 - 1.0;
    }
    v.set_linear(xd);
    v.set_angular(w);
    box->set_velocity(v);
    
    // set the first no K.E. time 
    double no_KE_time = 0.0;

    // setup the maximum violation
    double max_vio = 0.0;

    // reset steps
    unsigned step = 0;

    // integrate the box until it has no kinetic energy for 1/2 sec 
    while (true)
    {
      // execute the step
      sim->step(DT);
 
      // get the constraint violation
      const vector<PairwiseDistInfo>& pdi = sim->get_pairwise_distances();
      for (unsigned i=0; i< pdi.size(); i++)
        max_vio = std::min(max_vio, pdi[i].dist);

      // get the kinetic energy of the box
      double KE = box->calc_kinetic_energy();
      energy_out << KE << std::endl;

      // see whether there is no kinetic energy
      if (KE < 1e-6)
      {
        if (sim->current_time - no_KE_time > 0.5)
          break;
      }
      else
        no_KE_time = sim->current_time;
      if (step++ % 1000 == 0)
        std::cerr << "." << std::flush;
    }

    // only want to print out one message about violation
    EXPECT_GT(max_vio, -TOL);
    std::cerr << "+" << std::flush;
    energy_out.close();
  }
}

