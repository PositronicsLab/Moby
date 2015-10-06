/*****************************************************************************
 * Utility for producing regression testing data 
 *****************************************************************************/

#include <errno.h>
#include <sys/time.h>
#include <dlfcn.h>
#include <cmath>
#include <fstream>
#include <boost/foreach.hpp>
#include <Moby/XMLReader.h>
#include <Moby/Log.h>
#include <Moby/Simulator.h>
#include <Moby/RigidBody.h>
#include <Ravelin/DynamicBodyd.h>

using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using Ravelin::VectorNd;
using Ravelin::DynamicBodyd;
using namespace Moby;

/// Handle for dynamic library loading
std::vector<void*> handles;

/// The current simulation iteration
unsigned ITER = 0;

/// The time the simulation starts
clock_t start_time; 

/// The default simulation step size
const double DEFAULT_STEP_SIZE = .001;

/// The simulation step size
double STEP_SIZE = DEFAULT_STEP_SIZE;

/// Total (CPU) clock time used by the simulation
double TOTAL_TIME = 0.0;

/// The maximum number of iterations (default infinity)
unsigned MAX_ITER = std::numeric_limits<unsigned>::max(); 

/// The maximum time of the simulation (default infinity)
double MAX_TIME = std::numeric_limits<double>::max();

/// The output file
std::ofstream outfile;

/// Outputs to stdout
bool OUTPUT_ITER_NUM = false;
bool OUTPUT_SIM_RATE = false;

/// The map of objects read from the simulation XML file
std::map<std::string, BasePtr> READ_MAP;

/// Pointer to the controller's initializer, called once (if any)
typedef void (*init_t)(void*, const std::map<std::string, BasePtr>&, double);
std::list<init_t> INIT;

/// Gets the current time (as a floating-point number)
double get_current_time()
{
  const double MICROSEC = 1.0/1000000;
  timeval t;
  gettimeofday(&t, NULL);
  return (double) t.tv_sec + (double) t.tv_usec * MICROSEC;
}

bool compbody(ControlledBodyPtr b1, ControlledBodyPtr b2)
{
  return b1->id < b2->id;
}

/// runs the simulator and updates all transforms
bool step(void* arg)
{
  // get the simulator pointer
  boost::shared_ptr<Simulator> s = *(boost::shared_ptr<Simulator>*) arg;

  // get the generalized coordinates for all bodies in alphabetical order
  std::vector<ControlledBodyPtr> bodies = s->get_dynamic_bodies();
  std::sort(bodies.begin(), bodies.end(), compbody);
  VectorNd q;
  outfile << s->current_time;
  for (unsigned i=0; i< bodies.size(); i++)  
  {
    shared_ptr<DynamicBodyd> db = dynamic_pointer_cast<DynamicBodyd>(bodies[i]);
    db->get_generalized_coordinates_euler(q);
    for (unsigned j=0; j< q.size(); j++)
      outfile << " " << q[j];
  }
  outfile << std::endl;

  // output the iteration #
  if (OUTPUT_ITER_NUM)
    std::cout << "iteration: " << ITER << "  simulation time: " << s->current_time << std::endl;
  if (Log<OutputToFile>::reporting_level > 0)
    FILE_LOG(Log<OutputToFile>::reporting_level) << "iteration: " << ITER << "  simulation time: " << s->current_time << std::endl;

  // step the simulator and update visualization
  clock_t pre_sim_t = clock();
  s->step(STEP_SIZE);
  clock_t post_sim_t = clock();
  double total_t = (post_sim_t - pre_sim_t) / (double) CLOCKS_PER_SEC;
  TOTAL_TIME += total_t;

  // output the iteration / stepping rate
  if (OUTPUT_SIM_RATE)
    std::cout << "time to compute last iteration: " << total_t << " (" << TOTAL_TIME / ITER << "s/iter, " << TOTAL_TIME / s->current_time << "s/step)" << std::endl;

  // update the iteration #
  ITER++;

  // check that maximum number of iterations or maximum time not exceeded
  if (ITER >= MAX_ITER || s->current_time > MAX_TIME)
    return false;

  return true;
}

// attempts to read control code plugin
void read_plugin(const char* filename)
{
  // attempt to read the file
  void* plugin = dlopen(filename, RTLD_LAZY);
  if (!plugin)
  {
    // get the error string, in case we need it
    char* dlerror_str = dlerror();

    // attempt to use the plugin path
    char* plugin_path = getenv("MOBY_PLUGIN_PATH");
    if (plugin_path)
    {
      // get the plugin path and make sure it has a path string at the end
      std::string plugin_path_str(plugin_path);
      if (plugin_path_str.at(plugin_path_str.size()-1) != '/')
        plugin_path_str += '/';

      // concatenate
      plugin_path_str += filename;

      // attempt to re-open the plugin
      plugin = dlopen(plugin_path_str.c_str(), RTLD_LAZY);
    }

    // check whether the plugin was successfully loaded
    if (!plugin)
    { 
      std::cerr << "regress: failed to read plugin from " << filename << " (or using " << plugin_path << ")" << std::endl;
      std::cerr << "  " << dlerror_str << std::endl;
      exit(-1);
    }
  }
  handles.push_back(plugin);

  // attempt to load the initializer
  dlerror();
  INIT.push_back((init_t) dlsym(plugin, "init"));
  const char* dlsym_error = dlerror();
  if (dlsym_error)
  {
    std::cerr << "regress warning: cannot load symbol 'init' from " << filename << std::endl;
    std::cerr << "        error follows: " << std::endl << dlsym_error << std::endl;
    exit(-1);
  }
}

// where everything begins...
int main(int argc, char** argv)
{
  const unsigned ONECHAR_ARG = 3, TWOCHAR_ARG = 4;

  // check that syntax is ok
  if (argc < 3)
  {
    std::cerr << "syntax: regress [OPTIONS] <xml file> <output file>" << std::endl;
    return -1;
  }

  // get all options
  for (int i=1; i< argc-2; i++)
  {
    // get the option
    std::string option(argv[i]);

    // process options
    if (option.find("-oi") != std::string::npos)
      OUTPUT_ITER_NUM = true;
    else if (option.find("-or") != std::string::npos)
      OUTPUT_SIM_RATE = true;
    else if (option.find("-s=") != std::string::npos)
    {
      STEP_SIZE = std::atof(&argv[i][ONECHAR_ARG]);
      assert(STEP_SIZE >= 0.0 && STEP_SIZE < 1);
    }
    else if (option.find("-mi=") != std::string::npos)
    {
      MAX_ITER = std::atoi(&argv[i][TWOCHAR_ARG]);
      assert(MAX_ITER > 0);
    }
    else if (option.find("-mt=") != std::string::npos)
    {
      MAX_TIME = std::atof(&argv[i][TWOCHAR_ARG]);
      assert(MAX_TIME > 0);
    }
    else if (option.find("-p=") != std::string::npos)
      read_plugin(&argv[i][ONECHAR_ARG]);
  }

  // setup the simulation 
  READ_MAP = XMLReader::read(std::string(argv[argc-2]));

  // get the (only) simulation object
  boost::shared_ptr<Simulator> s;
  for (std::map<std::string, BasePtr>::const_iterator i = READ_MAP.begin(); i != READ_MAP.end(); i++)
  {
    s = boost::dynamic_pointer_cast<Simulator>(i->second);
    if (s)
      break;
  }

  // make sure that a simulator was found
  if (!s)
  {
    std::cerr << "regress: no simulator found in " << argv[argc-2] << std::endl;
    return -1;
  } 

  // setup the output file
  outfile.open(argv[argc-1]);

  // call the initializers, if any
  if (!INIT.empty())
  {
    BOOST_FOREACH(init_t i, INIT)
      (*i)(NULL, READ_MAP, STEP_SIZE);
  }

  // begin timing
  start_time = clock();

  // begin rendering
  while (step((void*) &s))
  {
  }

  // close the loaded library
  for(size_t i = 0; i < handles.size(); ++i){
    dlclose(handles[i]);
  }

  // write the number of clock ticks elapsed
  clock_t end_time = clock();
  double elapsed = (end_time - start_time) / (double) CLOCKS_PER_SEC;
  outfile << elapsed << std::endl;

  // close the output file
  outfile.close();
}

