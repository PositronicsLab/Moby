/*****************************************************************************
 * A program for updating kinematics given "control" loops 
 *****************************************************************************/

#include <errno.h>
#include <sys/time.h>
#include <dlfcn.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <boost/foreach.hpp>
#include <Moby/XMLReader.h>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Geode>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>
#include <Moby/Log.h>
#include <Moby/Simulator.h>
#include <Moby/RigidBody.h>
#include <Moby/EventDrivenSimulator.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Handle for dynamic library loading
void* HANDLE = NULL;

/// Beginning iteration for logging
unsigned LOG_START = 0;

/// Ending iteration for logging
unsigned LOG_STOP = std::numeric_limits<unsigned>::max();

/// The logging reporting level
unsigned LOG_REPORTING_LEVEL = 0;

/// The default simulation step size
const double DEFAULT_STEP_SIZE = .001;

/// The simulation step size
double STEP_SIZE = DEFAULT_STEP_SIZE;

/// The current simulation iteration
unsigned ITER = 0;

/// Determines whether to output timings
bool OUTPUT_TIMINGS = false;

/// The maximum number of iterations (default infinity)
unsigned MAX_ITER = std::numeric_limits<unsigned>::max(); 

/// The maximum time of the simulation (default infinity)
double MAX_TIME = std::numeric_limits<double>::max();

/// The total (CPU) clock time used by the simulation
double TOTAL_TIME = 0.0;

/// Last 3D output iteration and time output
unsigned LAST_3D_WRITTEN = -1;
double LAST_3D_WRITTEN_T = -std::numeric_limits<double>::max()/2.0;

/// Last image iteration output
unsigned LAST_IMG_WRITTEN = -1;
double LAST_IMG_WRITTEN_T = -std::numeric_limits<double>::max()/2.0;
      
/// Outputs to stdout
bool OUTPUT_FRAME_RATE = false;
bool OUTPUT_ITER_NUM = false;
bool OUTPUT_SIM_RATE = false;

/// Render Contact Points
bool RENDER_CONTACT_POINTS = false;

/// The map of objects read from the simulation XML file
std::map<std::string, BasePtr> READ_MAP;

/// The OpenInventor group node for Moby
osg::Group* MOBY_GROUP;

/// The OpenInventor root group node for this application
osg::Group* MAIN_GROUP;

/// Pointer to the viewer
osgViewer::Viewer* viewer_pointer;

/// Pointer to the controller's initializer, called once (if any)
typedef void (*init_t)(void*, const std::map<std::string, BasePtr>&, double);
std::list<init_t> INIT;

/// Checks whether was compiled with OpenSceneGraph support
bool check_osg()
{
  return true;
}

/// runs the simulator and updates all transforms
void step(void* arg)
{
  static double current_time = 0.0;

  // get the simulator pointer
  boost::shared_ptr<Simulator> s = *(boost::shared_ptr<Simulator>*) arg;

  // see whether to activate logging
  if (ITER >= LOG_START && ITER <= LOG_STOP)
    Log<OutputToFile>::reporting_level = LOG_REPORTING_LEVEL;
  else
    Log<OutputToFile>::reporting_level = 0;

  // output the iteration #
  if (OUTPUT_ITER_NUM)
    std::cout << "iteration: " << ITER << "  simulation time: " << current_time << std::endl;
  if (Log<OutputToFile>::reporting_level > 0)
    FILE_LOG(Log<OutputToFile>::reporting_level) << "iteration: " << ITER << "  simulation time: " << current_time << std::endl;

  // only update the graphics if it is necessary
  s->update_visualization();

  // step the dynamics
  BOOST_FOREACH(DynamicBodyPtr db, s->get_dynamic_bodies())
    if (db->controller)
      (*db->controller)(db, current_time, db->controller_arg);

  // update the current "time"
  current_time += STEP_SIZE;

  // update the iteration #
  ITER++;

  // check that maximum number of iterations or maximum time not exceeded
  if (ITER >= MAX_ITER || current_time > MAX_TIME)
    exit(0);
}

// attempts to read control code plugin
void read_plugin(const char* filename)
{
  // attempt to read the file
  HANDLE = dlopen(filename, RTLD_LAZY);
  if (!HANDLE)
  {
    std::cerr << "kinematics: failed to read plugin from " << filename << std::endl;
    std::cerr << "  " << dlerror() << std::endl;
    exit(-1);
  }

  // attempt to load the initializer
  dlerror();
  INIT.push_back((init_t) dlsym(HANDLE, "init"));
  const char* dlsym_error = dlerror();
  if (dlsym_error)
  {
    std::cerr << "kinematics warning: cannot load symbol 'init' from " << filename << std::endl;
    std::cerr << "        error follows: " << std::endl << dlsym_error << std::endl;
    INIT.pop_back();
  }
}

/// Adds lights to the scene when no scene background file specified
void add_lights()
{
}

/// Gets the XML sub-tree rooted at the specified tag
shared_ptr<const XMLTree> find_subtree(shared_ptr<const XMLTree> root, const std::string& name)
{
  // if we found the tree, return it
  if (strcasecmp(root->name.c_str(), name.c_str()) == 0)
    return root;

  // otherwise, look for it recursively
  const std::list<XMLTreePtr>& children = root->children;
  for (std::list<XMLTreePtr>::const_iterator i = children.begin(); i != children.end(); i++)
  {
    shared_ptr<const XMLTree> node = find_subtree(*i, name);
    if (node)
      return node;
  }

  // return NULL if we are here
  return shared_ptr<const XMLTree>();
}

// finds and processes given XML tags
void process_tag(const std::string& tag, shared_ptr<const XMLTree> root, void (*fn)(shared_ptr<const XMLTree>))
{
  // if this node is of the given type, process it 
  if (strcasecmp(root->name.c_str(), tag.c_str()) == 0)
    fn(root);
  else
  {
    const std::list<XMLTreePtr>& child_nodes = root->children;
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
      process_tag(tag, *i, fn);
  }
}

/// processes the 'camera' tag
void process_camera_tag(shared_ptr<const XMLTree> node)
{
  // read all attributes
  XMLAttrib* target_attr = node->get_attrib("target");
  XMLAttrib* position_attr = node->get_attrib("position");
  XMLAttrib* up_attr = node->get_attrib("up");
  if (!target_attr || !position_attr || !up_attr)
    return;

  // get the actual values
  Vector3d up;
  Point3d target, position;
  target_attr->get_vector_value(target);
  position_attr->get_vector_value(position);
  up_attr->get_vector_value(up);

  // setup osg vectors
  osg::Vec3d position_osg(position[0], position[1], position[2]);
  osg::Vec3d target_osg(target[0], target[1], target[2]);
  osg::Vec3d up_osg(up[0], up[1], up[2]);

  // set the camera view
  if (viewer_pointer && viewer_pointer->getCameraManipulator())
  {
    osg::Camera* camera = viewer_pointer->getCamera();
    camera->setViewMatrixAsLookAt(position_osg, target_osg, up_osg);

    // setup the manipulator using the camera, if necessary
    viewer_pointer->getCameraManipulator()->setHomePosition(position_osg, target_osg, up_osg); 
  }
}

/// processes the 'window' tag
void process_window_tag(shared_ptr<const XMLTree> node)
{
  // read window location
  XMLAttrib* loc_attr = node->get_attrib("location");

  // read window size
  XMLAttrib* size_attr = node->get_attrib("size");

  // get the actual values
  Vector2d loc(0,0), size(640,480);
  if (loc_attr)
    loc_attr->get_vector_value(loc);
  if (size_attr)
    size_attr->get_vector_value(size);

  // setup the window 
  viewer_pointer->setUpViewInWindow(loc[0], loc[1], size[0], size[1]);
}

/// processes all 'kinematics' options in the XML file
void process_xml_options(const std::string& xml_fname)
{
  // *************************************************************
  // going to remove any path from the argument and change to that
  // path; this is done so that all files referenced from the
  // local path of the XML file are found
  // *************************************************************

  // set the filename to use as the argument, by default
  std::string filename = xml_fname;

  // get the current pathname
  size_t BUFSIZE = 128;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "process_xml_options() - unable to allocate sufficient memory!" << std::endl;
      return;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = xml_fname.find_last_of('/');
  if (last_path_sep != std::string::npos)
  {
    // get the new working path
    std::string pathname = xml_fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = xml_fname.substr(last_path_sep+1,std::string::npos);
  }

  // read the XML Tree 
  shared_ptr<const XMLTree> kinematics_tree = XMLTree::read_from_xml(filename);
  if (!kinematics_tree)
  {
    std::cerr << "process_xml_options() - unable to open file " << xml_fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return;
  }
  
  // find the kinematics tree 
  kinematics_tree = find_subtree(kinematics_tree, "kinematics");

  // make sure that the kinematics node was found
  if (!kinematics_tree)
  {
    chdir(cwd.get());
    return;
  }

  // process tags
  process_tag("window", kinematics_tree, process_window_tag);
  process_tag("camera", kinematics_tree, process_camera_tag);

  // change back to current directory 
  chdir(cwd.get());
}

// where everything begins...
int main(int argc, char** argv)
{
  const unsigned ONECHAR_ARG = 3, TWOCHAR_ARG = 4;

  const double DYNAMICS_FREQ = 0.001;
  osgViewer::Viewer viewer;
  viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer_pointer = &viewer;

  // setup some default options
  std::string scene_path;

  // check that syntax is ok
  if (argc < 2)
  {
    std::cerr << "syntax: kinematics [OPTIONS] <xml file>" << std::endl;
    std::cerr << "        (see README for OPTIONS)" << std::endl;
    return -1;
  }

  // get all options
  for (int i=1; i< argc-1; i++)
  {
    // get the option
    std::string option(argv[i]);

    // process options
    if (option.find("-of") != std::string::npos)
      OUTPUT_FRAME_RATE = true;
    else if (option.find("-ot") != std::string::npos)
      OUTPUT_TIMINGS = true;
    else if (option.find("-oi") != std::string::npos)
      OUTPUT_ITER_NUM = true;
    else if (option.find("-or") != std::string::npos)
      OUTPUT_SIM_RATE = true;
    else if (option.find("-s=") != std::string::npos)
    {
      STEP_SIZE = std::atof(&argv[i][ONECHAR_ARG]);
      assert(STEP_SIZE >= 0.0 && STEP_SIZE < 1);
    }
    else if (option.find("-lf=") != std::string::npos)
    {
      std::string fname(&argv[i][TWOCHAR_ARG]);
      OutputToFile::stream.open(fname.c_str());
    }  
    else if (option.find("-l=") != std::string::npos)
    {
      LOG_REPORTING_LEVEL = std::atoi(&argv[i][ONECHAR_ARG]);
      Log<OutputToFile>::reporting_level = LOG_REPORTING_LEVEL;
    }
    else if (option.find("-lt=") != std::string::npos)
    {
      LOG_START = std::atoi(&argv[i][TWOCHAR_ARG]);
    }
    else if (option.find("-lp=") != std::string::npos)
    {
      LOG_STOP = std::atoi(&argv[i][TWOCHAR_ARG]);
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
    else if (option.find("-x=") != std::string::npos)
    {
      check_osg();
      scene_path = std::string(&argv[i][ONECHAR_ARG]);
    }
    else if (option.find("-p=") != std::string::npos)
      read_plugin(&argv[i][ONECHAR_ARG]);
  }

  // setup the simulation 
  READ_MAP = XMLReader::read(std::string(argv[argc-1]));

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
    std::cerr << "kinematics: no simulator found in " << argv[argc-1] << std::endl;
    return -1;
  } 

  // setup osg window
  viewer.setCameraManipulator(new osgGA::TrackballManipulator());
  viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
  viewer.addEventHandler(new osgViewer::WindowSizeHandler);
  viewer.addEventHandler(new osgViewer::StatsHandler);

  // init the main group
  MAIN_GROUP = new osg::Group;

  // call the initializers, if any
  if (!INIT.empty())
  {
    BOOST_FOREACH(init_t i, INIT)
      (*i)(MAIN_GROUP, READ_MAP, STEP_SIZE);
  }

  // look for a scene description file
  if (scene_path != "")
  {
    std::ifstream in(scene_path.c_str());
    if (in.fail())
    {
      std::cerr << "kinematics: unable to find scene description from " << scene_path << std::endl;
      add_lights();
    }
    else
    {
      in.close();
      osg::Node* node = osgDB::readNodeFile(scene_path);
      if (!node)
      {
        std::cerr << "kinematics: unable to open scene description file!" << std::endl;
        add_lights();
      }
      else
        MAIN_GROUP->addChild(node);
    }
  }
  else
    add_lights();

  // process XML options
  process_xml_options(std::string(argv[argc-1]));

  // get the simulator visualization
  MAIN_GROUP->addChild(s->get_persistent_vdata());
  MAIN_GROUP->addChild(s->get_transient_vdata());

  // prepare to render
  viewer.setSceneData(MAIN_GROUP);
  viewer.realize();

  // begin rendering
  while (true)
  {   
    if (viewer.done())
      break;
    viewer.frame();
    step((void*) &s);
    usleep(DYNAMICS_FREQ);
  }

  // close the loaded library
  if (HANDLE)
    dlclose(HANDLE);
}

