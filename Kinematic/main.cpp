#include <errno.h>
#include <pthread.h>
#include <boost/foreach.hpp>
#include <fstream>
#include <iostream>
#include <boost/foreach.hpp>
#include <Moby/XMLReader.h>
#include <Moby/XMLWriter.h>
#include "ObjectFormImpl.h"
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Geode>
#include <osg/Material>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>
#include <qapplication.h>
#include <qmutex.h>
#include <QTimer>

using namespace boost;
using namespace Moby;
using std::list;

// the mutex for updating graphics
ObjectFormImpl* g_oform = NULL;
QMutex _mutex;
osg::PolygonMode* polygon_mode = NULL;
osgViewer::Viewer* viewer_pointer = NULL;
osg::Node* g_persistent_vdata = NULL;
osg::Node* g_transient_vdata = NULL;
osg::Node* g_new_persistent_vdata = NULL;
osg::Node* g_new_transient_vdata = NULL;
osg::Group* g_moby_group = NULL;
osg::Group* g_main_group = NULL;
osg::Group* g_colliding_tris_group = NULL;
Vector2 g_loc(0,0), g_size(640,480);
Vector3 g_up(0,1,0), g_position(0,0,10), g_target(0,0,0);
bool g_change = false;
shared_ptr<Simulator> g_simulator;
list<list<Triangle> > g_colliding_tris;
list<Vector3> g_tri_colors;

// generates a uniformly distributed random value in [0,1]
Real randD()
{
  return (Real) rand() / RAND_MAX;
}

void* render(void* arg)
{
  const unsigned X = 0, Y = 1, Z = 2;
  typedef std::pair<Triangle, Triangle> TrianglePair;

  // keep running
  while (true) {

  // lock the mutex
  _mutex.lock();

  // see whether we are making changes to the visualization data (a new
  // XML file has been loaded)
  if (g_change)
  {
    assert(g_new_persistent_vdata);
    assert(g_new_transient_vdata);
    if (g_persistent_vdata)
      g_moby_group->removeChild(g_persistent_vdata);
    if (g_transient_vdata)
      g_moby_group->removeChild(g_transient_vdata);
    g_persistent_vdata = g_new_persistent_vdata;
    g_transient_vdata = g_new_transient_vdata;
    g_moby_group->addChild(g_persistent_vdata);
    g_moby_group->addChild(g_transient_vdata);
    g_change = false;
  }

  // only update visualization if we have a simulator; 
  if (g_simulator)
    g_simulator->update_visualization();

  // remove all colliding triangles
  g_colliding_tris_group->removeChildren(0, g_colliding_tris_group->getNumChildren());

  // ******
  // following code draws colliding triangles (if any)
  // ******
  if (g_colliding_tris.empty())
  {
    // draw everything as normal
    polygon_mode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::FILL);
  }
  else
  {
    // draw wireframe
    polygon_mode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE);

    // setup an iterator for the colors
    list<Vector3>::iterator color_it = g_tri_colors.begin();

    // process every set of colliding triangles
    BOOST_FOREACH(const list<Triangle>& tris, g_colliding_tris)
    {
      // create a new group
      osg::Group* subsep = new osg::Group;
      g_colliding_tris_group->addChild(subsep);

      // create a material (and a new color, if necessary)
      osg::Material* mat = new osg::Material;
      if (color_it == g_tri_colors.end())
      {
        Vector3 c(randD()*0.5+0.5,randD()*0.5+0.5,randD()*0.5+0.5);
        g_tri_colors.push_back(c);
      }
      Vector3& c = *color_it++;
      mat->setColorMode(osg::Material::DIFFUSE);
      mat->setDiffuse(osg::Material::FRONT, osg::Vec4(c[0], c[1], c[2], 1.0f));
      subsep->getOrCreateStateSet()->setAttribute(mat);

      // create a new coordinate object
      osg::Geode* geode = new osg::Geode;
      osg::Geometry* geom = new osg::Geometry;
      geode->addDrawable(geom);
      subsep->addChild(geode);

      // create the vertex array
      osg::Vec3Array* varray = new osg::Vec3Array(tris.size()*3);
      unsigned idx = 0;
      BOOST_FOREACH(const Triangle& tri, tris)
      {
        (*varray)[idx++] = osg::Vec3((float) tri.a[X], (float) tri.a[Y], (float) tri.a[Z]);
        (*varray)[idx++] = osg::Vec3((float) tri.b[X], (float) tri.b[Y], (float) tri.b[Z]);
        (*varray)[idx++] = osg::Vec3((float) tri.c[X], (float) tri.c[Y], (float) tri.c[Z]);
      }
      geom->setVertexArray(varray);

      // create the faces
      for (unsigned i=0; i< tris.size(); i++)
      {
        osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
        face->push_back(i*3);
        face->push_back(i*3+1);
        face->push_back(i*3+2);
        geom->addPrimitiveSet(face);
      }
    }
  }

  // release the mutex
  _mutex.unlock();

  if (viewer_pointer->done())
    break;
  viewer_pointer->frame();
  usleep(1000);   // sleep for 1ms
  } // end while loop

  return NULL;
}

/// Adds generic lights and camera to a scene
void add_lights_camera()
{
}

/// Reads simulator from fname into the g_simulator pointer; also loads visualization data 
void read_simulator(const std::string& fname)
{
  // read the simulator(s)
  std::map<std::string, BasePtr> read_map = XMLReader::read(fname);

  // get all simulator objects
  std::list<SimulatorPtr> simulators;
  for (std::map<std::string, BasePtr>::const_iterator i = read_map.begin(); i != read_map.end(); i++)
  {
    SimulatorPtr s = dynamic_pointer_cast<Simulator>(i->second);
    if (s)
      simulators.push_back(s);
  }

  // make sure that exactly one simulator was found
  if (simulators.empty())
    std::cerr << " --> no simulator found in XML file '" << fname << "'" << std::endl;
  else if (simulators.size() > 1)
    std::cerr << " --> multiple simulators found in XML file '" << fname << "'" << std::endl;
  else
  {
    // update the simulator in the ObjectForm (if the form has been created)
    if (g_oform)
      g_oform->set_simulator(simulators.front());

    // setup the new persistent and transient visualization data
    g_new_persistent_vdata = simulators.front()->get_persistent_vdata();
    g_new_transient_vdata = simulators.front()->get_transient_vdata();

    // save the simulator
    g_simulator = simulators.front();

    // indicate that we are making a change
    g_change = true;
  }
}

/// Gets the XML sub-tree rooted at the specified tag
XMLTreeConstPtr find_subtree(XMLTreeConstPtr root, const std::string& name)
{
  // if we found the tree, return it
  if (strcasecmp(root->name.c_str(), name.c_str()) == 0)
    return root;

  // otherwise, look for it recursively
  const std::list<XMLTreePtr>& children = root->children;
  for (std::list<XMLTreePtr>::const_iterator i = children.begin(); i != children.end(); i++)
  {
    XMLTreeConstPtr node = find_subtree(*i, name);
    if (node)
      return node;
  }

  // return NULL if we are here
  return XMLTreeConstPtr();
}

// finds and processes given XML tags
void process_tag(const std::string& tag, XMLTreeConstPtr root, void (*fn)(XMLTreeConstPtr))
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
void process_camera_tag(XMLTreeConstPtr node)
{
  // read all attributes
  const XMLAttrib* target_attr = node->get_attrib("target");
  const XMLAttrib* position_attr = node->get_attrib("position");
  const XMLAttrib* up_attr = node->get_attrib("up");
  if (!target_attr || !position_attr || !up_attr)
    return;

  // get the actual values
  target_attr->get_vector_value(g_target);
  position_attr->get_vector_value(g_position);
  up_attr->get_vector_value(g_up);

  // setup osg vectors
  osg::Vec3d position_osg(g_position[0], g_position[1], g_position[2]);
  osg::Vec3d target_osg(g_target[0], g_target[1], g_target[2]);
  osg::Vec3d up_osg(g_up[0], g_up[1], g_up[2]);

  // set the camera view
  osg::Camera* camera = viewer_pointer->getCamera();
  camera->setViewMatrixAsLookAt(position_osg, target_osg, up_osg);

  // setup the manipulator using the camera, if necessary
  if (viewer_pointer->getCameraManipulator())
    viewer_pointer->getCameraManipulator()->setHomePosition(position_osg, target_osg, up_osg); 
}

/// processes the 'window' tag
void process_window_tag(XMLTreeConstPtr node)
{
  // read window location
  const XMLAttrib* loc_attr = node->get_attrib("location");

  // read window size
  const XMLAttrib* size_attr = node->get_attrib("size");

  // get the actual values
  if (loc_attr)
    loc_attr->get_vector_value(g_loc);
  if (size_attr)
    size_attr->get_vector_value(g_size);

  // setup the window 
  viewer_pointer->setUpViewInWindow(g_loc[0], g_loc[1], g_size[0], g_size[1]);
}

/// processes all 'driver' options in the XML file
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
  XMLTreeConstPtr driver_tree = XMLTree::read_from_xml(filename);
  if (!driver_tree)
  {
    std::cerr << "process_xml_options() - unable to open file " << xml_fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return;
  }
  
  // find the driver tree 
  driver_tree = find_subtree(driver_tree, "driver");

  // make sure that the driver node was found
  if (!driver_tree)
  {
    chdir(cwd.get());
    return;
  }

  // process tags
  process_tag("window", driver_tree, process_window_tag);
  process_tag("camera", driver_tree, process_camera_tag);

  // change back to current directory 
  chdir(cwd.get());
}

int main(int argc, char ** argv )
{
  std::string scene_path;

  // setup the application
  QApplication qapp(argc, argv, true);

  // see whether an XML file was specified
  bool xml_spec = (argc > 1 && argv[argc-1][0] != '-');

  // read the simulator, if given
  if (xml_spec)
  {
    // note: we decrement argc for the options loop, below
    argc--;
    read_simulator(argv[argc]);
  }

  // setup the state set for all visualization
  osg::StateSet* state_set = new osg::StateSet;
  polygon_mode = new osg::PolygonMode;
  polygon_mode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::FILL);
  state_set->setAttributeAndModes(polygon_mode, osg::StateAttribute::ON);

  // create the group to hold all Moby data
  g_moby_group = new osg::Group;
  g_moby_group->setStateSet(state_set);

  // add a group for drawing colliding triangles
  g_colliding_tris_group = new osg::Group;

  // setup the state set for colliding triangles
  osg::StateSet* cstate_set = new osg::StateSet;
  osg::PolygonMode* pmode = new osg::PolygonMode;
  pmode->setMode(osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::FILL);
  cstate_set->setAttributeAndModes(pmode, osg::StateAttribute::ON);
  g_colliding_tris_group->setStateSet(cstate_set);

  // create the main group
  g_main_group = new osg::Group;

  // get all options
  for (int i=1; i< argc; i++)
  {
    // options *must* be at the beginning
    std::string option(argv[i]);
    if (option.find('-') == std::string::npos && xml_spec && i != argc-1)
    {
      std::cerr << "kinematic: options must be specified before XML file!" << std::endl;
      return -1;
    }
      
    // look for scene description file
    if (option.find("-x=") != std::string::npos)
      scene_path = std::string(&argv[i][3]);
    else
      std::cerr << "kinematic: unknown option '" << option << "'" << std::endl;
  }    

  // look for a scene description file, if specified
  if (scene_path != "")
  {
    std::ifstream in(scene_path.c_str());
    if (in.fail())
    {
      std::cerr << "kinematic: unable to find scene description from " << scene_path << std::endl;
      add_lights_camera();
    }
    else
    {    
      in.close();
      osg::Node* node = osgDB::readNodeFile(scene_path);
      if (!node)
      {
        std::cerr << "kinematic: unable to open scene description file" << std::endl;
        add_lights_camera();
      }
      else
        g_main_group->addChild(node); 
    }
  }
  else
    add_lights_camera();

  // now that lights and cameras are added to the main group; add lights
  // and camera
  g_main_group->addChild(g_colliding_tris_group);
  g_main_group->addChild(g_moby_group);

  // get the simulator visualization, if any
  if (g_simulator)
  {
    g_persistent_vdata = g_simulator->get_persistent_vdata();
    g_transient_vdata = g_simulator->get_transient_vdata();
    g_moby_group->addChild(g_persistent_vdata);
    g_moby_group->addChild(g_transient_vdata);
  }
  
  // create the object form
  g_oform = new ObjectFormImpl(g_simulator, &_mutex, &read_simulator, &g_colliding_tris);
  g_oform->show();

  // setup OSG
  osgViewer::Viewer viewer;
  viewer_pointer = &viewer;
  viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer.setCameraManipulator(new osgGA::TrackballManipulator());
  viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
  viewer.addEventHandler(new osgViewer::WindowSizeHandler);
  viewer.addEventHandler(new osgViewer::StatsHandler);

  // process XML options
  if (xml_spec)
    process_xml_options(std::string(argv[argc]));

  // start rendering
  viewer.setSceneData(g_main_group);
  viewer.realize();
  pthread_t thread;
  pthread_create(&thread, NULL, &render, NULL);

  // start processing GUI events
  qapp.exec();  

  // delete the object form
  delete g_oform;

  return 0;
}

