/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <dlfcn.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <stack>
#include <queue>

#ifdef USE_OSG
#include <Moby/OSGGroupWrapper.h>
#endif

#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/Constants.h>
#include <Moby/Simulator.h>
#include <Moby/EventDrivenSimulator.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/RigidBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/FixedJoint.h>
//#include <Moby/MCArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/UniversalJoint.h>
#include <Moby/BulirschStoerIntegrator.h>
#include <Moby/RungeKuttaIntegrator.h>
#include <Moby/RungeKuttaFehlbergIntegrator.h>
#include <Moby/RungeKuttaImplicitIntegrator.h>
#include <Moby/ODEPACKIntegrator.h>
#include <Moby/EulerIntegrator.h>
#include <Moby/VariableEulerIntegrator.h>
#include <Moby/GravityForce.h>
#include <Moby/StokesDragForce.h>
#include <Moby/DampingForce.h>
#include <Moby/XMLTree.h>
#include <Moby/XMLReader.h>

using boost::shared_ptr;
using namespace Moby;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
std::map<std::string, BasePtr> XMLReader::read(const std::string& fname)
{
  // setup the list of IDs
  std::map<std::string, BasePtr> id_map;
  
  // *************************************************************
  // going to remove any path from the argument and change to that
  // path; this is done so that all files referenced from the
  // local path of the XML file are found
  // *************************************************************

  // set the filename to use as the argument, by default
  std::string filename = fname;

  // get the current pathname
  size_t BUFSIZE = 8192;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "XMLReader::read() - unable to allocate sufficient memory!" << std::endl;
      return id_map;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != std::string::npos)
  {
    // get the new working path
    std::string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,std::string::npos);
  }

  // read the XML Tree 
  shared_ptr<const XMLTree> root_tree = XMLTree::read_from_xml(filename);
  if (!root_tree)
  {
    std::cerr << "XMLReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return id_map;
  }

 
  // find the moby tree 
  shared_ptr<XMLTree> moby_tree = boost::const_pointer_cast<XMLTree>(find_subtree(root_tree, "moby"));

  // mark the moby root as processed
  moby_tree->processed = true;

   // make sure that the Moby node was found
  if (!moby_tree)
  {
    std::cerr << "XMLReader::read() - no moby tag found!" << std::endl;
    chdir(cwd.get());
    return id_map;
  }

  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RigidBody is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // read and construct all primitives
  process_tag("Box", moby_tree, &read_box, id_map);
  process_tag("Sphere", moby_tree, &read_sphere, id_map);
  process_tag("Cylinder", moby_tree, &read_cylinder, id_map);
  process_tag("Cone", moby_tree, &read_cone, id_map);
  process_tag("Heightmap", moby_tree, &read_heightmap, id_map);
  process_tag("Plane", moby_tree, &read_plane, id_map);
/*
  process_tag("TriangleMesh", moby_tree, &read_trimesh, id_map);
  process_tag("TetraMesh", moby_tree, &read_tetramesh, id_map);
  process_tag("PrimitivePlugin", moby_tree, &read_primitive_plugin, id_map);
  process_tag("CSG", moby_tree, &read_CSG, id_map);
*/

  // read and construct all integrators
  process_tag("EulerIntegrator", moby_tree, &read_euler_integrator, id_map);
  process_tag("VariableEulerIntegrator", moby_tree, &read_variable_euler_integrator, id_map);
  process_tag("BulirschStoerIntegrator", moby_tree, &read_bulirsch_stoer_integrator, id_map);
  process_tag("RungeKuttaIntegrator", moby_tree, &read_rk4_integrator, id_map);
  process_tag("RungeKuttaFehlbergIntegrator", moby_tree, &read_rkf4_integrator, id_map);
  process_tag("RungeKuttaImplicitIntegrator", moby_tree, &read_rk4i_integrator, id_map);
  process_tag("ODEPACKIntegrator", moby_tree, &read_odepack_integrator, id_map);

  // read and construct all recurrent forces (except damping)
  process_tag("GravityForce", moby_tree, &read_gravity_force, id_map);
  process_tag("StokesDragForce", moby_tree, &read_stokes_drag_force, id_map);

  #ifdef USE_OSG
  // read and construct all OSGGroupWrapper objects
  process_tag("OSGGroup", moby_tree, &read_osg_group, id_map);
  #endif

  // read and construct all rigid bodies (including articulated body links)
  process_tag("RigidBody", moby_tree, &read_rigid_body, id_map);

  // read and construct all joints -- we do this after the links have been read
  process_tag("RevoluteJoint", moby_tree, &read_revolute_joint, id_map);
  process_tag("PrismaticJoint", moby_tree, &read_prismatic_joint, id_map);
  process_tag("SphericalJoint", moby_tree, &read_spherical_joint, id_map);
  process_tag("UniversalJoint", moby_tree, &read_universal_joint, id_map);
  process_tag("FixedJoint", moby_tree, &read_fixed_joint, id_map);
  process_tag("JointPlugin", moby_tree, &read_joint_plugin, id_map);

  // read and construct all articulated bodies
//  process_tag("MCArticulatedBody", moby_tree, &read_mc_abody, id_map);
  process_tag("RCArticulatedBody", moby_tree, &read_rc_abody, id_map);
  process_tag("RCArticulatedBodySymbolicPlugin", moby_tree, &read_rc_abody_symbolic, id_map);

  // damping forces must be constructed after bodies
  process_tag("DampingForce", moby_tree, &read_damping_force, id_map);

  // finally, read and construct the simulator objects -- must be done last
  process_tag("Simulator", moby_tree, &read_simulator, id_map);
  process_tag("EventDrivenSimulator", moby_tree, &read_event_driven_simulator, id_map);
  process_tag("TimeSteppingSimulator", moby_tree, &read_time_stepping_simulator, id_map);

  // change back to the initial working directory
  chdir(cwd.get());

  // output unprocessed tags / attributes
  std::queue<shared_ptr<const XMLTree> > q;
  q.push(moby_tree);
  while (!q.empty())
  {
    // get the node off the front of the queue
    shared_ptr<const XMLTree> node = q.front();
    q.pop();

    // check whether the tag was processed
    if (!node->processed)
    {
      std::cerr << "XMLReader::read() warning- tag '" << node->name << "' not processed" << std::endl;
      continue;
    }

    // verify that all attributes were processed
    BOOST_FOREACH(const XMLAttrib& a, node->attribs)
      if (!a.processed)
        std::cerr << "XMLReader::read() warning- attribute '" << a.name << "' in tag '" << node->name << "' not processed" << std::endl;

    // add all children to the queue
    BOOST_FOREACH(XMLTreePtr child, node->children)
      q.push(child);
  }

  return id_map;
}

/// Finds and processes given tags
void XMLReader::process_tag(const std::string& tag, shared_ptr<const XMLTree> root, void (*fn)(shared_ptr<const XMLTree>, std::map<std::string, BasePtr>&), std::map<std::string, BasePtr>& id_map)
{
  // NOTE: if a tag is encountered, we do not process its descendants: 
  // load_from_xml() is responsible for that

  // if this node is of the given type, process it 
  if (strcasecmp(root->name.c_str(), tag.c_str()) == 0)
    fn(root, id_map);
  else
  {
    const std::list<XMLTreePtr>& child_nodes = root->children;
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      (*i)->processed = true;
      process_tag(tag, *i, fn, id_map);
    }
  }
}

/// Reads and constructs a geometry plugin object
void XMLReader::read_primitive_plugin(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // get the name of the plugin to load
  XMLAttrib* plugin_attr = node->get_attrib("plugin");
  if (!plugin_attr)
  {
    std::cerr << "XMLReader::read_primitive_plugin() - no plugin attribute!" << std::endl;
    return;
  }
  std::string pluginname = plugin_attr->get_string_value();

  // verify that the plugin can be found
  struct stat filestatus;
  if (stat(pluginname.c_str(), &filestatus) != 0)
  {
    std::cerr << "XMLReader::read_primitive_plugin() - unable to find plugin '" << pluginname << "'" << std::endl;
    return;
  }

  // load the plugin
  void* plugin = dlopen(pluginname.c_str(), RTLD_LAZY);
  if (!plugin)
  {
    std::cerr << "XMLReader::read_primitive_plugin()- cannot load plugin: " << dlerror() << std::endl;
    return;
  }

  // load the factory symbol
  boost::shared_ptr<Primitive> (*factory)(void) = (boost::shared_ptr<Primitive> (*) (void)) dlsym(plugin, "factory");
  if (!factory)
  {
    std::cerr << "XMLReader::read_primitive_plugin()- factory() not found in " << pluginname << std::endl;
    return;
  }

  // create a new Base object
  boost::shared_ptr<Base> primitive_plugin = factory();
  
  // populate the object
  primitive_plugin->load_from_xml(node, id_map);
}

/// Reads and constructs the OSGGroupWrapper object
void XMLReader::read_osg_group(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "OSGGroup") == 0);

  #ifdef USE_OSG
  // create a new OSGGroupWrapper object
  OSGGroupWrapperPtr group(new OSGGroupWrapper());

  // populate the object
  group->load_from_xml(node, id_map);
  #endif
}

/// Reads and constructs the SpherePrimitive object
void XMLReader::read_sphere(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Sphere") == 0);

  // create a new Base object
  boost::shared_ptr<Base> b(new SpherePrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the ConePrimitive object
void XMLReader::read_cone(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Cone") == 0);

  // create a new ConePrimitive object
  boost::shared_ptr<Base> b(new ConePrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the CylinderPrimitive object
void XMLReader::read_cylinder(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Cylinder") == 0);

  // create a new CylinderPrimitive object
  boost::shared_ptr<Base> b(new CylinderPrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs a CSG object
void XMLReader::read_CSG(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "CSG") == 0);

  // create a new CSG object
//  boost::shared_ptr<Base> b(new CSG());

  // populate the object
//  b->load_from_xml(node, id_map);
}

/// Reads and constructs the TriangleMeshPrimitive object
void XMLReader::read_tetramesh(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "TetraMesh") == 0);

  // create a new IndexedTetraArray object
//  boost::shared_ptr<Base> b(new IndexedTetraArray());
  
  // populate the object
//  b->load_from_xml(node, id_map);
}

/// Reads and constructs the TriangleMeshPrimitive object
void XMLReader::read_trimesh(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "TriangleMesh") == 0);

  // create a new TriangleMeshPrimitive object
//  boost::shared_ptr<Base> b(new TriangleMeshPrimitive());
  
  // populate the object
//  b->load_from_xml(node, id_map);
}

/// Reads and constructs a plane object
void XMLReader::read_plane(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Plane") == 0);

  // create a new plane object
  boost::shared_ptr<Base> b(new PlanePrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs a heightmap object
void XMLReader::read_heightmap(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Heightmap") == 0);

  // create a new Heightmap object
  boost::shared_ptr<Base> b(new HeightmapPrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the BoxPrimitive object
void XMLReader::read_box(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Box") == 0);

  // create a new BoxPrimitive object
  boost::shared_ptr<Base> b(new BoxPrimitive());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the EulerIntegrator object
/**
 * \pre node name is EulerIntegrator
 */
void XMLReader::read_euler_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "EulerIntegrator") == 0);

  // only create VectorN integrators
  boost::shared_ptr<Base> b(new EulerIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the VariableEulerIntegrator object
/**
 * \pre node name is VariableEulerIntegrator
 */
void XMLReader::read_variable_euler_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "VariableEulerIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new VariableEulerIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and construct the RungaKuttaFehlbergIntegrator object
/**
 * \pre node name is BulirschStoerIntegrator
 */
void XMLReader::read_bulirsch_stoer_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "BulirschStoerIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new BulirschStoerIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and construct the RungaKuttaFehlbergIntegrator object
/**
 * \pre node name is RungeKuttaFehlbergIntegrator
 */
void XMLReader::read_rkf4_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RungeKuttaFehlbergIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new RungeKuttaFehlbergIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and construct the RungaKuttaIntegrator object
/**
 * \pre node name is RungeKuttaIntegrator
 */
void XMLReader::read_rk4_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RungeKuttaIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new RungeKuttaIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and construct the RungaKuttaImplicitIntegrator object
/**
 * \pre node name is RungeKuttaImplicitIntegrator
 */
void XMLReader::read_rk4i_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RungeKuttaImplicitIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new RungeKuttaImplicitIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and construct the ODEPACKIntegrator object
/**
 * \pre node name is ODEPACKIntegrator
 */
void XMLReader::read_odepack_integrator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "ODEPACKIntegrator") == 0);

  // only create VectorN type integrators
  boost::shared_ptr<Base> b(new ODEPACKIntegrator());

  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the TimeSteppingSimulator object
/**
 * \pre node is named TimeSteppingSimulator
 */
void XMLReader::read_time_stepping_simulator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "TimeSteppingSimulator") == 0);

  // create a new TimeSteppingSimulator object
  boost::shared_ptr<Base> b(new TimeSteppingSimulator());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the EventDrivenSimulator object
/**
 * \pre node is named EventDrivenSimulator
 */
void XMLReader::read_event_driven_simulator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "EventDrivenSimulator") == 0);

  // create a new EventDrivenSimulator object
  boost::shared_ptr<Base> b(new EventDrivenSimulator());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the Simulator object
/**
 * \pre node is named Simulator
 */
void XMLReader::read_simulator(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Simulator") == 0);

  // create a new Simulator object
  boost::shared_ptr<Base> b(new Simulator());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the RigidBody object
/**
 * \pre node is named RigidBody
 */
void XMLReader::read_rigid_body(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RigidBody") == 0);

  // create a new RigidBody object
  boost::shared_ptr<Base> b(new RigidBody());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the MCArticulatedBody object
/**
 * \pre node is named MCArticulatedBody
 */
void XMLReader::read_mc_abody(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "MCArticulatedBody") == 0);

  // create a new MCArticulatedBody object
//  boost::shared_ptr<MCArticulatedBody> link(new MCArticulatedBody());
  
  // populate the object
//  link->load_from_xml(node, id_map);
}

/// Reads and constructs the RCArticulatedBody object
/**
 * \pre node is named RCArticulatedBody
 */
void XMLReader::read_rc_abody(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RCArticulatedBody") == 0);

  // create a new RCArticulatedBody object
  boost::shared_ptr<RCArticulatedBody> link(new RCArticulatedBody());
  
  // populate the object
  link->load_from_xml(node, id_map);
}

/// Reads and constructs a symbolic RCArticulatedBody object
/**
 * \pre node is named RCArticulatedBodySymbolic
 */
void XMLReader::read_rc_abody_symbolic(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RCArticulatedBodySymbolicPlugin") == 0);

  // get the name of the plugin to load
  XMLAttrib* plugin_attr = node->get_attrib("plugin");
  if (!plugin_attr)
  {
    std::cerr << "XMLReader::read_rc_abody_symbolic() - no plugin attribute!" << std::endl;
    return;
  }
  std::string pluginname = plugin_attr->get_string_value();

  // verify that the plugin can be found
  struct stat filestatus;
  if (stat(pluginname.c_str(), &filestatus) != 0)
  {
    std::cerr << "XMLReader::read_rc_abody_symbolic() - unable to find plugin '" << pluginname << "'" << std::endl;
    return;
  }

  // load the plugin
  void* plugin = dlopen(pluginname.c_str(), RTLD_NOW);
  if (!plugin)
  {
    std::cerr << "XMLReader::read_rc_abody_symbolic()- cannot load plugin: " << dlerror() << std::endl;
    return;
  }

  // load the factory symbol
  boost::shared_ptr<RCArticulatedBody> (*factory)(void) = (boost::shared_ptr<RCArticulatedBody> (*) (void)) dlsym(plugin, "factory");
  if (!factory)
  {
    std::cerr << "XMLReader::read_rc_abody_symbolic()- factory() not found found in " << pluginname << std::endl;
    return;
  }

  // create a new RCArticulatedBody object
  boost::shared_ptr<RCArticulatedBody> body = factory();
  
  // populate the object
  body->load_from_xml(node, id_map);
}

/// Reads and constructs a plugin Joint object
/**
 * \pre node is named JointPlugin 
 */
void XMLReader::read_joint_plugin(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "JointPlugin") == 0);

  // get the name of the plugin to load
  XMLAttrib* plugin_attr = node->get_attrib("plugin");
  if (!plugin_attr)
  {
    std::cerr << "XMLReader::read_joint_plugin() - no plugin attribute!" << std::endl;
    return;
  }
  std::string pluginname = plugin_attr->get_string_value();

  // verify that the plugin can be found
  struct stat filestatus;
  if (stat(pluginname.c_str(), &filestatus) != 0)
  {
    std::cerr << "XMLReader::read_joint_plugin() - unable to find plugin '" << pluginname << "'" << std::endl;
    return;
  }

  // load the plugin
  void* plugin = dlopen(pluginname.c_str(), RTLD_LAZY);
  if (!plugin)
  {
    std::cerr << "XMLReader::read_joint_plugin()- cannot load plugin: " << dlerror() << std::endl;
    return;
  }

  // load the factory symbol
  boost::shared_ptr<Joint> (*factory)(void) = (boost::shared_ptr<Joint> (*) (void)) dlsym(plugin, "factory");
  if (!factory)
  {
    std::cerr << "XMLReader::read_joint_plugin()- factory() not found in " << pluginname << std::endl;
    return;
  }

  // create a new Joint object
  boost::shared_ptr<Joint> joint_plugin = factory();
  
  // populate the object
  joint_plugin->load_from_xml(node, id_map);
}

/// Reads and constructs the UniversalJoint object
void XMLReader::read_universal_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "UniversalJoint") == 0);

  // create a new UniversalJoint object
  boost::shared_ptr<UniversalJoint> uj(new UniversalJoint());
  
  // populate the object
  uj->load_from_xml(node, id_map);
}

/// Reads and constructs the SphericalJoint object
void XMLReader::read_spherical_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "SphericalJoint") == 0);

  // create a new SphericalJoint object
  boost::shared_ptr<SphericalJoint> sj(new SphericalJoint());
  
  // populate the object
  sj->load_from_xml(node, id_map);
}

/// Reads and constructs the FixedJoint object
void XMLReader::read_fixed_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "FixedJoint") == 0);

  // create a new FixedJoint object
  boost::shared_ptr<FixedJoint> fj(new FixedJoint());
  
  // populate the object
  fj->load_from_xml(node, id_map);
}

/// Reads and constructs the RevoluteJoint object
void XMLReader::read_revolute_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RevoluteJoint") == 0);

  // create a new RevoluteJoint object
  boost::shared_ptr<RevoluteJoint> rj(new RevoluteJoint());
  
  // populate the object
  rj->load_from_xml(node, id_map);
}

/// Reads and constructs the PrismaticJoint object
void XMLReader::read_prismatic_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "PrismaticJoint") == 0);

  // create a new RevoluteJoint object
  boost::shared_ptr<Base> b(new PrismaticJoint());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the StokesDragForce object
/**
 * \pre node is named StokesDragForce
 */
void XMLReader::read_stokes_drag_force(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "StokesDragForce") == 0);

  // create a new StokesDragForce object
  boost::shared_ptr<StokesDragForce> sdf(new StokesDragForce());
  
  // populate the object
  sdf->load_from_xml(node, id_map);
}

/// Reads and constructs the DampingForce object
/**
 * \pre node is named DampingForce
 */
void XMLReader::read_damping_force(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "DampingForce") == 0);

  // create a new DampingForce object
  boost::shared_ptr<DampingForce> df(new DampingForce());
  
  // populate the object
  df->load_from_xml(node, id_map);
}

/// Reads and constructs the GravityForce object
/**
 * \pre node is named GravityForce
 */
void XMLReader::read_gravity_force(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "GravityForce") == 0);

  // create a new GravityForce object
  boost::shared_ptr<GravityForce> gf(new GravityForce());
  
  // populate the object
  gf->load_from_xml(node, id_map);
}

/// Gets the sub-tree rooted at the specified tag
shared_ptr<const XMLTree> XMLReader::find_subtree(shared_ptr<const XMLTree> root, const std::string& name)
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

/// Gets the tuple type from a node
XMLReader::TupleType XMLReader::get_tuple(shared_ptr<const XMLTree> node)
{
  std::string type;

  // get the 'type' attribute
  XMLAttrib* type_attr = node->get_attrib("type");
  if (type_attr)
    type = type_attr->get_string_value();

  // look for possible tuple types
  if (strcasecmp(type.c_str(), "VectorN") == 0)
    return eVectorN;
  if (strcasecmp(type.c_str(), "Vector3") == 0)
    return eVector3;
  if (strcasecmp(type.c_str(), "Quat") == 0)
    return eQuat;

  // still here?  not one of the above...
  return eNone;
}

