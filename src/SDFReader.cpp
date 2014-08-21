/****************************************************************************
 * Copyright 2014 Evan Drumwright
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

// TODO:
// implement triangle mesh

#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/Constants.h>
#include <Moby/Simulator.h>
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RigidBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
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
#include <Moby/SDFReader.h>

using std::map;
using std::vector;
using std::list;
using boost::shared_ptr;
using namespace Moby;
using namespace Ravelin;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
shared_ptr<EventDrivenSimulator> SDFReader::read(const std::string& fname)
{
  vector<vector<DynamicBodyPtr> > models;

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
      std::cerr << "SDFReader::read() - unable to allocate sufficient memory!" << std::endl;
      return shared_ptr<EventDrivenSimulator>();
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
    std::cerr << "SDFReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return shared_ptr<EventDrivenSimulator>();
  }

  // find the SDF tree 
  shared_ptr<XMLTree> sdf_tree = boost::const_pointer_cast<XMLTree>(find_subtree(root_tree, "SDF"));

  // mark the root as processed
  sdf_tree->processed = true;

  // make sure that the SDF node was found
  if (!sdf_tree)
  {
    std::cerr << "SDFReader::read() - no SDF tag found!" << std::endl;
    chdir(cwd.get());
    return shared_ptr<EventDrivenSimulator>();
  }

  // read in all world tags
  std::list<shared_ptr<const XMLTree> > world_nodes = find_tag("world", sdf_tree);

  // read in worlds
  if (world_nodes.size() != 1)
    throw std::runtime_error("SDFReader::read() - there is not exactly one world!");
  shared_ptr<EventDrivenSimulator> sim = read_world(world_nodes.front()); 

  // change back to the initial working directory
  chdir(cwd.get());

  return sim;
}

/// Constructs the event-driven simulator using proper settings
shared_ptr<EventDrivenSimulator> SDFReader::read_world(shared_ptr<const XMLTree> world_tree)
{
  // create the simulator
  shared_ptr<EventDrivenSimulator> sim(new EventDrivenSimulator); 

  // read the models
  vector<DynamicBodyPtr> models = read_models(world_tree, sim);

  // these defaults will be replaced with specific settings from SDF
  sim->integrator = shared_ptr<BulirschStoerIntegrator>(new BulirschStoerIntegrator);
  sim->rel_err_tol = 1e-3;
  sim->abs_err_tol = 1e-3;
  sim->minimum_step = 1e-5;

  // find the physics node
  shared_ptr<const XMLTree> physics_node = find_one_tag("physics", world_tree);
  if (physics_node)
  {
/*
    // setup and read the maximum step size, if present
    shared_ptr<const XMLTree> max_step_sz_node = find_one_tag("max_step_size", physics_node);
    if (max_step_sz_node)
      STEP_SIZE = read_double(max_step_sz_node);
*/
    // read the gravity vector
    shared_ptr<const XMLTree> gravity_node = find_one_tag("gravity", physics_node);
    if (gravity_node)
    {
      // create and add the gravity force to all bodies
      shared_ptr<GravityForce> grav(new GravityForce);
      BOOST_FOREACH(DynamicBodyPtr db, models)
        db->get_recurrent_forces().push_back(grav); 

      // set the force
      grav->gravity = read_Vector3(gravity_node);
    } 

    // read the Moby tag
    shared_ptr<const XMLTree> moby_node = find_one_tag("moby", physics_node);
    if (moby_node)
    {
      // read the integrator node
      shared_ptr<const XMLTree> int_node = find_one_tag("integrator", moby_node);
      if (int_node)
      {
        // set a pointer to a variable step integrator
        shared_ptr<VariableStepIntegrator> vsi;
        XMLAttrib* type_attr = int_node->get_attrib("type");
        if (strcasecmp(type_attr->get_string_value().c_str(), "BulirschStoer") == 0)
        {
          shared_ptr<BulirschStoerIntegrator> bsi(new BulirschStoerIntegrator);
          sim->integrator = bsi;
          vsi = bsi;
        }
        else if (strcasecmp(type_attr->get_string_value().c_str(), "RKF") == 0)
        {
          shared_ptr<RungeKuttaFehlbergIntegrator> rkf(new RungeKuttaFehlbergIntegrator);
          sim->integrator = rkf;
          vsi = rkf;
        }
        else if (strcasecmp(type_attr->get_string_value().c_str(), "ODEPACK") == 0)
        {
          shared_ptr<ODEPACKIntegrator> ode(new ODEPACKIntegrator);
          sim->integrator = ode;
          vsi = ode;
        }
        else if (strcasecmp(type_attr->get_string_value().c_str(), "rk4") == 0)
        {
          shared_ptr<RungeKuttaIntegrator> rk4(new RungeKuttaIntegrator);
          sim->integrator = rk4;
        }
        else if (strcasecmp(type_attr->get_string_value().c_str(), "rki") == 0)
        {
          shared_ptr<RungeKuttaImplicitIntegrator> rki(new RungeKuttaImplicitIntegrator);
          sim->integrator = rki;
        }

        // read the minimum step size
        shared_ptr<const XMLTree> min_step_node = find_one_tag("min_step_size", moby_node);
        if (min_step_node && vsi)
          vsi->min_step_size = read_double(min_step_node);

        // read the absolute error tolerance
        shared_ptr<const XMLTree> ae_tol_node = find_one_tag("abs_err_tol", moby_node);
        if (ae_tol_node && vsi)
          vsi->aerr_tolerance = read_double(ae_tol_node);

        // read the relative error tolerance
        shared_ptr<const XMLTree> re_tol_node = find_one_tag("rel_err_tol", moby_node);
        if (re_tol_node && vsi)
          vsi->rerr_tolerance = read_double(re_tol_node);
      }
    } 
  }

  return sim;
}

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
vector<DynamicBodyPtr> SDFReader::read_models(shared_ptr<const XMLTree> world_tree, shared_ptr<EventDrivenSimulator> sim)
{
  vector<DynamicBodyPtr> models;
  map<RigidBodyPtr, shared_ptr<SurfaceData> > sdata;

  // get all model nodes
  std::list<shared_ptr<const XMLTree> > model_nodes = find_tag("model", world_tree); 
  BOOST_FOREACH(shared_ptr<const XMLTree> model_node, model_nodes)
    models.push_back(read_model(model_node, sdata));

  // add the models to the simulator
  BOOST_FOREACH(DynamicBodyPtr b, models)
    sim->add_dynamic_body(b); 

  // now attempt to add contact data
  for (map<RigidBodyPtr, shared_ptr<SurfaceData> >::const_iterator i = sdata.begin(); i != sdata.end(); i++)
  {
    map<RigidBodyPtr, shared_ptr<SurfaceData> >::const_iterator j = i;
    for (j++; j != sdata.end(); j++)
    {
      // create a ContactParameters object
      shared_ptr<ContactParameters> cp(new ContactParameters);

      // setup the restitution
      cp->epsilon = (i->second->epsilon + j->second->epsilon)*0.5;

      // setup the number of friction directions
      cp->NK = std::max(i->second->NK, j->second->NK);

      // setup the friction coefficients
      cp->mu_coulomb = (i->second->mu_c + j->second->mu_c)*0.5;
      cp->mu_viscous = (i->second->mu_v + j->second->mu_v)*0.5;

      // add the contact data
      sim->contact_params[make_sorted_pair(i->first, j->first)] = cp;  
    }
  }

  return models;
}

/// Find a particular tag (top of the recursive function)
list<shared_ptr<const XMLTree> > SDFReader::find_tag(const std::string& tag, shared_ptr<const XMLTree> root)
{
  // create the list
  list<shared_ptr<const XMLTree> > l;

  find_tag(tag, root, l);
  return l;
}

/// Find a particular tag
void SDFReader::find_tag(const std::string& tag, shared_ptr<const XMLTree> root, list<shared_ptr<const XMLTree> >& l)
{
  // process all children of the root
  const std::list<XMLTreePtr>& child_nodes = root->children;
  for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    if (strcasecmp((*i)->name.c_str(), tag.c_str()) == 0)
      l.push_back(*i);
}

/// Find a particular tag
shared_ptr<const XMLTree> SDFReader::find_one_tag(const std::string& tag, shared_ptr<const XMLTree> root)
{
  const std::list<XMLTreePtr>& child_nodes = root->children;
  for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    if (strcasecmp((*i)->name.c_str(), tag.c_str()) == 0)
      return *i;

  return shared_ptr<const XMLTree>();
}

/// Reads an unsigned int value
unsigned SDFReader::read_uint(shared_ptr<const XMLTree> node)
{
  // convert the string to a uint 
  return (unsigned) std::atoi(node->content.c_str()); 
}

/// Reads a double value
double SDFReader::read_double(shared_ptr<const XMLTree> node)
{
  // convert the string to a double
  return std::atof(node->content.c_str()); 
}

/// Reads a Boolean value
bool SDFReader::read_bool(shared_ptr<const XMLTree> node)
{
  std::string val = node->content;
  std::transform(val.begin(), val.end(), val.begin(), ::tolower);
  if (val == "false" || val == "0")
    return false;
  else if (val == "true" || val == "1")
    return true;
  else
    throw std::runtime_error("SDFReader::read_bool() - value is not a Boolean");
}

/// Reads a Vector3 value
Vector3d SDFReader::read_Vector3(shared_ptr<const XMLTree> node)
{
  VectorNd w = VectorNd::parse(node->content);
  if (w.size() != 3)
    throw MissizeException();
  return Vector3d(w[0], w[1], w[2]);
}

/// Reads a joint
JointPtr SDFReader::read_joint(shared_ptr<const XMLTree> node, const std::map<std::string, RigidBodyPtr>& link_map, RigidBodyPtr& base_link)
{
  std::string name, type;
  JointPtr joint;
  shared_ptr<Pose3d> P(new Pose3d);
  RigidBodyPtr parent, child;
  shared_ptr<RevoluteJoint> rj;
  shared_ptr<PrismaticJoint> pj;
  shared_ptr<UniversalJoint> uj;
  shared_ptr<SphericalJoint> sj;

  // get the id of the joint
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
    name = name_attr->get_string_value();

  // get the joint type
  XMLAttrib* type_attr = node->get_attrib("type");
  if (type_attr)
  {
    type = type_attr->get_string_value();
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "revolute")
    {
      rj = shared_ptr<RevoluteJoint>(new RevoluteJoint);
      joint = rj;
    }
    else if (type == "prismatic")
    {
      pj = shared_ptr<PrismaticJoint>(new PrismaticJoint);
      joint = pj;
    }
    else if (type == "universal")
    {
      uj = shared_ptr<UniversalJoint>(new UniversalJoint);
      joint = uj;
    }
    else if (type == "ball")
    {
      sj = shared_ptr<SphericalJoint>(new SphericalJoint);
      joint = sj;
    }
    else if (type == "revolute2")
      throw std::runtime_error("revolute2 joint unsupported");
    else if (type == "gearbox")
      throw std::runtime_error("gearbox joint unsupported");
    else if (type == "piston")
      throw std::runtime_error("piston joint unsupported");
    else
    {
      std::string err = type + " joint unsupported";
      throw std::runtime_error(err.c_str());
    }
  }

  // setup the generic components of the joint
  joint->id = name;

  // read in the name of the parent link
  shared_ptr<const XMLTree> parent_tag = find_one_tag("parent", node);
  assert(parent_tag);
  std::string parent_link = parent_tag->content;
  if (strcasecmp(parent_link.c_str(), "world") == 0)
  {
    // joint is attached to the world; create the body, if not already done
    if (!base_link)
    {
      base_link = RigidBodyPtr(new RigidBody);
      base_link->set_enabled(false);
    }
    parent = base_link;
  }
  else if (link_map.find(parent_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_joint(.)- parent link '" + parent_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  else 
    parent = link_map.find(parent_link)->second;

  // read in the name of the child link  
  shared_ptr<const XMLTree> child_tag = find_one_tag("child", node);
  assert(child_tag);
  std::string child_link = child_tag->content; 
  if (link_map.find(child_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_joint(.)- child link '" + child_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  child = link_map.find(child_link)->second;

  // set child and parent
  joint->set_inboard_link(parent, false);
  joint->set_outboard_link(child, false);

  // read the axis tag (contains limits, joint damping/friction)
  shared_ptr<const XMLTree> axis_node = find_one_tag("axis", node);
  if (axis_node)
  {
    // read whether to use the parent model frame
    shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis_node);

    // read the axis, if any
    shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis_node);
    if (xyz_node)
    {
      // get the axis
      Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose 
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P; 

      // set the axis
      if (rj)
        rj->set_axis(axis);
      else if (pj)
        pj->set_axis(axis);
      else if (uj)
        uj->set_axis(axis, UniversalJoint::eAxis1);
      else if (sj)
        sj->set_axis(axis, SphericalJoint::eAxis1);
    }

    // read viscous and Coulomb friction
    shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        joint->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        joint->mu_fv = read_double(damping_node);
    }

    // read joint limits
    shared_ptr<const XMLTree> limit_node = find_one_tag("limit", axis_node);
    if (limit_node)
    {
      // attempt to read the lower limit 
      shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        joint->lolimit[0] = read_double(llimit_node);

      // attempt to read the upper limit
      shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        joint->hilimit[0] = read_double(ulimit_node);

      // attempt to read the maximum force
      shared_ptr<const XMLTree> effort_node = find_one_tag("effort", limit_node);
      if (effort_node)
        joint->maxforce[0] = read_double(effort_node);
    }
  }

  // read the axis tag (contains limits, joint damping/friction)
  shared_ptr<const XMLTree> axis2_node = find_one_tag("axis2", node);
  if (axis2_node)
  {
    // read whether to use the parent model frame
    shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis2_node);

    // read the axis, if any
    shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis2_node);
    if (xyz_node)
    {
      // get the axis
      Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose 
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P; 

      // set the axis
      if (rj)
        rj->set_axis(axis);
      else if (pj)
        pj->set_axis(axis);
      else if (uj)
        uj->set_axis(axis, UniversalJoint::eAxis2);
      else if (sj)
        sj->set_axis(axis, SphericalJoint::eAxis2);
    }

    // read viscous and Coulomb friction
    shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis2_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        joint->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        joint->mu_fv = read_double(damping_node);
    }

    // read joint limits
    shared_ptr<const XMLTree> limit_node = find_one_tag("limits", axis2_node);
    if (limit_node)
    {
      // attempt to read the lower limit 
      shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        joint->lolimit[1] = read_double(llimit_node);

      // attempt to read the upper limit
      shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        joint->hilimit[1] = read_double(ulimit_node);

      // attempt to read the maximum force
      shared_ptr<const XMLTree> effort_node = find_one_tag("effort", limit_node);
      if (effort_node)
        joint->maxforce[1] = read_double(effort_node);
    }
  }

  return joint;
}

/// Reads and constructs the SpherePrimitive object
PrimitivePtr SDFReader::read_sphere(shared_ptr<const XMLTree> node)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Sphere") == 0);

  // create a new Base object
  boost::shared_ptr<SpherePrimitive> s(new SpherePrimitive());

  // get the radius 
  shared_ptr<const XMLTree> radius_node = find_one_tag("radius", node);
  if (radius_node)
    s->set_radius(read_double(radius_node));

  return s;  
}

/// Reads and constructs the CylinderPrimitive object
PrimitivePtr SDFReader::read_cylinder(shared_ptr<const XMLTree> node)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Cylinder") == 0);

  // create a new CylinderPrimitive object
  boost::shared_ptr<CylinderPrimitive> c(new CylinderPrimitive());
  
  // get the length and radius attributes
  XMLAttrib* radius_attr = node->get_attrib("radius");
  XMLAttrib* len_attr = node->get_attrib("length");

  // get the length and radius 
  shared_ptr<const XMLTree> radius_node = find_one_tag("radius", node);
  shared_ptr<const XMLTree> length_node = find_one_tag("length", node);

  // set values for the object
  if (radius_node && length_node)
  {
    c->set_radius(read_double(radius_node));
    c->set_height(read_double(length_node));
  }

  return c;
}

/// Reads and constructs the Plane object
PrimitivePtr SDFReader::read_plane(shared_ptr<const XMLTree> node)
{  
  const unsigned X = 0, Y = 1, Z = 2;

  // sanity check
  assert(strcasecmp(node->name.c_str(), "Plane") == 0);

  // create a new PlanePrimitive object
  boost::shared_ptr<PlanePrimitive> b(new PlanePrimitive());

  // get the plane normal here
  shared_ptr<const XMLTree> normal_node = find_one_tag("normal", node);
  if (normal_node)
  {
    Vector3d normal = read_Vector3(normal_node);

    // determine *a* rotation matrix that converts from [0 1 0] to the normal
    Vector3d tan1, tan2;
    Vector3d::determine_orthonormal_basis(normal, tan1, tan2);
    Matrix3d R;
    R.set_column(X, tan1);
    R.set_column(Y, normal);
    R.set_column(Z, -tan2);

    // setup the pose for the primitive
    Pose3d P;
    P.x.set_zero();
    P.q = R;
    b->set_pose(P);
  }

  return b;      
}

/// Reads and constructs the TriangleMeshPrimitive object
PrimitivePtr SDFReader::read_trimesh(shared_ptr<const XMLTree> node)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "TriangleMesh") == 0);

  // TODO: finish implementing this
  // create a new TriangleMeshPrimitive object
  boost::shared_ptr<TriangleMeshPrimitive> b(new TriangleMeshPrimitive());
  
  // populate the object
//  b->load_from_xml(node, id_map);
}

/// Reads and constructs the heightmap object
PrimitivePtr SDFReader::read_heightmap(shared_ptr<const XMLTree> node)
{  
  // sanity check
  assert(strcasecmp(node->name.c_str(), "heightmap") == 0);

  // TODO: finish implementing this
  // create a new HeightmapPrimitive object
  boost::shared_ptr<HeightmapPrimitive> b(new HeightmapPrimitive());
  
  return b;
}

/// Reads and constructs the BoxPrimitive object
PrimitivePtr SDFReader::read_box(shared_ptr<const XMLTree> node)
{  
  const double X = 0, Y = 1, Z = 2;

  // sanity check
  assert(strcasecmp(node->name.c_str(), "box") == 0);

  // create a new BoxPrimitive object
  boost::shared_ptr<BoxPrimitive> b(new BoxPrimitive());

  // get the length attributes, if specified
  shared_ptr<const XMLTree> size_node = find_one_tag("size", node);

  // get the lengths
  if (size_node) 
  {
    Vector3d len = read_Vector3(size_node);
    b->set_size(len[X], len[Y], len[Z]);
  }
  
  return b;
}

/// Reads and constructs a pointer to a DynamicBody object from a Model tag
/**
 * \pre node is named Model 
 */
DynamicBodyPtr SDFReader::read_model(shared_ptr<const XMLTree> node, map<RigidBodyPtr, shared_ptr<SDFReader::SurfaceData> >& sdata)
{
  vector<RigidBodyPtr> links;
  vector<JointPtr> joints;

  // sanity check
  assert(strcasecmp(node->name.c_str(), "Model") == 0);

  // get the model name
  XMLAttrib* name_attr = node->get_attrib("name");

  // read links and joints
  std::list<shared_ptr<const XMLTree> > link_nodes = find_tag("link", node);
  std::list<shared_ptr<const XMLTree> > joint_nodes = find_tag("joint", node);

  // get the pose for the body, if specified
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);

  // see whether there is an individual rigid body
  if (link_nodes.size() == 1 && joint_nodes.empty())
  {
    // create the rigid body
    shared_ptr<SurfaceData> sd;
    RigidBodyPtr rb = read_link(link_nodes.front(), sd); 
    if (sd)
      sdata[rb] = sd;    

    // set the name
    if (name_attr)
      rb->id = name_attr->get_string_value();

    // see whether the model is static 
    shared_ptr<const XMLTree> static_node = find_one_tag("static", node);
    if (static_node && read_bool(static_node))
      rb->set_enabled(false);

    // transform the body if desired
    if (pose_node)
    {
      // read the pose
      Pose3d P = read_pose(node);

      // apply the rotation first
      rb->rotate(P.q);

      // apply the translation
      rb->translate(P.x);
    }

    return rb;
  }
  else
  {
    // create an articulated body
    RCArticulatedBodyPtr rcab(new RCArticulatedBody);

    // read all of the links
    BOOST_FOREACH(shared_ptr<const XMLTree> link_node, link_nodes)
    {
      shared_ptr<SurfaceData> sd;
      links.push_back(read_link(link_node, sd));
      if (sd)
        sdata[links.back()] = sd;
    }

    // transform each rigid body if desired
    if (pose_node)
    {
      // read the pose
      Pose3d P = read_pose(node);

      BOOST_FOREACH(RigidBodyPtr link, links)
      {
        // apply the rotation first
        link->rotate(P.q);

        // apply the translation
        link->translate(P.x);
      }
    }

    // construct a mapping from link id's to links
    std::map<std::string, RigidBodyPtr> link_map;
    for (unsigned i=0; i< links.size(); i++)
      link_map[links[i]->id] = links[i];

    // read all of the joints
    RigidBodyPtr base_link;
    BOOST_FOREACH(shared_ptr<const XMLTree> joint_node, joint_nodes)
      joints.push_back(read_joint(joint_node, link_map, base_link));

    // if the base link is created, add it to the vector and make the rcab
    // have a fixed base
    if (base_link)
    {
      links.insert(links.begin(), base_link);
      rcab->set_floating_base(false);
    }

    // set the links and joints
    rcab->set_links_and_joints(links, joints);

    // set the name
    if (name_attr)
      rcab->id = name_attr->get_string_value();

 
    return rcab;
  }
}

/// Reads and constructs a RigidBody object from a Link tag
/**
 * \pre node is named Link 
 */
RigidBodyPtr SDFReader::read_link(shared_ptr<const XMLTree> node, shared_ptr<SDFReader::SurfaceData>& sd)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Link") == 0);

  // create a new RigidBody object
  boost::shared_ptr<RigidBody> rb(new RigidBody());

  // get the link name
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
    rb->id = name_attr->get_string_value();
  
  // get the pose attribute for the body, if specified
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
    rb->set_pose(read_pose(node));

  // get the inertial properties for the body, if specified
  shared_ptr<const XMLTree> inertial_node = find_one_tag("inertial", node);
  if (inertial_node)
    rb->set_inertia(read_inertial(inertial_node, rb));

  // read the Visual tag
  shared_ptr<const XMLTree> visual_node = find_one_tag("visual", node);
  if (visual_node)
    read_visual_node(visual_node, rb); 
    
  std::list<shared_ptr<const XMLTree> > collision_nodes = find_tag("collision", node);
  
  // read the Collision tag/tags
  if (!collision_nodes.empty())
  {    
     if (collision_nodes.size() == 1){
       read_collision_node(collision_nodes.front(), rb, sd);
     }
     else{
       BOOST_FOREACH( shared_ptr<const XMLTree> collision_node, collision_nodes){
         read_collision_node(collision_node, rb, sd);
       }
     }
  } 

  return rb;
}

/// Reads a visual node
void SDFReader::read_visual_node(shared_ptr<const XMLTree> node, RigidBodyPtr rb)
{
  // read the pose of the visualization
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
  {
    Pose3d P(read_pose(node));
    P.rpose = rb->get_visualization_pose()->rpose;
    rb->set_visualization_relative_pose(P);
  }

  // read the geometry type
  shared_ptr<const XMLTree> geom_node = find_one_tag("geometry", node);
  if (geom_node)
  {
    PrimitivePtr geom = read_geometry(geom_node);
    rb->set_visualization_data(geom->get_visualization());
  }
}

/// Reads a collision node
void SDFReader::read_collision_node(shared_ptr<const XMLTree> node, RigidBodyPtr rb, shared_ptr<SDFReader::SurfaceData>& sd)
{
  // setup a collision geometry
  CollisionGeometryPtr cg(new CollisionGeometry);

  // set the rigid body and add the geometry
  rb->geometries.push_back(cg);
  cg->set_single_body(rb);

  // get the link name
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
    cg->id = name_attr->get_string_value();

   // read the pose of the collision geometry
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
  {
    Pose3d P(read_pose(node));
    P.update_relative_pose(rb->get_pose());
    cg->set_relative_pose(P);
  }

  // read the geometry type
  shared_ptr<const XMLTree> geom_node = find_one_tag("geometry", node);
  if (geom_node)
    cg->set_geometry(read_geometry(geom_node));

  // read the surface data, if any
  read_surface(find_one_tag("surface", node), sd);
}

/// Reads the surface data
void SDFReader::read_surface(shared_ptr<const XMLTree> node, shared_ptr<SDFReader::SurfaceData>& sd)
{
  // setup the SurfaceData
  sd = shared_ptr<SurfaceData>(new SurfaceData);

  // setup reasonable defaults for surface data
  sd->epsilon = 0.0;
  sd->NK = 4;
  sd->mu_c = 0.2;
  sd->mu_v = 0.0;

  // exit if node is null
  if (!node)
    return;

  // read the bounce (coefficient of restitution), if specified
  shared_ptr<const XMLTree> bounce_node = find_one_tag("bounce", node);
  if (bounce_node)
  {
    shared_ptr<const XMLTree> cor_node = find_one_tag("restitution_coefficient", bounce_node);
    if (cor_node)
      sd->epsilon = read_double(cor_node);
  }

  // read the friction node, if specified
  shared_ptr<const XMLTree> friction_node = find_one_tag("friction", node);
  if (friction_node)
  {
    shared_ptr<const XMLTree> moby_node = find_one_tag("moby", friction_node);
    if (moby_node)
    {
      // attempt to read mu_Coulomb, mu_viscous, and # of friction cone edges
      shared_ptr<const XMLTree> muc_node = find_one_tag("mu_coulomb", moby_node);
      shared_ptr<const XMLTree> muv_node = find_one_tag("mu_viscous", moby_node);
      shared_ptr<const XMLTree> nk_node = find_one_tag("num_friction_edges", moby_node);

      // set appropriate parts of surface data
      if (muc_node) sd->mu_c = read_double(muc_node);
      if (muv_node) sd->mu_v = read_double(muv_node);
      if (nk_node) sd->NK = read_uint(nk_node);
    }
  }
}

/// Reads geometry
PrimitivePtr SDFReader::read_geometry(shared_ptr<const XMLTree> node)
{
  // look for a box
  shared_ptr<const XMLTree> box_node = find_one_tag("box", node);
  if (box_node)
    return read_box(box_node);

  // look for a cylinder
  shared_ptr<const XMLTree> cylinder_node = find_one_tag("cylinder", node);
  if (cylinder_node)
    return read_cylinder(cylinder_node);

  // look for a sphere
  shared_ptr<const XMLTree> sphere_node = find_one_tag("sphere", node);
  if (sphere_node)
    return read_sphere(sphere_node);

  // look for a heightmap
  shared_ptr<const XMLTree> heightmap_node = find_one_tag("heightmap", node);
  if (heightmap_node)
    return read_heightmap(heightmap_node);

  // look for a triangle mesh
  shared_ptr<const XMLTree> trimesh_node = find_one_tag("mesh", node);
  if (trimesh_node)
    return read_trimesh(trimesh_node);

  // look for a plane
  shared_ptr<const XMLTree> plane_node = find_one_tag("plane", node);
  if (plane_node)
    return read_plane(plane_node);

  // shouldn't still be here...
  throw std::runtime_error("Geometry tag found that we couldn't handle!"); 

  return PrimitivePtr();
}

/// Reads a pose
Pose3d SDFReader::read_pose(shared_ptr<const XMLTree> node)
{
  // look for the pose tag
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);

  // get the pose
  VectorNd pose = VectorNd::parse(pose_node->content);

  // setup the pose
  Pose3d P;
  P.x = Origin3d(pose[0], pose[1], pose[2]);
  P.q = Quatd::rpy(pose[3], pose[4], pose[5]);
  return P;
}

/// Reads the inertia from the inertial node
SpatialRBInertiad SDFReader::read_inertial(shared_ptr<const XMLTree> node, RigidBodyPtr rb)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the spatial rigid body inertia
  SpatialRBInertiad J;

  // get the mass
  shared_ptr<const XMLTree> mass_node = find_one_tag("mass", node);
  if (mass_node)
    J.m = read_double(mass_node); 

  // get the pose of the inertial frame and set it with respect to the link
  // reference frame
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
  {
    Pose3d P = read_pose(node);
    P.rpose = rb->get_pose();
    rb->set_inertial_pose(P);
  }

  // find the inertia
  shared_ptr<const XMLTree> inertia_node = find_one_tag("inertia", node);
  assert(inertia_node);  

  // get the xx inertia  
  shared_ptr<const XMLTree> ixx_node = find_one_tag("ixx", inertia_node);
  if (ixx_node)
    J.J(X,X) = read_double(ixx_node); 

  // get the xy inertia  
  shared_ptr<const XMLTree> ixy_node = find_one_tag("ixy", inertia_node);
  if (ixy_node)
    J.J(Y,X) = J.J(X,Y) = read_double(ixy_node); 

  // get the xz inertia  
  shared_ptr<const XMLTree> ixz_node = find_one_tag("ixz", inertia_node);
  if (ixz_node)
    J.J(Z,X) = J.J(X,Z) = read_double(ixz_node); 

  // get the yy inertia  
  shared_ptr<const XMLTree> iyy_node = find_one_tag("iyy", inertia_node);
  if (iyy_node)
    J.J(Y,Y) = read_double(iyy_node); 

  // get the yz inertia  
  shared_ptr<const XMLTree> iyz_node = find_one_tag("iyz", inertia_node);
  if (iyz_node)
    J.J(Y,Z) = J.J(Z,Y) = read_double(iyz_node); 

  // get the zz inertia  
  shared_ptr<const XMLTree> izz_node = find_one_tag("izz", inertia_node);
  if (izz_node)
    J.J(Z,Z) = read_double(izz_node); 

  return J;
}

/// Gets the sub-tree rooted at the specified tag
shared_ptr<const XMLTree> SDFReader::find_subtree(shared_ptr<const XMLTree> root, const std::string& name)
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


