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
#include <boost/foreach.hpp>

#ifdef USE_OSG
#include <Moby/OSGGroupWrapper.h>
#include <osg/MatrixTransform>
#include <osg/Matrixd>
#include "Color.h"
#endif

#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/Constants.h>
#include <Moby/Simulator.h>
#include <Moby/RigidBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/FixedJoint.h>
//#include <Moby/MCArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/UniversalJoint.h>
#include <Moby/GravityForce.h>
#include <Moby/StokesDragForce.h>
#include <Moby/DampingForce.h>
#include <Moby/XMLTree.h>
#include <Moby/RigidBody.h>
#include <Moby/SDFReader.h>

using std::map;
using std::vector;
using std::list;
using boost::shared_ptr;
using namespace Moby;
using namespace Ravelin;

std::vector<shared_ptr<OSGGroupWrapper> > SDFReader::_osg_wrappers;
std::vector<shared_ptr<Primitive> > SDFReader::_primitives;

#ifdef USE_OSG
/// Copies this matrix to an OpenSceneGraph Matrixd object
static void to_osg_matrix(const Pose3d& src, osg::Matrixd& tgt)
{
  // get the rotation matrix
  Matrix3d M = src.q;

  // setup the rotation components of tgt
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = M(i,j);

  // setup the translation components of tgt
  for (unsigned i=X; i<= Z; i++)
    tgt(W,i) = src.x[i];

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (double) 0.0;
  tgt(W,W) = (double) 1.0;
}
#endif

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
shared_ptr<TimeSteppingSimulator> SDFReader::read(const std::string& fname)
{
  vector<vector<ControlledBodyPtr> > models;

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
      return shared_ptr<TimeSteppingSimulator>();
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
    return shared_ptr<TimeSteppingSimulator>();
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
    return shared_ptr<TimeSteppingSimulator>();
  }

  // read in all world tags
  std::list<shared_ptr<const XMLTree> > world_nodes = find_tag("world", sdf_tree);

  // read in worlds
  if (world_nodes.size() != 1)
    throw std::runtime_error("SDFReader::read() - there is not exactly one world!");
  shared_ptr<TimeSteppingSimulator> sim = read_world(world_nodes.front());

  // change back to the initial working directory
  chdir(cwd.get());

  return sim;
}

/// Reads models only from SDF file
/**
 * \return a map of IDs to read objects
 */
std::map<std::string, ControlledBodyPtr> SDFReader::read_models(const std::string& fname)
{
  std::map<std::string, ControlledBodyPtr> model_map;
  vector<ControlledBodyPtr> models;

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
      std::cerr << "SDFReader::read_model() - unable to allocate sufficient memory!" << std::endl;
      return model_map;
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
    std::cerr << "SDFReader::read_model() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return model_map;
  }

  // find the SDF tree
  shared_ptr<XMLTree> sdf_tree = boost::const_pointer_cast<XMLTree>(find_subtree(root_tree, "SDF"));

  // mark the root as processed
  sdf_tree->processed = true;

  // make sure that the SDF node was found
  if (!sdf_tree)
  {
    std::cerr << "SDFReader::read_model() - no SDF tag found!" << std::endl;
    chdir(cwd.get());
    return model_map;
  }

  // read in all world tags
  std::list<shared_ptr<const XMLTree> > world_nodes = find_tag("world", sdf_tree);

  // read in worlds
  if (world_nodes.size() > 1)
    throw std::runtime_error("SDFReader::read() - more than one world found!");

  // create the simulator
  shared_ptr<TimeSteppingSimulator> sim(new TimeSteppingSimulator);

  // read the models
  if (world_nodes.empty())
    models = read_models(sdf_tree, sim);
  else
    models = read_models(world_nodes.front(), sim);

  // clear all models from the simulator and convert to a map
  const vector<ControlledBodyPtr>& bodies = sim->get_dynamic_bodies();
  BOOST_FOREACH(ControlledBodyPtr db, models)
  {
    sim->remove_dynamic_body(db);
    model_map[db->id] = db;
  }

  // change back to the initial working directory
  chdir(cwd.get());

  return model_map;
}

/// Constructs the event-driven simulator using proper settings
shared_ptr<TimeSteppingSimulator> SDFReader::read_world(shared_ptr<const XMLTree> world_tree)
{
  // create the simulator
  shared_ptr<TimeSteppingSimulator> sim(new TimeSteppingSimulator);

  // read the models
  vector<ControlledBodyPtr> models = read_models(world_tree, sim);

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
      BOOST_FOREACH(ControlledBodyPtr db, models)
        db->get_recurrent_forces().push_back(grav);

      // set the force
      grav->gravity = read_Vector3(gravity_node);
    }

    // read the Moby tag
    shared_ptr<const XMLTree> moby_node = find_one_tag("moby", physics_node);
    if (moby_node)
    {
    }
  }

  return sim;
}

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
vector<ControlledBodyPtr> SDFReader::read_models(shared_ptr<const XMLTree> world_tree, shared_ptr<TimeSteppingSimulator> sim)
{
  vector<ControlledBodyPtr> models;
  map<RigidBodyPtr, shared_ptr<SurfaceData> > sdata;

  // get all model nodes
  std::list<shared_ptr<const XMLTree> > model_nodes = find_tag("model", world_tree);
  BOOST_FOREACH(shared_ptr<const XMLTree> model_node, model_nodes)
    models.push_back(read_model(model_node, sdata));

  // add the models to the simulator
  BOOST_FOREACH(ControlledBodyPtr b, models)
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

/// Reads a string value
std::string SDFReader::read_string(shared_ptr<const XMLTree> node)
{
  // convert the string to a uint
  return node->content;
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

  // set the pose for the cylinder
  c->set_pose(Pose3d(Matrix3d::rot_X(M_PI_2)));

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

/// Reads and constructs an OSGGroupWrapper object for meshes
shared_ptr<OSGGroupWrapper> SDFReader::read_OSG_file(shared_ptr<const XMLTree> node)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "mesh") == 0);

  // get the filename
  shared_ptr<const XMLTree> uri_node = find_one_tag("uri", node);
  if (!uri_node)
    throw std::runtime_error("Expected a 'uri' subnode under 'mesh' node");
 
  // construct the OSGGroupWrapper
  std::string fname = read_string(uri_node);
  shared_ptr<OSGGroupWrapper> osgg(new OSGGroupWrapper(fname));
  return osgg; 
}

/// Reads and constructs the PolyhedralPrimitive object
PrimitivePtr SDFReader::read_polyhedron(shared_ptr<const XMLTree> node)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "mesh") == 0);

  // get the filename
  shared_ptr<const XMLTree> uri_node = find_one_tag("uri", node);
  if (!uri_node)
    throw std::runtime_error("Expected a 'uri' subnode under 'mesh' node");
  
  // ensure that the file is a Wavefront OBJ
  std::string fname = read_string(uri_node);
  std::string fname_lower = fname;
  std::transform(fname.begin(), fname.end(), fname_lower.begin(), ::tolower);
  unsigned st = fname_lower.find(".obj");
  if (st != fname_lower.size()-4)
    throw std::runtime_error("Expect 'uri' to be of type Wavefront .obj");

  // read in the file using an indexed triangle array
  IndexedTriArray ita = IndexedTriArray::read_from_obj(fname);

  // get all of the vertices and compute the convex hull (yielding a
  // tessellated polyhedron)
  const std::vector<Origin3d>& vertices = ita.get_vertices();
  TessellatedPolyhedronPtr tessellated_poly = CompGeom::calc_convex_hull(vertices.begin(), vertices.end());   

  // convert the tessellated polyhedron to a standard polyhedron
  Polyhedron poly;
  tessellated_poly->to_polyhedron(poly);

  // create the polyhedral primitive object 
  shared_ptr<PolyhedralPrimitive> p(new PolyhedralPrimitive);
  p->set_polyhedron(poly);

  return p;
}

/// Reads and constructs the trianglemesh primitive object
PrimitivePtr SDFReader::read_trimesh(shared_ptr<const XMLTree> node)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "mesh") == 0);

  // get the filename
  shared_ptr<const XMLTree> uri_node = find_one_tag("uri", node);
  if (!uri_node)
    throw std::runtime_error("Expected a 'uri' subnode under 'mesh' node");
  
  // ensure that the file is a Wavefront OBJ
  std::string fname = read_string(uri_node);
  std::string fname_lower = fname;
  std::transform(fname.begin(), fname.end(), fname_lower.begin(), ::tolower);
  unsigned st = fname_lower.find(".obj");
  if (st != fname_lower.size()-4)
    throw std::runtime_error("Expect 'uri' to be of type Wavefront .obj");

  // create a new TriangleMesh object
  boost::shared_ptr<TriangleMeshPrimitive> b(new TriangleMeshPrimitive(fname, false));

  return b;
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
ControlledBodyPtr SDFReader::read_model(shared_ptr<const XMLTree> node, map<RigidBodyPtr, shared_ptr<SDFReader::SurfaceData> >& sdata)
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
    {
      rb->id = name_attr->get_string_value();
      rb->body_id = rb->id;
    }

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

    // set algorithm and frame
    rcab->algorithm_type = RCArticulatedBody::eCRB;
    rcab->set_computation_frame_type(Ravelin::eLinkCOM);

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
      read_joint(joint_node, link_map, base_link, std::back_inserter(joints));

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
    {
      rcab->id = name_attr->get_string_value();
      rcab->body_id = rcab->id;
    }

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

  // set default inertial properties (according to SDF 1.5)
  SpatialRBInertiad J;
  J.m = 1.0;
  J.J.set_identity();

  // get the link name
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
  {
    rb->id = name_attr->get_string_value();
    rb->body_id = rb->id;
  }

  // get the pose attribute for the body, if specified
  shared_ptr<Pose3d> origin(new Pose3d);
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
    *origin = read_pose(node);

  // get the inertial properties for the body, if specified
  shared_ptr<Pose3d> inertial_pose(new Pose3d);
  shared_ptr<const XMLTree> inertial_node = find_one_tag("inertial", node);
  if (inertial_node)
    J = read_inertial(inertial_node, rb, inertial_pose);

  // read Visual tags
  std::list<shared_ptr<const XMLTree> > visual_nodes = find_tag("visual", node);
  BOOST_FOREACH(shared_ptr<const XMLTree> visual_node, visual_nodes)
    read_visual_node(visual_node, rb, origin);

  // read collision tags
  std::list<shared_ptr<const XMLTree> > collision_nodes = find_tag("collision", node);
  std::map<CollisionGeometryPtr, shared_ptr<Pose3d> > collision_poses;
  BOOST_FOREACH( shared_ptr<const XMLTree> collision_node, collision_nodes)
    read_collision_node(collision_node, rb, collision_poses, sd);

  // setup the inertial pose
  inertial_pose->rpose = origin;
  inertial_pose->update_relative_pose(rb->get_pose()->rpose);
  rb->set_pose(*origin);

  // setup the poses to point to the origin, and then update
  for (map<shared_ptr<CollisionGeometry>, shared_ptr<Pose3d> >::iterator i = collision_poses.begin(); i != collision_poses.end(); i++)
  {
    i->second->rpose = origin;
    i->second->update_relative_pose(rb->get_pose());
    i->first->set_relative_pose(*i->second);
  }

  // set inertia
  J.pose = rb->get_pose();
  rb->set_inertia(J);

  return rb;
}

/// Reads a visual node
void SDFReader::read_visual_node(shared_ptr<const XMLTree> node, RigidBodyPtr rb, shared_ptr<Pose3d> origin)
{
  // setup the pose
  Pose3d P;

  // read the pose of the visualization
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
    P = read_pose(node);

  // read the geometry type if OSG is used
#ifdef USE_OSG
  shared_ptr<const XMLTree> geom_node = find_one_tag("geometry", node);
  if (geom_node)
  {
    // create a new transform
    osg::MatrixTransform* tg = new osg::MatrixTransform;
    osg::Group* group = rb->get_visualization_data();
    group->addChild(tg);

    // add the primitive to the transform
    osg::Node* geom = read_visual_geometry(geom_node);

    shared_ptr<const XMLTree> material_node = find_one_tag("material", node);
    
    if(material_node){
      // Change the RGBA color of the link if provided
      shared_ptr<const XMLTree> color_node = find_one_tag("diffuse", material_node);
      if (color_node){
        // get the pose
        /// Color to add to the rigid body when rendered
        Ravelin::VectorNd color_rgba = Ravelin::VectorNd::parse(color_node->content);
        CcolorVisitor  newColor;
        newColor.setColor( color_rgba[0], color_rgba[1], color_rgba[2], color_rgba[3] );
        geom->accept( newColor );
      }
    }

    tg->addChild(geom);

    // set the transform
    P.update_relative_pose(GLOBAL);
    osg::Matrix m;
    to_osg_matrix(P, m);
    tg->setMatrix(m);
    }

  #endif
}

/// Reads a collision node
void SDFReader::read_collision_node(shared_ptr<const XMLTree> node, RigidBodyPtr rb, std::map<CollisionGeometryPtr, shared_ptr<Pose3d> >& collision_poses, shared_ptr<SDFReader::SurfaceData>& sd)
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
    shared_ptr<Pose3d> P(new Pose3d(read_pose(node)));
    collision_poses[cg] = P;
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

/// Reads visual geometry
osg::Node* SDFReader::read_visual_geometry(shared_ptr<const XMLTree> node)
{
  #ifdef USE_OSG
  // look for a box
  shared_ptr<const XMLTree> box_node = find_one_tag("box", node);
  if (box_node)
  {
    _primitives.push_back(read_box(box_node));
    return _primitives.back()->get_visualization();
  }

  // look for a cylinder
  shared_ptr<const XMLTree> cylinder_node = find_one_tag("cylinder", node);
  if (cylinder_node)
  {
    _primitives.push_back(read_cylinder(cylinder_node));
    return _primitives.back()->get_visualization();
  }

  // look for a sphere
  shared_ptr<const XMLTree> sphere_node = find_one_tag("sphere", node);
  if (sphere_node)
  {
    _primitives.push_back(read_sphere(sphere_node));
    return _primitives.back()->get_visualization();
  }

  // look for a heightmap
  shared_ptr<const XMLTree> heightmap_node = find_one_tag("heightmap", node);
  if (heightmap_node)
  {
    _primitives.push_back(read_heightmap(heightmap_node));
    return _primitives.back()->get_visualization();
  }

  // look for a triangle mesh
  shared_ptr<const XMLTree> trimesh_node = find_one_tag("mesh", node);
  if (trimesh_node)
  {
    _osg_wrappers.push_back(read_OSG_file(trimesh_node));
    return _osg_wrappers.back()->get_group();
  }

  // look for a plane
  shared_ptr<const XMLTree> plane_node = find_one_tag("plane", node);
  if (plane_node)
  {
    _primitives.push_back(read_plane(plane_node));
    return _primitives.back()->get_visualization();
  }

  // shouldn't still be here...
  throw std::runtime_error("Geometry tag found that we couldn't handle!");
  #endif
  return NULL;
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
    return read_polyhedron(trimesh_node);

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
SpatialRBInertiad SDFReader::read_inertial(shared_ptr<const XMLTree> node, RigidBodyPtr rb, shared_ptr<Pose3d>& inertial_pose)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the spatial rigid body inertia
  SpatialRBInertiad J = rb->get_inertia();

  // get the mass
  shared_ptr<const XMLTree> mass_node = find_one_tag("mass", node);
  if (mass_node)
    J.m = read_double(mass_node);

  // get the pose of the inertial frame and set it with respect to the link
  // reference frame
  shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  if (pose_node)
    *inertial_pose = read_pose(node);

  // find the inertia
  shared_ptr<const XMLTree> inertia_node = find_one_tag("inertia", node);
  if (inertia_node)
  {
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
  }

  // transform the inertia to the rigid body frame
  J = Pose3d::transform(rb->get_pose(), J);

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


