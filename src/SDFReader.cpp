/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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

using std::vector;
using std::list;
using boost::shared_ptr;
using namespace Moby;
using namespace Ravelin;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
void SDFReader::read(const std::string& fname, std::vector<std::vector<DynamicBodyPtr> >& models)
{
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
      return;
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
    return;
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
    return;
  }

  // read in all world tags
  std::list<shared_ptr<const XMLTree> > world_nodes = find_tag("world", sdf_tree);

  // process each world tag
  BOOST_FOREACH(shared_ptr<const XMLTree> world_node, world_nodes)
  {
    // create a vector of models for this world
    models.push_back(vector<DynamicBodyPtr>());

    // get all model nodes
    std::list<shared_ptr<const XMLTree> > model_nodes = find_tag("model", world_node); 
    BOOST_FOREACH(shared_ptr<const XMLTree> model_node, model_nodes)
      models.back().push_back(read_model(model_node));
  }

  // change back to the initial working directory
  chdir(cwd.get());
}

/// Find a particular tag (top of the recursive function)
list<shared_ptr<const XMLTree> > SDFReader::find_tag(const std::string& tag, shared_ptr<const XMLTree> root)
{
  // create the list
  list<shared_ptr<const XMLTree> > l;

  // if this node is of the given type, process it 
  if (strcasecmp(root->name.c_str(), tag.c_str()) == 0)
  {
    l.push_back(root);
    return l;
  }
  else
  {
    const std::list<XMLTreePtr>& child_nodes = root->children;
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
      find_tag(tag, *i, l);
  }
}

/// Find a particular tag (recursive function)
void SDFReader::find_tag(const std::string& tag, shared_ptr<const XMLTree> root, list<shared_ptr<const XMLTree> >& l)
{
  // if this node is of the given type, process it and attempt to go no further 
  if (strcasecmp(root->name.c_str(), tag.c_str()) == 0)
  {
    l.push_back(root);
    return;
  }
  else
  {
    const std::list<XMLTreePtr>& child_nodes = root->children;
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
      find_tag(tag, *i, l);
  }
}

/// Find a particular tag (recursive function)
shared_ptr<const XMLTree> SDFReader::find_one_tag(const std::string& tag, shared_ptr<const XMLTree> root)
{
  // if this node is of the given type, process it and attempt to go no further 
  if (strcasecmp(root->name.c_str(), tag.c_str()) == 0)
    return root;
  else
  {
    const std::list<XMLTreePtr>& child_nodes = root->children;
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      shared_ptr<const XMLTree> node = find_one_tag(tag, *i);
      if (node)
        return node;
    } 
  }

  return shared_ptr<const XMLTree>();
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
  if (val == "false")
    return false;
  else if (val == "true")
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
JointPtr SDFReader::read_joint(shared_ptr<const XMLTree> node, const std::map<std::string, RigidBodyPtr>& link_map)
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

  // read in the name of the parent link
  shared_ptr<const XMLTree> parent_tag = find_one_tag("parent", node);
  assert(parent_tag);
  std::string parent_link = parent_tag->content;
  assert(link_map.find(parent_link) != link_map.end());
  parent = link_map.find(parent_link)->second;

  // read in the name of the child link  
  shared_ptr<const XMLTree> child_tag = find_one_tag("child", node);
  assert(child_tag);
  std::string child_link = child_tag->content; 
  assert(link_map.find(child_link) != link_map.end());
  child = link_map.find(child_link)->second;

  // set child and parent
  joint->set_inboard_link(parent);
  joint->set_outboard_link(child);

  // read the pose (offset from child link to joint frame) in child link frame
  if (find_one_tag("pose", node))
  {
    *P = read_pose(node);
    P->rpose = child->get_pose();
    joint->set_pose(P); 
  }

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
    shared_ptr<const XMLTree> limit_node = find_one_tag("limits", axis_node);
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

  // setup the generic components of the joint
  joint->id = name;

  return joint;
}

/// Reads and constructs the SpherePrimitive object
PrimitivePtr SDFReader::read_sphere(shared_ptr<const XMLTree> node)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Sphere") == 0);

  // create a new Base object
  boost::shared_ptr<SpherePrimitive> s(new SpherePrimitive());

  // get the radius attribute
  XMLAttrib* radius_attr = node->get_attrib("radius");
  if (radius_attr)
    s->set_radius(radius_attr->get_real_value());

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

  // set values for the object
  if (radius_attr && len_attr)
  {
    c->set_radius(radius_attr->get_real_value());
    c->set_height(len_attr->get_real_value());
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
  XMLAttrib* normal_attr = node->get_attrib("normal");
  if (normal_attr)
  {
    Vector3d normal;
    normal_attr->get_vector_value(normal);

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
  XMLAttrib* size_attr = node->get_attrib("size");

  // get the lengths
  Vector3d len;
  if (size_attr) 
  {
    size_attr->get_vector_value(len);
    b->set_size(len[X], len[Y], len[Z]);
  }
  
  return b;
}

/// Reads and constructs a pointer to a DynamicBody object from a Model tag
/**
 * \pre node is named Model 
 */
DynamicBodyPtr SDFReader::read_model(shared_ptr<const XMLTree> node)
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
    RigidBodyPtr rb = read_link(link_nodes.front()); 

    // set the name
    if (name_attr)
      rb->id = name_attr->get_string_value();

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
      links.push_back(read_link(link_node));

    // construct a mapping from link id's to links
    std::map<std::string, RigidBodyPtr> link_map;
    for (unsigned i=0; i< links.size(); i++)
      link_map[links[i]->id] = links[i];

    // read all of the joints
    BOOST_FOREACH(shared_ptr<const XMLTree> joint_node, joint_nodes)
      joints.push_back(read_joint(joint_node, link_map));

    // set the links and joints
    rcab->set_links_and_joints(links, joints);

    // set the name
    if (name_attr)
      rcab->id = name_attr->get_string_value();

    // transform the body if desired
    if (pose_node)
    {
      // read the pose
      Pose3d P = read_pose(node);

      // apply the rotation first
      rcab->rotate(P.q);

      // apply the translation
      rcab->translate(P.x);
    }
 
    return rcab;
  }
}

/// Reads and constructs a RigidBody object from a Link tag
/**
 * \pre node is named Link 
 */
RigidBodyPtr SDFReader::read_link(shared_ptr<const XMLTree> node)
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
  shared_ptr<const XMLTree> inertia_node = find_one_tag("inertia", node);
  if (inertia_node)
    rb->set_inertia(read_inertia(inertia_node, rb));

  // read the Collision tag
  shared_ptr<const XMLTree> collision_node = find_one_tag("collision", node);
  if (collision_node)
    read_collision_node(collision_node, rb);   
}

/// Reads a collision node
void SDFReader::read_collision_node(shared_ptr<const XMLTree> node, RigidBodyPtr rb)
{
  // setup a collision geometry
  CollisionGeometryPtr cg(new CollisionGeometry);

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

  // add the collision geometry to the rigid body
  rb->geometries.push_back(cg);
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
SpatialRBInertiad SDFReader::read_inertia(shared_ptr<const XMLTree> node, RigidBodyPtr rb)
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
    rb->set_inertial_pose(read_pose(node));

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


