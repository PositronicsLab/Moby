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

using std::list;
using boost::shared_ptr;
using namespace Moby;
using namespace Ravelin;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
std::map<std::string, BasePtr> SDFReader::read(const std::string& fname)
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
      std::cerr << "SDFReader::read() - unable to allocate sufficient memory!" << std::endl;
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
    std::cerr << "SDFReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return id_map;
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
    return id_map;
  }

  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RigidBody is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // TODO: finish this

  // read and construct all rigid bodies (including articulated body links)
  //process_tag("Link", sdf_tree, &read_rigid_body, id_map);

  // change back to the initial working directory
  chdir(cwd.get());

  // output unprocessed tags / attributes
  std::queue<shared_ptr<const XMLTree> > q;
  q.push(sdf_tree);
  while (!q.empty())
  {
    // get the node off the front of the queue
    shared_ptr<const XMLTree> node = q.front();
    q.pop();

    // check whether the tag was processed
    if (!node->processed)
    {
      std::cerr << "SDFReader::read() warning- tag '" << node->name << "' not processed" << std::endl;
      continue;
    }

    // verify that all attributes were processed
    BOOST_FOREACH(const XMLAttrib& a, node->attribs)
      if (!a.processed)
        std::cerr << "SDFReader::read() warning- attribute '" << a.name << "' in tag '" << node->name << "' not processed" << std::endl;

    // add all children to the queue
    BOOST_FOREACH(XMLTreePtr child, node->children)
      q.push(child);
  }

  return id_map;
}

/// Finds and processes given tags
void SDFReader::process_tag(const std::string& tag, shared_ptr<const XMLTree> root, void (*fn)(shared_ptr<const XMLTree>, std::map<std::string, BasePtr>&), std::map<std::string, BasePtr>& id_map)
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


/// Reads and constructs the OSGGroupWrapper object
void SDFReader::read_osg_group(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
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
  XMLAttrib* parent_link_attr = parent_tag->get_attrib("link");
  if (parent_link_attr)
  {
    std::string parent_link = parent_link_attr->get_string_value();
    assert(link_map.find(parent_link) != link_map.end());
    parent = link_map.find(parent_link)->second;
  }

  // read in the name of the child link  
  shared_ptr<const XMLTree> child_tag = find_one_tag("child", node);
  assert(child_tag);
  XMLAttrib* child_link_attr = parent_tag->get_attrib("link");
  if (child_link_attr)
  {
    std::string child_link = child_link_attr->get_string_value();
    assert(link_map.find(child_link) != link_map.end());
    child = link_map.find(child_link)->second;
  }

  // read the pose (offset from child link to joint frame) in child link frame
  *P = read_pose(node);
  P->rpose = child->get_pose();
  joint->set_pose(P); 

  // TODO: set child and parent

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

/// Reads and constructs a RigidBody object from a Link tag
/**
 * \pre node is named Link 
 */
RigidBodyPtr SDFReader::read_rigid_body(shared_ptr<const XMLTree> node)
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
  shared_ptr<const XMLTree> pose_node = find_one_tag("origin", node);
  if (pose_node)
    rb->set_pose(read_pose(pose_node));

  // get the inertial properties for the body, if specified
  shared_ptr<const XMLTree> inertia_node = find_one_tag("inertia", node);
  if (inertia_node)
    rb->set_inertia(read_inertia(node, rb));

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
  shared_ptr<const XMLTree> pose_node = find_one_tag("origin", node);
  if (pose_node)
  {
    Pose3d P(read_pose(pose_node));
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
  XMLAttrib* pose_attrib = node->get_attrib("pose");
  assert(pose_attrib);

  // get the pose
  VectorNd pose;
  pose_attrib->get_vector_value(pose);

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
  XMLAttrib* mass_attrib = node->get_attrib("mass");
  if (mass_attrib)
    J.m = mass_attrib->get_real_value();

  // get the pose of the inertial frame and set it with respect to the link
  // reference frame
  shared_ptr<const XMLTree> origin_tag = find_one_tag("origin", node);
  if (origin_tag)
    rb->set_inertial_pose(read_pose(origin_tag));

  // find the inertia
  shared_ptr<const XMLTree> inertia_node = find_one_tag("inertia", node);
  assert(inertia_node);  

  // get the xx inertia  
  XMLAttrib* xx_attrib = inertia_node->get_attrib("ixx");
  if (xx_attrib)
    J.J(X,X) = xx_attrib->get_real_value();

  // get the xy inertia  
  XMLAttrib* xy_attrib = inertia_node->get_attrib("ixy");
  if (xy_attrib)
    J.J(Y,X) = J.J(X,Y) = xy_attrib->get_real_value();

  // get the xz inertia  
  XMLAttrib* xz_attrib = inertia_node->get_attrib("ixz");
  if (xz_attrib)
    J.J(Z,X) = J.J(X,Z) = xz_attrib->get_real_value();

  // get the yy inertia  
  XMLAttrib* yy_attrib = inertia_node->get_attrib("iyy");
  if (yy_attrib)
    J.J(Y,Y) = yy_attrib->get_real_value();

  // get the yz inertia  
  XMLAttrib* yz_attrib = inertia_node->get_attrib("iyz");
  if (yz_attrib)
    J.J(Y,Z) = J.J(Z,Y) = yz_attrib->get_real_value();

  // get the zz inertia  
  XMLAttrib* zz_attrib = inertia_node->get_attrib("izz");
  if (zz_attrib)
    J.J(Z,Z) = zz_attrib->get_real_value();

  return J;
}





/// Reads and constructs the MCArticulatedBody object
/**
 * \pre node is named MCArticulatedBody
 */
void SDFReader::read_mc_abody(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
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
void SDFReader::read_rc_abody(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RCArticulatedBody") == 0);

  // create a new RCArticulatedBody object
  boost::shared_ptr<RCArticulatedBody> link(new RCArticulatedBody());
  
  // populate the object
  link->load_from_xml(node, id_map);
}

/// Reads and constructs the UniversalJoint object
void SDFReader::read_universal_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "UniversalJoint") == 0);

  // create a new UniversalJoint object
  boost::shared_ptr<UniversalJoint> uj(new UniversalJoint());
  
  // populate the object
  uj->load_from_xml(node, id_map);
}

/// Reads and constructs the SphericalJoint object
void SDFReader::read_spherical_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "SphericalJoint") == 0);

  // create a new SphericalJoint object
  boost::shared_ptr<SphericalJoint> sj(new SphericalJoint());
  
  // populate the object
  sj->load_from_xml(node, id_map);
}

/// Reads and constructs the FixedJoint object
void SDFReader::read_fixed_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "FixedJoint") == 0);

  // create a new FixedJoint object
  boost::shared_ptr<FixedJoint> fj(new FixedJoint());
  
  // populate the object
  fj->load_from_xml(node, id_map);
}

/// Reads and constructs the RevoluteJoint object
void SDFReader::read_revolute_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "RevoluteJoint") == 0);

  // create a new RevoluteJoint object
  boost::shared_ptr<RevoluteJoint> rj(new RevoluteJoint());
  
  // populate the object
  rj->load_from_xml(node, id_map);
}

/// Reads and constructs the PrismaticJoint object
void SDFReader::read_prismatic_joint(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "PrismaticJoint") == 0);

  // create a new RevoluteJoint object
  boost::shared_ptr<Base> b(new PrismaticJoint());
  
  // populate the object
  b->load_from_xml(node, id_map);
}

/// Reads and constructs the DampingForce object
/**
 * \pre node is named DampingForce
 */
void SDFReader::read_damping_force(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "DampingForce") == 0);

  // create a new DampingForce object
  boost::shared_ptr<DampingForce> df(new DampingForce());
  
  // populate the object
  df->load_from_xml(node, id_map);
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

/// Gets the tuple type from a node
SDFReader::TupleType SDFReader::get_tuple(shared_ptr<const XMLTree> node)
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

