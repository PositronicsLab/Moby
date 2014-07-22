/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <dlfcn.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <stack>

#ifdef USE_OSG
#include <osg/Material>
#include <osg/Group>
#include <osg/Texture2D>
#include <osg/MatrixTransform>
#include <osgDB/ReadFile>
#endif

/*
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/CSG.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
*/
#include <Ravelin/AAngled.h>
#include <Moby/Primitive.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/FixedJoint.h>
#include <Moby/MCArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/UniversalJoint.h>
#include <Moby/XMLTree.h>
#include <Moby/URDFReader.h>

//#define DEBUG_URDF

using namespace Ravelin;
using namespace Moby;
using std::set;
using std::make_pair;
using std::queue;
using std::vector;
using std::pair;
using std::map;
using std::string;
using std::list;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
bool URDFReader::read(const string& fname, std::string& name, vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  // *************************************************************
  // going to remove any path from the argument and change to that
  // path; this is done so that all files referenced from the
  // local path of the XML file are found
  // *************************************************************

  // set the filename to use as the argument, by default
  string filename = fname;

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
      std::cerr << "URDFReader::read() - unable to allocate sufficient memory!" << std::endl;
      return false;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != string::npos)
  {
    // get the new working path
    string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,string::npos);
  }

  // read the XML Tree 
  shared_ptr<const XMLTree> tree = XMLTree::read_from_xml(filename);
  if (!tree)
  {
    std::cerr << "URDFReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return false;
  }
  
  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RigidBody is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // read and construct (single) robot 
  if (strcasecmp(tree->name.c_str(), "Robot") == 0)
  {
    URDFData data;
    if (!read_robot(tree, data, name, links, joints))
      return false;
  }
  else
  {
    std::cerr << "URDFReader::read() error - root element of URDF file is not a 'Robot' tag" << std::endl;
    return false;
  }

  // change back to the initial working directory
  chdir(cwd.get());

  return true;
}

/// Reads and constructs a robot object
bool URDFReader::read_robot(shared_ptr<const XMLTree> node, URDFData& data, string& name, vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  // read the robot name
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (name_attrib)
  {
    name = name_attrib->get_string_value();
    read_links(node, data, links);
    read_joints(node, data, links, joints);
  }
  else
  {
    std::cerr << "URDFReader::read_robot() - robot name not specified! not processing further..." << std::endl;
    return false;
  }

  return true;
}

/// Reads robot links
void URDFReader::read_links(shared_ptr<const XMLTree> node, URDFData& data, vector<RigidBodyPtr>& links)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_link(*i, data, links);
}

/// Reads robot joints 
void URDFReader::read_joints(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_joint(*i, data, links, joints);
}

/// Attempts to read a robot link from the given node
void URDFReader::read_link(shared_ptr<const XMLTree> node, URDFData& data, vector<RigidBodyPtr>& links)
{
  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Link") != 0)
    return;

  // link must have the name attribute
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (!name_attrib)
    std::cerr << "URDFReader::read_link() - link name not specified! not processing further..." << std::endl;

  // read and construct the link
  RigidBodyPtr link(new RigidBody);
  link->id = name_attrib->get_string_value();

  // read link properties
  read_inertial(node, data, link);
  read_visual(node, data, link);
  read_collision(node, data, link);

  // add the link to the set of links
  links.push_back(link);
}

/// Finds all children of the given link
void URDFReader::find_children(const URDFData& data, RigidBodyPtr link, queue<RigidBodyPtr>& q, map<RigidBodyPtr, RigidBodyPtr>& parents)
{
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
    if (i->second == link)
    {
      assert(data.joint_child.find(i->first) != data.joint_child.end());
      RigidBodyPtr child = data.joint_child.find(i->first)->second;
      q.push(child);
      parents[child] = link;
    } 
}

/// Finds the inner joint for a link
JointPtr URDFReader::find_joint(const URDFData& data, RigidBodyPtr outboard)
{
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    if (i->second == outboard)
      return i->first;

  return JointPtr();
}

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

/// Copies an OpenSceneGraph Matrixd object to this matrix 
static void from_osg_matrix(const osg::Matrixd& src, Pose3d& tgt)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  Matrix3d R;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      R(j,i) = src(i,j);

  tgt.q = R;
  tgt.x = Origin3d(src(W,X), src(W, Y), src(W, Z));
}
#endif

void URDFReader::find_outboards(const URDFData& data, RigidBodyPtr link, vector<pair<JointPtr, RigidBodyPtr> >& outboards, map<RigidBodyPtr, RigidBodyPtr>& parents)
{
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
    if (i->second == link)
    {
      assert(data.joint_child.find(i->first) != data.joint_child.end());
      RigidBodyPtr outboard = data.joint_child.find(i->first)->second;
      outboards.push_back(make_pair(i->first, outboard));
      parents[outboard] = link;
    } 
}

void URDFReader::output_data(const URDFData& data, RigidBodyPtr link)
{
  #ifdef DEBUG_URDF
  std::cout << "link id: " << link->id << std::endl;
  JointPtr joint = find_joint(data, link);
  if (joint)
  {
    std::cout << "  Moby joint position: " << joint->get_local() << std::endl;
    shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
    shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
    if (rj)
      std::cout << "  Moby joint axis: " << rj->get_axis() << std::endl;
    else if (pj)
      std::cout << "  Moby joint axis: " << pj->get_axis() << std::endl;
  }
  std::cout << "  Moby pose: " << std::endl << *link->get_pose();
  if (!link->geometries.empty())
  {
    std::cout << "  Moby collision pose: " << std::endl << *link->geometries.front()->get_pose();
  }
  #endif
}

/// Attempts to read a robot joint from the given node
void URDFReader::read_joint(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  JointPtr joint;
  RigidBodyPtr inboard, outboard;

  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Joint") != 0)
    return;

  // link must have the name attribute
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (!name_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint name not specified! not processing further..." << std::endl;
    return;
  }

  // link must have the name attribute
  XMLAttrib* type_attrib = node->get_attrib("type");
  if (!type_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint type not specified! not processing further..." << std::endl;
    return;
  }

  // read and construct the joint
  if (strcasecmp(type_attrib->get_string_value().c_str(), "revolute") == 0)
  {
    shared_ptr<RevoluteJoint> rj(new RevoluteJoint);
    rj->lolimit[0] = -M_PI_2;
    rj->hilimit[0] = M_PI_2;
    joint = rj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "continuous") == 0)
  {
    shared_ptr<RevoluteJoint> rj(new RevoluteJoint);
    rj->lolimit[0] = -10000.0;
    rj->hilimit[0] = 10000.0;
    joint = rj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "prismatic") == 0)
  {
    shared_ptr<PrismaticJoint> pj(new PrismaticJoint);
    pj->lolimit[0] = -10000.0;
    pj->hilimit[0] = 10000.0;
    joint = pj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "fixed") == 0)
    joint = shared_ptr<FixedJoint>(new FixedJoint);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "floating") == 0)
  {
    std::cerr << "URDFReader::read_joint() - [deprecated] floating joint type specified! not processing further..." << std::endl;
    return;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "planar") == 0)
  {
    std::cerr << "URDFReader::read_joint() - planar joint type currently unsupported in Moby! not processing further..." << std::endl;
    return;
  }
  else
  {
    std::cerr << "URDFReader::read_joint() - invalid joint type specified! not processing further..." << std::endl;
    return;
  }

  // read and verify required properties
  joint->id = name_attrib->get_string_value();
  if (!(inboard = read_parent(node, data, links)))
  {
    std::cerr << "URDFReader::read_joint() - failed to properly read parent link! not processing further..." << std::endl;
    return;
  }
  if (!(outboard = read_child(node, data, links)))
  {
    std::cerr << "URDFReader::read_joint() - failed to properly read child link! not processing further..." << std::endl;
    return;
  }

  // setup the appropriate pointers
  data.joint_parent[joint] = inboard;
  data.joint_child[joint] = outboard;

  // joint frame is defined relative to the parent link frame
  shared_ptr<Pose3d> origin(new Pose3d(read_origin(node, data)));
  origin->rpose = inboard->get_pose(); 
  Point3d location_origin(0.0, 0.0, 0.0, origin);
  Point3d location = Pose3d::transform_point(GLOBAL, location_origin);
  origin->update_relative_pose(outboard->get_pose()->rpose);

  // child frame is defined relative to the joint frame
  outboard->set_pose(*origin);

  // setup the inboard and outboard links for the joint
  joint->set_location(location, inboard, outboard);

  // read optional properties
  read_axis(node, data, joint);
  read_limits(node, data, joint);
  read_calibration(node, data, joint);
  read_dynamics(node, data, joint);
  read_safety_controller(node, data, joint);

  // add the joint to the set of joints 
  joints.push_back(joint);
}

/// Attempts to read the parent for the joint
RigidBodyPtr URDFReader::read_parent(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "parent") == 0)
    {
      // read the link attribute
      XMLAttrib* link_attrib = (*i)->get_attrib("link");
      if (!link_attrib)
        continue;
      string link_id = link_attrib->get_string_value();

      // find parent link
      for (unsigned j=0; j< links.size(); j++)
        if (links[j]->id == link_id)
          return links[j];
    }
  }

  return RigidBodyPtr();
}

/// Attempts to read the child for the joint
RigidBodyPtr URDFReader::read_child(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "child") == 0)
    {
      // read the link attribute
      XMLAttrib* link_attrib = (*i)->get_attrib("link");
      if (!link_attrib)
        continue;
      string link_id = link_attrib->get_string_value();

      // find child link
      for (unsigned j=0; j< links.size(); j++)
        if (links[j]->id == link_id)
          return links[j];
    }
  }

  return RigidBodyPtr();
}

/// Attempts to read axis for the joint
void URDFReader::read_axis(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
{
  Vector3d axis(1,0,0);
  bool axis_specified = false;

  // look for the axis tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "axis") == 0)
    {
      // read the attributes first
      XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (!xyz_attrib)
        continue;
      xyz_attrib->get_vector_value(axis);
      axis_specified = true;
    }
  }

  // get the outboard link - it's pose is identical to the joint pose
  RigidBodyPtr outboard = joint->get_outboard_link();

  // setup the axis frame
  axis.pose = joint->get_pose();

  // verify that joint is of the proper type
  shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
  shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
  if (rj || pj)
  {
    if (rj)
      rj->set_axis(axis);
    else if (pj)
      pj->set_axis(axis);
    else
      assert(false);
  }
  else if (axis_specified)
    std::cerr << "URDFReader::read_axis() - joint axis specified for joint w/o axis!" << std::endl;
}

/// Attempts to read dynamic properties for the joint
void URDFReader::read_dynamics(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
{
  // look for the dynamics tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "dynamics") == 0)
    {
      // read the attributes first
      XMLAttrib* damping_attrib = (*i)->get_attrib("damping");
      XMLAttrib* friction_attrib = (*i)->get_attrib("friction");

      // verify that joint is of the proper type
      shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
      shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
      if ((rj || pj) && (damping_attrib || friction_attrib))
      {
        // read the values
        joint->mu_fc = (friction_attrib) ? friction_attrib->get_real_value() : 0.0;
        joint->mu_fv = (damping_attrib) ? damping_attrib->get_real_value() : 0.0; 

        // multiple tags unsupported
        return;
      }

      // still here? unsupported type
      std::cerr << "URDFReader::read_dynamics() - unsupported joint type" << std::endl;
      return;  
    }
  }
}

/// Attempts to read limits for the joint
void URDFReader::read_limits(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
{
  // look for the limit tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "limit") == 0)
    {
      // read the attributes first
      XMLAttrib* lower_attrib = (*i)->get_attrib("lower");
      XMLAttrib* upper_attrib = (*i)->get_attrib("upper");
      XMLAttrib* effort_attrib = (*i)->get_attrib("effort");
      XMLAttrib* velocity_attrib = (*i)->get_attrib("velocity");

      // verify that joint is of the proper type
      shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
      shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
      if ((rj || pj) && (effort_attrib || lower_attrib || upper_attrib))
      {
        // if the velocity attribute was read, indicate that it is not used
        if (velocity_attrib)
          std::cerr << "URDFReader::read_limits() warning- 'velocity' attribute unsupported in Moby" << std::endl;

        // read the values
        if (lower_attrib)
          joint->lolimit[0] = lower_attrib->get_real_value();
        if (upper_attrib)
          joint->hilimit[0] = upper_attrib->get_real_value(); 
        if (effort_attrib)
          joint->maxforce[0] = effort_attrib->get_real_value();

        // multiple tags unsupported
        return;
      }
    }
  }
}

/// Attempts to read calibration for the joint
void URDFReader::read_calibration(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
{
  // look for the safety_controller tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "calibration") == 0)
    {
      std::cerr << "URDFReader::read_calibration() - calibration tag currently unused" << std::endl;
      return;
    }
  }
}

/// Attempts to reads the safety controller for the joint (currently unused)
void URDFReader::read_safety_controller(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
{
  // look for the safety_controller tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "safety_controller") == 0)
    {
      std::cerr << "URDFReader::read_safety_controller() - safety controller unsupported in Moby" << std::endl;
      return;
    }
  }
}

/// Attempts to read and set link inertial properties 
void URDFReader::read_inertial(shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link)
{
  // look for the inertial tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertial") == 0)
    {
      double mass = read_mass(*i, data);
      Matrix3d inertia = read_inertia(*i, data);

      // read the inertial frame
      shared_ptr<Pose3d> origin(new Pose3d(read_origin(*i, data)));

      // set the inertial frame relative to the link frame
      origin->rpose = link->get_pose();
      link->set_inertial_pose(*origin);

      // set inertial properties
      SpatialRBInertiad J(origin);
      J.m = mass;
      J.J = inertia;
      link->set_inertia(J);

      // reading inertial was a success, attempt to read no further...
      // (multiple inertial tags not supported)
      return;
    }
  }
}

/// Attempts to read an RGBA color
bool URDFReader::read_color(shared_ptr<const XMLTree> node, URDFData& data, VectorNd& color)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "color") == 0)
    {
      XMLAttrib* rgba_attrib = (*i)->get_attrib("rgba");
      if (rgba_attrib)
      {
        rgba_attrib->get_vector_value(color);
        if (color.size() == 4)
          return true;
      }
      else
        std::cerr << "URDFReader::read_color() warning - no rgba attribute found in color tag!" << std::endl;
    }
  }

  // no color read
  return false;
}

/// Attempts to read a texture filename
bool URDFReader::read_texture(shared_ptr<const XMLTree> node, URDFData& data, string& texture_fname)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "texture") == 0)
    {
      XMLAttrib* tfname_attrib = (*i)->get_attrib("filename");
      if (tfname_attrib)
      {
        texture_fname = tfname_attrib->get_string_value();
        return true;
      }
      else
        std::cerr << "URDFReader::read_texture() warning - no filename attribute found in color tag!" << std::endl;
    }
  }

  // no filename read
  return false;
}

/// Attempts to read an OSG Material
void URDFReader::read_material(shared_ptr<const XMLTree> node, URDFData& data, void* osg_node)
{
  #ifdef USE_OSG
  osg::Node* n = (osg::Node*) osg_node;

  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "material") == 0)
    {
      // read the name
      string material_name;
      XMLAttrib* name_attr = (*i)->get_attrib("name");
      if (name_attr)
        material_name = name_attr->get_string_value(); 

      // get the color and texture filename
      bool material_exists = (data.materials.find(material_name) != data.materials.end());
      pair<VectorNd, string>& material_pair = data.materials[material_name];
      VectorNd& color = material_pair.first;
      string& texture_fname = material_pair.second;

      // attempt to read color
      if (read_color(*i, data, color) || (material_exists && color.size() == 4))
      {
        const float R = (float) color[0];
        const float G = (float) color[1];
        const float B = (float) color[2];
        const float A = (float) color[3];
        osg::Material* material = new osg::Material;
        material->setColorMode(osg::Material::DIFFUSE);
        material->setDiffuse(osg::Material::FRONT, osg::Vec4(R, G, B, A));
        n->getOrCreateStateSet()->setAttribute(material);
      }

      // attempt to read texture
      if (read_texture(*i, data, texture_fname) || (material_exists && texture_fname.size() > 0))
      {
        osg::Texture2D* texture = new osg::Texture2D;
        texture->setDataVariance(osg::Object::DYNAMIC); 
        osg::Image* texture_image = osgDB::readImageFile(texture_fname.c_str());
        if (texture_image)
        {
          texture->setImage(texture_image);
          n->getOrCreateStateSet()->setTextureAttributeAndModes(0,texture,osg::StateAttribute::ON);
        }
        else
          std::cerr << "URDFReader::read_visual() - couldn't find texture filename " << texture_fname << std::endl;
      }
    }
  }      
  #endif
}

/// Attempts to read and set link visual properties
void URDFReader::read_visual(shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link)
{
  #ifdef USE_OSG
  VectorNd color;
  string texture_fname;

  // look for the visual tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "visual") == 0)
    {
      // read the primitive
      PrimitivePtr primitive = read_primitive(*i, data);
      if (!primitive)
        continue;

      // setup visualization 
      osg::MatrixTransform* group = new osg::MatrixTransform;
      group->ref();
      data.visual_transform_nodes[link] = group;
      link->set_visualization_data(group);
 
      // read the material
      read_material(*i, data, group);
      group->addChild(primitive->get_visualization());

      // setup the relative pose
      shared_ptr<Pose3d> origin(new Pose3d(read_origin(*i, data)));
      origin->rpose = link->get_pose();
      link->set_visualization_relative_pose(*origin);

      // reading visual was a success, attempt to read no further...
      // (multiple visual tags not supported)
      return;
    }
  }
  #endif
}

/// Attempts to read and set link collision properties
void URDFReader::read_collision(shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link)
{
  // look for the collision tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "collision") == 0)
    {
      // read the primitive
      PrimitivePtr primitive = read_primitive(*i, data);
      if (!primitive)
        continue;

      // construct collision geometry
      CollisionGeometryPtr cg(new CollisionGeometry);
      link->geometries.push_back(cg);
      cg->set_single_body(link);
      cg->set_geometry(primitive);

      // setup the relative pose
      shared_ptr<Pose3d> origin(new Pose3d(read_origin(*i, data)));
      origin->rpose = link->get_pose();
      cg->set_relative_pose(*origin);

      // reading collision geometry was a success, attempt to read no further...
      // (multiple collision tags not supported)
      return;
    }
  }
}

/// Reads Primitive from a geometry tag 
PrimitivePtr URDFReader::read_primitive(shared_ptr<const XMLTree> node, URDFData& data)
{
  PrimitivePtr primitive;

  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "geometry") == 0)
    {
      // read geometry 
      if ((primitive = read_box(*i, data)))
        return primitive;
//      else if ((primitive = read_cylinder(*i, data)))
//        return primitive;
      else if ((primitive = read_sphere(*i, data)))
        return primitive;
//      else if ((primitive = read_trimesh(*i, data)))
//        return primitive;
    }
  }

  return primitive;
}

/// Reads a box primitive
shared_ptr<BoxPrimitive> URDFReader::read_box(shared_ptr<const XMLTree> node, URDFData& data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "box") == 0)
    {
      XMLAttrib* size_attrib = (*i)->get_attrib("size");
      if (!size_attrib)
        continue;

      // get the size and construct the box
      Vector3d size;
      size_attrib->get_vector_value(size);
      return shared_ptr<BoxPrimitive>(new BoxPrimitive(size[X], size[Y], size[Z]));
    }
  }

  return shared_ptr<BoxPrimitive>();
}

/*
/// Reads a trimesh primitive
shared_ptr<TriangleMeshPrimitive> URDFReader::read_trimesh(shared_ptr<const XMLTree> node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "mesh") == 0)
    {
      XMLAttrib* filename_attrib = (*i)->get_attrib("filename");
      XMLAttrib* scale_attrib = (*i)->get_attrib("scale");
      if (!filename_attrib)
        continue;

      // warn that scale attribute is not used
      if (scale_attrib)
        std::cerr << "URDFReader::read_trimesh() warning- 'scale' attribute is not used" << std::endl;

      // construct the triangle mesh primitive 
      return shared_ptr<TriangleMeshPrimitive>(new TriangleMeshPrimitive(filename_attrib->get_string_value(),false));
    }
  }

  return shared_ptr<TriangleMeshPrimitive>();
}

/// Reads a cylinder primitive
shared_ptr<CylinderPrimitive> URDFReader::read_cylinder(shared_ptr<const XMLTree> node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "cylinder") == 0)
    {
      XMLAttrib* radius_attrib = (*i)->get_attrib("radius");
      XMLAttrib* length_attrib = (*i)->get_attrib("length");
      if (!radius_attrib || !length_attrib)
        continue;

      // cylinder must be rotated around x axis
      Matrix3d rx90 = Matrix3d::rot_X(M_PI_2);
      Pose3d T(rx90, Origin3d::zero());

      // construct the cylinder 
      double radius = radius_attrib->get_real_value();
      double length = length_attrib->get_real_value();
      return shared_ptr<CylinderPrimitive>(new CylinderPrimitive(radius, length, T));
    }
  }

  return shared_ptr<CylinderPrimitive>();
}
*/

/// Reads a sphere primitive
shared_ptr<SpherePrimitive> URDFReader::read_sphere(shared_ptr<const XMLTree> node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "sphere") == 0)
    {
      XMLAttrib* radius_attrib = (*i)->get_attrib("radius");
      if (!radius_attrib)
        continue;

      // get the size and construct the sphere
      double radius = radius_attrib->get_real_value();
      return shared_ptr<SpherePrimitive>(new SpherePrimitive(radius));
    }
  }

  return shared_ptr<SpherePrimitive>();
}


/// Attempts to read an "origin" tag
Pose3d URDFReader::read_origin(shared_ptr<const XMLTree> node, URDFData& data)
{
  Origin3d xyz;
  Vector3d rpy;

  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "origin") == 0)
    {
      // look for xyz attribute 
      XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (xyz_attrib)
        xyz = xyz_attrib->get_origin_value();

      // look for rpy attribute
      XMLAttrib* rpy_attrib = (*i)->get_attrib("rpy");
      if (rpy_attrib)
        rpy_attrib->get_vector_value(rpy);

      // reading tag was a success, attempt to read no further...
      // (multiple such tags not supported)
      break;
    }
  }

  Quatd rpy_quat = Quatd::rpy(rpy[0], rpy[1], rpy[2]);
  return Pose3d(rpy_quat, xyz);
}

/// Attempts to read a "mass" tag
double URDFReader::read_mass(shared_ptr<const XMLTree> node, URDFData& data)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "mass") == 0)
    {
      // look for the "value" attribute
      XMLAttrib* value_attrib = (*i)->get_attrib("value");
      if (value_attrib)
        // reading tag was a success, attempt to read no further...
        // (multiple such tags not supported)
        return value_attrib->get_real_value(); 
    }
  }

  // couldn't find the tag.. return 0
  return 0.0;
}

/// Attempts to read an "inertia" tag
Matrix3d URDFReader::read_inertia(shared_ptr<const XMLTree> node, URDFData& data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup J to zero initially
  Matrix3d J = Matrix3d::zero();

  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertia") == 0)
    {
      // look for the six attributes
      XMLAttrib* ixx_attrib = (*i)->get_attrib("ixx");
      XMLAttrib* ixy_attrib = (*i)->get_attrib("ixy");
      XMLAttrib* ixz_attrib = (*i)->get_attrib("ixz");
      XMLAttrib* iyy_attrib = (*i)->get_attrib("iyy");
      XMLAttrib* iyz_attrib = (*i)->get_attrib("iyz");
      XMLAttrib* izz_attrib = (*i)->get_attrib("izz");

      // set values from present attributes
      if (ixx_attrib)
        J(X,X) = ixx_attrib->get_real_value();
      if (iyy_attrib)
        J(Y,Y) = iyy_attrib->get_real_value();
      if (izz_attrib)
        J(Z,Z) = izz_attrib->get_real_value();
      if (ixy_attrib)
        J(X,Y) = J(Y,X) = ixy_attrib->get_real_value();
      if (ixz_attrib)
        J(X,Z) = J(Z,X) = ixz_attrib->get_real_value();
      if (iyz_attrib)
        J(Y,Z) = J(Z,Y) = iyz_attrib->get_real_value();

      // reading tag was a success, attempt to read no further...
      // (multiple such tags not supported)
      return J;
    }
  }

  // no inertia read.. return default J (0 matrix)
  return J;
}



