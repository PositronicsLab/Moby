/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <dlfcn.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <stack>

#ifdef USE_OSG
#include <osg/Material>
#include <osg/Group>
#include <osg/Texture2D>
#include <osg/MatrixTransform>
#include <osgDB/ReadFile>
#endif

#include <Moby/CSG.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/FixedJoint.h>
#include <Moby/MCArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/UniversalJoint.h>
#include <Moby/AAngle.h>
#include <Moby/DeformableCCD.h>
#include <Moby/XMLTree.h>
#include <Moby/URDFReader.h>

using namespace Moby;
using std::set;
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
map<string, pair<vector<RigidBodyPtr>, vector<JointPtr> > > URDFReader::read(const string& fname)
{
  // setup the list of robots 
  map<string, pair<vector<RigidBodyPtr>, vector<JointPtr> > > robot_map;
  
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
      return robot_map;
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
  XMLTreeConstPtr tree = XMLTree::read_from_xml(filename);
  if (!tree)
  {
    std::cerr << "URDFReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return robot_map;
  }
  
  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RigidBody is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // read and construct all robots 
  const list<XMLTreePtr>& child_nodes = tree->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    if (strcasecmp((*i)->name.c_str(), "Robot") == 0)
    {
      URDFData data;
      read_robot(*i, data, robot_map);
    }

  // change back to the initial working directory
  chdir(cwd.get());

  return robot_map;
}

/// Reads and constructs a robot object
void URDFReader::read_robot(XMLTreeConstPtr node, URDFData& data, map<string, pair<vector<RigidBodyPtr>, vector<JointPtr> > >& robot_map)
{
  // sanity check
  assert(strcasecmp(node->name.c_str(), "Robot") == 0);

  // read the robot name
  const XMLAttrib* name_attrib = node->get_attrib("name");
  if (name_attrib)
  {
    string name = name_attrib->get_string_value();
    pair<vector<RigidBodyPtr>, vector<JointPtr> >& robot = robot_map[name];
    read_links(node, data, robot.first);
    read_joints(node, data, robot.first, robot.second);
    if (!transform_frames(data, robot.first, robot.second))
      robot_map.erase(robot_map.find(name));
  }
  else
    std::cerr << "URDFReader::read_robot() - robot name not specified! not processing further..." << std::endl;
}

/// Reads robot links
void URDFReader::read_links(XMLTreeConstPtr node, URDFData& data, vector<RigidBodyPtr>& links)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_link(*i, data, links);
}

/// Attempts to read a robot link from the given node
void URDFReader::read_link(XMLTreeConstPtr node, URDFData& data, vector<RigidBodyPtr>& links)
{
  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Link") != 0)
    return;

  // link must have the name attribute
  const XMLAttrib* name_attrib = node->get_attrib("name");
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
void URDFReader::find_children(const URDFData& data, const vector<RigidBodyPtr>& links, RigidBodyPtr link, queue<RigidBodyPtr>& q, map<RigidBodyPtr, RigidBodyPtr>& parents)
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
static void to_osg_matrix(const Matrix4& src, osg::Matrixd& tgt)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= W; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(i,j) = src(j,i);

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (Real) 0.0;
  tgt(W,W) = (Real) 1.0;
}
#endif

/// Transforms frames to Moby format (link positions at c.o.m.)
/**
 *  URDF defines the reference frame for a link to be identical to the inner
 * joint frame. Inertial, collision, and visualization frames are defined
 * relative to this frame. Joint axes are defined relative to this frame
 * as well.
 */
bool URDFReader::transform_frames(const URDFData& data, const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  map<RigidBodyPtr, RigidBodyPtr> parents; 
  Matrix3 J_out;
  Vector3 com_out;

  // find the base (it will be the link that does not exist as a child)
  set<RigidBodyPtr> candidates;
  for (unsigned i=0; i< links.size(); i++)
    candidates.insert(links[i]);
  for (std::map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    candidates.erase(i->second);

  // must be *exactly* one candidate remaining
  if (candidates.size() != 1)
  {
    std::cerr << "URDFReader::transform_frames() - no bases or multiple bases specified! not reading further!" << std::endl;
    return false;
  }
  RigidBodyPtr base = *candidates.begin();

  // need to transform base inertia orientation to be relative to global frame
  Matrix4 baseT = data.inertia_transforms.find(base)->second; 
  Vector3 com_in = baseT.get_translation();
  baseT.set_translation(ZEROS_3);
  Primitive::transform_inertia(base->get_mass(), base->get_inertia(), com_in, baseT, J_out, com_out); 
  base->set_inertia(J_out);

  // update base position and orientation
  base->set_position(com_out);
  base->set_orientation(Quat(&IDENTITY_3x3));

  // add all children to the queue  
  queue<RigidBodyPtr> q;
  find_children(data, links, base, q, parents);
  while (!q.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = q.front();
    q.pop();

    // get the parent link
    assert(parents.find(link) != parents.end());
    RigidBodyPtr parent = parents.find(link)->second;

    // get the joint transform
    JointPtr joint = find_joint(data, link);
    assert(data.joint_transforms.find(joint) != data.joint_transforms.end());
    const Matrix4& Tj = data.joint_transforms.find(joint)->second;
    const Matrix3 Rj = Tj.get_rotation();

    // compute the *reference* frame (w.r.t. the global frame) for the child 
    // (and the joint)
    Matrix4 Tref = parent->get_transform() * Tj;

    // get the inertial reference frame for the child
    assert(data.inertia_transforms.find(link) != data.inertia_transforms.end());
    Matrix4 iF = data.inertia_transforms.find(link)->second;
    com_in = Tref.mult_point(iF.get_translation());
    iF.set_translation(ZEROS_3);

    // inertia for child is relative to reference frame; make it relative to
    // frame aligned with global and located at child c.o.m.
    Primitive::transform_inertia(link->get_mass(), link->get_inertia(), com_in, iF, J_out, com_out);

    // setup the link position and orientation
    link->set_position(com_out);
    link->set_orientation(Quat(&Rj));     

    // compute the relative collision transform
    if (data.collision_transforms.find(link) != data.collision_transforms.end())
    {
      // get the true frame of the collision transform
      const Matrix4& Trel = data.collision_transforms.find(link)->second;
      Matrix4 Tc = Tref * Trel;

      // now compute the collision transform relative to the link frame
      Matrix4 Tc_rel = link->get_transform().inverse_transform() * Tc;
      assert(!link->geometries.empty());
      link->geometries.front()->set_rel_transform(Tc_rel, false); 
    }

    // compute the relative visual transform
    if (data.visual_transforms.find(link) != data.visual_transforms.end())
    {
      // get the true frame of the visual transform
      const Matrix4& Trel = data.visual_transforms.find(link)->second;
      Matrix4 Tv = Tref * Trel;

      // now compute the collision transform relative to the link frame
      Matrix4 Tv_rel = link->get_transform().inverse_transform() * Tv;
      if (data.visual_transform_nodes.find(link) != data.visual_transform_nodes.end())
      {
        #ifdef USE_OSG
        osg::MatrixTransform* group = (osg::MatrixTransform*) data.visual_transform_nodes.find(link)->second;
        osg::Matrixd m;
        to_osg_matrix(Tv_rel, m);
        group->setMatrix(m);
        #endif
      }
    }

    // get the joint position and the two vectors we need (all in global frame) 
    Vector3 joint_position = Tref.get_translation();
    Vector3 joint_to_child_com = link->get_position() - joint_position;
    Vector3 parent_com_to_joint = joint_position - parent->get_position();

    // setup the pointers between the joints and links
    link->add_inner_joint(parent, joint, Tref.inverse_mult_point(joint_to_child_com), link->get_transform().inverse_mult_point(joint_to_child_com));
    parent->add_outer_joint(link, joint, parent->get_transform().inverse_mult_point(parent_com_to_joint));

    // setup the joint axis, if necessary
    if (data.joint_axes.find(joint) != data.joint_axes.end())
    {
      Vector3 axis = data.joint_axes.find(joint)->second;
      shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
      shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
      if (pj)
        pj->set_axis_global(axis);
      else if (rj)
        rj->set_axis_global(axis);
    }
    
    // add all children of link to the queue
    find_children(data, links, link, q, parents);
  }

  return true;
}

/// Attempts to read a robot joint from the given node
void URDFReader::read_joint(XMLTreeConstPtr node, URDFData& data, const vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  JointPtr joint;
  RigidBodyPtr inboard, outboard;

  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Joint") != 0)
    return;

  // link must have the name attribute
  const XMLAttrib* name_attrib = node->get_attrib("name");
  if (!name_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint name not specified! not processing further..." << std::endl;
    return;
  }

  // link must have the name attribute
  const XMLAttrib* type_attrib = node->get_attrib("type");
  if (!type_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint type not specified! not processing further..." << std::endl;
    return;
  }

  // read and construct the joint
  if (strcasecmp(type_attrib->get_string_value().c_str(), "revolute"))
    joint = shared_ptr<RevoluteJoint>(new RevoluteJoint);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "continuous"))
    joint = shared_ptr<RevoluteJoint>(new RevoluteJoint);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "prismatic"))
    joint = shared_ptr<PrismaticJoint>(new PrismaticJoint);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "fixed"))
    joint = shared_ptr<FixedJoint>(new FixedJoint);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "floating"))
  {
    std::cerr << "URDFReader::read_joint() - [deprecated] floating joint type specified! not processing further..." << std::endl;
    return;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "planar"))
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

  // read optional properties
  data.joint_transforms[joint] = read_origin(node, data);
  read_axis(node, data, joint);
  read_limits(node, data, joint);
  read_calibration(node, data, joint);
  read_dynamics(node, data, joint);
  read_safety_controller(node, data, joint);

  // setup the appropriate pointers
  data.joint_parent[joint] = inboard;
  data.joint_child[joint] = outboard;

  // add the joint to the set of joints 
  joints.push_back(joint);
}

/// Attempts to read the parent for the joint
RigidBodyPtr URDFReader::read_parent(XMLTreeConstPtr node, URDFData& data, const vector<RigidBodyPtr>& links)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "parent") == 0)
    {
      // read the link attribute
      const XMLAttrib* link_attrib = (*i)->get_attrib("link");
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
RigidBodyPtr URDFReader::read_child(XMLTreeConstPtr node, URDFData& data, const vector<RigidBodyPtr>& links)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "child") == 0)
    {
      // read the link attribute
      const XMLAttrib* link_attrib = (*i)->get_attrib("link");
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
void URDFReader::read_axis(XMLTreeConstPtr node, URDFData& data, JointPtr joint)
{
  Vector3 axis(1,0,0);
  bool axis_specified = false;

  // look for the axis tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "axis") == 0)
    {
      // read the attributes first
      const XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (!xyz_attrib)
        continue;
      xyz_attrib->get_vector_value(axis);
      axis_specified = true;
    }
  }

  // verify that joint is of the proper type
  shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
  shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
  if (rj || pj)
    data.joint_axes[joint] = axis;
  else if (axis_specified)
    std::cerr << "URDFReader::read_axis() - joint axis specified for joint w/o axis!" << std::endl;
}

/// Attempts to read dynamic properties for the joint
void URDFReader::read_dynamics(XMLTreeConstPtr node, URDFData& data, JointPtr joint)
{
  // look for the dynamics tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "dynamics") == 0)
    {
      // read the attributes first
      const XMLAttrib* damping_attrib = (*i)->get_attrib("damping");
      const XMLAttrib* friction_attrib = (*i)->get_attrib("friction");

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
      std::cerr << "URDFReader::read_limits() - unsupported joint type for limits" << std::endl;
      return;  
    }
  }
}

/// Attempts to read limits for the joint
void URDFReader::read_limits(XMLTreeConstPtr node, URDFData& data, JointPtr joint)
{
  // look for the limit tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "limit") == 0)
    {
      // read the attributes first
      const XMLAttrib* lower_attrib = (*i)->get_attrib("lower");
      const XMLAttrib* upper_attrib = (*i)->get_attrib("upper");
      const XMLAttrib* effort_attrib = (*i)->get_attrib("effort");
      const XMLAttrib* velocity_attrib = (*i)->get_attrib("velocity");

      // verify that joint is of the proper type
      shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
      shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
      if ((rj || pj) && (effort_attrib || lower_attrib || upper_attrib))
      {
        // if the velocity attribute was read, indicate that it is not used
        if (velocity_attrib)
          std::cerr << "URDFReader::read_limits() warning- 'velocity' attribute unsupported in Moby" << std::endl;

        // read the values
        joint->lolimit[0] = (lower_attrib) ? lower_attrib->get_real_value() : 0.0;
        joint->hilimit[0] = (upper_attrib) ? upper_attrib->get_real_value() : 0.0; 
        joint->maxforce[0] = (effort_attrib) ? effort_attrib->get_real_value() : 0.0;

        // multiple tags unsupported
        return;
      }
    }
  }
}

/// Attempts to read calibration for the joint
void URDFReader::read_calibration(XMLTreeConstPtr node, URDFData& data, JointPtr joint)
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
void URDFReader::read_safety_controller(XMLTreeConstPtr node, URDFData& data, JointPtr joint)
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
void URDFReader::read_inertial(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link)
{
  // look for the inertial tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertial") == 0)
    {
      Real mass = read_mass(*i, data);
      Matrix3 inertia = read_inertia(*i, data);

      // set inertial properties
      link->set_mass(mass);
      link->set_inertia(inertia);

      // store relative transform
      data.inertia_transforms[link] = read_origin(*i, data);

      // reading inertial was a success, attempt to read no further...
      // (multiple inertial tags not supported)
      return;
    }
  }
}

/// Attempts to read an RGBA color
bool URDFReader::read_color(XMLTreeConstPtr node, URDFData& data, VectorN& color)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "color") == 0)
    {
      const XMLAttrib* rgba_attrib = (*i)->get_attrib("rgba");
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
bool URDFReader::read_texture(XMLTreeConstPtr node, URDFData& data, string& texture_fname)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "texture") == 0)
    {
      const XMLAttrib* tfname_attrib = (*i)->get_attrib("filename");
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
void URDFReader::read_material(XMLTreeConstPtr node, URDFData& data, void* osg_node)
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
      const XMLAttrib* name_attr = (*i)->get_attrib("name");
      if (name_attr)
        material_name = name_attr->get_string_value(); 

      // get the color and texture filename
      bool material_exists = (data.materials.find(material_name) != data.materials.end());
      pair<VectorN, string>& material_pair = data.materials[material_name];
      VectorN& color = material_pair.first;
      string& texture_fname = material_pair.second;

      // attempt to read color
      if (read_color(*i, data, color) || (material_exists && color.size() == 4))
      {
        const Real R = color[0];
        const Real G = color[1];
        const Real B = color[2];
        const Real A = color[3];
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
void URDFReader::read_visual(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link)
{
  #ifdef USG_OSG
  VectorN color;
  string texture_fname;

  // look for the visual tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "visual") == 0)
    {
      // read the origin
      data.visual_transforms[link] = read_origin(*i, data);

      // read the primitive
      PrimitivePtr primitive = read_primitive(*i);
      if (!primitive)
        continue;

      // setup visualization 
      osg::MatrixTransform* group = new osg::MatrixTransform;
      visual_transform_nodes[link] = group;
      group->addChild(primitive->get_visualization());
      link->set_visualization_data(group);
 
      // read the material
      read_material(node, group);

      // reading visual was a success, attempt to read no further...
      // (multiple visual tags not supported)
      return;
    }
  }
  #endif
}

/// Attempts to read and set link collision properties
void URDFReader::read_collision(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link)
{
  // look for the collision tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "collision") == 0)
    {
      // read the origin
      data.collision_transforms[link] = read_origin(*i, data);

      // read the primitive
      PrimitivePtr primitive = read_primitive(*i, data);
      if (!primitive)
        continue;

      // construct collision geometry
      CollisionGeometryPtr cg(new CollisionGeometry);
      cg->set_geometry(primitive);
      link->geometries.push_back(cg);

      // reading collision geometry was a success, attempt to read no further...
      // (multiple collision tags not supported)
      return;
    }
  }
}

/// Reads Primitive from a geometry tag 
PrimitivePtr URDFReader::read_primitive(XMLTreeConstPtr node, URDFData& data)
{
  PrimitivePtr primitive;

  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "geometry") == 0)
    {
      // read geometry 
      if (primitive = read_box(*i, data))
        return primitive;
      else if (primitive = read_cylinder(*i, data))
        return primitive;
      else if (primitive = read_sphere(*i, data))
        return primitive;
      else if (primitive = read_trimesh(*i, data))
        return primitive;
    }
  }

  return primitive;
}

/// Reads a box primitive
shared_ptr<BoxPrimitive> URDFReader::read_box(XMLTreeConstPtr node, URDFData& data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "box") == 0)
    {
      const XMLAttrib* size_attrib = (*i)->get_attrib("size");
      if (!size_attrib)
        continue;

      // get the size and construct the box
      Vector3 size;
      size_attrib->get_vector_value(size);
      return shared_ptr<BoxPrimitive>(new BoxPrimitive(size[X], size[Y], size[Z]));
    }
  }

  return shared_ptr<BoxPrimitive>();
}

/// Reads a trimesh primitive
shared_ptr<TriangleMeshPrimitive> URDFReader::read_trimesh(XMLTreeConstPtr node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "mesh") == 0)
    {
      const XMLAttrib* filename_attrib = (*i)->get_attrib("filename");
      const XMLAttrib* scale_attrib = (*i)->get_attrib("scale");
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
shared_ptr<CylinderPrimitive> URDFReader::read_cylinder(XMLTreeConstPtr node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "cylinder") == 0)
    {
      const XMLAttrib* radius_attrib = (*i)->get_attrib("radius");
      const XMLAttrib* length_attrib = (*i)->get_attrib("length");
      if (!radius_attrib || !length_attrib)
        continue;

      // cylinder must be rotated around x axis
      Matrix3 rx90 = Matrix3::rot_X(M_PI_2);
      const Matrix4& T = Matrix4(&rx90, &ZEROS_3);

      // construct the cylinder 
      Real radius = radius_attrib->get_real_value();
      Real length = length_attrib->get_real_value();
      return shared_ptr<CylinderPrimitive>(new CylinderPrimitive(radius, length, T));
    }
  }

  return shared_ptr<CylinderPrimitive>();
}

/// Reads a sphere primitive
shared_ptr<SpherePrimitive> URDFReader::read_sphere(XMLTreeConstPtr node, URDFData& data)
{
  // look for the geometry tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "sphere") == 0)
    {
      const XMLAttrib* radius_attrib = (*i)->get_attrib("radius");
      if (!radius_attrib)
        continue;

      // get the size and construct the sphere
      Real radius = radius_attrib->get_real_value();
      return shared_ptr<SpherePrimitive>(new SpherePrimitive(radius));
    }
  }

  return shared_ptr<SpherePrimitive>();
}


/// Attempts to read an "origin" tag
Matrix4 URDFReader::read_origin(XMLTreeConstPtr node, URDFData& data)
{
  Vector3 xyz = ZEROS_3;
  Vector3 rpy = ZEROS_3;

  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "origin") == 0)
    {
      // look for xyz attribute 
      const XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (xyz_attrib)
        xyz_attrib->get_vector_value(xyz);

      // look for rpy attribute
      const XMLAttrib* rpy_attrib = (*i)->get_attrib("rpy");
      if (rpy_attrib)
        rpy_attrib->get_vector_value(rpy);

      // reading tag was a success, attempt to read no further...
      // (multiple such tags not supported)
      break;
    }
  }

  Matrix3 rpy_matrix = Matrix3::rpy(rpy[0], rpy[1], rpy[2]); 
  return Matrix4(&rpy_matrix, &xyz);
}

/// Attempts to read a "mass" tag
Real URDFReader::read_mass(XMLTreeConstPtr node, URDFData& data)
{
  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "mass") == 0)
    {
      // look for the "value" attribute
      const XMLAttrib* value_attrib = (*i)->get_attrib("value");
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
Matrix3 URDFReader::read_inertia(XMLTreeConstPtr node, URDFData& data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup J to zero initially
  Matrix3 J = ZEROS_3x3;

  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertia") == 0)
    {
      // look for the six attributes
      const XMLAttrib* ixx_attrib = (*i)->get_attrib("ixx");
      const XMLAttrib* ixy_attrib = (*i)->get_attrib("ixy");
      const XMLAttrib* ixz_attrib = (*i)->get_attrib("ixz");
      const XMLAttrib* iyy_attrib = (*i)->get_attrib("iyy");
      const XMLAttrib* iyz_attrib = (*i)->get_attrib("iyz");
      const XMLAttrib* izz_attrib = (*i)->get_attrib("izz");

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



