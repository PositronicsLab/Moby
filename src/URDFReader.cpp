/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
  XMLTreeConstPtr tree = XMLTree::read_from_xml(filename);
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
bool URDFReader::read_robot(XMLTreeConstPtr node, URDFData& data, string& name, vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  // read the robot name
  const XMLAttrib* name_attrib = node->get_attrib("name");
  if (name_attrib)
  {
    name = name_attrib->get_string_value();
    read_links(node, data, links);
    read_joints(node, data, links, joints);
    if (!transform_frames(data, links, joints))
      return false; 
  }
  else
  {
    std::cerr << "URDFReader::read_robot() - robot name not specified! not processing further..." << std::endl;
    return false;
  }

  return true;
}

/// Reads robot links
void URDFReader::read_links(XMLTreeConstPtr node, URDFData& data, vector<RigidBodyPtr>& links)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_link(*i, data, links);
}

/// Reads robot joints 
void URDFReader::read_joints(XMLTreeConstPtr node, URDFData& data, const vector<RigidBodyPtr>& links, vector<JointPtr>& joints)
{
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_joint(*i, data, links, joints);
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

/// Fix for Moby
void URDFReader::fix_Moby(URDFData& data, const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  // find the base (it will be the link that does not exist as a child)
  set<RigidBodyPtr> candidates;
  for (unsigned i=0; i< links.size(); i++)
    candidates.insert(links[i]);
  for (std::map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    candidates.erase(i->second);

  // must be *exactly* one candidate remaining
  if (candidates.size() != 1)
  {
    std::cerr << "URDFReader::fix_Moby() - no bases or multiple bases specified! not reading further!" << std::endl;
    return;
  }
  RigidBodyPtr base = *candidates.begin();

  // determine leaf links 
  set<RigidBodyPtr> leafs;
  for (unsigned i=0; i< links.size(); i++)
    leafs.insert(links[i]);
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
    leafs.erase(i->second);

  // add all leafs to a queue
  queue<RigidBodyPtr> q;
  BOOST_FOREACH(RigidBodyPtr link, leafs)
    q.push(link);

  // setup indicator for when a link has been processed
  set<RigidBodyPtr> processed;

  // process from leafs inward 
  while (!q.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = q.front();
    q.pop();

    // see whether the link has already been processed
    if (processed.find(link) != processed.end())
      continue;
    else
      processed.insert(link);

    // get the parent link and the joint
    JointPtr joint = find_joint(data, link);
    RigidBodyPtr parent = data.joint_parent.find(joint)->second;

    // see whether this relative transform is ok
    if (!valid_transformation(data, parent, joint, link))
    {
      // get the joint frame / outboard frame
      Matrix4 Tj = data.joint_transforms.find(joint)->second;

      // propagate transform
      propagate_transform(data, link, Tj.get_rotation());
    }

    // add parent to the queue for processing, unless the parent is the base 
    if (parent != base)
      q.push(parent);
  }
}

// Propagates a rotation through a link and all children
void URDFReader::propagate_transform(URDFData& data, RigidBodyPtr link, const Matrix3& Rx)
{
  // add the link to the queue
  std::queue<RigidBodyPtr> q;
  q.push(link);

  // process until queue is empty
  while (!q.empty())
  {
    // get the link off of the front of the queue
    link = q.front();
    q.pop();

    // get the inner joint
    JointPtr joint = find_joint(data, link);

    // determine Tx
    Matrix4 Tj = data.joint_transforms[joint];
    Vector3 xlat_des = Tj.get_translation();
    Matrix3 Rdes = Rx * Tj.get_rotation();
    Matrix4 Tdes(&IDENTITY_3x3, &xlat_des);
    Matrix4 Tx = Tdes * Matrix4::inverse_transform(Tj);

    // update inertial frame, visual frame, collision frame
    data.joint_transforms[joint] = Tx * data.joint_transforms[joint];
    data.inertia_transforms[link] = Tx * data.inertia_transforms[link];
    data.visual_transforms[link] = Tx * data.visual_transforms[link];
    data.collision_transforms[link] = Tx * data.collision_transforms[link];

    // update joint axis, if it is present
    if (data.joint_axes.find(joint) != data.joint_axes.end())
      data.joint_axes[joint] = Tx.mult_vector(data.joint_axes[joint]);

break;
    // add all children of the link to the queue
    for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
      if (i->second == link)
      {
        assert(data.joint_child.find(i->first) != data.joint_child.end());
        RigidBodyPtr child = data.joint_child.find(i->first)->second;
        q.push(child);
      } 
  }
}

/// Determine whether it is possible to get from one reference frame to another
bool URDFReader::valid_transformation(const URDFData& data, RigidBodyPtr inboard, JointPtr joint, RigidBodyPtr outboard)
{
  // get the inboard reference frame (if it's available- it won't be if the
  // inboard is the base)
  Matrix4 Tref_inboard = IDENTITY_4x4;
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    if (i->second == inboard)
    {
      Tref_inboard = data.joint_transforms.find(i->first)->second;
      break;
    } 

  // get the outboard reference frame
  assert(data.joint_transforms.find(joint) != data.joint_transforms.end());
  Matrix4 Tref_outboard = data.joint_transforms.find(joint)->second;

  // get the relative transform and relative rotation
  Matrix4 Trel = Matrix4::inverse_transform(Tref_inboard) * Tref_outboard;
  Matrix3 Rrel = Trel.get_rotation();

  // see whether we can get from Tref_inboard to Tref_outboard based on
  // joint movement alone
  if (dynamic_pointer_cast<RevoluteJoint>(joint))
  {
    // get the joint axis relative to inboard frame
    Vector3 axis_0 = Tref_outboard.mult_vector(data.joint_axes.find(joint)->second);
    Vector3 axis = Tref_inboard.transpose_mult_vector(axis_0);

    // verify that rotation around that axis gives the desired rotation
    AAngle aa(&Rrel, &axis);
    Matrix3 Raa(&aa);
    return std::fabs(Quat::calc_angle(Quat(&Rrel), Quat(&Raa))) < NEAR_ZERO;    
  }
  else if (dynamic_pointer_cast<PrismaticJoint>(joint) || dynamic_pointer_cast<FixedJoint>(joint))
  {
    // both orientations must be the same
    Quat q(&Rrel);
    Quat eye = Quat::identity();
    return std::fabs(Quat::calc_angle(q, eye)) < NEAR_ZERO;
  }
  else
    // should never be here!
    assert(false);

  // req'd by compiler...
  return false; 
}

void URDFReader::output_data(const URDFData& data, RigidBodyPtr link)
{
  std::cout << "link id: " << link->id << std::endl;
  JointPtr joint = find_joint(data, link);
  if (joint)
  {
    std::cout << "  modified URDF joint transform: " << std::endl << data.joint_transforms.find(joint)->second << std::endl;
    if (data.joint_axes.find(joint) != data.joint_axes.end())
      std::cout << "  modified URDF joint axis: " << data.joint_axes.find(joint)->second << std::endl;
    shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
    shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
    if (rj)
      std::cout << "  Moby joint axis: " << rj->get_axis_global() << std::endl;
    else if (pj)
      std::cout << "  Moby joint axis: " << pj->get_axis_global() << std::endl;
  }
  if (data.visual_transforms.find(link) != data.visual_transforms.end())
    std::cout << "  modified URDF visualization transform: " << std::endl << data.visual_transforms.find(link)->second;
  if (data.inertia_transforms.find(link) != data.inertia_transforms.end())
    std::cout << "  modified URDF inertia transform: " << std::endl << data.inertia_transforms.find(link)->second;
  if (data.collision_transforms.find(link) != data.collision_transforms.end())
    std::cout << "  modified URDF collision transform: " << std::endl << data.collision_transforms.find(link)->second;
  std::cout << "  Moby transform: " << std::endl << link->get_transform();
  if (!link->geometries.empty())
    std::cout << "  Moby relative collision transform: " << std::endl << link->geometries.front()->get_rel_transform();
}

/// Transforms frames to Moby format (link positions at c.o.m.)
/**
 *  URDF defines the reference frame for a link to be identical to the inner
 * joint frame. Inertial, collision, and visualization frames are defined
 * relative to this frame. Joint axes are defined relative to this frame
 * as well.
 */
bool URDFReader::transform_frames(URDFData& data, const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  map<RigidBodyPtr, RigidBodyPtr> parents; 
  Matrix3 J_out;
  Vector3 com_out;

  // output all data (before Moby fixes it)
  std::cout << "data (pre-fixes)" << std::endl;
  for (unsigned i=0; i< links.size(); i++)
    output_data(data, links[i]);
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "data (post-fixes)" << std::endl;

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

  // fix frame data for Moby 
  fix_Moby(data, links, joints);

  // need to transform base inertia orientation to be relative to global frame
  Matrix4 baseT = data.inertia_transforms.find(base)->second; 
  Vector3 com_in = baseT.get_translation();
  baseT.set_translation(ZEROS_3);
  Primitive::transform_inertia(base->get_mass(), base->get_inertia(), com_in, baseT, J_out, com_out); 
  base->set_inertia(J_out);

  // update base position and orientation
  base->set_position(com_out);
  base->set_orientation(Quat(&IDENTITY_3x3));

  // compute the relative collision transform
  if (data.collision_transforms.find(base) != data.collision_transforms.end())
  {
    // get the true frame of the collision transform
    const Matrix4& Trel = data.collision_transforms.find(base)->second;
    Matrix4 Tc = baseT * Trel;

    // now compute the collision transform relative to the base frame
    Matrix4 Tc_rel = base->get_transform().inverse_transform() * Tc;
    assert(!base->geometries.empty());
    base->geometries.front()->set_rel_transform(Tc_rel, false); 
  }

  // compute the relative visual transform
  if (data.visual_transforms.find(base) != data.visual_transforms.end())
  {
    // get the true frame of the visual transform
    const Matrix4& Trel = data.visual_transforms.find(base)->second;
    Matrix4 Tv = baseT * Trel;

    // now compute the collision transform relative to the base frame
    Matrix4 Tv_rel = base->get_transform().inverse_transform() * Tv;
    #ifdef USE_OSG
    osg::MatrixTransform* group = (osg::MatrixTransform*) data.visual_transform_nodes.find(base)->second;
    osg::Matrixd m;
    to_osg_matrix(Tv_rel, m);
    group->setMatrix(m);
    #endif
  }

  output_data(data, base);

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

    // get the joint transform *relative to parent link reference frame*
    JointPtr joint = find_joint(data, link);
    assert(data.joint_transforms.find(joint) != data.joint_transforms.end());
    const Matrix4& Tj = data.joint_transforms.find(joint)->second;

    // compute the reference frame for the parent
    Vector3 parent_ref_xlat = data.inertia_transforms.find(parent)->second.get_translation();
    Matrix4 Tref_parent = parent->get_transform();
    Vector3 com_to_ref = Tref_parent.mult_vector(-parent_ref_xlat);
    Tref_parent.set_translation(Tref_parent.get_translation() + com_to_ref); 

    // compute the *reference* frame (w.r.t. the global frame) for the child 
    // (and the joint)
    Matrix4 Tref = Tref_parent * Tj;
    std::cout << "link " << link->id << " reference frame: " << std::endl << Tref;

    // compute the desired global orientation for the child
    Matrix3 Rdes = Tref.get_rotation();

    // determine the transform from the reference frame to Moby frame
    assert(data.inertia_transforms.find(link) != data.inertia_transforms.end());
    Matrix4 iF = data.inertia_transforms.find(link)->second;
    Matrix3 iR = iF.get_rotation();
    Matrix4 mTr = IDENTITY_4x4;
    mTr.set_translation(iF.get_translation());

    // get the current center-of-mass 
    com_in = Tref.mult_point(iF.get_translation());

    // inertia for child is relative to reference frame; make it relative to
    // frame aligned with global and located at child c.o.m.
    Primitive::transform_inertia(link->get_mass(), link->get_inertia(), com_in, iR, J_out, com_out);

    // setup the link position and orientation
    link->set_position(com_out);
    link->set_orientation(Quat(&Rdes));     

    // compute the relative collision transform
    if (data.collision_transforms.find(link) != data.collision_transforms.end())
    {
      Matrix4 Tc_rel = mTr * data.collision_transforms.find(link)->second;
      assert(!link->geometries.empty());
      link->geometries.front()->set_rel_transform(Tc_rel, false); 
    }

    // compute the relative visual transform
    if (data.visual_transforms.find(link) != data.visual_transforms.end())
    {
      Matrix4 Tv_rel = mTr * data.visual_transforms.find(link)->second;
      #ifdef USE_OSG
      osg::MatrixTransform* group = (osg::MatrixTransform*) data.visual_transform_nodes.find(link)->second;
      osg::Matrixd m;
      to_osg_matrix(Tv_rel, m);
      group->setMatrix(m);
      #endif
    }

    // get the joint position and the two vectors we need (all in global frame) 
    Vector3 joint_position = Tref.get_translation();
    Vector3 joint_to_child_com = link->get_position() - joint_position;
    Vector3 parent_com_to_joint = joint_position - parent->get_position();

    // setup the pointers between the joints and links
    link->add_inner_joint(parent, joint, Tref.transpose_mult_vector(joint_to_child_com), link->get_transform().transpose_mult_vector(joint_to_child_com));
    parent->add_outer_joint(link, joint, parent->get_transform().transpose_mult_vector(parent_com_to_joint));

    // setup the joint axis, if necessary, and determine q_tare
    if (data.joint_axes.find(joint) != data.joint_axes.end())
    {
      Vector3 axis = Tref.mult_vector(data.joint_axes.find(joint)->second);
      shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
      shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
      if (pj)
      {
        pj->set_axis_global(axis);
        pj->determine_q_tare();
      }
      else if (rj)
      {
        rj->set_axis_global(axis);
        rj->determine_q_tare();
      }
    }
    
output_data(data, link);
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
void URDFReader::read_visual(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link)
{
  #ifdef USE_OSG
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
      link->geometries.push_back(cg);
      cg->set_single_body(link);
      cg->set_geometry(primitive);

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

  Quat rpy_quat = Quat::rpy(rpy[0], rpy[1], rpy[2]);
  Matrix3 rpy_matrix = Matrix3(&rpy_quat); 
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



