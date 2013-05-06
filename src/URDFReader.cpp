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

#include <Ravelin/AAngled.h>
#include <Moby/CSG.h>
#include <Moby/Primitive.h>
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
#include <Moby/DeformableCCD.h>
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
static void to_osg_matrix(const Opsd::Pose3d& src, osg::Matrixd& tgt)
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

/// Fix for Moby
void URDFReader::fix_Moby(URDFData& data, const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  std::map<RigidBodyPtr, RigidBodyPtr> parents;
  Matrix3d J_out;
  Point3d com_out;

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

  // add all children of the base link to the queue 
  queue<RigidBodyPtr> q;
  find_children(data, base, q, parents);

  // process from base outward 
  while (!q.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = q.front();
    q.pop();

    // get the parent link and the joint
    JointPtr joint = find_joint(data, link);
    RigidBodyPtr parent = data.joint_parent.find(joint)->second;

    // see whether this relative transform is ok
    if (!valid_transformation(data, parent, joint, link))
    {
      // correct the transform; first, determine the desired outboard transform
      // (relative to inboard)
      const Pose3d& To = data.joint_transforms.find(joint)->second;
      Pose3d Tx = To;
      Tx.q.set_identity();

      // now compute the necessary transformation and its inverse
      Pose3d oTx = Pose3d::inverse(To) * Tx;
      Pose3d xTo = Pose3d::inverse(oTx);

      std::cerr << "URDFReader warning! Changing transforms! " << std::endl;
      std::cerr << "  joint transform (" << joint->id << ") was: " << std::endl << To;
      std::cerr << "  now: " << std::endl << Tx;
      std::cerr << "  inertial transform (" << link->id << ") was: " << std::endl << data.inertia_transforms[link];
      std::cerr << "  now: " << std::endl << (data.inertia_transforms[link] * oTx);
      std::cerr << "  visual transform (" << link->id << ") was: " << std::endl << data.visual_transforms[link];
      std::cerr << "  now: " << std::endl << (data.visual_transforms[link] * oTx);
      std::cerr << "  collision transform (" << link->id << ") was: " << std::endl << data.collision_transforms[link];
      std::cerr << "  now: " << std::endl << (data.collision_transforms[link] * oTx);

      // update all transformations / axes of the outboard 
      data.joint_transforms[joint] = data.joint_transforms[joint] * oTx;
      data.inertia_transforms[link] = data.inertia_transforms[link] * oTx;
      data.visual_transforms[link] = data.visual_transforms[link] * oTx;
      data.collision_transforms[link] = data.collision_transforms[link] * oTx;
    
      // update joint axis, if it is present
      if (data.joint_axes.find(joint) != data.joint_axes.end())
        data.joint_axes[joint] = xTo.transform(data.joint_axes[joint]);

      // now update *just* the joint transformation
      vector<pair<JointPtr, RigidBodyPtr> > outboards;
      find_outboards(data, link, outboards, parents);
      for (unsigned i=0; i< outboards.size(); i++)
      {
        // get the i'th outboard and its inner joint
        JointPtr joint = outboards[i].first;
        RigidBodyPtr outboard = outboards[i].second;

        // update the transform
        std::cerr << "  child joint transform (" << joint->id << ") was: " << std::endl << data.joint_transforms[joint];
        std::cerr << "  now: " << std::endl << (xTo * data.joint_transforms[joint]);
        data.joint_transforms[joint] = xTo * data.joint_transforms[joint];
      }
    }

    // add all children to the queue for processing
    find_children(data, link, q, parents);
  }

  // ok, now transform all inertias so that their orientation is unchanged
  // from the reference frame
  // start with the base
  Pose3d& jTref = data.inertia_transforms.find(base)->second;
  Matrix3d jRrefT = Matrix3d::transpose(Matrix3d(jTref.q));
  Point3d jxref = jTref.x;
  Primitive::transform_inertia(base->get_mass(), base->get_inertia(), jxref, jRrefT, J_out, com_out); 
  jTref.x = com_out;
  jTref.q.set_identity();
  base->set_inertia(J_out);  

  // now process all children
  find_children(data, base, q, parents);

  // process from base outward 
  while (!q.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = q.front();
    q.pop();

    // update the inertia
    Pose3d& jTref = data.inertia_transforms.find(link)->second;
    Matrix3d jRrefT = Matrix3d::transpose(Matrix3d(jTref.q));
    Point3d jxref = jTref.x;
    Primitive::transform_inertia(link->get_mass(), link->get_inertia(), jxref, jRrefT, J_out, com_out); 
    jTref.q.set_identity();
    jTref.x = com_out;
    link->set_inertia(J_out);  

    // now process all children
    find_children(data, link, q, parents);
  }
}

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

/// Determine whether it is possible to get from one reference frame to another
bool URDFReader::valid_transformation(const URDFData& data, RigidBodyPtr inboard, JointPtr joint, RigidBodyPtr outboard)
{
  // get the inboard reference frame (if it's available- it won't be if the
  // inboard is the base)
  Pose3d Tref_inboard = IDENTITY_4x4;
  for (map<JointPtr, RigidBodyPtr>::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    if (i->second == inboard)
    {
      Tref_inboard = data.joint_transforms.find(i->first)->second;
      break;
    } 

  // get the outboard reference frame
  assert(data.joint_transforms.find(joint) != data.joint_transforms.end());
  Pose3d Tref_outboard = data.joint_transforms.find(joint)->second;

  // get the relative transform and relative rotation
  Pose3d Trel = Pose3d::inverse(Tref_inboard) * Tref_outboard;
  Matrix3d Rrel = Trel.q;

  // see whether we can get from Tref_inboard to Tref_outboard based on
  // joint movement alone
  if (dynamic_pointer_cast<RevoluteJoint>(joint))
  {
    // get the joint axis relative to inboard frame
    Vector3d axis_0 = Tref_outboard.transform(data.joint_axes.find(joint)->second);
    Vector3d axis = Tref_inboard.inverse_transform(axis_0);

    // verify that rotation around that axis gives the desired rotation
    AAngled aa(Rrel, axis);
    Matrix3d Raa(aa);
    return std::fabs(Quatd::calc_angle(Quatd(Rrel), Quatd(Raa))) < NEAR_ZERO;    
  }
  else if (dynamic_pointer_cast<PrismaticJoint>(joint) || dynamic_pointer_cast<FixedJoint>(joint))
  {
    // both orientations must be the same
    Quatd qd(Rrel);
    Quatd eye = Quatd::identity();
    return std::fabs(Quatd::calc_angle(qd, eye)) < NEAR_ZERO;
  }
  else
    // should never be here!
    assert(false);

  // req'd by compiler...
  return false; 
}

void URDFReader::output_data(const URDFData& data, RigidBodyPtr link)
{
  #ifdef DEBUG_URDF
  std::cout << "link id: " << link->id << std::endl;
  JointPtr joint = find_joint(data, link);
  if (joint)
  {
    std::cout << "  URDF joint transform: " << std::endl << data.joint_transforms.find(joint)->second << std::endl;
    if (data.joint_axes.find(joint) != data.joint_axes.end())
      std::cout << "  URDF joint axis: " << data.joint_axes.find(joint)->second << std::endl;
    std::cout << "  Moby joint position: " << joint->get_position_global() << std::endl;
    shared_ptr<RevoluteJoint> rj = dynamic_pointer_cast<RevoluteJoint>(joint);
    shared_ptr<PrismaticJoint> pj = dynamic_pointer_cast<PrismaticJoint>(joint);
    if (rj)
      std::cout << "  Moby joint axis: " << rj->get_axis_global() << std::endl;
    else if (pj)
      std::cout << "  Moby joint axis: " << pj->get_axis_global() << std::endl;
  }
  if (data.visual_transforms.find(link) != data.visual_transforms.end())
    std::cout << "  URDF visualization transform: " << std::endl << data.visual_transforms.find(link)->second;
  if (data.inertia_transforms.find(link) != data.inertia_transforms.end())
    std::cout << "  URDF inertia transform: " << std::endl << data.inertia_transforms.find(link)->second;
  if (data.collision_transforms.find(link) != data.collision_transforms.end())
    std::cout << "  URDF collision transform: " << std::endl << data.collision_transforms.find(link)->second;
  std::cout << "  Moby transform: " << std::endl << *link->get_transform();
  if (!link->geometries.empty())
  {
    std::cout << "  Moby relative collision transform: " << std::endl << link->geometries.front()->get_rel_transform();
    std::cout << "  Moby collision transform: " << std::endl << (*link->get_transform() * link->geometries.front()->get_rel_transform());
  }
  if (data.visual_transform_nodes.find(link) != data.visual_transform_nodes.end())
  {
    Pose3d Tv;
    #ifdef USE_OSG
    from_osg_matrix(((osg::MatrixTransform*) data.visual_transform_nodes.find(link)->second)->getMatrix(), Tv);
    #endif
    std::cout << "  Moby relative visual transform: " << std::endl << Tv;
    std::cout << "  Moby visual transform: " << std::endl << (*link->get_transform() * Tv);
  }
  #endif
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
  Matrix3d J_out;
  Point3d com_out;

  // output all data (before Moby fixes it)
  #ifdef DEBUG_URDF 
  std::cout << "data (pre-fixes)" << std::endl;
//  for (unsigned i=0; i< links.size(); i++)
//    output_data(data, links[i]);
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "data (post-fixes)" << std::endl;
  #endif

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
  Pose3d baseT = data.inertia_transforms.find(base)->second; 
  Point3d com_in = baseT.x;
  baseT.x.set_zero();
  Primitive::transform_inertia(base->get_mass(), base->get_inertia(), com_in, baseT, J_out, com_out); 
  base->set_inertia(J_out);

  // update base position and orientation
  base->set_position(com_out);
  base->set_orientation(Quatd::identity());

  // compute the relative collision transform
  if (data.collision_transforms.find(base) != data.collision_transforms.end())
  {
    // get the true frame of the collision transform
    const Pose3d& Trel = data.collision_transforms.find(base)->second;
    Pose3d Tc = baseT * Trel;

    // now compute the collision transform relative to the base frame
    Pose3d Tc_rel = Pose3d::inverse(*base->get_transform()) * Tc;
    assert(!base->geometries.empty());
    base->geometries.front()->set_rel_transform(Tc_rel, false); 
  }

  // compute the relative visual transform
  if (data.visual_transforms.find(base) != data.visual_transforms.end())
  {
    // get the true frame of the visual transform
    const Pose3d& Trel = data.visual_transforms.find(base)->second;
    Pose3d Tv = baseT * Trel;

    // now compute the collision transform relative to the base frame
    Pose3d Tv_rel = Pose3d::inverse(*base->get_transform()) * Tv;
    #ifdef USE_OSG
    osg::MatrixTransform* group = (osg::MatrixTransform*) data.visual_transform_nodes.find(base)->second;
    osg::Matrixd m;
    to_osg_matrix(Tv_rel, m);
    group->setMatrix(m);
    #endif
  }

  // output data for the base
  output_data(data, base);

  // add all children to the queue  
  queue<RigidBodyPtr> q;
  find_children(data, base, q, parents);
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
    const Pose3d& refJc = data.joint_transforms.find(joint)->second;

    // compute the reference frame for the parent
    Pose3d refPj = data.inertia_transforms.find(parent)->second;
    Pose3d oPj = *parent->get_transform();
    Pose3d oPref = oPj * Pose3d::inverse(refPj);

    // compute the *reference* frame (w.r.t. the global frame) for the child 
    // (and the joint)
    Pose3d oTref = oPref * refJc;
    #ifdef DEBUG_URDF
    std::cout << "link " << link->id << " reference frame: " << std::endl << oTref;
    #endif

    // compute the desired frame for the child
    Pose3d refTj = data.inertia_transforms.find(link)->second;
    Pose3d oTj = oTref * refTj; 

    // compute the relative collision transform
    if (data.collision_transforms.find(link) != data.collision_transforms.end())
    {
      Pose3d refTc = data.collision_transforms.find(link)->second;
      Pose3d jTc = Pose3d::inverse(refTj) * refTc;
      assert(!link->geometries.empty());
      link->geometries.front()->set_rel_transform(jTc, false); 
    }

    // compute the relative visual transform
    if (data.visual_transforms.find(link) != data.visual_transforms.end())
    {
      Pose3d refTv = data.visual_transforms.find(link)->second;
      Pose3d jTv = Pose3d::inverse(refTj) * refTv;
      #ifdef USE_OSG
      osg::MatrixTransform* group = (osg::MatrixTransform*) data.visual_transform_nodes.find(link)->second;
      osg::Matrixd m;
      to_osg_matrix(jTv, m);
      group->setMatrix(m);
      #endif
    }

    // setup the link position and orientation
    link->set_transform(oTj);

    // get the joint position and the two vectors we need (all in global frame) 
    Point3d joint_position = oTref.x;
    Vector3d joint_to_child_com = oTj.x - joint_position;
    Vector3d parent_com_to_joint = joint_position - parent->get_transform()->x;

    // setup the pointers between the joints and links
    link->add_inner_joint(parent, joint, joint_to_child_com, link->get_transform()->inverse_transform(joint_to_child_com));
    parent->add_outer_joint(link, joint, parent->get_transform()->inverse_transform(parent_com_to_joint));

    // setup the joint axis, if necessary, and determine q_tare
    if (data.joint_axes.find(joint) != data.joint_axes.end())
    {
      Vector3d axis = oTref.transform(data.joint_axes.find(joint)->second);
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

    // output data for the link    
    output_data(data, link);

    // add all children of link to the queue
    find_children(data, link, q, parents);
  }

  return true;
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
RigidBodyPtr URDFReader::read_parent(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links)
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
RigidBodyPtr URDFReader::read_child(shared_ptr<const XMLTree> node, URDFData& data, const vector<RigidBodyPtr>& links)
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
void URDFReader::read_dynamics(shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint)
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
bool URDFReader::read_color(shared_ptr<const XMLTree> node, URDFData& data, VectorNd& color)
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
bool URDFReader::read_texture(shared_ptr<const XMLTree> node, URDFData& data, string& texture_fname)
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
      const XMLAttrib* name_attr = (*i)->get_attrib("name");
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
void URDFReader::read_collision(shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link)
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
      else if ((primitive = read_cylinder(*i, data)))
        return primitive;
      else if ((primitive = read_sphere(*i, data)))
        return primitive;
      else if ((primitive = read_trimesh(*i, data)))
        return primitive;
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
      const XMLAttrib* size_attrib = (*i)->get_attrib("size");
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

/// Reads a trimesh primitive
shared_ptr<TriangleMeshPrimitive> URDFReader::read_trimesh(shared_ptr<const XMLTree> node, URDFData& data)
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
shared_ptr<CylinderPrimitive> URDFReader::read_cylinder(shared_ptr<const XMLTree> node, URDFData& data)
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

/// Reads a sphere primitive
shared_ptr<SpherePrimitive> URDFReader::read_sphere(shared_ptr<const XMLTree> node, URDFData& data)
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
      double radius = radius_attrib->get_real_value();
      return shared_ptr<SpherePrimitive>(new SpherePrimitive(radius));
    }
  }

  return shared_ptr<SpherePrimitive>();
}


/// Attempts to read an "origin" tag
Pose3d URDFReader::read_origin(shared_ptr<const XMLTree> node, URDFData& data)
{
  Point3d xyz = ZEROS_3;
  Vector3d rpy = ZEROS_3;

  // look for the tag
  const list<XMLTreePtr>& child_nodes = node->children;
  for (list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "origin") == 0)
    {
      // look for xyz attribute 
      const XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (xyz_attrib)
        xyz = xyz_attrib->get_point_value();

      // look for rpy attribute
      const XMLAttrib* rpy_attrib = (*i)->get_attrib("rpy");
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



