/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <stack>
#include <queue>
#include <Moby/Log.h>
#include <Moby/Joint.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/XMLReader.h>
#include <Moby/RCArticulatedBody.h>

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::queue;
using std::list;
using std::map;
using std::string;
using namespace Ravelin;
using namespace Moby;

/// Default constructor
/**
 * Constructs a reduced-coordinate articulated body with no joints and no links.
 */
RCArticulatedBody::RCArticulatedBody()
{
}

/// Clones this
RCArticulatedBodyPtr RCArticulatedBody::clone() const 
{
  // create a node for Moby
  XMLTreePtr node(new XMLTree("Moby"));

  // setup a list of shared objects 
  std::list<shared_ptr<const Base> > shared_objects;
  shared_objects.push_back(dynamic_pointer_cast<const Base>(RCArticulatedBodyd::shared_from_this()));

  // init a set of serialized objects
  std::set<shared_ptr<const Base> > serialized;

  // develop the XML tree until there is nothing more to serialize
  while (!shared_objects.empty())
  {
    // get the object off of the front of the queue
    shared_ptr<const Base> obj = shared_objects.front();
    assert(obj);
    shared_objects.pop_front();

    // if this object has already been serialized, skip it
    if (serialized.find(obj) != serialized.end())
      continue;

    // create a new node for this object under the parent
    XMLTreePtr new_node(new XMLTree(""));
    node->add_child(new_node);

    // serialize to this new node
    obj->save_to_xml(new_node, shared_objects);

    // indicate that the node has been serialized
    serialized.insert(obj);
  }

  // now read in all of the objects
  map<string, BasePtr> objects = XMLReader::construct_ID_map(node);

  // find the RCArticulatedBody object - verify that there is only one
  for (map<string, BasePtr>::const_iterator i = objects.begin(); i != objects.end(); i++)
  {
    RCArticulatedBodyPtr rcab = dynamic_pointer_cast<RCArticulatedBody>(i->second);

    if (rcab)
      return rcab;
  }

  // should not still be here
  assert(false);
  return RCArticulatedBodyPtr(); 
}

/// Compiles this body (updates the link transforms and velocities)
void RCArticulatedBody::compile()
{
  // call parent methods first
  ArticulatedBody::compile();
  RCArticulatedBodyd::compile();

  // verify all links are enabled
  if (!is_floating_base())
  {
    for (unsigned i=1; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("Only first link can be disabled in a reduced coordinate body with fixed-base!");
  }
  else
  {
    for (unsigned i=0; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("No links can be disabled in a reduced coordinate body with floating-base!");
  }

  // update processed vector size
  _processed.resize(_links.size());

  // setup explicit joint generalized coordinate and constraint indices
  const vector<shared_ptr<Jointd> >& ejoints = get_explicit_joints();
  for (unsigned i=0, cidx = 0, ridx = 0; i< ejoints.size(); i++)
  {
    ejoints[i]->set_coord_index(cidx);
    ejoints[i]->set_constraint_index(ridx);
    cidx += ejoints[i]->num_dof();
    ridx += ejoints[i]->num_constraint_eqns();
  }

  // setup explicit joint generalized coordinate and constraint indices
  const vector<shared_ptr<Jointd> >& ijoints = get_implicit_joints();
  for (unsigned i=0, cidx = 0, ridx=0; i< ijoints.size(); i++)
  {
    ijoints[i]->set_coord_index(cidx);
    ijoints[i]->set_constraint_index(ridx);
    cidx += ijoints[i]->num_dof();
    ridx += ijoints[i]->num_constraint_eqns();
  }

  // update link transforms and velocities
  update_link_poses();
  update_link_velocities();
}

/// Sets the vector of links and joints
void RCArticulatedBody::set_links_and_joints(const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  // call the parent method to update the link indices, etc.
  vector<shared_ptr<RigidBodyd> > d_links;
  vector<shared_ptr<Jointd> > d_joints;
  d_links.insert(d_links.end(), links.begin(), links.end());
  d_joints.insert(d_joints.end(), joints.begin(), joints.end());
  RCArticulatedBodyd::set_links_and_joints(d_links, d_joints);
}

/// Implements Base::load_from_xml()
/**
 * \pre all links and joints must be loaded using their respective
 *      serialization methods before this method is called
 */
void RCArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // load the parent data
  ArticulatedBody::load_from_xml(node, id_map);

  // don't verify the node name -- this class has derived classes
//  assert(strcasecmp(node->name.c_str(), "RCArticulatedBody") == 0);

  // get whether the body has a floating base
  XMLAttrib* fb_attr = node->get_attrib("floating-base");
  if (fb_attr)
    _floating_base = fb_attr->get_bool_value();

  // read the pointer to the forward dynamics algorithm, if provided
  XMLAttrib* fdyn_algo_attr = node->get_attrib("fdyn-algorithm");
  if (fdyn_algo_attr)
  {
    // get the ID
    string algo = fdyn_algo_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = algo.find_first_not_of(" \t\n\r");
    size_t last_nws_index = algo.find_last_not_of(" \t\n\r");
    algo = algo.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the algorithm type
    if (strcasecmp(algo.c_str(), "fsab") == 0)
      algorithm_type = eFeatherstone;
    else if (strcasecmp(algo.c_str(), "crb") == 0)
      algorithm_type = eCRB;
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "forward dynamics algorithm " << std::endl << "  type '";
      std::cerr << algo << "' -- valid types are 'fsab' and 'crb'";
      std::cerr << std::endl;
    }
  }

  // read the forward dynamics algorithm computation frame, if provided
  XMLAttrib* fdyn_frame_attr = node->get_attrib("fdyn-algorithm-frame");
  if (fdyn_frame_attr)
  {
    // get the computation reference frame as a string
    string frame = fdyn_frame_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = frame.find_first_not_of(" \t\n\r");
    size_t last_nws_index = frame.find_last_not_of(" \t\n\r");
    frame = frame.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the frame type
    if (strcasecmp(frame.c_str(), "global") == 0)
      set_computation_frame_type(eGlobal);
    else if (strcasecmp(frame.c_str(), "link") == 0)
      set_computation_frame_type(eLink);
    else if (strcasecmp(frame.c_str(), "linkinertia") == 0 || strcasecmp(frame.c_str(), "link-inertia") == 0)
      set_computation_frame_type(eLinkInertia);
    else if (strcasecmp(frame.c_str(), "linkcom") == 0 || strcasecmp(frame.c_str(), "link-com") == 0)
      set_computation_frame_type(eLinkCOM);
    else if (strcasecmp(frame.c_str(), "joint") == 0)
      set_computation_frame_type(eJoint);
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "computation reference " << std::endl << "  frame type '";
      std::cerr << frame << "' -- valid types are 'link' and 'global'";
      std::cerr << std::endl;
    }
  }

  // compile everything once again, for safe measure
  compile();

  // transform the body, if desired
  XMLAttrib* xlat_attr = node->get_attrib("translate");
  XMLAttrib* rotate_attr = node->get_attrib("rotate");

  if (rotate_attr)
  {
    shared_ptr<RigidBodyd> base = get_base_link();
    base->rotate(rotate_attr->get_rpy_value());
    update_link_poses();
  }
  if (xlat_attr)
  {
    shared_ptr<RigidBodyd> base = get_base_link();
    base->translate(xlat_attr->get_origin_value());
    update_link_poses();
  }
}

/// Implements Base::save_to_xml()
void RCArticulatedBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent method first
  ArticulatedBody::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "RCArticulatedBody";

  // write whether body has a floating base
  node->attribs.insert(XMLAttrib("floating-base", _floating_base));

  // write the forward dynamics algorithm type
  if (algorithm_type == eFeatherstone)
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("fsab")));
  else
  {
    assert(algorithm_type == eCRB);
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("crb")));
  }

  // write the forward dynamics algorithm frame -- note that the string()
  // is necessary on the second argument to XMLAttrib b/c the compiler
  // interprets a constant string as a bool, rather than as an string,
  // given a choice
  if (get_computation_frame_type() == eGlobal)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("global")));
  else if (get_computation_frame_type() == eLink)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("link")));
  else if (get_computation_frame_type() == eLinkCOM)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("linkcom")));
  else if (get_computation_frame_type() == eLinkInertia)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("linkinertia")));
  else
  {
    assert(get_computation_frame_type() == eJoint);
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("joint")));
  }
}

