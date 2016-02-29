/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iostream>
#include <boost/foreach.hpp>
#include <queue>
#include <Moby/XMLTree.h>
#include <Moby/Joint.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/NumericalException.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/InvalidStateException.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/FixedJoint.h>
#include <Moby/UniversalJoint.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/URDFReader.h>

using namespace Moby;
using namespace Ravelin;
using std::set;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::list;
using std::vector;
using std::map;
using std::string;
using std::queue;

ArticulatedBody::ArticulatedBody()
{
}

/// Integrates a dynamic body
/*
void ArticulatedBody::integrate(double t, double h, shared_ptr<Integrator> integrator)
{
  FILE_LOG(LOG_DYNAMICS) << "ArticulatedBody::integrate() - integrating from " << t << " by " << h << std::endl;

  // reset the acceleration events exceeded
  _vel_limit_exceeded = false;

  if (_kinematic_update)
  {
    FILE_LOG(LOG_DYNAMICS) << " -- body set to kinematic update --" << std::endl;
    if (controller)
      (*controller)(dynamic_pointer_cast<ArticulatedBody>(shared_from_this()), t, controller_arg);
    return;
  }

  shared_ptr<ArticulatedBody> shared_this = dynamic_pointer_cast<ArticulatedBody>(shared_from_this());

  get_generalized_coordinates(eEuler, gc);
  get_generalized_velocity(eSpatial, gv); // gv depends on gc
  gcgv.resize(gc.size()+gv.size());
  gcgv.set_sub_vec(0, gc);
  gcgv.set_sub_vec(gc.size(), gv);
  integrator->integrate(gcgv, &ArticulatedBody::ode_both, t, h, (void*) &shared_this);
  gcgv.get_sub_vec(0, gc.size(), gc);
  gcgv.get_sub_vec(gc.size(), gcgv.size(), gv);
  // NOTE: velocity must be set first (it's computed w.r.t. old frame)
  set_generalized_coordinates(eEuler, gc);
  set_generalized_velocity(eSpatial, gv);
}
*/

/// Prepares to compute the ODE  
void ArticulatedBody::prepare_to_calc_ode_sustained_constraints(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the shared pointer to this
  ArticulatedBodyPtr shared_this = dynamic_pointer_cast<ArticulatedBody>(ArticulatedBodyd::shared_from_this());

  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the coordinates and velocity from x
  SharedConstVectorNd gc = x.segment(0, NGC_EUL);
  SharedConstVectorNd gv = x.segment(NGC_EUL, x.size());

  // set the state
  set_generalized_coordinates_euler(gc);
  if (is_joint_constraint_violated())
    throw InvalidStateException();

  // set the velocity 
  set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
  {
    VectorNd tmp;

    // get the generalized forces
    (*controller)(get_this(), tmp, t, controller_arg);

    FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << id << std::endl;

    // apply the generalized forces
    add_generalized_force(tmp);
  }
}

/// Prepares to compute the ODE  
void ArticulatedBody::prepare_to_calc_ode(SharedConstVectorNd& x, double t, double dt, void* data)
{
  // get the shared pointer to this
  ArticulatedBodyPtr shared_this = dynamic_pointer_cast<ArticulatedBody>(ArticulatedBodyd::shared_from_this());

  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the coordinates and velocity from x
  SharedConstVectorNd gc = x.segment(0, NGC_EUL);
  SharedConstVectorNd gv = x.segment(NGC_EUL, x.size());

  // set the state
  set_generalized_coordinates_euler(gc);
  if (is_joint_constraint_violated())
    throw InvalidStateException();

  // set the velocity 
  set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // clear the force accumulators on the body
  reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(shared_this);

  // call the body's controller
  if (controller)
  {
    VectorNd tmp;

    // get the generalized forces
    (*controller)(get_this(), tmp, t, controller_arg);

    FILE_LOG(LOG_DYNAMICS) << "Computing controller forces for " << id << std::endl;

    // apply the generalized forces
    add_generalized_force(tmp);
  }
}

/// Returns the ODE's for position and velocity (concatenated into x)
void ArticulatedBody::ode(double t, double dt, void* data, SharedVectorNd& dx)
{
  // get the articulated body
  const unsigned NGC_EUL = num_generalized_coordinates(eEuler);

  // get the derivatives of coordinates and velocity from dx
  SharedVectorNd dgc = dx.segment(0, NGC_EUL);
  SharedVectorNd dgv = dx.segment(NGC_EUL, dx.size());

  // we need the generalized velocity as Rodrigues coordinates
  get_generalized_velocity(DynamicBodyd::eEuler, dgc);

  // calculate forward dynamics 
  get_generalized_acceleration(dgv);
}

// TODO: remove this
/// Updates joint constraint violation (after integration)
void ArticulatedBody::update_joint_constraint_violations()
{
  // update the size of the constraint violation vector
  _cvio.resize(num_joint_dof());

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    JointPtr joint = dynamic_pointer_cast<Joint>(_joints[i]);

    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      _cvio[k] = 0.0;
      if (joint->q[j] < joint->lolimit[j])
        _cvio[k] = joint->lolimit[j] - joint->q[j];
      else if (joint->q[j] > joint->hilimit[j])
        _cvio[k] = joint->q[j] - joint->hilimit[j]; 
    }
  }
}

/// Checks for a joint constraint violation
bool ArticulatedBody::is_joint_constraint_violated() const
{
  // obvious check
  if (_cvio.size() != num_joint_dof())
    return false;

  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    JointPtr joint = dynamic_pointer_cast<Joint>(_joints[i]);

    // loop over all DOF
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      if (joint->q[j] + _cvio[k] < joint->lolimit[j] ||
          joint->q[j] - _cvio[k] > joint->hilimit[j])
        return true; 
    }
  }

  // no violation if here
  return false;
}

/// Updates visualization for the body
void ArticulatedBody::update_visualization()
{
  BOOST_FOREACH(shared_ptr<RigidBodyd> rb, _links)
    dynamic_pointer_cast<RigidBody>(rb)->update_visualization(); 
  BOOST_FOREACH(shared_ptr<Jointd> joint, _joints)
    dynamic_pointer_cast<Joint>(joint)->update_visualization(); 
}

/// Loads a MCArticulatedBody object from an XML node
void ArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, std::map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // call parent method
  ControlledBody::load_from_xml(node, id_map);

  // save the id
  body_id = id;

  // don't verify the node name -- this class has derived classes
  // assert(strcasecmp(node->name().c_str(), "MCArticulatedBody") == 0);

  // see whether to load the model from a URDF file
  XMLAttrib* urdf_attr = node->get_attrib("urdf-filename");
  if (urdf_attr)
  {
    // get the URDF filename
    std::string urdf_fname = urdf_attr->get_string_value();

    // load robots from the URDF
    std::string robot_name;
    std::vector<shared_ptr<RigidBody> > links;
    std::vector<shared_ptr<Joint> > joints; 
    std::vector<shared_ptr<RigidBodyd> > linksd;
    std::vector<shared_ptr<Jointd> > jointsd; 
    if (URDFReader::read(urdf_fname, robot_name, links, joints))
    {
      linksd.insert(linksd.end(), links.begin(), links.end());
      jointsd.insert(jointsd.end(), joints.begin(), joints.end());
      set_links_and_joints(linksd, jointsd);
    }
    else
      std::cerr << "ArticulatedBody::load_from_xml()- unable to process URDF " << urdf_fname << std::endl;

    // do no more processing
    return;
  }

  // setup a list of joint nodes to find
  list<string> joint_node_names;
  joint_node_names.push_back("RevoluteJoint");
  joint_node_names.push_back("PrismaticJoint");
  joint_node_names.push_back("SphericalJoint");
  joint_node_names.push_back("UniversalJoint");
  joint_node_names.push_back("FixedJoint");
  joint_node_names.push_back("JointPlugin");
  joint_node_names.push_back("Gears");

  // read the set of joint nodes and concatenate them into a single list
  list<shared_ptr<const XMLTree> > joint_nodes = node->find_child_nodes(joint_node_names);
  
  // read the set of link nodes
  list<shared_ptr<const XMLTree> > link_nodes = node->find_child_nodes("RigidBody");

  // if there were links read or joints read, add them 
  if (!joint_nodes.empty() || !link_nodes.empty())
  {
    // setup new lists for joints and links
    list<JointPtr> joints;
    list<RigidBodyPtr> links;

    // process all link nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = link_nodes.begin(); i != link_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body links are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the link was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the link
      links.push_back(dynamic_pointer_cast<RigidBody>(id_iter->second));
    }

    // process all joint nodes in the same manner
    for (list<shared_ptr<const XMLTree> >::const_iterator i = joint_nodes.begin(); i != joint_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body joints are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the joint was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the joints
      joints.push_back(dynamic_pointer_cast<Joint>(id_iter->second));
    }

    // set the joints and links
    set_links_and_joints(vector<shared_ptr<RigidBodyd> >(links.begin(), links.end()), vector<shared_ptr<Jointd> >(joints.begin(), joints.end()));
  }
}

/// Saves this object to a XML tree
void ArticulatedBody::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call parent method
  ControlledBody::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "ArticulatedBody";

  // add all links
  for (unsigned i=0; i< _links.size(); i++)
  {
    // create the new node in the tree
    XMLTreePtr child_node(new XMLTree("RigidBody"));
    node->add_child(child_node);

    // write to this node
    dynamic_pointer_cast<RigidBody>(_links[i])->save_to_xml(child_node, shared_objects);
  }

  // add all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // create the new node in the tree -- note that save_to_xml() should set
    // the node name correctly
    XMLTreePtr child_node(new XMLTree("Joint"));
    node->add_child(child_node);

    // write to this node
    dynamic_pointer_cast<Joint>(_joints[i])->save_to_xml(child_node, shared_objects);
  }
}

