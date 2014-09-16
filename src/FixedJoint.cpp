/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/FixedJoint.h>

using namespace Ravelin;
using namespace Moby;
using std::vector;
using boost::shared_ptr;

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
FixedJoint::FixedJoint() : Joint(), FixedJointd()
{
  init_data();  // FixedJointd()'s constructor should call our virtual...
                // ... method but this is not happening. 
}

/// Initializes the joint with the specified inboard and outboard links
FixedJoint::FixedJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard), FixedJointd()
{
}  

/// Implements Base::load_from_xml()
void FixedJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "FixedJoint") == 0);

  // set the relative joint transformation
  FixedJointd::determine_q_tare();
}

/// Implements Base::save_to_xml()
void FixedJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "FixedJoint";
}

