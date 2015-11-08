/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/PlanarJoint.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using boost::dynamic_pointer_cast;
using namespace Moby;

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
PlanarJoint::PlanarJoint() : Joint(), PlanarJointd()
{
  init_data();  // PlanarJointd()'s constructor should call our virtual...
                // ... method but this is not happening. 
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
PlanarJoint::PlanarJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard), PlanarJointd()
{
}  

/// Implements Base::load_from_xml()
void PlanarJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "PlanarJoint") == 0);

  // read the plane normal, if given
  XMLAttrib* normal_attrib = node->get_attrib("normal");
  if (normal_attrib)
  {
    Vector3d normal;
    normal_attrib->get_vector_value(normal);
    set_normal(normal);  
  }

  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // compute _q_tare if necessary
  if (_determine_q_tare)
    PlanarJointd::determine_q_tare();
}

/// Implements Base::save_to_xml()
void PlanarJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "PlanarJoint";

  // save normal 
  node->attribs.insert(XMLAttrib("normal", _normal));
}

