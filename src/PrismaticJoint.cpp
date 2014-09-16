/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/RigidBody.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/PrismaticJoint.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The axis of rotation is set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
PrismaticJoint::PrismaticJoint() : Joint(), PrismaticJointd()
{
  init_data();  // PrismaticJointd()'s constructor should call our virtual...
                // ... method but this is not happening. 
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
PrismaticJoint::PrismaticJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
}  

/// Implements Base::load_from_xml()
void PrismaticJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "PrismaticJoint") == 0);

  // read the joint axis
  XMLAttrib* axis_attrib = node->get_attrib("axis");
  if (axis_attrib)
  {
    Vector3d axis;
    axis_attrib->get_vector_value(axis);
    set_axis(axis);
  }

  // set the joint tare
  if (_determine_q_tare)
    PrismaticJointd::determine_q_tare();
}

/// Implements Base::save_to_xml()
void PrismaticJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get the majority of the info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "PrismaticJoint";

  // save the joint axis (global coords)
  Vector3d u0 = Pose3d::transform_vector(shared_ptr<const Pose3d>(), _u);
  node->attribs.insert(XMLAttrib("axis", u0));
}

