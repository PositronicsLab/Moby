/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/SphericalJoint.h>

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
SphericalJoint::SphericalJoint() : Joint(), SphericalJointd()
{
  init_data();  // SphericalJointd()'s constructor should call our virtual...
                // ... method but this is not happening. 
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The axis of rotation is set to [0 0 0].
 */
SphericalJoint::SphericalJoint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard), SphericalJointd()
{
}  

/// Implements Base::load_from_xml()
void SphericalJoint::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "SphericalJoint") == 0);

  // read the global joint axes, if given
  XMLAttrib* axis_attrib = node->get_attrib("axis");
  if (axis_attrib)
  {
    Vector3d axis;
    axis_attrib->get_vector_value(axis);
    set_axis(axis, eAxis1);  
  }
  XMLAttrib* axis1_attrib = node->get_attrib("axis1");
  if (axis1_attrib)
  {
    Vector3d axis1;
    axis1_attrib->get_vector_value(axis1);
    set_axis(axis1, eAxis1);  
  }
  XMLAttrib* axis2_attrib = node->get_attrib("axis2");
  if (axis2_attrib)
  {
    Vector3d axis2;
    axis2_attrib->get_vector_value(axis2);
    set_axis(axis2, eAxis2);  
  }

  // compute _q_tare if necessary
  if (_determine_q_tare)
    SphericalJointd::determine_q_tare();
}

/// Implements Base::save_to_xml()
void SphericalJoint::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  Vector3d u0;

  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "SphericalJoint";

  // convert local axes to global axes and save
  u0 = Pose3d::transform_vector(shared_ptr<const Pose3d>(), _u[eAxis1]);
  node->attribs.insert(XMLAttrib("axis1", u0));
  u0 = Pose3d::transform_vector(shared_ptr<const Pose3d>(), _u[eAxis2]);
  node->attribs.insert(XMLAttrib("axis2", u0));
}

