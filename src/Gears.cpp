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
#include <Moby/Gears.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The gear ratio is set to 1:1 
 */
Gears::Gears() : Joint()
{
  init_data();  
  _ratio = 1.0;
}

/// Initializes the joint with the specified inboard and outboard links
/**
 * The gear ratio is set to 1:1 
 */
Gears::Gears(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard) : Joint(inboard, outboard)
{
  _ratio = 1.0;
}  

/// Implements Base::load_from_xml()
void Gears::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "Gears") == 0);

  // read the gear ratio, if given
  XMLAttrib* gear_ratio_attrib = node->get_attrib("gear-ratio");
  if (gear_ratio_attrib)
    gear_ratio_attrib->get_real_value(_ratio);
}

/// Implements Base::save_to_xml()
void Gears::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "Gears";

  // write the gear ratio 
  node->attribs.insert(XMLAttrib("gear-ratio", _ratio));
}

