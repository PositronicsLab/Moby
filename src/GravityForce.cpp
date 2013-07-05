/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/GravityForce.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// Constructs a default gravity vector of [0,0,0]
GravityForce::GravityForce()
{
  gravity.set_zero();
}

/// Copy constructor
GravityForce::GravityForce(const GravityForce& source)
{
  operator=(source);
}

/// Adds gravity to a body
void GravityForce::add_force(DynamicBodyPtr body)
{
  // check to see whether body is a single body first 
  shared_ptr<SingleBody> sb = dynamic_pointer_cast<SingleBody>(body);
  if (sb)
  {
    SForced w;
    Vector3d gx = Pose3d::transform(GLOBAL, sb->get_computation_frame(), gravity);
    w.set_force(gx * sb->get_mass());
    w.pose = gx.pose; 
    sb->add_force(w);        
  }
  else
  {
    // it's an articulated body, get it as such
    ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);
      
    // get the vector of links
    const std::vector<RigidBodyPtr>& links = ab->get_links();
      
    // apply gravity force to all links
    BOOST_FOREACH(RigidBodyPtr rb, links)
    {
      SForced w;
      Vector3d gx = Pose3d::transform(GLOBAL, rb->get_computation_frame(), gravity);
      w.set_force(gx * rb->get_mass());
      w.pose = gx.pose; 
      rb->add_force(w);        
    }
  }
}

/// Implements Base::load_from_xml()
void GravityForce::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // load XML data for the parent class
  RecurrentForce::load_from_xml(node, id_map);

  // verify that the name of this node is correct
  assert(strcasecmp(node->name.c_str(), "GravityForce") == 0);

  // read the acceleration due to gravity, if given
  const XMLAttrib* gravity_attrib = node->get_attrib("accel");
  if (gravity_attrib)
    gravity_attrib->get_vector_value(gravity);
}

/// Implements Base::save_to_xml()
void GravityForce::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const 
{
  // save XML data from the parent class
  RecurrentForce::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "GravityForce";

  // save the acceleration due to gravity 
  node->attribs.insert(XMLAttrib("accel", gravity));
}
 
