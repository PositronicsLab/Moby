/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/XMLTree.h>
#include <Moby/Integrator.h>

using namespace Moby;

/// Implements Base::save_to_xml()
void Integrator::save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // set the node's name
  node->name = "Integrator";
}


