/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/XMLTree.h>
#include <Moby/Integrator.h>

namespace Moby {

/// Implements Base::save_to_xml()
template <>
void Integrator<VectorN>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // set the node's name
  node->name = "Integrator";
}

/// Implements Base::save_to_xml()
template <>
void Integrator<Vector3>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // set the node's name
  node->name = "Integrator";
}

/// Implements Base::save_to_xml()
template <>
void Integrator<Quat>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // set the node's name
  node->name = "Integrator";
}

}

