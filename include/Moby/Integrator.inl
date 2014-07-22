/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Implements Base::load_from_xml()
void Integrator::load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  std::list<boost::shared_ptr<const XMLTree> > child_nodes;
  std::map<std::string, BasePtr>::const_iterator id_iter;

  // call parent method
  Base::load_from_xml(node, id_map);
}

/// Implements Base::save_to_xml()
void Integrator::save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // set the node's name
  node->name = "Integrator";
}


