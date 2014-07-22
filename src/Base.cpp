/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <sstream>
#include <iostream>
#include <Moby/XMLTree.h>
#include <Moby/Base.h>

using boost::shared_ptr;
using namespace Moby;

namespace Moby {
class BaseData
{
  public:
    BaseData(const Base* b);
    void populate(Base* b);
  
  private:
    boost::shared_ptr<void> _userdata;
}; // end class
} // end namespace

BaseData::BaseData(const Base* b)
{
  _userdata = b->userdata;
}

void BaseData::populate(Base* b)
{
  b->userdata = _userdata;
}

Base::Base() 
{
  std::ostringstream oss;
  oss << this;
  id = oss.str(); 
}

Base::Base(const Base* b) 
{ 
  userdata = b->userdata; 
  std::ostringstream oss;
  oss << this;
  id = oss.str(); 
}

/// Method for loading the data for this object from XML
/**
 * \param node the subtree under which all data necessary to load this
 *        object is stored
 * \param id_map a map from node IDs to read objects
 */
void Base::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the ID, if specified
  XMLAttrib
* id_attrib = node->get_attrib("id");
  if (id_attrib)
    id = id_attrib->get_string_value();

  // verify that the ID does not already exist in the map
  if (id_map.find(id) != id_map.end())
  {
    std::cerr << "Base::load_from_xml() - \"unique\" ID '" << id << "' ";
    std::cerr << "  already exists in the id_map!  Adding anyway...";
    std::cerr << std::endl;
  }

  // add the ID to the ID map
  id_map[id] = shared_from_this();
}

/// Method for saving this object to XML
/**
 * \param node the XML node to which this object should be serialized 
 * \param on output, a list of shared objects which should also be serialized
 */
void Base::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // set the ID for the node
  node->attribs.insert(XMLAttrib("id", id));
}

