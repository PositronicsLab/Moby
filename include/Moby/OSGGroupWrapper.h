/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _OSG_GROUP_WRAPPER_H_
#define _OSG_GROUP_WRAPPER_H_

#include <Moby/Base.h>

namespace osg
{
  class Node;
  class Group;
}

namespace Moby {

/// A wrapper for OpenInventor OSGGroup class, supporting serialization
class OSGGroupWrapper : public virtual Base
{
  public:
    OSGGroupWrapper();
    OSGGroupWrapper(osg::Node* n); 
    OSGGroupWrapper(const std::string& filename);
    ~OSGGroupWrapper();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    osg::Group* get_group() { return _group; }

  private:
    osg::Group* _group;
}; // end class

} // end namespace

#endif

