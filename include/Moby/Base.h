/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _BASE_H
#define _BASE_H

#include <boost/enable_shared_from_this.hpp>
#include <list>
#include <string>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/Types.h>

namespace Moby {

/// Class from which all Moby classes are derived
class Base : public boost::enable_shared_from_this<Base>
{
  public:
    Base();
    Base(const Base* b);
    virtual ~Base()  {}  
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);  

    /// Static method for cloning a shared pointer
    template <class T>
    static boost::shared_ptr<T> clone(boost::shared_ptr<T> x) { return (!x) ? x : boost::shared_ptr<T>(new T(*x)); }

    /// Any relevant userdata
    boost::shared_ptr<void> userdata;

    /// The unique ID for this object
    std::string id;
};

} // end namespace Moby

#endif
