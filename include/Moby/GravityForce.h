/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GRAVITY_FORCE_H
#define _GRAVITY_FORCE_H

#include <Ravelin/Vector3d.h>
#include <Moby/RecurrentForce.h>

namespace Moby {
class GravityForce : public RecurrentForce
{
  public:
    GravityForce();
    GravityForce(const GravityForce& source);
    virtual ~GravityForce() {}
    virtual void add_force(boost::shared_ptr<DynamicBody> body);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The gravity vector
    Ravelin::Vector3d gravity;
}; // end class
} // end namespace

#endif

