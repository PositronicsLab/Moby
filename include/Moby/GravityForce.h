/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _GRAVITY_FORCE_H
#define _GRAVITY_FORCE_H

#include <Moby/Vector3.h>
#include <Moby/RecurrentForce.h>

namespace Moby {
class GravityForce : public RecurrentForce
{
  public:
    GravityForce();
    GravityForce(const GravityForce& source);
    virtual ~GravityForce() {}
    virtual void add_force(boost::shared_ptr<DynamicBody> body);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// The gravity vector
    Vector3 gravity;
}; // end class
} // end namespace

#endif

