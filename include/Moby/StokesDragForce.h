/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _STOKES_DRAG_FORCE_H
#define _STOKES_DRAG_FORCE_H

#include <Moby/RecurrentForce.h>

namespace Moby {
class StokesDragForce : public RecurrentForce
{
  public:
    StokesDragForce();
    StokesDragForce(const StokesDragForce& source);
    virtual ~StokesDragForce() {}
    virtual void add_force(boost::shared_ptr<DynamicBody> body);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The drag coefficient 
    double b, b_ang;
}; // end class
} // end namespace

#endif

