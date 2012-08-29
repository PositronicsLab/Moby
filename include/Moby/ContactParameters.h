/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_CONTACT_PARAMETERS_H_
#define _MOBY_CONTACT_PARAMETERS_H_

#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/sorted_pair>

namespace Moby {

class ContactParameters : public Base
{
  public:
    ContactParameters();
    ContactParameters(BasePtr o1, BasePtr o2);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// The objects in contact
    sorted_pair<BasePtr> objects;  

    /// Coefficient of restitution for contact (default is 0.0)
    Real epsilon;

    /// Coefficient of Coulomb friction for contact (default is 0.0)
    Real mu_coulomb;

    /// Coefficient of viscous friction for contact (default is 0.0)
    Real mu_viscous;

    /// Number of edges in the polygon friction cone, minimum of 2 (default 2)
    unsigned NK;
}; // end class

} // end namespace

#endif  // #ifdef MOBY_CONTACT_PARAMETERS_H_
