/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The objects in contact
    sorted_pair<BasePtr> objects;  

    /// Coefficient of restitution for contact (default is 0.0)
    double epsilon;

    /// Coefficient of Coulomb friction for contact (default is 0.0)
    double mu_coulomb;

    /// Coefficient of viscous friction for contact (default is 0.0)
    double mu_viscous;

    /// Penalty Method Depth Penalty
    double penalty_Kp;
    
    /// Penalty Method Interpenetration Speed
    double penalty_Kv;

    /// Number of edges in the polygon friction cone, minimum of 4 (default 4)
    unsigned NK;
}; // end class

} // end namespace

#endif  // #ifdef MOBY_CONTACT_PARAMETERS_H_
