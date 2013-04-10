/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _DAMPING_FORCE_H
#define _DAMPING_FORCE_H

#include <Moby/RecurrentForce.h>

namespace Moby {

/// Defines a force that damps rigid bodies on any given step
class DampingForce : public RecurrentForce
{
  public:
    DampingForce() {}
    DampingForce(const DampingForce& source);
    virtual ~DampingForce() {}
    virtual void add_force(boost::shared_ptr<DynamicBody> body);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// The mapping from bodies to linear damping constants
    std::map<DynamicBodyPtr, double> kl;

    /// The mapping from bodies to angular damping constants
    std::map<DynamicBodyPtr, double> ka;

    /// The mapping from bodies to linear squared damping constants
    std::map<DynamicBodyPtr, double> klsq;

    /// The mapping from bodies to angular squared damping constants
    std::map<DynamicBodyPtr, double> kasq;
}; // end class
} // end namespace

#endif

