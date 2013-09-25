/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RESTING_CONTACT_FORCE_H
#define _RESTING_CONTACT_FORCE_H

#include <Ravelin/Vector3d.h>
#include <Moby/RecurrentForce.h>

namespace Moby {
class RestingContactForce : public RecurrentForce
{
  public:
    virtual ~RestingContactForce() {}

    /// Adds the resting contact force, if any, to the body
    virtual void add_force(boost::shared_ptr<DynamicBody> body)
    {
      // look for the resting contact force
      std::map<DynamicBodyPtr, Ravelin::VectorNd>::const_iterator i = resting_contact_forces.find(body);
      if (i == resting_contact_forces.end())
        return;
      else if (i->second.size() == 0)
        return;
      else
        i->first->add_generalized_force(i->second);
    }

    /// The mapping from bodies to (generalized) resting contact forces 
    std::map<DynamicBodyPtr, Ravelin::VectorNd> resting_contact_forces;
}; // end class
} // end namespace

#endif

