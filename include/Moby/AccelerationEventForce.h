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
class AccelerationEventForce : public RecurrentForce
{
  public:
    virtual ~AccelerationEventForce() {}

    /// Adds the resting contact force, if any, to the body
    virtual void add_force(boost::shared_ptr<DynamicBody> body)
    {
      // look for the resting contact force
      std::map<DynamicBodyPtr, Ravelin::VectorNd>::const_iterator i = accel_event_forces.find(body);
      if (i == accel_event_forces.end())
        return;
      else if (i->second.size() == 0)
        return;
      else
        i->first->add_generalized_force(i->second);
    }

    /// The mapping from bodies to (generalized) acceleration event forces 
    std::map<DynamicBodyPtr, Ravelin::VectorNd> accel_event_forces;
}; // end class
} // end namespace

#endif

