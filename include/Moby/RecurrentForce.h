/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RECURRENT_FORCE
#define _RECURRENT_FORCE

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/Base.h>

namespace Moby {

class DynamicBody;
  
/// Used for applying forces to the simulation on every time step 
/**
 * Recurrent forces may be constant (like gravity) or vary by the velocity
 * of the body (like wind resistance).  The recurrent force can be either an
 * actual force, or a torque, or both.
 */
class RecurrentForce : public virtual Base
{
  public:
    virtual ~RecurrentForce() { }
    
    /// Abstract method for applying this force/torque to a body
    virtual void add_force(DynamicBodyPtr body) = 0;  
}; // end class

} // end namespace Moby

#endif
