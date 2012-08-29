/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _IMPLICIT_STEP_COLLISION_METHOD_H
#define _IMPLICIT_STEP_COLLISION_METHOD_H

#include <Moby/CollisionMethod.h>

namespace Moby {

class Contact;

/// Defines an abstract collision method (i.e., a method for handling colliding, not resting, contact) that is able to step the simulation implicitly
class ImplicitStepCollisionMethod : public CollisionMethod
{
  public:
    ImplicitStepCollisionMethod() {}
    virtual ~ImplicitStepCollisionMethod() {}
  
    /// Implicit stepping methods are always global
    virtual bool is_global_method() const { return true; }

    /// Integrates the equations of motion implicitly for the bodies in the contacts
    virtual void integrate_implicit(const std::list<Contact>& contacts, Real dt) = 0;
}; // end class

} // end namespace

#endif

