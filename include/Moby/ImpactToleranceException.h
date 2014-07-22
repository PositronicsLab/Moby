/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _IMPACT_TOLERANCE_EXCEPTION_H_
#define _IMPACT_TOLERANCE_EXCEPTION_H_

#include <list>
#include <stdexcept>
#include <Moby/Event.h>

namespace Moby {

/// Exception thrown when trying to initialize a fixed size vector/matrix with the wrong size
class ImpactToleranceException : public std::runtime_error
{
  public:
    ImpactToleranceException(const std::list<Event*>& impacting_events) : std::runtime_error("Event velocity is below tolerance after treatment!") { events = impacting_events; }

    virtual ~ImpactToleranceException() throw() { }

  std::list<Event*> events;
}; // end class


} // end namespace

#endif

