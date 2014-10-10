/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _IMPACT_TOLERANCE_EXCEPTION_H_
#define _IMPACT_TOLERANCE_EXCEPTION_H_

#include <list>
#include <stdexcept>
#include <Moby/UnilateralConstraint.h>

namespace Moby {

/// Exception thrown when trying to initialize a fixed size vector/matrix with the wrong size
class ImpactToleranceException : public std::runtime_error
{
  public:
    ImpactToleranceException(const std::list<UnilateralConstraint*>& impacting_constraints) : std::runtime_error("Constraint velocity is below tolerance after treatment!") { constraints = impacting_constraints; }

    virtual ~ImpactToleranceException() throw() { }

  std::list<UnilateralConstraint*> constraints;
}; // end class


} // end namespace

#endif

