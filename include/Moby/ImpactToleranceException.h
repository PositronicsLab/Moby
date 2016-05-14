/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _IMPACT_TOLERANCE_EXCEPTION_H_
#define _IMPACT_TOLERANCE_EXCEPTION_H_

#include <list>
#include <stdexcept>
#include <Moby/Constraint.h>

namespace Moby {

/// Exception thrown when constraint violation is greater than a desired tolerance 
class ImpactToleranceException : public std::runtime_error
{
  public:
    ImpactToleranceException(const std::list<Constraint*>& impacting_constraints, double max_vio) : std::runtime_error("Constraint velocity is below tolerance after treatment!") { constraints = impacting_constraints; violation = max_vio; }

    virtual ~ImpactToleranceException() throw() { }

  std::list<Constraint*> constraints;
  double violation;  // the amount of constraint violation
}; // end class


} // end namespace

#endif

