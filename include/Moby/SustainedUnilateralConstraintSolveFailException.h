/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_RESTING_EVENT_FAIL_EXCEPTION_H_
#define _MOBY_RESTING_EVENT_FAIL_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when general numerical error occurs 
class SustainedUnilateralConstraintSolveFailException : public std::runtime_error
{
  public:
    SustainedUnilateralConstraintSolveFailException() : std::runtime_error("SustainedUnilateralConstraintSolveFailException error") {}
    SustainedUnilateralConstraintSolveFailException(const char* error) : std::runtime_error(error) {}
}; // end class


} // end namespace

#endif

