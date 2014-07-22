/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_INVALID_VELOCITY_EXCEPTION_H_
#define _MOBY_INVALID_VELOCITY_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when an integrator tries to evaluate a derivative at an invalid velocity 
class InvalidVelocityException : public std::runtime_error
{
  public:
    InvalidVelocityException(double t) : std::runtime_error("Integrator tries to evaluate derivative at invalid velocity") { evaluation_time = t; }
    double evaluation_time; 
}; // end class

} // end namespace

#endif

