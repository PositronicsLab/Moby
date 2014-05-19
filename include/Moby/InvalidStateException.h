/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_INVALID_STATE_EXCEPTION_H_
#define _MOBY_INVALID_STATE_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when an integrator tries to evaluate a derivative at an invalid state 
class InvalidStateException : public std::runtime_error
{
  public:
    InvalidStateException() : std::runtime_error("Integrator tries to evaluate derivative at invalid state") {}
    InvalidStateException(const char* error) : std::runtime_error(error) {}
}; // end class

} // end namespace

#endif

