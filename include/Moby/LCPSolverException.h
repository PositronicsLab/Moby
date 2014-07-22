/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_LCP_SOLVER_EXCEPTION_H_
#define _MOBY_LCP_SOLVER_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when LCP solver can't solve 
class LCPSolverException : public std::runtime_error
{
  public:
    LCPSolverException() : std::runtime_error("Unable to solve event LCP/QP!") {}
}; // end class


} // end namespace

#endif

