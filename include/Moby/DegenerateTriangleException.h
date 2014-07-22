/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _DEGENERATE_TRIANGLE_EXCEPTION_H_
#define _DEGENERATE_TRIANGLE_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when trying to perform an operation that requires a normal on a degenerate triangle 
class DegenerateTriangleException : public std::runtime_error
{
  public:
    DegenerateTriangleException() : std::runtime_error("Tried to perform an operation that requires a normal on a degenerate triangle") {}
}; // end class


} // end namespace

#endif

