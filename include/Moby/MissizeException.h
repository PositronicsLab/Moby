/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MISSIZE_EXCEPTION_H_
#define _MISSIZE_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when trying to initialize a fixed size vector/matrix with the wrong size
class MissizeException : public std::runtime_error
{
  public:
    MissizeException() : std::runtime_error("Matrix / vector size mismatch") {}
}; // end class


} // end namespace

#endif

