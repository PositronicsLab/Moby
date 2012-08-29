/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _NONCONVEXITY_EXCEPTION_H_
#define _NONCONVEXITY_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when a non-convex object is encountered and a convex object is expected 
class NonconvexityException : public std::runtime_error
{
  public:
    NonconvexityException(const std::string& error) : std::runtime_error(error) {}
}; // end class


} // end namespace

#endif

