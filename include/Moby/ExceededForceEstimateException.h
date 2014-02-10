/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_EXCEEDED_FORCE_ESTIMATE_EXCEPTION_H_
#define _MOBY_EXCEEDED_FORCE_ESTIMATE_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when an integrator discovers the magnitude of a force exceeds its estimate 
class ExceededForceEstimateException : public std::runtime_error
{
  public:
    ExceededForceEstimateException() : std::runtime_error("Force magnitude estimate exceeded") {}
    ExceededForceEstimateException(const char* error) : std::runtime_error(error) {}
}; // end class

} // end namespace

#endif

