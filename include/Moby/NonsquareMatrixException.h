/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _NONSQUARE_MATRIX_EXCEPTION_H_
#define _NONSQUARE_MATRIX_EXCEPTION_H_

#include <stdexcept>

namespace Moby {

/// Exception thrown when trying to initialize a fixed size vector/matrix with the wrong size
class NonsquareMatrixException : public std::runtime_error
{
  public:
    NonsquareMatrixException() : std::runtime_error("Matrix should be square!") {}
}; // end class


} // end namespace

#endif

