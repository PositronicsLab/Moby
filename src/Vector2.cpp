/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Matrix2.h>
#include <Moby/MatrixN.h>
#include <Moby/Vector2.h>

using namespace Moby;
using namespace boost;

/// Constructs this vector with the given values
Vector2::Vector2(Real x, Real y)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
Vector2::Vector2(const Real* array)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
}

/// Computes the outer product of two Vector2 objects
Matrix2* Vector2::outer_prod(const Vector2& v1, const Vector2& v2, Matrix2* M)
{
  const unsigned X = 0, Y = 1;
  Matrix2& MM = *M;
  MM(X,X) = v1[X]*v2[X];
  MM(X,Y) = v1[X]*v2[Y];
  MM(Y,X) = v1[Y]*v2[X];
  MM(Y,Y) = v1[Y]*v2[Y];
  return M;
}

