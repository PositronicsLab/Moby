/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_MEASURE_DATA_H
#define _MOBY_MEASURE_DATA_H

#include <Moby/Vector3.h>
#include <Moby/Types.h>

namespace Moby {

/// Encapsulates data for measurements used by global contact methods
/**
 * This class is used to encapsulate data for measurements used for global
 * contact methods.
 */
class MeasureData
{
  public:
    MeasureData() { point = NULL; }

    /// The point at which the measurement will be taken
    const Vector3* point;

    /// The unit-vector direction of the measurement
    Vector3 dir;

    /// The first body used in the measurement
    SingleBodyPtr sb1;

    /// The second body used in the measurement
    SingleBodyPtr sb2;

    /// A pointer to the relevant contact
    Contact* contact;
}; // end class

} // end namespace

#endif

