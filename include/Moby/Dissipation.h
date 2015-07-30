/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _DISSIPATION_H
#define _DISSIPATION_H

#include <Moby/RecurrentForce.h>

namespace Moby {
class Dissipation : public Base 
{
  public:
    virtual ~Dissipation() {}
    virtual void apply(const std::vector<DynamicBodyPtr>& bodies) = 0;
}; // end class
} // end namespace

#endif

