/****************************************************************************
 * Copyright 2014 Samuel Zapolsky 
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#ifndef _PENALTY_EVENT_HANDLER_H
#define _PENALTY_EVENT_HANDLER_H

#include <list>
#include <vector>
#include <map>
#include <Ravelin/LinAlgd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/UnilateralConstraintProblemData.h>

namespace Moby {

/// Defines the mechanism for handling Penalty constraints
class PenaltyConstraintHandler
{
  public:
    PenaltyConstraintHandler();
    void process_constraints(const std::vector<UnilateralConstraint>& constraints) const;
  private:
    void apply_model(const std::vector<UnilateralConstraint>& constraints) const;
    static double sqr(double x) { return x*x; }
}; // end class
} // end namespace

#endif

