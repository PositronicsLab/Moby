/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RNE_ALGO_H
#define _RNE_ALGO_H

#include <Moby/RCArticulatedBodyInvDynAlgo.h>

namespace Moby {

/// Implementation of the Recursive Newton-Euler algorithm for inverse dynamics
/**
 * Algorithm taken from Featherstone, 1987.
 */ 
class RNEAlgorithm : public RCArticulatedBodyInvDynAlgo
{
  public:
    std::map<JointPtr, Ravelin::VectorNd> calc_inv_dyn(RCArticulatedBodyPtr body, const std::map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data);
    void calc_constraint_forces(RCArticulatedBodyPtr body);

  private:
    std::map<JointPtr, Ravelin::VectorNd> calc_inv_dyn_fixed_base(RCArticulatedBodyPtr body, const std::map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const;
    std::map<JointPtr, Ravelin::VectorNd> calc_inv_dyn_floating_base(RCArticulatedBodyPtr body, const std::map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const;
};
} // end namespace

#endif

