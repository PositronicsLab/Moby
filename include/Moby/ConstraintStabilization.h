/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CONSTRAINT_STABILIZATION_H
#define _CONSTRAINT_STABILIZATION_H

#include <vector>
#include <Moby/DynamicBody.h>
#include <Moby/LCP.h>
#include <Moby/UnilateralConstraintProblemData.h>
#include <Moby/PairwiseDistInfo.h>

namespace Moby {

class EventDrivenSimulator;

class ConstraintStabilization 
{
  public:
    ConstraintStabilization(boost::shared_ptr<EventDrivenSimulator> sim);
    void stabilize(); 

    // tolerance to solve constraints to
    double eps;

  private:
    void get_body_configurations(Ravelin::VectorNd& q);
    void update_body_configurations(const Ravelin::VectorNd& q);
    void update_q(const Ravelin::VectorNd& dq, Ravelin::VectorNd& q);
    void compute_problem_data(std::vector<UnilateralConstraintProblemData>& pd);
    static void set_unilateral_constraint_data(UnilateralConstraintProblemData& pd);
    void determine_dq(const UnilateralConstraintProblemData& pd, Ravelin::VectorNd& dq);
    static double get_min_pairwise_dist(const std::vector<PairwiseDistInfo>& pdi); 
    static double compute_s(const std::vector<PairwiseDistInfo>& pdi);
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);

    boost::shared_ptr<EventDrivenSimulator> _sim;

    // the LCP solver
    LCP _lcp;
}
;

} // end namespace

#endif
