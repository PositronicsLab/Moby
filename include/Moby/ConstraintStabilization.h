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

class ConstraintSimulator;

class ConstraintStabilization 
{
  public:
    ConstraintStabilization();
    void stabilize(boost::shared_ptr<ConstraintSimulator> sim); 

    // tolerance to solve constraints to
    double eps;

  private:
    void get_body_configurations(Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void update_body_configurations(const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void update_q(const Ravelin::VectorNd& dq, Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void compute_problem_data(std::vector<UnilateralConstraintProblemData>& pd, boost::shared_ptr<ConstraintSimulator> sim);
    void add_contact_constraints(std::vector<UnilateralConstraint>& constraints, RigidBodyPtr rb1, RigidBodyPtr rb2, boost::shared_ptr<ConstraintSimulator> sim);
    void add_articulate_limit_constraint(std::vector<UnilateralConstraint>& constraints, ArticulatedBodyPtr ab);
    void generate_body_index_map(std::map<DynamicBodyPtr, unsigned>& body_index_map, boost::shared_ptr<ConstraintSimulator> sim);
    static void set_unilateral_constraint_data(UnilateralConstraintProblemData& pd);
    void determine_dq(UnilateralConstraintProblemData& pd, Ravelin::VectorNd& dqm, const std::map<DynamicBodyPtr, unsigned>& body_index_map);
    void update_from_stacked(const Ravelin::VectorNd& z, UnilateralConstraintProblemData& pd);
    void update_velocities(const UnilateralConstraintProblemData& pd);
    static double get_min_pairwise_dist(const std::vector<PairwiseDistInfo>& pdi); 
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);
    double ridders(double x1, double x2, double fx1, double fx2, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double eval(double t, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    static void save_velocities(boost::shared_ptr<ConstraintSimulator> sim, std::vector<Ravelin::VectorNd>& qd);
    static void restore_velocities(boost::shared_ptr<ConstraintSimulator> sim, const std::vector<Ravelin::VectorNd>& qd);

    // the LCP solver
    LCP _lcp;
}
;

} // end namespace

#endif
