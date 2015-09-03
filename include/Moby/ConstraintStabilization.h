/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CONSTRAINT_STABILIZATION_H
#define _CONSTRAINT_STABILIZATION_H

#include <vector>
#include <Ravelin/DynamicBodyd.h>
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

    // tolerance to solve unilateral constraints to
    double eps;

    // tolerance to solve bilateral constraints to
    double bilateral_eps;

  private:
    void get_body_configurations(Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void update_body_configurations(const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void update_q(const Ravelin::VectorNd& dq, Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void compute_problem_data(std::vector<UnilateralConstraintProblemData>& pd, boost::shared_ptr<ConstraintSimulator> sim);
    void add_contact_constraints(std::vector<UnilateralConstraint>& constraints, boost::shared_ptr<Ravelin::RigidBodyd> rb1, boost::shared_ptr<Ravelin::RigidBodyd> rb2, boost::shared_ptr<ConstraintSimulator> sim);
    void add_articulate_limit_constraint(std::vector<UnilateralConstraint>& constraints, boost::shared_ptr<Ravelin::ArticulatedBodyd> ab);
    void generate_body_index_map(std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& body_index_map, boost::shared_ptr<ConstraintSimulator> sim);
    static void set_unilateral_constraint_data(UnilateralConstraintProblemData& pd);
    void determine_dq(UnilateralConstraintProblemData& pd, Ravelin::VectorNd& dqm, const std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& body_index_map);
    void update_from_stacked(const Ravelin::VectorNd& z, UnilateralConstraintProblemData& pd);
    void update_velocities(const UnilateralConstraintProblemData& pd);
    static double get_min_pairwise_dist(const std::vector<PairwiseDistInfo>& pdi); 
    static boost::shared_ptr<Ravelin::DynamicBodyd> get_super_body(boost::shared_ptr<Ravelin::SingleBodyd> sb);
    double ridders_unilateral(double x1, double x2, double fx1, double fx2, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double ridders_bilateral(double x1, double x2, double fx1, double fx2, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double eval_unilateral(double t, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double eval_bilateral(double t, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    static void save_velocities(boost::shared_ptr<ConstraintSimulator> sim, std::vector<Ravelin::VectorNd>& qd);
    static void restore_velocities(boost::shared_ptr<ConstraintSimulator> sim, const std::vector<Ravelin::VectorNd>& qd);
    static void add_contact_to_Jacobian(const UnilateralConstraint& c, SparseJacobian& Cn, const std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& gc_map, unsigned contact_idx);
    static double evaluate_implicit_constraints(boost::shared_ptr<ConstraintSimulator> sim, std::vector<double>& C);

    // the LCP solver
    LCP _lcp;
}
;

} // end namespace

#endif
