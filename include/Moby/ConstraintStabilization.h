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
#include <Moby/Constraint.h>
#include <Moby/PairwiseDistInfo.h>

namespace Moby {

class ConstraintSimulator;

/// Projected constraint stabilization mechanism
class ConstraintStabilization 
{
  public:
    ConstraintStabilization();
    void stabilize(boost::shared_ptr<ConstraintSimulator> sim); 

    // tolerance to solve unilateral constraints to
    double eps;

    // tolerance to solve bilateral constraints to
    double bilateral_eps;

    // maximum number of iterations for constraint stabilization
    unsigned max_iterations;

    // the inequality tolerance for solving QPs (set to negative value to set programmatically)
    double inequality_tolerance;

  private:
    void get_body_configurations(Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void update_body_configurations(const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    bool update_q(const Ravelin::VectorNd& dq, Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    void add_contact_constraints(std::vector<Constraint>& constraints, CollisionGeometryPtr cg1, CollisionGeometryPtr cg2, boost::shared_ptr<ConstraintSimulator> sim);
    void add_joint_constraints(const std::vector<ControlledBodyPtr>& bodies, std::vector<Constraint>& constraints, boost::shared_ptr<ConstraintSimulator> sim);
    void add_limit_constraints(const std::vector<ControlledBodyPtr>& bodies, std::vector<Constraint>& constraints, boost::shared_ptr<ConstraintSimulator> sim);
    void generate_body_index_map(std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& body_index_map, boost::shared_ptr<ConstraintSimulator> sim);
    void determine_dq(std::vector<Constraint*>& pd, const std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& super_bodies, Ravelin::VectorNd& dqm, const std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, unsigned>& body_index_map);
    static double get_min_pairwise_dist(const std::vector<PairwiseDistInfo>& pdi); 
    static boost::shared_ptr<Ravelin::DynamicBodyd> get_super_body_from_rigid_body(boost::shared_ptr<Ravelin::RigidBodyd> sb);
    static boost::shared_ptr<Ravelin::DynamicBodyd> get_super_body(boost::shared_ptr<Ravelin::DynamicBodyd> sb);
    double ridders_unilateral(double x1, double x2, double fx1, double fx2, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double ridders_bilateral(double x1, double x2, double fx1, double fx2, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double eval_unilateral(double t, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    double eval_bilateral(double t, unsigned i, const Ravelin::VectorNd& dq, const Ravelin::VectorNd& q, boost::shared_ptr<ConstraintSimulator> sim);
    static void save_velocities(boost::shared_ptr<ConstraintSimulator> sim, std::vector<Ravelin::VectorNd>& qd);
    static void restore_velocities(boost::shared_ptr<ConstraintSimulator> sim, const std::vector<Ravelin::VectorNd>& qd);
    static double evaluate_unilateral_constraints(boost::shared_ptr<ConstraintSimulator> sim, std::vector<double>& uC);
    static double evaluate_bilateral_constraints(boost::shared_ptr<ConstraintSimulator> sim, std::vector<double>& C);
    void compute_problem_data(std::vector<std::vector<Constraint*> >& pd_vector, std::vector<std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > >& super_bodies, boost::shared_ptr<ConstraintSimulator> sim);

    // temporary variables
    Ravelin::MatrixNd _H, _A, _M, _MM, _Haug, _Aaug, _Maug;
    Ravelin::VectorNd _c, _b, _q, _qq, _lb, _ub, _lbaug, _ubaug, _caug;

    // the LCP solver
    LCP _lcp;

    // the unilateral constraints
    std::vector<Constraint> constraints;

    /// Geometric pairs that should be checked for unilateral constraints (according to broad phase collision detection)
    std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> > _pairs_to_check;
}
;

} // end namespace

#endif
