/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _SIGNED_DIST_DOT_H
#define _SIGNED_DIST_DOT_H

#include <list>
#include <vector>
#include <map>
#include <Ravelin/LinAlgd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/UnilateralConstraintProblemData.h>

namespace Moby {

class ConstraintSimulator;

/// Defines the mechanism for handling impact constraints
class SignedDistDot
{
  friend class ConstraintSimulator;
  friend class ConstraintStabilization;

  public:
    static void compute_signed_dist_dot_Jacobians(UnilateralConstraintProblemData& q, Ravelin::MatrixNd& CnT, Ravelin::MatrixNd& CsT, Ravelin::MatrixNd& CtT, Ravelin::MatrixNd& LT);

  private:
    static double calc_signed_dist(boost::shared_ptr<Ravelin::SingleBodyd> sb1, boost::shared_ptr<Ravelin::SingleBodyd> sb2);
    static void restore_coords_and_velocities(const std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& isect, std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, Ravelin::VectorNd>& gc_map, std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, Ravelin::VectorNd>& gv_map);
    static void integrate_positions(const std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& isect, double dt);
    static void apply_impulse(const UnilateralConstraint& contact_constraint, const Ravelin::Vector3d& dir);
    static void apply_impulse(const UnilateralConstraint& contact_constraint);

    // a pointer to the simulator
    boost::shared_ptr<ConstraintSimulator> _simulator;
}; // end class
} // end namespace

#endif

