#ifndef _MOBY_EVENT_PROBLEM_DATA_H
#define _MOBY_EVENT_PROBLEM_DATA_H

#include <boost/foreach.hpp>
#include <vector>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/sorted_pair>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/SparseJacobian.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/Types.h>

namespace Moby {

class Simulator;

struct UnilateralConstraintProblemData
{
  // setup reasonable defaults
  UnilateralConstraintProblemData()
  {
    reset();
  }

  // copies constraint problem data
  UnilateralConstraintProblemData& operator=(const UnilateralConstraintProblemData& q)
  {
    return copy_from(q);
  }

  // copies constraint problem data
  UnilateralConstraintProblemData& copy_from(const UnilateralConstraintProblemData& q)
  {
    simulator = q.simulator;
    N_K_TOTAL = q.N_K_TOTAL;
    N_LIN_CONE = q.N_LIN_CONE;
    N_TRUE_CONE = q.N_TRUE_CONE;
    N_CONTACTS = q.N_CONTACTS;
    N_CONSTRAINTS = q.N_CONSTRAINTS;
    N_CONSTRAINT_EQNS_IMP = q.N_CONSTRAINT_EQNS_IMP;
    N_GC = q.N_GC;

    // copy indices
    CN_IDX = q.CN_IDX;
    CS_IDX = q.CS_IDX;
    CT_IDX = q.CT_IDX;
    NCS_IDX = q.NCS_IDX;
    NCT_IDX = q.NCT_IDX;
    L_IDX = q.L_IDX;
    N_VARS = q.N_VARS;  

    // copy constraint velocities
    Cn_v = q.Cn_v;
    Cs_v = q.Cs_v;
    Ct_v = q.Ct_v;
    L_v = q.L_v;
    Jx_v = q.Jx_v;

    // copy the signed distances
    signed_distances = q.signed_distances;

    // the vector of "super" bodies
    super_bodies = q.super_bodies; 

    // the vectors of constraints
    constraints = q.constraints;
    contact_constraints = q.contact_constraints;
    limit_constraints = q.limit_constraints;

    // cross-constraint terms
    Cn_X_CnT = q.Cn_X_CnT;
    Cn_X_CsT = q.Cn_X_CsT;
    Cn_X_CtT = q.Cn_X_CtT;
    Cn_X_LT = q.Cn_X_LT;
    Cn_X_JxT = q.Cn_X_JxT;
    Cs_X_CsT = q.Cs_X_CsT;
    Cs_X_CtT = q.Cs_X_CtT;
    Cs_X_LT = q.Cs_X_LT;
    Cs_X_JxT = q.Cs_X_JxT;
    Ct_X_CtT = q.Ct_X_CtT;
    Ct_X_LT = q.Ct_X_LT;
    Ct_X_JxT = q.Ct_X_JxT;
    L_X_LT = q.L_X_LT;
    L_X_JxT = q.L_X_JxT;
    Jx_X_JxT = q.Jx_X_JxT;

    // copy Cdot Jacobians
    Cdot_iM_CnT = q.Cdot_iM_CnT; 
    Cdot_iM_CsT = q.Cdot_iM_CsT; 
    Cdot_iM_CtT = q.Cdot_iM_CtT; 
    Cdot_iM_LT = q.Cdot_iM_LT; 

    // copy implicit constraints, Jacobian, and related terms
    island_ijoints = q.island_ijoints;
    J = q.J;
    Jfull = q.Jfull;
    Jx_iM_JxT = q.Jx_iM_JxT;
    iM_JxT = q.iM_JxT;
    active = q.active;

    // copy limit indices
    limit_indices = q.limit_indices;

    // copy impulse magnitudes 
    cn = q.cn;
    cs = q.cs;
    ct = q.ct;
    l = q.l;
    lambda = q.lambda;

    return *this;
  }

  // resets all constraint problem data
  void reset()
  {
    simulator.reset();
    N_K_TOTAL = N_LIN_CONE = N_TRUE_CONE = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_EQNS_IMP = 0;
    N_GC = 0;

    // clear all indices
    N_VARS = 0;
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = 0;
    L_IDX = 0;

    // clear all vectors
    super_bodies.clear();
    constraints.clear();
    contact_constraints.clear();
    limit_constraints.clear();
    limit_indices.clear();

    // clear signed distances
    signed_distances.clear();

    // clear implicit constraint related stuff
    island_ijoints.clear();
    J.blocks.clear();
    Jfull.blocks.clear();
    Jx_iM_JxT.resize(0,0);
    iM_JxT.resize(0,0);
    active.clear();

    // reset all Ravelin::VectorNd sizes
    Cn_v.resize(0);
    Cs_v.resize(0);
    Ct_v.resize(0);
    L_v.resize(0);
    Jx_v.resize(0);
    cn.resize(0);
    cs.resize(0);
    ct.resize(0);
    l.resize(0);
    lambda.resize(0);

    // reset all MatrixN sizes
    Cn_X_CnT.resize(0,0);
    Cn_X_CsT.resize(0,0);
    Cn_X_CtT.resize(0,0);
    Cn_X_LT.resize(0,0);
    Cn_X_JxT.resize(0,0);
    Cs_X_CsT.resize(0,0);
    Cs_X_CtT.resize(0,0);
    Cs_X_LT.resize(0,0);
    Cs_X_JxT.resize(0,0);
    Ct_X_CtT.resize(0,0);
    Ct_X_LT.resize(0,0);
    Ct_X_JxT.resize(0,0);
    L_X_LT.resize(0,0);
    L_X_JxT.resize(0,0);
    Jx_X_JxT.resize(0,0);
    Cdot_iM_CnT.resize(0,0);
    Cdot_iM_CsT.resize(0,0);
    Cdot_iM_CtT.resize(0,0);
    Cdot_iM_LT.resize(0,0);
  }

  // sets up indices for a QP
  void set_qp_indices()
  {
    CN_IDX = 0;
    CS_IDX = CN_IDX + N_CONTACTS;
    CT_IDX = CS_IDX + N_CONTACTS;
    NCS_IDX = CT_IDX + N_CONTACTS;
    NCT_IDX = NCS_IDX + N_CONTACTS;
    L_IDX = NCT_IDX + N_CONTACTS;
    N_VARS = L_IDX + N_LIMITS;
  }

  // sets up indices for a nonlinear QP
  void set_nqp_indices()
  {
    CN_IDX = 0;
    CS_IDX = CN_IDX + N_CONTACTS;
    CT_IDX = CS_IDX + N_CONTACTS;
    NCS_IDX = CT_IDX + N_CONTACTS;
    NCT_IDX = NCS_IDX + 0;
    L_IDX = NCT_IDX + 0;
    N_VARS = L_IDX + N_LIMITS;
  }

  // sets cn, cs, etc. from stacked vector (NQP version)
  void update_from_stacked_nqp(const Ravelin::VectorNd& z)
  {
    cn = z.segment(CN_IDX, CS_IDX);
    cs = z.segment(CS_IDX, CT_IDX);
    ct = z.segment(CT_IDX, L_IDX);
    l = z.segment(L_IDX, N_VARS);
  }

  // sets cn, cs, etc. from stacked vector (QP version)
  void update_from_stacked_qp(const Ravelin::VectorNd& z)
  {
    cn = z.segment(CN_IDX, CS_IDX);
    l = z.segment(L_IDX, N_VARS);

    // setup cs/ct -- first determine linearized friction cone forces
    cs.segment(0, N_LIN_CONE) = z.segment(CS_IDX, CT_IDX);
    cs.segment(0, N_LIN_CONE) -= z.segment(NCS_IDX, NCT_IDX);
    ct.segment(0, N_LIN_CONE) = z.segment(CT_IDX, NCS_IDX);
    ct.segment(0, N_LIN_CONE) -= z.segment(NCT_IDX, L_IDX);
  }

  // partitions constraint vectors into contact and limit constraints
  void partition_constraints()
  {
    const unsigned UINF = std::numeric_limits<unsigned>::max();
    contact_constraints.clear();
    limit_constraints.clear();

    BOOST_FOREACH(UnilateralConstraint* e, constraints)
    {
      if (e->constraint_type == UnilateralConstraint::eContact)
        contact_constraints.push_back(e);
      else
      {
        assert(e->constraint_type == UnilateralConstraint::eLimit);
        limit_constraints.push_back(e);
      }
    }
    
    // now, sort the contact constraints such that constraints that use a true friction
    // cone are at the end
    if(!contact_constraints.empty()){
      for (unsigned i=0, j=contact_constraints.size()-1; i< j; )
      {
        if (contact_constraints[i]->contact_NK == UINF)
        {
          std::swap(contact_constraints[i], contact_constraints[j]);
          j--;
        } 
        else
          i++;
      }
    }
  }

  // sets stacked vector from cn, cs, etc.
  Ravelin::VectorNd& to_stacked(Ravelin::VectorNd& z)
  {
    z.set_sub_vec(CN_IDX, cn);
    z.set_sub_vec(CS_IDX, cs.segment(0, N_LIN_CONE));
    z.set_sub_vec(CT_IDX, ct.segment(0, N_LIN_CONE));
    z.set_sub_vec(NCS_IDX, cs.segment(N_LIN_CONE, cs.size()));
    z.set_sub_vec(NCT_IDX, ct.segment(N_LIN_CONE, ct.size()));
    for (unsigned i=CS_IDX, j=NCS_IDX; i< CT_IDX; i++, j++)
    {
      if (z[i] < 0.0)
      {
        z[j] = -z[i];
        z[i] = 0.0;
      }
      else
        z[j] = 0.0;
    }
    for (unsigned i=CT_IDX, j=NCT_IDX; i< NCS_IDX; i++, j++)
    {
      if (z[i] < 0.0)
      {
        z[j] = -z[i];
        z[i] = 0.0;
      }
      else
        z[j] = 0.0;
    }
    z.set_sub_vec(L_IDX, l);
    return z;
  }

  // starting index of cn in the stacked vector
  unsigned CN_IDX;

  // starting index of cs in the stacked vector
  unsigned CS_IDX;

  // starting index of ct in the stacked vector
  unsigned CT_IDX;

  // starting index of -cs in the stacked vector
  unsigned NCS_IDX;

  // starting index of -ct in the stacked vector
  unsigned NCT_IDX;

  // starting index of l in the stacked vector
  unsigned L_IDX;

  // total number of variables
  unsigned N_VARS;

  // the total number of linearized friction tangents for contact constraints
  unsigned N_K_TOTAL;

  // the number of contacts with linearized friction cones
  unsigned N_LIN_CONE;

  // the number of contacts with true friction cones
  unsigned N_TRUE_CONE;

  // the number of contacts (total)
  unsigned N_CONTACTS;

  // the number of limits
  unsigned N_LIMITS;

  // the total number of constraints
  unsigned N_CONSTRAINTS;

  // the total number of generalized coordinates
  unsigned N_GC;

  // the number of implicit joint constraint equations (total)
  unsigned N_CONSTRAINT_EQNS_IMP;

  // pairwise distances between rigid bodies
  std::vector<PairwiseDistInfo> signed_distances;

  // the vector of "super" bodies
  std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > super_bodies; 

  // the vectors of constraints
  std::vector<UnilateralConstraint*> constraints, contact_constraints, limit_constraints;

  // the vector indicating which contact constraints are in the linear constraint set 
  // cross-constraint terms
  Ravelin::MatrixNd Cn_X_CnT, Cn_X_CsT, Cn_X_CtT, Cn_X_LT,    Cn_X_JxT;
  Ravelin::MatrixNd            Cs_X_CsT, Cs_X_CtT, Cs_X_LT,   Cs_X_JxT;
  Ravelin::MatrixNd                       Ct_X_CtT, Ct_X_LT,  Ct_X_JxT;
  Ravelin::MatrixNd                                   L_X_LT, L_X_JxT;
  Ravelin::MatrixNd                                           Jx_X_JxT;

  // bilateral constraint inertia terms 
  Ravelin::MatrixNd iM_JxT, Jx_iM_JxT;

  // X times Jacobian transposes
  Ravelin::MatrixNd X_CnT, X_CsT, X_CtT, X_LT, X_JxT;

  // Cdot Jacobians
  Ravelin::MatrixNd Cdot_iM_CnT, Cdot_iM_CsT, Cdot_iM_CtT, Cdot_iM_LT;

  // vector-based terms
  Ravelin::VectorNd Cn_v, Cs_v, Ct_v, L_v, Jx_v;

  // impulse magnitudes determined by solve_qp()
  Ravelin::VectorNd cn, cs, ct, l;

  // bilateral constraint impulses
  Ravelin::VectorNd lambda;

  // active implicit bilateral constraints
  std::vector<bool> active;

  // implicit bilateral constraints in this island
  std::vector<JointPtr> island_ijoints;

  // the implicit constraint Jacobian (full)
  SparseJacobian Jfull;

  // the implicit constraint Jacobian (reduced)
  SparseJacobian J;

  // indices in the generalized coordinates of the various limits
  std::vector<unsigned> limit_indices;

  // a pointer to the simulator
  boost::shared_ptr<Simulator> simulator;
}; // end struct

} // end namespace Moby

#endif

