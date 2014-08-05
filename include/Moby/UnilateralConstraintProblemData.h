#ifndef _MOBY_EVENT_PROBLEM_DATA_H
#define _MOBY_EVENT_PROBLEM_DATA_H

#include <boost/foreach.hpp>
#include <vector>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/Types.h>

namespace Moby {

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
    N_K_TOTAL = q.N_K_TOTAL;
    N_LIN_CONE = q.N_LIN_CONE;
    N_TRUE_CONE = q.N_TRUE_CONE;
    N_CONTACTS = q.N_CONTACTS;
    N_CONSTRAINTS = q.N_CONSTRAINTS;
    N_CONSTRAINT_EQNS_IMP = q.N_CONSTRAINT_EQNS_IMP;
    N_ACT_K = q.N_ACT_K;
    N_ACT_CONTACTS = q.N_ACT_CONTACTS;
    N_CONTACT_CONSTRAINTS = q.N_CONTACT_CONSTRAINTS;
    N_GC = q.N_GC;
    kappa = q.kappa;

    // copy indices
    CN_IDX = q.CN_IDX;
    CS_IDX = q.CS_IDX;
    CT_IDX = q.CT_IDX;
    NCS_IDX = q.NCS_IDX;
    NCT_IDX = q.NCT_IDX;
    L_IDX = q.L_IDX;
    ALPHA_X_IDX = q.ALPHA_X_IDX;
    N_VARS = q.N_VARS;  

    // copy contact constraints
    contact_constraint_set = q.contact_constraint_set;

    // copy constraint velocities
    Cn_v = q.Cn_v;
    Cs_v = q.Cs_v;
    Ct_v = q.Ct_v;
    L_v = q.L_v;
    Jx_v = q.Jx_v;

    // the vector of "super" bodies
    super_bodies = q.super_bodies; 

    // the vectors of constraints
    constraints = q.constraints;
    contact_constraints = q.contact_constraints;
    limit_constraints = q.limit_constraints;

    // cross-constraint terms
    Cn_iM_CnT = q.Cn_iM_CnT;
    Cn_iM_CsT = q.Cn_iM_CsT;
    Cn_iM_CtT = q.Cn_iM_CtT;
    Cn_iM_LT = q.Cn_iM_LT;
    Cn_iM_JxT = q.Cn_iM_JxT;
    Cs_iM_CsT = q.Cs_iM_CsT;
    Cs_iM_CtT = q.Cs_iM_CtT;
    Cs_iM_LT = q.Cs_iM_LT;
    Cs_iM_JxT = q.Cs_iM_JxT;
    Ct_iM_CtT = q.Ct_iM_CtT;
    Ct_iM_LT = q.Ct_iM_LT;
    Ct_iM_JxT = q.Ct_iM_JxT;
    L_iM_LT = q.L_iM_LT;
    L_iM_JxT = q.L_iM_JxT;
    Jx_iM_JxT = q.Jx_iM_JxT;

    // copy impulse magnitudes 
    cn = q.cn;
    cs = q.cs;
    ct = q.ct;
    l = q.l;
    alpha_x = q.alpha_x;

    return *this;
  }

  // resets all constraint problem data
  void reset()
  {
    N_K_TOTAL = N_LIN_CONE = N_TRUE_CONE = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_EQNS_IMP = 0;
    N_ACT_K = N_ACT_CONTACTS = N_CONTACT_CONSTRAINTS = 0;
    N_GC = 0;
    kappa = 0.0;

    // clear all indices
    N_VARS = 0;
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = 0;
    L_IDX = 0;
    ALPHA_X_IDX = 0;

    // clear all vectors
    contact_constraint_set.clear();
    super_bodies.clear();
    constraints.clear();
    contact_constraints.clear();
    limit_constraints.clear();

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
    alpha_x.resize(0);

    // reset all MatrixN sizes
    Cn_iM_CnT.resize(0,0);
    Cn_iM_CsT.resize(0,0);
    Cn_iM_CtT.resize(0,0);
    Cn_iM_LT.resize(0,0);
    Cn_iM_JxT.resize(0,0);
    Cs_iM_CsT.resize(0,0);
    Cs_iM_CtT.resize(0,0);
    Cs_iM_LT.resize(0,0);
    Cs_iM_JxT.resize(0,0);
    Ct_iM_CtT.resize(0,0);
    Ct_iM_LT.resize(0,0);
    Ct_iM_JxT.resize(0,0);
    L_iM_LT.resize(0,0);
    L_iM_JxT.resize(0,0);
    Jx_iM_JxT.resize(0,0);
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
    ALPHA_X_IDX = L_IDX + N_LIMITS;
    N_VARS = ALPHA_X_IDX + N_CONSTRAINT_EQNS_IMP;
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
    ALPHA_X_IDX = L_IDX + N_LIMITS;
    N_VARS = ALPHA_X_IDX + N_CONSTRAINT_EQNS_IMP;
  }

  // sets cn, cs, etc. from stacked vector (NQP version)
  void update_from_stacked_nqp(const Ravelin::VectorNd& z)
  {
    cn = z.segment(CN_IDX, CS_IDX);
    cs = z.segment(CS_IDX, CT_IDX);
    ct = z.segment(CT_IDX, L_IDX);
    l = z.segment(L_IDX, ALPHA_X_IDX);
    alpha_x = z.segment(ALPHA_X_IDX, N_VARS);
  }

  // sets cn, cs, etc. from stacked vector (QP version)
  void update_from_stacked_qp(const Ravelin::VectorNd& z)
  {
    cn = z.segment(CN_IDX, CS_IDX);
    l = z.segment(L_IDX, ALPHA_X_IDX);
    alpha_x = z.segment(ALPHA_X_IDX, N_VARS);

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
    z.set_sub_vec(ALPHA_X_IDX, alpha_x);
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

  // starting index of alpha_x in the stacked vector
  unsigned ALPHA_X_IDX;

  // total number of variables
  unsigned N_VARS;

  // the total number of linearized friction tangents for *active* contacts
  unsigned N_ACT_K;

  // the total number of linearized friction tangents for contact constraints
  unsigned N_K_TOTAL;

  // the number of contacts with linearized friction cones
  unsigned N_LIN_CONE;

  // the number of contacts with true friction cones
  unsigned N_TRUE_CONE;

  // the number of contacts (total)
  unsigned N_CONTACTS;

  // the number of contacts (active)
  unsigned N_ACT_CONTACTS;

  // the number of contact constraints in the optimization problem
  unsigned N_CONTACT_CONSTRAINTS;

  // the number of limits
  unsigned N_LIMITS;

  // the total number of constraints
  unsigned N_CONSTRAINTS;

  // the total number of generalized coordinates
  unsigned N_GC;

  // the velocity in the normal direction for frictionless contacts
  double kappa;

  // the number of implicit joint constraint equations (total)
  unsigned N_CONSTRAINT_EQNS_IMP;

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies; 

  // the vectors of constraints
  std::vector<UnilateralConstraint*> constraints, contact_constraints, limit_constraints;

  // the vector indicating which contact constraints are in the linear constraint set 
  std::vector<bool> contact_constraint_set;

  // cross-constraint terms
  Ravelin::MatrixNd Cn_iM_CnT, Cn_iM_CsT, Cn_iM_CtT, Cn_iM_LT, Cn_iM_JxT;
  Ravelin::MatrixNd            Cs_iM_CsT, Cs_iM_CtT, Cs_iM_LT, Cs_iM_JxT;
  Ravelin::MatrixNd                       Ct_iM_CtT, Ct_iM_LT, Ct_iM_JxT;
  Ravelin::MatrixNd                                   L_iM_LT, L_iM_JxT;
  Ravelin::MatrixNd                                            Jx_iM_JxT;

  // vector-based terms
  Ravelin::VectorNd Cn_v, Cs_v, Ct_v, L_v, Jx_v;

  // impulse magnitudes determined by solve_qp()
  Ravelin::VectorNd cn, cs, ct, l, alpha_x;
}; // end struct

} // end namespace Moby

#endif

