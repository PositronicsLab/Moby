/****************************************************************************
 * Copyright 2013 Samuel Zapolsky
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ***************************************************************************/

#ifndef SUSTAINEDCONSTRAINTPROBLEMDATA_H
#define SUSTAINEDCONSTRAINTPROBLEMDATA_H

#include <vector>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/Types.h>

namespace Moby {

struct SustainedUnilateralConstraintProblemData
{
  // setup reasonable defaults
  SustainedUnilateralConstraintProblemData()
  {
    reset();
  }

  // copies contact problem data
  SustainedUnilateralConstraintProblemData& operator=(const SustainedUnilateralConstraintProblemData& q)
  {
    return copy_from(q);
  }

  // copies contact problem data
  SustainedUnilateralConstraintProblemData& copy_from(const SustainedUnilateralConstraintProblemData& q)
  {
    N_K_TOTAL = q.N_K_TOTAL;
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
    ALPHA_X_IDX = q.ALPHA_X_IDX;
    N_VARS = q.N_VARS;

    // copy contact constraints
    contact_constraint_set = q.contact_constraint_set;

    // copy contact accelerations
    Cn_a = q.Cn_a;
    Cs_a = q.Cs_a;
    Ct_a = q.Ct_a;
    L_a = q.L_a;
    Jx_a = q.Jx_a;

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

    // copy the working set
    contact_working_set = q.contact_working_set;

    return *this;
  }

  // resets all contact problem data
  void reset()
  {
    N_K_TOTAL = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_EQNS_IMP = 0;
    N_GC = 0;

    // clear all indices
    N_VARS = 0;
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = L_IDX = ALPHA_X_IDX = 0;

    // clear all vectors
    contact_constraint_set.clear();
    super_bodies.clear();
    constraints.clear();
    contact_constraints.clear();
    limit_constraints.clear();

    // reset all Ravelin::VectorNd sizes
    Cn_a.resize(0);
    Cs_a.resize(0);
    Ct_a.resize(0);
    L_a.resize(0);
    Jx_a.resize(0);
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

    // reset the working set
    contact_working_set.clear();
  }

  // sets up indices 
  void set_indices()
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

  // sets cn, cs, etc. from stacked vectors
  void update_from_stacked(const Ravelin::VectorNd& z)
  {
    cn = z.segment(CN_IDX, CS_IDX);
    l = z.segment(L_IDX, ALPHA_X_IDX);
    alpha_x = z.segment(ALPHA_X_IDX, N_VARS);

    // setup cs/ct -- first determine linearized friction cone forces
    cs.segment(0, N_STICKING) = z.segment(CS_IDX, CT_IDX);
    cs.segment(0, N_STICKING) -= z.segment(NCS_IDX, NCT_IDX);
    ct.segment(0, N_STICKING) = z.segment(CT_IDX, NCS_IDX);
    ct.segment(0, N_STICKING) -= z.segment(NCT_IDX, L_IDX);
  }

  // sets stacked vector from cn, cs, etc.
  Ravelin::VectorNd& to_stacked(Ravelin::VectorNd& z)
  {
    z.set_sub_vec(CN_IDX, cn);
    z.set_sub_vec(CS_IDX, cs.segment(0, N_STICKING));
    z.set_sub_vec(CT_IDX, ct.segment(0, N_STICKING));
    z.set_sub_vec(NCS_IDX, cs.segment(N_STICKING, cs.size()));
    z.set_sub_vec(NCT_IDX, ct.segment(N_STICKING, ct.size()));

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

  // total number of primal variables
  unsigned N_VARS;

  // the total number of linearized friction tangents for contacts
  unsigned N_K_TOTAL;

  // the total number of sticking contacts (out of N_CONTACTS)
  unsigned N_STICKING;

  // the number of contacts
  unsigned N_CONTACTS;

  // the number of limit constraints 
  unsigned N_LIMITS;

  // the total number of constraints
  unsigned N_CONSTRAINTS;

  // the total number of generalized coordinates
  unsigned N_GC;

  // the number of implicit joint constraints (total)
  unsigned N_CONSTRAINT_EQNS_IMP;

  // indication of contacts that the solver is actively considering
  std::vector<bool> contact_working_set;

  // indication if contact is sticking
  std::vector<bool> contact_sticking; // NOTE!:THis value is not calculated yet

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies;

  // the vectors of acceleration-level constraints 
  std::vector<UnilateralConstraint*> constraints, contact_constraints, limit_constraints;

  // the vector indicating which contact constraints are in the contact constraint set (not currently used)
  std::vector<bool> contact_constraint_set;

  // cross-constraint terms
  Ravelin::MatrixNd Cn_iM_CnT, Cn_iM_CsT, Cn_iM_CtT, Cn_iM_LT, Cn_iM_JxT;
  Ravelin::MatrixNd Cs_iM_CnT, Cs_iM_CsT, Cs_iM_CtT, Cs_iM_LT, Cs_iM_JxT;
  Ravelin::MatrixNd Ct_iM_CnT, Ct_iM_CsT, Ct_iM_CtT, Ct_iM_LT, Ct_iM_JxT;
  Ravelin::MatrixNd L_iM_CnT,  L_iM_CsT,  L_iM_CtT,  L_iM_LT, L_iM_JxT;
  Ravelin::MatrixNd Jx_iM_CnT, Jx_iM_CsT, Jx_iM_CtT, Jx_iM_LT, Jx_iM_JxT;

  // vector-based terms
  Ravelin::VectorNd Cn_a, Cs_a, Ct_a, L_a, Jx_a;

  // force magnitudes determined by solve_lcp()
  Ravelin::VectorNd cn, cs, ct, l, alpha_x;
}; // end struct

} // end namespace Moby

#endif // CONTACTPROBLEMDATA_H
