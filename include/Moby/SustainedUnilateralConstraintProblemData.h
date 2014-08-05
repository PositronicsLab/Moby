/****************************************************************************
 * Copyright 2013 Samuel Zapolsky
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ***************************************************************************/

#ifndef CONTACTPROBLEMDATA_H
#define CONTACTPROBLEMDATA_H

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

    // copy indices
    CN_IDX = q.CN_IDX;
    CS_IDX = q.CS_IDX;
    CT_IDX = q.CT_IDX;
    NCS_IDX = q.NCS_IDX;
    NCT_IDX = q.NCT_IDX;
    CS_U_IDX = q.CS_U_IDX;
    CT_U_IDX = q.CT_U_IDX;

    // copy contact accelerations
    Cn_a = q.Cn_a;
    Cs_a = q.Cs_a;
    Ct_a = q.Ct_a;

    // the vector of "super" bodies
    super_bodies = q.super_bodies;

    // the vectors of contacts 
    constraints = q.constraints;

    // cross-contact terms
    Cn_iM_CnT = q.Cn_iM_CnT;
    Cn_iM_CsT = q.Cn_iM_CsT;
    Cn_iM_CtT = q.Cn_iM_CtT;
    Cs_iM_CnT = q.Cs_iM_CnT;
    Cs_iM_CsT = q.Cs_iM_CsT;
    Cs_iM_CtT = q.Cs_iM_CtT;
    Ct_iM_CnT = q.Ct_iM_CnT;
    Ct_iM_CsT = q.Ct_iM_CsT;
    Ct_iM_CtT = q.Ct_iM_CtT;

    // copy force magnitudes
    cn = q.cn;
    cs = q.cs;
    ct = q.ct;

    // copy the working set
    contact_working_set = q.contact_working_set;

    return *this;
  }

  // resets all contact problem data
  void reset()
  {
    N_K_TOTAL = N_CONTACTS = 0;

    // clear all indices
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = CS_U_IDX = CT_U_IDX = 0;

    // clear all vectors
    super_bodies.clear();
    constraints.clear();

    // reset all Ravelin::VectorNd sizes
    cn.resize(0);
    cs.resize(0);
    ct.resize(0);
    // reset all MatrixN sizes
    Cn_iM_CnT.resize(0,0);
    Cn_iM_CsT.resize(0,0);
    Cn_iM_CtT.resize(0,0);
    Cs_iM_CnT.resize(0,0);
    Cs_iM_CsT.resize(0,0);
    Cs_iM_CtT.resize(0,0);
    Ct_iM_CnT.resize(0,0);
    Ct_iM_CsT.resize(0,0);
    Ct_iM_CtT.resize(0,0);
    // reset the working set
    contact_working_set.clear();
  }

  // sets cn, cs, etc. from stacked vectors
  void update_from_stacked(const Ravelin::VectorNd& z)
  {
    cn += z.segment(CN_IDX, CS_IDX);

    // setup cs/ct -- first determine linearized friction cone forces
    cs.segment(0, N_STICKING) = z.segment(CS_IDX, CT_IDX);
    cs.segment(0, N_STICKING) -= z.segment(NCS_IDX, NCT_IDX);
    ct.segment(0, N_STICKING) = z.segment(CT_IDX, NCS_IDX);
    ct.segment(0, N_STICKING) -= z.segment(NCT_IDX, CS_U_IDX);
  }

  // sets stacked vector from cn, cs, etc.
  Ravelin::VectorNd& to_stacked(Ravelin::VectorNd& z)
  {
    z.set_sub_vec(CN_IDX, cn);
    z.set_sub_vec(CS_IDX, cs.segment(0, N_STICKING));
    z.set_sub_vec(CT_IDX, ct.segment(0, N_STICKING));
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

  // starting index of cs (unbounded) in the stacked vector
  unsigned CS_U_IDX;

  // starting index of ct (unbounded) in the stacked vector
  unsigned CT_U_IDX;

  // the total number of linearized friction tangents for contacts
  unsigned N_K_TOTAL;

  // the total number of sticking contacts (out of N_CONTACTS)
  unsigned N_STICKING;

  // the number of contacts
  unsigned N_CONTACTS;

  // indication of contacts that the solver is actively considering
  std::vector<bool> contact_working_set;

  // indication if contact is sticking
  std::vector<bool> contact_sticking; // NOTE!:THis value is not calculated yet

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies;

  // the vectors of acceleration-level constraints 
  std::vector<UnilateralConstraint*> constraints;

  // cross-contact terms
  Ravelin::MatrixNd Cn_iM_CnT, Cn_iM_CsT, Cn_iM_CtT, Cn_iM_muS_CqT;
  Ravelin::MatrixNd Cs_iM_CnT, Cs_iM_CsT, Cs_iM_CtT, Cs_iM_muS_CqT;
  Ravelin::MatrixNd Ct_iM_CnT, Ct_iM_CsT, Ct_iM_CtT, Ct_iM_muS_CqT;

  // vector-based terms
  Ravelin::VectorNd Cn_a, Cs_a, Ct_a;

  // EMD: why use this instead of cn/cs/ct?
  // force magnitudes determined by solve_lcp()
  Ravelin::VectorNd cn, cs, ct;
}; // end struct

} // end namespace Moby

#endif // CONTACTPROBLEMDATA_H
