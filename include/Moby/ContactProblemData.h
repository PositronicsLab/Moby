#ifndef CONTACTPROBLEMDATA_H
#define CONTACTPROBLEMDATA_H

#include <vector>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/Event.h>
#include <Moby/Types.h>

namespace Moby {

struct ContactProblemData
{
  // setup reasonable defaults
  ContactProblemData()
  {
    reset();
  }

  // copies contact problem data
  ContactProblemData& operator=(const ContactProblemData& q)
  {
    return copy_from(q);
  }

  // copies contact problem data
  ContactProblemData& copy_from(const ContactProblemData& q)
  {
    N_K_TOTAL = q.N_K_TOTAL;
    N_LOOPS = q.N_LOOPS;
    N_CONTACTS = q.N_CONTACTS;

    // copy indices
    CN_IDX = q.CN_IDX;
    CS_IDX = q.CS_IDX;
    CT_IDX = q.CT_IDX;
    NCS_IDX = q.NCS_IDX;
    NCT_IDX = q.NCT_IDX;
    CS_U_IDX = q.CS_U_IDX;
    CT_U_IDX = q.CT_U_IDX;
    ALPHA_X_IDX = q.ALPHA_X_IDX;
    BETA_X_IDX = q.BETA_X_IDX;
    BETA_T_IDX = q.BETA_T_IDX;
    N_VARS = q.N_VARS;

    // copy contact accelerations
    Cn_a = q.Cn_a;
    Cs_a = q.Cs_a;
    Ct_a = q.Ct_a;

    // the vector of "super" bodies
    super_bodies = q.super_bodies;

    // the vectors of contacts 
    contacts = q.contacts;

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
    fn = q.fn;
    fs = q.fs;
    ft = q.ft;

    // copy the working set
    contact_working_set = q.contact_working_set;

    return *this;
  }

  // resets all contact problem data
  void reset()
  {
    N_K_TOTAL = N_LIN_CONE = N_TRUE_CONE = N_LOOPS = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_DOF_EXP = N_CONSTRAINT_EQNS_IMP = 0;
    N_CONSTRAINT_DOF_IMP = 0;

    // clear all indices
    N_VARS = 0;
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = CS_U_IDX = CT_U_IDX = 0;
    BETA_T_IDX = 0;
    ALPHA_X_IDX = BETA_X_IDX = 0;

    // clear all vectors
    super_bodies.clear();
    contacts.clear();

    // reset all Ravelin::VectorNd sizes
    Cn_v.resize(0);
    Cs_v.resize(0);
    Ct_v.resize(0);
    Jx_v.resize(0);
    Dx_v.resize(0);
    cn.resize(0);
    cs.resize(0);
    ct.resize(0);
    beta_t.resize(0);
    alpha_x.resize(0);
    beta_x.resize(0);

    // reset all MatrixN sizes
    Cn_iM_CnT.resize(0,0);
    Cn_iM_CsT.resize(0,0);
    Cn_iM_CtT.resize(0,0);
    Cn_iM_DtT.resize(0,0);
    Cn_iM_JxT.resize(0,0);
    Cn_iM_DxT.resize(0,0);
    Cs_iM_CnT.resize(0,0);
    Cs_iM_CsT.resize(0,0);
    Cs_iM_CtT.resize(0,0);
    Cs_iM_DtT.resize(0,0);
    Cs_iM_JxT.resize(0,0);
    Cs_iM_DxT.resize(0,0);
    Ct_iM_CnT.resize(0,0);
    Ct_iM_CsT.resize(0,0);
    Ct_iM_CtT.resize(0,0);
    Ct_iM_DtT.resize(0,0);
    Ct_iM_JxT.resize(0,0);
    Ct_iM_DxT.resize(0,0);
    Dt_iM_CnT.resize(0,0);
    Dt_iM_CsT.resize(0,0);
    Dt_iM_CtT.resize(0,0);
    Dt_iM_DtT.resize(0,0);
    Dt_iM_JxT.resize(0,0);
    Dt_iM_DxT.resize(0,0);
    Jx_iM_CnT.resize(0,0);
    Jx_iM_CsT.resize(0,0);
    Jx_iM_CtT.resize(0,0);
    Jx_iM_DtT.resize(0,0);
    Jx_iM_JxT.resize(0,0);
    Jx_iM_DxT.resize(0,0);
    Dx_iM_CnT.resize(0,0);
    Dx_iM_CsT.resize(0,0);
    Dx_iM_CtT.resize(0,0);
    Dx_iM_DtT.resize(0,0);
    Dx_iM_JxT.resize(0,0);
    Dx_iM_DxT.resize(0,0);

    // reset the working set
    contact_working_set.clear();
  }

  // sets cn, cs, etc. from stacked vectors
  void update_from_stacked(const Ravelin::VectorNd& z)
  {
    cn += z.segment(CN_IDX, CS_IDX);
    beta_t += z.segment(BETA_T_IDX, ALPHA_X_IDX);
    alpha_x += z.segment(ALPHA_X_IDX, BETA_X_IDX);
    beta_x += z.segment(BETA_X_IDX, N_VARS);

    // setup cs/ct -- first determine linearized friction cone forces
    cs.segment(0, N_LIN_CONE) = z.segment(CS_IDX, CT_IDX);
    cs.segment(0, N_LIN_CONE) -= z.segment(NCS_IDX, NCT_IDX);
    ct.segment(0, N_LIN_CONE) = z.segment(CT_IDX, NCS_IDX);
    ct.segment(0, N_LIN_CONE) -= z.segment(NCT_IDX, CS_U_IDX);
  }

  // sets stacked vector from cn, cs, etc.
  Ravelin::VectorNd& to_stacked(Ravelin::VectorNd& z)
  {
    z.set_sub_vec(CN_IDX, cn);
    z.set_sub_vec(CS_IDX, cs.segment(0, N_LIN_CONE));
    z.set_sub_vec(CT_IDX, ct.segment(0, N_LIN_CONE));
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
    z.set_sub_vec(BETA_T_IDX, beta_t);
    z.set_sub_vec(ALPHA_X_IDX, alpha_x);
    z.set_sub_vec(BETA_X_IDX, beta_x);
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
  unsigned N_STICKING; // NOTE!:THis value is not calculated yet

  // the number of contacts
  unsigned N_CONTACTS;

  // indication of contacts that the solver is actively considering
  std::vector<bool> contact_working_set;

  // indication if contact is sticking
  std::vector<bool> contact_sticking; // NOTE!:THis value is not calculated yet

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies;

  // the vectors of acceleration-level events 
  std::vector<Event*> events;

  // cross-contact terms
  Ravelin::MatrixNd Cn_iM_CnT, Cn_iM_CsT, Cn_iM_CtT, Cn_iM_muS_CqT;
  Ravelin::MatrixNd Cs_iM_CnT, Cs_iM_CsT, Cs_iM_CtT, Cs_iM_muS_CqT;
  Ravelin::MatrixNd Ct_iM_CnT, Ct_iM_CsT, Ct_iM_CtT, Ct_iM_muS_CqT;

  // the Sticking Friction matrix N_STICKING x N_CONTACTS
  Ravelin::MatrixNd muK;

  // EMD: it's not necessary to set these up explicitly, and they're sparse.
  //      if you look at my 'A' matrix in my QP, you'll see how I do this 
  // the Tangent Friction Polygon constraint matrix contact_NK x N_CONTACTS
  Ravelin::MatrixNd Xs, Xt;

  // vector-based terms
  Ravelin::VectorNd Cn_a, Cs_a, Ct_a

  // EMD: why use this instead of cn/cs/ct?
  // impulse magnitudes determined by solve_qp()
  Ravelin::VectorNd fn, fs, ft;
}; // end struct

} // end namespace Moby

#endif // CONTACTPROBLEMDATA_H
