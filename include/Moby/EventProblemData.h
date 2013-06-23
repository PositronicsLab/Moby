#ifndef _MOBY_EVENT_PROBLEM_DATA_H
#define _MOBY_EVENT_PROBLEM_DATA_H

#include <vector>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Moby/Event.h>
#include <Moby/Types.h>

namespace Moby {

struct EventProblemData
{
  // setup reasonable defaults
  EventProblemData()
  {
    reset();
  }

  // copies event problem data
  EventProblemData& copy_from(const EventProblemData& q)
  {
    N_K_TOTAL = q.N_K_TOTAL;
    N_LIN_CONE = q.N_LIN_CONE;
    N_TRUE_CONE = q.N_TRUE_CONE;
    N_LOOPS = q.N_LOOPS;
    N_CONTACTS = q.N_CONTACTS;
    N_CONSTRAINTS = q.N_CONSTRAINTS;
    N_CONSTRAINT_DOF_EXP = q.N_CONSTRAINT_DOF_EXP;
    N_CONSTRAINT_EQNS_EXP = q.N_CONSTRAINT_EQNS_EXP;
    N_CONSTRAINT_DOF_IMP = q.N_CONSTRAINT_DOF_IMP;
    use_kappa = q.use_kappa;
    kappa = q.kappa;

    // copy indices
    CN_IDX = q.CN_IDX;
    CS_IDX = q.CS_IDX;
    CT_IDX = q.CT_IDX;
    NCS_IDX = q.NCS_IDX;
    NCT_IDX = q.NCT_IDX;
    CS_U_IDX = q.CS_U_IDX;
    CT_U_IDX = q.CT_U_IDX;
    L_IDX = q.L_IDX;
    ALPHA_X_IDX = q.ALPHA_X_IDX;
    BETA_X_IDX = q.BETA_X_IDX;
    BETA_T_IDX = q.BETA_T_IDX;
    N_VARS = q.N_VARS;  

    // copy event velocities
    Cn_v = q.Cn_v;
    Cs_v = q.Cs_v;
    Ct_v = q.Ct_v;
    L_v = q.L_v;
    Jx_v = q.Jx_v;
    Dx_v = q.Dx_v;

    // the vector of "super" bodies
    super_bodies = q.super_bodies; 

    // the vectors of events
    events = q.events;
    contact_events = q.contact_events;
    limit_events = q.limit_events;

    // cross-event terms
    Cn_iM_CnT = q.Cn_iM_CnT;
    Cn_iM_CsT = q.Cn_iM_CsT;
    Cn_iM_CtT = q.Cn_iM_CtT;
    Cn_iM_LT = q.Cn_iM_LT;
    Cn_iM_DtT = q.Cn_iM_DtT;
    Cn_iM_JxT = q.Cn_iM_JxT;
    Cn_iM_DxT = q.Cn_iM_DxT;
    Cs_iM_CsT = q.Cs_iM_CsT;
    Cs_iM_CtT = q.Cs_iM_CtT;
    Cs_iM_LT = q.Cs_iM_LT;
    Cs_iM_DtT = q.Cs_iM_DtT;
    Cs_iM_JxT = q.Cs_iM_JxT;
    Cs_iM_DxT = q.Cs_iM_DxT;
    Ct_iM_CtT = q.Ct_iM_CtT;
    Ct_iM_LT = q.Ct_iM_LT;
    Ct_iM_DtT = q.Ct_iM_DtT;
    Ct_iM_JxT = q.Ct_iM_JxT;
    Ct_iM_DxT = q.Ct_iM_DxT;
    L_iM_LT = q.L_iM_LT;
    L_iM_DtT = q.L_iM_DtT;
    L_iM_JxT = q.L_iM_JxT;
    L_iM_DxT = q.L_iM_DxT;
    Dt_iM_DtT = q.Dt_iM_DtT;
    Dt_iM_JxT = q.Dt_iM_JxT;
    Dt_iM_DxT = q.Dt_iM_DxT;
    Jx_iM_JxT = q.Jx_iM_JxT;
    Jx_iM_DxT = q.Jx_iM_DxT;
    Dx_iM_DxT = q.Dx_iM_DxT;

    // copy impulse magnitudes 
    cn = q.cn;
    cs = q.cs;
    l = q.l;
    beta_t = q.beta_t;
    alpha_x = q.alpha_x;
    beta_x = q.beta_x;     

    // copy the working set
    contact_working_set = q.contact_working_set;

    return *this;
  }

  // resets all event problem data
  void reset()
  {
    N_K_TOTAL = N_LIN_CONE = N_TRUE_CONE = N_LOOPS = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_DOF_EXP = N_CONSTRAINT_EQNS_EXP = 0;
    N_CONSTRAINT_DOF_IMP = 0;
    use_kappa = false;
    kappa = 0.0;

    // clear all indices
    N_VARS = 0;
    CS_IDX = CT_IDX = NCS_IDX = NCT_IDX = CS_U_IDX = CT_U_IDX = 0;
    L_IDX = BETA_T_IDX = 0;
    ALPHA_X_IDX = BETA_X_IDX = 0;

    // clear all vectors
    super_bodies.clear();
    events.clear();
    contact_events.clear();
    limit_events.clear();

    // reset all Ravelin::VectorNd sizes
    Cn_v.resize(0);
    Cs_v.resize(0);
    Ct_v.resize(0);
    L_v.resize(0);
    Jx_v.resize(0);
    Dx_v.resize(0);
    cn.resize(0);
    cs.resize(0);
    ct.resize(0);
    l.resize(0);
    beta_t.resize(0);
    alpha_x.resize(0);
    beta_x.resize(0);

    // reset all MatrixN sizes
    Cn_iM_CnT.resize(0,0);
    Cn_iM_CsT.resize(0,0);
    Cn_iM_CtT.resize(0,0);
    Cn_iM_LT.resize(0,0);
    Cn_iM_DtT.resize(0,0);
    Cn_iM_JxT.resize(0,0);
    Cn_iM_DxT.resize(0,0);
    Cs_iM_CsT.resize(0,0);
    Cs_iM_CtT.resize(0,0);
    Cs_iM_LT.resize(0,0);
    Cs_iM_DtT.resize(0,0);
    Cs_iM_JxT.resize(0,0);
    Cs_iM_DxT.resize(0,0);
    Ct_iM_CtT.resize(0,0);
    Ct_iM_LT.resize(0,0);
    Ct_iM_DtT.resize(0,0);
    Ct_iM_JxT.resize(0,0);
    Ct_iM_DxT.resize(0,0);
    L_iM_LT.resize(0,0);
    L_iM_DtT.resize(0,0);
    L_iM_JxT.resize(0,0);
    L_iM_DxT.resize(0,0);
    Dt_iM_DtT.resize(0,0);
    Dt_iM_JxT.resize(0,0);
    Dt_iM_DxT.resize(0,0);
    Jx_iM_JxT.resize(0,0);
    Jx_iM_DxT.resize(0,0);
    Dx_iM_DxT.resize(0,0);

    // reset the working set
    contact_working_set.clear();
  }

  // sets cn, cs, etc. from stacked vectors
  void update_from_stacked(const Ravelin::VectorNd& z)
  {
    cn += z.segment(CN_IDX, CS_IDX);
    l += z.segment(L_IDX, BETA_T_IDX);
    beta_t += z.segment(BETA_T_IDX, ALPHA_X_IDX);
    alpha_x += z.segment(ALPHA_X_IDX, BETA_X_IDX);
    beta_x += z.segment(BETA_X_IDX, N_VARS);

    // setup cs/ct -- first determine linearized friction cone forces
    cs.segment(0, N_LIN_CONE) = z.segment(CS_IDX, CT_IDX);
    cs.segment(0, N_LIN_CONE) -= z.segment(NCS_IDX, NCT_IDX);
    ct.segment(0, N_LIN_CONE) = z.segment(CT_IDX, NCS_IDX);
    ct.segment(0, N_LIN_CONE) -= z.segment(NCT_IDX, CS_U_IDX);

    // now determine quadratic friction cone forces
    cs.segment(N_LIN_CONE, cs.size()) = z.segment(CS_U_IDX, CT_U_IDX);
    ct.segment(N_LIN_CONE, ct.size()) = z.segment(CT_U_IDX, L_IDX);
  }

  // partitions event vectors into contact and limit events
  void partition_events()
  {
    const unsigned UINF = std::numeric_limits<unsigned>::max();
    contact_events.clear();
    limit_events.clear();

    BOOST_FOREACH(Event* e, events)
    {
      if (e->event_type == Event::eContact)
        contact_events.push_back(e);
      else
      {
        assert(e->event_type == Event::eLimit);
        limit_events.push_back(e);
      }
    }

    // now, sort the contact events such that events that use a true friction
    // cone are at the end
    for (unsigned i=0, j=contact_events.size()-1; i< j; )
    {
      if (contact_events[i]->contact_NK == UINF)
      {
        std::swap(contact_events[i], contact_events[j]);
        j--;
      } 
      else
        i++;
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

  // starting index of l in the stacked vector
  unsigned L_IDX;

  // starting index of beta_t in the stacked vector
  unsigned BETA_T_IDX;

  // starting index of alpha_x in the stacked vector
  unsigned ALPHA_X_IDX;

  // starting index of beta_t in the stacked vector
  unsigned BETA_X_IDX;

  // total number of variables
  unsigned N_VARS;

  // the total number of linearized friction tangents for contact events
  unsigned N_K_TOTAL;

  // the number of contacts with linearized friction cones
  unsigned N_LIN_CONE;

  // the number of contacts with true friction cones
  unsigned N_TRUE_CONE;

  // the number of kinematic loops for articulated bodies (only relevant for 
  // advanced joint friction models)
  unsigned N_LOOPS;

  // the number of contacts
  unsigned N_CONTACTS;

  // the number of limits
  unsigned N_LIMITS;

  // the total number of constraints
  unsigned N_CONSTRAINTS;

  // the number of explicit joint constraint degrees-of-freedom used in joint friction computation
  unsigned N_CONSTRAINT_DOF_EXP;

  // the number of explicit joint constraint equations (total)
  unsigned N_CONSTRAINT_EQNS_EXP;

  // the number of implicit joint constraint degrees-of-freedom used in joint friction computation
  unsigned N_CONSTRAINT_DOF_IMP;

  // indication of contacts that the solver is actively considering
  std::vector<bool> contact_working_set;

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies; 

  // the vectors of events
  std::vector<Event*> events, contact_events, limit_events;

  // cross-event terms
  Ravelin::MatrixNd Cn_iM_CnT, Cn_iM_CsT, Cn_iM_CtT, Cn_iM_LT, Cn_iM_DtT, Cn_iM_JxT, Cn_iM_DxT;
  Ravelin::MatrixNd            Cs_iM_CsT, Cs_iM_CtT, Cs_iM_LT, Cs_iM_DtT, Cs_iM_JxT, Cs_iM_DxT;
  Ravelin::MatrixNd                       Ct_iM_CtT, Ct_iM_LT, Ct_iM_DtT, Ct_iM_JxT, Ct_iM_DxT;
  Ravelin::MatrixNd                                  L_iM_LT,  L_iM_DtT,  L_iM_JxT,  L_iM_DxT;
  Ravelin::MatrixNd                                            Dt_iM_DtT, Dt_iM_JxT, Dt_iM_DxT;
  Ravelin::MatrixNd                                            Jx_iM_JxT, Jx_iM_DxT;
  Ravelin::MatrixNd                                                       Dx_iM_DxT;

  // vector-based terms
  Ravelin::VectorNd Cn_v, Cs_v, Ct_v, L_v, Jx_v, Dx_v;

  // kappa term
  double kappa;

  // determines whether to use kappa term
  bool use_kappa;

  // impulse magnitudes determined by solve_qp()
  Ravelin::VectorNd cn, cs, ct, l, beta_t, alpha_x, beta_x;
}; // end struct

} // end namespace Moby

#endif

