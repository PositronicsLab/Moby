#ifndef _MOBY_EVENT_PROBLEM_DATA_H
#define _MOBY_EVENT_PROBLEM_DATA_H

#include <vector>
#include <Moby/MatrixN.h>
#include <Moby/VectorN.h>
#include <Moby/Types.h>

namespace Moby {

class Event;

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

    // copy event velocities
    Jc_v = q.Jc_v;
    Dc_v = q.Dc_v;
    Jl_v = q.Jl_v;
    Jx_v = q.Jx_v;
    Dx_v = q.Dx_v;

    // the vector of "super" bodies
    super_bodies = q.super_bodies; 

    // the vectors of events
    contact_events = q.contact_events;
    limit_events = q.limit_events;
    constraint_events = q.constraint_events;

    // cross-event terms
    Jc_iM_JcT = q.Jc_iM_JcT;
    Jc_iM_DcT = q.Jc_iM_DcT;
    Jc_iM_JlT = q.Jc_iM_JlT;
    Jc_iM_DtT = q.Jc_iM_DtT;
    Jc_iM_JxT = q.Jc_iM_JxT;
    Jc_iM_DxT = q.Jc_iM_DxT;
    Dc_iM_DcT = q.Dc_iM_DcT;
    Dc_iM_JlT = q.Dc_iM_JlT;
    Dc_iM_DtT = q.Dc_iM_DtT;
    Dc_iM_JxT = q.Dc_iM_JxT;
    Dc_iM_DxT = q.Dc_iM_DxT;
    Jl_iM_JlT = q.Jl_iM_JlT;
    Jl_iM_DtT = q.Jl_iM_DtT;
    Jl_iM_JxT = q.Jl_iM_JxT;
    Jl_iM_DxT = q.Jl_iM_DxT;
    Dt_iM_DtT = q.Dt_iM_DtT;
    Dt_iM_JxT = q.Dt_iM_JxT;
    Dt_iM_DxT = q.Dt_iM_DxT;
    Jx_iM_JxT = q.Jx_iM_JxT;
    Jx_iM_DxT = q.Jx_iM_DxT;
    Dx_iM_DxT = q.Dx_iM_DxT;

    // copy impulse magnitudes 
    alpha_c = q.alpha_c;
    beta_c = q.beta_c;
    alpha_l = q.alpha_l;
    beta_t = q.beta_t;
    alpha_x = q.alpha_x;
    beta_x = q.beta_x;     
  }

  // resets all event problem data
  void reset()
  {
    N_K_TOTAL = N_LIN_CONE = N_TRUE_CONE = N_LOOPS = N_CONTACTS = 0;
    N_CONSTRAINTS = N_CONSTRAINT_DOF_EXP = N_CONSTRAINT_EQNS_EXP = 0;
    N_CONSTRAINT_DOF_IMP = 0;
    use_kappa = false;
    kappa = (Real) 0.0;

    // clear all vectors
    super_bodies.clear();
    contact_events.clear();
    limit_events.clear();
    constraint_events.clear();

    // reset all VectorN sizes
    Jc_v.resize(0);
    Dc_v.resize(0);
    Jl_v.resize(0);
    Jx_v.resize(0);
    Dx_v.resize(0);
    alpha_c.resize(0);
    beta_c.resize(0);
    alpha_l.resize(0);
    beta_t.resize(0);
    alpha_x.resize(0);
    beta_x.resize(0);

    // reset all MatrixN sizes
    Jc_iM_JcT.resize(0,0);
    Jc_iM_DcT.resize(0,0);
    Jc_iM_JlT.resize(0,0);
    Jc_iM_DtT.resize(0,0);
    Jc_iM_JxT.resize(0,0);
    Jc_iM_DxT.resize(0,0);
    Dc_iM_DcT.resize(0,0);
    Dc_iM_JlT.resize(0,0);
    Dc_iM_DtT.resize(0,0);
    Dc_iM_JxT.resize(0,0);
    Dc_iM_DxT.resize(0,0);
    Jl_iM_JlT.resize(0,0);
    Jl_iM_DtT.resize(0,0);
    Jl_iM_JxT.resize(0,0);
    Jl_iM_DxT.resize(0,0);
    Dt_iM_DtT.resize(0,0);
    Dt_iM_JxT.resize(0,0);
    Dt_iM_DxT.resize(0,0);
    Jx_iM_JxT.resize(0,0);
    Jx_iM_DxT.resize(0,0);
    Dx_iM_DxT.resize(0,0);
  }

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

  // the vector of "super" bodies
  std::vector<DynamicBodyPtr> super_bodies; 

  // the vectors of events
  std::vector<Event*> contact_events, limit_events, constraint_events;

  // cross-event terms
  MatrixN Jc_iM_JcT, Jc_iM_DcT, Jc_iM_JlT, Jc_iM_DtT, Jc_iM_JxT, Jc_iM_DxT;
  MatrixN            Dc_iM_DcT, Dc_iM_JlT, Dc_iM_DtT, Dc_iM_JxT, Dc_iM_DxT;
  MatrixN                       Jl_iM_JlT, Jl_iM_DtT, Jl_iM_JxT, Jl_iM_DxT;
  MatrixN                                  Dt_iM_DtT, Dt_iM_JxT, Dt_iM_DxT;
  MatrixN                                             Jx_iM_JxT, Jx_iM_DxT;
  MatrixN                                                        Dx_iM_DxT;

  // vector-based terms
  VectorN Jc_v, Dc_v, Jl_v, Jx_v, Dx_v;

  // kappa term
  Real kappa;

  // determines whether to use kappa term
  bool use_kappa;

  // impulse magnitudes determined by solve_qp()
  VectorN alpha_c, beta_c, alpha_l, beta_t, alpha_x, beta_x;
}; // end struct

} // end namespace Moby

#endif

