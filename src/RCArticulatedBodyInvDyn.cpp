/****************************************************************************
 * Copyright 2016 Samuel Zapolsky and Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/Jointd.h>
#include <Moby/insertion_sort>
#include <Moby/RCArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/RCArticulatedBodyInvDyn.h>

using std::vector;
using std::map;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Gets the body indices and the total number of degrees of freedom
unsigned RCArticulatedBodyInvDyn::get_body_indices(const std::vector<Moby::UnilateralConstraint>& c, map<shared_ptr<DynamicBodyd>, unsigned>& index_start)
{
  unsigned NDOFS = 0;

  // clear the map
  index_start.clear();

  // populate the map
  for (unsigned i=0; i< c.size(); i++)
  {
    // verify that this is a contact constraint
    assert(c[i].constraint_type == UnilateralConstraint::eContact);

    // get the two bodies involved
    shared_ptr<SingleBodyd> b1 = c[i].contact_geom1->get_single_body(); 
    shared_ptr<SingleBodyd> b2 = c[i].contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> sb1 = ImpactConstraintHandler::get_super_body(b1);
    shared_ptr<DynamicBodyd> sb2 = ImpactConstraintHandler::get_super_body(b2);

    // store sb1's starting index, if necessary 
    if (index_start.find(sb1) == index_start.end())
    {
      index_start[sb1] = NDOFS;
      NDOFS += sb1->num_generalized_coordinates(DynamicBodyd::eSpatial);
    }

    // store sb2's starting index, if necessary 
    if (index_start.find(sb2) == index_start.end())
    {
      index_start[sb2] = NDOFS;
      NDOFS += sb2->num_generalized_coordinates(DynamicBodyd::eSpatial);
    }
  }

  return NDOFS;
}

/// Computes contact Jacobians from the unilateral contact constraints
void RCArticulatedBodyInvDyn::calc_contact_jacobians(const std::vector<Moby::UnilateralConstraint>& c, Ravelin::MatrixNd& N,Ravelin::MatrixNd& S, Ravelin::MatrixNd& T, const map<shared_ptr<DynamicBodyd>, unsigned>& index_start, unsigned NDOFS)
{
  // get the number of contact constraints
  unsigned NC = c.size();
  if(NC==0) return;

  // size N, S, and T
  N.set_zero(NDOFS,NC);
  S.set_zero(NDOFS,NC);
  T.set_zero(NDOFS,NC);
  
  // Contact Jacobian [GLOBAL frame]
  for(unsigned i=0;i<NC;i++){
   
    // get the two single bodies involved in the contact 
    shared_ptr<SingleBodyd> b1 = c[i].contact_geom1->get_single_body();
    shared_ptr<SingleBodyd> b2 = c[i].contact_geom2->get_single_body();

    // get the two super bodies
    shared_ptr<DynamicBodyd> sb1 = ImpactConstraintHandler::get_super_body(b1);
    shared_ptr<DynamicBodyd> sb2 = ImpactConstraintHandler::get_super_body(b2);

    // setup the impulse frame 
    boost::shared_ptr<const Ravelin::Pose3d>
    impulse_frame(new Ravelin::Pose3d(Ravelin::Quatd::identity(),c[i].contact_point.data(),Moby::GLOBAL));
    
    // setup the normal and two tangent directions 
    const Ravelin::Vector3d& normal = c[i].contact_normal;
    const Ravelin::Vector3d& tan1 = c[i].contact_tan1;
    const Ravelin::Vector3d& tan2 = c[i].contact_tan2;
 
    // get the starting and ending indices
    assert(index_start.find(sb1) != index_start.end());
    unsigned st_idx = index_start.find(sb1)->second;
    unsigned ed_idx = st_idx + sb1->num_generalized_coordinates(DynamicBodyd::eSpatial);
   
    // compute the Jacobian for the first body 
    sb1->calc_jacobian(sb1->get_gc_pose(), impulse_frame,dynamic_pointer_cast<DynamicBodyd>(b1),_workM);

    // get the linear components of the Jacobian
    SharedMatrixNd J1 = _workM.block(0,3,0,_workM.columns());

    // Normal direction
    J1.transpose_mult(normal,_workv);
    N.column(i).segment(st_idx,ed_idx) = _workv;
    
    // 1st tangent
    J1.transpose_mult(tan1,_workv);
    S.column(i).segment(st_idx,ed_idx) = _workv;
    
    // 2nd tangent
    J1.transpose_mult(tan2,_workv);
    T.column(i).segment(st_idx,ed_idx) = _workv;

    // compute the Jacobian for the second body
    sb2->calc_jacobian(sb2->get_gc_pose(),impulse_frame,dynamic_pointer_cast<DynamicBodyd>(b2),_workM);

    // get the linear components of the Jacobian
    MatrixNd J2 = _workM.block(0,3,0,_workM.columns());

    // get the starting and ending indices
    assert(index_start.find(sb2) != index_start.end());
    st_idx = index_start.find(sb2)->second;
    ed_idx = st_idx + sb2->num_generalized_coordinates(DynamicBodyd::eSpatial);

    // Normal direction
    J2.transpose_mult(-normal,_workv);
    N.column(i).segment(st_idx,ed_idx) = _workv;
    
    // 1st tangent
    J2.transpose_mult(-tan1,_workv);
    S.column(i).segment(st_idx,ed_idx) = _workv;
    
    // 2nd tangent
    J2.transpose_mult(-tan2,_workv);
    T.column(i).segment(st_idx,ed_idx) = _workv;
  }

  FILE_LOG(LOG_DYNAMICS) << "N: " << std::endl << N;
  FILE_LOG(LOG_DYNAMICS) << "S: " << std::endl << S;
  FILE_LOG(LOG_DYNAMICS) << "T: " << std::endl << T;
}

void RCArticulatedBodyInvDyn::jacobian_ST_to_D(const Ravelin::MatrixNd& S,const Ravelin::MatrixNd& T, Ravelin::MatrixNd& D)
{
  unsigned NC = S.columns();
  unsigned N = S.rows();
  D.resize(N,NC*4);
  D.set_sub_mat(0,0,S);
  D.set_sub_mat(0,NC,T);
  _workM = S;
  _workM.negate();
  D.set_sub_mat(0,NC*2,_workM);
  _workM = T;
  _workM.negate();
  D.set_sub_mat(0,NC*3,_workM);
}

// ============================================================================
// ================================ CONTROLLER ================================
// ============================================================================

//     void (*controller)(boost::shared_ptr<ControlledBody>, double, void*);

void RCArticulatedBodyInvDyn::calc_inv_dyn(RCArticulatedBodyPtr body, const RCArticulatedBodyInvDynData& idyn_data, double dt, VectorNd& u)
{
  // get the total number of degrees-of-freedom for all bodies involved
  // and setup a mapping for components in N, S, T
  map<shared_ptr<DynamicBodyd>, unsigned> index_start;
  const unsigned NDOFS = get_body_indices(idyn_data.constraints, index_start);

  // get the bodies in order
  std::vector<shared_ptr<DynamicBodyd> > body_ordering;
  for (map<shared_ptr<DynamicBodyd>, unsigned>::const_iterator i = index_start.begin(); i != index_start.end(); i++)
    body_ordering.push_back(i->first);

  // get the generalized velocities and forces
  VectorNd fext(NDOFS), v(NDOFS);
  fext.set_zero();
  v.set_zero();
  for (map<shared_ptr<DynamicBodyd>, unsigned>::const_iterator i = index_start.begin(); i != index_start.end(); i++)
  {
    const unsigned ST_IDX = i->second;
    const unsigned ED_IDX = ST_IDX + i->first->num_generalized_coordinates(DynamicBodyd::eSpatial);
    SharedVectorNd fext_sub = fext.segment(ST_IDX, ED_IDX);
    SharedVectorNd v_sub = v.segment(ST_IDX, ED_IDX);
    i->first->get_generalized_velocity(DynamicBodyd::eSpatial, v_sub);
    i->first->get_generalized_forces(fext_sub);
  }

  // get the generalized inertias
  std::vector<MatrixNd> M(body_ordering.size());
  for (unsigned i=0; i< body_ordering.size(); i++)
    body_ordering[i]->get_generalized_inertia(M[i]);

  // look for no constraints
  if (idyn_data.constraints.empty())
  {
    unsigned NGC = body->num_generalized_coordinates(DynamicBodyd::eSpatial);
    unsigned NQ = body->num_joint_dof_explicit();

    // get the generalized inertia
    M.push_back(MatrixNd());
    body->get_generalized_inertia(M.front());

    // get the controllable joint indices
    std::vector<unsigned> cj_indices;
    const vector<shared_ptr<Jointd> >& joints = body->get_joints();
    for (unsigned i=0; i< joints.size(); i++)
      cj_indices.push_back(joints[i]->get_coord_index());
    std::sort(cj_indices.begin(), cj_indices.end());

    // setup P and P'
    _P.set_zero(NQ,NGC);
    for (unsigned i=0; i< cj_indices.size(); i++)
      _P(i,cj_indices[i]) = 1.0;

    // compute inverse dynamics forces
    VectorNd tau;
    _P.transpose_mult(idyn_data.qdd_des, tau);
    M.front().mult(tau, u);
    return;
  }

  // form Jacobians
  Ravelin::MatrixNd N,S,T,D;
  calc_contact_jacobians(idyn_data.constraints,N,S,T,index_start,NDOFS);
  jacobian_ST_to_D(S,T,D);

  // setup MU
  std::vector<double> mu(idyn_data.constraints.size());
  for( unsigned i = 0 ; i < idyn_data.constraints.size() ; i++ )
    mu[i] = idyn_data.constraints[i].contact_mu_coulomb;

  // prepare to calculate contact forces and inverse dynamics forces
  VectorNd tau, cf;

  // call the appropriate inverse dynamics algorithm
  inverse_dynamics_comp(body, index_start, v, idyn_data.qdd_des, M, N, D, fext, dt, mu, tau, cf);

  // compute the generalized torques
  _P.transpose_mult(tau, u);
}

// implements a controller callback for Moby
/*
Ravelin::VectorNd& controller_callback(boost::shared_ptr<Moby::ControlledBody> cbp,Ravelin::VectorNd& control_force, double t, void*)
{
  NDOFS = abrobot->num_joint_dof();

  
  static unsigned long long ITER = 0;
  static double last_time = -0.001;
  double dt = t - last_time;
  last_time = t;
  
  /////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// Get State: ///////////////////////////////////
  
//  Ravelin::VectorNd generalized_q(NDOFS+6);
  Ravelin::VectorNd generalized_qd(NDOFS+6), generalized_fext(NDOFS+6);
  
    {
//      // Arm
//      abrobot->get_generalized_coordinates_euler(WORKV);
//      generalized_q.segment(0,NDOFS) = WORKV;
      
      abrobot->get_generalized_velocity(Ravelin::DynamicBodyd::eSpatial,WORKV);
      generalized_qd.segment(0,NDOFS) = WORKV;
      
      abrobot->get_generalized_forces(WORKV);
      generalized_fext.segment(0,NDOFS) = WORKV;
      
      // Object
//      object->get_generalized_coordinates_euler(WORKV);
//      generalized_q.segment(NDOFS,NDOFS+7) = WORKV;
      
      object->get_generalized_velocity(Ravelin::DynamicBodyd::eSpatial,WORKV);
      generalized_qd.segment(NDOFS,NDOFS+6) = WORKV;
      
      object->get_generalized_forces(WORKV);
      generalized_fext.segment(NDOFS,NDOFS+6) = WORKV;
    }
  
  /////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// Get contacts: ///////////////////////////////////

    // pointer to the simulator
    boost::shared_ptr<Moby::Simulator> sim = sim_weak_ptr.lock();
    
    boost::shared_ptr<Moby::ConstraintSimulator> csim;
    csim = boost::dynamic_pointer_cast<Moby::ConstraintSimulator>(sim);
    
    std::vector<Moby::UnilateralConstraint>& rigid_constraints = csim->get_rigid_constraints();
    std::vector<Moby::UnilateralConstraint>& compliant_constraints = csim->get_compliant_constraints();
    std::vector<Moby::UnilateralConstraint> e;
    std::map<std::string, Moby::UnilateralConstraint> contacts;

    e.insert(e.end(), rigid_constraints.begin(), rigid_constraints.end());
    e.insert(e.end(), compliant_constraints.begin(), compliant_constraints.end());
    for(unsigned i=0;i<e.size();i++){
      if (e[i].constraint_type == Moby::UnilateralConstraint::eContact)
      {
        boost::shared_ptr<Ravelin::SingleBodyd> b1 = e[i].contact_geom1->get_single_body();
        boost::shared_ptr<Ravelin::SingleBodyd> b2 = e[i].contact_geom2->get_single_body();
        
        Ravelin::Vector3d
        normal = e[i].contact_normal,
        tangent = e[i].contact_tan1,
        impulse = Ravelin::Vector3d(0,0,0);//e[i].contact_impulse.get_linear();
        impulse.pose = e[i].contact_point.pose;
        tangent.normalize();
        normal.normalize();
        
        bool compliant =
        (e[i].compliance == Moby::UnilateralConstraint::eCompliant)? true : false;
        
        
        FILE_LOG(LOG_DYNAMICS) << "MOBY: contact-ids: " << b1->body_id << " <-- " << b2->body_id ;
        FILE_LOG(LOG_DYNAMICS) << "MOBY: compliant: " << compliant;
        FILE_LOG(LOG_DYNAMICS) << "MOBY: normal: " << normal;
        FILE_LOG(LOG_DYNAMICS) << "MOBY: tangent: " << tangent;
        FILE_LOG(LOG_DYNAMICS) << "MOBY: point: " << e[i].contact_point;
        
        if(b2->body_id.compare("BLOCK") == 0){
          continue;
//          throw std::runtime_error("sb2 nees to be a FINGER_, not the BLOCK");
        }

        contacts[b2->body_id] = e[i];
        }
    }
  
  std::map<std::string, unsigned> finger_index;
  finger_index["FINGER_0"] = 0;
  finger_index["FINGER_1"] = 1;
  finger_index["FINGER_2"] = 2;
  finger_index["FINGER_3"] = 3;

  // Create active contacts vector
  std::vector<Moby::UnilateralConstraint> c;
  std::vector<unsigned> indices;
  for( std::map<std::string, Moby::UnilateralConstraint>::iterator it = contacts.begin(); it != contacts.end(); ++it ) {
    
    boost::shared_ptr<Ravelin::SingleBodyd> b1 = it->second.contact_geom1->get_single_body();
    boost::shared_ptr<Ravelin::SingleBodyd> b2 = it->second.contact_geom2->get_single_body();
    
    if(b2->body_id.compare("BLOCK") == 0){
//      continue;
      throw std::runtime_error("sb2 nees to be a FINGER_, not the BLOCK");
    }
    
    c.push_back( it->second );
    indices.push_back(finger_index[b2->body_id]);
  }
  
  /////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// Get JACOBIANS: ///////////////////////////////////

  unsigned NC = c.size();

  // Make Jacobians
  Ravelin::MatrixNd N,S,T,D,MU;
  calc_contact_jacobians(c,N,S,T);
  D = jacobian_ST_to_D(S,T);

  MU.set_zero(NC,2);
  for( unsigned i = 0 ; i < c.size() ; i++ ) {
    std::fill(MU.row(i).begin(),MU.row(i).end(),c[i].contact_mu_coulomb);
  }

  Ravelin::VectorNd qdd_des = Ravelin::VectorNd::zero(NDOFS);

  Ravelin::MatrixNd M = Ravelin::MatrixNd::zero(NDOFS+6,NDOFS+6);
  
  abrobot->get_generalized_inertia(WORKM);
  M.block(0,NDOFS,0,NDOFS) = WORKM;
  
  object->get_generalized_inertia(WORKM);

  M.block(NDOFS,NDOFS+6,NDOFS,NDOFS+6) = WORKM;

  assert(M.rows() == N.rows());
  assert(generalized_qd.rows() == N.rows());
  
  // reset control force.
  control_force.set_zero(generalized_qd.rows());
  
  
  ////////////////////////// simulator DT IDYN //////////////////////////////

  // IDYN MAXIMAL DISSIPATION MODEL
  std::vector<std::string> controller_name;
//  controller_name.push_back("CFQP1");
//  controller_name.push_back("SCFQP");
  controller_name.push_back("CFLCP");
  controller_name.push_back("NSQP");
//  controller_name.push_back("SNSQP");
  controller_name.push_back("NSLCP");
  controller_name.push_back("CFQP");
  controller_name.push_back("NOCP");

  std::map<std::string,Ravelin::VectorNd> cf_map,uff_map;
  
  // Index of foot check (continuity of contacts)
  static std::vector<unsigned> last_indices = indices;
  static unsigned last_NC = NC;
  bool SAME_INDICES = false;
  if(NC == last_NC){
    SAME_INDICES = true;
    for (unsigned i=0; i<indices.size(); i++) {
      if(last_indices[i] != indices[i]){
        SAME_INDICES = false;
        break;
      }
    }
  }
  last_NC = NC;
  last_indices = indices;

  unsigned NUM_JOINT_DOFS = NDOFS;
  bool USE_LAST_CFS = false;
  
  for (std::vector<std::string>::iterator it=controller_name.begin();
       it!=controller_name.end(); it++) {
    
    bool solve_flag = false;
    
    Ravelin::VectorNd id = Ravelin::VectorNd::zero(NUM_JOINT_DOFS);
    Ravelin::VectorNd cf;// = Ravelin::VectorNd::zero(NC*5);
    
    std::string& name = (*it);
    
    FILE_LOG(LOG_DYNAMICS) << "CONTROLLER: " << name;
    FILE_LOG(LOG_DYNAMICS) << "USE_LAST_CFS: " << USE_LAST_CFS;
    //
    if (NC == 0) {
      solve_flag = inverse_dynamics_no_contact(generalized_qd,qdd_des,M,generalized_fext,dt,id);
      cf /= dt;
      cf_map[name] = cf;
      uff_map[name] = id;
      control_force = id;
      return control_force;
    }
//    else if(USE_LAST_CFS){
//      cf = cf_init;
//      cf *= dt;
//      FILE_LOG(LOG_DYNAMICS) << "USING LAST CFS: " << cf;
//      //      solve_flag = inverse_dynamics(qdd_des,M,N,D,generalized_fext,dt,id,cf);
//      solve_flag = inverse_dynamics_one_stage(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,MU,id,cf,indices,NC,SAME_INDICES);
//      FILE_LOG(LOG_DYNAMICS) << "SAME: " << cf;
//    }
    else {
#ifdef TIMING
      struct timeval start_t;
      struct timeval end_t;
      gettimeofday(&start_t, NULL);
#endif
      try{
        if(name.compare("NSQP") == 0){
          solve_flag = inverse_dynamics_no_slip(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,id,cf,indices,NC,SAME_INDICES);
        } else if(name.compare("NSLCP") == 0){
          solve_flag = inverse_dynamics_no_slip_fast(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,id,cf,false,indices,NC,SAME_INDICES);
        } else if(name.compare("CFQP") == 0){    // IDYN QP
          solve_flag = inverse_dynamics_two_stage(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,MU,id,cf,indices,NC,SAME_INDICES);
        } else if(name.compare("CFQP1") == 0){    // IDYN QP
          solve_flag = inverse_dynamics_one_stage(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,MU,id,cf,indices,NC,SAME_INDICES);
        } else if(name.compare("CFLCP") == 0){
          solve_flag = inverse_dynamics_ap(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,MU,id,cf);
        }  else if(name.compare("SCFQP") == 0){
          double damping = 0;
          solve_flag = inverse_dynamics_two_stage_simple(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,MU,id,cf, damping);
        } else if(name.compare("SNSQP") == 0){
          double damping = 0;
          solve_flag = inverse_dynamics_two_stage_simple_no_slip(generalized_qd,qdd_des,M,N,D,generalized_fext,dt,id,cf, damping);
        } else {
          solve_flag = inverse_dynamics_no_contact(generalized_qd,qdd_des,M,generalized_fext,dt,id);
          control_force = id;
        }
      } catch (std::exception e){
        solve_flag = false;
      }
        

    }
    
#if 0
    //// Check for Normal direction torque chatter
    if (NC > 0 && solve_flag){
      static std::map<std::string,Ravelin::VectorNd> last_cf;
      std::map<std::string,Ravelin::VectorNd>::iterator it = last_cf.find(name);
      if(it != last_cf.end()){
        if(last_cf[name].rows() == NC){
          Ravelin::VectorNd diff_cf  = (*it).second;
          diff_cf -= cf.segment(0,NC);
          if(diff_cf.norm() > 0.01)
            FILE_LOG(LOG_DYNAMICS) << "-- Torque chatter detected " << t;
        }
      } else {
        last_cf[name] = cf.segment(0,NC);
      }
    }
    
    double normal_sum = 0;
#endif
    
    cf /= dt;
    cf_map[name] = cf;
    uff_map[name] = id;
  }
  
  for (std::map<std::string,Ravelin::VectorNd>::iterator it=uff_map.begin();it!=uff_map.end();it++){
    std::string name = it->first;
    Ravelin::VectorNd& cf = cf_map[name];
//    control_force = uff_map[name];
  }

  return control_force;
}
*/

/// Computes inverse dynamics with contact under the assumption of complementarity
bool RCArticulatedBodyInvDyn::inverse_dynamics_comp(RCArticulatedBodyPtr body, const map<shared_ptr<DynamicBodyd>, unsigned>& start_indices, const VectorNd& vel, const VectorNd& qdd, const std::vector<MatrixNd>& M,const MatrixNd& NT, const MatrixNd& D_, const VectorNd& fext, double dt, const std::vector<double>& mu, VectorNd& x, VectorNd& cf)
{
  FILE_LOG(LOG_DYNAMICS) << ">> inverse_dynamics_ap() entered" << std::endl;
  
  // get number of degrees of freedom and number of contact points
  unsigned NGC = fext.size();
  unsigned NQ = body->num_joint_dof_explicit();
  unsigned NC = NT.columns();
  unsigned NK = D_.columns()/NC;

  // get the body start index
  const unsigned ST_IDX = start_indices.find(body)->second; 

  // get the controllable joint indices
  std::vector<unsigned> cj_indices;
  const vector<shared_ptr<Jointd> >& joints = body->get_joints();
  for (unsigned i=0; i< joints.size(); i++)
    cj_indices.push_back(ST_IDX + joints[i]->get_coord_index());
  std::sort(cj_indices.begin(), cj_indices.end());

  // setup P and P'
  _P.set_zero(NQ,NGC);
  for (unsigned i=0; i< cj_indices.size(); i++)
    _P(i,cj_indices[i]) = 1.0;
  MatrixNd::transpose(_P, _PT);
  
  // setup R
  _R.resize(NGC, NC*5 );
  _DT = D_;
  _D.resize(NGC,NC*NK);
  _N = NT;
  MatrixNd::transpose(_DT, _D);
  MatrixNd::transpose(NT, _N);
  _R.block(0,NGC,0,NC) = NT;
  _R.block(0,NGC,NC,NC*5) = _DT;
  
  // compute D, E, and F
  bool pass = factor_chols(M, _M_chols);
  assert(pass);
  
  Ravelin::VectorNd vq(NQ);
  vel.get_sub_vec(0,NQ,vq);
  
  Ravelin::VectorNd v = fext;
  solve_chol_fast(_M_chols, v);
  v *= dt;
  v += vel;
  
  Ravelin::VectorNd vqstar;
  //  ((vqstar = qdd) *= dt) += vq;
  ((vqstar = qdd) *= dt) += v.segment(0,NQ);
  //////////////////////////////////////////////////////////////

 
  // S
  solve_chol_fast(_M_chols,_workM = _DT);
  _N.mult(_workM,_Cn_iM_CdT);
  _D.mult(_workM,_Cd_iM_CdT);
  
  // N
  solve_chol_fast(_M_chols,_workM = NT);
  _N.mult(_workM,_Cn_iM_CnT);
  _D.mult(_workM,_Cd_iM_CnT);
  
  // P
  solve_chol_fast(_M_chols,_workM = _PT);
  _N.mult(_workM,_Cn_iM_JxT);
  _D.mult(_workM,_Cd_iM_JxT);
  _P.mult(_workM,_Jx_iM_JxT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = DT);
  //  P.mult(_workM,Jx_iM_CdT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = NT);
  //  P.mult(_workM,Jx_iM_CnT);
  
  Ravelin::VectorNd Cd_v, Cn_v, Jx_v;
  _D.mult(v,Cd_v);
  _N.mult(v,Cn_v);
  _P.mult(v,Jx_v);
  
  /////////////////////////////////////////////////////////////////////////////
  /// Solve System
  ///
  const unsigned J_IDX = 0;
  const unsigned N_IDX = 0;
  const unsigned D_IDX = NC;
  const unsigned E_IDX = D_IDX + NC*NK;
  
  // we do this by solving the MLCP:
  // |  A  C  | | x | + | g | = | 0 |
  // |  D  B  | | y |   | h |   | 0 |
  
  // g = [-M*v'; 0; 0; -vq_err]
  // h = [ 0 ];
  // A = [ M -P';
  //       P  0 ];
  
  // B = [ 0  0  0
  //         0  0  E
  //        mu -E' 0 ];
  
  // C = [ -N' -D' 0 ;
  //        0   0  0 ];
  // D = -C'
  // x = [ v^+; tau ]
  // y = [ cN, cD, lambda]
  
  // lower left & upper right block of matrix
  if (NK > 4)
    _LL.set_zero(NC*NK,NC+NC*NK);
  else
    _LL.set_zero(NC,NC+NC*NK);
  _UR.set_zero(_LL.columns(),_LL.rows());
  for(unsigned i=0,j=0,r=0;i<NC;i++)
  {
    if (NK > 4)
    {
      assert(NK == NK);
      unsigned NK4 = (NK+4)/4;
      for(unsigned k=0;k<NK4;k++)
      {
        // muK
        _LL(r+k,i)         = mu[i];
        // Xe
        _LL(r+k,NC+j)      = -cos((M_PI*k)/(2.0*NK4));
        _LL(r+k,NC+NC+j)   = -cos((M_PI*k)/(2.0*NK4));
        // Xf
        _LL(r+k,NC+NC*2+j) = -sin((M_PI*k)/(2.0*NK4));
        _LL(r+k,NC+NC*3+j) = -sin((M_PI*k)/(2.0*NK4));
        // XeT
        _UR(NC+j,r+k)      =  cos((M_PI*k)/(2.0*NK4));
        _UR(NC+NC+j,r+k)   =  cos((M_PI*k)/(2.0*NK4));
        // XfT
        _UR(NC+NC*2+j,r+k) =  sin((M_PI*k)/(2.0*NK4));
        _UR(NC+NC*3+j,r+k) =  sin((M_PI*k)/(2.0*NK4));
      }
      r+=NK4;
      j++;
    }
    else
    {
      //      assert(NK == 1);
      // muK
      _LL(r,i)         = mu[i];
      // Xe
      _LL(r,NC+j)      = -1.0;
      _LL(r,NC+NC+j)   = -1.0;
      // Xf
      _LL(r,NC+NC*2+j) = -1.0;
      _LL(r,NC+NC*3+j) = -1.0;
      // XeT
      _UR(NC+j,r)      =  1.0;
      _UR(NC+NC+j,r)   =  1.0;
      // XfT
      _UR(NC+NC*2+j,r) =  1.0;
      _UR(NC+NC*3+j,r) =  1.0;
      r += 1;
      j++;
    }
  }
  _B.set_zero(_LL.rows()+_UR.rows(),_LL.columns()+_UR.columns());
  _B.set_sub_mat(_UR.rows(),0,_LL);
  _B.set_sub_mat(0,_LL.columns(),_UR);
  
  // Assuming that C is of full row rank (no dependent joint constraints)
  // A is invertible; then we just need to solve the LCP:
  
  // | B - D*inv(A)*C | | v | + | h - D*inv(A)*g | = | w |
  // and use the result to solve for u:
  // u = -inv(A)*(g + Cv)
  
  // A is the matrix | M X'|
  //                 | X 0 |  where X is [ S; T; P ]
  // blockwise inversion yields
  // inv(A) = [
  //    inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y ;
  //    Y*X*inv(M)                    -Y          ]
  // where Y = inv(X*inv(M)*X')
  
  // defining Q = C and using the result above yields following LCP:
  // matrix: Q*inv(A)*Q' = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  // vector: Q*inv(A)*a  = -(-Q*v + Q*inv(M)*X'*Y*X*v - Q*inv(M)*X'*Y*[0;0;vq*])
  //                     = Q*v - Q*inv(M)*X'*Y*X*v + Q*inv(M)*X'*Y*[0;0;vq*])
  
  _Y = _Jx_iM_JxT;
  
  // Invert (chol factorize) matrix _Y
  bool success = _LA.factor_chol(_Y);
  assert(success);
  
  // inv(A) = [
  //    inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y ;
  //    Y*X*inv(M)                    -Y          ]
  // defining Y = inv(X*inv(M)*X') and Q = [Cn]
  // matrix:B-D*inv(A)*C = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  //
  // vector:h-Q*inv(A)*a  = -(-Q*inv(M)*M*v - -Q*inv(M)*X'*Y*X*inv(M)*M*v + -Q*inv(M)*X'*Y*[0;0;vq*])
  //                     = -(-Q*v + Q*inv(M)*X'*Y*X*v - Q*inv(M)*X'*Y*[0;0;vq*])
  //                     =    Q*v - Q*inv(M)*X'*Y*X*v + Q*inv(M)*X'*Y*[0;0;vq*]
  
  // setup Q*inv(M)*Q'
  _MM.set_zero(NC+NC*NK+NC,NC+NC*NK+NC);
  _MM.set_sub_mat(0,0,_Cn_iM_CnT);
  _MM.set_sub_mat(NC,0,_Cd_iM_CnT);
  _MM.set_sub_mat(0,NC,_Cn_iM_CdT);
  _MM.set_sub_mat(NC,NC,_Cd_iM_CdT);
  
  // setup Q*inv(M)*X'
  // NOTE: with will break if NK != 4
  _Q_iM_XT.resize(NC+NC*NK+NC,NQ);
  _Q_iM_XT.set_sub_mat(N_IDX, J_IDX, _Cn_iM_JxT);
  _Q_iM_XT.set_sub_mat(D_IDX, J_IDX, _Cd_iM_JxT);
  _Q_iM_XT.set_sub_mat(E_IDX, J_IDX, Ravelin::MatrixNd::zero(NC,NQ));
  
  //  Ravelin::MatrixNd Y;
  //  Y = Ravelin::MatrixNd::identity(_Y.rows());
  //  _LA.solve_chol_fast(_Y, Y);
  
  // compute Y*X*inv(M)*Q'
  Ravelin::MatrixNd::transpose(_Q_iM_XT, _workM);
  _LA.solve_chol_fast(_Y, _workM);
  
  // compute Q*inv(M)*X'*Y*X*inv(M)*Q'
  _Q_iM_XT.mult(_workM, _workM2);
  _MM  -= _workM2;
  _MM += _B;
  
  // setup -Q*v
  _qq.set_zero(NC*(1+NK+1));
  _qq.set_sub_vec(0,Cn_v);
  _qq.set_sub_vec(NC,Cd_v);
  
  // setup X*v
  _Xv.resize(NQ);
  _Xv.set_sub_vec(J_IDX, Jx_v);
  
  // compute Y*X*v
  _YXv = _Xv;
  _LA.solve_chol_fast(_Y, _YXv);
  
  // compute Q*inv(M)*X' * Y*X*v
  // & subtract from LCP Vector
  _Q_iM_XT.mult(_YXv, _qq,-1.0,1.0);
  
  // compute Q*inv(M)*X'*Y*[vq*]
  _Yvqstar.set_zero(NQ);
  _Yvqstar.set_sub_vec(J_IDX, vqstar);
  _LA.solve_chol_fast(_Y, _Yvqstar);
  // & add to LCP Vector
  _Q_iM_XT.mult(_Yvqstar, _qq,1.0,1.0);
  
  // setup remainder of LCP vector
  
  static Ravelin::VectorNd _v;
  
  FILE_LOG(LOG_DYNAMICS) << "-- using: Moby::LCP::lcp_lemke" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "-- M: " << _MM << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "-- q: " << _qq << std::endl;
  
  if (!_lcp.lcp_lemke_regularized(_MM, _qq, _v,-20,4,1))
      throw std::runtime_error("Unable to solve constraint LCP!");
  FILE_LOG(LOG_DYNAMICS) << "-- z: " << _v << std::endl;
  
  Ravelin::VectorNd tau;
  
  // compute the friction forces
  // u = -inv(A)*(a + Cv)
  // u = inv(A)*(Q'*[cn])  [b/c we don't care about new velocity]
  // recalling that inv(A) =
  // | inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y | ngc x ngc,    ngc x sz(x)
  // | Y*X*inv(M)                    -Y          | sz(x) x ngc,  sz(x) x sz(x)
  // Q is nlcp x (ngc + sz(x))
  // [cs; ct] = -Y*X*v - Y*X*inv(M)*Q'*[cn]
  tau = _YXv;
  tau -= _Yvqstar;
  _Q_iM_XT.transpose_mult(_v, _workv);
  _LA.solve_chol_fast(_Y, _workv);
  tau += _workv;
  tau.negate();
  
  Ravelin::VectorNd cn,cd;
  // setup impulses
  cn = _v.segment(0, NC);
  cd = _v.segment(NC, NC*(NK+1));
  
  cf.set_zero(NC*(NK+1));
  cf.set_sub_vec(0,cn);
  cf.set_sub_vec(NC,cd);
  
  //   IDYN MLCP
  //            A             C      x         g
  //       | M −S' -T' -P'   -N' | | v+ |   |-M v |    | 0 |
  //       | S  0   0   0     0  | | cS |   |  0  |  = | 0 |
  //       | T  0   0   0     0  | | cT |   |  0  |    | 0 |
  //       | P  0   0   0     0  | | a  |   |-vq* |    | 0 |
  //                              *       +
  //       | N  0   0   0     0  | | cN |   |  0  | >= | 0 |
  //            D            B      y         h
  
  //   Constraints
  //   [M −S' -T' -N' -P'  P' ].[ ]' - M*v  = 0
  //    S v+        = 0
  //    T v+        = 0
  //    P v+ - vq*  = 0
  //    N v+       >= 0
  
  //    M v+ - (S'cs + T'cT + N'cN) - a = M v
  //    M(v+ - v) = (S'cs + T'cT + N'cN)  +  a
  
  //   Therefore
  //    M(v+ - v)           = f_total
  //    S'cs + T'cT + N'cN  = contact_force
  //    f = cf + a
  
  // Using M(dv) - z = tau
  x = tau;
  x /= dt;
  
  FILE_LOG(LOG_DYNAMICS) << "contact forces: " << cf << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "inverse dynamics forces: " << x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "<< inverse_dynamics_ap() exited" << std::endl;
  return true;
}

/*
void log_idyn_matrices(const Ravelin::VectorNd& v, const Ravelin::MatrixNd& M,const  Ravelin::MatrixNd& nT, const Ravelin::MatrixNd& D, const Ravelin::VectorNd& fext, double dt)
{
  Ravelin::MatrixNd NT = nT;
  FILE_LOG(LOG_DYNAMICS) << ">> log_idyn_matrices() entered" << std::endl;
  
  // get number of degrees of freedom and number of contact points
  unsigned N = M.rows();
  unsigned NQ = n - 6;
  unsigned NC = NT.columns();
  unsigned NK = 2;
  
  // setup R
  Ravelin::MatrixNd R(n, NC*5 ),ST(n,NC),TT(n,NC),N(NC,n),S(n,NC),T(n,NC);
  R.block(0,n,0,NC) = NT;
  R.block(0,n,NC,NC*5) = D;
  ST = D.block(0,n,0,NC);
  TT = D.block(0,n,NC,NC*2);
  N = NT;
  N.transpose();
  S = ST;
  T = TT;
  S.transpose();
  T.transpose();
  
  // compute D, E, and F
  _M_chol = M;
  bool pass = _LA.factor_chol(_M_chol);
  assert(pass);
  
  _iM.set_identity(n);
  _LA.solve_chol_fast(_M_chol,_iM);
  
  // | F E'|  =  inv(M)
  // | E D |
  Ravelin::MatrixNd E(6,NQ);
  _iM.get_sub_mat(NQ,n,0,NQ,E);
  Ravelin::MatrixNd F(NQ,NQ);
  _iM.get_sub_mat(0,NQ,0,NQ,F);
  _F_chol = F;
  pass = _LA.factor_chol(_F_chol);
  assert(pass);
  
  // compute j and k
  // [F,E']
  Ravelin::MatrixNd FET(F.rows(),F.columns()+E.rows()),
  ET = E;
  ET.transpose();
  FET.set_sub_mat(0,0,F);
  FET.set_sub_mat(0,F.columns(),ET);
  
  //////////////////////////////////////////////////////////////
  _P.set_zero(NQ,n);
  _P.set_sub_mat(0,0,Ravelin::MatrixNd::identity(NQ));
  MatrixNd::transpose(_P, _PT);
  
  Ravelin::MatrixNd Cs_iM_CsT, Cs_iM_CtT, Cs_iM_JxT,
  Ct_iM_CsT, Ct_iM_CtT, Ct_iM_JxT,
  _Cn_iM_CsT, _Cn_iM_CtT,   _Cn_iM_CnT,   _Cn_iM_JxT,
  _Jx_iM_JxT;
  
  // S
  _LA.solve_chol_fast(_M_chol,_workM = ST);
  S.mult(_workM,Cs_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = TT);
  S.mult(_workM,Cs_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM =  NT);
  //  S.mult(_workM,Cs_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  S.mult(_workM,Cs_iM_JxT);
  
  // T
  _LA.solve_chol_fast(_M_chol,_workM = ST);
  T.mult(_workM,Ct_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = TT);
  T.mult(_workM,Ct_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = NT);
  //  T.mult(_workM,Ct_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  T.mult(_workM,Ct_iM_JxT);
  
  // N
  _LA.solve_chol_fast(_M_chol,_workM = ST);
  N.mult(_workM,_Cn_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = TT);
  N.mult(_workM,_Cn_iM_CtT);
  
  _LA.solve_chol_fast(_M_chol,_workM = NT);
  N.mult(_workM,_Cn_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  N.mult(_workM,_Cn_iM_JxT);
  
  // P
  //  _LA.solve_chol_fast(_M_chol,_workM = ST);
  //  P.mult(_workM,Jx_iM_CsT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = TT);
  //  P.mult(_workM,Jx_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = NT);
  //  P.mult(_workM,Jx_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  _P.mult(_workM,_Jx_iM_JxT);
  
  Ravelin::VectorNd Cs_v, Ct_v, Cn_v, Jx_v;
  S.mult(v,Cs_v);
  T.mult(v,Ct_v);
  N.mult(v,Cn_v);
  _P.mult(v,Jx_v);
  
  FILE_LOG(LOG_DYNAMICS) << "<< log_idyn_matrices() exited" << std::endl;
}
*/

/** Maximal Dissipation Model (contact + inverse dynamics) -- Algebraic setup -- 1 stage
 *  min{z}  v+' M v+
 *  such that:
 *  v+ = v + inv(M)([ N, ST+, ST- ] z + h fext + P tau)
 *  N' v+ >= 0
 *  f_N >= 0
 *  mu * f_N{i} >= 1' [f_S{i}; f_T{i}] for i=1:NC
 *  P' v+ = v + h qdd
 **/
/*
bool RCArticulatedBodyInvDyn::inverse_dynamics_one_stage(const VectorNd& v, const VectorNd& qdd, const MatrixNd& M,const  MatrixNd& N, const MatrixNd& ST, const VectorNd& fext_, double h, const MatrixNd& MU, VectorNd& x, VectorNd& cf_final, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS)
{
  FILE_LOG(LOG_DYNAMICS) << ">> inverse_dynamics() entered" << std::endl;
  
  VectorNd fext = fext_;
  // get number of degrees of freedom and number of contact points
  unsigned NGC = M.rows();
  unsigned NQ = NGC - 6;
  unsigned NC = N.columns();
  
  //  if (NC==0) {
  //    return false;
  //  }
  
  VectorNd fID;
  
  VectorNd vq(NQ);
  v.get_sub_vec(0,NQ,vq);
  
  VectorNd vb(6);
  v.get_sub_vec(NQ,NGC,vb);
  
  VectorNd vqstar;
  ((vqstar = qdd) *= h) += vq;
  
  // Log these function variables
  
  // compute A, B, and C
  // | C B'| = M
  // | B A |
  SharedConstMatrixNd C = M.block(0,NQ,0,NQ);
  SharedConstMatrixNd B = M.block(NQ,NGC,0,NQ);
  SharedConstMatrixNd A = M.block(NQ,NGC,NQ,NGC);
  
  // compute D, E, and F
  _M_chol = M;
  bool pass = _LA.factor_chol(_M_chol);
  assert(pass);
  
  _iM.set_identity(NGC);
  _LA.solve_chol_fast(_M_chol,_iM);
  //  _LA.solve_fast(M,iM);
  
  // | F E'|  =  inv(M)
  // | E D |
  SharedMatrixNd D = _iM.block(NQ,NGC,NQ,NGC);
  SharedMatrixNd E = _iM.block(NQ,NGC,0,NQ);
  SharedMatrixNd F = _iM.block(0,NQ,0,NQ);
  _F_chol = F;
  pass = _LA.factor_chol(_F_chol);
  assert(pass);
  
  // if in mid-air only return ID forces solution
  // fID + fext = M qdd  ==> fID = M qdd - fext
  C.mult((_workv = qdd),fID) -= fext.get_sub_vec(0,NQ,_workv2);
  
  const unsigned NK = (NC == 0)? 0 : ST.columns()/NC;
  const unsigned NVARS = NC + NC*(NK);

  // setup R
  _R.resize(NGC, NC + (NC*NK) );
  _R.block(0,NGC,0,NC) = N;
  _R.block(0,NGC,NC,NC*NK+NC) = ST;
  
  // compute j and k
  // [E,D]
  MatrixNd ED(E.rows(),E.columns()+D.columns());
  ED.set_sub_mat(0,0,E);
  ED.set_sub_mat(0,E.columns(),D);
  // [F,E']
  MatrixNd FET(F.rows(),F.columns()+E.rows());
  MatrixNd::transpose(E, _ET);
  FET.set_sub_mat(0,0,F);
  FET.set_sub_mat(0,F.columns(),_ET);
  
  // Predict Contact forces
  // Incorporate fID into acting forces on robot, then find contact forces
  fext.segment(0,NQ) += fID;
  
  /// Stage 1 optimization energy minimization
  VectorNd z(NVARS),cf(NVARS);
  
  // j and k
  
  // j = [E,D]*fext*h + vb
  VectorNd j;
  // fext + [0;fID]
  ED.mult(fext,(j = vb),h,1);
  
  // k = [F,E']*fext*h  +  vq
  VectorNd k;
  FET.mult(fext,(k = vq),h,1);
  
  // compute Z and p
  // Z = ( [E,D] - E inv(F) [F,E'] ) R
  MatrixNd Z(ED.rows(), _R.columns());
  _workM = FET;
  _LA.solve_chol_fast(_F_chol,_workM);
  E.mult(_workM,_workM2);
  _workM2 -= ED;
  _workM2.mult(_R,Z,-1,0);
  
  // p = j + E inv(F) (vq* - k)
  VectorNd p = j;
  _workv = vqstar;
  _workv -= k;
  _LA.solve_chol_fast(_F_chol,_workv);
  E.mult(_workv,p,1,1);
  
  // H = Z'A Z
  MatrixNd H(Z.columns(),Z.columns());
  // compute objective function
  Z.transpose_mult(A,_workM);
  _workM.mult(Z,H);
  
  if(cf_final.rows() != 0){
    // return the inverse dynamics forces
    // x = iF*(vqstar - k - FET*R*(cf))/h
    (x = vqstar) -= k;
    FET.mult(_R,_workM);
    _workM.mult(cf_final,x,-1,1);
    _LA.solve_chol_fast(_F_chol,x);
    x /= h;
    x += fID;
    return true;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  ///////////////// Stage 1 optimization:  IDYN Energy Min ////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////// OBJECTIVE ///////////////////////////////////
  // set Hessian:
  // qG = Z'A Z = [H]
  MatrixNd qG = H;
  // set Gradient:
  // qc = Z'A p + Z'B vq*;
  VectorNd qc(Z.columns());
  // HINT: _workM = Z'*A
  
  // qc = Z'A p
  _workM.mult(p,qc);
  
  Z.transpose_mult(B,_workM2);
  // HINT: workM2 = Z'B
  
  // qc += Z'B vqstar
  _workM2.mult(vqstar,qc,1,1);
  
  ////////////////////////////// CONSTRAINTS ///////////////////////////////////
  
  // setup linear inequality constraints -- noninterpenetration
  // N'[zeros(NQ,:) ; Z] z + N'[vq* ; p] >= 0
  
  // [zeros(NQ,:) ; Z]
  _workM.set_zero(NGC,Z.columns());
  _workM.set_sub_mat(NQ,0,Z);
  // constraint Jacobain 1:
  // qM1 = N'[zeros(NQ,:) ; Z]
  MatrixNd qM1(N.columns(),Z.columns());
  N.transpose_mult(_workM,qM1);
  // [vq* ; p]
  VectorNd vqstar_p(NGC);
  vqstar_p.set_sub_vec(0,vqstar);
  vqstar_p.set_sub_vec(NQ,p);
  
  // constraint vector 1
  // qq1 = -N'[vq* ; p]
  VectorNd qq1(N.columns());
  N.transpose_mult(vqstar_p,qq1);
  qq1.negate();
  
  // setup linear inequality constraints -- coulomb friction
  // where : z = [cN  cS cT  -cS -cT]'
  // mu_i cN_i - cS_i - cT_i >= 0
  
  MatrixNd qM2;
  VectorNd qq2;
  // rhs ia zero
  
  
  // inscribe friction polygon in friction cone (scale by cos(pi/NK))
  if(NK == 4){
    qM2.set_zero(NC, NVARS);
    qq2.set_zero(NC);
    for (unsigned ii=0;ii < NC;ii++){
      // normal force
      qM2(ii,ii) = MU(ii,0);
      // tangent forces [polygonal]
      for(unsigned kk=NC+ii;kk<NC+NK*NC;kk+=NC)
        qM2(ii,kk) = -1.0;
    }
  } else {
    qM2.set_zero(NC*NK/2, NVARS);
    qq2.set_zero(NC*NK/2);
    double polygon_rad = cos(M_PI/NK);
    // for each Contact
    for (unsigned ii=0;ii < NC;ii++){
      // for each Friction Direction
      for(unsigned k=0;k<NK/2;k++){
        // normal force
        qM2(ii*NK/2+k,ii) = MU(ii,k)*polygon_rad;
        // tangent forces [polygonal]
        for(unsigned kk=NC+ii+NC*k;kk<NC+NK*NC;kk+=NC*NK/2)
          qM2(ii*NK/2+k,kk) = -1.0;
      }
    }
  }
  
  // combine all linear inequality constraints
  assert(qM1.columns() == qM2.columns());
  MatrixNd qM(qM1.rows()+qM2.rows(),qM1.columns());
  VectorNd qq(qq1.rows()+qq2.rows());
  qM.set_sub_mat(0,0,qM1);
  qM.set_sub_mat(qM1.rows(),0,qM2);
  qq.set_sub_vec(0,qq1);
  qq.set_sub_vec(qq1.rows(),qq2);
  
  static VectorNd _v;
  bool warm_start = true;
//  if(_v.rows() != (qq.rows() + z.rows()) || !SAME_AS_LAST_CONTACTS)
    warm_start = false;
  
  if(!solve_qp_pos(qG,qc,qM,qq,z,_v,warm_start)){
    FILE_LOG(LOG_DYNAMICS)  << "%ERROR: Unable to solve stage 1!";
    return false;
  }
  
  // measure feasibility of solution
  // qM z - qq >= 0
  VectorNd feas;
  qM.mult(z,feas) -= qq;
  
  
  // push z into output vector
  cf_final = z;
  // return the inverse dynamics forces
  // x = iF*(vqstar - k - FET*R*(cf))/h
  (x = vqstar) -= k;
  FET.mult(_R,_workM);
  _workM.mult(cf,x,-1,1);
  _LA.solve_chol_fast(_F_chol,x);
  x /= h;
  x += fID;
  
  // Some debugging dialogue
  FILE_LOG(LOG_DYNAMICS) << "<< inverse_dynamics() exited" << std::endl;
  return true;
}
*/

/** Maximal Dissipation Model (contact + inverse dynamics) -- Algebraic setup -- 2 stage
 *  min{z}  v+' M v+
 *  such that:
 *  v+ = v + inv(M)([ N, ST+, ST- ] z + h fext + P tau)
 *  N' v+ >= 0
 *  f_N >= 0
 *  mu * f_N{i} >= 1' [f_S{i}; f_T{i}] for i=1:nc
 *  P' v+ = v + h qdd
 **/
/*
bool RCArticulatedBodyInvDyn::inverse_dynamics_two_stage(const VectorNd& v, const VectorNd& qdd, const MatrixNd& M,const  MatrixNd& N, const MatrixNd& ST, const VectorNd& fext_, double h, const MatrixNd& MU, VectorNd& x, VectorNd& cf_final, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS)
{
  FILE_LOG(LOG_DYNAMICS) << ">> inverse_dynamics() entered" << std::endl;
  
  VectorNd fext = fext_;
  // get number of degrees of freedom and number of contact points
  unsigned NGC = M.rows();
  unsigned NQ = NGC - 6;
  unsigned NC = N.columns();
  
  //  if (NC==0) {
  //    return false;
  //  }
  
  VectorNd fID;
  
  VectorNd vq(NQ);
  v.get_sub_vec(0,NQ,vq);
  
  VectorNd vb(6);
  v.get_sub_vec(NQ,NGC,vb);
  
  VectorNd vqstar;
  ((vqstar = qdd) *= h) += vq;
  
  // Log these function variables
  
  // compute A, B, and C
  // | C B'| = M
  // | B A |
  SharedConstMatrixNd C = M.block(0,NQ,0,NQ);
  SharedConstMatrixNd B = M.block(NQ,NGC,0,NQ);
  SharedConstMatrixNd A = M.block(NQ,NGC,NQ,NGC);
  
  // compute D, E, and F
  _M_chol = M;
  bool pass = _LA.factor_chol(_M_chol);
  assert(pass);
  
  _iM.set_identity(NGC);
  _LA.solve_chol_fast(_M_chol,_iM);
  //  _LA.solve_fast(M,iM);
  
  // | F E'|  =  inv(M)
  // | E D |
  SharedConstMatrixNd D = _iM.block(NQ,NGC,NQ,NGC);
  SharedConstMatrixNd E = _iM.block(NQ,NGC,0,NQ);
  SharedConstMatrixNd F = _iM.block(0,NQ,0,NQ);
  _F_chol = F;
  pass = _LA.factor_chol(_F_chol);
  assert(pass);
  
  // if in mid-air only return ID forces solution
  // fID + fext = M qdd  ==> fID = M qdd - fext
  C.mult((_workv = qdd),fID) -= fext.get_sub_vec(0,NQ,_workv2);
  
  unsigned NK = (NC == 0)? 0 : ST.columns()/NC;
  unsigned NVARS = NC + NC*(NK);
  
  // setup R
  _R.resize(NGC, NC + (NC*NK) );
  _R.block(0,NGC,0,NC) = N;
  _R.block(0,NGC,NC,NC*NK+NC) = ST;
  
  // compute j and k
  // [E,D]
  MatrixNd ED(E.rows(),E.columns()+D.columns());
  ED.set_sub_mat(0,0,E);
  ED.set_sub_mat(0,E.columns(),D);
  // [F,E']
  MatrixNd FET(F.rows(),F.columns()+E.rows());
  MatrixNd::transpose(E, _ET);
  FET.set_sub_mat(0,0,F);
  FET.set_sub_mat(0,F.columns(),_ET);
  
  // Predict Contact forces
  // Incorporate fID into acting forces on robot, then find contact forces
  fext.segment(0,NQ) += fID;
  
  /// Stage 1 optimization energy minimization
  VectorNd z(NVARS),cf(NVARS);
  
  // j and k
  
  // j = [E,D]*fext*h + vb
  VectorNd j;
  // fext + [0;fID]
  ED.mult(fext,(j = vb),h,1);
  
  // k = [F,E']*fext*h  +  vq
  VectorNd k;
  FET.mult(fext,(k = vq),h,1);
  
  // compute Z and p
  // Z = ( [E,D] - E inv(F) [F,E'] ) R
  MatrixNd Z(ED.rows(), _R.columns());
  _workM = FET;
  _LA.solve_chol_fast(_F_chol,_workM);
  E.mult(_workM,_workM2);
  _workM2 -= ED;
  _workM2.mult(_R,Z,-1,0);
  
  // p = j + E inv(F) (vq* - k)
  VectorNd p = j;
  _workv = vqstar;
  _workv -= k;
  _LA.solve_chol_fast(_F_chol,_workv);
  E.mult(_workv,p,1,1);
  
  // H = Z'A Z
  MatrixNd H(Z.columns(),Z.columns());
  // compute objective function
  Z.transpose_mult(A,_workM);
  _workM.mult(Z,H);
  
  if(cf_final.rows() != 0){
    // return the inverse dynamics forces
    // x = iF*(vqstar - k - FET*R*(cf))/h
    (x = vqstar) -= k;
    FET.mult(_R,_workM);
    _workM.mult(cf_final,x,-1,1);
    _LA.solve_chol_fast(_F_chol,x);
    x /= h;
    x += fID;
    return true;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  ///////////////// Stage 1 optimization:  IDYN Energy Min ////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////// OBJECTIVE ///////////////////////////////////
  // set Hessian:
  // qG = Z'A Z = [H]
  MatrixNd qG = H;
  // set Gradient:
  // qc = Z'A p + Z'B vq*;
  VectorNd qc(Z.columns());
  // HINT: _workM = Z'*A
  
  // qc = Z'A p
  _workM.mult(p,qc);
  
  Z.transpose_mult(B,_workM2);
  // HINT: workM2 = Z'B
  
  // qc += Z'B vqstar
  _workM2.mult(vqstar,qc,1,1);
  
  ////////////////////////////// CONSTRAINTS ///////////////////////////////////
  
  // setup linear inequality constraints -- noninterpenetration
  // N'[zeros(NQ,:) ; Z] z + N'[vq* ; p] >= 0
  
  // [zeros(NQ,:) ; Z]
  _workM.set_zero(NGC,Z.columns());
  _workM.set_sub_mat(NQ,0,Z);
  // constraint Jacobain 1:
  // qM1 = N'[zeros(NQ,:) ; Z]
  MatrixNd qM1(N.columns(),Z.columns());
  N.transpose_mult(_workM,qM1);
  // [vq* ; p]
  VectorNd vqstar_p(NGC);
  vqstar_p.set_sub_vec(0,vqstar);
  vqstar_p.set_sub_vec(NQ,p);
  
  // constraint vector 1
  // qq1 = -N'[vq* ; p]
  VectorNd qq1(N.columns());
  N.transpose_mult(vqstar_p,qq1);
  qq1.negate();
  
  // setup linear inequality constraints -- coulomb friction
  // where : z = [cN  cS cT  -cS -cT]'
  // mu_i cN_i - cS_i - cT_i >= 0
  
  MatrixNd qM2;
  VectorNd qq2;
  // rhs ia zero
  
  
  // inscribe friction polygon in friction cone (scale by cos(pi/NK))
  if(NK == 4){
    qM2.set_zero(NC, NVARS);
    qq2.set_zero(NC);
    for (unsigned ii=0;ii < NC;ii++){
      // normal force
      qM2(ii,ii) = MU(ii,0);
      // tangent forces [polygonal]
      for(unsigned kk=NC+ii;kk<NC+NK*NC;kk+=NC)
        qM2(ii,kk) = -1.0;
    }
  } else {
    qM2.set_zero(NC*NK/2, NVARS);
    qq2.set_zero(NC*NK/2);
    double polygon_rad = cos(M_PI/NK);
    // for each Contact
    for (unsigned ii=0;ii < NC;ii++){
      // for each Friction Direction
      for(unsigned k=0;k<NK/2;k++){
        // normal force
        qM2(ii*NK/2+k,ii) = MU(ii,k)*polygon_rad;
        // tangent forces [polygonal]
        for(unsigned kk=NC+ii+NC*k;kk<NC+NK*NC;kk+=NC*NK/2)
          qM2(ii*NK/2+k,kk) = -1.0;
      }
    }
  }
  
  // combine all linear inequality constraints
  assert(qM1.columns() == qM2.columns());
  MatrixNd qM(qM1.rows()+qM2.rows(),qM1.columns());
  VectorNd qq(qq1.rows()+qq2.rows());
  qM.set_sub_mat(0,0,qM1);
  qM.set_sub_mat(qM1.rows(),0,qM2);
  qq.set_sub_vec(0,qq1);
  qq.set_sub_vec(qq1.rows(),qq2);
  
  static VectorNd _v;
  bool warm_start = true;
//  if(_v.rows() != (qq.rows() + z.rows()) || !SAME_AS_LAST_CONTACTS)
    warm_start = false;
  
  if(!solve_qp_pos(qG,qc,qM,qq,z,_v,warm_start)){
    FILE_LOG(LOG_DYNAMICS)  << "%ERROR: Unable to solve stage 1!";
    return false;
  }
  
  // measure feasibility of solution
  // qM z - qq >= 0
  VectorNd feas;
  qM.mult(z,feas) -= qq;
  
  
  // push z into output vector
  cf = z;
  
  /////////////////////////////////////////////////////////////////////////////
  ///////////////// Stage 2 optimization: command smoothing ///////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  // H = Z'AZ
  MatrixNd P;
  _LA.nullspace(H,P);
  unsigned size_null_space = P.columns();
  if(size_null_space != 0)
  {
    // second optimization is necessary if the previous Hessian was PSD:
    // size_null_space > 0
    
    // compute U = [F,E']*R
    MatrixNd U;
    FET.mult(_R,U);
    
    /////////////////////////////// OBJECTIVE //////////////////////////////////
    
    // Objective Hessian:
    // qG = P'*U'*iF'*iF*U*P;
    _workM = U;
    _LA.solve_chol_fast(_F_chol,_workM);
    _workM.mult(P,_workM2);
    _workM2.transpose_mult(_workM2,qG);
    
    // HINT: workM2 = iF*U*P
    // HINT: _workM = iF*U
    
    // Objective Gradient:
    // qc = z'*U'*iF'*iF*U*P - vqstar'*iF'*iF*U*P + k'*iF'*iF*U*P;
    
    // qc = (iF*U*P)'*iF*U*z
    _workM2.transpose_mult(_workM.mult(z,_workv),_workv2);
    qc = _workv2;
    
    // qc -= (iF*U*P)'* iF*vqstar
    _workv = vqstar;
    _workM2.transpose_mult(_LA.solve_chol_fast(_F_chol,_workv),_workv2);
    qc -= _workv2;
    
    // qc += (iF*U*P)'* iF*k
    _workv = k;
    _workM2.transpose_mult(_LA.solve_chol_fast(_F_chol,_workv),_workv2);
    qc += _workv2;
    
    ////////////////////////////// CONSTRAINTS /////////////////////////////////
    
    // Linear Inequality Constraints:
    
    // Compressive force constraint (& polygonal tangential forces):
    // z + Pw >= 0 || P*w >= -z
    
    // Constraint Jacobian 1:
    // qM1 = null(H) = P
    qM1 = P;
    
    // Constraint Vector 1:
    // qq1 = z (contact impulses from Stage 1)
    qq1 = z;
    qq1.negate();
    
    // Non-Interpenetration:
    // SRZ: P = null( Z'*H*Z ) --> P = null(Z) this means:
    //       noninterpenetration & linear energy constraints always = 0
    
    // qM2 = N'*[zeros(NQ,NVARS_null);Z*P];
    // qq2 = N'*[Z*z + p];
    
    // Coulomb Friction Polygon:
    NVARS = P.columns();
    MatrixNd qM3;
    VectorNd qq3;

    if(NK == 4){
      qM3.set_zero(NC, NVARS);
      qq3.set_zero(NC);
      for (unsigned ii=0;ii < NC;ii++){
        // normal direction
        //  qM3(ii,:) = P(ii,:)
        //  qq3(ii) = -z(ii)
        qM3.row(ii) = ((_workv = P.row(ii))*=MU(ii,0));
        qq3[ii] = -z[ii]*MU(ii,0);
        
        // tangent directions
        // kk indexes matrix, k refers to contact direction
        for(unsigned kk=NC+ii;kk<NC+NK*NC;kk+=NC){
          qM3.row(ii) -= P.row(kk);
          qq3[ii]     += z[kk];
        }
      }
    } else {
      qM3.set_zero(NC*NK/2, NVARS);
      qq3.set_zero(NC*NK/2);
      for (unsigned ii=0;ii < NC;ii++){
        // for each Friction Direction
        for(unsigned k=0;k<NK/2;k++){
          // normal force
          qM3.row(ii*NK/2+k) = ((_workv = P.row(ii))*=MU(ii,k));
          qq3[ii*NK/2+k] = -z[ii]*MU(ii,k);
          // tangent forces [polygonal]
          for(unsigned kk=NC+ii+NC*k;kk<NC+NK*NC;kk+=NC*NK/2){
            qM3.row(ii*NK/2+k) -= P.row(kk);
            qq3[ii*NK/2+k]     += z[kk];
          }
        }
      }
    }
    
    // Set up constarint matrix
    // SRZ: constraint 2 (interpenetration) is not here
    qM.set_zero(qM1.rows()+qM3.rows(),qM1.columns());
    qq.set_zero(qM1.rows()+qM3.rows());
    qM.set_sub_mat(0,0,qM1);
    qM.set_sub_mat(qM1.rows(),0,qM3);
    qq.set_sub_vec(0,qq1);
    qq.set_sub_vec(qq1.rows(),qq3);
    
    // optimize system
    VectorNd w(size_null_space);
    if(!solve_qp(qG,qc,qM,qq,w)){
      FILE_LOG(LOG_DYNAMICS)  << "ERROR: Unable to solve stage 2!";
      return false;
      // then skip to calculating x from stage 1 solution
    } else {
      
      // measure feasibility of solution
      // qM w - qq >= 0
      feas = qq;
      qM.mult(w,feas,1,-1);
      
      // return the solution (contact forces)
      // cf = z + P*w;
      
      P.mult(w,cf,1,1);
      
    }
  }
  
  //  Note compare contact force prediction to Moby contact force
  
  cf_final = cf;
  // return the inverse dynamics forces
  // x = iF*(vqstar - k - FET*R*(cf))/h
  (x = vqstar) -= k;
  FET.mult(_R,_workM);
  _workM.mult(cf,x,-1,1);
  _LA.solve_chol_fast(_F_chol,x);
  x /= h;
  x += fID;
  
  // Some debugging dialogue
  FILE_LOG(LOG_DYNAMICS) << "<< inverse_dynamics() exited" << std::endl;
  return true;
}
*/

bool RCArticulatedBodyInvDyn::inverse_dynamics_no_slip_fast(const VectorNd& vel, const VectorNd& qdd, const MatrixNd& M,const  MatrixNd& nT, const MatrixNd& D, const VectorNd& fext, double dt, VectorNd& x, VectorNd& cf, bool frictionless, std::vector<unsigned>& indices, unsigned active_eefs,bool SAME_AS_LAST_CONTACTS)
{
  const MatrixNd& NT = nT;
  FILE_LOG(LOG_DYNAMICS) << ">> inverse_dynamics_no_slip_fast() entered" << std::endl;
  
  // get number of degrees of freedom and number of contact points
  unsigned NGC = M.rows();
  unsigned NQ = NGC - 6;
  unsigned NC = NT.columns();
  unsigned NK = 2;
  
  // setup R
  _R.resize(NGC, NC*5);
  _ST.resize(NGC,NC);
  _TT.resize(NGC,NC),
  _N.resize(NC,NGC);
  _S.resize(NGC,NC);
  _T.resize(NGC,NC);
  _R.block(0,NGC,0,NC) = NT;
  _R.block(0,NGC,NC,NC*5) = D;
  _ST = D.block(0,NGC,0,NC);
  _TT = D.block(0,NGC,NC,NC*2);
  MatrixNd::transpose(NT, _N);
  MatrixNd::transpose(_ST, _S);
  MatrixNd::transpose(_TT, _T);
  
  // compute D, E, and F
  _M_chol = M;
  bool pass = _LA.factor_chol(_M_chol);
  assert(pass);
  
  _iM.set_identity(NGC);
  _LA.solve_chol_fast(_M_chol,_iM);
  
  // | F E'|  =  inv(M)
  // | E D |
  SharedMatrixNd E = _iM.block(NQ,NGC,0,NQ);
  SharedMatrixNd F = _iM.block(0,NQ,0,NQ);
  _F_chol = F;
  pass = _LA.factor_chol(_F_chol);
  assert(pass);
  
  // compute j and k
  // [F,E']
  MatrixNd FET(F.rows(),F.columns()+E.rows());
  MatrixNd::transpose(E, _ET);
  FET.set_sub_mat(0,0,F);
  FET.set_sub_mat(0,F.columns(),_ET);
  
  SharedConstVectorNd vq = vel.segment(0,NQ);
  VectorNd v = vel;
  _iM.mult(fext,v,dt,1);
  
  VectorNd vqstar;
  ((vqstar = qdd) *= dt) += v.segment(0,NQ);
  
  //////////////////////////////////////////////////////////////
  _P.set_zero(NQ,NGC);
  _P.block(0,0,NQ,NQ).set_identity();
  MatrixNd::transpose(_P, _PT);
  
  MatrixNd Cs_iM_CsT, Cs_iM_CtT, Cs_iM_JxT,
  Ct_iM_CsT, Ct_iM_CtT, Ct_iM_JxT,
  _Cn_iM_CsT, _Cn_iM_CtT,   _Cn_iM_CnT,   _Cn_iM_JxT,
  _Jx_iM_JxT;
  
  // S
  _LA.solve_chol_fast(_M_chol,_workM = _ST);
  _S.mult(_workM,Cs_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _TT);
  _S.mult(_workM,Cs_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM =  NT);
  //  S.mult(_workM,Cs_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  _S.mult(_workM,Cs_iM_JxT);
  
  // T
  _LA.solve_chol_fast(_M_chol,_workM = _ST);
  _T.mult(_workM,Ct_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _TT);
  _T.mult(_workM,Ct_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = NT);
  //  T.mult(_workM,Ct_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  _T.mult(_workM,Ct_iM_JxT);
  
  // N
  _LA.solve_chol_fast(_M_chol,_workM = _ST);
  _N.mult(_workM,_Cn_iM_CsT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _TT);
  _N.mult(_workM,_Cn_iM_CtT);
  
  _LA.solve_chol_fast(_M_chol,_workM = NT);
  _N.mult(_workM,_Cn_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  _N.mult(_workM,_Cn_iM_JxT);
  
  // P
  //  _LA.solve_chol_fast(_M_chol,_workM = ST);
  //  P.mult(_workM,Jx_iM_CsT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = TT);
  //  P.mult(_workM,Jx_iM_CtT);
  
  //  _LA.solve_chol_fast(_M_chol,_workM = NT);
  //  P.mult(_workM,Jx_iM_CnT);
  
  _LA.solve_chol_fast(_M_chol,_workM = _PT);
  _P.mult(_workM,_Jx_iM_JxT);
  
  VectorNd Cs_v, Ct_v, Cn_v, Jx_v;
  _S.mult(v,Cs_v);
  _T.mult(v,Ct_v);
  _N.mult(v,Cn_v);
  _P.mult(v,Jx_v);
  
  /////////////////////////////////////////////////////////////////////////////
  /// Solve System
  ///
  std::vector<unsigned> S_indices,T_indices,J_indices;
  const unsigned N_IDX = 0;
  
  // we do this by solving the MLCP:
  // |  A  C  | | x | + | g | = | 0 |
  // |  D  B  | | y |   | h |   | 0 |
  
  // g = [-M*v'; 0; 0; -vq_err]
  // h = [ 0 ];
  // A = [ M -S' -T' -P';
  //       S  0   0   0 ;
  //       T  0   0   0 ;
  //       P  0   0   0 ];
  
  // B = [ 0 ];
  
  // C = [ -N';
  //        0 ;
  //        0 ];
  // D = -C'
  // x = [ v^+; cs; ct; a ]
  // y = [ cN]
  
  // Assuming that C is of full row rank (no dependent joint constraints)
  // A is invertible; then we just need to solve the LCP:
  
  // | B - D*inv(A)*C | | v | + | h - D*inv(A)*g | = | w |
  // and use the result to solve for u:
  // u = -inv(A)*(g + Cv)
  
  // A is the matrix | M X'|
  //                 | X 0 |  where X is [ S; T; P ]
  // blockwise inversion yields
  // inv(A) = [
  //    inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y ;
  //    Y*X*inv(M)                    -Y          ]
  // where Y = inv(X*inv(M)*X')
  
  // defining Q = C and using the result above yields following LCP:
  // matrix: Q*inv(A)*Q' = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  // vector: Q*inv(A)*a  = -(-Q*v + Q*inv(M)*X'*Y*X*v - Q*inv(M)*X'*Y*[0;0;vq*])
  //                     = Q*v - Q*inv(M)*X'*Y*X*v + Q*inv(M)*X'*Y*[0;0;vq*])
  
  // ********************************************************
  // find largest non-singular set of J, S, and T indices
  // ********************************************************
  
  // loop through joint constraints, forming J*inv(M)*J' and checking condition
  for (unsigned i=0; i< NQ; i++)
  {
    // add the index tentatively to the set
    J_indices.push_back(i);
    
    // This is an [I,0] matrix, it has full row rank
    continue;
  }
  
  // get the reduced Jx*iM*Jx' matrix
  _rJx_iM_JxT = _Jx_iM_JxT;
  
  // loop through contacts, forming matrix below and checking its condition
  // | S*inv(M)*S'  S*inv(M)*T' S*inv(M)*J' |
  // | T*inv(M)*S'  T*inv(M)*T' T*inv(M)*J' |
  // | J*inv(M)*S'  J*inv(M)*T' J*inv(M)*J' |
  
  bool last_success = false;
  for (unsigned i=0; i< NC; i++)
  {
    // update S indices
    S_indices.push_back(i);
    
    // setup indices
    unsigned S_IDX;
    S_IDX = 0;
    unsigned T_IDX;
    T_IDX = S_indices.size();
    unsigned J_IDX;
    J_IDX = T_IDX + T_indices.size();
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());
    
    // add S/S, T/T, J/J components to 'check' matrix
    Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, S_IDX, _workM);
    Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, T_IDX, _workM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);
    
    // add S/T components to 'check' matrix
    Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, T_IDX, _workM);
    _Y.set_sub_mat(T_IDX, S_IDX, _workM, eTranspose);
    
    // add S/J components to check matrix
    Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, S_IDX, _workM, eTranspose);
    
    // add T/J components to check matrix
    Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, T_IDX, _workM, eTranspose);
    
    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< _Y.rows(); j++)
      _Y(j,j) -= NEAR_ZERO;
    
    // see whether check matrix can be Cholesky factorized
    if (!_LA.factor_chol(_Y) || frictionless)
      S_indices.pop_back();
    
    // add index for T
    T_indices.push_back(i);
    
    // resize the check matrix
    T_IDX = S_indices.size();
    J_IDX = T_IDX + T_indices.size();
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());
    
    // add S/S, T/T, J/J components to 'check' matrix
    Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, S_IDX, _workM);
    Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, T_IDX, _workM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);
    
    // add S/T components to 'check' matrix
    Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, T_IDX, _workM);
    _Y.set_sub_mat(T_IDX, S_IDX, _workM, eTranspose);
    
    // add S/J components to check matrix
    Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, S_IDX, _workM, eTranspose);
    
    // add T/J components to check matrix
    Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, T_IDX, _workM, eTranspose);
    
    // skew the matrix away from positive definiteness
    for (unsigned j=0; j< _Y.rows(); j++)
      _Y(j,j) -= NEAR_ZERO;
    
    // see whether check matrix can be Cholesky factorized
    last_success = _LA.factor_chol(_Y);
    if (!last_success && !frictionless)
      T_indices.pop_back();
  }
  
  last_success = false;
  // output indices
  if (true)
  {
    std::ostringstream oss;
    oss << "s indices:";
    for (unsigned i=0; i< S_indices.size(); i++)
      oss << " " << S_indices[i];
    oss << "  t indices:";
    for (unsigned i=0; i< T_indices.size(); i++)
      oss << " " << T_indices[i];
    oss << "  j indices:";
    for (unsigned i=0; i< J_indices.size(); i++)
      oss << " " << J_indices[i];
    
    FILE_LOG(LOG_DYNAMICS) << oss.str() << std::endl;
  }
  
  // ********************************************************
  // reform Y if necessary
  // ********************************************************
  
  // setup indices
  unsigned S_IDX;
  S_IDX = 0;
  unsigned T_IDX;
  T_IDX = S_indices.size();
  unsigned J_IDX;
  J_IDX = T_IDX + T_indices.size();
  
  if (!last_success)
  {
    _Y.resize(J_IDX + J_indices.size(), J_IDX + J_indices.size());
    
    // add S/S, T/T, J/J components to X
    Cs_iM_CsT.select_square(S_indices.begin(), S_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, S_IDX, _workM);
    Ct_iM_CtT.select_square(T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, T_IDX, _workM);
    _Y.set_sub_mat(J_IDX, J_IDX, _rJx_iM_JxT);
    
    // add S/T components to X
    Cs_iM_CtT.select(S_indices.begin(), S_indices.end(), T_indices.begin(), T_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, T_IDX, _workM);
    _Y.set_sub_mat(T_IDX, S_IDX, _workM, eTranspose);
    
    // add S/J components to X
    Cs_iM_JxT.select(S_indices.begin(), S_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(S_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, S_IDX, _workM, eTranspose);
    
    // add T/J components to X
    Ct_iM_JxT.select(T_indices.begin(), T_indices.end(), J_indices.begin(), J_indices.end(), _workM);
    _Y.set_sub_mat(T_IDX, J_IDX, _workM);
    _Y.set_sub_mat(J_IDX, T_IDX, _workM, eTranspose);
    
  }
  
  // do the Cholesky factorization (should not fail)
  
  // check the condition number on Y
  // TODO: remove this code when satisfied Cholesky factorization is not problem
  MatrixNd tmp = _Y;
  double cond = _LA.cond(tmp);
  if (cond > 1e6){
    FILE_LOG(LOG_DYNAMICS) << "Condition number *may* be high (check!): " << cond << std::endl;
    bool success = _LA.factor_chol(_Y);
    assert(success);
  } else {
    bool success = _LA.factor_chol(_Y);
    assert(success);
  }
  
  // inv(A) = [
  //    inv(M)-inv(M)*X'*Y*X*inv(M)   inv(M)*X'*Y ;
  //    Y*X*inv(M)                    -Y          ]
  // defining Y = inv(X*inv(M)*X') and Q = [Cn]
  // matrix:-D*inv(A)*C = Q*inv(M)*Q' - Q*inv(M)*X'*Y*X*inv(M)*Q'
  //
  // vector:-Q*inv(A)*a  = -(-Q*inv(M)*M*v - -Q*inv(M)*X'*Y*X*inv(M)*M*v + -Q*inv(M)*X'*Y*[0;0;vq*])
  //                     = -(-Q*v + Q*inv(M)*X'*Y*X*v - Q*inv(M)*X'*Y*[0;0;vq*])
  //                     =    Q*v - Q*inv(M)*X'*Y*X*v + Q*inv(M)*X'*Y*[0;0;vq*]
  
  // setup Q*inv(M)*Q'
  //  _MM = _Cn_iM_CnT;
  
  // setup Q*inv(M)*X'
  _Q_iM_XT.resize(NC, S_indices.size() + T_indices.size() + J_indices.size());
  _Cn_iM_CsT.select_columns(S_indices.begin(), S_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, S_IDX, _workM);
  _Cn_iM_CtT.select_columns(T_indices.begin(), T_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, T_IDX, _workM);
  _Cn_iM_JxT.select_columns(J_indices.begin(), J_indices.end(), _workM);
  _Q_iM_XT.set_sub_mat(N_IDX, J_IDX, _workM);
  
  //  MatrixNd Y;
  //  Y = MatrixNd::identity(_Y.rows());
  //  _LA.solve_chol_fast(_Y, Y);
  
  // compute Y*X*inv(M)*Q'
  MatrixNd::transpose(_Q_iM_XT, _workM);
  _LA.solve_chol_fast(_Y, _workM);
  //  Y.mult_transpose(_Q_iM_XT,_workM);
  
  // compute Q*inv(M)*X'*Y*X*inv(M)*Q'
  _Q_iM_XT.mult(_workM, _workM2);
  (_MM = _Cn_iM_CnT) -= _workM2;
  
  // setup -Q*v
  _qq = Cn_v;
  
  // setup X*v
  VectorNd _Xv,_YXv;
  _Xv.resize(S_indices.size() + T_indices.size() + J_indices.size());
  Cs_v.select(S_indices.begin(), S_indices.end(), _workv);
  _Xv.set_sub_vec(S_IDX, _workv);
  Ct_v.select(T_indices.begin(), T_indices.end(), _workv);
  _Xv.set_sub_vec(T_IDX, _workv);
  Jx_v.select(J_indices.begin(), J_indices.end(), _workv);
  _Xv.set_sub_vec(J_IDX, _workv);
  
  // compute Y*X*v
  _YXv = _Xv;
  _LA.solve_chol_fast(_Y, _YXv);
  
  // compute Q*inv(M)*X' * Y*X*v
  _Q_iM_XT.mult(_YXv, _workv);
  
  // Add to LCP Vector
  _qq -= _workv;
  
  // compute Q*inv(M)*X'*Y*[0;0;vq*]
  _Yvqstar.set_zero(S_indices.size() + T_indices.size() + J_indices.size());
  
  vqstar.select(J_indices.begin(), J_indices.end(), _workv);
  _Yvqstar.set_sub_vec(J_IDX, _workv);
  _LA.solve_chol_fast(_Y, _Yvqstar);
  _Q_iM_XT.mult(_Yvqstar,_workv);
  _qq += _workv;
  
  // setup remainder of LCP vector
  static VectorNd _v;
  if(_v.size() != _qq.size() || !SAME_AS_LAST_CONTACTS)
    _v.resize(0);
  
  if (active_eefs > 0) {
    // attempt to solve the LCP using the fast method
    if(active_eefs > 2){
      
      if (!lcp_fast(_MM, _qq,indices, _v,NEAR_ZERO))
      {
        FILE_LOG(LOG_DYNAMICS) << "-- Principal pivoting method LCP solver failed; falling back to regularized lemke solver" << std::endl;
        
        if (!_lcp.lcp_lemke_regularized(_MM, _qq, _v,-20,4,1))
          throw std::runtime_error("Unable to solve constraint LCP!");
      }
      
    } else {
      FILE_LOG(LOG_DYNAMICS) << "-- using: Moby::LCP::lcp_fast" << std::endl;
      
      if (!_lcp.lcp_fast(_MM, _qq, _v))
      {
        FILE_LOG(LOG_DYNAMICS) << "Principal pivoting method LCP solver failed; falling back to regularized lemke solver" << std::endl;
        
        if (!_lcp.lcp_lemke_regularized(_MM, _qq, _v,-20,4,1))
          throw std::runtime_error("Unable to solve constraint LCP!");
      }
    }
  }
  
  // u = inv(A)*(Q'*[cn])  [b/c we don't care about new velocity]
  
  _cs_ct_tau = _YXv;
  _cs_ct_tau -= _Yvqstar;

  if (active_eefs > 0) {
    _Q_iM_XT.transpose_mult(_v, _workv);
    _LA.solve_chol_fast(_Y, _workv);
    _cs_ct_tau += _workv;
  }
  _cs_ct_tau.negate();
  
  VectorNd cn,cs,ct,tau;
  // setup impulses
  cn = _v.segment(0, NC);
  cs.set_zero(NC);
  ct.set_zero(NC);
  tau.set_zero(NQ);
  SharedConstVectorNd cs_vec = _cs_ct_tau.segment(S_IDX, T_IDX);
  SharedConstVectorNd ct_vec = _cs_ct_tau.segment(T_IDX, J_IDX);
  SharedConstVectorNd a_vec = _cs_ct_tau.segment(J_IDX, _cs_ct_tau.rows());
  cs.set(S_indices.begin(), S_indices.end(), cs_vec);
  ct.set(T_indices.begin(), T_indices.end(), ct_vec);
  tau.set(J_indices.begin(), J_indices.end(), a_vec);
  
  cf.set_zero(NC*5);
  for(unsigned i=0;i< NC;i++){
    cf[i] = cn[i];
    if(ct[i] >= 0)
      cf[NC+i] = cs[i];
    else
      cf[NC+i+NC*2] = -cs[i];
    if(ct[i] >= 0)
      cf[NC+i+NC] = ct[i];
    else
      cf[NC+i+NC*3] = -ct[i];
  }
  
  //   IDYN MLCP
  //            A             C      x         g
  //       | M −S' -T' -P'   -N' | | v+ |   |-M v |    | 0 |
  //       | S  0   0   0     0  | | cS |   |  0  |  = | 0 |
  //       | T  0   0   0     0  | | cT |   |  0  |    | 0 |
  //       | P  0   0   0     0  | | a  |   |-vq* |    | 0 |
  //                              *       +
  //       | N  0   0   0     0  | | cN |   |  0  | >= | 0 |
  //            D            B      y         h
  
  //   Constraints
  //   [M −S' -T' -N' -P'  P' ].[ ]' - M*v  = 0
  //    S v+        = 0
  //    T v+        = 0
  //    P v+ - vq*  = 0
  //    N v+       >= 0
  
  //    M v+ - (S'cs + T'cT + N'cN) - a = M v
  //    M(v+ - v) = (S'cs + T'cT + N'cN)  +  a
  
  //   Therefore
  //    M(v+ - v)           = f_total
  //    S'cs + T'cT + N'cN  = contact_force
  //    f = cf + a
  
  // Using M(dv) - z = tau
  x = tau;
  x /= dt;
  
  FILE_LOG(LOG_DYNAMICS) << "<< inverse_dynamics_no_slip_fast() exited" << std::endl;
  return true;
}

bool RCArticulatedBodyInvDyn::inverse_dynamics_no_contact(const VectorNd& vel, const VectorNd& qdd_des, const MatrixNd& M, const VectorNd& fext, double dt, VectorNd& x)
{
  MatrixNd nT(vel.rows(),0),D(vel.rows(),0);
  std::vector<unsigned> indices;
  VectorNd cf;
  inverse_dynamics_no_slip_fast(vel,qdd_des,M,nT,D,fext, dt,x, cf, false, indices, 0, false);
  return true;
}

bool RCArticulatedBodyInvDyn::solve_qp_pos(const MatrixNd& Q, const VectorNd& c, VectorNd& x, bool warm_start)
{
  const int n = Q.rows();
  const int m = 0;
  
  // solve the LCP
  x.set_zero();
  bool SOLVE_FLAG = true;
  
 
  double zero_tol = Q.norm_inf()*Q.rows()*std::numeric_limits<double>::epsilon() * 1e4;
  if (warm_start) {
    if(!_lcp.lcp_fast(Q,c,x)){
      if(!_lcp.lcp_lemke_regularized(Q,c,x,-20,4,0,-1.0,zero_tol))
        SOLVE_FLAG = false;
      else
        SOLVE_FLAG = std::isfinite(x.norm());
    } else {
      SOLVE_FLAG = false;
    }
  } else {
    if(!_lcp.lcp_lemke_regularized(Q,c,x,-20,4,0,-1.0,zero_tol))
      SOLVE_FLAG = false;
    else
      SOLVE_FLAG = std::isfinite(x.norm());
  }
  
#ifndef NDEBUG
  FILE_LOG(LOG_DYNAMICS)  << "Solutions" ;
  FILE_LOG(LOG_DYNAMICS)  << " << solve qp positive" ;
#endif
  return SOLVE_FLAG;
}

bool RCArticulatedBodyInvDyn::solve_qp(const MatrixNd& Q, const VectorNd& c, const MatrixNd& A, const VectorNd& b, VectorNd& x)
{
  const int n = Q.rows();
  const int m = A.rows();
  
  // setup the LCP matrix
  // MMM = |  Q -Q -A' |
  //       | -Q  Q  A' |
  //       |  A -A  0  |
  Ravelin::MatrixNd MMM,AT = A, nA = A, nQ = Q;
  MMM.set_zero(Q.rows()*2 + A.rows(), Q.rows()*2 + A.rows());
  
  nA.negate();
  AT.transpose();
  nQ.negate();
  
  Ravelin::VectorNd zzz,qqq, nc=c, nb=b;
  qqq.set_zero(MMM.rows());
  zzz.set_zero(MMM.rows());
  
  nc.negate();
  nb.negate();
  
  MMM.set_sub_mat(0,0,Q);
  MMM.set_sub_mat(n,n,Q);
  MMM.set_sub_mat(0,n,nQ);
  MMM.set_sub_mat(n,0,nQ);
  
  // setup linear inequality constraints in LCP matrix
  MMM.set_sub_mat(n,n*2,AT);
  AT.negate();
  MMM.set_sub_mat(0,n*2,AT);
  
  MMM.set_sub_mat(n*2,0,A);
  MMM.set_sub_mat(n*2,n,nA);
  
  // setup LCP vector qqq = [c;-c;-b]
  qqq.set_sub_vec(0,c);
  qqq.set_sub_vec(n,nc);
  qqq.set_sub_vec(2*n,nb);
  
#ifndef NDEBUG
  //  FILE_LOG(LOG_DYNAMICS)  << "% >> solve qp" << std::endl;
  //  FILE_LOG(LOG_DYNAMICS)  << "%QP variables" << std::endl;
  //  OUTLOG(Q,"Q",logDEBUG1);
  //  OUTLOG(c,"c",logDEBUG1);
  //  OUTLOG(A,"AA",logDEBUG1);
  //  OUTLOG(b,"bb",logDEBUG1);
  //  FILE_LOG(LOG_DYNAMICS)  << "LCP variables" << std::endl;
  //  OUTLOG(MMM,"MMM",logDEBUG1);
  //  OUTLOG(qqq,"qqq",logDEBUG1);
#endif
  
  // solve the LCP
  bool SOLVE_FLAG = true;
#ifndef SPLITTING_METHOD
  double zero_tol = MMM.norm_inf()*MMM.rows()*std::numeric_limits<double>::epsilon() * 1e4;
  if(!_lcp.lcp_lemke_regularized(MMM,qqq,zzz,-20,4,0,-1.0,zero_tol))
    //  if(!_lcp.lcp_lemke(MMM,qqq,zzz))
    SOLVE_FLAG = false;
  else
    SOLVE_FLAG = std::isfinite(zzz.norm());
#else
  _lcpsymm_iter(MMM, qqq, zzz, 0.5, 1.0, MAX_ITER);
#endif
  
  // extract x
  for(int i=0;i<n;i++)
    x[i] = zzz[i] - zzz[n+i];
  
#ifndef NDEBUG
  //  FILE_LOG(LOG_DYNAMICS)  << "%Solutions" << std::endl;
  //  OUTLOG(zzz,"zzz",logDEBUG1);
  //  OUTLOG(x,"xx",logDEBUG1);
  //  FILE_LOG(LOG_DYNAMICS)  << "% << solve qp" << std::endl;
#endif
  return SOLVE_FLAG;
}

bool RCArticulatedBodyInvDyn::solve_qp_pos(const MatrixNd& Q, const VectorNd& c, const MatrixNd& A, const VectorNd& b, VectorNd& x, VectorNd& v, bool warm_start, bool regularize)
{
  const int n = Q.rows();
  const int m = A.rows();
  
  Ravelin::MatrixNd MMM;
  MMM.set_zero(n + m,n + m);
  
  Ravelin::VectorNd zzz(n + m),qqq(n + m);
  qqq.set_zero();
  zzz.set_zero();
  // init and setup MMM
  Ravelin::MatrixNd nAT = A;
  nAT.transpose();
  nAT.negate();
  MMM.set_sub_mat(0,0,Q);
  MMM.set_sub_mat(0,n,nAT);
  MMM.set_sub_mat(n,0,A);
  
  // setup qqq
  qqq.set_sub_vec(0,c);
  Ravelin::VectorNd nb = b;
  nb.negate();
  qqq.set_sub_vec(n,nb);
  
  // solve the LCP
  bool SOLVE_FLAG = true;
  
#ifndef NDEBUG
  FILE_LOG(LOG_DYNAMICS)  << "% >> solve qp positive" ;
  FILE_LOG(LOG_DYNAMICS)  << "%QP variables" ;
  //OUT_LOG(logINFO)  << "LCP variables" << std::endl;
  //OUTLOG(MMM,"MM",logDEBUG1);
  //OUTLOG(qqq,"qq",logDEBUG1);
#endif
  
#ifndef SPLITTING_METHOD
  double zero_tol = MMM.norm_inf()*MMM.rows()*std::numeric_limits<double>::epsilon() * 1e4;
  if (warm_start) {
    zzz = v;
    if(!_lcp.lcp_fast(MMM,qqq,zzz)){
      if(regularize){
        if(!_lcp.lcp_lemke_regularized(MMM,qqq,zzz,-20,4,0,-1.0,zero_tol))
          SOLVE_FLAG = false;
        else
          SOLVE_FLAG = std::isfinite(zzz.norm());
      } else {
        if(!_lcp.lcp_lemke(MMM,qqq,zzz,-1.0,zero_tol))
          SOLVE_FLAG = false;
      }
    } else {
      SOLVE_FLAG = false;
    }
  } else {
    if(regularize){
      if(!_lcp.lcp_lemke_regularized(MMM,qqq,zzz,-20,4,0,-1.0,zero_tol))
        SOLVE_FLAG = false;
      else
        SOLVE_FLAG = std::isfinite(zzz.norm());
    } else {
      if(!_lcp.lcp_lemke(MMM,qqq,zzz,-1.0,zero_tol))
        SOLVE_FLAG = false;
    }
  }
  
  
  if (SOLVE_FLAG) {
    v = zzz;
  }
#endif
  // extract x
  for(int i=0;i<n;i++)
    x[i] = zzz[i];
#ifndef NDEBUG
  FILE_LOG(LOG_DYNAMICS)  << "%Solutions" ;
  //  OUTLOG(zzz,"zz");
  FILE_LOG(LOG_DYNAMICS)  << "% << solve qp positive" ;
#endif
  return SOLVE_FLAG;
}

/// Performs Cholesky factorizations on a group of generalized inertia matrices
bool RCArticulatedBodyInvDyn::factor_chols(const std::vector<MatrixNd>& M, std::vector<MatrixNd>& M_chols)
{
  M_chols.resize(M.size());
  bool success = true;
  for (unsigned i=0; i< M_chols.size(); i++)
  {
    M_chols[i] = M[i];
    success = (LinAlgd::factor_chol(M_chols[i]) & success);
  }

  return success;
}

/// Solves Cholesky factorizations on a group of generalized inertia matrices
void RCArticulatedBodyInvDyn::solve_chol_fast(const std::vector<MatrixNd>& M_chols, VectorNd& v)
{
  unsigned gc_start = 0;
  for (unsigned i=0; i< M_chols.size(); i++)
  {
    unsigned gc_end = gc_start + M_chols[i].rows();
    SharedVectorNd vsub = v.segment(gc_start, gc_end);
    LinAlgd::solve_chol_fast(M_chols[i], vsub);
    gc_start = gc_end;
  }
}

/// Solves Cholesky factorizations on a group of generalized inertia matrices
void RCArticulatedBodyInvDyn::solve_chol_fast(const std::vector<MatrixNd>& M_chols, MatrixNd& M)
{
  unsigned gc_start = 0;

  for (unsigned i=0; i< M_chols.size(); i++)
  {
    unsigned gc_end = gc_start + M_chols[i].rows();
    SharedMatrixNd Msub = M.block(gc_start, gc_end, 0, M.columns());
    LinAlgd::solve_chol_fast(M_chols[i], Msub);
    gc_start = gc_end;
  }
}

/// Get the minimum index of vector v; if there are multiple minima (within zero_tol), returns one randomly
unsigned RCArticulatedBodyInvDyn::rand_min(const VectorNd& v, double zero_tol)
{
  static std::vector<unsigned> minima;
  minima.clear();
  unsigned minv = std::min_element(v.begin(), v.end()) - v.begin();
  minima.push_back(minv);
  for (unsigned i=0; i< v.rows(); i++)
    if (i != minv && v[i] < v[minv] + zero_tol)
      minima.push_back(i);
  return minima[rand() % minima.size()];
}

unsigned RCArticulatedBodyInvDyn::select_pivot(const VectorNd& znbas, const std::vector<unsigned>& nonbas, const std::vector<unsigned>& indices, double zero_tol)
{
  // if nonbas is empty, return INF
  if (nonbas.empty())
    return std::numeric_limits<unsigned>::max();

  // find all indices i of znbas for which znbas[i] < 0
  static std::vector<unsigned> neg;
  neg.clear();
  for (unsigned i=0; i< znbas.size(); i++)
    if (znbas[i] < -zero_tol)
      neg.push_back(i);

  // of all negative indices, find those which have a contact point with the
  // same link in nonbas
  static std::vector<unsigned> repeated;
  repeated.clear();
  for (unsigned i=0; i< neg.size(); i++)
  {
    // get the z index
    unsigned z_index = nonbas[neg[i]];

    // get the link index
    unsigned link_idx = indices[z_index];

    // loop through all other z indices
    for (unsigned j=0; j< nonbas.size(); j++)
      if (z_index != nonbas[j] && indices[nonbas[j]] == link_idx)
      {
        repeated.push_back(i);
        break;
      }
  }

  // if there are no such contact points, ...
  if (repeated.empty())
  {
    // pick the minimum entry of znbas
    return rand_min(znbas, zero_tol);
  }
  else if (repeated.size() > 1)
  {
/*
    // there are multiple such contact points, pick the one with the most
    // negative z value
    unsigned most_neg = repeated.front();
    for (unsigned i=1; i< repeated.size(); i++)
      if (znbas[repeated[i]] < znbas[repeated[most_neg]])
        most_neg = repeated[i];
*/
    // there are multiple such contact points, pick cardinally
    return repeated.front();
  }
  else
  {
    // exactly one entry, return it
    return repeated.front();
  }
}

/// Fast pivoting algorithm for denerate, monotone LCPs with few nonzero, nonbasic variables
/**
 * \param M the LCP matrix
 * \param q the LCP vector
 * \param indices a set of indices corresponding to each foot
 * \param z the solution is returned here; if z is an n-dimensional vector,
 *        attempts warm-starting
 * \param zero_tol the tolerance for solving to zero
 * \return true if the solver is successful
 */
bool RCArticulatedBodyInvDyn::lcp_fast(const MatrixNd& M, const VectorNd& q, const std::vector<unsigned>& indices, VectorNd& z, double zero_tol)
{
  const unsigned N = q.rows();
  const unsigned UINF = std::numeric_limits<unsigned>::max();

  // verify that indices are the right size
  assert(indices.size() == N);

  // determine the number of links represented
  _represented.clear();
  for (unsigned i=0; i< indices.size(); i++)
    if (indices[i] >= _represented.size())
      _represented.resize(indices[i]+1);

  // indicate that no link is represented in the non-basic set
  std::fill(_represented.begin(), _represented.end(), false);

  // look for trivial solution
  if (N == 0)
  {
    z.set_zero(0);
    return true;
  }

  // set zero tolerance if necessary
  if (zero_tol < 0.0)
    zero_tol = M.rows() * M.norm_inf() * std::numeric_limits<double>::epsilon();

  // get minimum element of q (really w)
  unsigned minw = std::min_element(q.begin(), q.end()) - q.begin();
  if (q[minw] > -zero_tol)
  {
    z.set_zero(N);
    return true;
  }

  // look for warm-start
  if (z.size() == N)
  {
    _nonbas.clear();
    for (unsigned i=0; i< N; i++)
      if (z[i] > zero_tol)
        _nonbas.push_back(i);
  }
  else
  {
    // set initial nonbasic set to have a negative variable from each link index
    _nonbas.clear();
    for (unsigned i=0; i< N; i++)
      if (q[i] < zero_tol && !_represented[indices[i]])
      {
        _nonbas.push_back(i);
        _represented[indices[i]] = true;
      }

    // now add contacts from all links not represented
    for (unsigned i=0; i< N; i++)
      if (std::binary_search(_nonbas.begin(), _nonbas.end(), i) &&
          !_represented[indices[i]])
      {
        _nonbas.push_back(i);
        insertion_sort(_nonbas.begin(), _nonbas.end());
        _represented[indices[i]] = true;
      }
  }

  // setup basic indices
  _bas.clear();
  for (unsigned i=0; i< N; i++)
    if (!std::binary_search(_nonbas.begin(), _nonbas.end(), i))
      _bas.push_back(i);

  // loop for maximum number of pivots
  const unsigned MAX_PIV = std::max(N*N, (unsigned) 1000);
  for (unsigned piv=0; piv < MAX_PIV; piv++)
  {
    // select nonbasic indices
    M.select_square(_nonbas.begin(), _nonbas.end(), _Msub);
    M.select(_bas.begin(), _bas.end(), _nonbas.begin(), _nonbas.end(), _Mmix);
    q.select(_nonbas.begin(), _nonbas.end(), _z);
    q.select(_bas.begin(), _bas.end(), _qbas);
    _z.negate();

    // solve for nonbasic z
    try
    {
      _LA.solve_fast(_Msub, _z);
    }
    catch (SingularException e)
    {
      return false;
    }

    // compute w and find minimum value
    _Mmix.mult(_z, _w) += _qbas;
    minw = (_w.rows() > 0) ? rand_min(_w, zero_tol) : UINF;

    // if w >= 0, check whether any component of z < 0
    if (minw == UINF || _w[minw] > -zero_tol)
    {
      // find the (a) minimum of z
      unsigned minz = select_pivot(_z, _nonbas, indices, zero_tol);
      if (minz < UINF && _z[minz] < -zero_tol)
      {
        // get the original index and remove it from the nonbasic set
        unsigned idx = _nonbas[minz];
        _nonbas.erase(_nonbas.begin()+minz);

        // move index to basic set and continue looping
        _bas.push_back(idx);
        insertion_sort(_bas.begin(), _bas.end());
      }
      else
      {
        // found the solution
        z.set_zero(N);

        // set values of z corresponding to _z
        for (unsigned i=0, j=0; j < _nonbas.size(); i++, j++)
          z[_nonbas[j]] = _z[i];

        return true;
      }
    }
    else
    {
      // one or more components of w violating w >= 0
      // move component of w from basic set to nonbasic set
      unsigned idx = _bas[minw];
      _bas.erase(_bas.begin()+minw);
      _nonbas.push_back(idx);
      insertion_sort(_nonbas.begin(), _nonbas.end());

      // look whether any component of z needs to move to basic set
      unsigned minz = select_pivot(_z, _nonbas, indices, zero_tol);
      if (minz < UINF &&_z[minz] < -zero_tol)
      {
        // move index to basic set and continue looping
        unsigned idx = _nonbas[minz];
        _nonbas.erase(_nonbas.begin()+minz);
        _bas.push_back(idx);
        insertion_sort(_bas.begin(), _bas.end());
      }
    }
  }

  // if we're here, then the maximum number of pivots has been exceeded
  return false;
}

