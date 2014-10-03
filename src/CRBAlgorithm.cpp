/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <queue>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/RNEAlgorithm.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/NumericalException.h>
#include <Moby/Spatial.h>
#include <Moby/CRBAlgorithm.h>

using namespace Ravelin;
using namespace Moby;
using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;
using boost::shared_array;
using boost::shared_ptr;

CRBAlgorithm::CRBAlgorithm()
{
}

/// Computes the parent array for sparse Cholesky factorization
void CRBAlgorithm::setup_parent_array()
{
  // get the number of generalized coordinates
  RCArticulatedBodyPtr body(_body);
  const unsigned N = body->num_generalized_coordinates(DynamicBody::eSpatial);

  // get explicit joints
  const vector<JointPtr>& ijoints = body->get_explicit_joints();

  // determine parent array (lambda)
  _lambda.resize(N);

  // set all values of lambda to inf initially
  for (unsigned i=0; i< N; i++)
    _lambda[i] = std::numeric_limits<unsigned>::max();
  
  // loop over all explicit joints
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the index of this joint
    unsigned idx = ijoints[i]->get_coord_index();

    // get the parent joint and its index
    RigidBodyPtr inboard = ijoints[i]->get_inboard_link();
    JointPtr parent = inboard->get_inner_joint_explicit();
    if (!parent)
      continue;
    unsigned pidx = parent->get_coord_index();

    // now set the elements of lambda
    for (unsigned j=0; j< ijoints[i]->num_dof(); j++)
    {
      _lambda[idx] = pidx;
      pidx = idx;
      idx++;
    }
  }
}

/// Factorizes (Cholesky) the generalized inertia matrix, exploiting sparsity
bool CRBAlgorithm::factorize_cholesky(MatrixNd& M)
{
  // get the number of degrees of freedom
  const unsigned n = M.rows();

  // loop
  for (unsigned kk=n; kk> 0; kk--)
  {
    unsigned k = kk - 1;
    if (M(k,k) < (double) 0.0)
    {
      assert(false);
      return false;
    }
    M(k,k) = std::sqrt(M(k,k));
    unsigned i = _lambda[k];
    while (i != std::numeric_limits<unsigned>::max())
    {
      M(k,i) /= M(k,k);
      i = _lambda[i];
    }
    i = _lambda[k];
    while (i != std::numeric_limits<unsigned>::max())
    {
      unsigned j=i;
      while (j != std::numeric_limits<unsigned>::max())
      {
        M(i,j) -= M(k,i)*M(k,j);
        j = _lambda[j];
      }
      i = _lambda[i];
    }
  }

  return true;
}

// Transforms (as necessary) and multiplies
void CRBAlgorithm::transform_and_mult(shared_ptr<const Pose3d> target, const SpatialRBInertiad& I, const vector<SVelocityd>& s, vector<SMomentumd>& Is)
{
  // transform s
  Pose3d::transform(I.pose, s, _sprime);

  // do the multiplication
  mult(I, _sprime, _Isprime);

  // do the transformation
  Pose3d::transform(target, _Isprime, Is);
}

/// Gets the frame in which kinematics and dynamics computations occur
/**
 * Only valid for bodies with floating bases.
 */
shared_ptr<const Pose3d> CRBAlgorithm::get_computation_frame(RCArticulatedBodyPtr body)
{
  assert(body->is_floating_base());

  // get the base 
  RigidBodyPtr base = body->get_base_link();

  switch (body->get_computation_frame_type())
  {
    case eLink:
    case eJoint:
      return base->get_pose();

    case eLinkInertia:
      return base->get_inertial_pose();

    case eLinkCOM:
      return base->get_gc_pose();

    case eGlobal:
      return shared_ptr<const Pose3d>();

    default:
      assert(false);
  }

  return shared_ptr<const Pose3d>();
}

/// Calculates the generalized inertia of this body
/**
 * Specialized function for use with the CRB algorithm
 */
void CRBAlgorithm::calc_generalized_inertia(RCArticulatedBodyPtr body)
{
  const unsigned SPATIAL_DIM = 6;

  // get the appropriate M
  MatrixNd& M = this->_M;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // get explicit joints
  const vector<JointPtr>& ijoints = body->get_explicit_joints();

  // compute the joint space inertia
  calc_joint_space_inertia(body, _H, _Ic);

  // get the number of base degrees-of-freedom
  const unsigned N_BASE_DOF = (body->is_floating_base()) ? 6 : 0;

  // resize M
  M.resize(N_BASE_DOF + body->num_joint_dof_explicit(), N_BASE_DOF + body->num_joint_dof_explicit());

  // set appropriate part of H
  M.set_sub_mat(0, 0, _H);

  // see whether we are done
  if (!body->is_floating_base())
    return;

  // ************************************************************************
  // floating base
  // ************************************************************************

  // setup the indices for the base
  const unsigned BASE_START = body->num_joint_dof_explicit();

  // get components of M
  SharedMatrixNd Ic0 = M.block(BASE_START, M.rows(), BASE_START, M.columns());
  SharedMatrixNd KS = M.block(0, BASE_START, BASE_START, M.columns()); 
  SharedMatrixNd K = M.block(BASE_START, M.rows(), 0, BASE_START); 

  // get composite inertia in desired frame
  shared_ptr<const Pose3d> P = get_computation_frame(body);
  Pose3d::transform(P, _Ic.front()).to_PD_matrix(Ic0);

  FILE_LOG(LOG_DYNAMICS) << "Ic0: " << std::endl << Ic0;

  // compute K
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the spatial axes for the joint
    JointPtr joint = ijoints[i];
    const std::vector<SVelocityd>& s = joint->get_spatial_axes();
    if (joint->num_dof() == 0)
      continue;

    // get the index for this joint
    unsigned jidx = joint->get_coord_index();

    // get the outboard link and link index
    RigidBodyPtr outboard = joint->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // transform and multiply
    transform_and_mult(P, _Ic[oidx], s, _Is);

    // compute the requisite columns of K
    SharedMatrixNd Kb = K.block(0, SPATIAL_DIM, jidx, jidx+joint->num_dof()); 
    SharedMatrixNd KSb = KS.block(jidx, jidx+joint->num_dof(), 0, SPATIAL_DIM); 
    to_matrix(_Is, Kb); 
    MatrixNd::transpose(Kb, KSb); 
  }

  FILE_LOG(LOG_DYNAMICS) << "[H K'; K Ic0] (permuted): " << std::endl << M;
}

/// Computes *just* the joint space inertia matrix
void CRBAlgorithm::calc_joint_space_inertia(RCArticulatedBodyPtr body, MatrixNd& H, vector<SpatialRBInertiad>& Ic)
{
  queue<RigidBodyPtr> link_queue;
  const unsigned SPATIAL_DIM = 6;

  // get the reference frame
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // get the sets of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_explicit_joints();
  const vector<JointPtr>& joints = body->get_joints();

  // set the composite inertias to the isolated inertias initially 
  Ic.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    Ic[links[i]->get_index()] = links[i]->get_inertia();

    // check for degenerate inertia
    #ifndef NDEBUG
    SpatialRBInertiad& J = Ic[links[i]->get_index()];
    if (links[i]->is_base() && body->is_floating_base() && 
        (J.m <= 0.0 || J.J.norm_inf() <= 0.0))
      throw std::runtime_error("Attempted to compute dynamics given degenerate inertia for a floating base body");
    #endif

    if (LOGGING(LOG_DYNAMICS))
    {
      MatrixNd X;
      Ic[links[i]->get_index()].to_matrix(X);
      FILE_LOG(LOG_DYNAMICS) << "isolated inertia for link " << links[i]->id << ": " << endl << X;
    }
  }

  // ************************************************************************
  // first, determine the supports for the joints and the number of joint DOF
  // ************************************************************************
 
  // create and initialize the supports array
  _supports.resize(joints.size());
  for (unsigned i=0; i< joints.size(); i++)
  {
    _supports[i].resize(links.size());
    for (unsigned j=0; j< links.size(); j++)
      _supports[i][j] = false;
  }

  // add all leaf links to the link queue
  for (unsigned i=1; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);

  // process until the queue is empty
  while (!link_queue.empty())
  {
    // get the element out of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // get the explicit inner joint for this link
    JointPtr joint(link->get_inner_joint_explicit());
    unsigned jidx = joint->get_index();
    assert(joint);

    // add this link to the support for the joint
    _supports[jidx][link->get_index()] = true;

    // add all supports from the outer joints of this link
    list<RigidBodyPtr> child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr child, child_links)
    {
      // don't process children with lower link indices (loops)
      if (child->get_index() < link->get_index())
        continue;

      // get the inner explicit joint
      JointPtr child_joint = child->get_inner_joint_explicit();
      unsigned jiidx = child_joint->get_index();

      // setup the supports
      for (unsigned i=0; i< links.size(); i++)
        if (_supports[jiidx][i])
          _supports[jidx][i] = true;
    }  

    // add the parent of this link to the queue, if it is not the base and
    // if the parent's link index is lower
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base() && parent->get_index() < link->get_index())
      link_queue.push(parent);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    FILE_LOG(LOG_DYNAMICS) << "supports: " << endl;
    for (unsigned i=0; i< _supports.size(); i++)
    {
      std::ostringstream str;
      for (unsigned j=0; j< _supports[i].size(); j++)
        str << _supports[i][j] << " ";
      FILE_LOG(LOG_DYNAMICS) << i << ": " << str.str() << endl;
    }
  }

  // resize H 
  H.set_zero(body->num_joint_dof_explicit(), body->num_joint_dof_explicit());

  // ************************************************************************
  // compute spatial composite inertias 
  // ************************************************************************

  // now determine the composite inertias
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // put all leaf links into a queue
  assert(link_queue.empty());
  for (unsigned i=1; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);

  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // get the index for this link
    unsigned i = link->get_index();
    
    // see whether this link has already been processed
    if (body->_processed[i])
      continue;

    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue;    

    // process the parent link, if possible
    RigidBodyPtr parent(link->get_parent_link());
    if (parent && i > parent->get_index())
    {
      // put the parent on the queue
      link_queue.push(parent);
    
      // get the parent index
      unsigned h = parent->get_index();
    
      // add this inertia to its parent
      Ic[h] += Pose3d::transform(Ic[h].pose, Ic[i]); 

      if (LOGGING(LOG_DYNAMICS))
      {
        MatrixNd X;
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (child) link " << link->id << ": " << std::endl << Ic[i].to_matrix(X);
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (child) link " << link->id << ": " << std::endl << Ic[i];
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (parent) link " << parent->id << ": " << std::endl << Ic[h].to_matrix(X);
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (parent) link " << parent->id << ": " << std::endl << Ic[h];
      }
    }

    // indicate that the link has been processed
    body->_processed[i] = true;
  }

  // ************************************************************************
  // compute H
  // ************************************************************************

  // compute the forces
  _momenta.resize(links.size());
  for (unsigned i=0; i < ijoints.size(); i++)
  {
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link(); 
    unsigned oidx = outboard->get_index();
    const std::vector<SVelocityd>& s = ijoints[i]->get_spatial_axes();
    Pose3d::transform(Ic[oidx].pose, s, _sprime);
    mult(Ic[oidx], _sprime, _momenta[oidx]);
    if (_sprime.size() > 0)
    {
      FILE_LOG(LOG_DYNAMICS) << "s: " << _sprime[0] << std::endl;
      FILE_LOG(LOG_DYNAMICS) << "Is[" << i << "]: " << _momenta[oidx][0] << std::endl;
    }
  } 

  // setup H
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the number of degrees of freedom for joint i
    const unsigned NiDOF = ijoints[i]->num_dof();

    // get the starting coordinate index for this joint
    unsigned iidx = ijoints[i]->get_coord_index();

    // get the outboard link for joint i
    RigidBodyPtr obi = ijoints[i]->get_outboard_link();
    unsigned oiidx = obi->get_index();

    // get the spatial axes for jointi
    const std::vector<SVelocityd>& s = ijoints[i]->get_spatial_axes();

    // get the appropriate submatrix of H
    SharedMatrixNd subi = H.block(iidx, iidx+NiDOF, iidx, iidx+NiDOF); 

    // compute the H term for i,i
    transform_and_transpose_mult(s, _momenta[oiidx], subi);

    // determine what will be the new value for m
    for (unsigned j=i+1; j< ijoints.size(); j++)
    {
      // get the number of degrees of freedom for joint j 
      const unsigned NjDOF = ijoints[j]->num_dof();

      // get the outboard link for joint j
      RigidBodyPtr obj = ijoints[j]->get_outboard_link();
      unsigned ojidx = obj->get_index();

      // if link j is not supported by joint i, contribution to H is zero
      if (!_supports[ijoints[i]->get_index()][ojidx])
        continue;

      // get the starting coordinate index for joint j
      unsigned jidx = ijoints[j]->get_coord_index();

      // get the appropriate submatrices of H
      SharedMatrixNd subj = H.block(iidx, iidx+NiDOF, jidx, jidx+NjDOF); 
      SharedMatrixNd subjT = H.block(jidx, jidx+NjDOF, iidx, iidx+NiDOF); 

      // compute the appropriate submatrix of H
      transform_and_transpose_mult(s, _momenta[ojidx], subj);

      // set the transposed part
      MatrixNd::transpose(subj, subjT);
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "joint space inertia: " << endl << H;
}

/// Calculates the generalized inertia matrix for the given representation
/**
 * Generic method provided for use with generalized coordinates.
 */
void CRBAlgorithm::calc_generalized_inertia(SharedMatrixNd& M)
{
  // do the precalculation
  RCArticulatedBodyPtr body(_body);
  precalc(body);

  // get the set of links
  ReferenceFrameType rftype = body->get_computation_frame_type();
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_explicit_joints();

// TODO: this should be able to be removed; this function should already be
//       computed using precalc(.)
  // get the joint space inertia and composite inertias
//  calc_joint_space_inertia(body, _H, _Ic);

  // get the number of base degrees-of-freedom
  const unsigned N_BASE_DOF = (body->is_floating_base()) ? 6 : 0;
  const unsigned SPATIAL_DIM = 6;

  // resize M and set H
  M.resize(N_BASE_DOF + body->num_joint_dof_explicit(), N_BASE_DOF + body->num_joint_dof_explicit());
  M.set_sub_mat(0, 0, _H);

  // look for simplest case
  if (!body->is_floating_base())
    return;

  // get the coordinates at which the base start 
  const unsigned BASE_START = M.rows() - SPATIAL_DIM;

  // ************************************************************************
  // floating base
  // ************************************************************************

  // get the number of explicit joint degrees-of-freedom
  const unsigned NjDOF = body->num_joint_dof_explicit();

  // get components of M
  SharedMatrixNd Ic0 = M.block(BASE_START, M.rows(), BASE_START, M.columns());
  SharedMatrixNd KS = M.block(0, BASE_START, BASE_START, M.columns()); 
  SharedMatrixNd K = M.block(BASE_START, M.rows(), 0, BASE_START); 

  // get composite inertia in desired frame
  RigidBodyPtr base = body->get_base_link();
  shared_ptr<const Pose3d> P = base->get_gc_pose();
  Pose3d::transform(P, _Ic[base->get_index()]).to_PD_matrix(Ic0); 

  // compute K
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the spatial axes for the joint
    JointPtr joint = ijoints[i];
    const std::vector<SVelocityd>& s = joint->get_spatial_axes();
    if (joint->num_dof() == 0)
      continue;

    // get the index for this joint
    unsigned jidx = joint->get_coord_index();

    // get the outboard link and link index
    RigidBodyPtr outboard = joint->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // transform and multiply
    transform_and_mult(P, _Ic[oidx], s, _Is);

    // compute the requisite columns of K
    SharedMatrixNd Kb = K.block(0, SPATIAL_DIM, jidx, jidx+joint->num_dof()); 
    SharedMatrixNd KSb = KS.block(jidx, jidx+joint->num_dof(), 0, SPATIAL_DIM); 
    to_matrix(_Is, Kb); 
    MatrixNd::transpose(Kb, KSb);
  }

  FILE_LOG(LOG_DYNAMICS) << "[H K'; K Ic0] (permuted): " << std::endl << M;
}

/// Performs necessary pre-computations for computing accelerations or applying impulses
void CRBAlgorithm::precalc(RCArticulatedBodyPtr body)
{
  // get the links and joints for the body
  const vector<JointPtr>& joints = body->get_explicit_joints();

  // compute spatial isolated inertias and generalized inertia matrix
  // do the calculations
  calc_generalized_inertia(body);

  // attempt to do a Cholesky factorization of M
  MatrixNd& fM = this->_fM;
  MatrixNd& M = this->_M;
  if ((_rank_deficient = !_LA->factor_chol(fM = M)))
  {
    fM = M;
    _LA->svd(fM, _uM, _sM, _vM);
  }
}

/// Executes the composite rigid-body method
void CRBAlgorithm::calc_fwd_dyn()
{
  // get the body
  RCArticulatedBodyPtr body(_body);

  // do necessary pre-calculations
  precalc(body);

  // execute the appropriate algorithm
  if (body->is_floating_base())
    calc_fwd_dyn_floating_base(body);
  else
    calc_fwd_dyn_fixed_base(body);
   
  // update the link accelerations
  update_link_accelerations(body);
}

/// Executes the composite rigid-body method without computing and factorizing inertia matrix
/**
 * This method is useful when the inertia matrix has already been computed-
 * considerable computation will then be avoided.
 */
void CRBAlgorithm::calc_fwd_dyn_special()
{
  // get the body and the reference frame
  RCArticulatedBodyPtr body(_body);

  // execute the appropriate algorithm
  if (body->is_floating_base())
    calc_fwd_dyn_floating_base(body);
  else
    calc_fwd_dyn_fixed_base(body);
   
  // update the link accelerations
  update_link_accelerations(body);
}

/// Solves for acceleration using the body inertia matrix
VectorNd& CRBAlgorithm::M_solve(VectorNd& xb) 
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  precalc(body);

  // setup xb
  SharedVectorNd xb_shared = xb.segment(0, xb.rows());

  M_solve_noprecalc(xb_shared); 
  return xb;
}

/// Solves for acceleration using the body inertia matrix
SharedVectorNd& CRBAlgorithm::M_solve(SharedVectorNd& xb) 
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  precalc(body);

  M_solve_noprecalc(xb); 
  return xb;
}

/// Solves for acceleration using the body inertia matrix
MatrixNd& CRBAlgorithm::M_solve(MatrixNd& XB)
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  precalc(body);

  // setup XB
  SharedMatrixNd XB_shared = XB.block(0, XB.rows(), 0, XB.columns());

  M_solve_noprecalc(XB_shared); 
  return XB;
}

/// Solves for acceleration using the body inertia matrix
SharedMatrixNd& CRBAlgorithm::M_solve(SharedMatrixNd& XB)
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  precalc(body);

  M_solve_noprecalc(XB); 
  return XB;
}

/// Solves for acceleration using the body inertia matrix
VectorNd& CRBAlgorithm::M_solve_noprecalc(VectorNd& xb)
{
  // setup xb
  SharedVectorNd xb_shared = xb.segment(0, xb.rows());

  // solve
  M_solve_noprecalc(xb_shared);

  return xb;
}

/// Solves for acceleration using the body inertia matrix
SharedVectorNd& CRBAlgorithm::M_solve_noprecalc(SharedVectorNd& xb)
{
  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
    _LA->solve_LS_fast(_uM, _sM, _vM, xb);
  else
    _LA->solve_chol_fast(_fM, xb);

  return xb;
}

/// Solves for acceleration using the body inertia matrix
MatrixNd& CRBAlgorithm::M_solve_noprecalc(MatrixNd& XB)
{
  // get a shared matrix
  SharedMatrixNd XB_shared = XB.block(0, XB.rows(), 0, XB.columns());

  // solve
  M_solve_noprecalc(XB_shared);

  return XB;
}

/// Solves for acceleration using the body inertia matrix
SharedMatrixNd& CRBAlgorithm::M_solve_noprecalc(SharedMatrixNd& XB)
{
  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
    _LA->solve_LS_fast(_uM, _sM, _vM, XB);
  else
    _LA->solve_chol_fast(_fM, XB);

  return XB;
}

/// Executes the composite rigid-body method on an articulated body with a fixed base
void CRBAlgorithm::calc_fwd_dyn_fixed_base(RCArticulatedBodyPtr body)
{
  // get the set of links and joints for the articulated body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_explicit_joints();

  // ***********************************************************************
  // first, calculate C
  // ***********************************************************************

  // call inverse dynamics to calculate C
  SForced f0;
  calc_generalized_forces(f0, _C);
  
  // get the number of degrees-of-freedom
  unsigned nDOF = _C.size();

  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************
  _Q.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    ijoints[i]->get_scaled_force(_Qi);
    _Qi += ijoints[i]->ff;
    unsigned j = ijoints[i]->get_coord_index();
    _Q.set_sub_vec(j, _Qi);
  }

  FILE_LOG(LOG_DYNAMICS) << "H: " << std::endl << _M;
  FILE_LOG(LOG_DYNAMICS) << "C: " << _C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "Q: " << _Q << std::endl;

  // subtract C from Q
  _Q -= _C;

  // get the pointer to the joint-space acceleration vector
  VectorNd& qdd = this->_qdd;

  // compute joint accelerations
  M_solve_noprecalc(qdd = _Q);

  FILE_LOG(LOG_DYNAMICS) << "qdd: " << qdd << std::endl;

  // set qdd
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned j = ijoints[i]->get_coord_index();
    qdd.get_sub_vec(j, j+ijoints[i]->num_dof(), ijoints[i]->qdd);
  }

  // set spatial acceleration of the base
  this->_a0.set_zero(); 
  SAcceld abase = links.front()->get_accel();
  links.front()->set_accel(SAcceld::zero(links.front()->get_computation_frame()));
}

/// Executes the composite rigid-body method on an articulated body with a floating base
/**
 * This algorithm is taken from [Featherstone, 1987], p. 123.  This is only
 * calculated in the global frame.
 */
void CRBAlgorithm::calc_fwd_dyn_floating_base(RCArticulatedBodyPtr body)
{
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() entered" << std::endl;

  // get the set of links and explicit joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_explicit_joints();

  // calculate C
  SForced f0;
  calc_generalized_forces(f0, _C);

  // compute generalized forces in proper frame
  SForced f0x = Pose3d::transform(get_computation_frame(body), f0);

  // get the number of degrees-of-freedom
  unsigned nDOF = _C.size();
  
  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************

  _Q.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    ijoints[i]->get_scaled_force(_Qi);
    _Qi += ijoints[i]->ff;
    unsigned j = ijoints[i]->get_coord_index();
    _Q.set_sub_vec(j, _Qi);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    Pose3d base_pose = links[0]->get_pose();
    base_pose.update_relative_pose(GLOBAL);
    FILE_LOG(LOG_DYNAMICS) << "base pose: " << base_pose << std::endl;
  }
  FILE_LOG(LOG_DYNAMICS) << "Q: " << _Q << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "C: " << _C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "M: " << std::endl << this->_M;

  // setup the simulataneous equations to solve: [Featherstone, 1987], eq. 7.24
  concat(_Q -= _C, -f0x, _b);
  FILE_LOG(LOG_DYNAMICS) << "b: " << _b << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "link + external forces on base: " << f0x << std::endl;

  // solve for accelerations
  M_solve_noprecalc(_augV = _b); 
  FILE_LOG(LOG_DYNAMICS) << "b: " << _b << std::endl;

  // get pointers to a0 and qdd vectors
  SAcceld& a0 = this->_a0;
  VectorNd& qdd = this->_qdd;

  // swap components of a0
  const unsigned SPATIAL_DIM = 6;
  const unsigned BASE_START = _augV.size() - SPATIAL_DIM;
  std::swap(_augV[BASE_START+0], _augV[BASE_START+3]);
  std::swap(_augV[BASE_START+1], _augV[BASE_START+4]);
  std::swap(_augV[BASE_START+2], _augV[BASE_START+5]);

  // get out a0, qdd
  qdd = _augV.segment(0, BASE_START);
  a0 = _augV.segment(BASE_START,_augV.size());
  a0.pose = f0x.pose;

  // set the base acceleration
  links.front()->set_accel(a0);

  FILE_LOG(LOG_DYNAMICS) << "base spatial acceleration: " << a0 << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "joint accelerations:";

  // set qdd
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned j = ijoints[i]->get_coord_index();
    qdd.get_sub_vec(j, j+ijoints[i]->num_dof(), ijoints[i]->qdd);
  }

  FILE_LOG(LOG_DYNAMICS) << qdd << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() exited" << std::endl;
}

/// Computes the vector "C", used for forward dynamics computation, using the recursive Newton-Euler algorithm
/**
 * \return the spatial vector of forces on the base, which can be ignored for forward dynamics for fixed bases
 */
void CRBAlgorithm::calc_generalized_forces(SForced& f0, VectorNd& C)
{
  const unsigned SPATIAL_DIM = 6;
  queue<RigidBodyPtr> link_queue;
  SForced w;

  // get the body and the reference frame
  RCArticulatedBodyPtr body(_body);

  // get the set of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_explicit_joints();
  if (links.empty())
  {
    C.resize(0);
    f0.set_zero();
    return;
  }

  // **************************************************************************
  // first, compute forward dynamics using RNE algorithm; we copy the algorithm
  // here, because we need some data out of it
  // **************************************************************************

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() entered" << std::endl;

  // ** STEP 1: compute accelerations

  // get the base link
  RigidBodyPtr base = links.front();

  // setup the map of link accelerations
  _a.resize(links.size());

  // setup the acceleration for the base
  _a[base->get_index()].set_zero(base->get_velocity().pose);
  
  // add all child links of the base to the processing queue
  list<RigidBodyPtr> children;
  base->get_child_links(std::back_inserter(children));
  BOOST_FOREACH(RigidBodyPtr rb, children)
    link_queue.push(rb);
    
  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
  body->_processed[base->get_index()] = true;

  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();  
    unsigned i = link->get_index();
    body->_processed[i] = true;

    // push all children of the link onto the queue, unless they were already
    // processed  
    list<RigidBodyPtr> child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the link's parent
    RigidBodyPtr parent(link->get_parent_link());
    unsigned h = parent->get_index();

    // get the joint for this link
    JointPtr joint(link->get_inner_joint_explicit());

    // get the spatial link velocity
    const SVelocityd& vx = link->get_velocity(); 

    // get spatial axes and derivative for this link's inner joint
    const std::vector<SVelocityd>& s = joint->get_spatial_axes();
    const std::vector<SVelocityd>& sdot = joint->get_spatial_axes_dot();

    // get the current joint velocity
    const VectorNd& qd = joint->qd;

    // **** compute acceleration

    // add this link's contribution
    Pose3d::transform(_a[i].pose, s, _sprime);
    if (_sprime.empty())
      _a[i].set_zero(vx.pose);
    else
    {
      SVelocityd sqd = mult(_sprime, qd);
      sqd.pose = vx.pose;
      _a[i] = vx.cross(sqd);
    }
    if (!sdot.empty())
      _a[i] += SAcceld(mult(Pose3d::transform(_a[i].pose, sdot, _sprime), qd)); 

    // now add parent's contribution
    _a[i] += Pose3d::transform(_a[i].pose, _a[h]);
    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->id << std::endl;
    if (s.size() > 0)
      FILE_LOG(LOG_DYNAMICS) << "  spatial joint velocity: " << (mult(s,qd)) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << link->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << _a[i] << std::endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  // use a map to determine which links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // setup a map of link forces, all set to zero initially
  _w.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    _w[i].set_zero();
    _w[i].pose = links[i]->get_computation_frame();
  }

  // add all leaf links to the queue
  for (unsigned i=0; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // if this link has already been processed, do not process it again
    if (body->_processed[i])
      continue;

    // verify all children have been processed
    if (!body->all_children_processed(link))
      continue;
   
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << _w[i] << std::endl;    
    if (LOGGING(LOG_DYNAMICS) && link != body->get_base_link())
      FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (link->get_inertia() * _a[i]) << std::endl;

    // add I*a to the link force and Euler torque components 
    const SVelocityd& vx = link->get_velocity(); 
    _w[i] += link->get_inertia() * _a[i];
    _w[i] += vx.cross(link->get_inertia() * vx);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << _w[i] << std::endl;

    // subtract external forces
    SForced wext = link->sum_forces(); 
    _w[i] -= wext;
    FILE_LOG(LOG_DYNAMICS) << "  external forces: " << wext << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << _w[i] << std::endl;

    // update the parent force and add parent for processing (if parent)
    RigidBodyPtr parent = link->get_parent_link();
    if (parent)
    {
      unsigned h = parent->get_index();
      _w[h] += Pose3d::transform(_w[h].pose, _w[i]);
      link_queue.push(parent);
    }

    // indicate that this link has been processed
    body->_processed[i] = true;
  }
  
  // ** STEP 3: compute actuator forces (C)

  // determine the length of the C vector
  const unsigned nDOF = body->num_joint_dof_explicit();

  // compute actuator forces (C)
  C.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get links, joints, etc.
    JointPtr joint = ijoints[i];
    RigidBodyPtr ob = joint->get_outboard_link();

    // get indices
    unsigned jidx = joint->get_coord_index();
    unsigned oidx = ob->get_index();

    // compute appropriate components of C
    SharedVectorNd Csub = C.segment(jidx, jidx+joint->num_dof()); 
    const std::vector<SVelocityd>& s = joint->get_spatial_axes();
    transform_and_transpose_mult(s, _w[oidx], Csub);

    FILE_LOG(LOG_DYNAMICS) << " -- computing C for link " << ob->id << std::endl;
    if (LOGGING(LOG_DYNAMICS))
    {
      SForced w_mixed = Pose3d::transform(ob->get_mixed_pose(), _w[oidx]);
//      SForced w_mixed = Pose3d::transform(_w[oidx].pose, _w[oidx]);
      std::vector<SVelocityd> s_mixed;
      Pose3d::transform(ob->get_mixed_pose(), s, s_mixed);
//      Pose3d::transform(_w[oidx].pose, s, s_mixed);
      if (!s.empty())
      {
        FILE_LOG(LOG_DYNAMICS) << " -- s' (mixed pose): " << s_mixed[0] << std::endl;
        FILE_LOG(LOG_DYNAMICS) << " -- force (mixed pose): " << w_mixed << std::endl;
      }
    }
    FILE_LOG(LOG_DYNAMICS) << "   -- forces: " << _w[oidx] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "   -- component of C: " << Csub << std::endl;
  }  

  FILE_LOG(LOG_DYNAMICS) << "------------------------------------------------" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "forces on base: " << _w[0] << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() exited" << std::endl;

  // store forces on base
  f0 = _w[0]; 
}

/// Updates all link accelerations (except the base)
void CRBAlgorithm::update_link_accelerations(RCArticulatedBodyPtr body)
{
  queue<RigidBodyPtr> link_queue;

  // get the set of links and their velocities
  const vector<RigidBodyPtr>& links = body->get_links();

  // if there are no links, there is nothing to do
  if (links.empty())
    return;

  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // get the base link and mark it as processed
  RigidBodyPtr base = links.front();
  body->_processed[base->get_index()] = true;
  
  // get the spatial acceleration of the base link (should have already been
  // computed)
  base->set_accel(this->_a0);

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::update_link_accelerations() entered" << std::endl;
  
  // add all children of the base to the link queue
  list<RigidBodyPtr> children;
  base->get_child_links(std::back_inserter(children));
  BOOST_FOREACH(RigidBodyPtr rb, children)
    link_queue.push(rb);
  
  // propagate link accelerations 
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();
    body->_processed[i] = true;

    // push all unprocessed children of the link onto the queue
    list<RigidBodyPtr> children;
    link->get_child_links(std::back_inserter(children));
    BOOST_FOREACH(RigidBodyPtr rb, children)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the parent link
    JointPtr joint(link->get_inner_joint_explicit());
    RigidBodyPtr parent(link->get_parent_link());
    unsigned h = parent->get_index();
 
    // set link acceleration
    const SAcceld& ah = parent->get_accel();
    SAcceld ai = Pose3d::transform(link->get_accel().pose, ah);

    // get the link spatial axis
    const std::vector<SVelocityd>& s = joint->get_spatial_axes(); 

    // determine the link accel
    Pose3d::transform(ai.pose, s, _sprime);
    if (!_sprime.empty())
    {
      SVelocityd sqd = mult(_sprime, joint->qd);
      ai += SAcceld(link->get_velocity().cross(sqd));
      ai += SAcceld(mult(_sprime, joint->qdd)); 
    }
    link->set_accel(ai);

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent acceleration: " << ah << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- velocity: " << link->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qdd: " << joint->qdd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- acceleration: " << ai << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::update_link_accelerations() exited" << std::endl;
}

/*
/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
void CRBAlgorithm::apply_impulse(const SForced& w, RigidBodyPtr link)
{
  // An alternative method for applying impulses using generalized coordinates
  // below...
  SMatrix6 Xcp = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, link->get_position());
  SMatrix6 X0p = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, ZEROS_3);
  SVector6 jx(jj[0], jj[1], jj[2], jk[0], jk[1], jk[2]);
 
  // convert the impulse to a generalized impulse
  VectorNd gf;
  body->convert_to_generalized_force(link, point, jj, jk, gf);

  // get the generalized inertia and invert it
  MatrixNd M;
  body->get_generalized_inertia(M);
  LinAlg::inverse_PD(M);

  // get the generalized velocity
  VectorNd gv;
  body->get_generalized_velocity(gv);

  // compute the change in velocity
  VectorNd dv;
  M.mult(gf, dv);

  SMatrix6 Xc0 = SMatrix6::calc_spatial_transform(IDENTITY_3x3, ZEROS_3, IDENTITY_3x3, link->get_position());
  SMatrix6 X0c = SMatrix6::calc_spatial_transform(IDENTITY_3x3, link->get_position(), IDENTITY_3x3, ZEROS_3);
  MatrixNd Mx;
  body->get_generalized_inertia(Mx, eGlobal);
  SMatrix6 M2(Mx.begin());

  FILE_LOG(LOG_DYNAMICS) << "generalized inertia (global frame): " << std::endl << Mx;
  FILE_LOG(LOG_DYNAMICS) << "spatial transform (global to centroidal): " << std::endl << Xc0;
  FILE_LOG(LOG_DYNAMICS) << "[permuted] spatial inertia (centroidal frame): " << std::endl << (Xc0 * M2 * X0c);

  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (global frame): " << (X0p * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (centroidal frame): " << (Xcp * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "base transform: " << std::endl << body->get_links().front()->get_pose();
  FILE_LOG(LOG_DYNAMICS) << "impulse: " << jj << " / " << (jk + Vector3::cross(jj, r)) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "generalized impulse: " << gf << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "inverse generalized inertia: " << std::endl << M;
  FILE_LOG(LOG_DYNAMICS) << "generalized v: " << gv << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "delta v: " << dv << std::endl;
  gv += dv;
  FILE_LOG(LOG_DYNAMICS) << "new v: " << gv << std::endl;
  body->set_generalized_velocity(gv);
  FILE_LOG(LOG_DYNAMICS) << "new base linear velocity: " << body->get_links().front()->get_lvel() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "new base angular velocity: " << body->get_links().front()->get_avel() << std::endl;
}
*/

/// TODO: fix this for bodies with kinematic loops
/// Applies an impulse to an articulated body with a floating base; complexity O(n^2)
void CRBAlgorithm::apply_impulse(const SMomentumd& w, RigidBodyPtr link)
{
  const unsigned OPSPACE_DIM = 6;
  SVelocityd dv0;
  VectorNd& b = _b;
  VectorNd& augV = _augV;
  VectorNd& workv = _workv;

  // get the body
  RCArticulatedBodyPtr body(_body);

  // do necessary pre-calculations
  precalc(body); 

  // does not work for bodies with kinematic loops
  assert(body->_ijoints.empty());

  // get the base link
  RigidBodyPtr base = body->get_base_link();

  // compute the Jacobian for the floating base w.r.t. the contact point
  _J.clear();

  // compute the Jacobian with respect to the contact point
  RigidBodyPtr l = link;
  JointPtr j;
  while ((j = l->get_inner_joint_explicit()))
  {
    // compute the Jacobian column(s) for the joint
    const std::vector<SVelocityd>& s = j->get_spatial_axes();

    // transform spatial axes to global frame
    Pose3d::transform(GLOBAL, s, _sprime);

    // set the column(s) of the Jacobian
    _J.insert(_J.end(), _sprime.begin(), _sprime.end());

    // set l to its parent
    l = RigidBodyPtr(l->get_parent_link());
  }

  // transform the impulse to the global frame
  SMomentumd w0 = Pose3d::transform(GLOBAL, w); 

  // compute the impulse applied to the joints
  transpose_mult(_J, w0, workv);

  FILE_LOG(LOG_DYNAMICS) << "  impulse (last frame): " << w0 << std::endl;

  // special case: floating base
  if (body->_floating_base)
  {
    // determine the index where the base starts
    const unsigned BASE_START = workv.size();

    // form vector to solve for b
    concat(workv, w0, b);
 
    // compute changes in base and joint velocities
    M_solve_noprecalc(workv = b);

    // swap base velocity change linear and angular components 
    const unsigned BASE_A = BASE_START + 0;
    const unsigned BASE_B = BASE_START + 1;
    const unsigned BASE_G = BASE_START + 2;
    const unsigned BASE_X = BASE_START + 3;
    const unsigned BASE_Y = BASE_START + 4;
    const unsigned BASE_Z = BASE_START + 5;
    std::swap(workv[BASE_A], workv[BASE_X]);
    std::swap(workv[BASE_B], workv[BASE_Y]);
    std::swap(workv[BASE_G], workv[BASE_Z]);
 
    // get change in base and change in joint velocities
    Vector3d dv0_angular(workv[BASE_A], workv[BASE_B], workv[BASE_G]);
    Vector3d dv0_linear(workv[BASE_X], workv[BASE_Y], workv[BASE_Z]);
    SVelocityd dv0(dv0_angular, dv0_linear, w0.pose);

    // update the base velocity
    SVelocityd basev = base->get_velocity();
    basev += Pose3d::transform(basev.pose, dv0);
    base->set_velocity(basev);

    FILE_LOG(LOG_DYNAMICS) << "  change in base velocity: " << dv0 << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << basev << std::endl;
  }

  // apply the change and update link velocities
  const vector<JointPtr>& ijoints = body->get_explicit_joints();
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned idx = ijoints[i]->get_coord_index();
    FILE_LOG(LOG_DYNAMICS) << " joint " << ijoints[i]->id << " qd: " << ijoints[i]->qd << "  dqd: " << workv.segment(idx, idx+ijoints[i]->num_dof()) << std::endl;  
    ijoints[i]->qd += workv.segment(idx, idx+ijoints[i]->num_dof());
  }  
  body->update_link_velocities();

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::apply_impulse() exited" << std::endl;
}


