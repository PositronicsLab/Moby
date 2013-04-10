/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <queue>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/RNEAlgorithm.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/NumericalException.h>
#include <Moby/CRBAlgorithm.h>

using namespace Moby;
using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;
using boost::shared_array;

CRBAlgorithm::CRBAlgorithm()
{
  _position_data_valid = false;
  _velocity_data_valid = false;
}

/// Computes the parent array for sparse Cholesky factorization
void CRBAlgorithm::setup_parent_array()
{
  // get the number of generalized coordinates
  RCArticulatedBodyPtr body(_body);
  const unsigned N = body->num_generalized_coordinates(DynamicBody::eAxisAngle);

  // get implicit joints
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // determine parent array (lambda)
  _lambda.resize(N);

  // set all values of lambda to inf initially
  for (unsigned i=0; i< N; i++)
    _lambda[i] = std::numeric_limits<unsigned>::max();
  
  // loop over all implicit joints
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the index of this joint
    unsigned idx = ijoints[i]->get_coord_index();

    // get the parent joint and its index
    RigidBodyPtr inboard = ijoints[i]->get_inboard_link();
    JointPtr parent = inboard->get_inner_joint_implicit();
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
bool CRBAlgorithm::factorize_cholesky(MatrixN& M)
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

/// Calculates the generalized inertia of this body
void CRBAlgorithm::calc_generalized_inertia(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  SAFESTATIC MatrixN H;
  SAFESTATIC SMatrix6N K, Is;
  SAFESTATIC vector<SpatialRBInertia> Ic;
  const unsigned SPATIAL_DIM = 6;

  // get the appropriate M
  MatrixN& M = this->_M;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // get implicit joints
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // compute the joint space inertia
  calc_joint_space_inertia(body, rftype, H, Ic);

  // get the number of base degrees-of-freedom
  const unsigned n_base_DOF = (body->is_floating_base()) ? 6 : 0;

  // resize M
  M.resize(n_base_DOF + body->num_joint_dof_implicit(), n_base_DOF + body->num_joint_dof_implicit());

  // set appropriate part of H
  M.set_sub_mat(n_base_DOF, n_base_DOF, H);

  // see whether we are done
  if (!body->is_floating_base())
    return;

  // ************************************************************************
  // floating base: compute K
  // ************************************************************************

  // compute K
  K.resize(6, body->num_joint_dof_implicit());
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the spatial axes for the joint
    const SMatrix6N& s = ijoints[i]->get_spatial_axes(rftype);

    // get the index for this joint
    unsigned jidx = ijoints[i]->get_coord_index() - SPATIAL_DIM;

    // get the outboard link for this joint and its index
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // compute the requisite columns of K
    Ic[oidx].mult(s, Is);
    K.set_sub_mat(0, jidx, Is);
  }

  // get composite inertia in matrix form
  double Ic0[SPATIAL_DIM*SPATIAL_DIM];
  Ic.front().to_matrix(Ic0);

  // setup the remainder of the augmented inertia matrix
  for (unsigned i=0, k=0; i< SPATIAL_DIM; i++)
    for (unsigned j=0; j< SPATIAL_DIM; j++)
      M(j,i) = Ic0[k++];
  M.set_sub_mat(0,6,K);
  MatrixN KT = K.transpose();
  M.set_sub_mat(6,0,KT);

  FILE_LOG(LOG_DYNAMICS) << "(unpermuted) [Ic0 K; K' H]: " << std::endl << M;

  // swap first and second three rows if body is a floating base
  if (body->is_floating_base())
    for (unsigned i=0; i< 3; i++)
      for (unsigned j=0; j< M.rows(); j++)
        std::swap(M(i,j), M(i+3, j));
}

/// Calculates the generalized inertia matrix for the given representation
void CRBAlgorithm::calc_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M)
{
  // do the precalculation
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;
  precalc(body, rftype);

  // calculate the generalized inertia
  if (gctype == DynamicBody::eAxisAngle)
    calc_generalized_inertia_axisangle(M);
  else
  {
    assert(gctype == DynamicBody::eRodrigues);
    calc_generalized_inertia_rodrigues(M);
  }
}

/// Computes *just* the joint space inertia matrix
void CRBAlgorithm::calc_joint_space_inertia(RCArticulatedBodyPtr body, ReferenceFrameType rftype, MatrixN& H, vector<SpatialRBInertia>& Ic) const
{
  SAFESTATIC MatrixN tmp, sub;
  SAFESTATIC vector<SMatrix6N> forces;
  SAFESTATIC vector<vector<bool> > supports;
  queue<RigidBodyPtr> link_queue;
  const unsigned SPATIAL_DIM = 6;

  // get the sets of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();
  const vector<JointPtr>& joints = body->get_joints();

  // set the composite inertias intially to the spatial isolated inertias
  Ic.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Ic[idx] = links[i]->get_spatial_iso_inertia(rftype);
  }

  // ************************************************************************
  // first, determine the supports for the joints and the number of joint DOF
  // ************************************************************************
 
  // create and initialize the supports array
  supports.resize(joints.size());
  for (unsigned i=0; i< joints.size(); i++)
  {
    supports[i].resize(links.size());
    for (unsigned j=0; j< links.size(); j++)
      supports[i][j] = false;
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

    // get the implicit inner joint for this link
    JointPtr joint(link->get_inner_joint_implicit());
    unsigned jidx = joint->get_index();
    assert(joint);

    // add this link to the support for the joint
    supports[jidx][link->get_index()] = true;

    // add all supports from the outer joints of this link
    list<RigidBodyPtr> child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr child, child_links)
    {
      // don't process children with lower link indices (loops)
      if (child->get_index() < link->get_index())
        continue;

      // get the inner implicit joint
      JointPtr child_joint = child->get_inner_joint_implicit();
      unsigned jiidx = child_joint->get_index();

      // setup the supports
      for (unsigned i=0; i< links.size(); i++)
        if (supports[jiidx][i])
          supports[jidx][i] = true;
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
    for (unsigned i=0; i< supports.size(); i++)
    {
      std::ostringstream str;
      for (unsigned j=0; j< supports[i].size(); j++)
        str << supports[i][j] << " ";
      FILE_LOG(LOG_DYNAMICS) << i << ": " << str.str() << endl;
    }
  }

  // resize H 
  H.set_zero(body->num_joint_dof_implicit(), body->num_joint_dof_implicit());

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
      unsigned im1 = parent->get_index();
    
      // add this inertia to its parent
      if (rftype == eGlobal)
        Ic[im1] += Ic[i];
      else
      {
        assert(rftype == eLink);
        SpatialTransform X_im1_i = link->get_spatial_transform_backward();
        Ic[im1] += X_im1_i.transform(Ic[i]);
      }

      FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (child) link " << link->id << ": " << std::endl << Ic[i];
      FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (parent) link " << parent->id << ": " << std::endl << Ic[im1];
    }

    // indicate that the link has been processed
    body->_processed[i] = true;
  }

  // ************************************************************************
  // compute H
  // ************************************************************************

  // compute the forces
  forces.resize(links.size());
  for (unsigned i=0; i < ijoints.size(); i++)
  {
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link(); 
    unsigned oidx = outboard->get_index();
    const SMatrix6N& si = ijoints[i]->get_spatial_axes(rftype);
    Ic[oidx].mult(si, forces[oidx]);
  } 

  // if this is a floating base, we must subtract from the coordinate indices
  unsigned add = (body->is_floating_base()) ? SPATIAL_DIM : 0;

  // setup H
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the starting coordinate index for this joint
    unsigned iidx = ijoints[i]->get_coord_index() - add;

    // get the outboard link for joint i
    RigidBodyPtr outboardi = ijoints[i]->get_outboard_link();
    unsigned oiidx = outboardi->get_index();

    // get the spatial axes for jointi
    const SMatrix6N& si = ijoints[i]->get_spatial_axes(rftype);

    // compute the H term for i,i
    si.transpose_mult(forces[oiidx], sub);

    // set the appropriate part of H
    H.set_sub_mat(iidx,iidx,sub);

    // determine what will be the new value for m
    for (unsigned j=i+1; j< ijoints.size(); j++)
    {
      // get the outboard link for joint j
      RigidBodyPtr outboardj = ijoints[j]->get_outboard_link();
      unsigned ojidx = outboardj->get_index();

      // if link j is not supported by joint i, contribution to H is zero
      if (!supports[ijoints[i]->get_index()][ojidx])
        continue;

      // get the starting coordinate index for joint j
      unsigned jidx = ijoints[j]->get_coord_index() - add;

      // compute the appropriate submatrix of H
      if (rftype == eGlobal)
        si.transpose_mult(forces[ojidx], sub);
      else
      {
        SpatialTransform X_i_j(outboardj->get_transform(), outboardi->get_transform());
        si.transpose_mult(X_i_j.transform(forces[ojidx]), sub);
      }

      // set the appropriate parts of H
      H.set_sub_mat(iidx,jidx,sub);
      H.set_sub_mat(jidx,iidx,sub, true);
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "joint space inertia: " << endl << H;
}

/// Calculates the generalized inertia matrix for the given representation
void CRBAlgorithm::calc_generalized_inertia_axisangle(MatrixN& M) const
{
  SAFESTATIC SMatrix6N K, sb, Is;
  SAFESTATIC MatrixN H;
  SAFESTATIC MatrixN KT;
  SAFESTATIC vector<SpatialRBInertia> Ic;
  const unsigned SPATIAL_DIM = 6;

  // get the set of links
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // get the joint space inertia and composite inertias
  calc_joint_space_inertia(body, rftype, H, Ic);

  // get the number of base degrees-of-freedom
  const unsigned n_base_DOF = (body->is_floating_base()) ? 6 : 0;

  // resize M and set H
  M.resize(n_base_DOF + body->num_joint_dof_implicit(), n_base_DOF + body->num_joint_dof_implicit());
  M.set_sub_mat(n_base_DOF, n_base_DOF, H);

  // look for simplest case
  if (!body->is_floating_base())
    return;
 
  // ************************************************************************
  // floating base: transform all inertias 
  // ************************************************************************
  if (rftype == eLink)
  {
    for (unsigned i=0; i< Ic.size(); i++)
    { 
      unsigned idx = links[i]->get_index();
      const Matrix4& linkX = links[i]->get_transform();
      Matrix4 baseX(&IDENTITY_3x3, &links.front()->get_position());
      SpatialTransform X_b_i(linkX, baseX);
      Ic[idx] = X_b_i.transform(Ic[idx]);
    }
  }
  else
  {
    assert(rftype == eGlobal);
    for (unsigned i=0; i< Ic.size(); i++)
    { 
      unsigned idx = links[i]->get_index();
      Matrix4 baseX(&IDENTITY_3x3, &links.front()->get_position());
      SpatialTransform X_b_0(IDENTITY_4x4, baseX);
      Ic[idx] = X_b_0.transform(Ic[idx]);
    }
  }

  // ************************************************************************
  // floating base: compute K
  // ************************************************************************

  // compute K
  K.resize(n_base_DOF, body->num_joint_dof_implicit());
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the joint index and spatial axes
    unsigned jidx = ijoints[i]->get_coord_index() - SPATIAL_DIM;
    const SMatrix6N& si = ijoints[i]->get_spatial_axes(rftype);
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link();
    unsigned idx = outboard->get_index(); 

    if (rftype == eLink)
    {
      const Matrix4& linkX = outboard->get_transform();
      Matrix4 baseX(&IDENTITY_3x3, &links.front()->get_position());
      SpatialTransform X_b_i(linkX, baseX);
      X_b_i.transform(si, sb);
      Ic[idx].mult(sb, Is);
    }
    else
    {
      assert(rftype == eGlobal);
      Matrix4 baseX(&IDENTITY_3x3, &links.front()->get_position());
      SpatialTransform X_b_0(IDENTITY_4x4, baseX);
      X_b_0.transform(si, sb);
      Ic[idx].mult(sb, Is);
    }

    K.set_sub_mat(0, jidx, Is);
  }

  // get composite inertia in matrix form
  double Ic0[SPATIAL_DIM*SPATIAL_DIM];
  Ic.front().to_matrix(Ic0);

  // setup the remainder of the augmented inertia matrix
  for (unsigned i=0, k=0; i< SPATIAL_DIM; i++)
    for (unsigned j=0; j< SPATIAL_DIM; j++)
      M(j,i) = Ic0[k++];
  M.set_sub_mat(0,6,K);
  SMatrix6N::transpose(K, KT);
  M.set_sub_mat(6,0,KT);

  FILE_LOG(LOG_DYNAMICS) << "(unpermuted) [Ic0 K; K' H]: " << std::endl << M;

  // swap first and second three columns if body is a floating base
  if (body->is_floating_base())
    for (unsigned i=0; i< 3; i++)
      for (unsigned j=0; j< M.rows(); j++)
        std::swap(M(j,i), M(j, i+3));

  FILE_LOG(LOG_DYNAMICS) << "(permuted) [Ic0 K; K' H]: " << std::endl << M;
}

/// Calculates the generalized inertia of this body with respect to Rodrigues parameters
/**
 * \note this method does not utilize cached data nor does it cache any data
 *       to speed repeated calculations.
 */
void CRBAlgorithm::calc_generalized_inertia_rodrigues(MatrixN& M) const
{
  SAFESTATIC MatrixN K, K2, L, tmp, tmp1, tmp2;
  SAFESTATIC vector<SpatialRBInertia> Ic;
  SAFESTATIC MatrixN Ic0_7, H;
  SAFESTATIC SMatrix6N Is;
  const unsigned SPATIAL_DIM = 6;

  // get the set of links
  RCArticulatedBodyPtr body(_body);
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // get the joint space inertia and composite inertias
  calc_joint_space_inertia(body, eLink, H, Ic);

// TODO: need to convert everything to the proper frame (see axis-angle version)
  throw std::runtime_error("Need conversion!");

  // get the number of base degrees-of-freedom
  const unsigned n_base_DOF = (body->is_floating_base()) ? 7 : 0;

  // resize M and set H
  M.resize(n_base_DOF + body->num_joint_dof_implicit(), n_base_DOF + body->num_joint_dof_implicit());
  M.set_sub_mat(n_base_DOF, n_base_DOF, H);

  // look for simplest case
  if (!body->is_floating_base())
    return;
 
  // convert six dimensional spatial inertia to seven dimensional
  to_spatial7_inertia(Ic.front(), links.front()->get_orientation(), Ic0_7);

  // ************************************************************************
  // floating base: compute K
  // ************************************************************************

  // compute K
  K.resize(7, body->num_joint_dof_implicit());
  for (unsigned i=0; i< K.columns(); i++)
    K(6,i) = (double) 0.0;
  K2.resize(body->num_joint_dof_implicit(), 7);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    // get the joint coordinate index
    unsigned jidx = ijoints[i]->get_coord_index() - SPATIAL_DIM;

    // get the outboard link and its index
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // compute the requisite column of K
    const SMatrix6N& si = ijoints[i]->get_spatial_axes(eLink);
    Ic[oidx].mult(si, Is);
    K.set_sub_mat(0, jidx, Is);

    // compute the requisite row of K2
    Is.get_sub_mat(0, 3, 0, Is.columns(), tmp1);
    K2.set_sub_mat(jidx, 4, tmp1, true);
    outboard->get_orientation().determine_L(L);
    L*= (double) 2.0;
    Is.get_sub_mat(3, 6, 0, Is.columns(), tmp1);
    L.transpose_mult(tmp1, tmp2);
    K2.set_sub_mat(jidx, 0, tmp2, true);
  }

  // setup the remainder of the augmented inertia matrix
  M.set_sub_mat(0,0,Ic0_7);
  M.set_sub_mat(0,7,K);
  M.set_sub_mat(7,0,K2);

  FILE_LOG(LOG_DYNAMICS) << "(unpermuted) [Ic0 K; K' H]: " << std::endl << M;

  // swap columns
  MatrixN subL = M.get_sub_mat(0, M.rows(), 0, 4);
  MatrixN subR = M.get_sub_mat(0, M.rows(), 4, 7);  
  M.set_sub_mat(0, 0, subR);
  M.set_sub_mat(0, 3, subL);
}

void CRBAlgorithm::to_spatial7_inertia(const SpatialRBInertia& I, const Quat& q, MatrixN& I7)
{
  SAFESTATIC MatrixN work;
  SAFESTATIC MatrixN work2, L;

  // first resize I7
  I7.resize(7,7);

  // now get the subcomponents of I 
  Matrix3 IUL = Matrix3::skew_symmetric(-I.h);
  Matrix3 IUR = IDENTITY_3x3*I.m;
  Matrix3 ILL = I.J;
  Matrix3 ILR = -IUL;

  // copy the invariant parts
  I7.set_sub_mat(0, 4, IUR);
  I7.set_sub_mat(3, 4, ILR);

  // transform and set the non-invariant parts
  q.determine_L(L) *= (double) 2.0;
  work.resize(3,3);
  work.set_sub_mat(0, 0, IUL);
  work.mult(L, work2);
  I7.set_sub_mat(0, 0, work2);
  work.set_sub_mat(0, 0, ILL);
  work.mult(L, work2);
  I7.set_sub_mat(3, 0, work2);

  // set the last row
  I7(6,0) = q.w;
  I7(6,1) = q.x;
  I7(6,2) = q.y;
  I7(6,3) = q.z;
  I7(6,4) = I7(6,5) = I7(6,6) = (double) 0.0;
}

/// Sets all spatial velocities
/**
 * This method sets spatial velocities- it does not compute them recursively
 * as done in Featherstone.  Rather, it sets the velocities from the link
 * linear and angular velocities.
 */
void CRBAlgorithm::set_spatial_velocities(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // get the vector of spatial velocities
  vector<SVector6>& v = this->_velocities;

  // init the vector of spatial velocities
  v.resize(links.size());

  // determine spatial velocities
  for (unsigned i=0; i< links.size(); i++)
    v[links[i]->get_index()] = links[i]->get_spatial_velocity(rftype);
}

/// Performs necessary pre-computations for computing accelerations or applying impulses
void CRBAlgorithm::precalc(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  // get the last reference frame type for the body
  if (_rftype != rftype)
  {
    _rftype = rftype;
    _position_data_valid = false;
    _velocity_data_valid = false;
  }

  // get the links and joints for the body
  const vector<JointPtr>& joints = body->get_implicit_joints();

  // invalidate spatial values for all joints, if necessary
  if (body->positions_invalidated())
  {
    for (unsigned i=0; i< joints.size(); i++)
      joints[i]->reset_spatial_axis();
    body->validate_positions();
  }

  // compute link velocities, if necessary
  if (!_velocity_data_valid)
  {
    set_spatial_velocities(body, rftype);
    _velocity_data_valid = true;
  }

  // compute spatial isolated inertias and generalized inertia matrix
  if (!_position_data_valid)
  {
    // do the calculations
    calc_generalized_inertia(body, rftype);

    // attempt to do a Cholesky factorization of M
    MatrixN& fM = this->_fM;
    MatrixN& M = this->_M;
    fM.copy_from(M);
    if (_rank_deficient = !LinAlg::factor_chol(fM))
    {
      fM.copy_from(M);
//    if (_rank_deficient = !factorize_cholesky(fM))
      try
      {
        LinAlg::pseudo_inverse(fM, LinAlg::svd1);
      }
      catch (NumericalException e)
      {
        fM.copy_from(M);
        LinAlg::pseudo_inverse(fM, LinAlg::svd2);
      }
    }

    // validate position data
    _position_data_valid = true;
  }
}

/// Executes the composite rigid-body method
void CRBAlgorithm::calc_fwd_dyn()
{
  // get the body and the reference frame
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // do necessary pre-calculations
  precalc(body, rftype);

  // execute the appropriate algorithm
  if (body->is_floating_base())
    calc_fwd_dyn_floating_base(body, rftype);
  else
    calc_fwd_dyn_fixed_base(body, rftype);
   
  // update the link accelerations
  update_link_accelerations(body, rftype);
}

/// Solves for acceleration using the body inertia matrix
VectorN& CRBAlgorithm::M_solve(const VectorN& v, VectorN& result) 
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;
  precalc(body, rftype);

  return M_solve_noprecalc(v, result); 
}

/// Solves for acceleration using the body inertia matrix
MatrixN& CRBAlgorithm::M_solve(const MatrixN& m, MatrixN& result)
{
  // do necessary pre-calculations
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;
  precalc(body, rftype);

  return M_solve_noprecalc(m, result); 
}

/// Solves for acceleration using the body inertia matrix
VectorN& CRBAlgorithm::M_solve_noprecalc(const VectorN& v, VectorN& result) const
{
  // get the factorized inertia matrix
  const MatrixN& fJ = this->_fM; 

  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
  {
    // matrix is rank deficient, use pseudo-inverse solution 
    fJ.mult(v, result);
  }
  else
  {
    // matrix is not rank deficient, use Cholesky factorization
    result.copy_from(v);
    LinAlg::solve_chol_fast(fJ, result);
  }

  return result;
}

/// Solves for acceleration using the body inertia matrix
MatrixN& CRBAlgorithm::M_solve_noprecalc(const MatrixN& m, MatrixN& result) const
{
  // get the factorized inertia matrix
  const MatrixN& fJ = this->_fM; 

  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
  {
    // matrix is rank deficient, use pseudo-inverse solution 
    fJ.mult(m, result);
  }
  else
  {
    // matrix is not rank deficient, use Cholesky factorization
    result.copy_from(m);
    LinAlg::solve_chol_fast(fJ, result);
  }

  return result;
}

/// Executes the composite rigid-body method on an articulated body with a fixed base
void CRBAlgorithm::calc_fwd_dyn_fixed_base(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  SAFESTATIC VectorN C, Q, Qi;

  // get the set of links and joints for the articulated body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // ***********************************************************************
  // first, calculate C
  // ***********************************************************************

  // call inverse dynamics to calculate C
  SVector6 f0;
  calc_generalized_forces(f0, C);
  
  // get the number of degrees-of-freedom
  unsigned nDOF = C.size();

  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************
  Q.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    ijoints[i]->get_scaled_force(Qi);
    Qi += ijoints[i]->ff;
    unsigned j = ijoints[i]->get_coord_index();
    Q.set_sub_vec(j, Qi);
  }

  FILE_LOG(LOG_DYNAMICS) << "H: " << std::endl << _M;
  FILE_LOG(LOG_DYNAMICS) << "C: " << C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "Q: " << Q << std::endl;

  // subtract C from Q
  Q -= C;

  // get the pointer to the joint-space acceleration vector
  VectorN& qdd = this->_qdd;

  // compute joint accelerations
  M_solve_noprecalc(Q, qdd);

  FILE_LOG(LOG_DYNAMICS) << "qdd: " << qdd << std::endl;

  // set qdd
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned j = ijoints[i]->get_coord_index();
    qdd.get_sub_vec(j, j+ijoints[i]->num_dof(), ijoints[i]->qdd);
  }

  // set spatial acceleration of the base
  this->_a0 = ZEROS_6; 
  links.front()->set_laccel(ZEROS_3);
  links.front()->set_aaccel(ZEROS_3);
}

/// Executes the composite rigid-body method on an articulated body with a floating base
/**
 * This algorithm is taken from [Featherstone, 1987], p. 123.  This is only
 * calculated in the global frame.
 */
void CRBAlgorithm::calc_fwd_dyn_floating_base(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  SAFESTATIC VectorN Qi, b, augV, C, Q;
  const unsigned SPATIAL_DIM = 6;

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() entered" << std::endl;

  // get the set of links and implicit joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();

  // calculate C
  SVector6 f0;
  calc_generalized_forces(f0, C);

  // get the number of degrees-of-freedom
  unsigned nDOF = C.size();
  
  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************

  Q.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    ijoints[i]->get_scaled_force(Qi);
    Qi += ijoints[i]->ff;
    unsigned j = ijoints[i]->get_coord_index() - SPATIAL_DIM;
    Q.set_sub_vec(j, Qi);
  }

  // setup the simulataneous equations to solve: [Featherstone, 1987], eq. 7.24
  VectorN::concat(-f0, Q-= C, b);

  // swap first three and next three elements of b
  for (unsigned i=0; i< 3; i++)
    std::swap(b[i], b[i+3]);

  FILE_LOG(LOG_DYNAMICS) << "link + external forces on base: " << f0 << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "Q: " << Q << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "C: " << C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "b: " << b << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "M: " << std::endl << this->_M;
  
  // solve for accelerations
  M_solve_noprecalc(b, augV); 

  // get pointers to a0 and qdd vectors
  SVector6& a0 = this->_a0;
  VectorN& qdd = this->_qdd;

  // get out a0, qdd
  a0 = augV.get_sub_vec(0,6);
  qdd.resize(augV.size()-6);
  if (augV.size() > 6)
    augV.get_sub_vec(6,augV.size(), qdd);

  // set the base acceleration
  links.front()->set_spatial_accel(this->_a0, rftype);

  FILE_LOG(LOG_DYNAMICS) << "base spatial acceleration: " << this->_a0 << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "joint accelerations:";

  // set qdd
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned j = ijoints[i]->get_coord_index() - SPATIAL_DIM;
    qdd.get_sub_vec(j, j+ijoints[i]->num_dof(), ijoints[i]->qdd);
  }

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() exited" << std::endl;
}

/// Computes the vector "C", used for forward dynamics computation, using the recursive Newton-Euler algorithm
/**
 * \return the spatial vector of forces on the base, which can be ignored for forward dynamics for fixed bases
 */
void CRBAlgorithm::calc_generalized_forces(SVector6& f0, VectorN& C)
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC VectorN Q;
  SAFESTATIC vector<SVector6> forces, a;
  SAFESTATIC vector<SpatialRBInertia> Iiso;
  queue<RigidBodyPtr> link_queue;
  SVector6 tmp1, tmp2;

  // get the body and the reference frame
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;

  // check existing data
  if (rftype != _rftype)
  {
    _position_data_valid = false;
    _velocity_data_valid = false;
    _rftype = rftype;
  }

  // get the set of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& ijoints = body->get_implicit_joints();
  if (links.empty())
  {
    C.resize(0);
    f0 = ZEROS_6;
    return;
  }

  // compute link velocities if necessary
  if (!_velocity_data_valid)
  {
    set_spatial_velocities(body, rftype);
    _velocity_data_valid = true;
  }

  // determine spatial isolated inertias
  Iiso.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Iiso[idx] = links[i]->get_spatial_iso_inertia(rftype);
  }

  // **************************************************************************
  // first, compute forward dynamics using RNE algorithm; we copy the algorithm
  // here, because we need some data out of it
  // **************************************************************************

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces entered" << std::endl;

  // ** STEP 1: compute accelerations

  // get the link velocities
  vector<SVector6>& v = this->_velocities;

  // get the base link
  RigidBodyPtr base = links.front();

  // setup the map of link accelerations
  a.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
    a[i] = ZEROS_6;
  
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
    unsigned idx = link->get_index();
    body->_processed[idx] = true;

    // push all children of the link onto the queue, unless they were already
    // processed  
    list<RigidBodyPtr> child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the link's parent
    RigidBodyPtr parent(link->get_parent_link());

    // get the joint for this link
    JointPtr joint(link->get_inner_joint_implicit());

    // get the spatial link velocity
    const SVector6& vx = v[idx]; 

    // get the reference to the spatial link acceleration
    SVector6& ax = a[idx];
 
    // get spatial axes for this link's inner joint
    const SMatrix6N& s = joint->get_spatial_axes(rftype);

    // get derivative of the spatial axes for this link's inner joint
    const SMatrix6N& s_dot = joint->get_spatial_axes_dot(rftype);

    // get the current joint velocity
    const VectorN& qd = joint->qd;

    // **** compute acceleration

    // add this link's contribution
    ax += SVector6::spatial_cross(vx, s.mult(qd)) + s_dot.mult(qd);

    // now add parent's contribution
    if (rftype == eGlobal)
      ax += a[parent->get_index()];
    else
    {
      SpatialTransform X_i_im1 = link->get_spatial_transform_forward();
      ax += X_i_im1.transform(a[parent->get_index()]);
    }

    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  spatial axis: " << s << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  spatial joint velocity: " << (s.mult(qd)) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << vx << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << ax << std::endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  // use a map to determine which links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // setup a map of link forces, all set to zero initially
  forces.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
    forces[i] = ZEROS_6;

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
    unsigned idx = link->get_index();

    // if this link has already been processed, do not process it again
    if (body->_processed[idx])
      continue;

    // verify all children have been processed
    if (!body->all_children_processed(link))
      continue;

    // get the forces for this link and this link's parent
    SVector6& fi = forces[idx];
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << fi << std::endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * a[idx]) << std::endl;

    // get the spatial velocity for this link
    const SVector6& vx = v[idx];

    // add I*a to the link force
    fi += Iiso[idx] * a[idx];

    // update the determined force to this link w/Coriolis + centrifugal terms
    fi += SVector6::spatial_cross(vx, Iiso[idx] * vx);

    // determine external forces in link frame
    const Vector3& fext = link->sum_forces();
    const Vector3& text = link->sum_torques();
    const Matrix4& T = link->get_transform();
    SVector6 fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << fi << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  sum of external forces / torques: " << fext << " / " << text << std::endl;

    // subtract external forces in the appropriate frame
    if (rftype == eGlobal)
    {
      SpatialTransform X_0_i = link->get_spatial_transform_link_to_global();
      fi -= X_0_i.transform(fx);
    }
    else
      fi -= fx;

    FILE_LOG(LOG_DYNAMICS) << "  external force in link frame: " << fx << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << fi << std::endl;

    // update the parent force and add parent for processing (if parent)
    RigidBodyPtr parent = link->get_parent_link();
    if (parent)
    {
      // put the parent on the queue
      SVector6& fim1 = forces[parent->get_index()];
      if (rftype == eGlobal)
        fim1 += fi;
      else
      {
        assert(rftype == eLink);
        fim1 += link->get_spatial_transform_backward().transform(fi);
      }
      link_queue.push(parent);
    }

    // indicate that this link has been processed
    body->_processed[idx] = true;
  }
  
  // ** STEP 3: compute actuator forces (C)

  // determine the length of the C vector
  const unsigned nDOF = body->num_joint_dof_implicit();

  // determine whether we need to offset coordinates
  const unsigned ADD = body->is_floating_base() ? SPATIAL_DIM : 0;

  // compute actuator forces (C)
  C.resize(nDOF);
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned jidx = ijoints[i]->get_coord_index() - ADD;
    const SMatrix6N& s = ijoints[i]->get_spatial_axes(rftype);
    RigidBodyPtr outboard = ijoints[i]->get_outboard_link();
    unsigned oidx = outboard->get_index();
    s.transpose_mult(forces[oidx], Q);
    C.set_sub_vec(jidx, Q);

    FILE_LOG(LOG_DYNAMICS) << " -- computing C for link " << outboard->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "   -- s: " << std::endl << s;
    FILE_LOG(LOG_DYNAMICS) << "   -- forces: " << forces[oidx] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "   -- component of C: " << Q << std::endl;
  }  

  FILE_LOG(LOG_DYNAMICS) << "------------------------------------------------" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "forces on base: " << forces[0] << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces exited" << std::endl;

  // set forces on base
  Matrix3 R = base->get_transform().get_rotation();
  if (rftype == eGlobal)
  {
    SpatialTransform X(IDENTITY_3x3, ZEROS_3, R, base->get_position());
    f0 = X.transform(forces[0]);
  }
  else
  {
    f0.set_upper(R.transpose_mult(forces[0].get_upper()));
    f0.set_lower(R.transpose_mult(forces[0].get_lower()));
  }
}

/// Updates all link accelerations (except the base)
void CRBAlgorithm::update_link_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype) const
{
  Matrix3 Ri, Rim1;
  queue<RigidBodyPtr> link_queue;
  SVector6 tmp1, tmp2;

  // get the set of links and their velocities
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<SVector6>& v = this->_velocities;

  // if there are no links, there is nothing to do
  if (links.empty())
    return;

  // setup the vector of accelerations
  SAFESTATIC vector<SVector6> a;
  a.resize(links.size());

  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // get the base link and mark it as processed
  RigidBodyPtr base = links.front();
  body->_processed[base->get_index()] = true;
  
  // get the spatial acceleration of the base link (should have already been
  // computed)
  a.front() = this->_a0;

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
    unsigned idx = link->get_index();
    body->_processed[idx] = true;

    // push all unprocessed children of the link onto the queue
    list<RigidBodyPtr> children;
    link->get_child_links(std::back_inserter(children));
    BOOST_FOREACH(RigidBodyPtr rb, children)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the parent link
    JointPtr joint(link->get_inner_joint_implicit());
    RigidBodyPtr parent(link->get_parent_link());
    unsigned pidx = parent->get_index();

    // get the parent acceleration (should be in the desired frame)
    const SVector6& aparent = a[pidx];
  
    // get the link spatial axis
    const SMatrix6N& s = joint->get_spatial_axes(rftype); 

    // get the link velocity
    const SVector6& vx = v[idx];

    // determine the link accel
    a[idx] = aparent;
    a[idx] += SVector6::spatial_cross(vx, s.mult(joint->qd)) + s.mult(joint->qdd);

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent acceleration: " << aparent << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- s: " << joint->get_spatial_axes(rftype) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- velocity: " << vx << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qdd: " << joint->qdd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- acceleration: " << a[idx] << std::endl;
  }

  // set all accelerations
  for (unsigned i=1; i< links.size(); i++)
    links[i]->set_spatial_accel(a[links[i]->get_index()], rftype);

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::update_link_accelerations() exited" << std::endl;
}

/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
void CRBAlgorithm::apply_impulse(const Vector3& jj, const Vector3& jk, const Vector3& point, RigidBodyPtr link)
{
  // get the body
  RCArticulatedBodyPtr body(_body);
  ReferenceFrameType rftype = body->computation_frame_type;

  // do necessary pre-calculations
  precalc(body, rftype); 

  // two different methods, depending on whether body has floating base
  // or not
  if (body->is_floating_base())
    apply_impulse_floating_base(body, jj, jk, point, link);
  else
    apply_impulse_fixed_base(body, jj, jk, point, link);

//*/
  // An alternative method for applying impulses using generalized coordinates
  // below...
/*
  SMatrix6 Xcp = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, link->get_position());
  SMatrix6 X0p = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, ZEROS_3);
  SVector6 jx(jj[0], jj[1], jj[2], jk[0], jk[1], jk[2]);
 
  // convert the impulse to a generalized impulse
  VectorN gf;
  body->convert_to_generalized_force(link, point, jj, jk, gf);

  // get the generalized inertia and invert it
  MatrixN M;
  body->get_generalized_inertia(M);
  LinAlg::inverse_PD(M);

  // get the generalized velocity
  VectorN gv;
  body->get_generalized_velocity(gv);

  // compute the change in velocity
  VectorN dv;
  M.mult(gf, dv);

  SMatrix6 Xc0 = SMatrix6::calc_spatial_transform(IDENTITY_3x3, ZEROS_3, IDENTITY_3x3, link->get_position());
  SMatrix6 X0c = SMatrix6::calc_spatial_transform(IDENTITY_3x3, link->get_position(), IDENTITY_3x3, ZEROS_3);
  MatrixN Mx;
  body->get_generalized_inertia(Mx, eGlobal);
  SMatrix6 M2(Mx.begin());

  FILE_LOG(LOG_DYNAMICS) << "generalized inertia (global frame): " << std::endl << Mx;
  FILE_LOG(LOG_DYNAMICS) << "spatial transform (global to centroidal): " << std::endl << Xc0;
  FILE_LOG(LOG_DYNAMICS) << "[permuted] spatial inertia (centroidal frame): " << std::endl << (Xc0 * M2 * X0c);

  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (global frame): " << (X0p * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (centroidal frame): " << (Xcp * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "base transform: " << std::endl << body->get_links().front()->get_transform();
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
*/
}

/// TODO: fix this for bodies with kinematic loops
/// Applies an impulse to an articulated body with a fixed base; complexity O(n^2)
void CRBAlgorithm::apply_impulse_fixed_base(RCArticulatedBodyPtr body, const Vector3& jj, const Vector3& jk, const Vector3& point, RigidBodyPtr link)
{
  const unsigned OPSPACE_DIM = 6;
  SAFESTATIC MatrixN J, col;
  SAFESTATIC VectorN p, dq, dq_sub, work;

  // get the vector of implicit joints and all joints
  const vector<JointPtr>& ijoints = body->get_implicit_joints();
  const vector<JointPtr>& joints = body->get_joints();

  // this does not work for bodies with kinematic loops
  assert(body->_ejoints.empty());

  // get # of joints
  unsigned n = this->_M.rows();

  // compute the Jacobian with respect to the application point
  J.set_zero(OPSPACE_DIM, n);
  RigidBodyPtr l = link;
  JointPtr j;
  while ((j = l->get_inner_joint_implicit()))
  {
    // compute the Jacobian column(s) for the joint
    body->calc_jacobian_column(j, point, col);

    // set the column(s) of the Jacobian
    const unsigned start = j->get_coord_index();
    J.set_sub_mat(0,start, col);

    // set l to its parent
    l = RigidBodyPtr(l->get_parent_link());
  }

  // form a vector from the impulse
  VectorN::concat(jj, jk, p);

  // compute the change in velocity at the joints
  J.transpose_mult(p,work);
  M_solve_noprecalc(work, dq);

  // apply the change and update link velocities
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned jidx = ijoints[i]->get_coord_index();
    dq.get_sub_vec(jidx, jidx + ijoints[i]->num_dof(), dq_sub);
    ijoints[i]->qd += dq_sub;
  }  
  body->update_link_velocities();

  // reset all force and torque accumulators -- impulses drive them to zero
  const vector<RigidBodyPtr>& links = body->get_links();
  for (unsigned i=0; i< joints.size(); i++)
    joints[i]->reset_force(); 
  for (unsigned i=0; i< links.size(); i++)
    links[i]->reset_accumulators();
}

/// TODO: fix this for bodies with kinematic loops
/// Applies an impulse to an articulated body with a floating base; complexity O(n^2)
void CRBAlgorithm::apply_impulse_floating_base(RCArticulatedBodyPtr body, const Vector3& jj, const Vector3& jk, const Vector3& point, RigidBodyPtr link)
{
  const unsigned OPSPACE_DIM = 6;
  SAFESTATIC SMatrix6N J;
  SAFESTATIC VectorN augV, b, dqd, work;
  SVector6 dv0;

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::apply_impulse_floating_base() entered" << std::endl;

  // does not work for bodies with kinematic loops
  assert(body->_ejoints.empty());

  // get the computation frame type
  ReferenceFrameType rftype = _rftype;

  // get the base link
  RigidBodyPtr base = body->get_base_link();

  // compute the Jacobian for the floating base w.r.t. the contact point
  J.set_zero(OPSPACE_DIM, body->num_joint_dof_implicit());

  // compute the Jacobian with respect to the contact point
  RigidBodyPtr l = link;
  JointPtr j;
  while ((j = l->get_inner_joint_implicit()))
  {
    // compute the Jacobian column(s) for the joint
    const SMatrix6N& columns = j->get_spatial_axes(rftype);

    // set the column(s) of the Jacobian
    const unsigned start = j->get_coord_index();
    J.set_sub_mat(0,start, columns);

    // set l to its parent
    l = RigidBodyPtr(l->get_parent_link());
  }

  // form a vector from the impulse
  SVector6 p(jj, jk);

  // transform p into the desired frame
  SpatialTransform X_f_pt;
  if (rftype == eGlobal)
    X_f_pt = SpatialTransform(IDENTITY_3x3, point, IDENTITY_3x3, ZEROS_3);
  else
  {
    assert(rftype == eLink);
    RigidBodyPtr base = body->get_base_link();
    const Matrix3 R = base->get_transform().get_rotation();
    const Vector3 x = base->get_position();
    X_f_pt = SpatialTransform(IDENTITY_3x3, point, R, x);
  }
  SVector6 pf = X_f_pt.transform(p);

  // form vector to solve for b
  VectorN::concat(pf, J.transpose_mult(pf, work), b);

  // swap top three rows and next three rows of b
  std::swap(b[0], b[3]);
  std::swap(b[1], b[4]);
  std::swap(b[2], b[5]);

  // compute changes in base and joint velocities
  M_solve_noprecalc(b, augV);

  // get change in base and change in joint velocities
  augV.get_sub_vec(0,OPSPACE_DIM, dv0);

  // update the base velocity
  SVector6 basev = base->get_spatial_velocity(rftype);
  base->set_spatial_velocity(basev + SVector6(dv0.begin()), rftype);

  FILE_LOG(LOG_DYNAMICS) << "  Jacobian: " << std::endl << J;
  FILE_LOG(LOG_DYNAMICS) << "  impulse (impulse frame): " << p << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  impulse (last frame): " << pf << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  change in base velocity: " << dv0 << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  new base linear velocity: " << base->get_lvel() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  new base angular velocity: " << base->get_avel() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  change in link velocities: " << augV.get_sub_vec(OPSPACE_DIM, augV.size()) << std::endl;

  // apply the change and update link velocities
  const vector<JointPtr>& ijoints = body->get_implicit_joints();
  for (unsigned i=0; i< ijoints.size(); i++)
  {
    unsigned idx = ijoints[i]->get_coord_index();
    augV.get_sub_vec(idx, idx+ijoints[i]->num_dof(), dqd);
    FILE_LOG(LOG_DYNAMICS) << " joint " << ijoints[i]->id << " qd: " << ijoints[i]->qd << "  dqd: " << dqd << std::endl;  
    ijoints[i]->qd += dqd;
  }  
  body->update_link_velocities();

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::apply_impulse_floating_base() exited" << std::endl;
}


