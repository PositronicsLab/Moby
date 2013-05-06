/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <queue>
#include <map>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/LinAlg.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/NumericalException.h>
#include <Moby/FSABAlgorithm.h>
#include <Moby/FastThreadable.h>

using namespace Ravelin;
using namespace Moby;
using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;

FSABAlgorithm::FSABAlgorithm()
{
  _position_data_valid = false;
}

/// Calculates the inverse generalized inertia matrix
void FSABAlgorithm::calc_inverse_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixNd& iM)
{
  vector<Twistd> sprime;

  // we currently only handle axis-angle generalized coordinates...
  assert(gctype == DynamicBody::eAxisAngle);

  // get the body
  RCArticulatedBodyPtr body(_body);

  // prepare to calculate the impulses
  calc_impulse_dyn(body, eGlobal);

  // now, compute s'I
  const vector<RigidBodyPtr>& links = body->get_links();
  vector<Wrenchd> sTI(links.size());
  for (unsigned i=1; i< links.size(); i++)
  {
    JointPtr joint = links[i]->get_inner_joint_implicit();
    if (links[i]->get_computation_type() == eJoint)
    {
      const vector<Twistd>& s = joint->get_spatial_axes();
      transform(s.front().pose, links[i]->get_computation_pose(), sprime); 
    }
    else
      sprime = joint->get_spatial_axes();
    transpose_mult(sprime, _I[i], sTI[i]);
  }

  // get the number of generalized coords
  const unsigned NGC = body->num_generalized_coordinates(gctype);

  // resize inv(M)
  iM.resize(NGC, NGC);
  _workv.set_zero(NGC);

  // compute the inverse inertia
  #pragma omp parallel for
  for (unsigned i=0; i< NGC; i++)
  {
    // setup the generalized impulse
    _workv.set_zero(NGC);
    _workv[i] = (double) 1.0;

    // apply the generalized impulse
    apply_generalized_impulse(i, sTI, _workv);

    // set the appropriate column of inv(M)
    iM.set_column(i, _workv);
  } 

  // make inverse(M) symmetric
  const unsigned LD = NGC;
  double* data = iM.data();
  for (unsigned i=0, ii=0; i< NGC; i++, ii+= LD)
    for (unsigned j=0, jj=0; j< i; j++, jj+= LD)
      data[ii] = data[jj];
  
  FILE_LOG(LOG_DYNAMICS) << "inverse M: " << std::endl << iM;
}

/// Applies a generalized impulse using the algorithm of Drumwright
void FSABAlgorithm::apply_generalized_impulse(unsigned index, const vector<Wrenchd>& sTI, VectorNd& vgj) const
{
  queue<RigidBodyPtr> link_queue;
  const unsigned SPATIAL_DIM = 6;
  const unsigned BASE_IDX = 0;

  // get the body as well as sets of links and joints
  RCArticulatedBodyPtr body(_body);
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_implicit_joints();

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = (body->is_floating_base()) ? SPATIAL_DIM : 0;

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << vgj << endl;

  // TODO: is this necessary?
  // get the computation reference frame
  assert(_rftype == eGlobal);
  // TODO: can we verify that things are all in the same frame?

  // clear values for vectors
  _Y.resize(links.size());
  _processed.resize(links.size());

  // reset Y and processed
  for (unsigned i=0; i< links.size(); i++)
  {
    _Y[i].set_zero();
    _processed[i] = false;
  }

  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (i > BASE_IDX && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    unsigned idx = link_pqueue.top();
    RigidBodyPtr link = links[idx];
    link_pqueue.pop();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (_processed[idx])
      continue;
  
    // indicate that this link has been processed
    _processed[idx] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->get_index() == BASE_IDX)
    {
      link_pqueue.push(parent->get_index());
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    unsigned pidx = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_implicit());
    const vector<Twistd>& s = joint->get_spatial_axes();

    // get I, Y, and mu
    const SpatialABInertiad& I = _I[idx];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && parent->get_index() == BASE_IDX)
      continue;
 
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    vgj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu << endl;
    
    // compute the qm subexpression
    _mu -= transpose_mult(s, _Y[idx], _workv2);

    // get Is
    const vector<Wrenchd>& Is = _Is[idx];

    // prepare to update parent Y
    solve_sIs(idx, _mu, _sIsmu);
    _Y[pidx] += _Y[idx] + mult(Is, _sIsmu);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Is: " << Is << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[idx] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    parent Y is now: " << endl << _Y[pidx] << endl;
  }

  // TODO: fix this
  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // get the linear and angular components (global-aligned frame centered
    // at base)
    RigidBodyPtr base = links.front();
    const Twistd& v = base->velocity();
    Vector3 linear(vgj[0], vgj[1], vgj[2]);
    Vector3 angular(vgj[3], vgj[4], vgj[5]);

    // compute Y[0] update (global frame)
    Vector3 cross = Vector3::cross(-base->get_position(), -linear);
    _Y.front() += SVector6(-linear, -angular - cross);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
    FILE_LOG(LOG_DYNAMICS) << endl;
    for (unsigned i=st_idx; i< links.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << "Articulated zero-velocity delta vector for " << links[i]->id << ": " << _Y[links[i]->get_index()] << endl;  
    FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_velocity_deltas() ended" << endl;
  }

  // ************************************************************
  // determine the changes in joint and link velocities
  // ************************************************************

  // get the base link
  RigidBodyPtr base = links.front();
  const unsigned NUM_LINKS = links.size();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();

    // determine change in generalized velocities
    for (unsigned i=0; i< N_BASE_GC; i++)
      vgj[i] = _dv.front()[i];
  }
  else 
    _dv.front().set_zero();

  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // get the parent link
    RigidBodyPtr parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    JointPtr joint(link->get_inner_joint_implicit());

    // check the starting index of this joint
    const unsigned CSTART = joint->get_coord_index();
    if (CSTART > index)
      continue; 
    
    // get spatial axes of the inner joint for link i
    const vector<Twistd>& s = joint->get_spatial_axes();
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
//    FILE_LOG(LOG_DYNAMICS) << "    -- s: " << s << endl;
    
    // determine appropriate components of gj
    vgj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu);

    // transform s to Y's frame
    transform(s.front().pose, _Y[i].pose, s, sprime); 

    // compute s'Y
    transpose_mult(sprime, _Y[i], _sTY);

    // determine the joint and link velocity updates
    sTI[i].mult(_dv[h], _workv3) += _sTY;
    _workv2.negate() += _mu;
    solve_sIs(i, _workv2, _qd_delta);
    _dv[i] = _dv[h] + mult(s, _qd_delta);

    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- I: " << endl << _I[i];
    FILE_LOG(LOG_DYNAMICS) << "    -- s'I: " << endl << sTI[i];
    FILE_LOG(LOG_DYNAMICS) << "    -- Qi: " << _mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- Qi - (s'I*dv + s'Y): " << endl << _workv2;
    FILE_LOG(LOG_DYNAMICS) << "    -- [Qi - (s'I*dv + s'Y)]/s'Is = qd_delta: " << _qd_delta << endl; 
    FILE_LOG(LOG_DYNAMICS) << "    -- dv[parent]: " << _dv[h] << endl;    
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);

    // update vgj
    for (unsigned k=0; k< joint->num_dof(); k++)
      vgj[CSTART+k] = _qd_delta[k];
  }

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() exited" << endl;
}

/// Applies a generalized impulse using the algorithm of Drumwright
void FSABAlgorithm::apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gj)
{
  SAFESTATIC VectorNd tmp, tmp2;
  SAFESTATIC vector<VectorNd> mu;
  queue<RigidBodyPtr> link_queue;

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = (gctype == DynamicBody::eAxisAngle) ? 6 : 7;

  // get the body
  RCArticulatedBodyPtr body(_body);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << gj << endl;

   // get the computation reference frame
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // prepare to calculate the impulse
  calc_impulse_dyn(body, rftype);

  // get the sets of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_implicit_joints();

  // init mu and Y 
  _mu.resize(links.size());
  _Y.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
    _Y[i].set_zero();

  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
 
  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    unsigned idx = link_pqueue.top();
    RigidBodyPtr link = links[idx];
    link_pqueue.pop();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[idx])
      continue;
 
    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue;
 
    // indicate that this link has been processed
    body->_processed[idx] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    unsigned pidx = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_implicit());
    unsigned jidx = joint->get_index();
    const vector<Twistd>& s = joint->get_spatial_axes();

    // get I
    const SpatialABInertiad& I = _I[idx];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[idx]);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu[idx] << endl;
    
    // compute the qm subexpression
    _mu[idx] -= transpose_mult(s, _Y[idx], tmp2);

    // get Is
    const vector<Wrenchd>& Is = _Is[idx];

    // prepare to update parent Y
    solve_sIs(idx, _mu[idx], _sIsmu);
    Wrenchd uY = _Y[idx] + mult(Is, _sIsmu);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Is: " << Is << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu[idx] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[idx] << endl;

    // get the parent current zero velocity delta and inertia
    Wrenchd& Yim1 = Y[pidx];

    // handle cases of global and link frames separately
    if (rftype == eGlobal)
    {
      Yim1 += uY;
      FILE_LOG(LOG_DYNAMICS) << "    parent Y is now: " << endl << Yim1 << endl;
    }
    else
    {
      // get the backward and forward spatial transforms
      SpatialTransform X_im1_i = link->get_spatial_transform_backward();
    
      // compute transformed Y
      SVector6 tY = X_im1_i.transform(uY);
        
      // update ZVD
      Yim1 += tY;

      FILE_LOG(LOG_DYNAMICS) << "    untransformed Y: " << uY << endl;
      FILE_LOG(LOG_DYNAMICS) << "    transformed Y: " << tY << endl;
      FILE_LOG(LOG_DYNAMICS) << "    parent Y is now: " << endl << Yim1 << endl;
    }
  }

  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // get the linear and angular components
    RigidBodyPtr base = links.front();
    Vector3 linear(gj[0], gj[1], gj[2]);
    Vector3 angular(gj[3], gj[4], gj[5]);

    // compute Y[0] update in link frame
    SVector6 Y0_delta(-linear, -angular);

    // see whether we need to transform Y0
    if (gctype == DynamicBody::eAxisAngle)
    {
      SpatialTransform X_b_0(IDENTITY_3x3, base->get_position(), base->get_pose().get_rotation(), base->get_position());
      Y0_delta = X_b_0.transform(Y0_delta);
    }

    // update Y[0]
    if (rftype == eLink)
      _Y.front() += Y0_delta;
    else
    {
      assert(rftype == eGlobal);
      SpatialTransform X0i = base->get_spatial_transform_link_to_global();
      _Y.front() += X0i.transform(Y0_delta);
    }
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
    FILE_LOG(LOG_DYNAMICS) << endl;
    for (unsigned i=st_idx; i< links.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << "Articulated zero-velocity delta vector for " << links[i]->id << ": " << Y[links[i]->get_index()] << endl;  
    FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_velocity_deltas() ended" << endl;
  }

  // ************************************************************
  // determine the new joint and link velocities
  // ************************************************************

  // get the base link
  RigidBodyPtr base = links.front();
  const unsigned NUM_LINKS = links.size();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    // update the base velocity
    _v.front() += _dv.front();

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base linear velocity: " << base->get_lvel() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  current base angular velocity: " << base->get_avel() << endl;

    // set the base linear and angular velocities
    base->set_spatial_velocity(_v.front(), rftype);

    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << Y.front() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base linear velocity: " << base->get_lvel() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base angular velocity: " << base->get_avel() << endl;
  }
  else 
    _dv.front().set_zero();
  
  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // get the parent link
    RigidBodyPtr parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    JointPtr joint(link->get_inner_joint_implicit());
    const unsigned jidx = joint->get_index();
    
    // get spatial axes of the inner joint for link i
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);
    
    // get articulated body inertia for link i
    const SpatialABInertiad& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- s: " << s << endl;
    
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _Qi);

    // determine the joint and link velocity updates
    if (rftype == eGlobal)
    {
      transpose_mult(s, (I * _dv[h]) + _Y[i], tmp2).negate();
      tmp2 += _Qi;
      solve_sIs(i, tmp2, _qd_delta);
      _dv[i] = _dv[h] + mult(s, _qd_delta);

      FILE_LOG(LOG_DYNAMICS) << "    -- I * dv[parent]: " << _dv[h] << endl;
    }
    else
    {
      SpatialTransform X_i_im1 = link->get_spatial_transform_forward();
      SVector6 xdv = X_i_im1.transform(_dv[h]);
      transpose_mult(s, (I * xdv) + _Y[i], tmp2).negate();
      tmp2 += _Qi;
      solve_sIs(i, tmp2, _qd_delta);
      _dv[i] = xdv + mult(s, _qd_delta);

      FILE_LOG(LOG_DYNAMICS) << "    -- I * X_i_im1 * dv[parent]: " << (I * X_i_im1.transform(_dv[h])) << endl;
    }
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    _v[i] += _dv[i];
    link->set_spatial_velocity(_v[i], rftype);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << "  v: " << _v[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  for (unsigned i=0; i< joints.size(); i++)
    joints[i]->reset_force(); 
  for (unsigned i=0; i< NUM_LINKS; i++)
    links[i]->reset_accumulators();

  // recompute the coriolis vectors so we can validate velocity data
  calc_spatial_coriolis_vectors(body, rftype);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() exited" << endl;
}

/// Sets all spatial isolated inertias
void FSABAlgorithm::set_spatial_iso_inertias(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // clear spatial values for all links
  _Iiso.resize(links.size());

  // determine spatial isolated inertias
  #pragma omp parallel for
  for (unsigned i=0; i< links.size(); i++)
    _Iiso[i] = links[i]->get_spatial_iso_inertia(rftype);
}

/// Computes the combined spatial coriolis / centrifugal forces vectors for the body
void FSABAlgorithm::calc_spatial_coriolis_vectors(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  SVector6 tmp;
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_coriolis_vectors() entered" << endl;

  // get the set of links
  const vector<RigidBodyPtr >& links = body->get_links();

  // clear spatial values for all links
  _c.resize(links.size());

  // process all links except the base
  for (unsigned i=1; i< links.size(); i++)
  {  
    // get the link
    RigidBodyPtr link = links[i];
    unsigned idx = link->get_index();

    // get the link's joint
    JointPtr joint(link->get_inner_joint_implicit());

    // get the spatial velocity in spatial cross-product form
    const Twistd& v = _v[idx];

    // get the spatial axis
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);

    // compute the coriolis vector
    _c[idx] = SVector6::spatial_cross(v, mult(s, joint->qd));

    FILE_LOG(LOG_DYNAMICS) << "processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "v: " << v << endl;
    FILE_LOG(LOG_DYNAMICS) << "s * qdot: " << mult(s, joint->qd) << endl;
    FILE_LOG(LOG_DYNAMICS) << "c: " << _c[idx] << endl;
  }
}

/// Computes articulated body zero acceleration forces used for computing forward dynamics
void FSABAlgorithm::calc_spatial_zero_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  VectorNd tmp, tmp2;
  queue<RigidBodyPtr> link_queue;

  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() entered" << endl;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // clear spatial values for all links
  _Z.resize(links.size());
  _mu.resize(links.size());

  // ************************************************
  // form the isolated z.a. vectors
  // ************************************************

  // process all links (including the base, if it's floating)
  unsigned start_idx = (body->is_floating_base()) ? 0 : 1;

  // process all links except the base
  for (unsigned i=start_idx; i< links.size(); i++)
  {
    // get the link
    RigidBodyPtr link = links[i];
    unsigned idx = link->get_index();

    // get the spatial link velocity
    const Twistd& v = _v[idx];

    // get the isolated spatial inertia for this link
    const SpatialRBInertiad& Iiso = _Iiso[idx]; 

    // set 6-dimensional spatial isolated zero-acceleration vector of link  
    Twistd Z = SVector6::spatial_cross(v, Iiso * v);
    
    // determine external forces in link frame; note that external forces
    // are always applied to the link c.o.m., so at most, they will need to
    // be rotated to coincide with the link
    const Matrix4& T = link->get_pose();
    SVector6 f(T.transpose_mult_vector(link->sum_forces()), T.transpose_mult_vector(link->sum_torques()));

    // subtract external forces in the appropriate frame
    if (rftype == eGlobal)
    {
      SpatialTransform X_0_i = link->get_spatial_transform_link_to_global();
      Z -= X_0_i.transform(f);
    }
    else
      Z -= f;

    // set the spatial zero acceleration (isolated)
    _Z[idx] = Z;

    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso ZA: " << endl << _Z[idx] << endl;
  }
 
  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false; 

  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    unsigned idx = link_pqueue.top();
    RigidBodyPtr link = links[idx];
    link_pqueue.pop();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[idx])
      continue;
  
    // verify that all children have already been processed
    if (!body->all_children_processed(link))
      continue;

    // indicate that this link has been processed
    body->_processed[idx] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
      link_pqueue.push(parent->get_index());
    unsigned pidx = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_implicit());
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);

    // get I, c, and Z
    const SpatialABInertiad& I = _I[idx];
    const Wrenchd& c = _c[idx];
    const Twistd& Z = _Z[idx];
    
    // compute the qm subexpression
    joint->get_scaled_force(_Q);
    _Q += joint->ff;
    _mu[idx] = _Q;
    _mu[idx] -= _Is[idx].transpose_mult(_c[idx], tmp);
    _mu[idx] -= transpose_mult(s, Z, tmp2);

    // get Is
    const SMatrix6N& Is = _Is[idx];

    // get the qm subexpression
    const VectorNd& mu = _mu[idx];
  
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    c: " << c << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Is: " << Is << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Z: " << Z << endl;

    // don't update Z for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // compute a couple of necessary matrices
    solve_sIs(idx, mu, sIsmu);
    Twistd uZ = Z + (I*c) + mult(Is, sIsmu);

    // get the parent current zero acceleration and inertia
    Twistd& Zim1 = _Z[pidx];

    // handle cases of global and link frames separately
    if (rftype == eGlobal)
    {
      Zim1 += uZ;
      FILE_LOG(LOG_DYNAMICS) << "    parent Z is now: " << endl << Zim1 << endl;
    }
    else
    {
      // get the backward and forward spatial transforms
      SpatialTransform X_im1_i = link->get_spatial_transform_backward();
    
      // compute transformed Z
      Twistd tZ = X_im1_i.transform(uZ);
        
      // update ZA
      Zim1 += tZ;

      FILE_LOG(LOG_DYNAMICS) << "    untransformed Z: " << uZ << endl;
      FILE_LOG(LOG_DYNAMICS) << "    transformed Z: " << tZ << endl;
      FILE_LOG(LOG_DYNAMICS) << "    parent Z is now: " << endl << Zim1 << endl;
    }
  }

  FILE_LOG(LOG_DYNAMICS) << endl;
  unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
  for (unsigned i=st_idx; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << "Articulated zero-acceleration vector for " << links[i]->id << ": " << _Z[links[i]->get_index()] << endl;  
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() ended" << endl;
}

/// Computes articulated body inertias used for computing forward dynamics
void FSABAlgorithm::calc_spatial_inertias(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  SMatrix6N tmp;

  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() entered" << endl;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // clear spatial values for all links
  _rank_deficient.resize(links.size());
  _I.resize(links.size());
  _Is.resize(links.size());
  _sIs.resize(links.size());

  // ************************************************
  // form the isolated spatial inertia vectors
  // ************************************************

  // process all links (including the base, if it's floating)
  unsigned start_idx = (body->is_floating_base()) ? 0 : 1;

  // process all links except the base
  for (unsigned i=start_idx; i< links.size(); i++)
  {
    // get the link
    RigidBodyPtr link = links[i];
    unsigned idx = link->get_index();

    // get the isolated spatial inertia for this link
    const SpatialRBInertiad& Iiso = _Iiso[idx]; 

    // set the articulated body inertia for this link to be its isolated
    // spatial inertia (this will be updated in the phase below)
    _I[idx] = Iiso;

    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso inertia: " << endl << _I[idx];
  }
 
  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
 
  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    unsigned idx = link_pqueue.top();
    RigidBodyPtr link = links[idx];
    link_pqueue.pop();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[idx])
      continue;

    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue; 
  
    // indicate that this link has been processed
    body->_processed[idx] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    unsigned pidx = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_implicit());
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);

    // get I
    const SpatialABInertiad& I = _I[idx];
    
    // compute Is
    I.mult(s, _Is[idx]);

    // compute sIs
    s.transpose_mult(_Is[idx], _sIs[idx]);

    // get whether s is rank deficient
    _rank_deficient[idx] = joint->is_singular_config();

    // if the joint is not rank deficient, compute a Cholesky factorization 
    // of sIs
    if (!_rank_deficient[idx])
      LinAlg::factor_chol(_sIs[idx]);
    else
    {
      try
      {
        LinAlg::pseudo_inverse(_sIs[idx], LinAlg::svd1);
      }
      catch (NumericalException e)
      {
        s.transpose_mult(_Is[idx], _sIs[idx]);
        LinAlg::pseudo_inverse(_sIs[idx], LinAlg::svd2);
      }
    }

    // get Is
    const SMatrix6N& Is = _Is[idx];
 
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Is: " << Is << endl;

    // don't update I for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // compute a couple of necessary matrices
    transpose_solve_sIs(idx, s, _sIss);
    SpatialABInertiad uI = I - SpatialABInertia::mult(mult(Is, _sIss, tmp), I);

    // get the parent current zero acceleration and inertia
    SpatialABInertiad& Iim1 = _I[pidx];

    // handle cases of global and link frames separately
    if (rftype == eGlobal)
    {
      Iim1 += uI;
      FILE_LOG(LOG_DYNAMICS) << "    parent I is now: " << endl << Iim1;
    }
    else
    {
      // get the backward spatial transform and compute transformed I
      SpatialTransform X_im1_i = link->get_spatial_transform_backward();
      SpatialABInertiad tI = X_im1_i.transform(uI);
        
      // update I
      Iim1 += tI;

      FILE_LOG(LOG_DYNAMICS) << "    untransformed I: " << endl << uI;
      FILE_LOG(LOG_DYNAMICS) << "    backward transform: " << endl << X_im1_i;
      FILE_LOG(LOG_DYNAMICS) << "    transformed I: " << endl << tI;
      FILE_LOG(LOG_DYNAMICS) << "    parent I is now: " << endl << Iim1;
    }
  }
}

/// Computes joint and spatial link accelerations 
void FSABAlgorithm::calc_spatial_accelerations(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  queue<RigidBodyPtr> link_queue;
  VectorNd result;
  SVector6 tmp, tmp2;

  // get the links
  const vector<RigidBodyPtr>& links = body->get_links();

  // clear spatial values for all links
  _a.resize(links.size());

  // get the base link
  RigidBodyPtr base = links.front();

  // verify that base is index 0
  assert(base->get_index() == 0);

  // set spatial acceleration of base
  if (!body->is_floating_base())
  {
    _a.front().set_zero(); 
    FILE_LOG(LOG_DYNAMICS) << "  negated base Z: " << (-_Z.front()) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base acceleration: " << _a.front() << endl;
  }
  else
  {
    _a.front() = _I.front().inverse_mult(-_Z.front());
    base->set_spatial_accel(_a.front(), rftype);
  }
  
  // *****************************************************************
  // compute joint accelerations
  // *****************************************************************

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // process all links for forward recursion
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();
    unsigned idx = link->get_index();

    // push all children of the link onto the queue
    push_children(link, link_queue);
      
    // get the parent link and the inner joint
    RigidBodyPtr parent(link->get_parent_link());
    JointPtr joint(link->get_inner_joint_implicit());

    // get the parent index
    unsigned pidx = parent->get_index();
    
    // compute transformed parent link acceleration
    Twistd aim1;
    if (rftype == eGlobal) 
      aim1 = _a[pidx];
    else
      aim1 = link->get_spatial_transform_forward().transform(_a[pidx]);

    // get the spatial axis and its derivative
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);
    const vector<Twistd>& s_dot = joint->get_spatial_axes_dot(rftype);

    // get the Is and qm subexpressions
    const SMatrix6N& Is = _Is[idx];
    const VectorNd& mu = _mu[idx];    
    const SVector6& c = _c[idx];

    // compute joint i acceleration
    transpose_mult(Is, aim1, result);
    result.negate();
    result += mu;
    solve_sIs(idx, result, joint->qdd);
    
    // compute link i spatial acceleration
    _a[idx] = aim1 + c + mult(s_dot, joint->qd) + mult(s, joint->qdd);
    link->set_spatial_accel(_a[idx], rftype);

    FILE_LOG(LOG_DYNAMICS) << endl << endl << "  *** Forward recursion processing link " << link->id << endl;  
    FILE_LOG(LOG_DYNAMICS) << "    aim1: " << aim1 << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm(subexp): " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qdd: " << joint->qdd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    spatial acceleration: " << _a[idx] << endl;
  }
  
  FILE_LOG(LOG_DYNAMICS) << "    joint accel: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << JointPtr(links[i]->get_inner_joint_implicit())->qdd << " ";
  FILE_LOG(LOG_DYNAMICS) << "joint q: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << JointPtr(links[i]->get_inner_joint_implicit())->q << " ";
  FILE_LOG(LOG_DYNAMICS) << endl;
  FILE_LOG(LOG_DYNAMICS) << "calc_joint_accels() exited" << endl;

  FILE_LOG(LOG_DYNAMICS) << endl << "Outputting link data" << endl;
  for (unsigned i=0; i< links.size(); i++)
  {
    FILE_LOG(LOG_DYNAMICS) << "*** Link: " << links[i]->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  v: " << _v[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  a: " << _a[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  c: " << _c[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  Z: " << _Z[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  I: " << endl << _I[links[i]->get_index()];
  }
}

/// Computes the joint accelerations (forward dynamics) for a  manipulator
/**
 * Featherstone Algorithm taken from Mirtich's thesis (p. 113).  Note that Mirtich's numbering is a little funny;
 * I decrement his joint indices by one, while leaving his link indices intact.
 */
void FSABAlgorithm::calc_fwd_dyn()
{
  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() entered" << endl;

  // get the body and the reference frame
  RCArticulatedBodyPtr body(_body);
  if (!body->_ejoints.empty())
    throw std::runtime_error("FSABAlgorithm cannot process bodies with kinematic loops!");
  ReferenceFrameType rftype = body->computation_frame_type;

  // save the reference frame
  if (_rftype != rftype)
    _position_data_valid = false;
  _rftype = rftype;

  // get the links and joints for the body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_implicit_joints();

  // store reference frame
  _rftype = rftype;

  // get the base link
  RigidBodyPtr base = links.front();
 
  // invalidate spatial values for all joints, if necessary 
  if (body->positions_invalidated())
  {
    for (unsigned i=0; i< joints.size(); i++)
      joints[i]->reset_spatial_axis();
    body->validate_positions();
  }

  // compute spatial isolated inertias, if necessary
  if (!_position_data_valid)
    set_spatial_iso_inertias(body, rftype);

  // compute spatial coriolis vectors, if necessary
  if (!_velocity_data_valid)
    calc_spatial_coriolis_vectors(body, rftype);

  // compute spatial articulated body inertias, if necessary
  if (!_position_data_valid)
    calc_spatial_inertias(body, rftype);

  // compute spatial ZAs
  calc_spatial_zero_accelerations(body, rftype);
    
  // compute spatial accelerations
  calc_spatial_accelerations(body, rftype);

  // validate position and velocity data
  _position_data_valid = true;

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() exited" << endl;
}

/// Signum function with NEAR_ZERO cutoff
double FSABAlgorithm::sgn(double x)
{
  if (x > NEAR_ZERO)
    return (double) 1.0;
  else if (x < -NEAR_ZERO)
    return (double) -1.0;
  else
    return (double) 0.0;
}

/// Pushes all children of the given link onto the given queue
void FSABAlgorithm::push_children(RigidBodyPtr link, queue<RigidBodyPtr>& q)
{
  const list<RigidBody::OuterJointData>& ojd_list = link->get_outer_joints_data();
  BOOST_FOREACH(const RigidBody::OuterJointData& ojd, ojd_list)
  {
    RigidBodyPtr child(ojd.child);
    q.push(child);
  }
}

/// Computes necessary vectors to apply impulses to an articulated body
/**
 * Featherstone Algorithm taken from Mirtich's thesis (p. 113).  Note that Mirtich's numbering is a little funny;
 * I decrement his joint indices by one, while leaving his link indices intact.
 */
void FSABAlgorithm::calc_impulse_dyn(RCArticulatedBodyPtr body, ReferenceFrameType rftype)
{
  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_impulse_dyn() entered" << endl;

  // save the reference frame
  if (_rftype != rftype)
  {
    _position_data_valid = false;
    _velocity_data_valid = false;
  }
  _rftype = rftype;

  // get the links and joints for the body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_implicit_joints();

  // invalidate spatial values for all joints
  if (body->positions_invalidated())
  {
    for (unsigned i=0; i< joints.size(); i++)
      joints[i]->reset_spatial_axis();
    body->validate_positions();
  }

  // store reference frame
  _rftype = rftype;

  // get the base link
  RigidBodyPtr base = links.front();

  // initialize base velocities to zero, if not a floating base
  if (!body->is_floating_base())
  {
    base->set_avel(ZEROS_3);
    base->set_lvel(ZEROS_3);
  }
  
  if (!_position_data_valid)
  {
    // compute everything else we need
    set_spatial_iso_inertias(body, rftype);
    calc_spatial_inertias(body, rftype);
    _position_data_valid = true;
  }

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_impulse_dyn() exited" << endl;
}

/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
/**
 * \pre forward dynamics have been computed already by calc_fwd_dyn()
 */
void FSABAlgorithm::apply_impulse(const Vector3& j, const Vector3& k, const Vector3& point, RigidBodyPtr link)
{
  SAFESTATIC MatrixNd tmp, sT;
  SAFESTATIC VectorNd tmp2;

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << " -- applying linear impulse " << j << ", angular impulse " << k << " to " << point << endl;

  // get the computation reference frame
  RCArticulatedBodyPtr body(_body);   
  if (!body->_ejoints.empty())
    throw std::runtime_error("FSABAlgorithm cannot process bodies with kinematic loops!");
  ReferenceFrameType rftype = body->computation_frame_type;

  // prepare to compute impulses
  calc_impulse_dyn(body, rftype);

  // initialize spatial zero velocity deltas to zeros
  const unsigned NUM_LINKS = body->get_links().size();
  _Y.resize(NUM_LINKS);
  for (unsigned i=0; i< NUM_LINKS; i++)
    _Y[i].set_zero();

  // **********************************************************************
  // NOTE: uses articulated body inertias and spatial axes already computed 
  // **********************************************************************

  // determine the orientation of the target spatial transform
  Matrix3 R;
  if (rftype == eGlobal)
    R = Matrix3::identity();
  else
    link->get_pose().get_rotation(&R);

  // determine the origin of the target spatial transform
  const Vector3& o = (rftype == eGlobal) ? ZEROS_3 : link->get_position();

  // form the spatial impulse from frame with no rotation, origin at contact point
  // to frame with no rotation (global reference frame) / link rotation (link
  // reference frame), origin at link center
  SVector6 pcoll = SpatialTransform(IDENTITY_3x3, point, R, o).transform(SVector6(j, k));
 
  // get the index for this link
  unsigned i = link->get_index();
  
  // save the spatial impulse
  _Y[i] = -pcoll;
  
  FILE_LOG(LOG_DYNAMICS) << "  -- impulse applied to link " << link->id << " = " << pcoll << endl;
  FILE_LOG(LOG_DYNAMICS) << "  -- recursing backward" << endl;

  // do not recurse backward if this is the base
  if (!link->is_base())
  {
    // recurse backward
    while (!link->is_base())
    {
      // get the parent link
      RigidBodyPtr parent(link->get_parent_link());
    
      // get the parent index
      unsigned h = parent->get_index();

      // get spatial axes of the inner joint for link i
      JointPtr joint(link->get_inner_joint_implicit());
      const vector<Twistd>& s = joint->get_spatial_axes(rftype);
    
      // get Is for link i
      const SMatrix6N& Is = _Is[i];
    
      // compute Is * inv(sIs) * s'
      transpose_solve_sIs(i, s, _sIss);
      MatrixNd::mult(Is, _sIss, tmp); 
/*
      SMatrix6 IssIss(tmp.begin()); 

      // compute impulse for h in i's frame (or global frame)
      SVector6 _Yi = (SMatrix6::identity() - IssIss) * _Y[i];
*/
      Wrenchd Yi = _Y[i] - (SpatialABInertiad(tmp) * _Y[i]);

      // transform the spatial impulse, if necessary
      if (rftype == eGlobal)
        _Y[h] = Yi;
      else
        _Y[h] = link->get_spatial_transform_backward().transform(Yi);
   
      FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- this transformed impulse is: " << _Y[i] << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- parent is link: " << h << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- transformed spatial impulse for parent: " << _Y[h] << endl; 
    
      // update the link to be the parent
      link = parent;
      i = h;
    }
  }

  // ************************************************************
  // determine the new joint and link velocities
  // ************************************************************
  
  // setup a queue for recursing forward through the links
  queue<RigidBodyPtr> link_queue;

  // get the base link
  const vector<RigidBodyPtr>& links = body->get_links();
  RigidBodyPtr base = links.front();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    // update the base velocity
    _v.front() += _dv.front();

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base linear velocity: " << base->get_lvel() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  current base angular velocity: " << base->get_avel() << endl;

    // set the base linear and angular velocities
    base->set_spatial_velocity(_v.front(), rftype);

    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << _Y.front() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base linear velocity: " << base->get_lvel() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base angular velocity: " << base->get_avel() << endl;
  }
  else 
    _dv.front().set_zero();
  
  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // get the parent link
    RigidBodyPtr parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    JointPtr joint(link->get_inner_joint_implicit());
    
    // get spatial axes of the inner joint for link i
    const vector<Twistd>& s = joint->get_spatial_axes(rftype);
    
    // get articulated body inertia for link i
    const SpatialABInertiad& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- s: " << s << endl;
    
    // determine the joint and link velocity updates
    if (rftype == eGlobal)
    {
      s.transpose_mult((I * _dv[h]) + _Y[i], tmp2);
      solve_sIs(i, tmp2, _qd_delta).negate();
      _dv[i] = _dv[h] + s.mult(_qd_delta);

      FILE_LOG(LOG_DYNAMICS) << "    -- I * dv[parent]: " << _dv[h] << endl;
    }
    else
    {
      SpatialTransform X_i_im1 = link->get_spatial_transform_forward();
      SVector6 xdv = X_i_im1.transform(_dv[h]);
      s.transpose_mult((I * xdv) + _Y[i], tmp2).negate();
      solve_sIs(i, tmp2, _qd_delta);
      _dv[i] = xdv + s.mult(_qd_delta);

      FILE_LOG(LOG_DYNAMICS) << "    -- I * X_i_im1 * dv[parent]: " << (I * X_i_im1.transform(_dv[h])) << endl;
    }
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    _v[i] += _dv[i];
    link->set_spatial_velocity(_v[i], rftype);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << "  v: " << _v[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  const vector<JointPtr>& joints = body->get_implicit_joints();
  for (unsigned i=0; i< joints.size(); i++)
    joints[i]->reset_force(); 
  for (unsigned i=0; i< NUM_LINKS; i++)
    links[i]->reset_accumulators();

  // recompute the coriolis vectors so we can validate velocity data
  calc_spatial_coriolis_vectors(body, rftype);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_impulse() exited" << endl;
}

/// Solves a system for sIs*x = m' using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MatrixNd& FSABAlgorithm::transpose_solve_sIs(unsigned idx, const SMatrix6N& m, MatrixNd& result) const
{
  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[idx])
    return SMatrix6N::mult_transpose(_sIs[idx],m,result);
  else
  {
    _tmpM = _sIs[idx];
    LinAlg::inverse_chol(_tmpM);
    return SMatrix6N::mult_transpose(_tmpM,m,result);
  }
}

/// Solves a system for sIs*x = m using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MatrixNd& FSABAlgorithm::solve_sIs(unsigned idx, const MatrixNd& m, MatrixNd& result) const
{
  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[idx])
    return _sIs[idx].mult(m, result);
  else
  {
    result = m;
    LinAlg::solve_chol_fast(_sIs[idx], result);
    return result;
  }
}

/// Solves a system for sIs*x = v using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
VectorNd& FSABAlgorithm::solve_sIs(unsigned idx, const VectorNd& v, VectorNd& result) const
{
  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[idx])
    return _sIs[idx].mult(v, result);
  else
  {
    result = v;
    LinAlg::solve_chol_fast(_sIs[idx], result);
    return result;
  }
}

