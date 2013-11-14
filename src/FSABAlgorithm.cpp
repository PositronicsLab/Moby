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
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/NumericalException.h>
#include <Moby/Spatial.h>
#include <Moby/FSABAlgorithm.h>

using namespace Ravelin;
using namespace Moby;
using std::set;
using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;

FSABAlgorithm::FSABAlgorithm()
{
}

/*
/// Calculates the inverse generalized inertia matrix
void FSABAlgorithm::calc_inverse_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixNd& iM)
{
  vector<SVelocityd> sprime;

  // we currently only handle axis-angle generalized coordinates...
  assert(gctype == DynamicBody::eSpatial);

  // get the body
  RCArticulatedBodyPtr body(_body);

  // prepare to calculate the impulses
  calc_impulse_dyn(body, eGlobal);

  // now, compute s'I
  const vector<RigidBodyPtr>& links = body->get_links();
  vector<SForced> sTI(links.size());
  for (unsigned i=1; i< links.size(); i++)
  {
    JointPtr joint = links[i]->get_inner_joint_explicit();
    if (links[i]->get_computation_type() == eJoint)
    {
      const vector<SVelocityd>& s = joint->get_spatial_axes();
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
*/

/// Applies a generalized impulse using the algorithm of Drumwright
void FSABAlgorithm::apply_generalized_impulse(unsigned index, const vector<vector<SForced> >& sTI, VectorNd& vgj)
{
  queue<RigidBodyPtr> link_queue;
  const unsigned SPATIAL_DIM = 6;
  const unsigned BASE_IDX = 0;
  vector<SVelocityd> sprime;

  // get the body as well as sets of links and joints
  RCArticulatedBodyPtr body(_body);
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_explicit_joints();

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = (body->is_floating_base()) ? SPATIAL_DIM : 0;

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << vgj << endl;

  // clear values for vectors
  _Y.resize(links.size());
  _processed.resize(links.size());

  // reset Y and processed
  for (unsigned j=0; j< links.size(); j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
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
    unsigned i = link_pqueue.top();
    RigidBodyPtr link = links[i];
    link_pqueue.pop();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (_processed[i])
      continue;
  
    // indicate that this link has been processed
    _processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->get_index() == BASE_IDX)
    {
      link_pqueue.push(parent->get_index());
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_explicit());
    const vector<SVelocityd>& s = joint->get_spatial_axes();

    // get I, Y, and mu
    const SpatialABInertiad& I = _I[i];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && h == BASE_IDX)
      continue;
 
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    vgj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[i]);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu[i] << endl;
    
    // compute the qm subexpression
    Pose3d::transform(_Y[i].pose, s, sprime);
    _mu[i] -= transpose_mult(sprime, _Y[i], _workv2);

    // get Is and sIsmu
    const vector<SMomentumd>& Is = _Is[i];

    // update parent impulsive force 
    solve_sIs(i, _mu[i], _sIsmu);
    SMomentumd Y = _Y[i] + SMomentumd::from_vector(_sIsmu, _Y[i].pose);
    _Y[h] += Pose3d::transform(_Y[h].pose, Y);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    parent Y is now: " << endl << _Y[h] << endl;
  }

  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // update base components 
    RigidBodyPtr base = links.front();

    // generalized impulse on base will be in global frame, by Moby convention 
    SMomentumd w(vgj[0], vgj[1], vgj[2], vgj[3], vgj[4], vgj[5], GLOBAL);
    _Y.front() += Pose3d::transform(_Y.front().pose, w);
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
    JointPtr joint(link->get_inner_joint_explicit());

    // check the starting index of this joint
    const unsigned CSTART = joint->get_coord_index();
    if (CSTART > index)
      continue; 
    
    // get spatial axes of the inner joint for link i
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
    
    // determine appropriate components of gj
    vgj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[i]);

    // transform s to Y's frame
    Pose3d::transform(_Y[i].pose, s, sprime); 

    // update _Yi
    SMomentumd Y = _Y[i] + _I[i]*_dv[h]; 

    // compute -s'Y
    transpose_mult(sprime, Y, _workv2);
    _workv2.negate();

    // solve using s'Is
    solve_sIs(i, _workv2, _qd_delta);

    // update the joint velocity   
    _dv[i] = Pose3d::transform(_dv[i].pose, _dv[h]);
    _dv[i] += mult(sprime, _qd_delta);

/*
    // compute s'Y
    transpose_mult(sprime, _Y[i], _sTY);

    // determine the joint and link velocity updates
    transpose_mult(sTI[i], _dv[h], _workv2) += _sTY;
    _workv2.negate() += _mu[i];
    solve_sIs(i, _workv2, _qd_delta);
    _dv[i] = Pose3d::transform(_dv[h].pose, _dv[i].pose, _dv[h]);
    _dv[i] += mult(sprime, _qd_delta);
*/
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- I: " << endl << _I[i];
    FILE_LOG(LOG_DYNAMICS) << "    -- Qi: " << _mu[i] << endl;
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
void FSABAlgorithm::apply_generalized_impulse(const VectorNd& gj)
{
  SAFESTATIC VectorNd tmp, tmp2;
  SAFESTATIC vector<VectorNd> mu;
  queue<RigidBodyPtr> link_queue;
  vector<SVelocityd> sprime;

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = 6;

  // get the body
  RCArticulatedBodyPtr body(_body);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << gj << endl;

  // prepare to calculate the impulse
  calc_impulse_dyn(body);

  // get the sets of links and joints
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_explicit_joints();

  // init mu and Y 
  _mu.resize(links.size());
  _Y.resize(links.size());
  for (unsigned j=0; j< links.size(); j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
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
    RigidBodyPtr link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index(); 
 
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;
 
    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue;
 
    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_explicit());
    const vector<SVelocityd>& s = joint->get_spatial_axes();

    // get I
    const SpatialABInertiad& I = _I[i];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[i]);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu[i] << endl;
    
    // compute the qm subexpression
    Pose3d::transform(_Y[i].pose, s, sprime);
    _mu[i] -= transpose_mult(sprime, _Y[i], tmp2);

    // get Is
    const vector<SMomentumd>& Is = _Is[i];

    // prepare to update parent Y
    solve_sIs(i, _mu[i], _sIsmu);
    SMomentumd uY = _Y[i] + SMomentumd::from_vector(mult(Is, _sIsmu, tmp), _Y[i].pose);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[i] << endl;

    // update parent force
    _Y[h] += Pose3d::transform(_Y[h].pose, _Y[i]); 
  }

  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // update base components 
    RigidBodyPtr base = links.front();

    // momentum is in global frame by Moby convention 
    SMomentumd basew(gj[0], gj[1], gj[2], gj[3], gj[4], gj[5], GLOBAL);
    _Y.front() += Pose3d::transform(_Y.front().pose, basew);
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
    base->set_velocity(base->get_velocity() + _dv.front());

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base velocity: " << base->get_velocity() << endl;

    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << _Y.front() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << base->get_velocity() << endl;
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
    JointPtr joint(link->get_inner_joint_explicit());
    
    // get spatial axes of the inner joint for link i
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    
    // get articulated body inertia for link i
    const SpatialABInertiad& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
    
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _Qi);

    // determine the joint and link velocity updates
    Pose3d::transform(_Y[i].pose, s, sprime);
    SMomentumd w1 = _I[i] * Pose3d::transform(_Y[i].pose, _dv[h]);
    transpose_mult(sprime, w1 + _Y[i], tmp2).negate();
    tmp2 += _Qi;
    solve_sIs(i, tmp2, _qd_delta);
    _dv[i] = Pose3d::transform(_Y[i].pose, _dv[h]);
    if (joint->num_dof() > 0)
      _dv[i] += mult(sprime, _qd_delta);
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    link->set_velocity(link->get_velocity() + _dv[i]);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_generalized_impulse() exited" << endl;
}

/// Computes the combined spatial coriolis / centrifugal forces vectors for the body
void FSABAlgorithm::calc_spatial_coriolis_vectors(RCArticulatedBodyPtr body)
{
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_coriolis_vectors() entered" << endl;
  vector<SVelocityd> sprime;

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
    JointPtr joint(link->get_inner_joint_explicit());

    // get the spatial axis and transform it
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    Pose3d::transform(link->get_computation_frame(), s, sprime);

    // compute the coriolis vector
    SVelocityd sqd = mult(sprime, joint->qd);
    _c[idx] = link->get_velocity().cross(sqd);

    FILE_LOG(LOG_DYNAMICS) << "processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "v: " << link->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "s * qdot: " << mult(sprime, joint->qd) << endl;
    FILE_LOG(LOG_DYNAMICS) << "c: " << _c[idx] << endl;
  }
}

/// Computes articulated body zero acceleration forces used for computing forward dynamics
void FSABAlgorithm::calc_spatial_zero_accelerations(RCArticulatedBodyPtr body)
{
  VectorNd tmp, workv;
  MatrixNd workM;
  queue<RigidBodyPtr> link_queue;
  vector<SVelocityd> sprime;

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
  for (unsigned j=start_idx; j< links.size(); j++)
  {
    // get the link
    RigidBodyPtr link = links[j];
    const unsigned i = link->get_index();

    // get the spatial link velocity
    const SVelocityd& v = link->get_velocity();

    // set 6-dimensional spatial isolated zero-acceleration vector of link  
    _Z[i] = v.cross(link->get_inertia() * v) - link->sum_forces();
    
    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso ZA: " << endl << _Z[i] << endl;
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
    RigidBodyPtr link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index(); 
 
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;
  
    // verify that all children have already been processed
    if (!body->all_children_processed(link))
      continue;

    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
      link_pqueue.push(parent->get_index());
    const unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_explicit());
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    Pose3d::transform(link->get_computation_frame(), s, sprime);

    // get I, c, and Z
    const SpatialABInertiad& I = _I[i];
    const SAcceld& c = _c[i];
    const SForced& Z = _Z[i];
    
    // compute the qm subexpression
    joint->get_scaled_force(_Q);
    _Q += joint->ff;
    _mu[i] = _Q;
    _mu[i] -= transpose_mult(sprime, Z + I*c, workv);

    // get Is
    const vector<SMomentumd>& Is = _Is[i];

    // get the qm subexpression
    const VectorNd& mu = _mu[i];
  
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    parent link: " << parent->id << endl;
    if (LOGGING(LOG_DYNAMICS) && !sprime.empty()) 
    FILE_LOG(LOG_DYNAMICS) << "    s': " << sprime.front() << endl;  
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    c: " << c << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Z: " << Z << endl;

    // don't update Z for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // compute a couple of necessary matrices
    solve_sIs(i, mu, _sIsmu);
    SForced uZ = Z + (I*c) + SForced::from_vector(mult(Is, _sIsmu, workv), Z.pose);

    // update the parent zero acceleration
    _Z[h] += Pose3d::transform(_Z[h].pose, uZ);
  }

  FILE_LOG(LOG_DYNAMICS) << endl;
  unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
  for (unsigned i=st_idx; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << "Articulated zero-acceleration vector for " << links[i]->id << ": " << _Z[links[i]->get_index()] << endl;  
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() ended" << endl;
}

/// Computes articulated body inertias used for computing forward dynamics
void FSABAlgorithm::calc_spatial_inertias(RCArticulatedBodyPtr body)
{
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() entered" << endl;
  vector<SVelocityd> sprime;
  MatrixNd tmp, tmp2, tmp3;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // clear spatial values for all links
  _rank_deficient.resize(links.size());
  _I.resize(links.size());
  _Is.resize(links.size());
  _sIs.resize(links.size());
  _usIs.resize(links.size());
  _ssIs.resize(links.size());
  _vsIs.resize(links.size());

  // ************************************************
  // form the isolated spatial inertia vectors
  // ************************************************

  // process all links (including the base, if it's floating)
  unsigned start_idx = (body->is_floating_base()) ? 0 : 1;

  // process all links except the base
  for (unsigned j=start_idx; j< links.size(); j++)
  {
    // get the link
    RigidBodyPtr link = links[j];
    const unsigned i = link->get_index();

    // set the articulated body inertia for this link to be its isolated
    // spatial inertia (this will be updated in the phase below)
   _I[i] = link->get_inertia();

    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso inertia: " << endl << _I[i];
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
    RigidBodyPtr link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;

    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue; 
  
    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    RigidBodyPtr parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
FILE_LOG(LOG_DYNAMICS) << "added link " << parent->id << " to queue for processing" << endl;
    }
    const unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    JointPtr joint(link->get_inner_joint_explicit());
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    Pose3d::transform(_I[i].pose, s, sprime);

    // get I
    const SpatialABInertiad& I = _I[i];
    
    // compute Is
    mult(I, sprime, _Is[i]);

    // compute sIs
    transpose_mult(sprime, _Is[i], _sIs[i]);

    // get whether s is rank deficient
    _rank_deficient[i] = joint->is_singular_config();

    // if the joint is not rank deficient, compute a Cholesky factorization 
    // of sIs
    if (_sIs[i].rows() == 1)
      _sIs[i].data()[0] = 1.0/_sIs[i].data()[0];
    else
    { 
      if (!_rank_deficient[i])
        _LA->factor_chol(_sIs[i]);
      else
        _LA->svd(_sIs[i], _usIs[i], _ssIs[i], _vsIs[i]);
    }

    // get Is
    const vector<SMomentumd>& Is = _Is[i];
 
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;

    // don't update I for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
// NOTE: if we allow a spatial inertia to be formed from a general matrix, we
// need to check structure
    // compute a couple of necessary matrices
    transpose_solve_sIs(i, sprime, _sIss);
    mult(Is, _sIss, tmp);
    I.to_matrix(tmp2);
    MatrixNd::mult(tmp, tmp2, tmp3);
    SpatialABInertiad uI = I - SpatialABInertiad::from_matrix(tmp3, I.pose);

    // output the updates
    if (LOGGING(LOG_DYNAMICS) && _Is[i].size() > 0)
      FILE_LOG(LOG_DYNAMICS) << "  Is: " << _Is[i][0] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  s/(s'Is): " << _sIss << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  Is*s/(s'Is): " << std::endl << tmp;
    FILE_LOG(LOG_DYNAMICS) << "  Is*s/(s'Is)*I: " << std::endl << tmp3;
    FILE_LOG(LOG_DYNAMICS) << "  inertial update: " << uI << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  transformed I: " << Pose3d::transform(_I[h].pose, uI) << std::endl;

    // update the parent inertia
    _I[h] += Pose3d::transform(_I[h].pose, uI);
  }
}

/// Computes joint and spatial link accelerations 
void FSABAlgorithm::calc_spatial_accelerations(RCArticulatedBodyPtr body)
{
  queue<RigidBodyPtr> link_queue;
  VectorNd result;
  vector<SVelocityd> sprime, sdotprime;

  // get the links
  const vector<RigidBodyPtr>& links = body->get_links();

  // get the base link
  RigidBodyPtr base = links.front();

  // verify that base is index 0
  assert(base->get_index() == 0);

  // set spatial acceleration of base
  if (!body->is_floating_base())
  {
    base->set_accel(SAcceld::zero(base->get_computation_frame()));
    FILE_LOG(LOG_DYNAMICS) << "  base acceleration: (zero)" << endl;
  }
  else
  {
    SAcceld a0 = _I.front().inverse_mult(-_Z.front());
    base->set_accel(a0);
    FILE_LOG(LOG_DYNAMICS) << "  articulated base inertia: " << Pose3d::transform(base->get_mixed_pose(), _I.front()) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  negated base Z: " << Pose3d::transform(base->get_mixed_pose(), -_Z.front()) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base acceleration: " << Pose3d::transform(base->get_mixed_pose(), a0) << endl;
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
    const unsigned i = link->get_index();

    // push all children of the link onto the queue
    push_children(link, link_queue);
      
    // get the parent link and the inner joint
    RigidBodyPtr parent(link->get_parent_link());
    JointPtr joint(link->get_inner_joint_explicit());

    // get the parent index
    const unsigned h = parent->get_index();
    
    // compute transformed parent link acceleration
    SAcceld ah = Pose3d::transform(link->get_computation_frame(), parent->get_accel()); 

    // get the spatial axis and its derivative
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    const vector<SVelocityd>& sdot = joint->get_spatial_axes_dot();

    // transform spatial axes
    Pose3d::transform(link->get_computation_frame(), s, sprime);
    Pose3d::transform(link->get_computation_frame(), sdot, sdotprime);

    // get the Is and qm subexpressions
    const VectorNd& mu = _mu[i];    
    const SAcceld& c = _c[i];

    // compute joint i acceleration
    SForced w = _I[i] * ah;
    transpose_mult(sprime, w, result);
    result.negate();
    result += mu;
    solve_sIs(i, result, joint->qdd);
    
    // compute link i spatial acceleration
    SAcceld ai = ah + c;
    if (!sprime.empty())
      ai += SAcceld(mult(sprime, joint->qdd));
    if (!sdotprime.empty())
      ai += SAcceld(mult(sdotprime, joint->qd));
    link->set_accel(ai);

    FILE_LOG(LOG_DYNAMICS) << endl << endl << "  *** Forward recursion processing link " << link->id << endl;  
    FILE_LOG(LOG_DYNAMICS) << "    a[h]: " << ah << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm(subexp): " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qdd: " << joint->qdd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    spatial acceleration: " << link->get_accel() << endl;
  }
  
  FILE_LOG(LOG_DYNAMICS) << "    joint accel: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << JointPtr(links[i]->get_inner_joint_explicit())->qdd << " ";
  FILE_LOG(LOG_DYNAMICS) << "joint q: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << JointPtr(links[i]->get_inner_joint_explicit())->q << " ";
  FILE_LOG(LOG_DYNAMICS) << endl;
  FILE_LOG(LOG_DYNAMICS) << "calc_joint_accels() exited" << endl;

  FILE_LOG(LOG_DYNAMICS) << endl << "Outputting link data" << endl;
  for (unsigned i=0; i< links.size(); i++)
  {
    FILE_LOG(LOG_DYNAMICS) << "*** Link: " << links[i]->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  v: " << links[i]->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  a: " << links[i]->get_accel() << endl;
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
  if (!body->_ijoints.empty())
    throw std::runtime_error("FSABAlgorithm cannot process bodies with kinematic loops!");

  // get the links and joints for the body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_explicit_joints();

  // get the base link
  RigidBodyPtr base = links.front();
 
  // compute spatial coriolis vectors
  calc_spatial_coriolis_vectors(body);

  // compute spatial articulated body inertias, if necessary
  calc_spatial_inertias(body);

  // compute spatial ZAs
  calc_spatial_zero_accelerations(body);
    
  // compute spatial accelerations
  calc_spatial_accelerations(body);

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
  const set<JointPtr>& ojs = link->get_outer_joints();
  BOOST_FOREACH(JointPtr j, ojs)
  {
    RigidBodyPtr child(j->get_outboard_link());
    q.push(child);
  }
}

/// Computes necessary vectors to apply impulses to an articulated body
/**
 * Featherstone Algorithm taken from Mirtich's thesis (p. 113).  Note that Mirtich's numbering is a little funny;
 * I decrement his joint indices by one, while leaving his link indices intact.
 */
void FSABAlgorithm::calc_impulse_dyn(RCArticulatedBodyPtr body)
{
  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_impulse_dyn() entered" << endl;

  // get the links and joints for the body
  const vector<RigidBodyPtr>& links = body->get_links();
  const vector<JointPtr>& joints = body->get_explicit_joints();

  // get the base link
  RigidBodyPtr base = links.front();

  // initialize base velocities to zero, if not a floating base
  if (!body->is_floating_base())
    base->set_velocity(SVelocityd::zero(base->get_computation_frame()));
  
  // compute everything else we need
  calc_spatial_inertias(body);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_impulse_dyn() exited" << endl;
}

/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
/**
 * \pre forward dynamics have been computed already by calc_fwd_dyn()
 */
void FSABAlgorithm::apply_impulse(const SMomentumd& w, RigidBodyPtr link)
{
  SAFESTATIC MatrixNd tmp;
  SAFESTATIC VectorNd tmp2, workv;
  vector<SVelocityd> sprime;

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::apply_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << " -- applying impulse " << w << endl;

  // get the computation reference frame
  RCArticulatedBodyPtr body(_body);   
  if (!body->_ijoints.empty())
    throw std::runtime_error("FSABAlgorithm cannot process bodies with kinematic loops!");

  // prepare to compute impulses
  calc_impulse_dyn(body);

  // initialize spatial zero velocity deltas to zeros
  const vector<RigidBodyPtr>& links = body->get_links();
  const unsigned NUM_LINKS = links.size();
  _Y.resize(NUM_LINKS);
  for (unsigned j=0; j< NUM_LINKS; j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
  } 

  // **********************************************************************
  // NOTE: uses articulated body inertias and spatial axes already computed 
  // **********************************************************************


  // get the index for this link
  unsigned i = link->get_index();
  
  // transform the impulse
  _Y[i] = -Pose3d::transform(link->get_computation_frame(), w);
  
  FILE_LOG(LOG_DYNAMICS) << "  -- impulse applied to link " << link->id << " = " << w << endl;
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
      const unsigned h = parent->get_index();

      // get spatial axes of the inner joint for link i
      JointPtr joint(link->get_inner_joint_explicit());
      const vector<SVelocityd>& s = joint->get_spatial_axes();
      Pose3d::transform(link->get_computation_frame(), s, sprime);
    
      // get Is for link i
      const vector<SMomentumd>& Is = _Is[i];
    
      // compute Is * inv(sIs) * s'
      transpose_solve_sIs(i, sprime, _sIss);
      mult(Is, _sIss, tmp); 
      tmp.mult(_Y[i], workv);

      // compute impulse for h in i's frame (or global frame)
      SMomentumd Yi = _Y[i] - SMomentumd::from_vector(workv,  _Y[i].pose);

      // transform the spatial impulse, if necessary
      _Y[h] = Pose3d::transform(_Y[h].pose, Yi);
   
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

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base velocity: " << base->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << _Y.front() << endl;

    // update the base velocity
    base->set_velocity(base->get_velocity() + _dv.front());
    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << base->get_velocity() << endl;
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
    JointPtr joint(link->get_inner_joint_explicit());
    
    // get spatial axes of the inner joint for link i
    const vector<SVelocityd>& s = joint->get_spatial_axes();
    Pose3d::transform(link->get_computation_frame(), s, sprime);
    
    // get articulated body inertia for link i
    const SpatialABInertiad& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent->id << endl;
    
    // determine the joint and link velocity updates
    SVelocityd dvh = Pose3d::transform(_dv[i].pose, _dv[h]);
    SMomentumd f = (I * dvh) + _Y[i];
    transpose_mult(sprime, (I * dvh) + _Y[i], tmp2);
    solve_sIs(i, tmp2, _qd_delta).negate();
    _dv[i] = dvh + mult(s, _qd_delta);

    FILE_LOG(LOG_DYNAMICS) << "    -- I * dv[parent]: " << _dv[h] << endl;
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    link->set_velocity(link->get_velocity() + _dv[i]);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorithm::appl_impulse() exited" << endl;
}

/// Solves a system for sIs*x = m' using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MatrixNd& FSABAlgorithm::transpose_solve_sIs(unsigned i, const vector<SVelocityd>& m, MatrixNd& result) const
{
  // transpose m
  transpose_to_matrix(m, result);   

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);
 
  return result;
}

/// Solves a system for sIs*x = m using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MatrixNd& FSABAlgorithm::solve_sIs(unsigned i, const MatrixNd& m, MatrixNd& result) const
{
  result = m;

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);
    
  return result;
}

/// Solves a system for sIs*x = v using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
VectorNd& FSABAlgorithm::solve_sIs(unsigned i, const VectorNd& v, VectorNd& result) const
{
  result = v;

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);

  return result;
}

