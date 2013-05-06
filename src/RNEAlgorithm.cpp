/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <queue>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/Spatial.h>
#include <Moby/RNEAlgorithm.h>

using Ravelin::VectorNd;
using Ravelin::MatrixNd;
using Ravelin::Wrenchd;
using Ravelin::Twistd;
using Ravelin::SpatialRBInertiad;
using Ravelin::Pose3d;
using namespace Moby;
using std::vector;
using std::map;
using std::endl;
using std::queue;
using std::list;

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics
/**
 * Computed joint actuator forces are stored in inv_dyn_data.
 */
map<JointPtr, VectorNd> RNEAlgorithm::calc_inv_dyn(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data)
{
  if (!body->is_floating_base())
    return calc_inv_dyn_fixed_base(body, inv_dyn_data);
  else
    return calc_inv_dyn_floating_base(body, inv_dyn_data);
}

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a fixed base
/**
 * Computed joint actuator forces are stored in inv_dyn_data.
 */
map<JointPtr, VectorNd> RNEAlgorithm::calc_inv_dyn_fixed_base(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<RigidBodyPtr> link_queue;
  map<RigidBodyPtr, RCArticulatedBodyInvDynData>::const_iterator idd_iter;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_fixed_base() entered" << endl;

  // get the reference frame for computation
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // ** STEP 0: compute isolated inertias 

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // ** STEP 1: compute velocities and accelerations

  // get the base link
  RigidBodyPtr base = links.front();

  // setup the vector of link accelerations
  vector<Twistd> accels(links.size(), Twistd::zero());

  // all accelerations should be in the same frame as the individual links
  for (unsigned i=0; i< links.size(); i++)
    accels[i].pose = links[i]->velocity().pose; 
  
  // add all child links of the base to the processing queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);
  
  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // push all children of the link onto the queue
    child_links.clear();
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      link_queue.push(rb);

    // get the link's parent
    RigidBodyPtr parent(link->get_parent_link());
    unsigned h = parent->get_index();

    // get the joint for this link
    JointPtr joint(link->get_inner_joint_implicit());

    // get the spatial link velocity
    const Twistd& v = link->velocity(); 

    // get the reference to the spatial link acceleration
    Twistd& a = accels[i];
 
    // get the current joint velocity
    const VectorNd& qd = joint->qd;

    // **** compute acceleration

    // get the desired joint acceleration
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VectorNd& qdd_des = idd_iter->second.qdd;  

    // put s into the proper frame (that of v/a) (if necessary)
    transform(s.front().pose, v.pose, s, sprime);  

    // add this link's contribution
    a += spatial_cross(v, mult(sprime, qd));
    a += mult(s, qdd_des);

    // put s into the proper frame (that of v/a) (if necessary)
    if (!s_dot.empty())
    {
      transform(s.front().pose, v.pose, s_dot, sprime);  
      a += mult(sprime, qd);
    }

    // now add parent's contribution
    a += Pose3d::transform(accels[h].pose, a.pose, accels[h]);

    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->id << endl;
//    FILE_LOG(LOG_DYNAMICS) << "  spatial axes: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << v << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << a << endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<Wrenchd> forces(links.size(), Wrenchd::zero());

  // all wrenches should be in the same frame as the individual links
  for (unsigned i=0; i< links.size(); i++)
    forces[i].pose = links[i]->velocity().pose;

   // add all leaf links to the queue
  for (unsigned i=1; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // if this link has already been processed, do not process it again
    if (processed[i])
      continue;

    // if the link is the base, continue the loop
    if (link->is_base())
      continue;
    
    // link is not the base; add the parent to the queue for processing
    RigidBodyPtr parent(link->get_parent_link());
    link_queue.push(parent);
    unsigned h = parent->get_index();

    // get the forces for this link and this link's parent
    Wrenchd& wi = forces[i];
    Wrenchd& wim1 = forces[h];
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << wi << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * accels[idx]) << endl;

    // get the spatial velocity for this link
    const Twistd& v = link->velocity();

    // add I*a to the link force
    wi += Iiso[idx] * accels[idx];

    // update the determined force to this link w/Coriolis + centrifugal terms
    wi += spatial_cross(v, Iiso[idx] * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << wi << endl;    

    // determine external forces in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const Wrenchd& wx = idd_iter->second.wext;  

    // convert wrench to appropriate frame, if necessary
    if (wx.pose == wi.pose)
      wi -= wx; 
    else
      wi -= Pose3d::transform(wx.pose, wi.pose, wx);

    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << wi << endl;

    // indicate that this link has been processed
    processed[idx] = true;

    // update the parent force, if the parent is not the base
    if (parent->is_base())
      continue;
    else 
      if (wi.pose == wim1.pose)
        wim1 += wi;
      else
        wim1 += Pose3d::transform(wi.pose, wim1.pose, wi); 
  }
  
  // ** STEP 3: compute joint forces

  // setup a map from joints to actuator forces
  map<JointPtr, VectorNd> actuator_forces;

  // compute actuator forces
  for (unsigned i=1; i< links.size(); i++)
  {
    RigidBodyPtr link = links[i];
    JointPtr joint(link->get_inner_joint_implicit());
    const vector<Twistd>& s = s_prime[joint->get_index()];
    VectorNd& Q = actuator_forces[joint]; 
    Moby::transpose_mult(s, forces[link->get_index()], Q);
  
    FILE_LOG(LOG_DYNAMICS) << "joint " << joint->id << " inner joint force: " << actuator_forces[joint] << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_fixed_base() exited" << endl;

  return actuator_forces;
}

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a fixed base
/**
 * \pre Uses joint accelerations computed by forward dynamics, so the 
 *      appropriate forward dynamics function must be run first.
 */
/*
void RNEAlgorithm::calc_constraint_forces(RCArticulatedBodyPtr body)
{
  queue<RigidBodyPtr> link_queue;
  vector<SpatialRBInertiad> Iiso;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_constraint_forces() entered" << endl;

  // get the reference frame for computation
  ReferenceFrameType rftype = body->computation_frame_type;

  // ** STEP 0: compute isolated inertias 

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // get the isolated inertiae
  Iiso.resize(links.size());
  for (unsigned i=1; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Iiso[idx] = links[i]->get_spatial_iso_inertia(rftype); 
  }

   // ** STEP 1: compute velocities and accelerations

  // get the base link
  RigidBodyPtr base = links.front();

  // setup the vector of link accelerations
  vector<Twistd> accels(links.size(), Twistd::zero());
  
  // add all child links of the base to the processing queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);

  // ** STEP 1: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<Wrenchd> forces(links.size(), Wrench::zero());

  // add all leaf links to the queue
  for (unsigned i=1; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue and push all children of the link onto the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();    
    unsigned idx = link->get_index();

    // if this link has already been processed, do not process it again
    if (processed[idx])
      continue;

    // if the link is the base, continue the loop
    if (link->is_base())
      continue;
    
    // link is not the base; add the parent to the queue for processing
    RigidBodyPtr parent(link->get_parent_link());
    link_queue.push(parent);
    unsigned pidx = parent->get_index();

    // get the forces for this link and this link's parent
    Wrenchd& fi = forces[idx];
    Wrenchd& fim1 = forces[pidx];

    // get this link's acceleration
    Twistd& a = link->get_accel(rftype);
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << fi << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * a) << endl;

    // get the spatial velocity for this link
    const Twistd& v = link->velocity();

    // add I*a to the link force
    fi += Iiso[idx] * a;

    // update the determined force to this link w/Coriolis + centrifugal terms
    fi += SVector6d::spatial_cross(v, Iiso[idx] * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << fi << endl;    

    // determine external forces in link frame
    const Wrenchd& wext = link->sum_wrench(); 
    shared_ptr<const Pose3d>& T = link->get_pose();
    Wrenchd fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

    // subtract external forces in the appropriate frame
    if (rftype == eGlobal)
    {
      SpatialTransform X_0_i = link->get_spatial_transform_link_to_global();
      fi -= X_0_i.transform(fx);
    }
    else
      fi -= fx;

    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << fi << endl;

    // indicate that this link has been processed
    processed[idx] = true;

    // update the parent force, if the parent is not the base
    if (parent->is_base())
      continue;
    else
    { 
      if (rftype == eGlobal)
        fim1 += fi;
      else
        fim1 += link->get_spatial_transform_backward().transform(fi);
    }
  }
  
  // ** STEP 2: compute constraint forces

  // compute actuator forces
  for (unsigned i=1; i< links.size(); i++)
  {
    RigidBodyPtr link = links[i];
    JointPtr joint(link->get_inner_joint_implicit());
    joint->get_spatial_constraints(rftype, s);
    transpose_mult(s, forces[link->get_index()], joint->lambda);
  
    FILE_LOG(LOG_DYNAMICS) << "joint " << joint->id << " constraint force: " << joint->lambda << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_constraint_forces() exited" << endl;
}
*/

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a floating base
/**
 * \param inv_dyn_data a mapping from links to the external forces (and
 *        torques) applied to them and to the desired inner joint
 *        accelerations; note that all links in the robot should be included
 *        in this map (even the base link, although inner joint acceleration
 *        is not applicable in that case and will be ignored for it)
 * \return a mapping from joints to actuator forces
 */
map<JointPtr, VectorNd> RNEAlgorithm::calc_inv_dyn_floating_base(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<RigidBodyPtr> link_queue;
  map<RigidBodyPtr, RCArticulatedBodyInvDynData>::const_iterator idd_iter;
  vector<SpatialRBInertiad> Iiso, I;
  vector<Twistd> v, a;
  vector<Wrenchd> Z;
  vector<vector<Twistd> > s_prime, s_dot_prime;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() entered" << endl;

  // get the reference frame type
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // ** STEP 0: compute isolated inertias 

  // get the isolated inertiae
  Iiso.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Iiso[idx] = links[i]->get_inertia(); 
  }

  // resize s_prime / s_prime_dot
  s_prime.resize(body->get_joints().size()); 
  s_dot_prime.resize(body->get_joints().size()); 

  // ** STEP 1: compute velocities and relative accelerations

  // set all desired velocities and accelerations (latter relative to the base)
  // to zero initially
  v.resize(links.size());
  a.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    v[i] = a[i] = Twistd::zero();
    v[i].pose = a[i].pose = links[i]->velocity().pose;
  }  

  // get the base link
  RigidBodyPtr base = links.front();
  
  // set velocity for the base
  v.front() = base->velocity();

  // add all child links of the base to the processing queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);
    
  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();
    
    // add all child links of the link to the processing queue
    child_links.clear();
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      link_queue.push(rb);
    
    // get the parent link and inner joint
    RigidBodyPtr parent(link->get_parent_link());
    JointPtr joint(link->get_inner_joint_implicit());
    
    // get the index of this link and its parent
    unsigned i = link->get_index();
    unsigned im1 = parent->get_index();

    // put s and s_dot into the proper frame (that of v/a) (if necessary)
    if (rftype != eJoint)
    {
      const vector<Twistd>& s = joint->get_spatial_axes();
      const vector<Twistd>& s_dot = joint->get_spatial_axes_dot();
      transform(s.front().pose, v[i].pose, s, s_prime[joint->get_index()]);
      transform(s_dot.front().pose, v[i].pose, s_dot, s_dot_prime[joint->get_index()]);
    }
    else
    {
      const vector<Twistd>& s = joint->get_spatial_axes();
      const vector<Twistd>& s_dot = joint->get_spatial_axes_dot();
      s_prime[joint->get_index()] = s;
      s_dot_prime[joint->get_index()] = s_dot;
    }

    // get spatial axes and time derivatives for this link's inner joint
    const vector<Twistd>& s = s_prime[joint->get_index()];
    const vector<Twistd>& s_dot = s_dot_prime[joint->get_index()];

    // compute s * qdot
    Twistd sqd = mult(s, joint->qd);
    
    // get the desired acceleration for the current link
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VectorNd& qdd_des = idd_iter->second.qdd;

    // compute velocity and relative acceleration
    v[i] = v[im1] + sqd;
    a[i] = a[im1] + mult(s, qdd_des) + mult(s_dot, joint->qd) + spatial_cross(v[i], sqd);

//    FILE_LOG(LOG_DYNAMICS) << "  s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  velocity for link " << links[i]->id << ": " << v[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  s * qdd: " << mult(s, qdd_des) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  v x s * qd: " << Twistd(Moby::spatial_cross(v[i], sqd)) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  relative accel for link " << links[i]->id << ": " << a[i] << endl;
  }
  
  // ** STEP 2: compute composite inertias and Z.A. forces

  // resize vectors of spatial inertias and Z.A. vectors
  I.resize(links.size());
  Z.resize(links.size());

  // zero out I and Z and set poses
  for (unsigned i=0; i< links.size(); i++)
  {
    I[i].set_zero();
    Z[i].set_zero();
    I[i].pose = Z[i].pose = links[i]->velocity().pose;;
  }

  // set all spatial isolated inertias and Z.A. forces
  for (unsigned i=0; i< links.size(); i++)
  {
    // get the i'th link
    RigidBodyPtr link = links[i];
    unsigned idx = link->get_index();

    // add the spatial isolated inertia for this link to the composite inertia
    I[idx] += Iiso[idx];

    // setup forces due to (relative) acceleration on link
    Z[idx] = Iiso[idx] * a[idx];

    // add coriolis and centrifugal forces on link
    Z[idx] += spatial_cross(v[i], Iiso[idx] * v[idx]);

    // determine external forces on the link in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const Wrenchd& wx = idd_iter->second.wext;

    // transform external forces and subtract from Z.A. vector
    if (wx.pose == Z[idx].pose)
      Z[idx] -= wx;
    else
      Z[idx] -= Pose3d::transform(wx.pose, Z[idx].pose, wx);

    FILE_LOG(LOG_DYNAMICS) << " -- processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "   external wrench: " << wx << endl;
    FILE_LOG(LOG_DYNAMICS) << "   ZA vector: " << Z[idx] << endl;
    FILE_LOG(LOG_DYNAMICS) << "   isolated spatial-inertia: " << endl << Iiso[idx];
  }
  
  // *** compute composite inertias and zero acceleration vectors

  // setup vector that indicates when links have been processed
  vector<bool> processed(links.size(), false);

  // put all leaf links into the queue
  for (unsigned i=0; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
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
    if (processed[i])
      continue;
    
    // process the parent link, if possible
    RigidBodyPtr parent(link->get_parent_link());
    if (parent)
    {
      // put the parent on the queue
      link_queue.push(parent);
    
      // get the parent index
      unsigned im1 = parent->get_index();
    
      // add this inertia and Z.A. force to its parent
      if (I[i].pose == I[im1].pose)
      {
        I[im1] += I[i];
        Z[im1] += Z[i];
      }
      else
      {
        I[im1] += Pose3d::transform(I[i].pose, I[im1].pose, I[i]);
        Z[im1] += Pose3d::transform(Z[i].pose, Z[im1].pose, Z[i]);
      }

      // indicate that the link has been processed
      processed[i] = true;
    }
  }

  // ** STEP 3: compute base acceleration
  a.front() = I.front().inverse_mult(-Z.front());
  
  FILE_LOG(LOG_DYNAMICS) << "  Composite inertia for the base: " << endl << I.front();
  FILE_LOG(LOG_DYNAMICS) << "  ZA vector for the base: " << Z.front() << endl;
  FILE_LOG(LOG_DYNAMICS) << "  Determined base acceleration: " << a.front() << endl;

  // ** STEP 4: compute joint forces
  
  // setup the map of actuator forces
  map<JointPtr, VectorNd> actuator_forces;

  // compute the forces
  for (unsigned i=1; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    JointPtr joint(links[i]->get_inner_joint_implicit());
    const vector<Twistd>& s = s_prime[joint->get_index()];
    VectorNd& Q = actuator_forces[joint];
    Moby::transpose_mult(s, (I[idx] * a.front() + Z[idx]), Q);

    FILE_LOG(LOG_DYNAMICS) << "  processing link: " << links[i]->id << endl;
//    FILE_LOG(LOG_DYNAMICS) << "    spatial axis: " << endl << s;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << endl << I[idx];
    FILE_LOG(LOG_DYNAMICS) << "    Z: " << endl << Z[idx];
    FILE_LOG(LOG_DYNAMICS) << "    actuator force: " << actuator_forces[joint] << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() exited" << endl;

  return actuator_forces;
}

