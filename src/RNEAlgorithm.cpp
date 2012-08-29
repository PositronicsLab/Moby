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
#include <Moby/RNEAlgorithm.h>

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
map<JointPtr, VectorN> RNEAlgorithm::calc_inv_dyn(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data)
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
map<JointPtr, VectorN> RNEAlgorithm::calc_inv_dyn_fixed_base(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<RigidBodyPtr> link_queue;
  map<RigidBodyPtr, RCArticulatedBodyInvDynData>::const_iterator idd_iter;
  vector<SpatialRBInertia> Iiso;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_fixed_base() entered" << endl;

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
  vector<SVector6> accels(links.size(), ZEROS_6);
  
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
    unsigned idx = link->get_index();

    // push all children of the link onto the queue
    child_links.clear();
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      link_queue.push(rb);

    // get the link's parent
    RigidBodyPtr parent(link->get_parent_link());
    unsigned pidx = parent->get_index();

    // get the joint for this link
    JointPtr joint(link->get_inner_joint_implicit());

    // get the spatial link velocity
    const SVector6& v = link->get_spatial_velocity(rftype); 

    // get the reference to the spatial link acceleration
    SVector6& a = accels[idx];
 
    // get spatial axes for this link's inner joint
    const SMatrix6N& s = joint->get_spatial_axes(rftype);

    // get derivative of the spatial axes for this link's inner joint
    const SMatrix6N& s_dot = joint->get_spatial_axes_dot(rftype);

    // get the current joint velocity
    const VectorN& qd = joint->qd;

    // **** compute acceleration

    // get the desired joint acceleration
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VectorN& qdd_des = idd_iter->second.qdd;  

    // add this link's contribution
    a += SVector6::spatial_cross(v, s.mult(qd)) + s.mult(qdd_des) + s_dot.mult(qd);

    // now add parent's contribution
    if (rftype == eGlobal)
      a += accels[pidx];
    else
    {
      SpatialTransform X_i_im1 = link->get_spatial_transform_forward();
      a += X_i_im1.transform(accels[pidx]);
    }

    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  spatial axis: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  spatial joint velocity: " << s.mult(qd) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << v << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << a << endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<SVector6> forces(links.size(), SVector6(0,0,0,0,0,0));

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
    SVector6& fi = forces[idx];
    SVector6& fim1 = forces[pidx];
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << fi << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * accels[idx]) << endl;

    // get the spatial velocity for this link
    const SVector6& v = link->get_spatial_velocity(rftype);

    // add I*a to the link force
    fi += Iiso[idx] * accels[idx];

    // update the determined force to this link w/Coriolis + centrifugal terms
    fi += SVector6::spatial_cross(v, Iiso[idx] * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << fi << endl;    

    // determine external forces in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const Vector3& fext = idd_iter->second.fext;  
    const Vector3& text = idd_iter->second.text;  
    const Matrix4& T = link->get_transform();
    SVector6 fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

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
      if (rftype == eGlobal)
        fim1 += fi;
      else
        fim1 += link->get_spatial_transform_backward().transform(fi);
  }
  
  // ** STEP 3: compute joint forces

  // setup a map from joints to actuator forces
  map<JointPtr, VectorN> actuator_forces;

  // compute actuator forces
  for (unsigned i=1; i< links.size(); i++)
  {
    RigidBodyPtr link = links[i];
    JointPtr joint(link->get_inner_joint_implicit());
    const SMatrix6N& s = joint->get_spatial_axes(rftype);
    VectorN& Q = actuator_forces[joint];  
    s.transpose_mult(forces[link->get_index()], Q);
  
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
void RNEAlgorithm::calc_constraint_forces(RCArticulatedBodyPtr body)
{
  queue<RigidBodyPtr> link_queue;
  SMatrix6N s;
  vector<SpatialRBInertia> Iiso;

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
  vector<SVector6> accels(links.size(), ZEROS_6);
  
  // add all child links of the base to the processing queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);

  // ** STEP 1: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<SVector6> forces(links.size(), ZEROS_6);

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
    SVector6& fi = forces[idx];
    SVector6& fim1 = forces[pidx];

    // get this link's acceleration
    SVector6 a = link->get_spatial_accel(rftype);
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << fi << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * a) << endl;

    // get the spatial velocity for this link
    const SVector6& v = link->get_spatial_velocity(rftype);

    // add I*a to the link force
    fi += Iiso[idx] * a;

    // update the determined force to this link w/Coriolis + centrifugal terms
    fi += SVector6::spatial_cross(v, Iiso[idx] * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << fi << endl;    

    // determine external forces in link frame
    const Vector3& fext = link->sum_forces();  
    const Vector3& text = link->sum_torques();  
    const Matrix4& T = link->get_transform();
    SVector6 fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

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
    s.transpose_mult(forces[link->get_index()], joint->lambda);
  
    FILE_LOG(LOG_DYNAMICS) << "joint " << joint->id << " constraint force: " << joint->lambda << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_constraint_forces() exited" << endl;
}

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a floating base
/**
 * \param inv_dyn_data a mapping from links to the external forces (and
 *        torques) applied to them and to the desired inner joint
 *        accelerations; note that all links in the robot should be included
 *        in this map (even the base link, although inner joint acceleration
 *        is not applicable in that case and will be ignored for it)
 * \return a mapping from joints to actuator forces
 */
map<JointPtr, VectorN> RNEAlgorithm::calc_inv_dyn_floating_base(RCArticulatedBodyPtr body, const map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<RigidBodyPtr> link_queue;
  map<RigidBodyPtr, RCArticulatedBodyInvDynData>::const_iterator idd_iter;
  vector<SpatialRBInertia> Iiso, I;
  vector<SVector6> Z, v, a;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() entered" << endl;

  // get the reference frame type
  ReferenceFrameType rftype = body->computation_frame_type;

  // get the set of links
  const vector<RigidBodyPtr>& links = body->get_links();

  // ** STEP 0: compute isolated inertias 

  // get the isolated inertiae
  Iiso.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Iiso[idx] = links[i]->get_spatial_iso_inertia(rftype); 
  }

  // ** STEP 1: compute velocities and relative accelerations

  // set all desired velocities and accelerations (latter relative to the base)
  // to zero initially
  v.resize(links.size());
  a.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
    v[i] = a[i] = ZEROS_6;
  
  // get the base link
  RigidBodyPtr base = links.front();
  
  // set velocity for the base
  v.front() = base->get_spatial_velocity(rftype);

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
    
    // get the parent link
    RigidBodyPtr parent(link->get_parent_link());
    
    // get the index of this link and its parent
    unsigned i = link->get_index();
    unsigned im1 = parent->get_index();

    // get the spatial axes (and derivative) of this link's inner joint
    JointPtr joint(link->get_inner_joint_implicit());
    const SMatrix6N& s = joint->get_spatial_axes(rftype);
    const SMatrix6N& s_dot = joint->get_spatial_axes_dot(rftype);

    // compute s * qdot
    SVector6 sqd = s.mult(joint->qd);
    
    // get the desired acceleration for the current link
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VectorN& qdd_des = idd_iter->second.qdd;

    // compute velocity and relative acceleration
    v[i] = v[im1] + sqd;
    a[i] = a[im1] + s.mult(qdd_des) + s_dot.mult(joint->qd) + SVector6::spatial_cross(v[i], sqd);

    FILE_LOG(LOG_DYNAMICS) << "  s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  velocity for link " << links[i]->id << ": " << v[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  s * qdd: " << s.mult(qdd_des) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  v x s * qd: " << SVector6::spatial_cross(v[i], sqd) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  relative accel for link " << links[i]->id << ": " << a[i] << endl;
  }
  
  // ** STEP 2: compute composite inertias and Z.A. forces

  // resize vectors of spatial inertias and Z.A. vectors
  I.resize(links.size());
  Z.resize(links.size());

  // zero out I and Z
  for (unsigned i=0; i< links.size(); i++)
  {
    I[i].set_zero();
    Z[i] = ZEROS_6;
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
    Z[idx] += SVector6::spatial_cross(v[i], Iiso[idx] * v[idx]);

    // determine external forces on the link in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const Vector3& fext = idd_iter->second.fext;  
    const Vector3& text = idd_iter->second.text;  
    const Matrix4& T = link->get_transform();
    SVector6 fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

    // transform external forces and subtract from Z.A. vector
    SpatialTransform X_0_i = link->get_spatial_transform_link_to_global();
    Z[idx] -= X_0_i.transform(fx);

    FILE_LOG(LOG_DYNAMICS) << " -- processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "   external force / torque: " << fext << " / " << text << endl;
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
      I[im1] += I[i];
      Z[im1] += Z[i];

      // indicate that the link has been processed
      processed[i] = true;
    }
  }

  // ** STEP 3: compute base acceleration
  a.front() = I.front().inverse_mult(-Z.front());
  
  SpatialTransform X_i_0 = base->get_spatial_transform_global_to_link(); 
  FILE_LOG(LOG_DYNAMICS) << "  Composite inertia for the base: " << endl << I.front();
  FILE_LOG(LOG_DYNAMICS) << "  ZA vector for the base (link frame): " << X_i_0.transform(Z.front()) << endl;
  FILE_LOG(LOG_DYNAMICS) << "  Determined base acceleration (link frame): " << X_i_0.transform(a.front()) << endl;

  // ** STEP 4: compute joint forces
  
  // setup the map of actuator forces
  map<JointPtr, VectorN> actuator_forces;

  // compute the forces
  for (unsigned i=1; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    JointPtr joint(links[i]->get_inner_joint_implicit());
    const SMatrix6N& s = joint->get_spatial_axes(rftype);
    VectorN& Q = actuator_forces[joint];
    s.transpose_mult((I[idx] * a.front()) + Z[idx], Q);

    FILE_LOG(LOG_DYNAMICS) << "  processing link: " << links[i]->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "    spatial axis: " << endl << s;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << endl << I[idx];
    FILE_LOG(LOG_DYNAMICS) << "    Z: " << endl << Z[idx];
    FILE_LOG(LOG_DYNAMICS) << "    actuator force: " << actuator_forces[joint] << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() exited" << endl;

  return actuator_forces;
}

