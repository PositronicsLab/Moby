/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <stack>
#include <queue>
#include <Moby/LinAlg.h>
#include <Moby/Log.h>
#include <Moby/AAngle.h>
#include <Moby/Joint.h>
#include <Moby/RigidBody.h>
#include <Moby/CRBAlgorithm.h>
#include <Moby/XMLTree.h>
#include <Moby/FSABAlgorithm.h>
#include <Moby/NumericalException.h>
#include <Moby/cblas.h>
#include <Moby/RCArticulatedBody.h>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::queue;
using std::list;
using std::map;
using std::string;

/// Default constructor
/**
 * Constructs a reduced-coordinate articulated body with no joints and no links.
 */
RCArticulatedBody::RCArticulatedBody()
{
  _floating_base = false;
  _n_joint_DOF_implicit = 0;
  
  // set default algorithm to FSAB and computation frame to global
  algorithm_type = eFeatherstone;
  computation_frame_type = eGlobal;

  // setup baumgarte parameters
  b_alpha = (Real) 0.0;
  b_beta = (Real) 0.0;

  // mark inverse inertia as invalid
  _invM_valid = false;
}

/// Determines whether all of the children of a link have been processed
bool RCArticulatedBody::all_children_processed(RigidBodyPtr link) const
{
  for (unsigned k=0; k< link->num_child_links(); k++)
    if (!_processed[link->get_child_link(k)->get_index()])
      return false;

  return true;
}

/// Gets the number of generalized coordinates for this body
unsigned RCArticulatedBody::num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const
{
  // look for trivial cases
  if (_links.empty() || _joints.empty() || !is_enabled())
    return 0;

  if (!_floating_base)
    return _n_joint_DOF_implicit;
  else
    return _n_joint_DOF_implicit + _links.front()->num_generalized_coordinates(gctype);
}

/// Multiplies the generalized inertia matrix by a matrix
MatrixN& RCArticulatedBody::generalized_inertia_mult(const MatrixN& M, MatrixN& result)
{
  if (_floating_base)
    generalized_inertia_mult_floating(M, result);
  else
    generalized_inertia_mult_fixed(M, result);

  return result;
}

/// Multiplies the generalized inertia matrix by a matrix
void RCArticulatedBody::generalized_inertia_mult_fixed(const MatrixN& M, MatrixN& result)
{
  // get the number of links
  const unsigned NLINKS = _links.size();

  // determine the parent for each link and the coordinate indices for the
  // inner joint
  _pidx.resize(NLINKS);
  _coord_indices.resize(NLINKS);
  for (unsigned i=0; i< NLINKS; i++)
  {
    unsigned li = _links[i]->get_index();
    RigidBodyPtr parent = _links[i]->get_parent_link();
    if (parent)
    {
      unsigned lh = parent->get_index();
      _pidx[li] = lh; 
      JointPtr joint = _links[i]->get_inner_joint_implicit();
      _coord_indices[li] = joint->get_coord_index();
    } 
    else
      _pidx[li] = std::numeric_limits<unsigned>::max();
  }

  // determine m and n
  const unsigned n = M.rows();
  const unsigned m = M.columns();

  // resize the result matrix
  result.resize(n,m);

  // compute the spatial axes
  _sx.resize(NLINKS);
  for (unsigned i=1; i< _links.size(); i++)
  {
    JointPtr joint = _links[i]->get_inner_joint_implicit();
    unsigned idx = _links[i]->get_index();
    _sx[idx] = &joint->get_spatial_axes(eGlobal);
  }

  // first, compute composite inertia matrices
  _Ic.resize(NLINKS);
  for (unsigned i=1; i< NLINKS; i++)
    _Ic[_links[i]->get_index()] = _links[i]->get_spatial_iso_inertia(eGlobal);
  for (unsigned i=NLINKS, idx=NLINKS-1; i> 0; i--,idx--)
  {
    unsigned li = _links[idx]->get_index();
    unsigned lh = _pidx[li];
    if (lh == std::numeric_limits<unsigned>::max() || lh == 0)
      continue;
    _Ic[lh] += _Ic[li];
  }

  // compute _dv and f for each column of M
  _dv.resize(NLINKS);
  _f.resize(NLINKS);

  // setup velocity change on the base
  _dv[0] = ZEROS_6;

  // iterate
  for (unsigned j=0; j< m; j++)
  {
    // get access to the desired column of M
    const Real* Mj = &M(0,j);

    // loop over each row
    for (unsigned i=1; i< n; i++)
    {
      // compute _dv
      unsigned li = _links[i]->get_index();
      unsigned lh = _pidx[li];
      if (lh != std::numeric_limits<unsigned>::max())
      {
        const unsigned NDOF = _sx[li]->columns();
        _workv.resize(NDOF);
        std::copy(Mj+_coord_indices[li], Mj+_coord_indices[li]+NDOF, _workv.begin());
        _dv[li] = _sx[li]->mult(_workv);
        _dv[li] += _dv[lh];
      }

      // compute f
      _f[li] = _Ic[li] * _dv[li];
    }

    // get the start of the column of the result
    Real* Uj = &result(0,j);

    // now compute the result column while updating f
    for (unsigned i=NLINKS-1; i>0; i--)
    {
      // do the multiplication and set appropriate part of u
      unsigned li = _links[i]->get_index();
      unsigned lh = _pidx[li];
      _sx[li]->transpose_mult(_f[li], _workv);
      std::copy(_workv.begin(), _workv.end(), Uj+_coord_indices[li]);
      if (lh != std::numeric_limits<unsigned>::max() && lh > 0)
        _f[lh] += _f[li];
    }
  }
}

/// Multiplies the generalized inertia matrix by a matrix
void RCArticulatedBody::generalized_inertia_mult_floating(const MatrixN& M, MatrixN& result)
{
  // get the number of links
  const unsigned NLINKS = _links.size();

  // determine the parent for each link and the coordinate indices for the
  // inner joint
  _pidx.resize(NLINKS);
  _coord_indices.resize(NLINKS);
  for (unsigned i=0; i< NLINKS; i++)
  {
    unsigned li = _links[i]->get_index();
    RigidBodyPtr parent = _links[i]->get_parent_link();
    if (parent)
    {
      unsigned lh = parent->get_index();
      _pidx[li] = lh; 
      JointPtr joint = _links[i]->get_inner_joint_implicit();
      _coord_indices[li] = joint->get_coord_index();
    } 
    else
      _pidx[li] = std::numeric_limits<unsigned>::max();
  }

  // determine m and n
  const unsigned n = M.rows();
  const unsigned m = M.columns();

  // resize the result matrix
  result.resize(n,m);

  // compute the spatial axes
  _sx.resize(NLINKS);
  for (unsigned i=1; i< _links.size(); i++)
  {
    JointPtr joint = _links[i]->get_inner_joint_implicit();
    unsigned idx = _links[i]->get_index();
    _sx[idx] = &joint->get_spatial_axes(eGlobal);
  }

  // first, compute composite inertia matrices
  _Ic.resize(NLINKS);
  for (unsigned i=0; i< NLINKS; i++)
    _Ic[_links[i]->get_index()] = _links[i]->get_spatial_iso_inertia(eGlobal);
  for (unsigned i=NLINKS, idx=NLINKS-1; i> 0; i--,idx--)
  {
    unsigned li = _links[idx]->get_index();
    unsigned lh = _pidx[li];
    if (lh == std::numeric_limits<unsigned>::max())
      continue;
    _Ic[lh] += _Ic[li];
  }

  // compute _dv and f for each column of M
  _dv.resize(NLINKS);
  _f.resize(NLINKS);

  // setup velocity change on the base
  _dv[0] = ZEROS_6;
  _f[0] = ZEROS_6;

  // iterate
  for (unsigned j=0; j< m; j++)
  {
    // get access to the desired column of M
    const Real* Mj = &M(0,j);

    // loop over each row
    for (unsigned i=1; i< n; i++)
    {
      // compute _dv
      unsigned li = _links[i]->get_index();
      unsigned lh = _pidx[li];
      if (lh != std::numeric_limits<unsigned>::max())
      {
        const unsigned NDOF = _sx[li]->columns();
        _workv.resize(NDOF);
        std::copy(Mj+_coord_indices[li], Mj+_coord_indices[li]+NDOF, _workv.begin());
        _dv[li] = _sx[li]->mult(_workv);
        _dv[li] += _dv[lh];
      }

      // compute f
      _f[li] = _Ic[li] * _dv[li];
    }

    // get the start of the column of the result
    Real* Uj = &result(0,j);

    // now compute the result column while updating f
    for (unsigned i=NLINKS, idx=NLINKS-1; i>0; i--,idx--)
    {
      // do the multiplication and set appropriate part of u
      unsigned li = _links[idx]->get_index();
      unsigned lh = _pidx[li];
      if (lh != std::numeric_limits<unsigned>::max())
        _f[lh] += _f[li];
    }

    // compute floating base velocity change 
    _dv[0] = -_Ic[0].inverse_mult(_f[0]);
    std::copy(_dv[0].begin(), _dv[0].end(), Uj);
 
    // compute joint accelerations
    for (unsigned i=1; i< NLINKS; i++)
    { 
      unsigned li = _links[i]->get_index();
      _sx[li]->transpose_mult(_Ic[li]*_dv[0] + _f[li], _workv);
      std::copy(_workv.begin(), _workv.end(), Uj+_coord_indices[li]);
    }
  }
}

/// Updates inverse generalized inertia matrix, as necessary
void RCArticulatedBody::update_factorized_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype)
{
  // see whether we need to update
  if (!is_enabled() || (_invM_valid && _invM_type == gctype))
    return;

  // we do need to  update; mark inverse(M) as non rank-deficient for now
  _invM_rankdef = false;

  // do the update
  if (algorithm_type == eFeatherstone)
    _fsab.calc_inverse_generalized_inertia(gctype, _invM);

  // indicate inverse inertia is valid    
  _invM_valid = true;
  _invM_type = gctype;
}

/// Gets the inverse of the generalized inertia matrix
MatrixN& RCArticulatedBody::get_generalized_inertia_inverse(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eRodrigues || _invM_rankdef)
    M.copy_from(_invM);
  else
  {
    assert(algorithm_type == eCRB && gctype == eAxisAngle);
    _workM.set_identity(num_generalized_coordinates(eAxisAngle));
    _crb.M_solve(_workM, M);
  }
 
  return M;
}

/// Solves using a generalized inertia matrix
VectorN& RCArticulatedBody::solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& v, VectorN& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eRodrigues || _invM_rankdef)
    _invM.mult(v, result);
  else
  {
    assert(algorithm_type == eCRB && gctype == eAxisAngle);
    _crb.M_solve(v, result);
  }
  
  return result;
}

/// Solves using a generalized inertia matrix
MatrixN& RCArticulatedBody::solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const MatrixN& m, MatrixN& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eRodrigues || _invM_rankdef)
    _invM.mult(m, result);
  else
  {
    assert(algorithm_type == eCRB && gctype == eAxisAngle);
    _crb.M_solve(m, result);
  }
  
  return result;
}

/// Solves using the transpose of the generalized inertia matrix
MatrixN& RCArticulatedBody::solve_generalized_inertia_transpose(DynamicBody::GeneralizedCoordinateType gctype, const MatrixN& m, MatrixN& result)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eRodrigues || _invM_rankdef)
    _invM.mult_transpose(m, result);
  else
  {
    assert(algorithm_type == eCRB && gctype == eAxisAngle);
    MatrixN::transpose(m, _workM);
    _crb.M_solve(_workM, result);
  }
  
  return result;
}

/// Applies a generalized impulse to the articulated body
void RCArticulatedBody::apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gj)
{
  if (!is_enabled())
    return;

  if (algorithm_type == eFeatherstone)
    _fsab.apply_generalized_impulse(gctype, gj);
  else
  {
    assert(algorithm_type == eCRB);

    // get the current generalized velocity
    get_generalized_velocity(gctype, _workv);

    // we'll solve for the change in generalized velocity
    solve_generalized_inertia(gctype, gj, _workv2);

    // apply the change in generalized velocity
    _workv += _workv2;
    set_generalized_velocity(gctype, _workv);
  }

  // reset the force and torque accumulators
  reset_accumulators();
}

/// Adds a generalized force to the articulated body
void RCArticulatedBody::add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf)
{
  if (!is_enabled())
    return;

  unsigned index = 0;

  if (_floating_base)
  {
    // get the base
    RigidBodyPtr base = _links.front();

    // first, get the force on the base link
    index = _links.front()->num_generalized_coordinates(gctype);
    gf.get_sub_vec(0, index, _workv);

    // handle Rodrigues parameter case 
    if (gctype == DynamicBody::eRodrigues)
    {
      // transform the linear force from the base frame
      Vector3 basef(_workv[0], _workv[1], _workv[2]);
      basef = base->get_transform().mult_vector(basef);
      _workv.set_sub_vec(0, basef);

      // transform the torque from the base frame, and then to a quaternion
      Vector3 baset(_workv[3], _workv[4], _workv[5]);
      baset = base->get_transform().mult_vector(baset);
      Quat q = base->get_orientation().G_transpose_mult(baset);
      q *= (Real) 2.0;
      _workv[3] = q.w;
      _workv[4] = q.x;
      _workv[5] = q.y;
      _workv[6] = q.z;
    }

    // add the generalized force to the base
    base->add_generalized_force(gctype, _workv);
  }

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // add to joint forces
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    gf.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _workv);
    _ijoints[i]->force += _workv;
  } 
}

/// Determines whether the link is effectively a leaf link
bool RCArticulatedBody::treat_link_as_leaf(RigidBodyPtr link) const
{
  // if all children have lower link indices, link treatable as end-effector
  BOOST_FOREACH(const RigidBody::OuterJointData& ojd, link->get_outer_joints_data())
  {
    RigidBodyPtr child(ojd.child);
    if (child->get_index() > link->get_index())
      return false;
  }

  return true;
}

/// Sets whether the base of this body is "floating" (or fixed)
void RCArticulatedBody::set_floating_base(bool flag)
{
  _floating_base = flag;
  compile();
}

/// Compiles this body (updates the link transforms and velocities)
void RCArticulatedBody::compile()
{
  // verify all links are enabled
  if (!is_floating_base())
  {
    for (unsigned i=1; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("Only first link can be disabled in a reduced coordinate body with fixed-base!");
  }
  else
  {
    for (unsigned i=0; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("No links can be disabled in a reduced coordinate body with floating-base!");
  }

  // update processed vector size
  _processed.resize(_links.size());

  // setup implicit joint generalized coordinate and constraint indices
  const unsigned SPATIAL_DIM = 6;
  const unsigned START_GC = (is_floating_base()) ? SPATIAL_DIM : 0;
  for (unsigned i=0, cidx = START_GC, ridx=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->set_coord_index(cidx);
    _ijoints[i]->set_constraint_index(ridx);
    cidx += _ijoints[i]->num_dof();
    ridx += _ijoints[i]->num_constraint_eqns();
  }

  // setup explicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx=0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->set_coord_index(cidx);
    _ejoints[i]->set_constraint_index(ridx);
    cidx += _ejoints[i]->num_dof();
    ridx += _ejoints[i]->num_constraint_eqns();
  }

  // point both algorithms to this body
  _crb.set_body(get_this());
  _fsab.set_body(get_this());

  // update link transforms and velocities
  update_link_transforms();
  update_link_velocities();
}

/// Sets the vector links 
void RCArticulatedBody::set_links(const vector<RigidBodyPtr>& links)
{
  // call the parent method to update the link indices, etc.  
  ArticulatedBody::set_links(links);

  // setup the processed vector
  _processed.resize(links.size());

  // check to see whether user's numbering scheme is acceptable
  for (unsigned i=1; i< _links.size(); i++)
  {
    // look for an unknown constraint
    BOOST_FOREACH(const RigidBody::InnerJointData& ijd, _links[i]->get_inner_joints_data())
    {
      JointPtr joint(ijd.inner_joint);

      // look for an unknown constraint type -- that means joints have not been processed
      if (joint->get_constraint_type() == Joint::eUnknown)
        return;
    }

    // no unknown constraint; look for an implicit constraint
    if (!_links[i]->get_inner_joint_implicit())
      throw std::runtime_error("Nonzero link does not have an inner implicit joint!");
  }
}

/// Appends a joint to the list of joints
void RCArticulatedBody::set_joints(const vector<JointPtr>& joints)
{
  // clear the vectors of joints
  _ijoints.clear();
  _ejoints.clear();

  // don't do anything if there are no joints
  if (joints.empty())
  {
    _joints.clear();
    return;
  }

  // find the base link and setup a processed map
  map<RigidBodyPtr, bool> processed;
  RigidBodyPtr base;
  for (unsigned i=0; i< joints.size(); i++)
  {
    RigidBodyPtr inboard = joints[i]->get_inboard_link();
    RigidBodyPtr outboard = joints[i]->get_outboard_link();
    RigidBodyPtr parent = inboard->get_parent_link();
    if (!parent)
    {
      if (base && inboard != base)
        throw std::runtime_error("Multiple base links detected!");
      base = inboard;
    }
    processed[inboard] = false;
    processed[outboard] = false;
  }

  // if there is no clearly defined base link, there are no links *and* there
  // are joints, we can't do anything
  if (!base)
  {
    if (_links.empty())
      throw std::runtime_error("Unable to construct RCArticulatedBody with loops w/o links set first!");
    else 
      base = _links.front();
  }

  // verify that base is first link
  if (base->get_index() != 0)
    throw std::runtime_error("Base is not first link!");

  // start processed at the base link
  queue<RigidBodyPtr> link_queue;
  link_queue.push(base);
  while (!link_queue.empty())
  {
    // get the link at the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();

    // if the link has already been processed, no need to process it again
    if (processed[link])
      continue;

    // get all outer joints for this link
    const list<RigidBody::OuterJointData>& outer = link->get_outer_joints_data();
    BOOST_FOREACH(const RigidBody::OuterJointData& ojd, outer)
    {
      // get the joint
      JointPtr joint(ojd.outer_joint);

      // see whether the child has already been processed
      RigidBodyPtr child(ojd.child);
      if (processed[child])
        _ejoints.push_back(joint);
      else
      {
        link_queue.push(child);
        _ijoints.push_back(joint);
      } 
    }

    // indicate that the link has been processed
    processed[link] = true;
  }

  // recalculate the implicit joint degrees-of-freedom of this body
  _n_joint_DOF_implicit = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    _n_joint_DOF_implicit += _ijoints[i]->num_dof();

  // mark joints as the correct type
  for (unsigned i=0; i< _ijoints.size(); i++)
    _ijoints[i]->set_constraint_type(Joint::eImplicit);
  for (unsigned i=0; i< _ejoints.size(); i++)
    _ejoints[i]->set_constraint_type(Joint::eExplicit);

  // call the articulated body method
  ArticulatedBody::set_joints(joints);
}

/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as 
 * well as the base momentum; therefore, the derivative of the state vector is 
 * composed of the joint-space accelerations and base forces (and torques).
 */
VectorN& RCArticulatedBody::get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, VectorN& ga)
{
  const unsigned GC_AA_DIM = 6, GC_ROD_DIM = 7;

  // look for easy case
  if (!is_enabled())
    return ga.resize(0);

  // resize the state-derivative vector
  ga.resize(num_generalized_coordinates(gctype));
  
  // NOTE: we do not lock the mutex here, because forward dynamics computation
  // and inverse dynamics computation will not modify any of the velocities or accelerations utilized

  // setup the generalized acceleration for the base (if any) 
  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    base->get_generalized_acceleration(gctype, _workv);
    if (gctype == DynamicBody::eRodrigues)
    {
      // note: we have to put base linear acceleration in the base frame
      Vector3 base_xdd(_workv[0], _workv[1], _workv[2]);
      base_xdd = base->get_transform().transpose_mult_vector(base_xdd);
      _workv.set_sub_vec(0, base_xdd);
    }
    ga.set_sub_vec(0, _workv);
  }

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // setup the state for the joints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the starting index for this joint
    unsigned idx = _ijoints[i]->get_coord_index() + add;

    // get the qdd vectors for this joint
    ga.set_sub_vec(idx, _ijoints[i]->qdd);
  }

  return ga;
}

/// Updates the transforms of the links based on the current joint positions
/**
 * \note all transforms for the body are recalculated, not just joint values that have changed
 */
void RCArticulatedBody::update_link_transforms()
{
  if (_links.empty() || _joints.empty())
    return;

  // get the base link
  RigidBodyPtr base = _links.front();
  
  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_transforms() entered" << std::endl;
  
  // add all joints of children of the base link to the transform queue for processed
  queue<JointPtr> transform_queue;
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links));
  BOOST_FOREACH(RigidBodyPtr rb, child_links) 
    transform_queue.push(rb->get_inner_joint_implicit());

  // reset processed vector
  for (unsigned i=0; i< _links.size(); i++)
    _processed[i] = false;
  _processed[base->get_index()] = true;
  
  // now set the transforms for all links beside the base
  while (!transform_queue.empty())
  {
    // get the joint
    JointPtr joint = transform_queue.front();
    transform_queue.pop();
  
    // get the outboard link (and its parent) for this joint
    RigidBodyPtr outboard(joint->get_outboard_link());
    RigidBodyPtr parent(joint->get_inboard_link());

    // we need to compute the transform oTp * pTj * jTx * xTc = oTc
    // where o = global frame
    //       p = parent link frame
    //       j = joint frame
    //       x = induced joint transformation
    //       c = child link frame

    // compute T_parent = oTp * pTj
    const Matrix4& Px = parent->get_transform();
    const Vector3& v_parent_to_joint = parent->get_outer_joint_data(outboard).com_to_joint_vec;
    Vector3 Tparent_x = Px.mult_point(v_parent_to_joint);
    Matrix4 Tparent = Px;
    Tparent.set_translation(Tparent_x);

    // compute T_local = jTx * xTc 
    const Matrix4& Jx = joint->get_transform();
    const Vector3& v_joint_to_child = outboard->get_inner_joint_data(parent).joint_to_com_vec_jf;
    Vector3 Tlink_x = Jx.mult_point(v_joint_to_child);
    Matrix4 Tlink = Jx;
    Tlink.set_translation(Tlink_x);

    // set the link transform and mark it as processed
    outboard->set_transform(Tparent * Tlink);
    _processed[outboard->get_index()] = true;

    // add joints of all child links to the queue for processed
    child_links.clear();
    outboard->get_child_links(std::back_inserter(child_links));
    BOOST_FOREACH(RigidBodyPtr rb, child_links) 
      if (!_processed[rb->get_index()])
        transform_queue.push(rb->get_inner_joint_implicit());

    FILE_LOG(LOG_DYNAMICS) << "_transform for link " << outboard->id << std::endl << outboard->get_transform() << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_transforms() exited" << std::endl;
}

/// Updates the link velocities
void RCArticulatedBody::update_link_velocities()
{
  Matrix3 Ri, Rim1;
  queue<RigidBodyPtr> link_queue;
  vector<SVector6> v(_links.size());

  // look for easy exit
  if (_links.empty() || _joints.empty() || !is_enabled())
    return;

  // get the base link
  RigidBodyPtr base = _links.front();

  // get the spatial velocity (in link coordinates) of the base link
  v.front() = base->get_spatial_velocity(eLink);
  
  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_velocities() entered" << std::endl;
  
  // add all children of the base to the link queue
  list<RigidBodyPtr> child_links;
  base->get_child_links(std::back_inserter(child_links));
  BOOST_FOREACH(RigidBodyPtr rb, child_links)
    link_queue.push(rb);
 
  // reset processed vector
  for (unsigned i=0; i< _links.size(); i++)
    _processed[i] = false;
  _processed[base->get_index()] = true;
 
  // propagate link velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    RigidBodyPtr link = link_queue.front();
    link_queue.pop();  
    unsigned idx = link->get_index();

    // and push all children of the link onto the queue
    list<RigidBodyPtr> child_links;
    link->get_child_links(std::back_inserter(child_links));
    BOOST_FOREACH(RigidBodyPtr rb, child_links)
      if (!_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the parent link
    JointPtr joint(link->get_inner_joint_implicit());
    RigidBodyPtr parent(joint->get_inboard_link());
    unsigned pidx = parent->get_index();

    // get the transform from the parent link to this link
    SpatialTransform X_i_im1 = link->get_spatial_transform_forward();

    // get the parent velocity in the link frame
    const SVector6& vparent = X_i_im1.transform(v[pidx]);
  
    // get the link spatial axis
    const SMatrix6N& s = joint->get_spatial_axes(eLink); 

    // determine the link velocity due to the parent velocity + joint velocity
    v[idx] = s.mult(joint->qd) + vparent;

    // set the link linear and angular velocities
    link->set_spatial_velocity(v[idx], eLink);
    _processed[link->get_index()] = true;

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << link->id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- s (local): " << joint->get_spatial_axes(eLink) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent velocity: " << vparent << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- velocity (w/o parent): " << (v[idx] - vparent) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- velocity (w/parent): " << v[idx] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- new linear velocity: " << link->get_lvel() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- new angular velocity: " << link->get_avel() << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::update_link_velocities() exited" << std::endl;
}

/// Calculates the column of a Jacobian matrix for a floating base with respect to a given point
/**
 * \param point the point in 3D space which the Jacobian is calculated against
 * \return a pointer to a 6x6 matrix; top three dimensions will be linear
 *         velocity and bottom three dimensions will be angular velocity
 */
MatrixN& RCArticulatedBody::calc_jacobian_floating_base(const Vector3& point, MatrixN& J)
{
  const unsigned SIX_D = 6, A = 0, B = 1, C = 2, D = 3;

  // init the base Jacobian
  J.resize(SIX_D, SIX_D);

  // get the base transform
  const Matrix4& base_transform = _links.front()->get_transform();

  // get the offset from the base center-of-mass
  Vector3 offset = point - base_transform.get_translation();

  // compute the column-wise cross product of the base orientation and the offset
  Matrix3 RB;
  base_transform.get_rotation(&RB);
  Vector3 RBa, RBb, RBc;
  RB.get_column(A, RBa.begin());
  RB.get_column(B, RBb.begin());
  RB.get_column(C, RBc.begin());
  Vector3 cross1 = Vector3::cross(RBa, offset);
  Vector3 cross2 = Vector3::cross(RBb, offset);
  Vector3 cross3 = Vector3::cross(RBc, offset);

  // make a 3x3 matrix from the cross products
  Matrix3 cross;
  cross.set_column(A, cross1);
  cross.set_column(B, cross2);
  cross.set_column(C, cross3);

  // set the upper left submatrix of the Jacobian
  J.set_sub_mat(A,A, IDENTITY_3x3);

  // set the lower left submatrix of the Jacobian
  J.set_sub_mat(D,A, ZEROS_3x3);

  // set the upper right submatrix of the Jacobian
  J.set_sub_mat(A,D, cross);  

  // set the lower right submatrix of the Jacobian
  J.set_sub_mat(D,D, RB);

  return J;
}

/// Calculates the Jacobian for the current robot configuration at a given point and with respect to a given link
MatrixN& RCArticulatedBody::calc_jacobian(const Vector3& p, RigidBodyPtr link, MatrixN& J)
{
  const unsigned NSPATIAL = 6;

  // resize the Jacobian
  J.set_zero(NSPATIAL, num_generalized_coordinates(DynamicBody::eAxisAngle));

  if (is_floating_base())
  {
    // calculate the floating base
    calc_jacobian_floating_base(p, _workM2);

    // setup the floating base
    J.set_sub_mat(0,0,_workM);
  }

  // calculate all relevant columns
  RigidBodyPtr base = get_base_link();
  while (link != base)
  {
    JointPtr joint = link->get_inner_joint_implicit();
    calc_jacobian_column(joint, p, _workM2); 
    J.set_sub_mat(0,joint->get_coord_index(), _workM2);
    link = link->get_parent_link();
  }

  return J;
}

/// Calculates the column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \param base_transform the transform to use for the base
 * \param q a mapping from joints to joint positions to use; joints without
 *        positions specified will be set to their current values 
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 * \note this method works by changing the joint values, recomputing all link
 *       transforms, computing the Jacobian, restoring the old joint values and
 *       link transforms; thus, it will be slower than a special purpose method
 */
MatrixN& RCArticulatedBody::calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const map<JointPtr, VectorN>& q, MatrixN& Jc)
{
  // store current joint values
  map<JointPtr, VectorN> currentQ;
  for (unsigned i=0; i< _ijoints.size(); i++)
    currentQ[_ijoints[i]].copy_from(_ijoints[i]->q);

  // overwrite current joint values
  for (map<JointPtr, VectorN>::const_iterator i = q.begin(); i != q.end(); i++)
    i->first->q = i->second;
  
  // update link transforms
  update_link_transforms();

  // compute and store the Jacobian
  calc_jacobian_column(joint, point, Jc);

  // restore joint values
  for (map<JointPtr, VectorN>::const_iterator i = currentQ.begin(); i != currentQ.end(); i++)
    i->first->q.copy_from(i->second);

  // restore transforms
  update_link_transforms();

  return Jc;
}

/// Calculates the column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \param base_transform the transform to use for the base
 * \param q a mapping from joints to joint positions to use; joints without
 *        positions specified will be set to their current values 
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 * \note this method works by changing the joint values, recomputing all link
 *       transforms, computing the Jacobian, restoring the old joint values and
 *       link transforms; thus, it will be slower than a special purpose method
 */
MatrixN RCArticulatedBody::calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const map<JointPtr, VectorN>& q)
{
  // compute and store the Jacobian
  MatrixN Jc;
  calc_jacobian_column(joint, point, base_transform, q, Jc);

  return Jc;
}

/// Calculates column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 */
MatrixN RCArticulatedBody::calc_jacobian_column(JointPtr joint, const Vector3& point)
{
  MatrixN Jc;
  calc_jacobian_column(joint, point, Jc);
  return Jc;
}

/// Calculates column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 */
MatrixN& RCArticulatedBody::calc_jacobian_column(JointPtr joint, const Vector3& point, MatrixN& Jc)
{
  const unsigned SPATIAL_DIM = 6;

  // NOTE: spatial algebra provides us with a simple means to compute the
  // Jacobian of a joint with respect to a point.  The spatial axis of the
  // joint transforms joint velocity to spatial velocity. The spatial
  // transform is then used to transform the spatial velocity into a
  // convenient frame: the orientation of that frame is aligned with the
  // global frame, and the origin of the frame is aligned with the desired
  // point.  Applying the spatial transform from the spatial axis (either in
  // link or global frame) to the new frame will give us the desired vector.

  // get the spatial axis of the joint 
  const SMatrix6N& sa = joint->get_spatial_axes(eLink);
  Jc.resize(SPATIAL_DIM, sa.columns());

  // compute a spatial transformation using the point as the target frame
  Matrix4 target(&IDENTITY_3x3, &point);
  RigidBodyPtr outboard(joint->get_outboard_link());
  SpatialTransform X_0i_i(outboard->get_transform(), target);

  // calculate the Jacobian column
  for (unsigned i=0; i< sa.columns(); i++)
  {
    SVector6 sx = X_0i_i.transform(sa.get_column(i));
    Vector3 top = sx.get_upper();
    Vector3 bot = sx.get_lower();
    Jc(0,i) = bot[0];
    Jc(1,i) = bot[1];
    Jc(2,i) = bot[2];
    Jc(3,i) = top[0];
    Jc(4,i) = top[1];
    Jc(5,i) = top[2];
  }

  return Jc; 
}

/// Invalidates the position-computed variables
void RCArticulatedBody::invalidate_positions()
{
  ArticulatedBody::invalidate_positions();
  _fsab.invalidate_position_data();
  _crb.invalidate_position_data();
  _invM_valid = false;
}

/// Invalidates the position-computed variables
void RCArticulatedBody::invalidate_velocities()
{
  ArticulatedBody::invalidate_velocities();
  _fsab.invalidate_velocity_data();
  _crb.invalidate_velocity_data();
}

/// The signum function
Real RCArticulatedBody::sgn(Real x)
{
  if (x < -NEAR_ZERO)
    return (Real) -1.0;
  else if (x > NEAR_ZERO)
    return (Real) 1.0;
  else
    return (Real) 0.0;
}

/// Computes the forward dynamics
/**
 * Given the joint positions and velocities, joint forces, and external 
 * forces on the links, determines the joint and link accelerations as well as
 * floating base accelerations (if applicable).  The joint 
 * accelerations are stored in the individual joints and the link accelerations
 * are stored in the individual links.
 * \note only computes the forward dynamics if the state-derivative is no longer valid
 */
void RCArticulatedBody::calc_fwd_dyn(Real dt)
{
  if (!is_enabled())
    return;

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() entered" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  computing forward dynamics in ";
  if (computation_frame_type == eGlobal)
    FILE_LOG(LOG_DYNAMICS) << "global ";
  else
    FILE_LOG(LOG_DYNAMICS) << "link ";
  FILE_LOG(LOG_DYNAMICS) << "coordinate system" << std::endl;

  // see whether to use advanced friction model 
  if (use_advanced_friction_model)
  {
    calc_fwd_dyn_advanced_friction(dt);
    return;
  }

  // do forward dynamics with old Coulomb model
  for (unsigned i=0; i< _joints.size(); i++)
  {
    _workv.resize(_joints[i]->num_dof());
    for (unsigned j=0; j< _joints[i]->num_dof(); j++)
    {
      _workv[j] = (Real) -1.0;
      if (_joints[i]->qd[j] < (Real) 0.0)
        _workv[j] = (Real) 1.0;
      _workv[j] *= _joints[i]->mu_fc;
      _workv[j] -= _joints[i]->qd[j]*_joints[i]->mu_fv;
    }
    _joints[i]->add_force(_workv);
  }

  // if there are explicit joints, we must do the model with loops
  if (!_ejoints.empty())
  {
    calc_fwd_dyn_loops();
    return;
  } 

  // use the proper dynamics algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      _fsab.calc_fwd_dyn();
      break;

    case eCRB:
      _crb.calc_fwd_dyn();
      break;

    default:
      assert(false);
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() exited" << std::endl;
}

/// Computes the forward dynamics with loops
void RCArticulatedBody::calc_fwd_dyn_loops()
{
  Real Cx[6];

  // look for easy exit
  if (!is_enabled())
    return;

  // get the generalized velocity, generalized forces, and inverse generalized 
  // inertia matrix
  get_generalized_velocity(eAxisAngle, _v);
  get_generalized_forces(eAxisAngle, _fext);

  // determine how many explicit constraint equations
  unsigned N_EXPLICIT_CONSTRAINT_EQNS = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    N_EXPLICIT_CONSTRAINT_EQNS = _ejoints[i]->num_constraint_eqns();

  // evaluate explicit constraints
  _C.resize(N_EXPLICIT_CONSTRAINT_EQNS);
  for (unsigned i=0, r=0, s=0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->evaluate_constraints(Cx);
    for (unsigned j=0; j< _ejoints[i]->num_constraint_eqns(); j++)
      _C[r++] = Cx[j];
  }

  // get the explicit constraint Jacobian and its time derivative
  determine_explicit_constraint_jacobian(_Jx);
  determine_explicit_constraint_jacobian_dot(_Jx_dot);

  // compute Jx * v and Jx_dot * v
  _Jx.mult(_v, _Jx_v);
  _Jx_dot.mult(_v, _Jx_dot_v);

  // get movement Jacobian for explicit constraints and compute velocities
  determine_explicit_constraint_movement_jacobian(_Dx);
  _Dx.mult(_v, _Dx_v);
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    _Dx_v.get_sub_vec(k, k+_ejoints[i]->num_dof(), _ejoints[i]->qd);
    k += _ejoints[i]->num_dof();
  }

  // add in explicit actuator forces
  _beta_x.resize(_Dx.rows());
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->get_scaled_force(_workv);
    for (unsigned j=0; j< _ejoints[i]->num_dof(); j++)
      _beta_x[k++] = _workv[j];
  }
  _Dx.transpose_mult(_beta_x, _workv);
  _fext += _workv;

  // compute the constraint forces 
  solve_generalized_inertia(DynamicBody::eAxisAngle, _fext, _iM_fext);
  _Jx.mult(_iM_fext, _alpha_x) += _Jx_dot_v;
  _Jx.mult(_v, _workv) *= ((Real) 2.0 * b_alpha);
  _alpha_x += _workv;
  _C *= (b_beta*b_beta);
  _alpha_x += _C;
  solve_generalized_inertia_transpose(DynamicBody::eAxisAngle, _Jx, _workM);
  _Jx.mult(_workM, _Jx_iM_JxT);
  try
  {
    LinAlg::solve_LS_fast1(_Jx_iM_JxT, _alpha_x);
  }
  catch (NumericalException e)
  {
    _Jx.mult(_workM, _Jx_iM_JxT);
    LinAlg::solve_LS_fast2(_Jx_iM_JxT, _alpha_x);
  }

  // compute generalized acceleration
  _fext -= _Jx.transpose_mult(_alpha_x, _workv);
  solve_generalized_inertia(DynamicBody::eAxisAngle, _fext, _workv2);
  set_generalized_acceleration(DynamicBody::eAxisAngle, _workv2);
}

/// Computes the forward dynamics
/**
 * Given the joint positions and velocities, joint forces, and external 
 * forces on the links, determines the joint and link accelerations as well as
 * floating base accelerations (if applicable).  The joint 
 * accelerations are stored in the individual joints and the link accelerations
 * are stored in the individual links.
 * \note only computes the forward dynamics if the state-derivative is no longer valid
 */
void RCArticulatedBody::calc_fwd_dyn_advanced_friction(Real dt)
{
  Real Cx[6];
  const Real QP_TOL = std::sqrt(NEAR_ZERO);  // tolerance for solving the QP
  SAFESTATIC MatrixN iM, X, Jx_X_JxT, Dx_X_DxT, Jx_iM_JxT;
  SAFESTATIC MatrixN Jx_dot, Y, X_JxT, X_DxT, Dx_X_JxT, Jx_iM, Jx_Y, A, RG;
  SAFESTATIC MatrixN Jx_iM_DxT;
  SAFESTATIC VectorN v, fext, C, alpha_x, beta_x, Dx_v, Jx_v, Jx_dot_v, _workv;
  SAFESTATIC VectorN vddt_plus_Y_fext, Y_fext, x, ff, delta, iM_fext, a;

  // look for quick exit
  if (!is_enabled())
    return;

  // get the generalized velocity, generalized forces, and inverse generalized 
  // inertia matrix
  get_generalized_velocity(eAxisAngle, v);
  get_generalized_forces(eAxisAngle, fext);
  get_generalized_inertia_inverse(eAxisAngle, iM);

  // setup X and Y
  const unsigned SPATIAL_DIM = 6;
  const unsigned START_GC = (_floating_base) ? SPATIAL_DIM : 0;
  iM.get_sub_mat(START_GC, iM.rows(), START_GC, iM.columns(), X);
  iM.get_sub_mat(0, iM.rows(), START_GC, iM.columns(), Y);

  // determine how many explicit constraint equations
  unsigned N_EXPLICIT_CONSTRAINT_EQNS = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    N_EXPLICIT_CONSTRAINT_EQNS = _ejoints[i]->num_constraint_eqns();

  // evaluate explicit constraints
  C.resize(N_EXPLICIT_CONSTRAINT_EQNS);
  for (unsigned i=0, r=0, s=0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->evaluate_constraints(Cx);
    for (unsigned j=0; j< _ejoints[i]->num_constraint_eqns(); j++)
      C[r++] = Cx[j];
  }

  // get the explicit constraint Jacobian and its time derivative
  determine_explicit_constraint_jacobian(_Jx);
  determine_explicit_constraint_jacobian_dot(Jx_dot);

  // compute Jx * v and Jx_dot * v
  _Jx.mult(v, Jx_v);
  Jx_dot.mult(v, Jx_dot_v);

  // get movement Jacobian for explicit constraints and compute velocities
  determine_explicit_constraint_movement_jacobian(_Dx);
  _Dx.mult(v, Dx_v);
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    Dx_v.get_sub_vec(k, k+_ejoints[i]->num_dof(), _ejoints[i]->qd);
    k += _ejoints[i]->num_dof();
  }

  // add in explicit actuator forces
  beta_x.resize(_Dx.rows());
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->get_scaled_force(_workv);
    for (unsigned j=0; j< _ejoints[i]->num_dof(); j++)
      beta_x[k++] = _workv[j];
  }
  _Dx.transpose_mult(beta_x, _workv);
  fext += _workv;

  // determine the number of kinematic loops and number of explicit joint DOF
  const unsigned N_LOOPS = _ejoints.size();
  unsigned N_EXPLICIT_DOF = 0, N_IMPLICIT_DOF = 0;
  for (unsigned i=0; i< _joints.size(); i++)
  {
    if (_joints[i]->get_constraint_type() == Joint::eExplicit)
      N_EXPLICIT_DOF += _joints[i]->num_dof();
    else if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      N_IMPLICIT_DOF += _joints[i]->num_dof();
    else
      assert(false);
  }

  // setup some constants
  const unsigned N_JOINT_DOF = N_IMPLICIT_DOF + N_EXPLICIT_DOF;
  const unsigned FF_START = 0;
  const unsigned ALPHAX_START = N_IMPLICIT_DOF;
  const unsigned BETA_START = N_IMPLICIT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;

  // precompute some things
  X.mult_transpose(_Jx, X_JxT);
  X.mult_transpose(_Dx, X_DxT);
  _Jx.mult(X_JxT, Jx_X_JxT);
  _Dx.mult(X_JxT, Dx_X_JxT);
  _Dx.mult(X_DxT, Dx_X_DxT);

  // setup data needed for convex optimization
  ABFwdDynOptData copt_data;
  copt_data.N_IMPLICIT_DOF = N_IMPLICIT_DOF;
  copt_data.N_EXPLICIT_CONSTRAINT_EQNS = N_EXPLICIT_CONSTRAINT_EQNS;
  copt_data.N_JOINT_DOF = N_JOINT_DOF;
  copt_data.N_LOOPS = N_LOOPS;
  copt_data.Dx.copy_from(_Dx);
  copt_data.fext.copy_from(fext);

  // first compute the nullspace and the homogeneous solution
  const unsigned NN = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS + N_LOOPS;
  MatrixN& R = copt_data.R;
  VectorN& z = copt_data.z;
  iM.mult(fext, iM_fext);
  if (N_EXPLICIT_CONSTRAINT_EQNS > 0)
  {
    // setup the right hand side
    _Jx.mult(iM_fext, z) *= dt;
    z += Jx_dot_v;
    _Jx.mult(v, _workv) *= ((Real) 2.0 * b_alpha);
    z += _workv;
    C *= (b_beta*b_beta);
    z += C;

    // compute the nullspace
    _Jx.mult(Y, Jx_Y);
    _Jx.mult(iM, Jx_iM);
    Jx_iM.mult_transpose(_Jx, Jx_iM_JxT);
    Jx_iM.mult_transpose(_Dx, Jx_iM_DxT);
    A.resize(N_EXPLICIT_CONSTRAINT_EQNS, N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS);
    A.set_sub_mat(0, 0, Jx_Y);
    A.set_sub_mat(0, N_IMPLICIT_DOF, Jx_iM_JxT);
    A.set_sub_mat(0, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, Jx_iM_DxT);
    A *= dt;
    LinAlg::nullspace(A, RG);

    // now add in parts for delta (unrelated to nullspace)
    R.set_zero(NN, RG.columns()+N_LOOPS);
    R.set_sub_mat(0,0,RG);
    for (unsigned i=0; i< N_LOOPS; i++)
      R(N_JOINT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+i,RG.columns()+i) = (Real) 1.0;

    // solve for explicit constraint forces 
    A.copy_from(Jx_iM_JxT);
    A *= dt;
    _workv.copy_from(z);
    try
    {
      LinAlg::solve_LS_fast1(A, _workv);
    }
    catch (NumericalException e)
    {
      A.copy_from(Jx_iM_JxT);
      A *= dt;
      LinAlg::solve_LS_fast2(A, _workv);
    } 
    z.set_zero(NN);
    z.set_sub_vec(ALPHAX_START, _workv);
  }
  else
  {
    z.set_zero(N_JOINT_DOF);
    R.set_zero(N_JOINT_DOF, N_JOINT_DOF);
    for (unsigned i=0; i< N_JOINT_DOF; i++)
      R(i,i) = (Real) 1.0;
  }

  // setup components of z and R to make things faster for gradient and Hessian
  // calculations
  z.get_sub_vec(FF_START, FF_START+N_IMPLICIT_DOF, copt_data.zff);
  z.get_sub_vec(BETA_START, BETA_START+N_EXPLICIT_DOF, _workv);
  _Dx.transpose_mult(_workv, copt_data.zbetax);
  R.get_sub_mat(FF_START, FF_START+N_IMPLICIT_DOF, 0, R.columns(), copt_data.Rff);
  R.get_sub_mat(BETA_START, BETA_START+N_EXPLICIT_DOF, 0, R.columns(), A);
  _Dx.transpose_mult(A, copt_data.DxTRbetax);

  // setup the data for the optimization problem
  const unsigned N = R.columns(); 
  const unsigned M = N_JOINT_DOF + N_LOOPS*2;
  OptParams copt_params(N, M, 0, calc_fwd_dyn_f0, calc_fwd_dyn_fx, NULL, calc_fwd_dyn_grad0, calc_fwd_dyn_cJac, NULL, calc_fwd_dyn_hess);
  copt_params.max_iterations = N*N*N;
  copt_params.data = (void*) &copt_data;

  // determine loops and compute Z matrices 
  vector<unsigned>& loop_indices = copt_data.loop_indices;
  vector<vector<unsigned> > loop_links;
  ArticulatedBody::find_loops(loop_indices, loop_links); 
  ArticulatedBody::compute_Z_matrices(loop_indices, loop_links, copt_data.Zd, copt_data.Z1d, copt_data.Z);

  // verify that joint indices are as we expect
  #ifndef NDEBUG
  for (unsigned i=0; i< _joints.size(); i++)
    assert(_joints[i]->get_index() == i);
  #endif

  // setup true indices: converts from an [implicit explicit] DOF to the
  // joint's index that contains that DOF
  vector<unsigned>& true_indices = copt_data.true_indices;
  true_indices.resize(N_JOINT_DOF);
  for (unsigned i=0, ki=0, ke=0; i< _joints.size(); i++)
  {
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
     for (unsigned j=0; j< _joints[i]->num_dof(); j++)
       true_indices[ki++] = i;
    else
    {
      assert(_joints[i]->get_constraint_type() == Joint::eExplicit);
      for (unsigned j=0; j< _joints[i]->num_dof(); j++)
        true_indices[N_IMPLICIT_DOF + ke++] = i;
    }
  }

  // initialize mu_c and viscous force vectors
  vector<Real>& mu_c = copt_data.mu_c;
  vector<Real>& visc = copt_data.visc;
  mu_c.resize(N_JOINT_DOF);
  visc.resize(N_JOINT_DOF);

  // setup mu_c and viscous force vectors for implicit joints
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
    for (unsigned j=0; j< _ijoints[i]->num_dof(); j++, k++)
    {
      mu_c[k] = _ijoints[i]->mu_fc*_ijoints[i]->mu_fc;
      Real tmp = _ijoints[i]->mu_fv * std::fabs(_ijoints[i]->qd[j]);
      visc[k] = tmp*tmp;
    }

  // setup mu_c and viscous force vectors for explicit joints
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
    for (unsigned j=0; j< _ejoints[i]->num_dof(); j++, k++)
    {
      mu_c[k+N_IMPLICIT_DOF] = _ejoints[i]->mu_fc * _ejoints[i]->mu_fc;
      Real tmp = _ejoints[i]->mu_fv * std::fabs(_ejoints[i]->qd[j]);
      visc[k+N_IMPLICIT_DOF] = tmp*tmp; 
    }

  // original (w/o nullspace) optimization problem is:
  // 0.5 * (v + h*inv(M)*tau)' * M * (v + h*inv(M)*tau)
  // ==> 0.5 * h^2 tau' inv(M) tau + v' * h * tau + 0.5 * v' * M * v'
  // ==> 0.5 v' M v +
  //     0.5 (h^2 fext' inv(M) fext + h^2 fext' inv(M) ff + 
  //     h^2 fext' inv(M) Jx'ax + h^2 fext' inv(M) Dx'bx) + v' * h * fext +
  //     0.5 (h^2 ff' inv(M) fext + h^2 ff' inv(M) ff + h^2 ff' inv(M) Jx' ax +
  //          h^2 ff' inv(M) Dx' bx) + v' * h * ff +
  //     0.5 (h^2 ax'Jx inv(M) fext + h^2 ax'Jx inv(M) ff + 
  //          h^2 ax'Jx inv(M) Jx' ax + h^2 ax'Jx inv(M) Dx' bx) + 
  //          v' * h * Jx' * ax +
  //     0.5 (h^2 bx'Dx inv(M) fext + h^2 bx'Dx inv(M) ff + 
  //          h^2 bx'Dx inv(M) Jx' ax + h^2 bx'Dx inv(M) Dx' bx) + 
  //          v' * h * Dx' * bx
  // G = | inv(M)     inv(M) Jx'     inv(M) Dx'    0 |
  //     | Jx inv(M)  Jx inv(M) Jx'  Jx inv(M) Dx' 0 |
  //     | Dx inv(M)  Dx inv(M) Jx'  Dx inv(M) Dx' 0 |
  //     | 0          0              0             0 |
  // c = | inv(M) fext + v/h       |
  //     | Jx inv(M) fext + Jx*v/h |
  //     | Dx inv(M) fext + Dx*v/h |
  //     | 0                       |

  // compute the quadratic matrix
  MatrixN& G = copt_data.G;
  G.set_zero(NN, NN);
  G.set_sub_mat(0,0,X); 
  G.set_sub_mat(0,N_IMPLICIT_DOF,X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF,0,X_JxT,true);
  G.set_sub_mat(0,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,X_DxT);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,0,X_DxT,true);
  G.set_sub_mat(N_IMPLICIT_DOF,N_IMPLICIT_DOF,Jx_X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,Dx_X_JxT, true);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,N_IMPLICIT_DOF,Dx_X_JxT);
  G.set_sub_mat(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS,Dx_X_DxT);
  FILE_LOG(LOG_DYNAMICS) << "G (before nullspace): " << std::endl << G;
  R.transpose_mult(G, RG);
  RG.mult(R, G);

  // (z + Ry)G(z + Ry) + (z' + y'R')c = 0.5 y'R'G'Ry + y'(R'Gz + R'c)
  // compute the linear optimization vector
  VectorN& c = copt_data.c;
  c.set_zero(NN);
  Y.mult(fext, Y_fext);
  vddt_plus_Y_fext.copy_from(v) /= dt;
  vddt_plus_Y_fext += Y_fext;
  c.set_sub_vec(0, vddt_plus_Y_fext);
  _Jx.mult(vddt_plus_Y_fext, _workv);
  c.set_sub_vec(N_IMPLICIT_DOF, _workv);
  _Dx.mult(vddt_plus_Y_fext, _workv);
  c.set_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, _workv);
  FILE_LOG(LOG_DYNAMICS) << "c (before nullspace): " << c << std::endl;
  R.transpose_mult(c, _workv);
  c.copy_from(_workv);
  RG.mult(z, _workv);
  c += _workv;
  FILE_LOG(LOG_DYNAMICS) << "G (after nullspace): " << std::endl << G;
  FILE_LOG(LOG_DYNAMICS) << "c (after nullspace): " << std::endl << c;
  
  // setup an optimization vector
  x.set_zero(N);  

  // set delta to a middle value
  for (unsigned i=0,j=R.columns()-1; i< N_LOOPS; i++, j--)
    x[j] = (Real) 0.5;

  // optimize
  if (G.norm_inf() > std::sqrt(NEAR_ZERO))
  {
//  Optimization::optimize_convex_pd(copt_params, x);
    Optimization::sqp(copt_params, x);
  }

  // output debugging data
  FILE_LOG(LOG_DYNAMICS) << "R: " << std::endl << R;
  FILE_LOG(LOG_DYNAMICS) << "z: " << z << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "x: " << x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "v: " << v << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "dt: " << dt << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "fext: " << fext << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "inv(M): " << std::endl << iM;
  FILE_LOG(LOG_DYNAMICS) << "Jx: " << std::endl << _Jx;
  FILE_LOG(LOG_DYNAMICS) << "\\dot{Jx}: " << std::endl << Jx_dot;
  FILE_LOG(LOG_DYNAMICS) << "Dx: " << std::endl << _Dx;

  // now retrieve the necessary variables
  R.mult(x, _workv);
  z += _workv;
  z.get_sub_vec(0, N_IMPLICIT_DOF, ff);
  z.get_sub_vec(N_IMPLICIT_DOF, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, alpha_x);
  z.get_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS, N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+N_EXPLICIT_DOF, beta_x);
  z.get_sub_vec(N_IMPLICIT_DOF+N_EXPLICIT_CONSTRAINT_EQNS+N_EXPLICIT_DOF, z.size(), delta);

  // output results
  FILE_LOG(LOG_DYNAMICS) << " implicit friction forces: " << ff << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " explicit constraint forces: " << alpha_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " explicit friction forces: " << beta_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " delta: " << delta << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " constraint evaluations: " << C << std::endl;

  // setup implicit joint friction forces
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    ff.get_sub_vec(k, k+_ijoints[i]->num_dof(), _ijoints[i]->ff);
    k += _ijoints[i]->num_dof();
  }

  // setup explicit joint friction forces
  for (unsigned i=0, k=0; i< _ejoints.size(); i++)
  {
    beta_x.get_sub_vec(k, k+_ejoints[i]->num_dof(), _ejoints[i]->ff);
    k += _ejoints[i]->num_dof();
  }

  // compute joint constraint forces
  fext += ff;
  fext += _Dx.transpose_mult(beta_x, _workv);
  ArticulatedBody::calc_joint_constraint_forces(loop_indices, delta, copt_data.Zd, copt_data.Z1d, copt_data.Z, fext);

  // compute generalized acceleration
  fext -= _Jx.transpose_mult(alpha_x, _workv);
  iM.mult(fext, a);
  set_generalized_acceleration(DynamicBody::eAxisAngle, a);
  if (LOGGING(LOG_DYNAMICS))
    for (unsigned i=0; i< _joints.size(); i++)
    {
      FILE_LOG(LOG_DYNAMICS) << "lambda for joint " << _joints[i]->id << ": " << _joints[i]->lambda << std::endl;
      FILE_LOG(LOG_DYNAMICS) << "Z: " << std::endl << copt_data.Z[i];
      FILE_LOG(LOG_DYNAMICS) << "Zd: " << std::endl << copt_data.Zd[i];
      FILE_LOG(LOG_DYNAMICS) << "Z1d: " << std::endl << copt_data.Z1d[i];
    }
  FILE_LOG(LOG_DYNAMICS) << " generalized acceleration: " << a << std::endl;

  // determine predicted energy
  if (LOGGING(LOG_DYNAMICS))
  {
    MatrixN M;
    get_generalized_inertia(eAxisAngle, M);
    M.mult(v, _workv);
    FILE_LOG(LOG_DYNAMICS) << " current energy: " << (0.5 * v.dot(_workv)) << std::endl;
    a *= dt;
    v += a;
    M.mult(v, _workv);
    FILE_LOG(LOG_DYNAMICS) << " new energy: " << (0.5 * v.dot(_workv)) << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::calc_fwd_dyn() exited" << std::endl;
}

/// Determines the constraint Jacobian for explicit constraints via impulse application (considerably slower than alternative method)
/*
void RCArticulatedBody::determine_explicit_constraint_movement_jacobian(MatrixN& D)
{
  // determine the number of explicit DOFs 
  unsigned NDOF = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    NDOF += _ejoints[i]->num_dof();
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // resize J 
  D.set_zero(NDOF, NGC);
  MatrixNN q(NGC);

  // init gj and dC
  VectorN gj = VectorN::zero(NGC);
  MatrixN qx(NDOF, NGC);
  VectorN gv(NGC);

  // loop
  for (unsigned i=0; i< NGC; i++)
  {
    // apply the generalized impulse
    if (i > 0)
      gj[i-1] = (Real) 0.0;
    gj[i] = (Real) 1.0;
    apply_generalized_impulse(DynamicBody::eAxisAngle, gj);

    // get the velocity
    get_generalized_velocity(DynamicBody::eAxisAngle, gv);
    q.set_column(i, gv);
    
    // compute the joint velocity
    for (unsigned r=0; r< _ejoints.size(); r++)
    {
      // get the coordinate index
      unsigned cidx = _ejoints[r]->get_coord_index();

      // get the spatial axis
      const SMatrix6N& si = _ejoints[r]->get_spatial_axes(eGlobal);
      RigidBodyPtr inboard = _ejoints[r]->get_inboard_link();
      RigidBodyPtr outboard = _ejoints[r]->get_outboard_link();
      SVector6 vh = inboard->get_spatial_velocity(eGlobal);
      SVector6 vi = outboard->get_spatial_velocity(eGlobal);
// vi = vh + si*qdot, except these are all loops
      SVector6 dv = vh - vi;
      MatrixN si_inv = LinAlg::pseudo_inverse(si);
      VectorN qdot = si_inv.mult(dv);
      for (unsigned j=0; j< _ejoints[r]->num_dof(); j++)
        qx(cidx + j, i) = qdot[j];
    }
  }

  // compute D 
  MatrixNN qinv;
  LinAlg::pseudo_inverse(q, qinv);
  qx.mult(qinv, D);

  // remove all impulses
  gj.set_one().negate();
  apply_generalized_impulse(DynamicBody::eAxisAngle, gj);
}
*/

/// Determines the ndof x ngc Jacobian for explicit constraint movement (ndof is the number of degrees of freedom of the explicit constraints)
void RCArticulatedBody::determine_explicit_constraint_jacobians(const EventProblemData& q, MatrixN& Jx, MatrixN& Dx)
{
  Real Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // determine the total number of explicit constraint DOF 
  unsigned NDOF = 0, NEQ = 0;
  for (unsigned i=0; i< q.constraint_events.size(); i++)
  {
    JointPtr joint = q.constraint_events[i]->constraint_joint;
    if (joint->get_constraint_type() == Joint::eExplicit)
    {
      NDOF += joint->num_dof();
      NEQ += joint->num_constraint_eqns();
    }
  }

  // resize Jx and Dx
  Jx.set_zero(NEQ, NGC);
  Dx.set_zero(NDOF, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the Jacobian for all explicit constraints
  for (unsigned i=0, ii=0, jj=0; i< q.constraint_events.size(); i++)
  {
    // get the joint
    JointPtr joint = q.constraint_events[i]->constraint_joint;

    // if it's not explicit, skip it
    if (joint->get_constraint_type() != Joint::eExplicit)
      continue;

    // it is explicit; see whether it belongs to this body
    if (joint->get_articulated_body() != get_this())
    {
      ii += joint->num_dof();
      jj += joint->num_constraint_eqns();
      continue;
    }

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = joint->get_inboard_link();
    RigidBodyPtr rbo = joint->get_outboard_link();

    // get the outboard body's transform
    const Matrix4& To = rbo->get_transform();

    // get the spatial movement (axes) of this joint (rbo frame)
    const SMatrix6N& si = joint->get_spatial_axes(eLink);
    _workM2.copy_from(si);
    try
    {
      LinAlg::pseudo_inverse(_workM2, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      _workM2.copy_from(si);
      LinAlg::pseudo_inverse(_workM2, LinAlg::svd2);
    }
/*
    // special case: floating base
    if (_floating_base)
    {
      SpatialTransform X(To, IDENTITY_3x3, base->get_position());
      MatrixN::transpose(_workM2, sx);
      X.transform(sx, so);
      D.set_sub_mat(r, 0, so, true); 
    }
*/

    // now, loop over implicit joints
    for (unsigned k=0; k< _ijoints.size(); k++)
    {
      // get the links
      RigidBodyPtr rba = _ijoints[k]->get_inboard_link();
      RigidBodyPtr rbb = _ijoints[k]->get_outboard_link();

      // get the coordinate index of the joint
      unsigned cidx = _ijoints[k]->get_coord_index();

      // get transform for the outboard link
      const Matrix4& rbbT = rbb->get_transform();

      // get the spatial axes for the joint (in rbb's frame)
      const SMatrix6N& sk = _ijoints[k]->get_spatial_axes(eLink);

      // compute the spatial transformation to the outboard link
      SpatialTransform X(rbbT, To);
      X.transform(sk, _so);
      _workM2.mult(_so, _workM); 

      // update Dx
      Dx.set_sub_mat(ii,cidx,_workM);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, so);
//        so.transpose_mult(svo, dot);
//        sub -= dot;

      // iterate over the constraint equations
      for (unsigned j=0; j< joint->num_constraint_eqns(); j++, jj++)
      {
        // get constraint equations for inner link
        joint->calc_constraint_jacobian(DynamicBody::eAxisAngle, rbi, j, Cqi); 

        // setup spatial vectors in rbi-centered frame and rbo-centered frame
        SVector6 svi(Cqi);

        // special case: floating base
        if (_floating_base)
        {
          // setup spatial vector in rbo-centered frame
          joint->calc_constraint_jacobian(DynamicBody::eAxisAngle,rbo, j, Cqo); 
          SVector6 svo(Cqo);

          // setup spatial transform for svi/svo
          SpatialTransform X;
          X.r = base->get_position() - rbi->get_position();
          SVector6 svix = X.transform(svi);
          X.r = base->get_position() - rbo->get_position();
          SVector6 svox = X.transform(svo);

          // update Jx corresponding to the base
          for (unsigned k=0; k< 6; k++) 
            Jx(jj,k) += svix[k] - svox[k];
        } 

        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(sk, _so);
        _so.transpose_mult(svi, _workv);

        // update Jx
        Jx.set_sub_mat(jj,cidx,_workv,true);
      }
    }
  
    // update ii
    ii += joint->num_dof();
  }
}

/// Determines the ndof x ngc Jacobian for explicit constraint movement (ndof is the number of degrees of freedom of the explicit constraints)
void RCArticulatedBody::determine_explicit_constraint_movement_jacobian(MatrixN& D)
{
  Real Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // determine the number of explicit constraint DOF 
  unsigned NDOF = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    NDOF += _ejoints[i]->num_dof();

  // resize D
  D.set_zero(NDOF, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the Jacobian for all explicit constraints
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // get the coordinate index of the joint
    unsigned r = _ejoints[i]->get_coord_index();

    // get the outboard body of this joint and its transform
    RigidBodyPtr rbo = _ejoints[i]->get_outboard_link();
    const Matrix4& To = rbo->get_transform();

    // get the spatial movement (axes) of this joint (rbo frame)
    const SMatrix6N& si = _ejoints[i]->get_spatial_axes(eLink);
    _workM2.copy_from(si);
    try
    {
      LinAlg::pseudo_inverse(_workM2, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      _workM2.copy_from(si);
      LinAlg::pseudo_inverse(_workM2, LinAlg::svd2);
    }
/*
    // special case: floating base
    if (_floating_base)
    {
      SpatialTransform X(To, IDENTITY_3x3, base->get_position());
      MatrixN::transpose(_workM2, sx);
      X.transform(sx, so);
      D.set_sub_mat(r, 0, so, true); 
    }
*/
    // now, loop over implicit joints
    for (unsigned k=0; k< _ijoints.size(); k++)
    {
      // get the outboard link
      RigidBodyPtr rbb = _ijoints[k]->get_outboard_link();
      unsigned idx = _ijoints[k]->get_coord_index();
      const Matrix4& Tb = rbb->get_transform();

      // get the spatial axes for the joint (rbb frame)
      const SMatrix6N& sk = _ijoints[k]->get_spatial_axes(eLink);

      // compute the spatial transformation
      SpatialTransform X(Tb, To);
      X.transform(sk, _so);
      _workM2.mult(_so, _workM); 
      
      // update D
      D.set_sub_mat(r,idx,_workM);
    }
  }
}

/// Determines the constraint Jacobian for explicit constraints via impulse application (considerably slower than alternative method)
/*
void RCArticulatedBody::determine_explicit_constraint_jacobian(MatrixN& J)
{
  Real Cx[6];

  // determine the number of explicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    NEQ += _ejoints[i]->num_constraint_eqns();
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // resize J 
  J.set_zero(NEQ, NGC);
  MatrixNN q(NGC);

  // init gj and dC
  VectorN gj = VectorN::zero(NGC);
  MatrixN dC(NEQ, NGC);
  VectorN gv(NGC);

  // loop
  for (unsigned i=0; i< NGC; i++)
  {
    // apply the generalized impulse
    if (i > 0)
      gj[i-1] = (Real) 0.0;
    gj[i] = (Real) 1.0;
    apply_generalized_impulse(DynamicBody::eAxisAngle, gj);

    // get the velocity
    get_generalized_velocity(DynamicBody::eAxisAngle, gv);
    q.set_column(i, gv);
    
    // compute the constraint dot
    for (unsigned r=0; r< _ejoints.size(); r++)
    {
      unsigned cidx = _ejoints[r]->get_constraint_index();
      _ejoints[r]->evaluate_constraints_dot(Cx);
      for (unsigned j=0; j< _ejoints[r]->num_constraint_eqns(); j++)
        dC(j+cidx, i) = Cx[j];
    }
  }

  // compute J
  MatrixNN qinv;
  LinAlg::pseudo_inverse(q, qinv);
  dC.mult(qinv, J);

  // remove all impulses
  gj.set_one().negate();
  apply_generalized_impulse(DynamicBody::eAxisAngle, gj);
}
*/

/// Determines the constraint Jacobian for explicit constraints
void RCArticulatedBody::determine_explicit_constraint_jacobian(MatrixN& J)
{
  Real Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // determine the number of explicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    NEQ += _ejoints[i]->num_constraint_eqns();

  // resize J 
  J.set_zero(NEQ, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // compute the constraint Jacobian for all explicit constraints
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // get the constraint index for this joint
    unsigned ridx = _ejoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ejoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ejoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ejoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ejoints[i]->calc_constraint_jacobian(DynamicBody::eAxisAngle, rbi, j, Cqi); 
      _ejoints[i]->calc_constraint_jacobian(DynamicBody::eAxisAngle, rbo, j, Cqo); 

      // setup spatial vectors in rbi-centered frame and rbo-centered frame
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++) 
          J(ridx,k) += svix[k] - svox[k];
      }        

      // now, loop over implicit joints
      for (unsigned k=0; k< _ijoints.size(); k++)
      {
        // get the links
        RigidBodyPtr rba = _ijoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ijoints[k]->get_outboard_link();

        // get the coordinate index of the joint
        unsigned cidx = _ijoints[k]->get_coord_index();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_transform();

        // get the spatial axes for the joint (in rbb's frame)
        const SMatrix6N& si = _ijoints[k]->get_spatial_axes(eLink);

        // setup components corresponding to rbi
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si, _so);
        _so.transpose_mult(svi, _workv);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, so);
//        so.transpose_mult(svo, dot);
//        _workM -= dot;

        // update J
        J.set_sub_mat(ridx,cidx,_workv,true);
      }
    }
  }
}

/// Determines the time derivative of the constraint Jacobian for explicit constraints
/**
 * Because this matrix is the product of two matrices, each of which is a 
 * function of time, we have to add the results together.
 */
void RCArticulatedBody::determine_explicit_constraint_jacobian_dot(MatrixN& J) 
{
  Real Cqi[6], Cqo[6];
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);

  // determine the number of explicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    NEQ += _ejoints[i]->num_constraint_eqns();

  // resize J 
  J.set_zero(NEQ, NGC);

  // get the base link
  RigidBodyPtr base = get_base_link();

  // Jx * \dot{Jy}
  // compute the constraint Jacobian for all explicit constraints
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // get the starting constraint index for this joint
    unsigned ridx = _ejoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ejoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ejoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ejoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ejoints[i]->calc_constraint_jacobian(DynamicBody::eAxisAngle, rbi, j, Cqi); 
      _ejoints[i]->calc_constraint_jacobian(DynamicBody::eAxisAngle, rbo, j, Cqo); 

      // setup spatial vectors
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++) 
          J(ridx,k) += svix[k] - svox[k];
      }        

      // now, loop over implicit joints
      for (unsigned k=0; k< _ijoints.size(); k++)
      {
        // get the coordinate index for the joint
        unsigned cidx = _ijoints[k]->get_coord_index();

        // get the links
        RigidBodyPtr rba = _ijoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ijoints[k]->get_outboard_link();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_transform();

        // get the time derivative of the spatial axes for the joint
        const SMatrix6N& si_dot = _ijoints[k]->get_spatial_axes_dot(eLink);

        // setup components corresponding to rbi 
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si_dot, _so);
        _so.transpose_mult(svi, _workv);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si_dot, _so);
//        _so.transpose_mult(svo, dot);
//        sub -= dot;

        // update J
        J.set_sub_mat(ridx,cidx,_workv,true);
      }
    }
  }

  // \dot{Jx} * Jy
  // compute the constraint Jacobian for all explicit constraints
  for (unsigned i=0, r=0; i< _ejoints.size(); i++)
  {
    // get the starting constraint index for this joint
    unsigned ridx = _ejoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    RigidBodyPtr rbi = _ejoints[i]->get_inboard_link();
    RigidBodyPtr rbo = _ejoints[i]->get_outboard_link();

    // get the constraint equations
    for (unsigned j=0; j< _ejoints[i]->num_constraint_eqns(); j++, ridx++)
    {
      // get constraint equations for inner and outer links
      _ejoints[i]->calc_constraint_jacobian_dot(DynamicBody::eAxisAngle, rbi, j, Cqi); 
      _ejoints[i]->calc_constraint_jacobian_dot(DynamicBody::eAxisAngle, rbo, j, Cqo); 

      // setup spatial vectors
      SVector6 svi(Cqi), svo(Cqo);

      // special case: floating base
      if (_floating_base)
      {
        RigidBodyPtr base = get_base_link();
        SpatialTransform X;
        X.r = base->get_position() - rbi->get_position();
        SVector6 svix = X.transform(svi);
        X.r = base->get_position() - rbo->get_position();
        SVector6 svox = X.transform(svo);
        for (unsigned k=0; k< 6; k++) 
          J(ridx,k) += svix[k] - svox[k];
      }        

      // now, loop over implicit joints
      for (unsigned k=0; k< _ijoints.size(); k++)
      {
        // get the coordinate index
        unsigned cidx = _ijoints[k]->get_coord_index();

        // get the links
        RigidBodyPtr rba = _ijoints[k]->get_inboard_link();
        RigidBodyPtr rbb = _ijoints[k]->get_outboard_link();

        // get transform for the outboard link
        const Matrix4& rbbT = rbb->get_transform();

        // get the spatial axes for the joint
        const SMatrix6N& si = _ijoints[k]->get_spatial_axes(eLink);

        // setup components corresponding to rbi
        // transform si from rbb's frame to rbi centered frame
        SpatialTransform X(rbbT, IDENTITY_3x3, rbi->get_position());
        X.transform(si, _so);
        _so.transpose_mult(svi, _workv);

        // setup components corresponding to rbo
        // transform si from rbb's frame to rbo centered frame
//        SpatialTransform Y(rbbT, IDENTITY_3x3, rbo->get_position());
//        Y.transform(si, _so);
//        _so.transpose_mult(svo, dot);
//        sub -= dot;

        // update J
        for (unsigned b=0; b< _workv.size(); b++)
            J(ridx,cidx+b) += _workv[b];
      }
    }
  }
}

/// Sets the generalized acceleration for this body
void RCArticulatedBody::set_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& a)
{
  if (!is_enabled())
    return;

  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    if (gctype == DynamicBody::eRodrigues)
    {
      const Quat& q = base->get_orientation();
      base->set_laccel(base->get_transform().mult_vector(Vector3(a[0], a[1], a[2])));
      base->set_aaccel(q.G_mult(a[3], a[4], a[5], a[6]) * (Real) 2.0);
    }
    else
    {
      assert(gctype == DynamicBody::eAxisAngle);
      base->set_laccel(Vector3(a[0], a[1], a[2]));
      base->set_aaccel(Vector3(a[3], a[4], a[5]));
    }
  }

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // set joint accelerations
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    a.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _ijoints[i]->qdd);
  }
}

/// Applies a linear impulse at the specified point and propagates it up the articulated body chain
/**
 * \param j the linear impulse
 * \param k the angular impulse
 * \param contact_point the point of collision
 * \param link the link that the impulse is applied to
 * \param fsab_algo if non-null, already computed quanities from dynamics
 *         algorithm will be used; otherwise, forward dynamics algorithm will
 *         be called using FSAB algorithm
 */
void RCArticulatedBody::apply_impulse(const Vector3& j, const Vector3& k, const Vector3& contact_point, RigidBodyPtr link)
{
  if (!is_enabled())
    return;

  // see whether we can skip this
  if (j.norm() < NEAR_ZERO && k.norm() < NEAR_ZERO)
    return;

  // compute the forward dynamics, given the algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      _fsab.apply_impulse(j, k, contact_point, link);
      break;

    case eCRB:
      _crb.apply_impulse(j, k, contact_point, link);
      break;

    default:
      assert(false);
  }
}

/// Synchronizes the visualization for all links and joints of this body to their corresponding transforms
void RCArticulatedBody::update_visualization()
{
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->update_visualization();
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->update_visualization();
}

/// Gets the generalized coordinates of this body
VectorN& RCArticulatedBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gc)
{
  // resize gc
  gc.resize(num_generalized_coordinates(gctype));

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // get the joint positions of all implicit joints
  for (unsigned i=0; i < _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    gc.set_sub_vec(idx, _ijoints[i]->q);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return gc;

  // get the generalized coordinates for the base -- NOTE: we put the base
  // acceleration in the base frame
  assert(!_links.empty());
  RigidBodyPtr base = _links.front();
  if (gctype == DynamicBody::eAxisAngle)
  {
    gc.set_sub_vec(0, base->get_lvel());
    gc.set_sub_vec(3, base->get_avel());
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);
    VectorN base_gc;
    base->get_generalized_coordinates(gctype, base_gc);
    Vector3 base_x(base_gc[0], base_gc[1], base_gc[2]);
    base_x = base->get_transform().transpose_mult_vector(base_x);
    base_gc.set_sub_vec(0, base_x);
    gc.set_sub_vec(0, base_gc);
  }
    
  return gc;
}

/// Sets the generalized position of this body
void RCArticulatedBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gc)
{
  assert(num_generalized_coordinates(gctype) == gc.size());

  // see whether the body has a floating base  
  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    gc.get_sub_vec(0, gc.size()-_n_joint_DOF_implicit, _workv);
    if (gctype == DynamicBody::eRodrigues)
    {
      Vector3 base_x(_workv[0], _workv[1], _workv[2]);
      base_x = base->get_transform().mult_vector(base_x);
      _workv.set_sub_vec(0, base_x);
    }
    base->set_generalized_coordinates(gctype, _workv);
  }

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // set the generalized coordinates for the implicit joints
  for (unsigned i=0; i < _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    gc.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _ijoints[i]->q);
  }

  // link transforms must now be updated
  update_link_transforms();
}

/// Sets the generalized velocity of this body
void RCArticulatedBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gv)
{
  assert(num_generalized_coordinates(gctype) == gv.size());

  if (!is_enabled())
    return;

  // see whether the body has a floating base  
  if (_floating_base)
  {
    RigidBodyPtr base = _links.front();
    gv.get_sub_vec(0, gv.size()-_n_joint_DOF_implicit, _workv);
    if (gctype == DynamicBody::eRodrigues)
    {
      Vector3 base_xd(_workv[0], _workv[1], _workv[2]);
      base_xd = base->get_transform().mult_vector(base_xd);
      _workv.set_sub_vec(0, base_xd);
    }
    base->set_generalized_velocity(gctype, _workv);
  }

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // set the generalized velocities for the implicit joints
  for (unsigned i=0; i < _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    gv.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _ijoints[i]->qd);
  }

  // compute explicit constraint velocities
  if (!_ejoints.empty()) 
  {
    determine_explicit_constraint_movement_jacobian(_Dx);
    _Dx.mult(gv, _workv2);
    for (unsigned i=0; i< _ejoints.size(); i++)
    {
      unsigned idx = _ejoints[i]->get_coord_index();
      _workv2.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->qd);
    }
  }

  // link velocities must now be updated
  update_link_velocities();
}

/// Gets the generalized velocity of this body
VectorN& RCArticulatedBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv)
{
  // look for easy exit
  if (!is_enabled())
    return gv.resize(0);

  // resize gv
  gv.resize(num_generalized_coordinates(gctype));

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // get the joint velocities of all joints
  for (unsigned i=0; i < _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index() + add;
    gv.set_sub_vec(idx, _ijoints[i]->qd);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return gv;

  // get the generalized velocities for the base
  RigidBodyPtr base = _links.front();
  if (gctype == DynamicBody::eAxisAngle)
  {
    gv.set_sub_vec(0, base->get_lvel());
    gv.set_sub_vec(3, base->get_avel());
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);
    VectorN base_gv;
    base->get_generalized_velocity(gctype, base_gv);
    Vector3 base_xd(base_gv[0], base_gv[1], base_gv[2]);
    base_xd = base->get_transform().transpose_mult_vector(base_xd);
    base_gv.set_sub_vec(0, base_xd);
    gv.set_sub_vec(0, base_gv);
  }

  return gv;
}

/// Gets the generalized inertia of this body
MatrixN& RCArticulatedBody::get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M) 
{
  // look for easy exit
  if (!is_enabled())
    return M.resize(0,0);

  // calculate the generalized inertia matrix
  _crb.calc_generalized_inertia(gctype, M);

  return M;
}

/// Gets the number of degrees of freedom for explicit constraints
unsigned RCArticulatedBody::num_joint_dof_explicit() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    ndof += _ejoints[i]->num_dof();
  return ndof;
}

/// Determines contact Jacobians
void RCArticulatedBody::determine_contact_jacobians(const EventProblemData& q, 
const VectorN& v, const MatrixN& M, MatrixN& Jc, MatrixN& Dc)
{
/*  SAFESTATIC VectorN vnew, vcurrent;

  // see whether velocity data is valid
  bool vvalid = !velocities_invalidated();

  // get # of generalized coordinates (axis angle representation)
  const unsigned NGC = v.size();

  // resize dJc, dDc
  _iM_JcT.set_zero(NGC, q.N_CONTACTS);
  _iM_DcT.set_zero(NGC, q.N_CONTACTS*2);

  // copy generalized velocity
  vcurrent.copy_from(v);

  // loop over all contact events
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++)
  {
    // don't process if contact event is inactive
    if (!q.contact_working_set[i])
      continue;

    // get the articulated bodies of the contacts
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
    
    // see whether we can skip contact event
    if (ab1.get() != this && ab2.get() != this)
    {
      ii += 2;
      continue;
    }

    // get the rigid bodies corresponding to sb1 and sb2
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(sb1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(sb2);

    // get the contact normal and contact point
    const Vector3& p = q.contact_events[i]->contact_point;
    const Vector3& normal = q.contact_events[i]->contact_normal;

    // apply the impulse at the contact point
    if (ab1 == get_this())
      apply_impulse(normal, ZEROS_3, p, rb1);
    if (ab2 == get_this())
      apply_impulse(-normal, ZEROS_3, p, rb2);

    // measure the change in velocity
    get_generalized_velocity(DynamicBody::eAxisAngle, vnew);
    vcurrent.negate() += vnew;
    _iM_JcT.set_column(i, vcurrent);
    vcurrent.copy_from(vnew);

    // do the first tangent direction
    if (ab1 == get_this())
      apply_impulse(q.contact_events[i]->contact_tan1, ZEROS_3, p, rb1);
    if (ab2 == get_this())
      apply_impulse(-q.contact_events[i]->contact_tan1, ZEROS_3, p, rb2);

    // measure the change in velocity
    get_generalized_velocity(DynamicBody::eAxisAngle, vnew);
    vcurrent.negate() += vnew;
    _iM_DcT.set_column(ii++, vcurrent);
    vcurrent.copy_from(vnew);     

    // do the second tangent direction
    if (ab1 == get_this())
      apply_impulse(q.contact_events[i]->contact_tan2, ZEROS_3, p, rb1);
    if (ab2 == get_this())
      apply_impulse(-q.contact_events[i]->contact_tan2, ZEROS_3, p, rb2);

    // measure the change in velocity
    get_generalized_velocity(DynamicBody::eAxisAngle, vnew);
    vcurrent.negate() += vnew;
    _iM_DcT.set_column(ii++, vcurrent);
    vcurrent.copy_from(vnew);     
  }

  // restore generalized velocity
  set_generalized_velocity(DynamicBody::eAxisAngle, v);

  // see whether to re-validate velocities
  if (vvalid)
    validate_velocities();

  // compute Jc and Dc
  _iM_JcT.transpose_mult(M, Jc);
  _iM_DcT.transpose_mult(M, Dc);
*/

  // resize Jc and Dc
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eAxisAngle);
  Jc.set_zero(q.N_CONTACTS, NGC);
  Dc.set_zero(q.N_CONTACTS*2, NGC);

  // determine the Jacobian over all contact events 
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+= 2)
  {
    // don't process if contact event is inactive
    if (!q.contact_working_set[i])
      continue;

    // get the articulated bodies of the contacts
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
    
    // see whether we can skip contact event
    if (ab1.get() != this && ab2.get() != this)
      continue;

    // get the rigid bodies corresponding to sb1 and sb2
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(sb1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(sb2);

    // get the contact point and orthonormal basis 
    const Vector3& p = q.contact_events[i]->contact_point;
    const Vector3& normal = q.contact_events[i]->contact_normal;
    const Vector3& tan1 = q.contact_events[i]->contact_tan1;
    const Vector3& tan2 = q.contact_events[i]->contact_tan2;

    // compute R
    _R.resize(3,3);
    _R.set_row(0, normal);
    _R.set_row(1, tan1);
    _R.set_row(2, tan2); 

    // compute the Jacobian at the contact point
    if (ab1 == get_this())
      calc_jacobian(p, rb1, _workM);
    else
    {
      calc_jacobian(p, rb2, _workM);
      _R.negate();
    }

    // transform the Jacobian
    _workM.get_sub_mat(0,3,0,_workM.columns(), _workM2);
    _R.mult(_workM2, _workM); 

    // set the appropriate blocks in the Jacobians
    _workM.get_sub_mat(0,1,0,_workM.columns(), _workM2);
    Jc.set_sub_mat(i, 0, _workM2);
    _workM.get_sub_mat(1,3,0,_workM.columns(), _workM2);
    Dc.set_sub_mat(ii, 0, _workM2);
  }

  // setup iM_JcT and iM_DcT
  solve_generalized_inertia_transpose(DynamicBody::eAxisAngle, Jc, _iM_JcT);
  solve_generalized_inertia_transpose(DynamicBody::eAxisAngle, Dc, _iM_DcT);

  FILE_LOG(LOG_EVENT) << "Jc:" << std::endl << Jc;
  FILE_LOG(LOG_EVENT) << "Dc:" << std::endl << Dc;
}

/// Updates the event data
void RCArticulatedBody::update_event_data(EventProblemData& q)
{
  const unsigned SPATIAL_DIM = 6;

  // get the generalized velocity (axis angle)
  get_generalized_velocity(DynamicBody::eAxisAngle, _workv);

  // get the generalized inertia matrix
  _crb.calc_generalized_inertia(DynamicBody::eAxisAngle, _workM);

  // determine contact normal and tangent Jacobians
  determine_contact_jacobians(q, _workv, _workM, _Jc, _Dc);

  // setup Jx (neqx x ngc) and Dx (nedof x ngc)
  determine_explicit_constraint_jacobians(q, _Jx, _Dx);

  // setup Jl (nl x ngc)
  const unsigned NGC = _workv.size();
  _Jl.set_zero(q.N_LIMITS, NGC);
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    // get the limit joint
    JointPtr joint = q.limit_events[i]->limit_joint;

    // see whether to skip this joint
    if (joint->get_articulated_body() != get_this())
      continue;

    // determine the degree-of-freedom
    unsigned dof = joint->get_coord_index() + q.limit_events[i]->limit_dof;

    // set the component of Jl appropriately
    _Jl(i, dof) = (q.limit_events[i]->limit_upper) ? (Real) -1.0 : (Real) 1.0;
  }

  // setup Dt (nt x ngc)
  _Dt.set_zero(q.N_CONSTRAINT_DOF_IMP, NGC);
  if (use_advanced_friction_model)
  {
    for (unsigned i=0, ii=0; i< q.constraint_events.size(); i++)
    {
      // get the joint
      JointPtr joint = q.constraint_events[i]->constraint_joint;

      // look for a fast skip
      if (joint->get_constraint_type() != Joint::eImplicit)
        continue;

      // determine whether the joint belongs to this body
      if (joint->get_articulated_body() != get_this())
      {
        ii += joint->num_dof();
        continue;
      }
 
      // get the starting index for this joint
      const unsigned ST_IDX = joint->get_coord_index();

      // setup entries of Dt
      for (unsigned j=0; j< joint->num_dof(); j++, ii++)
        _Dt(ii, ST_IDX+j) = (Real) 1.0;
    }
  }

  // setup some temporary matrices
  solve_generalized_inertia_transpose(eAxisAngle, _Jl, _iM_JlT);
  solve_generalized_inertia_transpose(eAxisAngle, _Jx, _iM_JxT);
  solve_generalized_inertia_transpose(eAxisAngle, _Dx, _iM_DxT);
  solve_generalized_inertia_transpose(eAxisAngle, _Dt, _iM_DtT);

  // update all matrices
  q.Jc_iM_JcT += _Jc.mult(_iM_JcT, _workM2);
  q.Jc_iM_DcT += _Jc.mult(_iM_DcT, _workM2);
  q.Jc_iM_JlT += _Jc.mult(_iM_JlT, _workM2);
  q.Jc_iM_DtT += _Jc.mult(_iM_DtT, _workM2);
  q.Jc_iM_JxT += _Jc.mult(_iM_JxT, _workM2);
  q.Jc_iM_DxT += _Jc.mult(_iM_DxT, _workM2);
  q.Dc_iM_DcT += _Dc.mult(_iM_DcT, _workM2);
  q.Dc_iM_JlT += _Dc.mult(_iM_JlT, _workM2);
  q.Dc_iM_DtT += _Dc.mult(_iM_DtT, _workM2);
  q.Dc_iM_JxT += _Dc.mult(_iM_JxT, _workM2);
  q.Dc_iM_DxT += _Dc.mult(_iM_DxT, _workM2);
  q.Jl_iM_JlT += _Jl.mult(_iM_JlT, _workM2);
  q.Jl_iM_DtT += _Jl.mult(_iM_DtT, _workM2);
  q.Jl_iM_JxT += _Jl.mult(_iM_JxT, _workM2);
  q.Jl_iM_DxT += _Jl.mult(_iM_DxT, _workM2);
  q.Dt_iM_DtT += _Dt.mult(_iM_DtT, _workM2);
  q.Dt_iM_JxT += _Dt.mult(_iM_JxT, _workM2);
  q.Dt_iM_DxT += _Dt.mult(_iM_DxT, _workM2);
  q.Jx_iM_JxT += _Jx.mult(_iM_JxT, _workM2);
  q.Jx_iM_DxT += _Jx.mult(_iM_DxT, _workM2);
  q.Dx_iM_DxT += _Dx.mult(_iM_DxT, _workM2);

  // update velocity vectors
  q.Jc_v += _Jc.mult(_workv, _workv2);
  q.Dc_v += _Dc.mult(_workv, _workv2);
  q.Jl_v += _Jl.mult(_workv, _workv2);
  q.Jx_v += _Jx.mult(_workv, _workv2);
  q.Dx_v += _Dx.mult(_workv, _workv2);
}

/// Updates the body velocities
void RCArticulatedBody::update_velocity(const EventProblemData& q)
{
  // get the current spatial velocity
  get_generalized_velocity(DynamicBody::eAxisAngle, _v);

  if (LOGGING(LOG_EVENT))
  {
    MatrixN M;
    FILE_LOG(LOG_EVENT) << "alpha_c: " << q.alpha_c << std::endl;
    FILE_LOG(LOG_EVENT) << "beta_c: " << q.beta_c << std::endl;
    FILE_LOG(LOG_EVENT) << "alpha_l: " << q.alpha_l << std::endl;
    FILE_LOG(LOG_EVENT) << "beta_t: " << q.beta_t << std::endl;
    FILE_LOG(LOG_EVENT) << "velocity before: " << _v << std::endl;
    FILE_LOG(LOG_EVENT) << "kinetic energy before: " << calc_kinetic_energy() << std::endl;
    FILE_LOG(LOG_EVENT) << "kinetic energy before (calculation 2): " << _v.dot(get_generalized_inertia(DynamicBody::eAxisAngle, M).mult(_v*0.5)) << std::endl;
  }

  // check whether we are using a subset of the contacts
  if (std::find(q.contact_working_set.begin(), q.contact_working_set.end(), false) != q.contact_working_set.end())
  {
    // setup normal and tangent contact indices
    _nrm_c_indices.clear();
    _tan_c_indices.clear();
    for (unsigned i=0, j=0; i< q.contact_working_set.size(); i++, j+=2)
      if (q.contact_working_set[i])
      {
        _nrm_c_indices.push_back(i);
        _tan_c_indices.push_back(j);
        _tan_c_indices.push_back(j+1);
      }

    // compute change in velocities
    _iM_JcT.select_columns(_nrm_c_indices.begin(), _nrm_c_indices.end(), _workM);
     _v += _workM.mult(q.alpha_c, _workv2);
    _iM_DcT.select_columns(_tan_c_indices.begin(), _tan_c_indices.end(), _workM);
     _v += _workM.mult(q.beta_c, _workv2);
  }
  else
  {
    // compute change in velocities
     _v += _iM_JcT.mult(q.alpha_c, _workv2);
     _v += _iM_DcT.mult(q.beta_c, _workv2);
  }
  _v += _iM_JlT.mult(q.alpha_l, _workv2);
  _v += _iM_JxT.mult(q.alpha_x, _workv2);
  _v += _iM_DxT.mult(q.beta_x, _workv2);
  _v += _iM_DtT.mult(q.beta_t, _workv2);

  // set the spatial velocity
  set_generalized_velocity(DynamicBody::eAxisAngle, _v);

  FILE_LOG(LOG_EVENT) << "velocity after: " << _v << std::endl;
  FILE_LOG(LOG_EVENT) << "kinetic energy after: " << calc_kinetic_energy() << std::endl;
}

/// Gets the generalized forces on this body
/**
 * \note does not add forces from explicit joint constraints!
 */
VectorN& RCArticulatedBody::get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, VectorN& f) 
{
  const unsigned SPATIAL_DIM = 6;

  // look for fast exit
  if (!is_enabled())
    return f.resize(0);

  // resize f
  f.resize(num_generalized_coordinates(gctype));

  // compute the generalized forces
  SVector6 f0;
  _crb.calc_generalized_forces(f0, _workv);

  // see whether to add one to joint coords
  unsigned add = (_floating_base && gctype == DynamicBody::eRodrigues) ? 1 : 0;

  // determine the vector of joint forces
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    unsigned idx = _ijoints[i]->get_coord_index();
    unsigned cidx = (_floating_base) ? idx - SPATIAL_DIM : idx;
    idx += add;
    for (unsigned m=0; m< _ijoints[i]->num_dof(); m++, idx++, cidx++)
      f[idx] = _ijoints[i]->force[m] - _workv[cidx]; 
  }

  // setup joint space part of f
  if (_floating_base)
  {
    // determine external and inertial forces on base 
    RigidBodyPtr base = _links.front();
    if (gctype == DynamicBody::eAxisAngle)
      f.set_sub_vec(0, -f0);
    else
    {
      assert(gctype == DynamicBody::eRodrigues);
      const Quat& q = base->get_orientation();
      Vector3 base_force = -f0.get_upper();
      Vector3 base_torque = -f0.get_lower();
      f.set_sub_vec(0, base_force);
      f.set_sub_vec(3, base_torque);
      f[6] = Quat::deriv(q, base->get_avel()).norm_sq();
    }
  }

  return f;
}

/// Converts a force to a generalized force
VectorN& RCArticulatedBody::convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& p, const Vector3& f, const Vector3& t, VectorN& gf)
{
  const unsigned SPATIAL_DIM = 6;

  // look for fast exit
  if (!is_enabled())
    return gf.resize(0);

  // get the body as a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(body);
  assert(rb);

  // make f/t a spatial vector in the frame
  SpatialTransform X = SpatialTransform(IDENTITY_3x3, p, IDENTITY_4x4);
  SVector6 sf(f, t);
  SVector6 ft = X.transform(sf);

  // determine the Jacobian in the global frame
  const unsigned BASE_DIM = (_floating_base) ? SPATIAL_DIM : 0;
  _workM.set_zero(SPATIAL_DIM, _n_joint_DOF_implicit);
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the gc index for the joint
    unsigned idx = _ijoints[i]->get_coord_index() - BASE_DIM;

    // only add in Jacobian contribution if joint i supports the link
    if (supports(_ijoints[i], rb))
      _workM.set_sub_mat(0, idx, _ijoints[i]->get_spatial_axes(eGlobal));
  }

  // get the torque on the joints
  if (!_floating_base)
  {
    _workM.transpose_mult(ft, gf);
    return gf;
  }
  else
    _workM.transpose_mult(ft, _workv);

  FILE_LOG(LOG_DYNAMICS) << "RCArticulatedBody::convert_to_generalized_force() - converting " << std::endl << "  " << f << " / " << t << " to generalized force" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  -- spatial force (global frame): " << ft << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  -- Jacobian: " << std::endl << _workM.transpose();
  FILE_LOG(LOG_DYNAMICS) << "  -- joint torques: " << _workv << std::endl;

  // resize the generalized force vector
  gf.resize(num_generalized_coordinates(gctype));
  gf.set_sub_vec(gf.size() - _n_joint_DOF_implicit, _workv);

  // determine the generalized force on the base
  RigidBodyPtr base = _links.front();
  if (gctype == DynamicBody::eAxisAngle)
  {
    gf.set_sub_vec(0, base->sum_forces());
    gf.set_sub_vec(3, base->sum_torques() - base->calc_inertial_forces());
  }
  else
  {
    Vector3 fb = base->get_transform().transpose_mult_vector(f);
    base->convert_to_generalized_force(gctype, base, p, fb, t, _workv);
    FILE_LOG(LOG_DYNAMICS) << "  -- generalized force on base: " << _workv << std::endl;
    gf.set_sub_vec(0, _workv);
  }

  // resize the generalized force vector
  return gf;
}

/// Determines whether a joint supports a link
bool RCArticulatedBody::supports(JointPtr joint, RigidBodyPtr link)
{
  // save the original link
  RigidBodyPtr l = link; 

  do
  {
    JointPtr j = l->get_inner_joint_implicit();
    // check for l is base
    if (!j)
      return false;
    if (j == joint)
      return true;

    // proceed up the chain
    l = j->get_inboard_link();
  }
  while (link != l);

  return false;
}

/// Implements Base::load_from_xml()
/**
 * \pre all links and joints must be loaded using their respective
 *      serialization methods before this method is called
 */
void RCArticulatedBody::load_from_xml(XMLTreeConstPtr node, map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // load the parent data
  ArticulatedBody::load_from_xml(node, id_map);

  // don't verify the node name -- this class has derived classes
//  assert(strcasecmp(node->name.c_str(), "RCArticulatedBody") == 0);

  // get whether the body has a floating base
  const XMLAttrib* fb_attr = node->get_attrib("floating-base");
  if (fb_attr)
    _floating_base = fb_attr->get_bool_value();

  // read the pointer to the forward dynamics algorithm, if provided
  const XMLAttrib* fdyn_algo_attr = node->get_attrib("fdyn-algorithm");
  if (fdyn_algo_attr)
  {
    // get the ID
    string algo = fdyn_algo_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = algo.find_first_not_of(" \t\n\r");
    size_t last_nws_index = algo.find_last_not_of(" \t\n\r");
    algo = algo.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the algorithm type
    if (strcasecmp(algo.c_str(), "fsab") == 0)
      algorithm_type = eFeatherstone;
    else if (strcasecmp(algo.c_str(), "crb") == 0)
      algorithm_type = eCRB;
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "forward dynamics algorithm " << std::endl << "  type '";
      std::cerr << algo << "' -- valid types are 'fsab' and 'crb'";
      std::cerr << std::endl;
    }
  }

  // read the forward dynamics algorithm computation frame, if provided
  const XMLAttrib* fdyn_frame_attr = node->get_attrib("fdyn-algorithm-frame");
  if (fdyn_frame_attr)
  {
    // get the computation reference frame as a string
    string frame = fdyn_frame_attr->get_string_value();

    // remove leading and trailing spaces from the string name
    size_t first_nws_index = frame.find_first_not_of(" \t\n\r");
    size_t last_nws_index = frame.find_last_not_of(" \t\n\r");
    frame = frame.substr(first_nws_index, last_nws_index-first_nws_index+1);

    // get the frame type
    if (strcasecmp(frame.c_str(), "global") == 0)
      computation_frame_type = eGlobal;
    else if (strcasecmp(frame.c_str(), "link") == 0)
      computation_frame_type = eLink;
    else
    {
      std::cerr << "RCArticulatedBody::load_from_xml() - unknown ";
      std::cerr << "computation reference " << std::endl << "  frame type '";
      std::cerr << frame << "' -- valid types are 'link' and 'global'";
      std::cerr << std::endl;
    }
  }

  // get baumgarte paramters
  const XMLAttrib* balpha_attr = node->get_attrib("baumgarte-alpha");
  if (balpha_attr)
    b_alpha = balpha_attr->get_real_value();
  const XMLAttrib* bbeta_attr = node->get_attrib("baumgarte-beta");
  if (bbeta_attr)
    b_beta = bbeta_attr->get_real_value();

  // compile everything once again, for safe measure
  compile();

  // transform the body, if desired
  const XMLAttrib* xlat_attr = node->get_attrib("translation");
  const XMLAttrib* transform_attr = node->get_attrib("transform");

  // verify both weren't used
  if (xlat_attr && transform_attr)
    std::cerr << "RCArticulatedBody::load_from_xml() warning- 'translation' and 'transform' attributes both specified; using neither" << std::endl;
  else if (xlat_attr)
  {
    RigidBodyPtr base = get_base_link();
    Vector3 x;
    xlat_attr->get_vector_value(x);
    Matrix4 T = base->get_transform();
    T.set_translation(x);
    base->set_transform(T);
    update_link_transforms();
  }
  else if (transform_attr)
  {
    Matrix4 T;
    transform_attr->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
    {
      std::cerr << "RCArticulatedBody::load_from_xml() warning: invalid transform ";
      std::cerr << std::endl << T << " when reading node " << std::endl;
      std::cerr << *node << std::endl;
      std::cerr << "  --> possibly a floating-point error..." << std::endl;
    }
    else
    {
      RigidBodyPtr base = get_base_link();
      base->set_transform(T * base->get_transform());
      update_link_transforms();
    }
  }
}

/// Implements Base::save_to_xml()
void RCArticulatedBody::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // call the parent method first
  ArticulatedBody::save_to_xml(node, shared_objects);  

  // (re)set the name of this node
  node->name = "RCArticulatedBody";

  // write whether body has a floating base
  node->attribs.insert(XMLAttrib("floating-base", _floating_base));

  // write the forward dynamics algorithm type
  if (algorithm_type == eFeatherstone)
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("fsab")));
  else
  {
    assert(algorithm_type == eCRB);
    node->attribs.insert(XMLAttrib("fdyn-algorithm", string("crb")));
  }

  // write the forward dynamics algorithm frame -- note that the string()
  // is necessary on the second argument to XMLAttrib b/c the compiler
  // interprets a constant string as a bool, rather than as an string,
  // given a choice
  if (computation_frame_type == eGlobal)
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("global")));
  else
  {
    assert(computation_frame_type == eLink);
    node->attribs.insert(XMLAttrib("fdyn-algorithm-frame", string("link")));
  }

  // save baumgarte parameters
  node->attribs.insert(XMLAttrib("baumgarte-alpha", b_alpha));
  node->attribs.insert(XMLAttrib("baumgarte-beta", b_beta));
}

