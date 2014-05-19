/*
template <class M>
M& RCArticulatedBody::solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const M& B, M& X)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eEuler || _invM_rankdef)
    _invM.mult(B, X);
  else
  {
    assert(algorithm_type == eCRB && gctype == eSpatial);
    _crb.M_solve(B, X);
  }

  return X;
}

template <class M>
M& RCArticulatedBody::transpose_solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const M& B, M& X)
{
  // update the inverse / factorized inertia (if necessary)
  update_factorized_generalized_inertia(gctype);

  if (algorithm_type == eFeatherstone || gctype == eEuler || _invM_rankdef)
    _invM.mult_transpose(B, X);
  else
  {
    assert(algorithm_type == eCRB && gctype == eSpatial);
    Ravelin::MatrixNd BT = Ravelin::MatrixNd::transpose(B);
    _crb.M_solve(BT, X);
  }

  return X;
}
*/

/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as 
 * well as the base momentum; therefore, the derivative of the state vector is 
 * composed of the joint-space accelerations and base forces (and torques).
 */
template <class V>
void RCArticulatedBody::get_generalized_acceleration_generic(V& ga)
{
  const unsigned GC_AA_DIM = 6, GC_ROD_DIM = 7;

  // resize the state-derivative vector
  ga.resize(num_generalized_coordinates(DynamicBody::eSpatial));
  
  // setup the generalized acceleration for the base (if any) 
  if (_floating_base)
  {
    Ravelin::SharedVectorNd base_ga = ga.segment(num_joint_dof_explicit(), ga.size());
    RigidBodyPtr base = _links.front();
    base->get_generalized_acceleration_generic(base_ga);
  }

  // setup the state for the joints
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    ga.set_sub_vec(idx, _ejoints[i]->qdd);
  }
}

/// Gets the generalized coordinates of this body
template <class V>
void RCArticulatedBody::get_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc)
{
  // resize gc
  gc.resize(num_generalized_coordinates(gctype));

  // get the joint positions of all explicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gc.set_sub_vec(idx, _ejoints[i]->q);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return;

  // get the generalized coordinates for the base -- NOTE: we put the base
  // acceleration in the base frame
  assert(!_links.empty());
  RigidBodyPtr base = _links.front();
  Ravelin::SharedVectorNd base_gc = gc.segment(num_joint_dof_explicit(), gc.size());
  base->get_generalized_coordinates_generic(gctype, base_gc);
}

/// Sets the generalized position of this body
template <class V>
void RCArticulatedBody::set_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc)
{
  assert(num_generalized_coordinates(gctype) == gc.size());

  // set the generalized coordinates for the implicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gc.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->q);
  }

  // update base gc, if necessary
  if (_floating_base)
  {
    assert(!_links.empty());
    RigidBodyPtr base = _links.front();
    Ravelin::SharedConstVectorNd base_gc = gc.segment(num_joint_dof_explicit(), gc.size());
    base->set_generalized_coordinates_generic(gctype, base_gc);
  }

  // link transforms must now be updated
  update_link_poses();
}

/// Sets the generalized velocity of this body
template <class V>
void RCArticulatedBody::set_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv)
{
  SAFESTATIC Ravelin::VectorNd Dx_qd;
  assert(num_generalized_coordinates(gctype) == gv.size());

  // set the generalized velocities for the explicit joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gv.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->qd);
  }

  // see whether the body has a floating base  
  if (_floating_base)
  {
    assert(!_links.empty());
    RigidBodyPtr base = _links.front();
    Ravelin::SharedConstVectorNd base_gv = gv.segment(num_joint_dof_explicit(), gv.size());
    base->set_generalized_velocity_generic(gctype, base_gv);
  }

  // compute implicit constraint velocities
  if (!_ijoints.empty()) 
  {
    determine_implicit_constraint_movement_jacobian(_Dx);
    _Dx.mult(gv, Dx_qd);
    for (unsigned i=0; i< _ijoints.size(); i++)
    {
      unsigned idx = _ijoints[i]->get_coord_index();
      Dx_qd.get_sub_vec(idx, idx+_ijoints[i]->num_dof(), _ijoints[i]->qd);
    }
  }

  // link velocities must now be updated
  update_link_velocities();
}

/// Gets the generalized velocity of this body
template <class V>
void RCArticulatedBody::get_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv)
{
  // resize gv
  gv.resize(num_generalized_coordinates(gctype));

  // get the joint velocities of all joints
  for (unsigned i=0; i < _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    gv.set_sub_vec(idx, _ejoints[i]->qd);
  }

  // see whether the body has a floating base
  if (!_floating_base)
    return;

  // get the generalized velocity for the base
  assert(!_links.empty());
  RigidBodyPtr base = _links.front();
  Ravelin::SharedVectorNd base_gv = gv.segment(num_joint_dof_explicit(), gv.size());
  base->get_generalized_velocity_generic(gctype, base_gv);
  gv.set_sub_vec(num_joint_dof_explicit(), base_gv);
}


