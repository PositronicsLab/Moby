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

