#include <sstream>
#include <Moby/RigidBody.h>
#include <Moby/Joint.h>
#include <Moby/XMLTree.h>
#include <Moby/Log.h>
#include <Moby/select>
#include <Moby/SingularException.h>
#include <Moby/NumericalException.h>
#include <Moby/MCArticulatedBody.h>

using std::vector;
using std::list;
using std::map;
using std::string;
using std::endl;
using boost::shared_ptr;
using boost::static_pointer_cast;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Sets Baumgarte stabilization factors to alpha = 0.8, beta = 0.9
MCArticulatedBody::MCArticulatedBody()
{
  b_alpha = (double) 0.8;
  b_beta = (double) 0.9;
}

/// Applies a generalized impulse to the articulated body
void MCArticulatedBody::apply_generalized_impulse(const VectorNd& gj)
{
  SAFESTATIC VectorNd j;
  j.resize(num_generalized_coordinates(DynamicBody::eSpatial));

  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    const unsigned ngc = _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
    gj.get_sub_vec(k, k+ngc, j);
    _links[i]->apply_generalized_impulse_single(j);
    k += ngc;
  }
}

/// Adds a generalized force to the articulated body
void MCArticulatedBody::add_generalized_force( const VectorNd& gf)
{
  SAFESTATIC VectorNd f;
  f.resize(num_generalized_coordinates(DynamicBody::eSpatial));

  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    const unsigned NGC = _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
    gf.get_sub_vec(k, k+NGC, f);
    _links[i]->add_generalized_force_single(f);
    k += NGC;
  }
}

/// Gets the generalized coordinates of the body 
VectorNd& MCArticulatedBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gc)
{
  // necessary temporary vector
  SAFESTATIC VectorNd gc_body;

  // concatenate all link gc's into a single vector
  const unsigned NGC = num_generalized_coordinates(gctype);
  gc.resize(NGC);

  // get the generalized coordinates
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    _links[i]->get_generalized_coordinates_single(gctype, gc_body);
    gc.set_sub_vec(k, gc_body);
    k += _links[i]->num_generalized_coordinates_single(gctype);
  }

  // evaluate all joints constraints
  if (LOGGING(LOG_DYNAMICS))
  {
    double C[7];
    for (unsigned i=0; i< _joints.size(); i++)
    {
      _joints[i]->evaluate_constraints(C);
      std::ostringstream oss;
      oss << "joint " << _joints[i]->id << " constraint evaluations: ";
      for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
        oss << C[j] << " ";
      FILE_LOG(LOG_DYNAMICS) << oss.str() << std::endl;
    }
  }

  return gc; 
}

/// Gets the generalized velocity for this body
VectorNd& MCArticulatedBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gv) 
{
  // necessary temporary vector
  SAFESTATIC VectorNd gv_body;

  // determine number of generalized coordinates 
  const unsigned NGC = num_generalized_coordinates(gctype);
  gv.resize(NGC);

  // get the generalized velocities 
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    _links[i]->get_generalized_velocity_single(gctype, gv_body);
    gv.set_sub_vec(k, gv_body);
    k += _links[i]->num_generalized_coordinates_single(gctype);
  }

  return gv; 
}

/// Gets the generalized velocity for this body
VectorNd& MCArticulatedBody::get_generalized_acceleration( VectorNd& ga) 
{
  // necessary temporary vector
  SAFESTATIC VectorNd ga_body;

  // determine number of generalized coordinates 
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  ga.resize(NGC);

  // get the generalized accelerations
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    _links[i]->get_generalized_acceleration_single(ga_body);
    ga.set_sub_vec(k, ga_body);
    k += _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
  }

  return ga; 
}

/// Compiles this body; compilation must occur before the body is used
void MCArticulatedBody::compile()
{
  // call parent method first
  ArticulatedBody::compile();

  // verify ordering correctness on the disabled link (if any)
  for (unsigned i=1; i< _links.size(); i++)
    if (!_links[i]->is_enabled())
      throw std::runtime_error("Only first link can be disabled in a maximal-coordinate articulated body!");

  // determine all (axis-angle) gc indices
  _gc_indices.clear();
  _gc_indices.push_back(0);
  for (unsigned i=1; i<= _links.size(); i++)
    _gc_indices.push_back(_gc_indices.back() + _links[i-1]->num_generalized_coordinates_single(DynamicBody::eSpatial));

  // setup implicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx=0; i< _joints.size(); i++)
  {
    _joints[i]->set_coord_index(cidx);
    _joints[i]->set_constraint_index(ridx);
    cidx += _joints[i]->num_dof();
    ridx += _joints[i]->num_constraint_eqns();
  }

  // mark all joints as implicit 
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->set_constraint_type(Joint::eImplicit);
}

/// Gets the individual body inertias
void MCArticulatedBody::determine_inertias()
{
  // setup the sparse inverse inertia matrix
  _iM.resize(_links.size());

  // get the generalized inertia
  for (unsigned i=0; i< _links.size(); i++)
  {
    // if the link is disabled, we do something special
    if (!_links[i]->is_enabled())
    {
      _iM[i].inv_mass = (double) 0.0;
      _iM[i].inv_inertia .set_zero()

x3;
    }
    else
    {
      _iM[i].inv_mass = _links[i]->get_inv_mass();
      Matrix3 R(&_links[i]->get_orientation());
      _iM[i].inv_inertia = R * _links[i]->get_inv_inertia() * Matrix3::transpose(R);
    }
  }
}

/// Gets the generalized inertia matrix
MatrixNd& MCArticulatedBody::get_generalized_inertia( MatrixNd& M)
{
  // setup the sparse inertia matrix and inverse inertia matrix
  MatrixNd gI;

  // get the generalized inertia
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  M.set_zero(NGC, NGC);
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    // set the proper component of the generalized inertia matrix
    _links[i]->get_generalized_inertia_single(gI);
    M.set_sub_mat(k, k, gI);
    k += _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
  }

  return M;
}

/// Performs necessary pre-computations for computing accelerations for applying impulses
void MCArticulatedBody::precalc()
{
  SAFESTATIC VectorNd jv;

  // determine joint positions
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->determine_q(_joints[i]->q);

  // get the constraint Jacobian and the constraint force transform
  get_constraint_jacobian(_Jx);
  get_constraint_jacobian_dot(_Jx_dot);

  // get the mechanism Jacobian
  get_mechanism_jacobian(_Dx, _Dx_dot);

  // determine the sparse inertias 
  determine_inertias();

  // now form Jx_iM_Jx' and its inverse
  _rank_def = false;
  calc_Jx_iM_JyT(_Jx, _Jx, _Jx_iM_JxT);
  _inv_Jx_iM_JxT = _Jx_iM_JxT;
  if (!_LA.factor_chol(_inv_Jx_iM_JxT))
  {  
    _inv_Jx_iM_JxT = _Jx_iM_JxT;
    try
    {
      LinAlg::pseudo_inverse(_inv_Jx_iM_JxT, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      _inv_Jx_iM_JxT = _Jx_iM_JxT;
      LinAlg::pseudo_inverse(_inv_Jx_iM_JxT, LinAlg::svd2);
    }
    _rank_def = true;
  }

  // determine the generalized velocity; note that this must be
  // done after determining Jacobians, which are position dependent
  // get the generalized velocities
  get_generalized_velocities(DynamicBody::eSpatial, _xd);

  // setup joint velocities
  const unsigned NDOF = _Dx.rows();
  mult_sparse(_Dx, _xd, jv);
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned idx = _joints[i]->get_coord_index();
    jv.get_sub_vec(idx, idx+_joints[i]->num_dof(), _joints[i]->qd);
  }
}

/// Solves J*iM*J'*x = rhs for x
VectorNd& MCArticulatedBody::solve_Jx_iM_JxT(const VectorNd& rhs, VectorNd& x) const
{
  if (!_rank_def)
  {
    x = rhs;
    _LA.solve_chol_fast(_inv_Jx_iM_JxT, x);
  }
  else
    _inv_Jx_iM_JxT.mult(rhs, x);

  return x;
}

/// Gets the velocity state-derivative of this articulated body
void MCArticulatedBody::integrate(double t, double h, boost::shared_ptr<Integrator<VectorNd> > integrator)
{
return DynamicBody::integrate(t, h, integrator);

  // do preliminary calculations
  precalc();

  // setup forces
  VectorNd fext;
  get_generalized_forces(DynamicBody::eSpatial, fext);

  // convert forces to impulses
  fext *= h;

  // setup generalized coordinates (for semi-implicit integration) 
  VectorNd q;
  get_generalized_coordinates(DynamicBody::eEuler, q);

  // evaluate constraints
  VectorNd C;
  get_constraint_evals(C);

  FILE_LOG(LOG_DYNAMICS) << "J: " << std::endl << dense_J(_Jx);
  FILE_LOG(LOG_DYNAMICS) << "fext: " << fext << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "C: " << C << std::endl;

  // we use the Schur complement method to solve for the unknown accelerations
  // and generalized constraint forces (see [Nocedal 2006, p. 455])
  // first we solve J*inv(M)*K*lambda = J*inv(M)*(fext*h) 

  // determine inv(M)*(fext*h)
  VectorNd iM_fext;
  iM_mult(fext, iM_fext);
  FILE_LOG(LOG_DYNAMICS) << "inv(M)*fext: " << iM_fext << std::endl;

  // solve Jx*inv(M)*Jx' * lambda = Jx*inv(M)*(fext*h) - b
  // where b is 2*b_alpha*Jx*xd + b_beta^2*C (Baumgarte stabilization terms)
  VectorNd rhs, xd, lambda;
  (xd = _xd) *= (2*b_alpha);
  iM_fext += xd; 
  mult_sparse(_Jx, iM_fext, rhs);
  C *= (b_beta*b_beta);
  rhs += C;
  FILE_LOG(LOG_DYNAMICS) << "rhs for solve: " << rhs << endl;
  solve_Jx_iM_JxT(rhs, lambda);

  // ok, now that we've solved for lambda, solving for xd_delta is trivial
  VectorNd xd_delta;
  mult_transpose_sparse(_Jx, lambda, rhs) -= fext;
  rhs.negate();
  iM_mult(rhs, xd_delta);
  FILE_LOG(LOG_DYNAMICS) << "Jx'*lambda: " << mult_transpose_sparse(_Jx, lambda, rhs) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "xd: " << _xd << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "xd_delta: " << xd_delta << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "|| Jx*xd_delta ||: " << mult_sparse(_Jx, xd_delta, rhs).norm() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "|| Jx*xd ||: " << mult_sparse(_Jx, _xd, rhs).norm() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "lambda: " << lambda << std::endl;

  // update generalized coordinates
  for (unsigned i=0, j=0, k=0; i< _links.size(); i++)
  {
    if (!_links[i]->is_enabled())
      continue;
    Vector3d xd(_xd[j+0], _xd[j+1], _xd[j+2]);
    Vector3d omega(_xd[j+3], _xd[j+4], _xd[j+5]);
    Quat qd = _links[i]->get_orientation().G_transpose_mult(omega) * 0.5;
    q[k+0] += xd[0]*h; 
    q[k+1] += xd[1]*h;
    q[k+2] += xd[2]*h;
    q[k+3] += qd.w*h;
    q[k+4] += qd.x*h;
    q[k+5] += qd.y*h;
    q[k+6] += qd.z*h;
    j += 6;
    k += 7;
  }

  // update velocities
  _xd += xd_delta;

  // update the generalized velocities
  FILE_LOG(LOG_DYNAMICS) << "new generalized coordinates (Rodrigues parameters): " << q << endl;
  FILE_LOG(LOG_DYNAMICS) << "new generalized velocities (axis-angle parameters): " << _xd << endl;

  // set the generalized coordinates / velocities
  set_generalized_coordinates(DynamicBody::eEuler, q);
  update_link_velocities();
}

/// Computes the velocity state-derivative of this articulated body
void MCArticulatedBody::calc_fwd_dyn(double dt)
{
  // work variables
  SAFESTATIC VectorNd ff, fext, C, tmpv, x, alpha_x, beta_x, delta, iM_fext, z;
  SAFESTATIC VectorNd Jx_dot_xd, Jx_xd, vddt_plus_iM_fext, lambda;
  SAFESTATIC MatrixNd RG, Jx_iM_DxT, A;
  SAFESTATIC MatrixNd Dx_iM_DxT;

  // do pre-calculations first
  precalc();

  // evaluate constraints
  get_constraint_evals(C);

  FILE_LOG(LOG_DYNAMICS) << "MCArticulatedBody::calc_fwd_dyn() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << " constraint evals: " << C << endl;

  // compute Jx * xd and \dot{Jx} * xd
  mult_sparse(_Jx_dot, _xd, Jx_dot_xd);
  mult_sparse(_Jx, _xd, Jx_xd);

  // compute actuator forces
  const unsigned NDOF = _Dx.rows();
  ff.resize(NDOF);
  for (unsigned i=0, k=0; i< _joints.size(); i++)
    for (unsigned j=0; j< _joints[i]->num_dof(); j++)
      ff[k++] = _joints[i]->force[j];

  // setup spatial forces
  get_generalized_forces(DynamicBody::eSpatial, fext);

  // transform joint forces to absolute coords 
  mult_transpose_sparse(_Dx, ff, tmpv);
  fext += tmpv;

  // see whether to use advanced friction model
  if (!use_advanced_friction_model)
  {
    // use old Coulomb model
    for (unsigned i=0; i< _joints.size(); i++)
    {
      ff.resize(_joints[i]->num_dof());
      for (unsigned j=0; j< _joints[i]->num_dof(); j++)
      {
        ff[j] = (double) -1.0;
        if (_joints[i]->qd[j] < (double) 0.0)
          ff[j] = (double) 1.0;
        ff[j] *= _joints[i]->mu_fc;
        ff[j] -= _joints[i]->qd[j]*_joints[i]->mu_fv;
      }
      _joints[i]->add_force(ff);
    }

    // compute constraint forces (will be stored in 'z')
    iM_mult(fext, iM_fext);
    mult_sparse(_Jx, iM_fext, z) += Jx_dot_xd;
    (tmpv = Jx_xd) *= ((double) 2.0 * b_alpha);
    z += tmpv;
    C *= (b_beta*b_beta);
    z += C;
    z.negate(); 
    solve_Jx_iM_JxT(z, lambda);
    FILE_LOG(LOG_DYNAMICS) << " constraint forces: " << lambda << endl;

    // compute accelerations
    mult_transpose_sparse(_Jx, lambda, tmpv);
    fext += tmpv;
    iM_mult(fext, _xdd);   
    FILE_LOG(LOG_DYNAMICS) << " generalized accelerations: " << _xdd << endl;
    update_link_accelerations();

    return;
  }

  // if we're still here, then we're using the full joint friction model
  // still here? full-on joint friction model... 

  // get the number of generalized coordinates
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // determine how many constraint equations
  unsigned N_EXPLICIT_CONSTRAINT_EQNS = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    N_EXPLICIT_CONSTRAINT_EQNS = _joints[i]->num_constraint_eqns();

  // setup some constants 
  const unsigned N_JOINT_DOF = num_joint_dof();
  const unsigned ALPHAX_START = 0;
  const unsigned BETA_START = N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned N_EXPLICIT_DOF = N_JOINT_DOF;

  // compute Jx*inv(M)*Dx', etc.
  calc_Jx_iM_JyT(_Jx, _Dx, Jx_iM_DxT);
  calc_Jx_iM_JyT(_Dx, _Dx, Dx_iM_DxT);

  // setup data needed for convex optimization
  ABFwdDynOptData copt_data;
  copt_data.N_IMPLICIT_DOF = 0;
  copt_data.N_EXPLICIT_CONSTRAINT_EQNS = N_EXPLICIT_CONSTRAINT_EQNS;
  copt_data.N_JOINT_DOF = N_JOINT_DOF;
  dense_J(_Dx, copt_data.Dx);
  copt_data.fext = fext;
  vector<unsigned>& loop_indices = copt_data.loop_indices;
  vector<vector<unsigned> > loop_links;
  ArticulatedBody::find_loops(loop_indices, loop_links); 
  ArticulatedBody::compute_Z_matrices(loop_indices, loop_links, copt_data.Zd, copt_data.Z1d, copt_data.Z);
  unsigned N_LOOPS = 0;
  for (unsigned i=0; i< loop_indices.size(); i++)
    if (loop_indices[i] < std::numeric_limits<unsigned>::max())
      N_LOOPS = std::max(N_LOOPS, loop_indices[i]); 
  copt_data.N_LOOPS = N_LOOPS;

  // first compute the homogeneous solution
  const unsigned NN = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS + N_LOOPS;
  MatrixNd& R = copt_data.R;
  VectorNd& zz = copt_data.z;
  iM_mult(fext, iM_fext);
  mult_sparse(_Jx, iM_fext, zz) += Jx_dot_xd;
  (tmpv = Jx_xd)  *= ((double) 2.0 * b_alpha);
  zz += tmpv;
  C *= (b_beta*b_beta);
  zz += C;
  zz.negate(); 
  A.set_zero(N_EXPLICIT_CONSTRAINT_EQNS, NN);
  A.set_sub_mat(0, 0, _Jx_iM_JxT);
  A.set_sub_mat(0, N_EXPLICIT_CONSTRAINT_EQNS, Jx_iM_DxT);

  // verify there *is* a nullspace 
  if (N_EXPLICIT_CONSTRAINT_EQNS > 0)
  {
    // compute nullspace
    _LA.nullspace(A, R);

    // this makes frictional forces and deltas zero
    A = _Jx_iM_JxT;
    tmpv = zz;
    _LA.solve_LS_fast(A, tmpv);

    zz.set_zero(NN);
    zz.set_sub_vec(ALPHAX_START, tmpv);
  }
  else
  {
    // hack together a nullspace
    zz.set_zero(N_JOINT_DOF);
    R.set_zero(N_JOINT_DOF, N_JOINT_DOF);
    for (unsigned i=0; i< N_JOINT_DOF; i++)
      R(i,i) = (double) 1.0;
  }

  // setup components of z and R to make things faster for gradient and Hessian
  // calculations
  copt_data.zff.resize(0);
  zz.get_sub_vec(BETA_START, BETA_START+N_EXPLICIT_DOF, copt_data.zbetax);
  copt_data.Rff.resize(0, R.columns());
  R.get_sub_mat(BETA_START, BETA_START+N_EXPLICIT_DOF, 0, R.columns(), A);
  mult_transpose_sparse(_Dx, A, copt_data.DxTRbetax);

  // setup the data for the optimization problem
  const unsigned N = R.columns(); 
  const unsigned M = N_JOINT_DOF + N_LOOPS*2;
  OptParams copt_params(N, M, 0, ArticulatedBody::calc_fwd_dyn_f0, ArticulatedBody::calc_fwd_dyn_fx, NULL, ArticulatedBody::calc_fwd_dyn_grad0, ArticulatedBody::calc_fwd_dyn_cJac, NULL, ArticulatedBody::calc_fwd_dyn_hess);
  copt_params.max_iterations = N*N*N;
  copt_params.data = (void*) &copt_data;

  // verify that joint indices are as we expect
  #ifndef NDEBUG
  for (unsigned i=0; i< _joints.size(); i++)
    assert(_joints[i]->get_index() == i);
  #endif

  // setup true indices: converts from an [explicit implicit] DOF to the
  // joint's index that contains that DOF
  vector<unsigned>& true_indices = copt_data.true_indices;
  true_indices.resize(N_JOINT_DOF);
  for (unsigned i=0, ke=0; i< _joints.size(); i++)
  {
    assert(_joints[i]->get_constraint_type() == Joint::eImplicit);
    for (unsigned j=0; j< _joints[i]->num_dof(); j++)
      true_indices[ke++] = i;
  }

  // initialize mu_c and viscous force vectors
  vector<double>& mu_c = copt_data.mu_c;
  vector<double>& visc = copt_data.visc;
  mu_c.resize(N_JOINT_DOF);
  visc.resize(N_JOINT_DOF);

  // setup mu_c and viscous force vectors
  for (unsigned i=0, k=0; i< _joints.size(); i++)
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      mu_c[k] = _joints[i]->mu_fc;
      double tmp = _joints[i]->mu_fv * std::fabs(_joints[i]->qd[j]);
      visc[k] = tmp*tmp;
    }

  // compute the quadratic matrix
  MatrixNd& G = copt_data.G;
  G.set_zero(NN, NN);
  G.set_sub_mat(0,0,_Jx_iM_JxT);
  G.set_sub_mat(0,N_EXPLICIT_CONSTRAINT_EQNS,Jx_iM_DxT);
  G.set_sub_mat(N_EXPLICIT_CONSTRAINT_EQNS,0,Jx_iM_DxT, true);
  G.set_sub_mat(N_EXPLICIT_CONSTRAINT_EQNS,N_EXPLICIT_CONSTRAINT_EQNS,Dx_iM_DxT);
  R.transpose_mult(G, RG);
  RG.mult(R, G);

  // compute the linear optimization vector
  VectorNd& c = copt_data.c;
  c.set_zero(NN);
  (vddt_plus_iM_fext = _xd) /= dt;
  vddt_plus_iM_fext += iM_fext;
  mult_sparse(_Jx, vddt_plus_iM_fext, tmpv);
  c.set_sub_vec(0, tmpv);
  mult_sparse(_Dx, vddt_plus_iM_fext, tmpv);
  c.set_sub_vec(N_EXPLICIT_CONSTRAINT_EQNS, tmpv);
  R.transpose_mult(c, tmpv);
  c = tmpv;
  RG.mult(z, tmpv);
  c += tmpv;
  
  // setup a feasible optimization vector
  x.set_zero(N);  

  // optimize
  Optimization::optimize_convex_pd(copt_params, x);

  // (z + Ry)G(Ry + z) + y'R'c = 0.5 y'R'G'Ry + (R'Gz + R'c)
  R.mult(x, tmpv);
  z += tmpv;

  // now retrieve the necessary variables
  z.get_sub_vec(0, N_EXPLICIT_CONSTRAINT_EQNS, alpha_x);
  z.get_sub_vec(N_EXPLICIT_CONSTRAINT_EQNS, N_EXPLICIT_CONSTRAINT_EQNS+N_JOINT_DOF, beta_x);
  z.get_sub_vec(N_EXPLICIT_CONSTRAINT_EQNS+N_JOINT_DOF, z.size(), delta);

  // output results
  FILE_LOG(LOG_DYNAMICS) << " external forces: " << fext << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " implicit constraint forces: " << alpha_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " implicit friction forces: " << beta_x << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " delta: " << delta << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " constraint evaluations: " << C << std::endl;

  // setup joint friction forces
  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    beta_x.get_sub_vec(k, k+_joints[i]->num_dof(), _joints[i]->ff);
    k += _joints[i]->num_dof();
  }

  // compute joint constraint forces
  fext += mult_transpose_sparse(_Dx, beta_x, tmpv);
  ArticulatedBody::calc_joint_constraint_forces(loop_indices, delta, copt_data.Zd, copt_data.Z1d, copt_data.Z, fext);

  // compute generalized acceleration and link accelerations
  fext += mult_transpose_sparse(_Jx, alpha_x, tmpv);
  iM_mult(fext, _xdd);
  update_link_accelerations();
}

/// Computes the matrix Dx * inv(M) for use with computing sticking joint friction forces
void MCArticulatedBody::calc_Dx_iM(SparseJacobian& Dx_iM) const
{
  const unsigned SPATIAL_DIM = 6;

  // resize Dx
  Dx_iM.set_zero(_Dx.rows(), _Dx.columns());
  Dx_iM.indices = _Dx.indices;

  // loop over all links 
  for (unsigned k=0; k< _links.size(); k++)
  {
    // skip disabled bodies
    if (!_links[k]->is_enabled())
      continue;

    // get the inverse mass and inertia (global frame)
    double inv_mass = _links[k]->get_inv_mass();
    Matrix3 R(&_links[k]->get_orientation());
    Matrix3 invJ = R * _links[k]->get_inv_inertia() * Matrix3::transpose(R);

    // loop over all rows of Dx
    for (unsigned rowDx=0; rowDx < _Dx.rows(); rowDx++)
    {
      for (unsigned m=0, r=0; m< _Dx.indices[rowDx].size(); m++, r+= SPATIAL_DIM)
      {
        // get the rigid body responsible for this block of Dx
        if (_Dx.indices[rowDx][m] != k)
          continue;

        // get the linear and angular vectors
        Vector3d lv(_Dx(rowDx, r+0), _Dx(rowDx, r+1), _Dx(rowDx, r+2));
        Vector3d av(_Dx(rowDx, r+3), _Dx(rowDx, r+4), _Dx(rowDx, r+5));

        // scale the spatial vectors using the inverse inertias
        lv *= inv_mass;
        av = invJ * av;

        // place the linear and angular vectors
        Dx_iM(rowDx,r+0) = lv[0]; 
        Dx_iM(rowDx,r+1) = lv[1]; 
        Dx_iM(rowDx,r+2) = lv[2];
        Dx_iM(rowDx,r+3) = av[0]; 
        Dx_iM(rowDx,r+4) = av[1]; 
        Dx_iM(rowDx,r+5) = av[2];
      }
    }
  }
}

/// Computes the matrix Dx * inverse(M) * Dx' for use with computing sticking joint friction forces
void MCArticulatedBody::calc_Dx_iM_DxT(MatrixNd& Dx_iM_DxT) const
{
  const unsigned SPATIAL_DIM = 6;

  // resize Dx_iM_DxT
  Dx_iM_DxT.set_zero(_Dx.rows(), _Dx.rows());

  // loop over all links 
  for (unsigned k=0; k< _links.size(); k++)
  {
    // skip disabled bodies
    if (!_links[k]->is_enabled())
      continue;

    // get the inverse mass and inertia (global frame)
    double inv_mass = _links[k]->get_inv_mass();
    Matrix3 R(&_links[k]->get_orientation());
    Matrix3 invJ = R * _links[k]->get_inv_inertia() * Matrix3::transpose(R);

    // loop over all rows of Dx (Dx) 
    for (unsigned rowDx=0; rowDx < _Dx.rows(); rowDx++)
    {
      // loop over all blocks in this row of Dx
      for (unsigned m=0, r=0; m< _Dx.indices[rowDx].size(); m++, r+= SPATIAL_DIM)
      {
        // get the rigid body responsible for this block of Dx
        if (_Dx.indices[rowDx][m] != k)
          continue;

        // get the velocity vectors
        Vector3d l1(_Dx(rowDx, r+0), _Dx(rowDx, r+1), _Dx(rowDx, r+2));
        Vector3d a1(_Dx(rowDx, r+3), _Dx(rowDx, r+4), _Dx(rowDx, r+5));

        // scale the spatial vectors using the inverse inertias
        l1 *= inv_mass;
        a1 = invJ * a1;

        // loop over all rows of Dy (Dx) 
        for (unsigned rowDy=0; rowDy < _Dx.rows(); rowDy++)
        {
          // loop over all blocks in this row of Dy
          for (unsigned n=0, s=0; n< _Dx.indices[rowDy].size(); n++, s+= SPATIAL_DIM)
          {
            // get the rigid body responsible for this row of Dy'
            if (_Dx.indices[rowDy][n] != k)
              continue;

            // get the velocity vectors
            Vector3d l2(_Dx(rowDy, s+0), _Dx(rowDy, s+1), _Dx(rowDy, s+2));
            Vector3d a2(_Dx(rowDy, s+3), _Dx(rowDy, s+4), _Dx(rowDy, s+5));

            // compute the dot product
            Dx_iM_DxT(rowDx, rowDy) += l1.dot(l2) + a1.dot(a2);
          }
        }
      }
    }
  }
}

/// Sets link velocities from member variable _xd
void MCArticulatedBody::update_link_velocities() const
{
  const unsigned SPATIAL_DIM = 6;

  for (unsigned i=0, idx=0; i< _links.size(); i++)
  {
    // skip disabled links
    if (!_links[i]->is_enabled())
      continue;

    // get the spatial accelerations
    Vector3d xd(_xd[idx], _xd[idx+1], _xd[idx+2]);
    Vector3d omega(_xd[idx+3], _xd[idx+4], _xd[idx+5]);

    // set the link spatial velocities 
    _links[i]->set_lvel(xd);
    _links[i]->set_avel(omega);

    // update idx
    idx += SPATIAL_DIM;
  }
}

/// Sets link accelerations from member variable _xdd
void MCArticulatedBody::update_link_accelerations() const
{
  const unsigned SPATIAL_DIM = 6;

  for (unsigned i=0, idx=0; i< _links.size(); i++)
  {
    // skip disabled links
    if (!_links[i]->is_enabled())
      continue;

    // get the spatial accelerations
    Vector3d xdd(_xdd[idx], _xdd[idx+1], _xdd[idx+2]);
    Vector3d alpha(_xdd[idx+3], _xdd[idx+4], _xdd[idx+5]);

    // set the link spatial accelerations 
    _links[i]->set_laccel(xdd);
    _links[i]->set_aaccel(alpha);

    // update idx
    idx += SPATIAL_DIM;
  }
}

/// Calculates joint accelerations
/**
 * Joint acceleration is equal to Dx*\ddot{x} + \dot{Dx}*\dot{x},
 * where Dx is the mechanism jacobian, x is the generalized coordinates of the
 * body.
 */
void MCArticulatedBody::calc_joint_accelerations()
{
  // calculate the joint accelerations
  SAFESTATIC VectorNd qdd, tmp;
  mult_sparse(_Dx_dot, _xd, tmp); 
  mult_sparse(_Dx, _xdd, qdd) += tmp;

  // set the joint accelerations
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned idx = _joints[i]->get_coord_index();
    qdd.get_sub_vec(idx, idx+_joints[i]->num_dof(), _joints[i]->qdd);
  }
}

/// Applies an impulse to the articulated body
void MCArticulatedBody::apply_impulse(const SForced& w, RigidBodyPtr link)
{
  SAFESTATIC VectorNd iM_Qj, rhs, delta_xd, Qj, gj, lambda;
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  FILE_LOG(LOG_DYNAMICS) << "MCArticulatedBody::apply_impulse() entered" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " applying impulses: " << w << " to " << link->id << std::endl;

  // do preliminary calculations
  precalc();

  // ensure that link belongs to this body
  #ifndef NDEBUG
  if (std::find(_links.begin(), _links.end(), link) == _links.end())
    throw std::runtime_error("Link does not belong to this body!");
  #endif

  // setup generalized impulses 
  link->convert_to_generalized_force_single(DynamicBody::eSpatial, link, w, gj);
  Qj.set_zero(NGC);
  unsigned gc_idx_start = _gc_indices[link->get_index()];
  Qj.set_sub_vec(gc_idx_start, gj); 

  // we use the Schur complement method to solve for the unknown accelerations
  // and generalized constraint forces (see [Nocedal 2006, p. 455])

  // solve M*x = Qj
  iM_mult(Qj, iM_Qj);

  // solve J*inv(M)*K * lambda = J*inv(M)*Qj
  mult_sparse(_Jx, iM_Qj, rhs);
  solve_Jx_iM_JxT(rhs, lambda);

  // ok, now that we've solved for lambda, solving for delta_xd is trivial
  mult_transpose_sparse(_Jx, lambda, rhs) -= Qj;
  rhs.negate();
  iM_mult(rhs, delta_xd);

  FILE_LOG(LOG_DYNAMICS) << " generalized force: " << gj << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " total generalized forces: " << Qj << std::endl;

  // update the spatial velocities
  FILE_LOG(LOG_DYNAMICS) << " constraint impulses: " << lambda << std::endl;
  FILE_LOG(LOG_DYNAMICS) << " old generalized velocity (axis-angle): " << _xd << std::endl;
  _xd += delta_xd;  
  FILE_LOG(LOG_DYNAMICS) << " new generalized velocity (axis-angle): " << _xd << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "MCArticulatedBody::apply_impulse() exiting" << std::endl;

  // set the link velocities
  update_link_velocities();
}

/// Sets up the generalized force vector
VectorNd& MCArticulatedBody::get_generalized_forces( VectorNd& gf) 
{
  SAFESTATIC VectorNd f;

  // resize the gf vector
  gf.resize(num_generalized_coordinates(DynamicBody::eSpatial)); 

  // get the generalized forces 
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    _links[i]->get_generalized_forces_single(f);
    gf.set_sub_vec(k, f);
    k += _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
  }

  return gf;
}

/// Sets up the vector of constraint equation evaluations
VectorNd& MCArticulatedBody::get_constraint_evals(VectorNd& C) const
{
  const unsigned GC_ROD_DIM = 7;
  double C_array[GC_ROD_DIM];

  // resize C vector
  unsigned nc = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    nc += _joints[i]->num_constraint_eqns();
  C.resize(nc);

  // process all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    _joints[i]->evaluate_constraints(C_array);
    unsigned idx = _joints[i]->get_constraint_index();
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
      C[idx+j] = C_array[j];
  }

  return C;
}

/// Gets the number of generalized coordinates in the articulated body
unsigned MCArticulatedBody::num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const
{
  unsigned ngc = 0;
  for (unsigned i=0; i< _links.size(); i++)
    ngc += _links[i]->num_generalized_coordinates_single(gctype);

  return ngc;
}

/// Gets a submatrix of the given sparse Jacobian
void MCArticulatedBody::get_sub_jacobian(const vector<unsigned>& rows, const MCArticulatedBody::SparseJacobian& J, MCArticulatedBody::SparseJacobian& Jx)
{
  SAFESTATIC VectorNd row;

  // first, setup the indices
  Jx.indices.clear();
  for (unsigned i=0; i< rows.size(); i++)
    Jx.indices.push_back(J.indices[rows[i]]);

  // resize Jx
  Jx.resize(rows.size(), J.columns());

  // now, get the specified rows
  for (unsigned i=0, j=0; i< rows.size(); i++)
  {
    J.get_row(rows[i], row);
    Jx.set_row(j++, row);
  }
}

/// Solves using the transpose of the generalized inertia matrix
MatrixNd& MCArticulatedBody::transpose_solve_generalized_inertia(const MatrixNd& B, MatrixNd& X)
{
  X.resize(num_generalized_coordinates(DynamicBody::eSpatial), B.columns()); 
  VectorNd workv;

  // do the multiplication
  for (unsigned j=0; j< B.rows(); j++)
    for (unsigned i=0; i< _iM.size(); i++)
    {
      B.get_sub_mat(j, j+1, _gc_indices[i], _gc_indices[i+1], workv);
      scale_inverse_inertia(i, workv);
      X.set_sub_mat(_gc_indices[i], j, workv);
    }

  return X;
}

/// Solves using the generalized inertia matrix
MatrixNd& MCArticulatedBody::solve_generalized_inertia(const MatrixNd& B, MatrixNd& X)
{
  X.resize(num_generalized_coordinates(DynamicBody::eSpatial), B.columns()); 
  VectorNd workv;

  // do the multiplication
  for (unsigned j=0; j< B.columns(); j++)
    for (unsigned i=0; i< _iM.size(); i++)
    {
      B.get_sub_mat(_gc_indices[i], _gc_indices[i+1], j, j+1, workv);
      scale_inverse_inertia(i, workv);
      X.set_sub_mat(_gc_indices[i], j, workv);
    }

  return X;
}

/// Multiplies the inverse of the inertia matrix by a vector
VectorNd& MCArticulatedBody::iM_mult(const VectorNd& v, VectorNd& result) const
{
  VectorNd sub;

  // resize the result vector
  result.resize(_gc_indices.back());

  // do the multiplication
  for (unsigned i=0; i< _iM.size(); i++)
  {
    v.get_sub_vec(_gc_indices[i], _gc_indices[i+1], sub);
    scale_inverse_inertia(i, sub);
    result.set_sub_vec(_gc_indices[i], sub);
  }

  return result;
}

/// Scales a spatial six (or zero) dimensional vector using the inverse inertia
VectorNd& MCArticulatedBody::scale_inverse_inertia(unsigned i, VectorNd& v) const
{
  if (v.size() == 0 || !_links[i]->is_enabled())
    return v.set_zero();

  // TODO: fix this as necessary to account for different computation frames
  // multiply by the inverse
  v[0] *= _iM[i].inv_mass;
  v[1] *= _iM[i].inv_mass;
  v[2] *= _iM[i].inv_mass;

  // get the inverse inertia
  Vector3d w = _iM[i].inv_inertia * Vector3d(v[3], v[4], v[5]);

  // set bottom part of v
  v[3] = w[0];
  v[4] = w[1];
  v[5] = w[2];

  return v;
}

/// Returns the dense version of J
MatrixNd& MCArticulatedBody::dense_J(const SparseJacobian& J, MatrixNd& dJ) const
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);

  // resize dJ
  dJ.set_zero(J.rows(), NGC);

  // setup dJ
  for (unsigned i=0; i< J.rows(); i++)
  {
    for (unsigned j=0, k=0; j< J.indices[i].size(); j++)
    {
      // get the rigid body corresponding to the index
      RigidBodyPtr link = _links[J.indices[i][j]];
      if (!link->is_enabled())
      {
        k += SPATIAL_DIM;
        continue;
      }

      // get the absolute coord index of the rigid body
      unsigned r = _gc_indices[J.indices[i][j]];

      // set appropriate part of dJ
      dJ(i,r++) = J(i,k++);
      dJ(i,r++) = J(i,k++);
      dJ(i,r++) = J(i,k++);
      dJ(i,r++) = J(i,k++);
      dJ(i,r++) = J(i,k++);
      dJ(i,r++) = J(i,k++);
    }
  }

  return dJ;
}

/// Returns the dense version of J
MatrixNd MCArticulatedBody::dense_J(const SparseJacobian& J) const
{
  const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
  const unsigned NROWS = J.rows();
  VectorNd v(NGC), result(NROWS);

  // set v to zero initially
  v.set_zero();
 
  // setup the matrix
  MatrixNd dJ(NROWS, NGC);
  for (unsigned i=0; i< NGC; i++)
  {
    v[i] = (double) 1.0;
    mult_sparse(J, v, result);
    v[i] = (double) 0.0;
    dJ.set_column(i, result);
  }

  return dJ;
}

void to_matrix(const SpatialTransform& X, MatrixNd& m)
{
  m.set_zero(6,6);
  m.set_sub_mat(0,0,X.E);
  m.set_sub_mat(3,3,X.E);
  m.set_sub_mat(3,0,-Matrix3::skew_symmetric(X.r)*X.E);
}

/// Computes the mechanism Jacobian (debugging purposes only!)
/**
 * Jacobian is ndof x ngc where ndof is number of joint degrees-of-freedom and 
 * ngc is number of generalized (absolute, axis angle) coordinates. Also sets up
 * the time derivative of the mechanism Jacobian.
 */
/*
void MCArticulatedBody::get_mechanism_jacobian()
{
  // setup the size of the mechanism Jacobian
  unsigned ndof = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    ndof += _joints[i]->num_dof();

  // the dof to start setting up the jacobian
  unsigned dof = 0;

  // save the current link velocities
  vector<SVector6> v(_links.size());
  for (unsigned i=0; i< v.size(); i++)
    v[i] = _links[i]->get_spatial_velocity(eLink);

  // reset all generalized velocities
  VectorNd gv = VectorNd::zero(_gc_indices.back());

  // setup dq
  MatrixNd dq(ndof, _gc_indices.back());

  for (unsigned i=0; i< gv.size(); i++)
  {
    // set the generalized velocities
    if (i > 0)
      gv[i-1] = 0.0;
    gv[i] = 1.0;
    set_generalized_velocity(DynamicBody::eSpatial, gv);

    // store the new joint velocities
    VectorNd vcol(ndof);
    for (unsigned j=0; j< _joints.size(); j++)
    {
      _joints[j]->determine_q_dot();
      unsigned idx = _joints[j]->get_coord_index();
      vcol.set_sub_vec(idx, _joints[j]->qd);
    }

    // setup column of dq
    dq.set_column(i, vcol);
  }

  // reset the current link velocities
  for (unsigned i=0; i< v.size(); i++)
    _links[i]->set_spatial_velocity(v[i], eLink);
}
*/

/// Computes the "special" spatial transform S, such that spatial_transpose(S) = transpose(X) * spatial_transpose(v), for any vector v
MatrixNd& MCArticulatedBody::reverse_transform(const SpatialTransform& X, const MatrixNd& pinv_s, MatrixNd& sx)
{
  const unsigned THREE_D = 3, SPATIAL_DIM = 6;

  // setup sx
  sx.resize(pinv_s.rows(), SPATIAL_DIM);
  
  // [st sb] [ E   0 ]
  //         [ -rE E ]

  // compute skew(r) * E
  Matrix3 LL = Matrix3::skew_symmetric(-X.r) * X.E;

  // compute st * E
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, pinv_s.rows(), THREE_D, THREE_D, (double) 1.0, pinv_s.begin(), pinv_s.rows(), X.E.begin(), THREE_D, (double) 0.0, sx.begin(), sx.rows());

  // compute sb * LL
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, pinv_s.rows(), THREE_D, THREE_D, (double) 1.0, pinv_s.begin()+pinv_s.rows()*THREE_D, pinv_s.rows(), LL.begin(), THREE_D, (double) 1.0, sx.begin(), sx.rows());

  // compute sb * E 
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, pinv_s.rows(), THREE_D, THREE_D, (double) 1.0, pinv_s.begin()+pinv_s.rows()*THREE_D, pinv_s.rows(), X.E.begin(), THREE_D, (double) 0.0, sx.begin()+pinv_s.rows()*THREE_D, sx.rows());

  return sx;
}

/// Sets up the mechanism Jacobian
/**
 * Jacobian is ndof x ngc where ndof is number of joint degrees-of-freedom and 
 * ngc is number of generalized (absolute, axis angle) coordinates. Also sets up
 * the time derivative of the mechanism Jacobian.
 */
void MCArticulatedBody::get_mechanism_jacobian(MCArticulatedBody::SparseJacobian& J, MCArticulatedBody::SparseJacobian& J_dot) const
{
  const unsigned SPATIAL_DIM = 6;
  SAFESTATIC MatrixNd left, right, pinv_s, pinv_s_dot;
  SAFESTATIC MatrixNd sx;

  // setup the size of the mechanism Jacobian
  unsigned ndof = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    ndof += _joints[i]->num_dof();
  J.set_zero(ndof, SPATIAL_DIM*2);
  J_dot.set_zero(ndof, SPATIAL_DIM*2);

  // setup the indices for the Jacobians
  J.indices.resize(ndof);
  J_dot.indices.resize(ndof);
  for (unsigned i=0, k=0; i< _joints.size(); i++)
  {
    RigidBodyPtr in = _joints[i]->get_inboard_link();
    RigidBodyPtr out = _joints[i]->get_outboard_link();
    for (unsigned j=0; j< _joints[i]->num_dof(); j++, k++)
    {
      J.indices[k].resize(2);
      J_dot.indices[k].resize(2);
      J.indices[k][0] = J_dot.indices[k][0] = in->get_index();
      J.indices[k][1] = J_dot.indices[k][1] = out->get_index();
    }
  }

  // process all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // get the coordinate index start
    unsigned dof = _joints[i]->get_coord_index();

    // get the inboard and outboard link of the joint
    RigidBodyPtr rbi = _joints[i]->get_inboard_link();
    RigidBodyPtr rbo = _joints[i]->get_outboard_link();
    assert(rbi && rbo);

    // get the spatial axes- and time derivatives- for joint in global frame
    const SMatrix6N& s = _joints[i]->get_spatial_axes(eGlobal);
    const SMatrix6N& s_dot = _joints[i]->get_spatial_axes_dot(eGlobal);

    // compute the pseudo inverse of s and \dot{s}
    pinv_s.copy_from(s);
    pinv_s_dot.copy_from(s_dot);
    try
    {
      LinAlg::pseudo_inverse(pinv_s, LinAlg::svd1);
    }
    catch (NumericalException e)
    {
      pinv_s.copy_from(s);
      LinAlg::pseudo_inverse(pinv_s, LinAlg::svd2);
    }

    // compute special spatial transformations 
    SpatialTransform X0i(IDENTITY_3x3, rbi->get_position(), IDENTITY_4x4);
    SpatialTransform X0o(IDENTITY_3x3, rbo->get_position(), IDENTITY_4x4);

    // process rbi if non-NULL
    if (rbi->is_enabled())
    {
      // setup components of J
      reverse_transform(X0i, pinv_s, sx);
      sx.negate();      
      sx.get_sub_mat(0, sx.rows(), 0, 3, right);
      sx.get_sub_mat(0, sx.rows(), 3, 6, left);
      J.set_sub_mat(dof, 0, left);
      J.set_sub_mat(dof, 3, right);

      // setup components of \dot{J}
      reverse_transform(X0i, pinv_s_dot, sx);
      sx.negate();      
      sx.get_sub_mat(0, sx.rows(), 0, 3, right);
      sx.get_sub_mat(0, sx.rows(), 3, 6, left);
      J_dot.set_sub_mat(dof, 0, left);
      J_dot.set_sub_mat(dof, 3, right);
    }

    // process rbo if non-null
    if (rbo->is_enabled())
    {
      // setup components of J
      reverse_transform(X0o, pinv_s, sx);
      sx.get_sub_mat(0, sx.rows(), 0, 3, right);
      sx.get_sub_mat(0, sx.rows(), 3, 6, left);
      J.set_sub_mat(dof, 6, left);
      J.set_sub_mat(dof, 9, right);

      // setup components of \dot{J}
      reverse_transform(X0o, pinv_s_dot, sx);
      sx.get_sub_mat(0, sx.rows(), 0, 3, right);
      sx.get_sub_mat(0, sx.rows(), 3, 6, left);
      J_dot.set_sub_mat(dof, 6, left);
      J_dot.set_sub_mat(dof, 9, right);
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "sparse J: " << std::endl << J;
  FILE_LOG(LOG_DYNAMICS) << "sparse Jdot: " << std::endl << J_dot;
}

/// Sets up the constraint Jacobian (fast, analytic version)
/**
 * Jacobian is nc x ngc where nc is number of constraints and ngc is number
 * of generalized coordinates (axis angle).
 */
void MCArticulatedBody::get_constraint_jacobian(MCArticulatedBody::SparseJacobian& J) const
{
  const unsigned SPATIAL_DIM = 6;
  double CJrb1[SPATIAL_DIM], CJrb2[SPATIAL_DIM];
  unsigned gc1 = std::numeric_limits<unsigned>::max();
  unsigned gc2 = std::numeric_limits<unsigned>::max();

  // setup the size of the constraint Jacobian
  unsigned neqns = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    neqns += _joints[i]->num_constraint_eqns();
  J.set_zero(neqns, SPATIAL_DIM*2);

  // setup the indices for the Jacobian
  J.indices.resize(neqns);
  for (unsigned i=0; i< _joints.size(); i++)
  {
    RigidBodyPtr in = _joints[i]->get_inboard_link();
    RigidBodyPtr out = _joints[i]->get_outboard_link();
    unsigned idx = _joints[i]->get_constraint_index();
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
    {
      J.indices[idx+j].resize(2);
      J.indices[idx+j][0] = in->get_index();
      J.indices[idx+j][1] = out->get_index();
    }
  }

  // the row to start setting up the jacobian
  unsigned row = 0;

  // process all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // get the inboard and outboard link of the joint
    RigidBodyPtr rb1 = _joints[i]->get_inboard_link();
    RigidBodyPtr rb2 = _joints[i]->get_outboard_link();
    assert(rb1 && rb2);

    // get the index of the constraint
    unsigned cidx = _joints[i]->get_constraint_index();

    // iterate over all joint constraint equations
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
    {
      // process rb1 if non-NULL
      if (rb1->is_enabled())
      {
        // get the starting spatial index
        unsigned rb1i = rb1->get_index(); 
        gc1 = _gc_indices[rb1i];
        assert(gc1 != std::numeric_limits<unsigned>::max());

        // get the constraint jacobian
        _joints[i]->calc_constraint_jacobian(DynamicBody::eSpatial, rb1,j,CJrb1);

        // set the appropriate part of the constraint Jacobian
        for (unsigned k=0; k< SPATIAL_DIM; k++)
          J(cidx+j, k) = CJrb1[k];

        FILE_LOG(LOG_DYNAMICS) << "constraint Jacobian for joint " << _joints[i]->id << " w.r.t. body " << rb1->id << endl;
        if (LOGGING(LOG_DYNAMICS))
          FILE_LOG(LOG_DYNAMICS) << "  eqn " << j << ": " << CJrb1[0] << " " << CJrb1[1] << " " << CJrb1[2] << " " << CJrb1[3] << " " << CJrb1[4] << " " << CJrb1[5] << endl;
      }

      // process rb2 if non-null
      if (rb2->is_enabled())
      {
        // get the starting spatial index
        unsigned rb2i = rb2->get_index();
        gc2 = _gc_indices[rb2i];
        assert(gc2 != std::numeric_limits<unsigned>::max());

        // get the constraint jacobian components for this body
        _joints[i]->calc_constraint_jacobian(DynamicBody::eSpatial,rb2,j,CJrb2);

        // set the appropriate part of the constraint Jacobian
        for (unsigned k=0; k< SPATIAL_DIM; k++)
          J(cidx+j, k+SPATIAL_DIM) = CJrb2[k];

        FILE_LOG(LOG_DYNAMICS) << "constraint Jacobian for joint " << _joints[i]->id << " w.r.t. body " << rb2->id << endl;
        if (LOGGING(LOG_DYNAMICS))
          FILE_LOG(LOG_DYNAMICS) << "  eqn " << j << ": " << CJrb2[0] << " " << CJrb2[1] << " " << CJrb2[2] << " " << CJrb2[3] << " " << CJrb2[4] << " " << CJrb2[5] << endl;
      }
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "sparse J: " << std::endl << J;
}

/// Sets up the time derivative of the constraint Jacobian (fast, analytic version)
/**
 * Jacobian is nc x ngc where nc is number of constraints and ngc is number
 * of generalized coordinates (axis angle).
 */
void MCArticulatedBody::get_constraint_jacobian_dot(MCArticulatedBody::SparseJacobian& J) const
{
  const unsigned SPATIAL_DIM = 6;
  double CJrb1[SPATIAL_DIM], CJrb2[SPATIAL_DIM];
  unsigned gc1 = std::numeric_limits<unsigned>::max();
  unsigned gc2 = std::numeric_limits<unsigned>::max();

  // setup the size of the constraint Jacobian
  unsigned neqns = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    neqns += _joints[i]->num_constraint_eqns();
  J.set_zero(neqns, SPATIAL_DIM*2);

  // setup the indices for the Jacobian
  J.indices.resize(neqns);
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned cidx = _joints[i]->get_constraint_index();
    RigidBodyPtr in = _joints[i]->get_inboard_link();
    RigidBodyPtr out = _joints[i]->get_outboard_link();
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
    {
      J.indices[cidx+j].resize(2);
      J.indices[cidx+j][0] = in->get_index();
      J.indices[cidx+j][1] = out->get_index();
    }
  }

  // process all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // get the inboard and outboard link of the joint
    RigidBodyPtr rb1 = _joints[i]->get_inboard_link();
    RigidBodyPtr rb2 = _joints[i]->get_outboard_link();
    assert(rb1 && rb2);

    // get the constraint index
    unsigned idx = _joints[i]->get_constraint_index();

    // iterate over all joint constraint equations
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
    {
      // process rb1 if non-NULL
      if (rb1->is_enabled())
      {
        // get the starting spatial index
        unsigned rb1i = rb1->get_index(); 
        gc1 = _gc_indices[rb1i];
        assert(gc1 != std::numeric_limits<unsigned>::max());

        // get the constraint jacobian
        _joints[i]->calc_constraint_jacobian_dot(DynamicBody::eSpatial, rb1,j,CJrb1);

        // set the appropriate part of the constraint Jacobian
        for (unsigned k=0; k< SPATIAL_DIM; k++)
          J(idx+j, k) = CJrb1[k];

        FILE_LOG(LOG_DYNAMICS) << "time deriv. of constraint Jacobian for joint " << _joints[i]->id << " w.r.t. body " << rb1->id << endl;
        if (LOGGING(LOG_DYNAMICS))
          FILE_LOG(LOG_DYNAMICS) << "  eqn " << j << ": " << CJrb1[0] << " " << CJrb1[1] << " " << CJrb1[2] << " " << CJrb1[3] << " " << CJrb1[4] << " " << CJrb1[5] << endl;
      }

      // process rb2 if non-null
      if (rb2->is_enabled())
      {
        // get the starting spatial index
        unsigned rb2i = rb2->get_index();
        gc2 = _gc_indices[rb2i];
        assert(gc2 != std::numeric_limits<unsigned>::max());

        // get the constraint jacobian components for this body
        _joints[i]->calc_constraint_jacobian_dot(DynamicBody::eSpatial,rb2,j,CJrb2);

        // set the appropriate part of the constraint Jacobian
        for (unsigned k=0; k< SPATIAL_DIM; k++)
          J(idx+j, k+SPATIAL_DIM) = CJrb2[k];

        FILE_LOG(LOG_DYNAMICS) << "time deriv of constraint Jacobian for joint " << _joints[i]->id << " w.r.t. body " << rb2->id << endl;
        if (LOGGING(LOG_DYNAMICS))
          FILE_LOG(LOG_DYNAMICS) << "  eqn " << j << ": " << CJrb2[0] << " " << CJrb2[1] << " " << CJrb2[2] << " " << CJrb2[3] << " " << CJrb2[4] << " " << CJrb2[5] << endl;
      }
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "sparse J: " << std::endl << J;
}

/// Utility function for numerically setting up constraint Jacobian
void MCArticulatedBody::increment_dof(RigidBodyPtr rb1, RigidBodyPtr rb2, unsigned k, double h)
{
  if (k > 5)
  {
    // dealing with rb2...  check whether it is disabled
    if (!rb2->is_enabled())
      return;

    if (k < 9)
    {
      Point3d x = rb2->get_position();
      x[k-6] += h;
      rb2->set_position(x);
    }
    else
    {
      Point3d omega .set_zero()

;
      omega[k-9] = h;
      Matrix3 omega_hat = Matrix3::skew_symmetric(omega);
      Matrix3 R = rb2->get_transform().get_rotation();
      R += omega_hat * R;
      R.orthonormalize();
      Matrix4 T = rb2->get_transform();
      T.set_rotation(&R);
      rb2->set_transform(T);
    }
  }
  else
  {
    // dealing with rb1... check whether it is disabled
    if (!rb1->is_enabled())
      return;

    if (k < 3)
    {
      Point3d x = rb1->get_position();
      x[k] += h;
      rb1->set_position(x);
    }
    else
    {
      Vector3d omega .set_zero()

;
      omega[k-3] = h;
      Matrix3 omega_hat = Matrix3::skew_symmetric(omega);
      Matrix3 R = rb1->get_transform().get_rotation();
      R += omega_hat * R;
      R.orthonormalize();
      Matrix4 T = rb1->get_transform();
      T.set_rotation(&R);
      rb1->set_transform(T);
    }
  }
}

/// Sets up the constraint Jacobian (numerical version for debugging)
/**
 * Jacobian is nc x ngc where nc is number of constraints and ngc is number
 * of generalized coordinates (axis angle).
 */
void MCArticulatedBody::get_constraint_jacobian_numerically(MCArticulatedBody::SparseJacobian& J) const
{
  double C1[6], C2[6];
  const double H = 1e-6;

  // get the dynamic state of the body
  VectorNd q;
  MCArticulatedBody* this_nc = (MCArticulatedBody*) this;
  this_nc->get_generalized_coordinates(DynamicBody::eEuler, q);

  // setup the size of the constraint Jacobian
  unsigned neqns = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    neqns += _joints[i]->num_constraint_eqns();
  J.set_zero(neqns, 12);

  // setup the indices for the Jacobian
  J.indices.resize(neqns);
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned idx = _joints[i]->get_constraint_index();
    RigidBodyPtr in = _joints[i]->get_inboard_link();
    RigidBodyPtr out = _joints[i]->get_outboard_link();
    for (unsigned j=0; j< _joints[i]->num_constraint_eqns(); j++)
    {
      J.indices[idx+j].resize(2);
      J.indices[idx+j][0] = in->get_index();
      J.indices[idx+j][1] = out->get_index();
    }
  }

  // process all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned idx = _joints[i]->get_constraint_index();

    // get the inboard and outboard link of the joint
    RigidBodyPtr rb1 = _joints[i]->get_inboard_link();
    RigidBodyPtr rb2 = _joints[i]->get_outboard_link();
    assert(rb1 && rb2);

    // iterate over all DOF
    for (unsigned k=0; k< 12; k++)
    {
      // measure the constraint values after
      _joints[i]->evaluate_constraints(C1);

      // increment the k'th degree of freedom of the two bodies
      increment_dof(rb1, rb2, k, H);

      // measure the constraint values after
      _joints[i]->evaluate_constraints(C2);

      // compute the difference
      for (unsigned m=0; m< 6; m++)
      {
        C2[m] -= C1[m];
        C2[m] /= H;
      }

      // setup the Jacobian
      for (unsigned m=0; m< _joints[i]->num_constraint_eqns(); m++)
        J(idx+m, k) = C2[m];
    }
  }

  // reset the dynamic state of the body
  this_nc->set_generalized_coordinates(DynamicBody::eEuler, q);

  FILE_LOG(LOG_DYNAMICS) << "sparse J: " << std::endl << J;
}

/// sets the generalized coordinates for this body 
void MCArticulatedBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gc) 
{
  SAFESTATIC VectorNd sub_vec;

  // get the generalized coordinates
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    const unsigned NGC = _links[i]->num_generalized_coordinates_single(gctype);
    gc.get_sub_vec(k, k+NGC, sub_vec);
    _links[i]->set_generalized_coordinates_single(gctype, sub_vec);
    k += NGC;
  }
}

/// Sets the generalized velocity for this articulated body
void MCArticulatedBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gv) 
{
  SAFESTATIC VectorNd tmpv;

  // get the generalized velocities 
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    const unsigned NGC = _links[i]->num_generalized_coordinates_single(gctype);
    gv.get_sub_vec(k, k+NGC, tmpv);
    _links[i]->set_generalized_velocity_single(gctype, tmpv);
    k += NGC;
  }

  // compute the joint velocities
  get_mechanism_jacobian(_Dx, _Dx_dot);
  mult_sparse(_Dx, gv, tmpv);
  for (unsigned i=0; i< _joints.size(); i++)
  {
    unsigned idx = _joints[i]->get_coord_index();
    tmpv.get_sub_vec(idx, idx+_joints[i]->num_dof(), _joints[i]->qd);
  }
}

/// Resets the accumulators on all links
void MCArticulatedBody::reset_accumulators()
{
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->reset_accumulators();
}

/// Determines a generalized force on the body
VectorNd& MCArticulatedBody::convert_to_generalized_force(SingleBodyPtr link, const SForced& w, const Point3d& p, VectorNd& gf)
{
  SAFESTATIC VectorNd link_gf;

  // setup a zero vector of all generalized forces
  gf.set_zero(num_generalized_coordinates(DynamicBody::eSpatial));

  // get the link as a rigid body
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(link);
  assert(rb);

  // get the generalized force for the link
  link->convert_to_generalized_force_single(link, w, link_gf);

  // set the appropriate part in the vector
  for (unsigned i=0, k=0; i< _links.size(); i++)
  {
    if (_links[i] == rb)
    {
      gf.set_sub_vec(k, link_gf);
      return gf;
    }
    k += _links[i]->num_generalized_coordinates_single(DynamicBody::eSpatial);
  }

  assert(false);
  return gf;
}

/// Loads a MCArticulatedBody object from an XML node
void MCArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // load the parent data
  ArticulatedBody::load_from_xml(node, id_map);

  // verify the node name -- this class has derived classes
  assert(strcasecmp(node->name.c_str(), "MCArticulatedBody") == 0);

  // get baumgarte paramters
  const XMLAttrib* balpha_attr = node->get_attrib("baumgarte-alpha");
  if (balpha_attr)
    b_alpha = balpha_attr->get_real_value();
  const XMLAttrib* bbeta_attr = node->get_attrib("baumgarte-beta");
  if (bbeta_attr)
    b_beta = bbeta_attr->get_real_value();
}

/// Saves this object to a XML tree
void MCArticulatedBody::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent method first
  ArticulatedBody::save_to_xml(node, shared_objects);  

  // (re)set the name of this node
  node->name = "MCArticulatedBody";

  // save baumgarte parameters
  node->attribs.insert(XMLAttrib("baumgarte-alpha", b_alpha));
  node->attribs.insert(XMLAttrib("baumgarte-beta", b_beta));
}

/// Updates the event data for Jx_v()
void MCArticulatedBody::update_Jx_v(EventProblemData& q)
{
  const unsigned NGC_ROD = 7;
  double Cq[NGC_ROD];
  Vector3d top, bottom;

  // update Jx_v
  for (unsigned i=0, ii=0; i< q.constraint_events.size(); i++)
  {
     // get the i'th constraint articulated body -- see whether to skip 
    // constraint i
    JointPtr joint = q.constraint_events[i]->constraint_joint;
    RigidBodyPtr outboard = joint->get_outboard_link();
    ArticulatedBodyPtr ab = outboard->get_articulated_body();
    assert(ab);
    if (ab.get() != this)
    {
      ii += joint->num_constraint_eqns();
      continue;
    }

    // get the inboard body of this constraint
    RigidBodyPtr inboard = joint->get_inboard_link();

    // get the linear and angular velocities of the inboard and outboard links
    const Vector3d& vi = inboard->get_lvel();
    const Vector3d& omegai = inboard->get_avel();
    const Vector3d& vo = outboard->get_lvel();
    const Vector3d& omegao = outboard->get_avel();

    // do the actual updating 
    for (unsigned m=0; m< joint->num_constraint_eqns(); m++)
    {
      joint->calc_constraint_jacobian(DynamicBody::eSpatial, inboard, m, Cq);
      q.Jx_v[ii+m] += vi.dot(Vector3d(Cq[0],Cq[1],Cq[2]));
      q.Jx_v[ii+m] += omegai.dot(Vector3d(Cq[3], Cq[4], Cq[5]));
      joint->calc_constraint_jacobian(DynamicBody::eSpatial, outboard, m, Cq);
      q.Jx_v[ii+m] += vo.dot(Vector3d(Cq[0],Cq[1],Cq[2]));
      q.Jx_v[ii+m] += omegao.dot(Vector3d(Cq[3], Cq[4], Cq[5]));
    }

    // update ii
    ii += joint->num_constraint_eqns();
  }
}

// update the joint limit event velocities
void MCArticulatedBody::update_Jl_v(EventProblemData& q)
{
  for (unsigned i=0; i< q.limit_events.size(); i++)
  {
    // get the i'th limit articulated body -- see whether to skip limit i
    JointPtr joint = q.limit_events[i]->limit_joint;
    ArticulatedBodyPtr ab = joint->get_articulated_body();
    assert(ab);
    if (ab.get() != this)
      continue;

    // update Jl_v
    if (q.limit_events[i]->limit_upper)
      q.Jl_v[i] -= joint->qd[q.limit_events[i]->limit_dof];
    else
      q.Jl_v[i] += joint->qd[q.limit_events[i]->limit_dof];
  }
}

// update the joint friction event velocities
void MCArticulatedBody::update_Dx_v(EventProblemData& q)
{
  for (unsigned i=0; i< q.constraint_events.size(); i++)
  {
    // get the i'th constraint articulated body -- see whether to skip event i
    JointPtr joint = q.constraint_events[i]->constraint_joint;
    ArticulatedBodyPtr ab = joint->get_articulated_body();
    assert(ab);
    if (ab.get() != this)
      continue;

    // update Dx_v for all constraint equations
    const unsigned NDOF = q.constraint_events[i]->constraint_joint->num_dof();
    for (unsigned j=0; j< NDOF; j++)
      q.Dx_v[i] += joint->qd[q.limit_events[i]->limit_dof];
  }
}

/// Gets whether the event affects the rigid body
bool MCArticulatedBody::affects(RigidBodyPtr rb, Event* e)
{
  if (e->event_type == Event::eNone)
    return false;
  else if (e->event_type == Event::eContact)
  {
    SingleBodyPtr sb1 = e->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e->contact_geom2->get_single_body();
    return (rb == sb1 || rb == sb2);
  }
  else if (e->event_type == Event::eLimit)
  {
    JointPtr joint = e->limit_joint;
    RigidBodyPtr inboard = joint->get_inboard_link();
    RigidBodyPtr outboard = joint->get_outboard_link();
    return (rb == inboard || rb == outboard);
  }
  else if (e->event_type == Event::eConstraint)
  {
    JointPtr joint = e->constraint_joint;
    RigidBodyPtr inboard = joint->get_inboard_link();
    RigidBodyPtr outboard = joint->get_outboard_link();
    return (rb == inboard || rb == outboard);
  }
  else
    assert(false);

  // make compiler happy
  return false;
}

/// Gets the number of subevents out of an event
unsigned MCArticulatedBody::num_sub_events(JacobianType jt, Event* e)
{
  switch (jt)
  {
    case MCArticulatedBody::eNone:
      return 0;

    case MCArticulatedBody::eContactNormal:
      return 1;

    case MCArticulatedBody::eContactTangent:
      return 2;

    case MCArticulatedBody::eLimit:
      return 1;

    case MCArticulatedBody::eJointFriction:
      return e->constraint_joint->num_dof();

    case MCArticulatedBody::eConstraint:
      if (e->constraint_joint->get_constraint_type() == Joint::eEXPLICIT)
        return 0;
      else
      {
        assert(e->constraint_joint->get_constraint_type() == Joint::eImplicit);
        return e->constraint_joint->num_constraint_eqns();
      }
  }

  // shouldn't still be here; nonetheless, keep compiler happy
  return 0;
}

/// Gets the particular data out of an event
void MCArticulatedBody::get_event_data(JacobianType jt, Event* e, RigidBodyPtr rb, unsigned subidx, Vector3d& tx, Vector3d& rx)
{
  const unsigned NGC_ROD = 7;
  double Cq[NGC_ROD];

  if (jt == MCArticulatedBody::eNone)
  {
    tx .set_zero()

;
    rx .set_zero()

;
    return;
  }
  else if (jt == MCArticulatedBody::eContactNormal)
  {
    assert(subidx == 0);

    // get the contact point and normal
    const Point3d& cp = e->contact_point;
    const Vector3d& normal = e->contact_normal;    
    Vector3d r = cp - rb->get_position();
    Vector3d cross = Vector3d::cross(r, normal); 

    // figure out which body maps to rb
    SingleBodyPtr sb1 = e->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e->contact_geom2->get_single_body();
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(sb1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(sb2);
    if (rb1 == rb)
    {
      tx = normal;
      rx = cross;
    }
    else if (rb2 == rb)
    {
      tx = -normal;
      rx = -cross;
    }
    else
      assert(false);
  }
  else if (jt == MCArticulatedBody::eContactTangent)
  {
    assert(subidx < 2);

    // get the contact point and normal
    const Point3d& cp = e->contact_point;
    const Vector3d& ctan = (subidx == 0) ? e->contact_tan1 : e->contact_tan2;    
    Vector3d r = cp - rb->get_position();
    Vector3d cross = Vector3d::cross(r, ctan); 

    // figure out which body maps to rb
    SingleBodyPtr sb1 = e->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e->contact_geom2->get_single_body();
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(sb1);
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(sb2);
    if (rb1 == rb)
    {
      tx = ctan;
      rx = cross;
    }
    else if (rb2 == rb)
    {
      tx = -ctan;
      rx = -cross;
    }
    else
      assert(false);
  }
  else if (jt == MCArticulatedBody::eLimit)
  {
    assert(subidx == 0); 

    // determine whether to negate
    bool negate = e->limit_upper;

    // get the joint motion
    const SMatrix6N& si = e->limit_joint->get_spatial_axes(eLink);

    // get the rotation and translation
    SVector6 Jcol = si.get_column(e->limit_dof);
    tx = Jcol.get_lower();
    rx = Jcol.get_upper();

    // transform them to the proper frame
    RigidBodyPtr outboard = e->limit_joint->get_outboard_link(); 
    const Matrix4& To = outboard->get_transform();
    tx = To.mult_vector(tx);
    rx = To.mult_vector(rx);

    // figure out which body maps to rb
    RigidBodyPtr rb1 = e->limit_joint->get_inboard_link(); 
    RigidBodyPtr rb2 = outboard; 
    if (rb2 == rb)
      negate = !negate;
    else
      assert(rb1 == rb);

    // see whether to negate
    if (negate)
    {
      tx = -tx;
      rx = -rx;
    }
  }
  else if (jt == MCArticulatedBody::eJointFriction)
  {
    assert(subidx <= e->constraint_joint->num_dof());

    // get the joint motion
    const SMatrix6N& si = e->constraint_joint->get_spatial_axes(eLink);

    // get the rotation and translation
    SVector6 Jcol = si.get_column(subidx);
    tx = Jcol.get_lower();
    rx = Jcol.get_upper();

    // transform them to the proper frame
    RigidBodyPtr outboard = e->constraint_joint->get_outboard_link(); 
    const Matrix4& To = outboard->get_transform();
    tx = To.mult_vector(tx);
    rx = To.mult_vector(rx);

    // figure out which body maps to rb
    RigidBodyPtr rb1 = e->constraint_joint->get_inboard_link(); 
    RigidBodyPtr rb2 = outboard; 
    if (rb2 == rb)
    {
      tx = -tx;
      rx = -rx;
    }
    else
      assert(rb1 == rb);
  }
  else if (jt == MCArticulatedBody::eConstraint)
  {
    assert(subidx <= e->constraint_joint->num_constraint_eqns());

    // get the column of the constraint Jacobian
    e->constraint_joint->calc_constraint_jacobian(DynamicBody::eSpatial, rb, subidx, Cq);

    // copy to vectors 
    tx[0] = Cq[0];
    tx[1] = Cq[1];
    tx[2] = Cq[2];
    rx[0] = Cq[3];
    rx[1] = Cq[4];
    rx[2] = Cq[5];
  }
  else
    assert(false);
}

/// Gets the proper vector of events
const vector<Event*>& MCArticulatedBody::get_events_vector(const EventProblemData& q, JacobianType jt)
{
  switch (jt)
  {
    case MCArticulatedBody::eContactNormal:
    case MCArticulatedBody::eContactTangent:
      return q.contact_events;

    case MCArticulatedBody::eLimit:
      return q.limit_events;

    case MCArticulatedBody::eJointFriction:
    case MCArticulatedBody::eConstraint:
      return q.constraint_events;

    // return some nonsense as we were given some nonsense
    default:
    case MCArticulatedBody::eNone:
      return q.contact_events;
  }
}

// updates any Jacobian matrix times inverse inertia times any other Jacobian matrix
void MCArticulatedBody::update_Jx_iM_JyT(EventProblemData& q, MatrixNd& Jx_iM_JyT, JacobianType j1t, JacobianType j2t)
{
  for (unsigned i=0; i< _links.size(); i++)
    update_Jx_iM_JyT(_links[i], q, Jx_iM_JyT, j1t, j2t);
}

// updates any Jacobian matrix times inverse inertia times any other Jacobian matrix (helper method)
void MCArticulatedBody::update_Jx_iM_JyT(RigidBodyPtr rb, EventProblemData& q, MatrixNd& Jx_iM_JyT, JacobianType j1t, JacobianType j2t)
{
  Vector3d t1, t2, r1, r2;

  // get the rigid body and its necessary properties
  Matrix3 R(&rb->get_orientation());
  Matrix3 Jinv = R * rb->get_inv_inertia() * Matrix3::transpose(R);
  const double inv_mass = rb->get_inv_mass();

  // get the first vector and second vector of events
  const vector<Event*>& events1 = get_events_vector(q, j1t);
  const vector<Event*>& events2 = get_events_vector(q, j2t);
  
  // loop over all of event #1 type
  for (unsigned i=0, ii=0; i< events1.size(); i++)
  {
    // get the number of subevents
    const unsigned SUB1 = num_sub_events(j1t, events1[i]);

    // see whether rb is coincident on the event
    if (!affects(rb, events1[i]))
    {
      ii += SUB1;
      continue;
    }

    // loop over all subevents
    for (unsigned k=0; k< SUB1; k++)
    {
      // get the data for the event
      get_event_data(j1t, events1[i], rb, k, t1, r1);

      // scale t1 and r1 appropriately
      t1 *= inv_mass;
      r1 = Jinv * r1;

      // loop over all of event #2 type
      for (unsigned j=0, jj=0; j< events2.size(); j++)
      {
        // get the number of subevents
        const unsigned SUB2 = num_sub_events(j2t, events2[j]);

        // see whether rb is coincident on the event
        if (!affects(rb, events2[j]))
        {
          jj += SUB2;
          continue;
        }

        // loop over all subevents
        for (unsigned m=0; m< SUB2; m++)
        {
          // get the data for the event
          get_event_data(j2t, events2[j], rb, m, t2, r2);

          // set the appropriate part of the matrix
          Jx_iM_JyT(ii+k, jj+m) += t1.dot(t2) + r1.dot(r2);
        }

        // update jj
        jj += SUB2;
      }
    }

    // update ii
    ii += SUB1;
  }
}

/// Updates the event data
void MCArticulatedBody::update_event_data(EventProblemData& q)
{
  // NOTE: we don't touch q.Ji

  // first, update contact data matrices (Jc_iM_JcT, Jc_iM_DcT, Dc_iM_DcT)
  // and vectors Jc_v and Dc_v
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->update_event_data(q);

  // do joint limit updates
  update_Jl_v(q);
  update_Jx_iM_JyT(q, q.Jc_iM_JlT, eContactNormal, eLimit);
  update_Jx_iM_JyT(q, q.Dc_iM_JlT, eContactTangent, eLimit);
  update_Jx_iM_JyT(q, q.Jl_iM_JlT, eLimit, eLimit);
  update_Jx_iM_JyT(q, q.Jl_iM_DxT, eLimit, eJointFriction);
  update_Jx_iM_JyT(q, q.Jl_iM_JxT, eLimit, eConstraint);

  // do joint constraint updates
  update_Jx_v(q);
  update_Jx_iM_JyT(q, q.Jc_iM_JxT, eContactNormal, eConstraint);
  update_Jx_iM_JyT(q, q.Dc_iM_JxT, eContactTangent, eConstraint);
  update_Jx_iM_JyT(q, q.Jx_iM_JxT, eConstraint, eConstraint);
  update_Jx_iM_JyT(q, q.Jx_iM_DxT, eConstraint, eJointFriction);

  // do joint friction updates
  update_Dx_v(q);
  update_Jx_iM_JyT(q, q.Jc_iM_DxT, eContactNormal, eJointFriction);
  update_Jx_iM_JyT(q, q.Dc_iM_JlT, eContactTangent, eJointFriction);
  update_Jx_iM_JyT(q, q.Dx_iM_DxT, eJointFriction, eJointFriction);
}

/// Selects a subset from sparse Jacobians
void MCArticulatedBody::select_sub_contact_Jacobians(const EventProblemData& q, SparseJacobian& Jc_sub, SparseJacobian& Dc_sub) const
{
  SAFESTATIC vector<unsigned> normal_indices, tangent_indices;

  // setup normal and tangent contact indices
  normal_indices.clear();
  tangent_indices.clear();
  for (unsigned i=0, j=0; i< q.contact_working_set.size(); i++, j+=2)
    if (q.contact_working_set[i])
    {
      normal_indices.push_back(i);
      tangent_indices.push_back(j);
      tangent_indices.push_back(j+1);
    }

  // select the appropriate rows 
  _Jc.select_rows(normal_indices.begin(), normal_indices.end(), Jc_sub); 
  _Dc.select_rows(tangent_indices.begin(), tangent_indices.end(), Dc_sub); 

  // setup the row indices
  select(_Jc.indices.begin(), normal_indices.begin(), normal_indices.end(), std::back_inserter(Jc_sub.indices));
  select(_Dc.indices.begin(), tangent_indices.begin(), tangent_indices.end(), std::back_inserter(Dc_sub.indices));
}

// Updates the body velocities
void MCArticulatedBody::update_velocity(const EventProblemData& q)
{
  SAFESTATIC VectorNd workv, workv2;
  SAFESTATIC SparseJacobian Jc_sub, Dc_sub;

  // get the current absolute velocities
  get_generalized_velocity(DynamicBody::eSpatial, _xd);

  // determine whether we are using a subset of the contacts
  if (std::find(q.contact_working_set.begin(), q.contact_working_set.end(), false) != q.contact_working_set.end())
  {
    select_sub_contact_Jacobians(q, Jc_sub, Dc_sub);
    mult_transpose_sparse(Jc_sub, q.alpha_c, workv);
    workv += mult_transpose_sparse(Dc_sub, q.beta_c, workv2); 
  }
  else
  {
    mult_transpose_sparse(_Jc, q.alpha_c, workv);
    workv += mult_transpose_sparse(_Dc, q.beta_c, workv2);
  }

  // compute change in velocities for events not reduced
  workv += mult_transpose_sparse(_Jl, q.alpha_l, workv2);
  workv += mult_transpose_sparse(_Jx, q.alpha_x, workv2);
  workv += mult_transpose_sparse(_Dx, q.beta_x, workv2);
  _xd += iM_mult(workv, workv2);

  // update the link velocities
  update_link_velocities(); 
}

/// The signum function, modified to use NEAR_ZERO instead of +/- 0.0
double MCArticulatedBody::sgn(double x)
{
  if (x > NEAR_ZERO)
    return (double) 1.0;
  else if (x < -NEAR_ZERO)
    return (double) -1.0;
  else
    return (double) 0.0;
}

/// Sets up a sparse Jacobian multiplied by the inverse inertia matrix multiplied by the transpose of a sparse Jacobian
MatrixNd& MCArticulatedBody::calc_Jx_iM_JyT(const SparseJacobian& Jx, const SparseJacobian& Jy, MatrixNd& Jx_iM_JyT) const
{
  const unsigned SPATIAL_DIM = 6;

  // setup the matrix Jx_iM_Jx'
  Jx_iM_JyT.set_zero(Jx.rows(), Jy.rows());

  // loop over all links 
  for (unsigned k=0; k< _links.size(); k++)
  {
    // skip disabled bodies
    if (!_links[k]->is_enabled())
      continue;

    // get the inverse mass and inverse inertia (global frame)
    double inv_mass = _links[k]->get_inv_mass();
    Matrix3 R(&_links[k]->get_orientation());
    Matrix3 J_inv = R * _links[k]->get_inv_inertia() * Matrix3::transpose(R);

    // loop over all rows of Jx and Jy 
    for (unsigned rowJx=0; rowJx < Jx.rows(); rowJx++)
    {
      // loop over all blocks in this row of Jx
      for (unsigned m=0, r=0; m< Jx.indices[rowJx].size(); m++, r+= SPATIAL_DIM)
      {
        // get the rigid body responsible for this block of Jx
        if (Jx.indices[rowJx][m] != k)
          continue;

        // get the velocity vectors
        Vector3d l1(Jx(rowJx, r+0), Jx(rowJx, r+1), Jx(rowJx, r+2));
        Vector3d a1(Jx(rowJx, r+3), Jx(rowJx, r+4), Jx(rowJx, r+5));

        // scale the spatial vectors using the inverse inertias
        l1 *= inv_mass;
        a1 = J_inv * a1;

        // loop over all rows of Jy 
        for (unsigned rowJy=0; rowJy < Jy.rows(); rowJy++)
        {
          // loop over all blocks in this row of Jy
          for (unsigned n=0, s=0; n< Jy.indices[rowJy].size(); n++, s+= SPATIAL_DIM)
          {
            // get the rigid body responsible for this row of Jy'
            if (Jy.indices[rowJy][n] != k)
              continue;

            // get the velocity vectors
            Vector3d l2(Jy(rowJy, s+0), Jy(rowJy, s+1), Jy(rowJy, s+2));
            Vector3d a2(Jy(rowJy, s+3), Jy(rowJy, s+4), Jy(rowJy, s+5));

            // compute the dot product
            Jx_iM_JyT(rowJx, rowJy) += l1.dot(l2) + a1.dot(a2);
          }
        }
      }
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "J*inv(M)*J': " << std::endl << Jx_iM_JyT;

  return Jx_iM_JyT;
}

/// Multiplies a sparse jacobian by a vector
VectorNd& MCArticulatedBody::mult_sparse(const SparseJacobian& J, const VectorNd& v, VectorNd& result) const
{
  const unsigned SPATIAL_DIM = 6;

  // resize the result vector
  result.set_zero(J.rows());
  assert(J.rows() == J.indices.size());

  // carry out the multiplication
  for (unsigned i=0; i< J.indices.size(); i++)
    for (unsigned j=0; j< J.indices[i].size(); j++)
    {
      // get the rigid body corresponding to the index
      RigidBodyPtr link = _links[J.indices[i][j]];
      if (!link->is_enabled())
        continue;

      // get the gc indices corresponding to the link
      unsigned gcidx = _gc_indices[link->get_index()];
      unsigned Jidx = j*SPATIAL_DIM;
      for (unsigned k=0; k< SPATIAL_DIM; k++)
        result[i] += J(i, Jidx++) * v[gcidx++];
    }

  return result; 
}

/// Multiplies the transpose of a sparse jacobian by a vector
MatrixNd& MCArticulatedBody::mult_transpose_sparse(const SparseJacobian& J, const MatrixNd& v, MatrixNd& result) const
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned NGC = _gc_indices.back();

  // resize the result matrix
  result.set_zero(NGC, v.columns());
  assert(J.rows() == J.indices.size());

  // loop over all links
  for (unsigned r=0; r< v.columns(); r++)
    for (unsigned i=0; i< _links.size(); i++)
    {
      if (!_links[i]->is_enabled())
        continue;

      // get the starting gc index of the link
      unsigned gcidx = _gc_indices[i];

      // loop over all indices of the sparse Jacobian, looking for this link
      for (unsigned j=0; j< J.indices.size(); j++)
        for (unsigned k=0; k< J.indices[j].size(); k++)
        {
          if (J.indices[j][k] != i)
            continue;
 
          // do the multiplication
          const unsigned Jidx = k*SPATIAL_DIM;
          for (unsigned m=0; m< SPATIAL_DIM; m++)
            result(gcidx+m,r) += J(j, Jidx+m) * result(j,r);
      }
  }

  return result; 
}

/// Multiplies the transpose of a sparse jacobian by a vector
VectorNd& MCArticulatedBody::mult_transpose_sparse(const SparseJacobian& J, const VectorNd& v, VectorNd& result) const
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned NGC = _gc_indices.back();

  // resize the result vector
  result.set_zero(NGC);
  assert(J.rows() == J.indices.size());

  // loop over all links
  for (unsigned i=0; i< _links.size(); i++)
  {
    if (!_links[i]->is_enabled())
      continue;

    // get the starting gc index of the link
    unsigned gcidx = _gc_indices[i];

    // loop over all indices of the sparse Jacobian, looking for this link
    for (unsigned j=0; j< J.indices.size(); j++)
      for (unsigned k=0; k< J.indices[j].size(); k++)
      {
        if (J.indices[j][k] != i)
          continue;

        // do the multiplication
        const unsigned Jidx = k*SPATIAL_DIM;
        for (unsigned m=0; m< SPATIAL_DIM; m++)
          result[gcidx+m] += J(j, Jidx+m) * v[j];
      }
  }

  return result; 
}


