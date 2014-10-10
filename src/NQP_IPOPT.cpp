/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <numeric>
#include <Moby/NQP_IPOPT.h>

extern "C"
{
void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);
void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
}

using boost::shared_array;
using namespace Moby;
using namespace Ravelin;
using namespace Ipopt;

static double sqr(double x) { return x*x; }

NQP_IPOPT::NQP_IPOPT()
{
  _ZERO_EPS = NEAR_ZERO;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)   
{
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_ACT_CONTACTS = epd->N_ACT_CONTACTS;
  const unsigned N_LIMITS = epd->N_LIMITS;
  n = N_ACT_CONTACTS*3 + N_LIMITS;
  m = N_CONTACTS + N_ACT_CONTACTS + N_LIMITS + 1;
  MatrixNd t1, t2;

  // compute nonzeros for objective part of Hessian 
  nnz_h_lag = 0;
  _nnz_h_obj = 0;
  _nnz_h_con.clear();
  _nnz_h_con.resize(N_ACT_CONTACTS, 0);

  // clear constraint h indices
  _h_con_indices.clear();
  _h_con_indices.resize(N_ACT_CONTACTS);
  _nnz_h_con.resize(N_ACT_CONTACTS);
  _h_con.resize(N_ACT_CONTACTS);
  std::map<std::pair<unsigned, unsigned>, unsigned> h_obj_nz_indices;

  // determine whether the problem is dense or not
  _dense = (R.columns() > 0);

  // compute nonzeros for (friction) constraint part of Hessian
  if (!_dense)
  {
    // compute nonzeros for objective part of Hessian
    for (unsigned i=0, k=0; i< n; i++)
      for (unsigned j=i; j< n; j++)
        if (std::fabs(H(i,j)) > _ZERO_EPS)
        {
          h_obj_nz_indices[std::make_pair(i,j)] = k++;
          _nnz_h_obj++;
          nnz_h_lag++;
        }

    // compute nonzeros for Coulomb friction part of Hessian
    for (unsigned i=0; i< N_ACT_CONTACTS; i++)
    {
      const unsigned N_IDX = i;
      const unsigned S_IDX = i + N_ACT_CONTACTS;
      const unsigned T_IDX = i + N_ACT_CONTACTS*2;

      // setup nnz_h_con
      _nnz_h_con[i] = 3;

      // setup constraint indices
      _h_con_indices[i].resize(3);
      _h_con_indices[i][0] = h_obj_nz_indices[std::make_pair(N_IDX,N_IDX)];
      _h_con_indices[i][1] = h_obj_nz_indices[std::make_pair(S_IDX,S_IDX)];
      _h_con_indices[i][2] = h_obj_nz_indices[std::make_pair(T_IDX,T_IDX)];

      // setup values for constraint Hessian
      _h_con[i] = shared_array<double>(new double[3]);
      _h_con[i][0] = 2.0*mu_c[i];
      _h_con[i][1] = -2.0;
      _h_con[i][2] = -2.0;
    }

    // setup objective component of Hessian and indices
    _h_obj = shared_array<double>(new double[_nnz_h_obj]);
    _h_iRow = shared_array<unsigned>(new unsigned[nnz_h_lag]);  
    _h_jCol = shared_array<unsigned>(new unsigned[nnz_h_lag]);
    for (unsigned i=0, k=0; i< n; i++)
      for (unsigned j=i; j< n; j++)
        if (std::fabs(H(i,j)) > _ZERO_EPS)
        {
          _h_obj[k] = H(i,j);
          _h_iRow[k] = i;
          _h_jCol[k] = j;
          k++;
        }
  }
  else
  {
    // setup the constraint indices
    for (unsigned k=0; k< N_ACT_CONTACTS; k++)
    {
      _h_con_indices[k].clear();
      for (unsigned i=0, r=0; i< n; i++)
        for (unsigned j=i; j< n; j++, r++)
          _h_con_indices[k][r] = r;
    }

    // setup objective component of Hessian and indices
    // NOTE we assume that the Hessian is dense
    _nnz_h_obj = n*(n-1);
    _h_obj = shared_array<double>(new double[_nnz_h_obj]);
    _h_iRow = shared_array<unsigned>(new unsigned[nnz_h_lag]);  
    _h_jCol = shared_array<unsigned>(new unsigned[nnz_h_lag]);
    for (unsigned i=0, k=0; i< n; i++)
      for (unsigned j=i; j< n; j++, k++)
      {
        _h_iRow[k] = i;
        _h_jCol[k] = j;
        _h_obj[k] = H(i,j);
      }

    // copy the constraint Hessians
    for (unsigned i=0; i< N_ACT_CONTACTS; i++)
    {
      const unsigned N_IDX = i;
      const unsigned S_IDX = i + N_ACT_CONTACTS;
      const unsigned T_IDX = i + N_ACT_CONTACTS*2;

      // compute the constraint Hessian
      SharedConstVectorNd rn = R.row(N_IDX);
      SharedConstVectorNd rs = R.row(S_IDX);
      SharedConstVectorNd rt = R.row(T_IDX);
      Opsd::outer_prod(rs, rs, t1);
      Opsd::outer_prod(rt, rt, t2);
      t1 += t2;
      t1 *= -2.0;
      Opsd::outer_prod(rn, rn, t2);
      t2 *= (2.0 * mu_c[N_IDX]);
      t1 += t2;

      // copy the constraint Hessian
      _h_con[i] = shared_array<double>(new double[n*(n-1)]);
      for (unsigned i=0, k=0; i< t1.rows(); i++)
        for (unsigned j=i; j< t1.columns(); j++, k++)
          _h_con[i][k] = t1(i,j);
    }
  }

  // setup values for kappa constraint
  _cJac_constant = n;

  // count remaining non-zeros
  if (!_dense)
  {
    // count number of non-zeros in the Cn_block 
    for (RowIteratord_const i = Cn_block.row_iterator_begin(); i != i.end(); i++)
      if (std::fabs(*i) > _ZERO_EPS)
        _cJac_constant++; 

    // determine number of non-zeros in L_block
    for (RowIteratord_const i = L_block.row_iterator_begin(); i != i.end(); i++)
      if (std::fabs(*i) > _ZERO_EPS)
        _cJac_constant++; 
  }
  else
  {
    // count number of non-zeros in the Cn_block 
    Cn_block.mult(R, _M);
    for (RowIteratord_const i = _M.row_iterator_begin(); i != i.end(); i++)
      if (std::fabs(*i) > _ZERO_EPS)
        _cJac_constant++; 

    // determine number of non-zeros in L_block
    L_block.mult(R, _M);
    for (RowIteratord_const i = _M.row_iterator_begin(); i != i.end(); i++)
      if (std::fabs(*i) > _ZERO_EPS)
        _cJac_constant++; 
  }

  // setup total number of nonzero Jacobian values
  nnz_jac_g = _cJac_constant + ((!_dense) ? N_ACT_CONTACTS*3 : N_ACT_CONTACTS*n);

  // setup constant Jacobian values
  _cJac = shared_array<double>(new double[_cJac_constant]);
  _cJac_iRow = shared_array<unsigned>(new unsigned[nnz_jac_g]);
  _cJac_jCol = shared_array<unsigned>(new unsigned[nnz_jac_g]);
  unsigned nv = 0;

  // count remaining non-zeros
  if (!_dense)
  {
    // count number of non-zeros in the Cn_block 
    for (unsigned i = 0; i < Cn_block.rows(); i++)
      for (unsigned j = 0; j < Cn_block.columns(); j++)
        if (std::fabs(Cn_block(i,j)) > _ZERO_EPS)
        {
          _cJac[nv] = Cn_block(i,j);
          _cJac_iRow[nv] = i; 
          _cJac_jCol[nv] = j;
          nv++;
        }

    // determine number of non-zeros in L_block
    for (unsigned i = 0; i < L_block.rows(); i++)
      for (unsigned j = 0; j < L_block.columns(); j++)
        if (std::fabs(L_block(i,j)) > _ZERO_EPS)
        {
          _cJac[nv] = L_block(i,j);
          _cJac_iRow[nv] = i + Cn_block.rows(); 
          _cJac_jCol[nv] = j;
          nv++;
        }

    // setup -1'*N constraint
    for (unsigned i=0; i< Cn_block.columns(); i++)
    {
      RowIteratord_const ci = Cn_block.column(i).row_iterator_begin();
      _cJac[nv] = -std::accumulate(ci, ci.end(), 0.0);
      _cJac_iRow[nv] = Cn_block.rows() + L_block.rows();
      _cJac_jCol[nv] = i;
      nv++;
    }

    // setup *only indices* for Coulomb friction constraint
    for (unsigned i=0; i< N_ACT_CONTACTS; i++)
    {
      const unsigned N_IDX = i;
      const unsigned S_IDX = i + N_ACT_CONTACTS;
      const unsigned T_IDX = i + N_ACT_CONTACTS*2;
      _cJac_iRow[nv] = Cn_block.rows() + L_block.rows() + 1 + i;
      _cJac_jCol[nv] = N_IDX;
      nv++;       
      _cJac_iRow[nv] = Cn_block.rows() + L_block.rows() + 1 + i;
      _cJac_jCol[nv] = S_IDX;
      nv++;       
      _cJac_iRow[nv] = Cn_block.rows() + L_block.rows() + 1 + i;
      _cJac_jCol[nv] = T_IDX;
      nv++;       
    }
  }  
  else
  {
    // count number of non-zeros in the Cn_block 
    Cn_block.mult(R, _M);
    for (unsigned i = 0; i < _M.rows(); i++)
      for (unsigned j = 0; j < _M.columns(); j++)
        if (std::fabs(_M(i,j)) > _ZERO_EPS)
        {
          _cJac[nv] = _M(i,j);
          _cJac_iRow[nv] = i; 
          _cJac_jCol[nv] = j;
          nv++;
        }

    // determine number of non-zeros in L_block
    L_block.mult(R, _M);
    for (unsigned i = 0; i < _M.rows(); i++)
      for (unsigned j = 0; j < _M.columns(); j++)
        if (std::fabs(_M(i,j)) > _ZERO_EPS)
        {
          _cJac[nv] = _M(i,j);
          _cJac_iRow[nv] = i + Cn_block.rows(); 
          _cJac_jCol[nv] = j;
          nv++;
        }

    // setup -1'*N constraint
    Cn_block.mult(R, _M);  
    for (unsigned i=0; i< N_CONTACTS; i++)
    {
      ColumnIteratord_const ci = _M.row(i).column_iterator_begin();
      _cJac[nv] = -std::accumulate(ci, ci.end(), 0.0);
      _cJac_iRow[nv] = Cn_block.rows() + L_block.rows();
      _cJac_jCol[nv] = i;
      nv++;
    }

    // setup *only* indices for Coulomb friction constraint
    // NOTE: we are assuming that R is dense
    for (unsigned i=0; i< N_ACT_CONTACTS; i++)
      for (unsigned j=0; j< n; j++)
      {
        _cJac_iRow[nv] = N_CONTACTS + N_LIMITS + 1 + i;
        _cJac_jCol[nv] = j;
        nv++;
      }
  }

  // C-indexing style
  index_style = TNLP::C_STYLE;

  return true;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
  const double INF = std::numeric_limits<double>::max();

  // get info
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_ACT_CONTACTS = epd->N_ACT_CONTACTS;
  const unsigned N_LIMITS = epd->N_LIMITS;

  // TODO: fix bounds for dense problem

  // set variable bounds
  for (unsigned i=0; i< n; i++)
  {
    // see whether this is a normal contact force or limit force
    if (i >= N_ACT_CONTACTS && i < N_ACT_CONTACTS*3)
      x_l[i] = -INF;
    else
      x_l[i] = 0.0;
    x_u[i] = INF;
  }

  // set bounds on inequality constraints 
  for (unsigned i=0; i< N_CONTACTS+N_ACT_CONTACTS+N_LIMITS+1; i++)
  {
    g_l[i] = 0.0;
    g_u[i] = INF;
  }

  return true;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lam)
{
  if (init_x)
    std::fill_n(x, n, 0.0);

  if (init_z)
    return false;

  if (init_lambda)
    std::fill_n(lam, m, 0.0);

  return true;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // get temporary data
  _w.resize(n);
  cblas_dcopy(n, c.data(), c.inc(), _w.data(), _w.inc());
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 0.5, H.data(), H.leading_dim(), x, 1, 1.0, _w.data(), _w.inc());

  // call the objective function
  obj_value = cblas_ddot(n, _w.data(), _w.inc(), x, 1);

  return true;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // compute the gradient 
  cblas_dcopy(n, c.data(), c.inc(), grad_f, 1);
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, H.data(), H.leading_dim(), x, 1, 1.0, grad_f, 1);

  return true;
}

/// User should not call this method (for IPOPT)
bool NQP_IPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_ACT_CONTACTS = epd->N_ACT_CONTACTS;
  const unsigned N_LIMITS = epd->N_LIMITS;

  // copy x
  _x.resize(n);
  std::copy(x, x+n, _x.begin());

  // compute w
  if (_dense)
    R.mult(_x, _w) += z;
  else
    _w = _x;

  // evaluate the kappa constraint 
  Cn_block.mult(_w, _workv);
  _workv += Cn_v;
  g[N_CONTACTS+N_LIMITS] = epd->kappa - std::accumulate(_workv.begin(), _workv.end(), 0.0) + _tol;

  // evaluate the non-interpenetration constraints
  std::copy(_workv.begin(), _workv.end(), g);

  // evaluate the joint limit constraints
  L_block.mult(_w, _workv) -= epd->L_v;
  std::copy(_workv.begin(), _workv.end(), g+N_ACT_CONTACTS);

  // evaluate the Coulomb friction constraint
  for (unsigned i=0, j=N_CONTACTS+N_LIMITS+1; i< N_ACT_CONTACTS; i++, j++)
  {
    const unsigned N_IDX = i;
    const unsigned S_IDX = i+N_ACT_CONTACTS;
    const unsigned T_IDX = i+N_ACT_CONTACTS*2;
    const double muc_sq = mu_c[N_IDX];
    const double mu_visc_sq = mu_visc[N_IDX];
    g[j] = muc_sq * sqr(_w[N_IDX]) + mu_visc_sq - sqr(_w[S_IDX]) - sqr(_w[T_IDX]) + _tol;
  }  

  return true;
}

bool NQP_IPOPT::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values)
{
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_ACT_CONTACTS = epd->N_ACT_CONTACTS;

  // only do these computations if 'values' is non-null
  if (values)
  {
    // setup invariant part of Jacobian 
    std::copy(_cJac.get(), _cJac.get() + _cJac_constant, values);

    // setup gradient of Coulomb friction constraint
    if (!_dense)
    {
      for (unsigned i=0, j=_cJac_constant; i< N_ACT_CONTACTS; i++, j+= 3)
      {
        const unsigned N_IDX = i;
        const unsigned S_IDX = i+N_ACT_CONTACTS;
        const unsigned T_IDX = i+N_ACT_CONTACTS*2;
        const double muc_sq = mu_c[N_IDX];
        const double mu_visc_sq = mu_visc[N_IDX];
        values[j+0] = 2.0*muc_sq*x[N_IDX];
        values[j+1] = -2.0*x[S_IDX];
        values[j+2] = -2.0*x[T_IDX];
      }
    }
    else
    {
      // copy x
      _x.resize(n);
      std::copy(x, x+n, _x.begin());

      // compute w
      R.mult(_x, _w) += z;

      for (unsigned i=0; i< N_ACT_CONTACTS; i++)
      {
        // get the contact event index
        const unsigned N_IDX = 0;
        const unsigned S_IDX = i + N_ACT_CONTACTS; 
        const unsigned T_IDX = i + N_ACT_CONTACTS*2;

        // compute contact friction
        R.get_row(N_IDX, _workv2) *= ((double) 2.0 * _w[N_IDX] * mu_c[N_IDX]);
        R.get_row(S_IDX, _workv) *= ((double) 2.0 * _w[S_IDX]);
        _workv2 -= _workv;
        R.get_row(T_IDX, _workv) *= ((double) 2.0 * _w[T_IDX]); 
        _workv2 -= _workv;
        std::copy(_workv2.begin(), _workv2.end(), values+_cJac_constant+n*i);
      }
    }
  }

  // setup the indices
  if (iRow && jCol)
  {
    std::copy(_cJac_iRow.get(), _cJac_iRow.get()+nele_jac, iRow);
    std::copy(_cJac_jCol.get(), _cJac_jCol.get()+nele_jac, jCol);
  }

  return true;
}

// sets the final solution 
void NQP_IPOPT::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index, const Number* g, const Number* lam, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ipcq)
{
  // copy x
  _x.resize(n);
  std::copy(x, x+n, _x.begin());

  // account for nullspace
  if (_dense)
  {
    R.mult(_x, _w) += z;
    z = _w;
  }
  else
    z = _x;
}

/// The Hessian
bool NQP_IPOPT::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lam, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_ACT_CONTACTS = epd->N_ACT_CONTACTS;
  const unsigned N_LIMITS = epd->N_LIMITS;

  // NOTE: lambda is of size N_CONTACTS + N_ACT_CONTACTS + N_LIMITS + 1

  // get pointer values
  double* hobj = _h_obj.get();

  // setup Hessian, if desired
  if (values)
  {
    // zero the values
    std::fill_n(values, nele_hess, 0.0);

    // scale the objective function part of the Hessian
    cblas_dcopy(_nnz_h_obj, hobj, 1, values, 1); 
    cblas_dscal(_nnz_h_obj, obj_factor, values, 1); 

    // get the starting index for lambda
    unsigned lambda_start = N_CONTACTS + N_LIMITS + 1;
    for (unsigned k=0; k< N_ACT_CONTACTS; k++)
    {
      // scale the values in the constraint Hessian
      _workv.set_zero(_nnz_h_con[k]);
      cblas_dcopy(_nnz_h_con[k], _h_con[k].get(), 1, _workv.data(), _workv.inc());
      cblas_dscal(_nnz_h_con[k], lam[lambda_start+k], _workv.data(), _workv.inc());

      // get the workv iterator
      ColumnIteratord_const workv_iter = _workv.column_iterator_begin();

      // update the combined Hessian
      for (unsigned i=0; i< _h_con_indices[k].size(); i++)
        values[_h_con_indices[k][i]] += workv_iter[i];
    }
  }      

  // setup indices
  if (iRow)
  {
    std::copy(_h_iRow.get(), _h_iRow.get()+nele_hess, iRow);
    std::copy(_h_jCol.get(), _h_jCol.get()+nele_hess, jCol);
  }

  return true;
}


