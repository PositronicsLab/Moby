/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <numeric>
#include <Moby/LCP_IPOPT.h>

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

LCP_IPOPT::LCP_IPOPT()
{
  _ZERO_EPS = NEAR_ZERO;
}

/// User should not call this method (for IPOPT)
bool LCP_IPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)   
{
  n = q.rows();
  m = q.rows();

  // compute nonzeros for objective part of Hessian 
  nnz_h_lag = 0;

  // compute nonzeros for objective part of Hessian
  for (unsigned i=0, k=0; i< n; i++)
    for (unsigned j=i; j< n; j++)
      if (std::fabs(M(i,j)) > _ZERO_EPS)
        nnz_h_lag++;

  // setup objective component of Hessian and indices
  _h_obj = shared_array<double>(new double[nnz_h_lag]);
  _h_iRow = shared_array<unsigned>(new unsigned[nnz_h_lag]);  
  _h_jCol = shared_array<unsigned>(new unsigned[nnz_h_lag]);
  for (unsigned i=0, k=0; i< n; i++)
    for (unsigned j=i; j< n; j++)
      if (std::fabs(M(i,j)) > _ZERO_EPS)
      {
        _h_obj[k] = M(i,j);
        _h_iRow[k] = i;
        _h_jCol[k] = j;
        k++;
      }

  // initialize number of zeros in the constraint 
  nnz_jac_g = nnz_h_lag;
  
  // C-indexing style
  index_style = TNLP::C_STYLE;

  return true;
}

/// User should not call this method (for IPOPT)
bool LCP_IPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
  const double INF = std::numeric_limits<double>::max();

  // get info
  const unsigned N_CONTACTS = epd->N_CONTACTS;
  const unsigned N_LIMITS = epd->N_LIMITS;

  // set variable bounds
  for (unsigned i=0; i< n; i++)
  {
    x_l[i] = 0.0;
    x_u[i] = INF;
  }

  // set bounds on inequality constraints 
  for (unsigned i=0; i< N_CONTACTS*2+N_LIMITS+1; i++)
  {
    g_l[i] = -q[i];
    g_u[i] = INF;
  }

  return true;
}

/// User should not call this method (for IPOPT)
bool LCP_IPOPT::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lam)
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
bool LCP_IPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // get temporary data
  _workv.resize(n);
  cblas_dcopy(n, q.data(), q.inc(), _workv.data(), _workv.inc());
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, M.data(), M.leading_dim(), x, 1, 1.0, _workv.data(), _workv.inc());

  // call the objective function
  obj_value = cblas_ddot(n, _workv.data(), _workv.inc(), x, 1);

  return true;
}

/// User should not call this method (for IPOPT)
bool LCP_IPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // compute the gradient 
  cblas_dcopy(n, q.data(), q.inc(), grad_f, 1);
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 2.0, M.data(), M.leading_dim(), x, 1, 1.0, grad_f, 1);

  return true;
}

/// User should not call this method (for IPOPT)
bool LCP_IPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // copy x
  z.resize(n);
  std::copy(x, x+n, z.begin());

  // evaluate M*x >= -q 
  M.mult(z, _workv);

  // evaluate the non-interpenetration constraints
  std::copy(_workv.begin(), _workv.end(), g);

  return true;
}

bool LCP_IPOPT::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values)
{
  const unsigned N_CONTACTS = epd->N_CONTACTS;

  // only do these computations if 'values' is non-null
  if (values)
    std::copy(_h_obj.get(), _h_obj.get() + nele_jac, values);

  // setup the indices
  if (iRow && jCol)
  {
    std::copy(_h_iRow.get(), _h_iRow.get()+nele_jac, iRow);
    std::copy(_h_jCol.get(), _h_jCol.get()+nele_jac, jCol);
  }

  return true;
}

// sets the final solution 
void LCP_IPOPT::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index, const Number* g, const Number* lam, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ipcq)
{
  // copy x
  z.resize(n);
  std::copy(x, x+n, z.begin());
}

/// The Hessian
bool LCP_IPOPT::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lam, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
  // get pointer values
  double* hobj = _h_obj.get();

  // setup Hessian, if desired
  if (values)
  {
    // zero the values
    std::fill_n(values, nele_hess, 0.0);

    // scale the objective function part of the Hessian
    cblas_dcopy(nele_hess, hobj, 1, values, 1); 
    cblas_dscal(nele_hess, obj_factor, values, 1); 
  }      

  // setup indices
  if (iRow)
  {
    std::copy(_h_iRow.get(), _h_iRow.get()+nele_hess, iRow);
    std::copy(_h_jCol.get(), _h_jCol.get()+nele_hess, jCol);
  }

  return true;
}


