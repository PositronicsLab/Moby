/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_OPT_H
#define _MOBY_OPT_H

#include <Moby/Constants.h>
#include <Moby/MatrixNN.h>
#include <Moby/VectorN.h>

namespace Moby {

/// Structure for performing linear programming
struct LPParams
{
  LPParams() {}
  LPParams(const VectorN& c, const MatrixN& A, const VectorN& b, const MatrixN& M, const VectorN& q)
  {
    n = c.size();
    this->c = c;
    this->A = A;
    this->b = b;
    this->M = M;
    this->q = q;
  }

  /// The number of variables in the problem
  unsigned n;
  
  /// The objective function coefficients
  VectorN c;

  /// A the coefficients of the equality constraint matrix (i.e., Ax = b after the linear program is solved)
  MatrixN A;

  /// b the right hand side of the equality constraints (i.e., Ax = b after the linear program is solved)
  VectorN b;

  /// M the coefficients of the inequality constraint matrix (i.e., Mx >= q after the linear program is solved)
  MatrixN M;

  /// q the right hand side of the inequality constraints (i.e., Mx >= q after the linear program is solved)
  VectorN q;
};

/// Structure for performing unconstrained optimization using BFGS
struct BFGSParams
{
  BFGSParams()
  {
    n = 0;
    fx = 0;
    grad = NULL;
    eps = NEAR_ZERO;
    max_iterations = std::numeric_limits<unsigned>::max();
    max_fn_evals = std::numeric_limits<unsigned>::max();
    min_alpha_step = std::numeric_limits<Real>::epsilon();
    max_alpha_step = (Real) 10.0;
    f_tol = NEAR_ZERO;
    grad_tol = NEAR_ZERO;
    x_rtol = NEAR_ZERO;
    data = NULL;
    tcheck = NULL;
    iterations = 0;
  }

  /// The number of iterations performed by BFGS (on return from BFGS())
  unsigned iterations;

  /// The dimensionality of the optimization problem
  unsigned n;

  /// The tolerance to which the gradient may be zero to continue
  Real eps;

  /// The maximum number of BFGS iterations
  unsigned max_iterations;

  /// The maximum number of function evaluations per BFGS iteration
  unsigned max_fn_evals;

  /// The minimum step size for determining alpha (machine epsilon by default)
  Real min_alpha_step;

  /// The maximum step size for determining alpha (10.0 by default)
  Real max_alpha_step;

  /// The sufficient decrease condition (decreases below this number lead to termination) for determining alpha
  Real f_tol;

  /// The directional derivative condition (norms of gradients below this number lead to termination) for determining alpha
  Real grad_tol;

  /// Minimum relative width of the interval of uncertainty (widths below this number lead to termination) for determining alpha
  Real x_rtol;

  /// The objective function to be minimized
  Real (*fx)(const VectorN&, void*);

  /// The gradient of the objective function
  void (*grad)(const VectorN&, void*, VectorN&);

  /// Function for checking for early termination (should return <b>true</b> if/when early termination is desired)
  bool (*tcheck)(const VectorN&, void*);

  /// Data passed to fx and grad
  void* data;
};

/// Structure for performing convex optimization
struct OptParams
{
  OptParams()
  {
    iterations = 0;
    n = m = r = 0;
    fx = NULL;
    grad = NULL;
    hess = NULL;
    tcheck = NULL;
    max_iterations = std::numeric_limits<unsigned>::max();
    mu = 10.0;
    alpha = 0.05;
    beta = 0.5;
    eps = NEAR_ZERO;
    zero_tol = NEAR_ZERO;
    lb.resize(0);
    ub.resize(0);
  }

  OptParams(unsigned n, unsigned m, unsigned r, Real (*fx)(const VectorN&, unsigned, void*), void (*grad)(const VectorN&, unsigned, VectorN&, void*), bool (*hess)(const VectorN&, unsigned, MatrixNN&, void*))
  {
    data = NULL;

    // set number of iterations
    this->iterations = 0;

    // store variable size
    this->n = n;

    // store number of nonlinear inequality constraints
    this->m = m;

    // store number of nonlinear equality constraints
    this->r = r;

    // setup maximum number of iterations
    max_iterations = std::numeric_limits<unsigned>::max();

    // setup function pointers
    this->fx = fx;
    this->grad = grad;
    this->hess = hess;
    tcheck = NULL;

    // make A, b, M, q zero by default
    A.resize(0, n);
    b.resize(0);
    M.resize(0, n);
    q.resize(0);

    // no bounds on the variables
    lb.resize(0);
    ub.resize(0);

    // setup defaults
    mu = 10.0;
    alpha = 0.05;
    beta = 0.5;
    eps = NEAR_ZERO;
    zero_tol = NEAR_ZERO;
  }

  OptParams(const OptParams& c) { operator=(c); }
  
  OptParams& operator=(const OptParams& c)
  {
    iterations = c.iterations;
    max_iterations = c.max_iterations;
    n = c.n;
    m = c.m;
    r = c.r;
    A = c.A;
    b = c.b;
    q = c.q;
    M = c.M;
    mu = c.mu;
    lb = c.lb;
    ub = c.ub;
    alpha = c.alpha;
    beta = c.beta;
    zero_tol = c.zero_tol;
    eps = c.eps;
    data = c.data;
    fx = c.fx;
    grad = c.grad;
    hess = c.hess;
    tcheck = c.tcheck;
    return *this;
  }

  /// Number of Newton/Quasi-Newton iterations performed by convex optimization method (on return)
  unsigned iterations;

  /// maximum number of iterations
  unsigned max_iterations;

  /// dimensionality of problem
  unsigned n;

  /// number of nonlinear inequality constraints
  unsigned m;

  /// number of nonlinear equality constraints
  unsigned r;

  /// Matrix A used in linear equality constraints Ax = b
  MatrixN A;

  /// Vector b used in linear equality constraints Ax = b
  VectorN b;

  /// Matrix M used in linear inequality constraints Mx >= q
  MatrixN M;

  /// Vector q used in linear inequality constraints Mx >= q
  VectorN q;

  /// Lower bounds on x
  VectorN lb;

  /// Upper bounds on x
  VectorN ub;

  /// mu parameter for increasing t for interior-point method (mu > 1), default = 10.0
  Real mu;

  /// Alpha parameter for backtracking line search (0 < alpha < 1), default = 0.05
  Real alpha;

  /// Beta parameter for backtracking line search (0 < beta < 1), default = 0.5
  Real beta;

  /// Tolerance for checking for zero
  Real zero_tol;

  /// interior point method terminates only after duality gap < eps
  Real eps;

  /// primal-dual optimization terminates only after surrogate duality gap < eps and primal-dual feasibility < eps_feas 
  Real eps_feas;

  /// optional data passed to fx(), grad(), and hess()
  void* data;

  /// fx(x,f,data) pointer to a function that evaluates the objective (f[0]) and inequality constraint functions (f[1..m]) 
  Real (*fx)(const VectorN&, unsigned, void*);

  /// grad(x,g,data) pointer to a function for computing the gradients of the objective (g[0]) and inequality constraint functions (g[1..m])
  void (*grad)(const VectorN&, unsigned, VectorN&, void*);

  /// hess(x,h,data) pointer to a function for computing the Hessians of the objective (h[0]) and inequality constraint functions (h[1..m]) -- unsigned parameter is the one-index of the constraint (index 0 corresponds to the objective); returns <b>true</b> if the Hessian is non-zero, <b>false</b> otherwise
  bool (*hess)(const VectorN&, unsigned, MatrixNN&, void*);

  /// tcheck(x,data) pointer to a function for determining whether or not the optimization should terminate at the given x
  bool (*tcheck)(const VectorN&, void*);
};

/// Routines for performing mathematical optimization 
class Optimization 
{
  public:
    static void sqp(OptParams& oparams, VectorN& x);
    static VectorN& find_cauchy_point(const MatrixNN& G, const VectorN& c, const VectorN& l, const VectorN& u, const VectorN& gradient, VectorN& x);
    static unsigned qp_gradproj(const MatrixNN& G, const VectorN& c, const VectorN& l, const VectorN& u, unsigned max_iter, VectorN& x, Real tol);
    static bool brent(Real x_lower, Real x_upper, Real& x, Real& fx, Real (*f)(Real, void*), Real tol, void* params);
    static bool optimize_convex(OptParams& cparams, VectorN& x);
    static bool make_feasible_convex(OptParams& cparams, VectorN& x, VectorN (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const OptParams&) = NULL);
    static bool make_feasible_convex2(OptParams& cparams, VectorN& x, VectorN (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const OptParams&) = NULL);
    static bool optimize_convex_pd(OptParams& cparams, VectorN& x, VectorN (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const OptParams&) = NULL);
    static bool mlcp(VectorN& y, VectorN& z, const MatrixNN& M, const VectorN& q, bool (*lcp_solver)(const MatrixNN&, const VectorN&, VectorN&) = NULL);
    static Real log_vec(const std::vector<Real>& v, unsigned st);
    static Real log_fn(Real (*f)(const VectorN&, unsigned, void*), unsigned m, const VectorN& x, void* data);
    static bool lp(const MatrixN& A, const VectorN& b, const VectorN& c, const VectorN& l, const VectorN& u, VectorN& x);
    static bool qp_convex_ip(const MatrixNN& G, const VectorN& c, OptParams& oparams, VectorN& x);
    static bool qp_convex_activeset_infeas_tcheck(const VectorN& x, void* data);
    static void qp_convex_activeset_infeas(const MatrixNN& G, const VectorN& c, Real upsilon, OptParams& qparams, VectorN& x, bool hot_start = false);
    static void qp_convex_activeset(const MatrixNN& G, const VectorN& c, OptParams& qparams, VectorN& x, bool hot_start = false);
    static void qp_to_lcp1(const MatrixNN& G, const VectorN& c, const MatrixN& M, const VectorN& q, MatrixNN& MM, VectorN& qq);
    static void qp_to_lcp2(const MatrixNN& G, const VectorN& c, const MatrixN& M, const VectorN& q, MatrixNN& MM, VectorN& qq);
    static bool lcp_gradproj(const MatrixNN& M, const VectorN& q, VectorN& x, Real tol, unsigned max_iterations = std::numeric_limits<unsigned>::max());
    static bool polyak(const MatrixNN& A, const VectorN& b, const VectorN& c, const VectorN& d, unsigned max_iter, VectorN& x);
    static bool make_feasible_qp(const MatrixN& A, const VectorN& b, const MatrixN& M, const VectorN& q, VectorN& x, Real tol = NEAR_ZERO);
    static bool lp_simplex(const LPParams& lpparams, VectorN& x);
    static void lcp_enum(const MatrixNN& M, const VectorN& q, std::vector<VectorN>& z);
    static bool lcp_lemke(const MatrixNN& M, const VectorN& q, VectorN& z, Real piv_tol = NEAR_ZERO, Real zero_tol = NEAR_ZERO);
    static bool lcp_lemke_regularized(const MatrixNN& M, const VectorN& q, VectorN& z, int min_exp = -20, unsigned step_exp = 4, int max_exp = 20, Real piv_tol = NEAR_ZERO, Real zero_tol = NEAR_ZERO);
    static bool lcp_convex_ip(const MatrixNN& M, const VectorN& q, VectorN& z, Real tol=NEAR_ZERO, Real eps=NEAR_ZERO, Real eps_feas=NEAR_ZERO, unsigned max_iterations = std::numeric_limits<unsigned>::max());
    static bool lcp_iter_PD(const MatrixNN& M, const VectorN& q, VectorN& z, Real tol = NEAR_ZERO, const unsigned iter = std::numeric_limits<unsigned>::max());
    static bool lcp_iter_symm(const MatrixNN& M, const VectorN& q, VectorN& z, Real tol = NEAR_ZERO, const unsigned iter = std::numeric_limits<unsigned>::max());
    static void solve_CG(const MatrixNN& A, const VectorN& b, VectorN& x, Real err_tol = 1.0);
    static void eliminate_redundant_constraints(MatrixN& A, VectorN& b);
    static void eliminate_redundant_constraints(MatrixN& A);
    static Real search_line(const VectorN& x, Real (*fn)(const VectorN&, void*), void (*grad)(const VectorN&, void*, VectorN&), const VectorN& sdir, void* data = NULL, Real f_tol = (Real) 0.0, Real grad_tol = (Real) 0.0, Real x_rtol = (Real) 0.0, Real step_min = NEAR_ZERO, Real step_max = 1e6, unsigned max_fn_evals = 10000);
    static Real search_line_backtrack(const VectorN& x, Real (*fn)(const VectorN&, void*), void (*grad)(const VectorN&, void*, VectorN&), const VectorN& dx, Real alpha, Real beta, void* data);
    static VectorN& BFGS(BFGSParams& params, VectorN& x);
    static bool optimize_convex_BFGS(OptParams& cparams, VectorN& x);
    static bool make_feasible_convex_BFGS(OptParams& cparams, VectorN& x);
    static bool lcp_ip_BFGS(const MatrixNN& M, const VectorN& q, VectorN& z, Real eps=NEAR_ZERO, Real eps_feas=NEAR_ZERO);

  private:
    static void add_to_working_set(const VectorN& vec, Real bi, MatrixN& AMr, VectorN& br);
    static bool add_to_working_set(unsigned i, const std::vector<unsigned>& bounded, const std::vector<unsigned>& free, MatrixN& AMr, MatrixN& AMr_bounded, MatrixN& Z, MatrixN& Y);
    static bool add_to_working_set(const VectorN& vec, const std::vector<unsigned>& bounded, const std::vector<unsigned>& free, MatrixN& AMr, MatrixN& AMr_bounded, MatrixN& Z, MatrixN& Y);
    static void reform_working_set(const MatrixN& Ar, const MatrixN& M, const std::vector<bool>& working, const std::vector<unsigned>& bounded, const std::vector<unsigned>& free, MatrixN& AMr, MatrixN& AMr_bounded, MatrixN& Z, MatrixN& Y);
    static void determine_selection(const std::vector<bool>& bworking, std::vector<unsigned>& bounded, std::vector<unsigned>& free);
    static void condition_and_factor_PD(MatrixNN& H);
    static void condition_hessian(MatrixNN& H);
    static bool tcheck_cvx_opt_BFGS(const VectorN& x, void* data);
    static Real fx_cvx_opt_BFGS(const VectorN& x, void* data);
    static void grad_cvx_opt_BFGS(const VectorN& x, void* data, VectorN& g);
    static Real absmax(Real a, Real b, Real c);
    static int cstep(Real stmin, Real stmax, Real& stx, Real& fx, Real& dx, Real& sty, Real& fy, Real& dy, Real& stp, Real& fp, Real& dp, bool& brackt);
    static VectorN calc_residual(const VectorN& x, const VectorN& lambda, const VectorN& nu, Real t, const OptParams& cparams);
    static VectorN solve_KKT_pd(const VectorN& x, const VectorN& lambda, const VectorN& r, const OptParams& cparams);
    static VectorN solve_KKT_lcp_ip_feas(const VectorN& x, const VectorN& lambda, const VectorN& r, const OptParams& cparams);
    static bool make_feasible_tcheck(const VectorN& y, void* data);
    static Real make_feasible_f(const VectorN& y, unsigned idx, void* data);
    static void make_feasible_grad(const VectorN& y, unsigned idx, VectorN& g, void* data);
    static bool make_feasible_hess(const VectorN& y, unsigned, MatrixNN& H, void* data);
    static bool feasible(Real (*f)(const VectorN&, unsigned, void*), unsigned m, const VectorN& x, Real infeas_tol, void* data);
    static bool feasible(const std::vector<Real>& f);
    static bool feasible(const VectorN& f);
    static VectorN ngrad(const VectorN& x, Real t, Real h, void* data, Real (*ofn)(const VectorN&, Real, void*));
    static MatrixNN nhess(const VectorN& x, Real t, Real h, void* data, Real (*ofn)(const VectorN&, Real, void*));
    static Real finitize(Real x);
    static VectorN remove_component(const VectorN& v, unsigned k);
    static VectorN insert_component(const VectorN& v, unsigned k);
    static bool lcp_ip_tcheck(const VectorN& x, void* data);
    static Real lcp_ip_fx(const VectorN& x, unsigned idx, void* data);
    static void lcp_ip_grad(const VectorN& x, unsigned idx, VectorN& g, void* data);
    static bool lcp_ip_hess(const VectorN& x, unsigned idx, MatrixNN& H, void* data);
    static Real qp_ip_fx(const VectorN& x, unsigned idx, void* data);
    static void qp_ip_grad(const VectorN& x, unsigned idx, VectorN& g, void* data);
    static bool qp_ip_hess(const VectorN& x, unsigned idx, MatrixNN& H, void* data);
    static void setup_C(OptParams& oparams, const VectorN& x, VectorN& C);
    static void setup_A(OptParams& oparams, const VectorN& x, MatrixN& A);

    /// Computes an integral power of a base
    template <class T>
    static T intpow(T base, unsigned pow)
    {
      T result = 1;
      for (unsigned i=pow; i > 0; i--)
        result *= base;
      return result;
    }
}; // end class

} // end namespace

#endif

