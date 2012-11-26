/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef USE_GLPK
extern "C" 
{ 
#include <glpk.h> 
}
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <numeric>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <iostream>
#include <set>
#include <boost/algorithm/minmax.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_array.hpp>
#include <Moby/cblas.h>
#include <Moby/select>
#include <Moby/LinAlg.h>
#include <Moby/FastThreadable.h>
#include <Moby/Log.h>
#include <Moby/SingularException.h>
#include <Moby/NonsquareMatrixException.h>
#include <Moby/NumericalException.h>
#include <Moby/Optimization.h>
#ifdef USE_PATH
#include <Moby/PathLCPSolver.h>
#endif

using namespace Moby;
using boost::shared_ptr;
using boost::shared_array;
using std::endl;
using std::cerr;
using std::vector;
using std::list;
using std::set;
using std::pair;

template <class T>
struct nlog : public std::unary_function<T, void>
{
  nlog() {}
  T operator()(T x) { return std::log(x); }
};

/// For SQP solver
void SLSQP(unsigned M, unsigned MEQ, unsigned LA, unsigned N, Real* X, const Real* XL, const Real* XU, Real F, Real* C, Real* G, Real* A, Real ACC, unsigned& ITER, int& MODE, Real* W, unsigned L_W, int* JW, unsigned L_JW, shared_ptr<void>& state);

/// Structure for computing QP data using primal/dual method
struct QPData
{
  MatrixN G;
  VectorN c;
  void* data;
};

/// Structure for making infeasible points feasible (for convex optimization)
struct FeasibilityData
{
  unsigned m;          // number of constraint functions in the original problem
  unsigned n;          // size of y in the original problem
  OptParams* oparams;  // original optimization parameters
  void* data;          // data in the original problem

  // functions in the original problem
  Real (*f0)(const VectorN&, void*);
  void (*fx)(const VectorN&, VectorN&, void*);
  void (*gx)(const VectorN&, VectorN&, void*);
  void (*grad0)(const VectorN&, VectorN&, void*);
  void (*cJac_f)(const VectorN&, MatrixN&, void*);
  void (*cJac_g)(const VectorN&, MatrixN&, void*);
  void (*hess)(const VectorN&, Real objscal, const VectorN& lambda, const VectorN& nu, MatrixN&, void*);
};

/// The signum function
static Real sign(Real x)
{
  if (x > (Real) 0.0)
    return (Real) 1.0;
  else if (x < (Real) 0.0)
    return (Real) -1.0;
  else
    return (Real) 0.0;
}

/// Sets up the C vector for SQP
void Optimization::setup_C(OptParams& oparams, const VectorN& x, VectorN& C)
{
  SAFESTATIC VectorN workv;
  SAFESTATIC vector<Real> ff, fg;

  // evaluate all nonlinear constraints
  evaluate_inequality_constraints(x, oparams, ff);
  evaluate_equality_constraints(x, oparams, fg);

  // resize C just in case
  C.resize(ff.size() + fg.size());


  // transform inequality constraints to f(x) >= 0 form
  std::transform(ff.begin(), ff.end(), C.begin()+fg.size(), std::negate<Real>()); 

  // copy to C
  std::copy(fg.begin(), fg.end(), C.begin());
}

/// Sets up the A matrix for SQP
void Optimization::setup_A(OptParams& oparams, const VectorN& x, MatrixN& A)
{
  SAFESTATIC VectorN workv;
  SAFESTATIC MatrixN Jf, Jg, sub;

  const unsigned N_LININEQ = oparams.q.size();
  const unsigned N_LINEQ = oparams.b.size();

  // setup the Jacobians 
  calc_Df(x, oparams, Jf);
  calc_Dg(x, oparams, Jg); 

  // reverse nonlinear inequalities so that they obey f(x) >= 0
  Jf.negate();

  // setup A
  A.set_sub_mat(0,0,Jg);
  A.set_sub_mat(Jg.rows(),0,Jf);
}

/// Does sequential quadratic programming
void Optimization::sqp(OptParams& oparams, VectorN& x)
{
  VectorN C, G, W, xlb, xub;
  MatrixN A;
  vector<int> JW;
  const Real INF = std::numeric_limits<Real>::max();

  // if there are no variables, quit now
  const unsigned N = oparams.n;
  if (N == 0)
    return;

  // setup data for slsqp
  const unsigned MEQ = oparams.b.size() + oparams.r;
  const unsigned M = MEQ + oparams.q.size() + oparams.m;
  const unsigned LA = std::max(M, (unsigned) 1);
  const Real ACC = oparams.eps;
  unsigned ITER = oparams.max_iterations;
  const unsigned N1 = N+1;
  const unsigned MINEQ = M - MEQ + 2*N1;
  unsigned L_W = (3*N1+M)*(N1+1) + (N1-MEQ+1)*(MINEQ+2) + 2*MINEQ+ 
                  (N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1 + N1*N/2 + 2*M+3*N + 3*N1+1;
  unsigned L_JW = MINEQ;

  // setup lower bounds on x
  if (oparams.lb.size() > 0)
    xlb.copy_from(oparams.lb);
  else
    xlb.set_one(N) *= (Real) -INF;

  // setup upper bounds on x
  if (oparams.ub.size() > 0)
    xub.copy_from(oparams.ub);
  else
    xub.set_one(N) *= (Real) INF;

  // setup x
  x.resize(N);

  // setup work vectors / matrices
  C.resize(M);
  G.resize(N);
  A.resize(LA,N1);
  W.resize(L_W);
  JW.resize(L_JW);

restart:
  // begin looping
  int MODE = 0;

  // do initial evaluations
  Real F = oparams.f0(x, oparams.data);
  oparams.grad0(x, G, oparams.data);
  setup_C(oparams, x, C);
  setup_A(oparams, x, A);

  // setup state pointer
  shared_ptr<void> state;
 
  do
  {
    // check for gradient evaluation
    if (MODE < 0)
    {
      // compute gradient and A matrix
      oparams.grad0(x, G, oparams.data);
      setup_A(oparams, x, A);
    }
    else if (MODE == 1)
    {
      // do function evaluations
      F = oparams.f0(x, oparams.data);
      setup_C(oparams, x, C);
    }

    // call the method
    SLSQP(M, MEQ, LA, N, x.begin(), xlb.begin(), xub.begin(), F, C.begin(), G.begin(), A.begin(), ACC, ITER, MODE, W.begin(), L_W, &JW.front(), L_JW, state);
  }
  while (MODE < 0 || MODE == 1);

  // check for error condition
  if (MODE > 1)
  {
    switch (MODE)
    {
      case 2:  
        throw std::runtime_error("Optimization::sqp()- too many runtime constraints");

      case 3:  
        FILE_LOG(LOG_OPT) << "Optimization::sqp()- too many least squares iterations" << endl;
        break;

      case 4:  
        FILE_LOG(LOG_OPT) << "Optimization::sqp()- inequality constraints incompatible" << endl;
        break;

      case 5:  
      case 6:  
        FILE_LOG(LOG_OPT) << "Optimization::sqp()- singular matrix in least square subproblem" << endl;
        break;

      case 7:  
        FILE_LOG(LOG_OPT) << "Optimization::sqp()- rank-deficient equality constraint subproblem" << endl;
        break;

      case 8:  
        FILE_LOG(LOG_OPT) << "Optimization::sqp()- positive directional derivative for line search" << endl;
        break;

      case 9:
        FILE_LOG(LOG_OPT) << "Optimization::sqp() - iteration max reached" << endl;
        break;
    }
  }

  // save iteration count
  oparams.iterations = ITER;

  // check for restart
  if (MODE >= 10)
  {
    L_W = MODE/1000;
    L_JW = MODE-1000*L_W;
    W.resize(L_W);
    JW.resize(L_JW);
    goto restart;
  }
}

/// Evaluates all equality constraints (nonlinear and linear)
void Optimization::evaluate_equality_constraints(const VectorN& x, const OptParams& cparams, vector<Real>& fc)
{
  SAFESTATIC VectorN workv;

  // determine constraint sizes
  const unsigned r = cparams.r;
  const unsigned N_LINEQ = cparams.b.size();
  fc.resize(r+N_LINEQ);

  // evaluate nonlinear inequalities first
  workv.resize(r);
  if (cparams.r > 0 && cparams.gx)
    cparams.gx(x, workv, cparams.data);
  std::copy(workv.begin(), workv.end(), fc.begin());

  // evaluate nonlinear constraints
  cparams.A.mult(x, workv) -= cparams.b;
  std::copy(workv.begin(), workv.end(), fc.begin()+r);
}

/// Evaluates all inequality constraints (nonlinear, bounds, and linear)
void Optimization::evaluate_inequality_constraints(const VectorN& x, const OptParams& cparams, vector<Real>& fc)
{
  SAFESTATIC VectorN workv;

  // determine constraint sizes
  const unsigned m = cparams.m;
  const unsigned N_LB = cparams.lb.size();
  const unsigned N_UB = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();
  fc.resize(m+N_LB+N_UB+N_LININEQ);

  // evaluate nonlinear inequalities first
  workv.resize(m);
  if (cparams.m > 0 && cparams.fx)
    cparams.fx(x, workv, cparams.data);
  std::copy(workv.begin(), workv.end(), fc.begin());

  // evaluate lower bounds constraints
  for (unsigned i=0, j=m; i< N_LB; i++)
    fc[j++] = cparams.lb[i] - x[i];

  // evaluate upper bounds constraints
  for (unsigned i=0, j=m+N_LB; i< N_UB; i++)
    fc[j++] = x[i] - cparams.ub[i];

  // evaluate nonlinear constraints
  cparams.M.mult(x, workv) -= cparams.q;
  workv.negate();
  std::copy(workv.begin(), workv.end(), fc.begin()+m+N_LB+N_UB);
}

/// Makes a point feasible for a convex optimization problem using BFGS
/**
 * Operates by solving the convex optimization problem:
 * minimize s
 *   subject to Ax = b
 *              Gx + s >= m
 * using the Quasi-Newton BFGS algorithm.
 */
bool Optimization::make_feasible_convex_BFGS(OptParams& cparams, VectorN& x)
{
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex_BFGS() entered" << endl;

  // setup the feasibility data
  const unsigned n = x.size();
  const unsigned m = cparams.m;
  FeasibilityData fdata;
  fdata.oparams = &cparams;
  fdata.n = n;
  fdata.m = m;
  fdata.f0 = cparams.f0;
  fdata.fx = cparams.fx;
  fdata.gx = cparams.gx;
  fdata.grad0 = cparams.grad0;
  fdata.cJac_f = cparams.cJac_f;
  fdata.cJac_g = cparams.cJac_g;
  fdata.hess = NULL; 
  fdata.data = cparams.data;

  // determine new constraints
  const unsigned N_LB = cparams.lb.size();
  const unsigned N_UB = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;

  // evaluate f at current x
  vector<Real> f;
  evaluate_inequality_constraints(x, cparams, f);

  // verify that x is not already feasible (setup v simultaneously)
  if (f.empty())
  {
    FILE_LOG(LOG_OPT) << " -- no inequality constraints" << endl;
    FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex_BFGS() exiting!" << endl;
    return true;
  }

  Real v = *std::max_element(f.begin(), f.end());
  if (v < (Real) 0.0)
  {
    FILE_LOG(LOG_OPT) << " -- point is already feasible (v=" << v << ")!" << endl;
    FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex_BFGS() exiting!" << endl;
    return true;
  }

  // setup y
  VectorN y(n+1);
  y.set_sub_vec(0, x);
  y[n] = v + 1.0; 
  FILE_LOG(LOG_OPT) << "initial v: " << v << endl;

  // setup convex optimization params
  OptParams cp;
  cp = cparams;
  cp.n++;
  cp.m = mm;
  cp.f0 = &make_feasible_f0;
  cp.fx = &make_feasible_fx;
  cp.grad0 = &make_feasible_grad0;
  cp.cJac_f = &make_feasible_cJac_f;
  cp.hess = NULL;
  cp.tcheck = &make_feasible_tcheck;
  cp.data = &fdata;
  cp.max_iterations = std::max(1000 + (unsigned) std::log(cp.n), cparams.max_iterations);

  // do convex optimization with primal/dual method
  optimize_convex_BFGS(cp, y);  

  // get x out
  x = y.get_sub_vec(0, n);

  // determine s
  evaluate_inequality_constraints(x, cparams, f);
  v = *std::max_element(f.begin()+1, f.end()); 
  bool result = v < NEAR_ZERO;

  FILE_LOG(LOG_OPT) << "   x: " << x << endl;
  FILE_LOG(LOG_OPT) << "   v: " << v << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex_BFGS() successful? " << result << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex_BFGS() exiting" << endl;

  return result;
}

/// Termination checking function for convex optimization using BFGS
bool Optimization::tcheck_cvx_opt_BFGS(const VectorN& x, void* data)
{
  // get the data
 pair<Real, OptParams>& cvx_opt_BFGS_params = *((pair<Real, OptParams>*) data); 

  // get the original termination check function, if any
  return cvx_opt_BFGS_params.second.tcheck && cvx_opt_BFGS_params.second.tcheck(x, cvx_opt_BFGS_params.second.data);
}

/// Objective function for convex optimization using BFGS
Real Optimization::f_cvx_opt_BFGS(const VectorN& x, void* data)
{
  SAFESTATIC vector<Real> fc;

  // get the data
  pair<Real, OptParams>& cvx_opt_BFGS_params = *((pair<Real, OptParams>*) data); 
  OptParams& cparams = cvx_opt_BFGS_params.second;

  // get t
  Real t = cvx_opt_BFGS_params.first;

  // evaluate all constraints
  evaluate_inequality_constraints(x, cparams, fc);

  // compute the objective function
  Real value = t * cparams.f0(x, cparams.data);

  // add in the logarithmically transformed constraint functions
  for (unsigned i=0; i< fc.size(); i++)
    value -= std::log(-fc[i]);

  if (std::isnan(value))
    value = std::numeric_limits<Real>::max(); 

  return value;
}

/// Gradient function for convex optimization using BFGS
void Optimization::grad_cvx_opt_BFGS(const VectorN& x, void* data, VectorN& g)
{
  SAFESTATIC VectorN workv; 
  SAFESTATIC MatrixN cJac;
  SAFESTATIC vector<Real> fc;

  // get the data
  pair<Real, OptParams>& cvx_opt_BFGS_params = *((pair<Real, OptParams>*) data); 
  OptParams& cparams = cvx_opt_BFGS_params.second;

  // setup constraint sizes
  const unsigned N_LB = cparams.lb.size();
  const unsigned N_UB = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();

  // get t
  Real t = cvx_opt_BFGS_params.first;

  // evaluate all inequality constraints
  evaluate_inequality_constraints(x, cparams, fc);

  // call the original gradient function
  cparams.grad0(x, g, cvx_opt_BFGS_params.second.data);
  g *= t;

  // get the constraint Jacobian
  cJac.resize(cparams.m, cparams.n);
  if (cparams.cJac_f)
    cparams.cJac_f(x, cJac, cparams.data);

  // add in the gradients of the nonlinear inequality constraint functions
  for (unsigned i=0; i< cparams.m; i++)
  {
    cJac.get_row(i, workv) /= fc[i];
    g -= workv;
  }

  // NOTE: we are accounting for future presence of r here, though it is not
  //       being utilized
  // add in lower bounds constraints
  for (unsigned i=0, j=cparams.m+cparams.r; i< N_LB; i++)
    g[i] += (Real) 1.0/fc[j++];

  // add in upper bounds constraints
  for (unsigned i=0, j=cparams.m+cparams.r+N_LB; i< N_UB; i++)
    g[i] -= (Real) 1.0/fc[j++];

  // add in linear inequality constraints
  for (unsigned i=0, j=cparams.m+cparams.r+N_LB+N_UB; i< N_LININEQ; i++)
  {
    cparams.M.get_row(i, workv);
    workv /= fc[j++];
    g += workv;
  }
}

/// Performs convex optimization using the logarithmic barrier method with BFGS
/**
 * Uses the Quasi-Newton BFGS algorithm for the inner iteration.
 * \param cparams parameters for the optimization
 * \param x the feasible value from which to start; contains the optimal value
 *        on return
 * \note does not currently handle equality constraints
 * \return <b>true</b> if optimization successful; <b>false</b> otherwise
 */
bool Optimization::optimize_convex_BFGS(OptParams& cparams, VectorN& x)
{
  const Real T_INIT = 1e-3;
  SAFESTATIC vector<Real> fx;

  FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_BFGS() entered" << endl;

  // init t
  Real t = T_INIT;

  // determine n and m
  const unsigned n = x.size();

  if (cparams.gx || cparams.cJac_g || cparams.r > 0 || cparams.b.size() > 0)
    throw std::runtime_error("Equality constraints not currently supported for BFGS!");

  // verify that constraint functions are met
  evaluate_inequality_constraints(x, cparams, fx);
  vector<Real>::const_iterator maxele = std::max_element(fx.begin(), fx.end());
  if (maxele != fx.end() && *maxele >= (Real) 0.0)
  {
    FILE_LOG(LOG_OPT) << "initial point is not feasible!" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_BFGS() exited" << endl;
    return false;
  }

  // check whether we can terminate early
  if (cparams.tcheck && cparams.tcheck(x, cparams.data))
  {
    FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_BFGS() exited" << endl;
    return true;
  }

  // setup data to be passed to BFGS
  pair<Real, OptParams> cvx_opt_BFGS_params(t, cparams);

  // setup BFGS parameters
  BFGSParams BFGS_params;
  BFGS_params.n = n;
  BFGS_params.eps = cparams.eps;
  BFGS_params.max_iterations = std::numeric_limits<unsigned>::max();
  BFGS_params.f = &f_cvx_opt_BFGS;
  BFGS_params.grad = &grad_cvx_opt_BFGS;
  BFGS_params.data = (void*) &cvx_opt_BFGS_params;
  BFGS_params.tcheck = &tcheck_cvx_opt_BFGS;

  // evaluate objective function at current value of x to obtain best f0
  Real f0_best = cparams.f0(x, cparams.data);
  VectorN x_best;
  x_best.copy_from(x);

  // setup iteration count
  cparams.iterations = 0;

  // do outer iteration
  while (true)
  {
    FILE_LOG(LOG_OPT) << "outer iteration begun, t = " << t << endl;

    // **************************************************
    // centering step and update
    // **************************************************

    // setup the maximum number of iterations
    BFGS_params.max_iterations = cparams.max_iterations - cparams.iterations;
    
    // set t in the data parameters
    cvx_opt_BFGS_params.first = t;

    // call BFGS
    BFGS(BFGS_params, x);

    // update and store best_x and best_f0
    Real f0 = cparams.f0(x, cparams.data);
    if (f0 < f0_best)
    {
      x_best.copy_from(x);
      f0_best = f0;
    }

    // update the iteration count
    cparams.iterations += BFGS_params.iterations;

    // check whether we can terminate early
    if (cparams.tcheck && cparams.tcheck(x, cparams.data))
    {
      FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;
      return true;
    }

    // check stopping condition
    if (n/t < cparams.eps)
      break;

    // increase t
    t *= cparams.mu;
  }

  // use best value of x obtained
  x = x_best;
  FILE_LOG(LOG_OPT) << "  optimum: " << x << endl;
  FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;

  return true;
}

/// The BFGS algorithm for unconstrained minimization
VectorN& Optimization::BFGS(BFGSParams& params, VectorN& x)
{
  VectorN p, delta_x, delta_g, gnew, Hy, g;
  MatrixN ssrho, ss1, ss, Hrys, inv_hess;
  const Real BETA = (Real) 1e-8;

  // indicate that we are starting
  bool restart = true;

  FILE_LOG(LOG_OPT) << "BFGS entered" << endl;

  // setup the inverse Hessian approximation
  inv_hess.set_identity(params.n);

  // setup the number of BFGS iterations (for output)
  params.iterations = 0;

  // get the gradient
  (*params.grad)(x, params.data, g);

  // loop until gradient norm smaller than epsilon
  for (unsigned i=0; i< params.max_iterations && g.norm() > params.eps; i++)
  {
    // check whether we can terminate now
    if (params.tcheck && (*params.tcheck)(x, params.data))
    {
      FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::BGFS() exited" << endl;
      return x;
    }

    // compute search direction
    inv_hess.mult(g, p);
    p.negate();

    // verify that we have a descent direction
    if (p.dot(g) >= (Real) 0.0)
    {
      FILE_LOG(LOG_OPT) << " -- p is not a descent direction; restarting with identity Hessian" << endl;
      restart = true;
      inv_hess.set_identity();
      continue;
    }

    FILE_LOG(LOG_OPT) << "BFGS iteration: " << i << endl;
    FILE_LOG(LOG_OPT) << "  (inv) Hessian: " << endl << inv_hess;
    FILE_LOG(LOG_OPT) << "  x: " << x << endl;
    FILE_LOG(LOG_OPT) << "  search direction: " << p << endl;

    // determine alpha using backtracking line search
    Real alpha = search_line_backtrack(x, params.f, params.grad, p, (Real) 0.05, (Real) 0.5, params.data); 

    FILE_LOG(LOG_OPT) << "  alpha: " << alpha << endl;

    // see whether we are done
    if (alpha == (Real) 0.0)
    {
      // if we're on the first iteration, Hessian is probably not scaled
      // correctly
      if (restart)
      {
        FILE_LOG(LOG_OPT) << " -- alpha determined via line search is 0.0; rescaling inverse Hessian by " << endl;
        FILE_LOG(LOG_OPT) << "    " << BETA << " and trying again..." << endl;
        inv_hess *= BETA;
        continue;
      }

      FILE_LOG(LOG_OPT) << " -- alpha determined via line search is 0.0; terminating BFGS now" << endl;
      FILE_LOG(LOG_OPT) << "BFGS exited" << endl;
      return x;
    }

    // determine xnew
    delta_x.copy_from(p);
    delta_x *= alpha;
    x += delta_x;

    // compute new gradient and determine delta_g
    (*params.grad)(x, params.data, gnew);
    delta_g.copy_from(gnew);
    delta_g -= g;
    g.copy_from(gnew);

    // if deltag'*deltax <= 0, do a more sophisticated line search to preserve
    // curvature condition
    if (delta_x.dot(delta_g) <= (Real) 0.0)
    {
      FILE_LOG(LOG_OPT) << "  -- delta_x'*delta_g = " << delta_x.dot(delta_g) << " <= 0; doing more sophisticated line search" << endl;

      // determine alpha using line search
      alpha = search_line(x, params.f, params.grad, p, params.data, 
                          params.f_tol, params.grad_tol, params.x_rtol, 
                          params.min_alpha_step, params.max_alpha_step, 
                          params.max_fn_evals);

      FILE_LOG(LOG_OPT) << "    ++ determined alpha: " << alpha << endl;

      // redetermine xnew
      delta_x.copy_from(p);
      delta_x *= alpha;
      x += delta_x;

      // recompute new gradient and redetermine delta_g
      (*params.grad)(x, params.data, gnew);
      delta_g.copy_from(gnew);
      delta_g -= g;
      g.copy_from(gnew);
    }

    // see whether we are done
    if (delta_g.norm() == (Real) 0.0)
    {
      FILE_LOG(LOG_OPT) << " -- delta g is 0.0; terminating BFGS now" << endl;
      FILE_LOG(LOG_OPT) << "BFGS exited" << endl;
      return x;
    }

    // determine rho
    Real rho = (Real) 1.0 / delta_x.dot(delta_g);

    // determine inverse Hessian update
    // H_{k+1} = (I - rsy')H_k(I - rys') + rss'
    // H_{k+1} = (H_k - rsy' H_k)(I - rys') + rss'
    // H_{k+1} = H_k - H_k rys' - rsy' H_k + rsy' H_k rys' + rss'
    // \Delta H_{k+1} = rsy' H_k rys' + rss' - H_k rys' - rsy' H_k

    // compute s*s'
    VectorN::outer_prod(delta_x, delta_x, &ss);

    // compute r*y'*H_k*y*r
    Real ryHyr = delta_g.dot(inv_hess * delta_g) * rho * rho;

    // compute 2*H_k rys
    inv_hess.mult(delta_g, Hy);
    VectorN::outer_prod(Hy, delta_x, &Hrys);
    Hrys *= rho;

    // compute r*s*s'
    ssrho.copy_from(ss);
    ssrho *= rho;

    // compute rsy' H_k rys'
    ss1.copy_from(ss);
    ss1 *= ryHyr;

    // update inverse Hessian
    inv_hess += ss1;
    inv_hess += ssrho;
    inv_hess -= Hrys;
    Hrys.transpose();
    inv_hess -= Hrys;

    // indicate we are not restarting
    restart = false;

    // update the number of BFGS iterations
    params.iterations++;
  } 

  FILE_LOG(LOG_OPT) << "BFGS exited" << endl;

  return x; 
}

/// Does a backtracking line search
Real Optimization::search_line_backtrack(const VectorN& x, Real (*fn)(const VectorN&, void*), void (*grad)(const VectorN&, void*, VectorN&), const VectorN& dx, Real alpha, Real beta, void* data) 
{
  // precompute some things
  Real s = (Real) 1.0;
  Real f_knot = (*fn)(x, data);
  SAFESTATIC VectorN g, xp;
  (*grad)(x, data, g);
  const Real grad_dot_dx = g.dot(dx);
  const Real dx_norm = dx.norm();

  // setup x+dx*s
  xp = dx;
  xp *= s;
  xp += x;

  // now do the line search
  Real f_prime = (*fn)(xp, data);
  while (f_prime > f_knot + alpha * s * grad_dot_dx)
  {
    if (dx_norm*s < std::numeric_limits<Real>::epsilon())
      return 0.0;
    s *= beta;
    xp = dx;
    xp *= s;
    xp += x;
    f_prime = (*fn)(xp, data);
  }

  return s;
}

/// Brent's method for univariation minimization
/**
 * \param x_lower the left end of the interval to search
 * \param x_supper the right end of the interval to search
 * \param x the initial estimate (f(x) < f(x_lower) and f(x) < f(x_upper))
 *        and the input that yields the optimal value on return
 * \param fx the optimal value on return
 * \param f the function to optimize
 * \param params parameters to pass to f
 * \param eps the tolerance
 * \return <b>true</b> if successful, <b>false</b> if does not converge
 */
bool Optimization::brent(Real x_lower, Real x_upper, Real& x, Real& fx, Real (*f)(Real, void*), Real eps, void* params)
{
  const unsigned MAX_ITER = 100;
  const Real GOLDEN = 0.3819660;

  Real v = GOLDEN * x_upper;
  Real w = v;
  Real d = 0;
  Real e = 0;
  Real f_vw = (*f)(v, params);
  Real f_v = f_vw;
  Real f_w = f_vw;

  // set the minimum
  Real f_lower = (*f)(x_lower, params);
  Real f_upper = (*f)(x_upper, params);
  fx = (*f)(x, params);
  if (fx > f_lower+std::numeric_limits<Real>::epsilon() || 
      fx > f_upper+std::numeric_limits<Real>::epsilon())
    return false;

  for (unsigned iter=0; iter< MAX_ITER; iter++)
  {
    const Real x_left = x_lower;
    const Real x_right = x_upper;
    const Real z = x;
    const Real f_z = fx;
    std::swap(d, e);
    Real u, f_u;

    Real w_lower = (z - x_left);
    Real w_upper = (x_right - z);

    const Real tolerance =  1.4901161193847656e-08 * std::fabs (z);

    Real p = 0, q = 0, r = 0;

    const Real midpoint = 0.5 * (x_left + x_right);
    if (std::fabs(x-midpoint) <= 2.0*(eps*std::fabs(x)+1e-10) - 0.5*(x_right-x_left))
      return true;

    if (std::fabs(e) > tolerance)
    {
      /* fit parabola */
      r = (z - w) * (f_z - f_v);
      q = (z - v) * (f_z - f_w);
      p = (z - v) * q - (z - w) * r;
      q = 2 * (q - r);

      if (q > 0)
        p = -p;
      else
        q = -q;

      r = e;
      e  = d;
    }

    if (std::fabs (p) < std::fabs (0.5 * q * r) && p < q * w_lower && p < q * w_upper)
    {
      Real t2 = 2 * tolerance ;
      d = p / q;
      u = z + d;

      if ((u - x_left) < t2 || (x_right - u) < t2)
        d = (z < midpoint) ? tolerance : -tolerance ;
    }
    else
    {
      e = (z < midpoint) ? x_right - z : -(z - x_left) ;
      d = GOLDEN * e;
    }

    if (std::fabs (d) >= tolerance)
      u = z + d;
    else
      u = z + ((d > 0) ? tolerance : -tolerance) ;

    f_u = (*f)(u, params);

    if (f_u <= f_z)
    {
      if (u < z)
      {
        x_upper = z;
        f_upper = f_z;
      }
      else
      {
        x_lower = z;
        f_lower = f_z;
      }

      v = w;
      f_v = f_w;
      w = z;
      f_w = f_z;
      x = u;
      fx = f_u;
    }
    else
    {
      if (u < z)
      {
        x_lower = u;
        f_lower = f_u;
        continue;
      }
      else
      {
        x_upper = u;
        f_upper = f_u;
        continue;
      }

      if (f_u <= f_w || w == z)
      {
        v = w;
        f_v = f_w;
        w = u;
        f_w = f_u;
        continue;
      }
      else if (f_u <= f_v || v == z || v == w)
      {
        v = u;
        f_v = f_u;
        continue;
      }
    }
  }

  return false;
}

/// Does line search using algorithm from More and Thuente
/**
 * \param x the point from which to search
 * \param fn pointer to the objective function
 * \param grad pointer to the gradient function
 * \param sdir direction to search
 * \param data data (if any) passed to fn and grad
 * \param f_tol sufficient decrease condition
 * \param grad_tol directional derivative condition
 * \param x_rtol minimum relative width of the interval of uncertainty
 * \param step_min the minimum step size (0.0 by default)
 * \param step_max the maximum step size (1e6 by default)
 * \param max_fn_evals the maximum number of function evaluations
 * \return the best determined value of alpha, where f(x + sdir*alpha) is the 
           function minimum
 */
Real Optimization::search_line(const VectorN& x, Real (*fn)(const VectorN&, void*), void (*grad)(const VectorN&, void*, VectorN&), const VectorN& sdir, void* data, Real f_tol, Real grad_tol, Real x_rtol, Real step_min, Real step_max, unsigned max_fn_evals)
{
  FILE_LOG(LOG_OPT) << "Optimization::search_line() entered" << endl;

  // make a copy of x
  SAFESTATIC VectorN xprime;
  xprime = x;

  // setup some constants
  const Real p5 = (Real) 0.5;
  const Real p66 = (Real) 0.66;
  const Real xtrapf = (Real) 4.0;

  // setup initial guess for step
  Real step = (Real) 1.0;

  // setup initial value for infoc
  unsigned infoc = 1;

  // check the input parameters for errors
  assert(step > 0.0);
  assert(f_tol >= 0.0);
  assert(grad_tol >= 0.0);
  assert(x_rtol >= 0.0);
  assert(step_min >= 0.0);
  assert(step_max > step_min);

  // compute the initial gradient in the search direction and check that s is
  // a descent direction
  VectorN g;
  (*grad)(x, data, g);
  Real dginit = g.dot(sdir);
  if (dginit >= 0.0)
  {
    FILE_LOG(LOG_OPT) << " -- initial gradient in search direction is not a descent direction!" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::search_line() terminating..." << endl;
    return 0.0;
  }

  // initialize local variables
  bool brackt = false;
  bool stage1 = true;
  unsigned nfev = 0;
  Real finit = (*fn)(x, data);
  Real dgtest = f_tol * dginit;
  Real width = step_max - step_min;
  Real width1 = 2.0*width;
  VectorN wa = x;

  // the variables stx, fx, dgx contain the values of the step, function, and
  // directional derivative at the best step.  The variables sty, fy, dgy
  // contain the value of the step, function, and derivative at the other
  // endpoint of the interval of uncertainty.  The variables step, f, dg contain
  // the values of the step, function, and derivative at the current step.
  Real stx = (Real) 0.0;
  Real fx = finit;
  Real dgx = dginit;
  Real sty = (Real) 0.0;
  Real fy = finit;
  Real dgy = dginit;

  // start of iteration
  while (true)
  {
    // set the minimum and maximum steps to correspond to the present level of
    // uncertainty
    Real stmin, stmax;
    if (brackt)
    {
      stmin = std::min(stx, sty);
      stmax = std::max(stx, sty);
    }
    else
    {
      stmin = stx;
      stmax = step + xtrapf*(step-stx);
    }

    FILE_LOG(LOG_OPT) << "  min step: " << stmin << "  max step: " << stmax << endl;

    // force the step to be within the bounds step_max and step_min
    step = std::max(step, step_min);
    step = std::min(step, step_max);

    // if an unusual termination is to occur then let step be the lowest point
    // obtained thus far
    if ((brackt && (step <= stmin || step >= stmax)) || nfev > max_fn_evals ||
        infoc == 0 || (brackt && stmax-stmin <= x_rtol*stmax))
      step = stx;

    // evaluate the function and gradient at step and compute the directional
    // derivative
    xprime = wa + step*sdir;
    Real f = (*fn)(xprime, data);
    nfev++;
    (*grad)(xprime, data, g);
    Real dg = g.dot(sdir);
    Real ftest1 = finit + step*dgtest;

    // test for convergence
    if ((brackt && (step <= stmin || step >= stmax)) || infoc == 0)
    {
      FILE_LOG(LOG_OPT) << " -- rounding errors prevent further progress" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }
    else if (step == step_max && f <= ftest1 && dg <= dgtest)
    {
      FILE_LOG(LOG_OPT) << " -- step is at the upper bound (step_max)" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }
    else if (step == step_min && (f > ftest1 || dg >= dgtest))
    {
      FILE_LOG(LOG_OPT) << " -- step is at the lower bound (step_min)" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }
    else if (nfev >= max_fn_evals)
    {
      FILE_LOG(LOG_OPT) << " -- number of calls has reached maximum (max_fn_evals)" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }
    else if (brackt && step_max-step_min <= x_rtol*step_max)
    {
      FILE_LOG(LOG_OPT) << " -- relative width of the interval of uncertainty is at most x_rtol" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }
    else if (f <= ftest1 && std::fabs(dg) <= grad_tol*(-dginit))
    {
      FILE_LOG(LOG_OPT) << " -- sufficient decrease and directional derivative conditions hold" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::search_line() exiting" << endl;
      return stx;
    }

    // in the first stage we seek a step for which the modified function has
    // a nonpositive value and nonnegative derivative
    if (stage1 && f <= ftest1 && dg >= std::min(f_tol, grad_tol)*dginit)
      stage1 = false;

    // a modified function is used to predict the step only if we have not
    // obtained a step for which the modified function has a nonpositive
    // function value and nonnegative derivative, and if a lower function value
    // has been obtained but the decrease is not sufficient
    if (stage1 && f <= fx && f > ftest1)
    {
      // define the modified function and derivative values
      Real fm = f - step*dgtest;
      Real fxm = fx - stx*dgtest;
      Real fym = fy - sty*dgtest;
      Real dgm = dg - dgtest;
      Real dgxm = dgx - dgtest;
      Real dgym = dgy - dgtest;

      // call cstep to update the interval of uncertainty and to compute the
      // new step
      infoc = cstep(step_min, step_max, stx, fxm, dgxm, sty, fym, dgym, step, fm, dgm, brackt);

      // reset the function and gradient values for f
      fx = fxm + stx*dgtest;
      fy = fym + sty*dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    }
    else
    {
      // call cstep to update the interval of uncertainty and to compute the
      // new step
      infoc = cstep(step_min, step_max, stx, fx, dgx, sty, fy, dgy, step, f, dg, brackt);
    }

    // force a significant decrease in the size of the interval of uncertainty
    if (brackt)
    {
      if (std::fabs(sty-stx) >= p66*width1)
        step = stx + p5*(sty - stx);
      width1 = width;
      width = std::fabs(sty-stx);
    }
  }
} 

/// Helper function for line search
/**
 * \param stmin the lower bound for the step
 * \param stmax the upper bound for the step
 * \param stx the step at the best step obtained so far
 * \param fx the function value at the best step obtained so far
 * \param dx the derivative at the best step obtained so far
 * \param sty the step at the other endpoint of interval of uncertainty
 * \param fy the function value at the other endpoint of interval of uncertainty
 * \param dy the derivative at the other endpoint of interval of uncertainty
 * \param stp the current step; on return the new step
 * \param fp the function value at the current step
 * \param dy the derivative at the current step
 * \param brackt if true, the minimizer has been bracketed; if brackt is set
 *        to true on input, then on input stp must be between stx and sty
 */
int Optimization::cstep(Real stmin, Real stmax, Real& stx, Real& fx, Real& dx, Real& sty, Real& fy, Real& dy, Real& stp, Real& fp, Real& dp, bool& brackt)
{
  int info = 0;

  FILE_LOG(LOG_OPT) << " -- cstep() entered; minimium bracketed? " << brackt << endl;
  FILE_LOG(LOG_OPT) << "   stx:" << stx << "  fx:" << fx << "  dx:" << dx << endl;
  FILE_LOG(LOG_OPT) << "   sty:" << sty << "  fy:" << fy << "  dy:" << dy << endl;
  FILE_LOG(LOG_OPT) << "   stp:" << stp << "  fp:" << fp << "  dp:" << dp << endl;

  // check the input parameters for errors
  assert(!brackt || (stp > std::min(stx,sty) && stp < std::max(stx,sty)));
  assert(dx*(stp-stx) < (Real) 0.0);
  assert(stmax >= stmin);

  // determine whether derivatives have opposite sign
  Real sgnd = dp * (dx / std::fabs(dx));

  // first case: a higher function value.  the minimum is bracketed.  if the
  // cubic step is closer to stx than the quadratic step, the cubic step is
  // taken, else the average of the cubic and quadratic steps is taken.
  bool bound;
  Real theta;
  Real s;
  Real gamma;
  Real p,q,r;
  Real stpc, stpq, stpf;

  if (fp > fx)
  {
    FILE_LOG(LOG_OPT) << " -- cstep() case 1 followed" << endl;
    info = 1;
    bound = true;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    s = absmax(theta, dx, dp);
    gamma = s * std::sqrt(((theta/s)*(theta/s))-(dx/s)*(dp/s));
    if (stp < stx)
      gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r * (stp-stx);
    stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2) * (stp-stx);
    if (std::fabs(stpc-stx) < std::fabs(stpq-stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / 2;

    brackt = true;
  }

  // second case: a lower function value and derivatives of opposite sign.
  // the minimum is bracketed.  if the cubic step is closer to stx than the 
  // quadratic (secant) step, the cubic step is taken, else the quadratic step
  // is taken.
  else if (sgnd < (Real) 0.0)
  {
    FILE_LOG(LOG_OPT) << " -- cstep() case 2 followed" << endl;
    info = 2;
    bound = false;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    s = absmax(theta,dx,dp);
    gamma = s * std::sqrt(((theta/s) * (theta/s)) - (dx/s) * (dp/s));
    if (stp > stx)
      gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (std::fabs(stpc - stp) > std::fabs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  }

  // third case: a lower function value, derivatives of the same sign, and the
  // magnitude of the derivative decreases.  The cubic step is only used if
  // the cubic tends to infinity in the direction of the step or if the min
  // of the cubic is beyond stp.  Otherwise the cubic step is defined to be
  // either stmin or stmax.  The quadratic (secant) step is also computed and
  // if the minimum is bracketed then the step closest to stx is taken, else
  // the step farthest away is taken

  else if (std::fabs(dp) < std::fabs(dx))
  {
    FILE_LOG(LOG_OPT) << " -- cstep() case 3 followed" << endl;
    info = 3;
    bound = true;
    theta = 3 * (fx-fp) / (stp - stx) + dx + dp;
    s = absmax(theta, dx, dp);

    // thet case gamma = 0 only arises if the cubic does not tend to inf
    // in the direction of the step
    gamma = s * std::sqrt(std::max((Real) 0.0,(theta/s) * (theta/s) - (dx/s) * (dp/s)));
    if (stp > stx)
      gamma = -gamma;

    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p / q;
    if ((r < (Real) 0.0) && (gamma != (Real) 0.0))
      stpc = stp + r * (stx - stp);
    else if (stp > stx)
      stpc = stmax;
    else
      stpc = stmin;

    stpq = stp + (dp/ (dp - dx)) * (stx - stp);
    if (brackt)
    {
      if (std::fabs(stp - stpc) < std::fabs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    }
    else
    {
      if (std::fabs(stp - stpc) > std::fabs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    }
  }

  // 4th case: a lower function value, derivatives of the same sign, and the
  // magnitude of the derivative does not decrease.  If the minimum is not
  // bracketed, the step is either stmin or stmax, else the cubic step is taken
  else
  {
    FILE_LOG(LOG_OPT) << " -- cstep() case 4 followed" << endl;
    info = 4;
    bound = false;
    if (brackt)
    {
      theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      s = absmax(theta, dy, dp);
      gamma = s * std::sqrt(((theta/s)*(theta/s)) - (dy/s)*(dp/s));
      if (stp > sty)
        gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    }
    else if (stp > stx)
      stpf = stmax; 
    else
      stpf = stmin;
  }

  // update the interval of uncertainty.  this update does not depend on the
  // new step or the case analysis above
  if (fp > fx)
  {
    sty = stp;
    fy = fp;
    dy = dp;
  }
  else
  {
    if (sgnd < (Real) 0.0)
    {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }

  // compute the new step and safeguard it
  stpf = std::min(stmax, stpf);
  stpf = std::max(stmin, stpf);
  stp = stpf;
  if (brackt && bound)
  {
    if (sty > stx)
      stp = std::min(stx + (Real) 0.66 * (sty - stx), stp);
    else
      stp = std::max(stx + (Real) 0.66 * (sty - stx), stp);
  }

  FILE_LOG(LOG_OPT) << " -- cstep() about to exit" << endl;
  FILE_LOG(LOG_OPT) << "   stx:" << stx << "  fx:" << fx << "  dx:" << dx << endl;
  FILE_LOG(LOG_OPT) << "   sty:" << sty << "  fy:" << fy << "  dy:" << dy << endl;
  FILE_LOG(LOG_OPT) << "   stp:" << stp << "  fp:" << fp << "  dp:" << dp << endl;
  FILE_LOG(LOG_OPT) << " -- cstep() exited; minimium bracketed? " << brackt << endl;

  return info;
}

/// Helper function for line search
Real Optimization::absmax(Real a, Real b, Real c)
{
  a = std::fabs(a);
  b = std::fabs(b);
  c = std::fabs(c);
  if (a > b)
    return (a > c) ? a : c;
  else
    return (b > c) ? b : c;
}

/// Makes a point feasible for a convex optimization problem
/**
 * We do this by solving the convex optimization problem:
 * minimize 1's
 *   subject to     Ax  = b
 *              Gx + s >= m
 *                   s >= 0
 */
bool Optimization::make_feasible_convex2(OptParams& cparams, VectorN& x, void (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const VectorN&, const MatrixN&, const MatrixN&, const OptParams&, VectorN&))
{
  // necessary b/c of the max_element we do below
  if (x.size() == 0)
    return true;

  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex2() entered" << endl;

  // setup the feasibility data
  const unsigned n = x.size();
  const unsigned m = cparams.m;
  FeasibilityData fdata;
  fdata.oparams = &cparams;
  fdata.n = n;
  fdata.m = m;
  fdata.f0 = cparams.f0;
  fdata.fx = cparams.fx;
  fdata.gx = cparams.gx;
  fdata.grad0 = cparams.grad0;
  fdata.cJac_f = cparams.cJac_f;
  fdata.cJac_g = cparams.cJac_g;
  fdata.hess = cparams.hess;
  fdata.data = cparams.data;

  // determine new constraints
  const unsigned N_LB = cparams.lb.size();
  const unsigned N_UB = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();  

  // evaluate f at current x
  vector<Real> f, g;
  evaluate_inequality_constraints(x, cparams, f);
  evaluate_equality_constraints(x, cparams, g);

  // setup total inequality and equality constraints
  const unsigned mm = f.size(); 
  const unsigned rr = g.size();

  // get max element and minimum/maximum elements of g
  typedef vector<Real>::const_iterator VIter;
  VIter max_f = std::max_element(f.begin(), f.end());
  pair<VIter, VIter> mmax_g = boost::minmax_element(g.begin(), g.end());
    
  // verify that x is not already feasible
  if ((max_f == f.end() || *max_f < (Real) 0.0) &&
      (mmax_g.first == g.end() || 
       std::max(std::fabs(*mmax_g.first), *mmax_g.second) < cparams.eps_feas))
  {
    FILE_LOG(LOG_OPT) << " -- point is already feasible!" << endl;
    FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex2() exiting!" << endl;
    return true;
  }

  // setup y
  VectorN y(n+mm+rr*2);
  y.set_sub_vec(0, x);

  // setup inequality components of y
  for (unsigned i=0; i< mm; i++)
    y[n+i] = (f[i] >= (Real) 0.0) ? f[i]+NEAR_ZERO : NEAR_ZERO;

  // setup equality components of y
  for (unsigned i=0, j=n+mm; i< rr; i++, j+= 2)
  {
    if (g[i] > (Real) 0.0)
    {
      y[j] = g[i] + NEAR_ZERO;
      y[j+1] = (Real) 0.0;
    }
    else
    {
      y[j] = (Real) 0.0;
      y[j+1] = -g[i] + NEAR_ZERO;
    }
  }

  FILE_LOG(LOG_OPT) << " -- initial y: " << y << endl;

  // setup INF
  const Real INF = std::sqrt(std::numeric_limits<Real>::max());

  // setup convex optimization params
  OptParams cp;
  cp = cparams;
  cp.n = y.size();
  cp.m = mm + rr*2;
  cp.r = 0;
  cp.fx = &make_feasible_fx;
  cp.f0 = &make_feasible_f0;
  cp.grad0 = &make_feasible_grad0;
  cp.cJac_f = &make_feasible_cJac_f;
  cp.cJac_g = &make_feasible_cJac_g;
  cp.hess = &make_feasible_hess;
  cp.tcheck = &make_feasible_tcheck;
  cp.data = &fdata;
  cp.A.resize(0, y.size());
  cp.b.resize(0);
  cp.max_iterations = std::max(1000 + (unsigned) std::log(cp.n), cparams.max_iterations);
  cp.M.resize(0, y.size());
  cp.q.resize(0);
  cp.ub.resize(0);

  // setup lower bounds on variables
  cp.lb.resize(cp.n);
  for (unsigned i=0; i< cparams.n; i++)
    cp.lb[i] = -INF;
  for (unsigned i=cparams.n; i< cp.n; i++)
    cp.lb[i] = (Real) 0.0;

  // do convex optimization with primal/dual method
  optimize_convex_pd(cp, y, solve_KKT);  

  // get x out
  y.get_sub_vec(0, n, x);

  // evaluate f and g at current x
  evaluate_inequality_constraints(x, cparams, f);
  evaluate_equality_constraints(x, cparams, g);

  // get max element and minimum/maximum elements of g
  max_f = std::max_element(f.begin(), f.end());
  mmax_g = boost::minmax_element(g.begin(), g.end());
   
  // redetermine v
  Real v = (Real) 0.0;
  if (max_f != f.end())
    v = std::max(v, -*max_f);

  // now use equality evaluations
  if (mmax_g.first != g.end())
    v = std::max(v, std::max(std::fabs(*mmax_g.first), *mmax_g.second));
 
  FILE_LOG(LOG_OPT) << "   x: " << x << endl;
  FILE_LOG(LOG_OPT) << "   v: " << v << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex2() successful? " << (v <= cparams.eps_feas) << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex2() exiting" << endl;

  // check feasibility 
  return (v <= cparams.eps_feas);
}

/// Makes a point feasible for a convex optimization problem
/**
 * We do this by solving the convex optimization problem:
 * minimize s
 *   subject to Ax = b
 *              Gx + s >= m
 */
bool Optimization::make_feasible_convex(OptParams& cparams, VectorN& x, void (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const VectorN&, const MatrixN&, const MatrixN&, const OptParams&, VectorN&))
{
  // necessary b/c of the max_element we do below
  if (x.size() == 0 || cparams.m + cparams.q.size() == 0)
    return true;

  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex() entered" << endl;

  // verify that there are no nonlinear inequality constraints
  if (cparams.r > 0)
    throw std::runtime_error("make_feasible_convex() does not support nonlinear equalities - use make_feasible_convex2() instead");

  // setup the feasibility data
  const unsigned n = x.size();
  const unsigned m = cparams.m;
  FeasibilityData fdata;
  fdata.oparams = &cparams;
  fdata.n = n;
  fdata.m = m;
  fdata.f0 = cparams.f0;
  fdata.fx = cparams.fx;
  fdata.grad0 = cparams.grad0;
  fdata.cJac_f = cparams.cJac_f;
  fdata.hess = cparams.hess;
  fdata.data = cparams.data;

  // evaluate f at current x
  vector<Real> f;
  evaluate_inequality_constraints(x, cparams, f);

  // *must* be at least one inequality constraint (from test above) 
  assert(!f.empty());

  // setup total inequality constraints
  const unsigned mm = f.size(); 

  // verify that x is not already feasible and setup v
  Real v = *std::max_element(f.begin(), f.end());
  if (v < (Real) 0.0)
  {
    FILE_LOG(LOG_OPT) << " -- point is already feasible!" << endl;
    FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex() exiting!" << endl;
    return true;
  }

  // setup y
  VectorN y(n+1);
  y.set_sub_vec(0, x);
  y[n] = v + (Real) 1.0; 
  FILE_LOG(LOG_OPT) << "initial v: " << v << endl;

  // make new matrix A to account for increased size of y
  MatrixN Anew(cparams.A.rows(), y.size());
  Anew.set_sub_mat(0,0,cparams.A);
  for (unsigned i=0; i< cparams.A.rows(); i++)
    Anew(i,n) = (Real) 0.0;

  // setup convex optimization params
  OptParams cp;
  cp = cparams;
  cp.n++;
  cp.m = mm;
  cp.fx = &make_feasible_fx;
  cp.f0 = &make_feasible_f0;
  cp.grad0 = &make_feasible_grad0;
  cp.cJac_f = &make_feasible_cJac_f;
  cp.hess = &make_feasible_hess;
  cp.tcheck = &make_feasible_tcheck;
  cp.data = &fdata;
  cp.A = Anew;
  cp.max_iterations = std::max(1000 + (unsigned) std::log(cp.n), cparams.max_iterations);
  cp.M.resize(0, y.size());
  cp.q.resize(0);
  cp.lb.resize(0);
  cp.ub.resize(0);

  // do convex optimization with primal/dual method
  optimize_convex_pd(cp, y, solve_KKT);  

  // get x out
  y.get_sub_vec(0, n, x);

  // redetermine v
  v = *std::max_element(f.begin(), f.end());
  evaluate_inequality_constraints(x, cparams, f);

  FILE_LOG(LOG_OPT) << "   x: " << x << endl;
  FILE_LOG(LOG_OPT) << "   v: " << v << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex() successful? " << (v <= (Real) 0.0) << endl;
  FILE_LOG(LOG_OPT) << " -- Optimization::make_feasible_convex() exiting" << endl;

  // check feasibility 
  return (v <= (Real) 0.0);
}

/// Termination check function for determining whether the convex optimization procedure can terminate early
bool Optimization::make_feasible_tcheck(const VectorN& y, void* data)
{
  SAFESTATIC VectorN x, work;

  // get the feasibility data
  FeasibilityData& fdata = *((FeasibilityData*) data);

  // determine how many constraints
  unsigned m = fdata.m;
  unsigned n = fdata.n;

  // get the feasibility tolerance (for equality constraints)
  const OptParams& oparams = *fdata.oparams;  
  Real eps_feas = oparams.eps_feas;

  // check whether lower bounds are met
  for (unsigned i=0; i< oparams.lb.size(); i++)
    if (y[i] < oparams.lb[i])
      return false;

  // check whether upper bounds are met
  for (unsigned i=0; i< oparams.ub.size(); i++)
    if (y[i] > oparams.ub[i])
      return false;

  // create new vector with 's' component(s) removed
  y.get_sub_vec(0, n, x);

  // check linear inequality constraints
  oparams.M.mult(x, work) -= oparams.q;
  for (unsigned i=0; i< work.size(); i++)
    if (work[i] < (double) 0.0)
      return false;

  // check linear equality constraints
  oparams.A.mult(x, work) -= oparams.b;
  for (unsigned i=0; i< work.size(); i++)
    if (std::fabs(work[i]) > eps_feas)
      return true;

  // determine whether any element of f is >= 0
  if (fdata.fx)
  {
    fdata.fx(x, work, fdata.data);
    Real* maxele = std::max_element(work.begin(), work.end());
    if (maxele != work.end() && *maxele >= (Real) 0.0)
      return false;
  }

  // determine whether any element of g is > tol 
  if (fdata.gx)
  {
    fdata.gx(x, work, fdata.data);
    pair<Real*, Real*> mmax = boost::minmax_element(work.begin(), work.end());
    if (mmax.first != work.end() && 
        std::max(std::fabs(*mmax.first), *mmax.second) > oparams.eps_feas)
      return false;
  }

  // still here?  all constraint functions are met 
  return true;
}

/// Evaluation function for making a point feasible
Real Optimization::make_feasible_f0(const VectorN& y, void* data)
{
  // get the feasibility data
  FeasibilityData& fdata = *((FeasibilityData*) data);

  // determine how many constraints
  unsigned n = fdata.n;
 
  // ** using method I
  if (y.size() == n+1)
    return y[n];
  else
    // ** using method II
    return std::accumulate(y.begin()+n, y.end(), (Real) 0.0);
} 

/// Evaluation function for making a point feasible
void Optimization::make_feasible_fx(const VectorN& y, VectorN& f, void* data)
{
  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);
  const OptParams& oparams = *fdata.oparams;

  // determine how many constraints
  const unsigned m = fdata.m;
  const unsigned n = fdata.n;
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;

  // create new vector y with v removed
  SAFESTATIC VectorN x;
  SAFESTATIC vector<Real> fx;
  y.get_sub_vec(0, n, x);

  // evaluate original constraints
  evaluate_inequality_constraints(x, oparams, fx);

  // determine which method we are using
  if (y.size() == n+1)
  {
    // ** using method I
    // get v
    Real v = y[n];

    // subtract v from fx
    std::transform(fx.begin(), fx.end(), fx.begin(), std::bind2nd(std::minus<Real>(), v));
  }
  else
  {
    // ** using method II

    // subtract v from fx
    std::transform(fx.begin(), fx.end(), y.begin()+n, fx.begin(), std::minus<Real>());
  } 
    
  // copy fx to f
  f.resize(fx.size());
  std::copy(fx.begin(), fx.end(), f.begin());
} 

/// Equality constraint evaluation function for making a point feasible
void Optimization::make_feasible_gx(const VectorN& y, VectorN& g, void* data)
{
  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);
  const OptParams& oparams = *fdata.oparams;

  // determine how many constraints
  const unsigned m = fdata.m;
  const unsigned n = fdata.n;
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;

  // create new vector y with v removed
  SAFESTATIC VectorN x;
  SAFESTATIC vector<Real> gx;
  y.get_sub_vec(0, n, x);
  const unsigned rr = gx.size();

  // evaluate original constraints
  evaluate_equality_constraints(x, oparams, gx);

  // must be using method II
  assert(y.size() == n + mm + rr*2);

  // evaluate the equality constraint
  for (unsigned i=0, j=n+mm; i< gx.size(); i++, j+=2)
    g[i] = gx[i] - y[j] + y[j+1];
}

/// Gradient function for making a point feasible
/**
 * \note this function needs to be tested
 */
void Optimization::make_feasible_grad0(const VectorN& y, VectorN& g, void* data)
{
  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);
  const OptParams& oparams = *fdata.oparams;

  // determine how many constraints
  const unsigned m = fdata.m;
  const unsigned n = fdata.n;
  const unsigned nn = y.size();
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;

  // compute gradient - works the same for both methods
  g.set_zero(nn);
  for (unsigned i=n; i< nn; i++)
    g[i] = (Real) 1.0;
}

/// Constraint gradient function for making a point feasible
void Optimization::make_feasible_cJac_f(const VectorN& y, MatrixN& J, void* data)
{
  SAFESTATIC VectorN x;
  SAFESTATIC MatrixN Jorig;

  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);
  const OptParams& oparams = *fdata.oparams;

  // determine how many constraints
  unsigned m = fdata.m;
  unsigned n = fdata.n;
  unsigned nn = y.size();
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;

  // resize the augmented Jacobian
  J.resize(mm,nn);
  Real* Jdata = J.data();

  // get the original vector
  y.get_sub_vec(0, n, x);

  // get the origin Jacobian and copy it
  Jorig.resize(m,n);
  if (m > 0)
  {
    assert(fdata.cJac_f);
    fdata.cJac_f(x, Jorig, data);
  }
  J.set_sub_mat(0,0,Jorig);

  // add lower bounds -- not incorporating v yet
  for (unsigned i=0, j=m; i< N_LB; i++, j++)
  {
    // first, zero the row
    BlockIterator bstart = J.block_start(j, j+1, 0, nn);
    std::fill_n(bstart, nn, (Real) 0.0);
    
    // now, set the gradient portion
    J(j, i) = (Real) -1.0;
  }

  // add upper bounds -- not incorporating v yet
  for (unsigned i=0, j=m+N_LB; i< N_UB; i++, j++)
  {
    // first, zero the row
    BlockIterator bstart = J.block_start(j, j+1, 0, nn);
    std::fill_n(bstart, nn, (Real) 0.0);

    // now, set the gradient portion
    J(j, i) = (Real) 1.0;
  }

  // add linear inequalities -- not incorporating v yet
  J.set_sub_mat(m+N_LB+N_UB, 0, oparams.M);
  for (unsigned i=0, j=m+N_LB+N_UB; i< N_LININEQ; i++, j++)
  {
    // negate the first part of the row 
    BlockIterator bstart = J.block_start(j, j+1, 0, n);
    BlockIterator bend = J.block_end(j, j+1, 0, n);
    std::transform(bstart, bend, bstart, std::negate<Real>());

    // zero the second part of the row 
    bstart = J.block_start(j, j+1, n, nn);
    std::fill_n(bstart, nn-n, (Real) 0.0);
  }

  // now add in the constraints for v 
  if (nn - n == 1)
  {
    // method 1 - need to subtract 1.0 from each constraint row to account for v
    // f(x) <= v
    BlockIterator bstart = J.block_start(0, mm, n, n+1);
    std::fill_n(bstart, mm, (Real) -1.0);
  }
  else
  {
    // setup the proper gradient terms for "v" inequalities
    for (unsigned i=0, j=n; i< mm; i++, j++)
      J(i, j) = (Real) -1.0;
  }
}

/// Constraint gradient function for making a point feasible
void Optimization::make_feasible_cJac_g(const VectorN& y, MatrixN& J, void* data)
{
  SAFESTATIC VectorN x;
  SAFESTATIC MatrixN Jred;

  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);
  const OptParams& oparams = *fdata.oparams;

  // determine how many constraints
  unsigned m = fdata.m;
  unsigned n = fdata.n;
  unsigned nn = y.size();
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();
  const unsigned mm = m + N_LB + N_UB + N_LININEQ;
  const unsigned rr = oparams.r + oparams.b.size();

  // get the original vector
  y.get_sub_vec(0, n, x);

  // get the equality gradients
  Jred.resize(rr, n);
  calc_Dg(x, oparams, Jred);

  // setup the new equality constrained Jacobian
  J.resize(rr, nn);

  // first, zero out columns
  BlockIterator bstart = J.block_start(0, rr, n+mm, n+mm+rr*2);
  std::fill_n(bstart, rr*rr*2, (Real) 0.0);

  // now, set "v" variables for equality constraints
  for (unsigned i=0, j=n+mm; i< rr; i++)
  {
    J(i, j++) = (Real) -1.0;
    J(i, j++) = (Real) 1.0;
  }
}

/// Hessian function for making a point feasible
void Optimization::make_feasible_hess(const VectorN& y, Real objscal, const VectorN& lambda, const VectorN& nu, MatrixN& H, void* data)
{
  // get the feasibility data
  const FeasibilityData& fdata = *((FeasibilityData*) data);

  // verify that nu is empty
  assert(nu.size() == 0);

  // determine how many constraints
  const unsigned m = fdata.m;
  const unsigned n = fdata.n;
  const unsigned nn = y.size();

  // setup original Hessian
  SAFESTATIC MatrixN Horig;
  SAFESTATIC VectorN x;
  y.get_sub_vec(0, n, x);

  // call the original Hessian functions
  Horig.resize(n,n);
  fdata.hess(x, (Real) 0.0, lambda, nu, Horig, fdata.data);

  // setup the new Hessian
  H.resize(nn,nn);
  H.set_sub_mat(0,0,Horig);

  // set the columns [n+1..nn] to zero
  BlockIterator colstart = H.block_start(0,nn,n,nn);
  std::fill_n(colstart, (nn-n)*nn, (Real) 0.0);

  // set the rows [n+1..nn] to zero
  BlockIterator rowstart = H.block_start(n,nn,0,nn);
  std::fill_n(rowstart, (nn-n)*nn, (Real) 0.0);
}

/// Performs convex optimization using the logarithmic barrier method
/**
 * \param cparams parameters for the optimization
 * \param x the feasible value from which to start; contains the optimal value
 *        on return
 * \return <b>true</b> if optimization successful; <b>false</b> otherwise
 */
bool Optimization::optimize_convex(OptParams& cparams, VectorN& x)
{
  const Real T_INIT = 1e-3;
  SAFESTATIC VectorN dx, g, gx, inv_f, workv, lambda;
  SAFESTATIC MatrixN H, outer_prod;
  SAFESTATIC MatrixN cJac;
  SAFESTATIC vector<Real> fx;

  FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() entered" << endl;

  // does not currently use equality constraints
  if (cparams.A.rows() > 0 || cparams.r > 0)
    throw std::runtime_error("optimize_convex() does not currently respect equality constraints!");

  // eliminate redundant constraints
  eliminate_redundant_constraints(cparams.A, cparams.b);

  // setup number of iterations for output
  cparams.iterations = 0;

  // init t
  Real t = T_INIT;

  // determine n and m
  const unsigned n = x.size();
  const unsigned m = cparams.m;
  const unsigned r = cparams.r;
  const unsigned N_LB = cparams.lb.size();
  const unsigned N_UB = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();

  // setup gradient 
  dx.resize(n);
  g.resize(n);
  H.resize(n,n);

  // verify that constraint functions are met
  evaluate_inequality_constraints(x, cparams, fx);
  std::vector<Real>::const_iterator minelm = std::min_element(fx.begin(), fx.end());
  if (minelm != fx.end() && *minelm >= (Real) 0.0)
  {
    FILE_LOG(LOG_OPT) << "initial point is not feasible!" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;
    return false;
  }

  // check whether we can terminate early
  if (cparams.tcheck && cparams.tcheck(x, cparams.data))
  {
    FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;
    return true;
  }

  // do outer iteration
  while (true)
  {
    FILE_LOG(LOG_OPT) << "outer iteration begun, t = " << t << endl;

    // **************************************************
    // centering step and update
    // **************************************************
    while (true)
    {
      FILE_LOG(LOG_OPT) << "iteration: " << cparams.iterations << endl;
      if (cparams.iterations++ > cparams.max_iterations)
      {
        FILE_LOG(LOG_OPT) << "maximum iterations reached..." << endl;
        FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;
        return false;
      }
      FILE_LOG(LOG_OPT) << "inner iteration begun" << endl;
      FILE_LOG(LOG_OPT) << " x: " << x << endl;

      // compute current function gradient of objective
      cparams.grad0(x, g, cparams.data);

      // compute the Jacobian of the nonlinear constraints
      if (m > 0)
        cparams.cJac_f(x, cJac, cparams.data);

      // compute inverses of f
      evaluate_inequality_constraints(x, cparams, fx);
      std::transform(fx.begin(), fx.end(), fx.begin(), std::bind1st(std::divides<Real>(), (Real) -1.0));

      // compute the Hessians
      inv_f.get_sub_vec(0, m, lambda);
      cparams.hess(x, t, lambda, EMPTY_VEC, H, cparams.data); 

      // update the Hessian using gradients
      for (unsigned i=0; i< m; i++)
      {
        const Real sq_inv_f = fx[i]*fx[i];
        cJac.get_row(i, gx);
        workv.copy_from(gx) *= sq_inv_f;
        H += *VectorN::outer_prod(gx, workv, &outer_prod);
      }

      // compute the gradient contributions of the nonlinear inequalities
      g *= t;
      for (unsigned i=0; i< m; i++)
      {
        cJac.get_row(i, gx);
        gx *= inv_f[i];
        g += gx;
      }

      // compute the gradient contributions of the lower bounds
      for (unsigned i=0, j=m+r; i< N_LB; i++)
        g[i] -= inv_f[j++];

      // NOTE: we are accounting for future inclusion of r here
      // compute the gradient contributions of the upper bounds
      for (unsigned i=0, j=m+r+N_LB; i< N_UB; i++)
        g[i] += inv_f[j++];

      // compute the gradient contributions of the linear inequalities
      for (unsigned i=0, j=m+r+N_LB+N_UB; i< N_LININEQ; i++)
      {
        cparams.M.get_row(i, workv) *= -inv_f[j++];
        g += workv;
      }

      FILE_LOG(LOG_OPT) << " -- Hessian: " << endl << H;
      FILE_LOG(LOG_OPT) << " -- gradient: " << g << endl;

      // compute the step
      dx.copy_from(g);
      dx.negate();

      // condition the Hessian matrix
      condition_and_factor_PD(H);
      LinAlg::solve_chol_fast(H, dx);

      // compute lambda^2 -- must be >= 0, b/c hessian PD (or PSD, at worst)
      Real lambdasq = -VectorN::dot(g, dx);

      FILE_LOG(LOG_OPT) << " -- dx: " << dx << endl;
      FILE_LOG(LOG_OPT) << " -- lambdasq: " << lambdasq << endl;

      // check the stopping criterion
      if (lambdasq/2.0 <= cparams.eps)
      {
        FILE_LOG(LOG_OPT) << "lambda^2 too small (" << lambdasq << "); terminating inner iteration" << endl;
        break;
      }

      // do backtracking line search
      Real s = 1.0;
      Real f_knot = cparams.f0(x, cparams.data)*t;
      evaluate_inequality_constraints(x, cparams, fx);
      std::transform(fx.begin(), fx.end(), fx.begin(), std::negate<Real>());
      std::transform(fx.begin(), fx.end(), fx.begin(), nlog<Real>());
      f_knot -= std::accumulate(fx.begin(), fx.end(), (Real) 0.0);      
      const Real grad_dot_dx = -lambdasq;
      (workv.copy_from(dx) *= s) += x;
      Real f_prime = cparams.f0(workv, cparams.data)*t;
      evaluate_inequality_constraints(workv, cparams, fx);
      std::transform(fx.begin(), fx.end(), fx.begin(), std::negate<Real>());
      std::transform(fx.begin(), fx.end(), fx.begin(), nlog<Real>());
      f_prime -= std::accumulate(fx.begin(), fx.end(), (Real) 0.0);
      while (f_prime > f_knot + cparams.alpha * s * grad_dot_dx)
      {
        s *= cparams.beta;
        (workv.copy_from(dx) *= s) += x;
        f_prime = cparams.f0(workv, cparams.data)*t;
        evaluate_inequality_constraints(workv, cparams, fx);
        std::transform(fx.begin(), fx.end(), fx.begin(), std::negate<Real>());
        std::transform(fx.begin(), fx.end(), fx.begin(), nlog<Real>());
        f_prime -= std::accumulate(fx.begin(), fx.end(), (Real) 0.0);
      }

      FILE_LOG(LOG_OPT) << " -- true objective function at x: " << cparams.f0(x, cparams.data) << endl;
      FILE_LOG(LOG_OPT) << " -- barrier objective function at x: " << f_knot << endl;
      FILE_LOG(LOG_OPT) << " -- determined s: " << s << endl;
      FILE_LOG(LOG_OPT) << " -- objective function at x+dx*s: " << f_prime << endl; 
      FILE_LOG(LOG_OPT) << " -- f' - f0: " << (f_prime - f_knot) << endl;

      if (std::fabs(f_prime - f_knot) < cparams.eps)
      {
        FILE_LOG(LOG_OPT) << "f' approx equals f0; breaking out of inner iterations" << endl;
        break;
      }
/*
      // if there is not much change to the objective function, and this is
      // not the last outer iteration, break out of inner iteration
      if (std::fabs(f_knot - f_prime) < eps && 2.0*n/t < eps)
      {
        FILE_LOG(LOG_OPT) << "not much difference between f' and f0; breaking out of inner iterations" << endl;
        break;
      }

      // if dx * s is close to zero, we can't continue any further...
      if (dx.norm() * s < NEAR_ZERO)
      {
        FILE_LOG(LOG_OPT) << "||dx||*s approx 0; breaking out of inner iterations" << endl;
        break;
      }
*/

      // update x
      dx *= s;
      x += dx;
      if (cparams.tcheck && (*cparams.tcheck)(x, cparams.data))
      {
        FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
        FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;
        return true;
      }
    } // end Newton iteration

    // check stopping condition
    if (n/t < cparams.eps)
      break;

    // increase t
    t *= cparams.mu;
  }

  FILE_LOG(LOG_OPT) << "  optimum: " << x << endl;
  FILE_LOG(LOG_OPT) << "Optimization::optimize_convex() exited" << endl;

  return true;
}

/// Conditions (makes positive-definite) but does not factorize the Hessian
void Optimization::condition_hessian(MatrixN& H)
{
  if (H.rows() != H.columns())
    throw NonsquareMatrixException();

  // if Hessian empty, quit
  if (H.rows() == 0)
    return;

  FILE_LOG(LOG_OPT) << "Optimization::condition_hessian() entered" << endl;

  // make a copy of H
  SAFESTATIC MatrixN Hcopy;
  Hcopy.copy_from(H);

  // Hessian _must_ be positive-definite for the optimization to work
  const Real BETA = 1e-3;
  Real tau = Hcopy(0,0);
  for (unsigned i=1; i< Hcopy.rows(); i++)
    tau = std::min(tau, Hcopy(i,i));
  tau = (tau > (Real) 0.0) ? BETA : -tau + BETA;
  if (!LinAlg::factor_chol(Hcopy))
    while (true)
    {
      // update H += I*tau
      for (unsigned i=0; i< H.rows(); i++)
        Hcopy(i,i) += tau;

      // try Cholesky factorization
      if (LinAlg::factor_chol(Hcopy))
        break;

      // remove I*tau from H
      for (unsigned i=0; i< H.rows(); i++)
        Hcopy(i,i) -= tau;

      // update tau
      tau = std::max(tau*(Real) 2.0, BETA);
    }

  // update H with tau
  for (unsigned i=0; i< Hcopy.rows(); i++)
    H(i,i) += tau;

  FILE_LOG(LOG_OPT) << "  -- regularization factor: " << tau << endl;
  FILE_LOG(LOG_OPT) << "Optimization::condition_hessian() exited" << endl;
}

/// Conditions (makes positive-definite) and factorizes (Cholesky) the Hessian
void Optimization::condition_and_factor_PD(MatrixN& H)
{
  if (H.rows() != H.columns())
    throw NonsquareMatrixException();

  // if Hessian empty, quit
  if (H.rows() == 0)
    return;

  FILE_LOG(LOG_OPT) << "Optimization::condition_and_factor_PD() entered" << endl;

  // make a copy of H
  SAFESTATIC MatrixN Hcopy;
  Hcopy.copy_from(H);

  // try factorization
  if (LinAlg::factor_chol(H))
    return;

  // Hessian _must_ be positive-definite for the optimization to work
  H.copy_from(Hcopy);
  const Real BETA = 1e-3;
  Real tau = H(0,0);
  for (unsigned i=1; i< H.rows(); i++)
    tau = std::min(tau, H(i,i));
  tau = (tau > (Real) 0.0) ? BETA : -tau + BETA;
  while (true)
  {
    // update H += I*tau
    for (unsigned i=0; i< H.rows(); i++)
      H(i,i) += tau;

    // try Cholesky factorization
    if (LinAlg::factor_chol(H))
      break;

    // update tau
    tau = std::max(tau*(Real) 2.0, BETA);

    // recopy H
    H.copy_from(Hcopy);
  }

  FILE_LOG(LOG_OPT) << "  -- regularization factor: " << tau << endl;
  FILE_LOG(LOG_OPT) << "Optimization::condition_and_factor_PD() exited" << endl;
}

/// Performs convex optimization using the primal-dual interior point method [Boyd, 2004]
/**
 *  Minimizes   g(x)
 *  subject to  f(x) <= 0
 *              A*x = b
 * where f is vector valued, convex, and twice differentiable and A has full
 * row rank.
 * \param x the initial, feasible point and the optimal value on return
 * \param solve_KKT optional KKT solver that takes; 
 *        solve_KKT(x, lambda, nu, r, Df, Dg, cparams, dy) solves equation (11.54) for (11.53)
 *        in [Boyd, 2004]
 * \param solver function pointer to linear equation solver
 * \return <b>true</b> if optimization successful; <b>false</b> otherwise
 */
bool Optimization::optimize_convex_pd(OptParams& cparams, VectorN& x, void (*solve_KKT)(const VectorN&, const VectorN&, const VectorN&, const VectorN&, const MatrixN&, const MatrixN&, const OptParams&, VectorN&))
{
  const unsigned n = x.size();
  const unsigned m = cparams.m; 
  const unsigned cr = cparams.r;
  Real rdual_nrm, rpri_nrm, f0_best;
  VectorN x_best, q_m_Mx, fx, gx, r, rplus, dy;
  MatrixN Df, Dg;
  VectorN& nu = cparams.nu;
  VectorN& lambda = cparams.lambda;

  // setup infeasibility on inequality constraints to zero initially 
  Real infeas_tol = (Real) 0.0;

  // resize A if necessary
  if (cparams.A.columns() != x.size())
  {
    if (cparams.A.rows() != 0)
      throw MissizeException();
    cparams.A.resize(0, x.size());
    cparams.b.resize(0);
  }

  // resize M if necessary
  if (cparams.M.columns() != x.size())
  {
    if (cparams.M.rows() != 0)
      throw MissizeException();
    cparams.M.resize(0, x.size());
    cparams.q.resize(0);
  }

  // see whether we can quit immediately
  if (cparams.tcheck && (*cparams.tcheck)(x, cparams.data))
  {
    FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_pd() exited" << endl;
    return true;
  }

  // determine number of box and inequality constraints
  const unsigned N_LBOX = cparams.lb.size(); 
  const unsigned N_UBOX = cparams.ub.size(); 
  const unsigned N_LININEQ = cparams.q.size();
  const unsigned mm = m + N_LBOX + N_UBOX + N_LININEQ;

  // get tolerances
  const Real EPS = cparams.eps;
  const Real EPS_FEAS = cparams.eps_feas;

  // indicate we'll warm start until we know otherwise
  bool warm_start = true;

  // setup nu
  unsigned nu_len = cparams.A.rows() + cparams.r;
  if (nu_len != nu.size())
  {
    nu.set_zero(nu_len);
    warm_start = false;
  }

  // setup mu*m
  const Real mum = cparams.mu*mm;

  // eliminate redundant equality constraints
  eliminate_redundant_constraints(cparams.A, cparams.b);

  // determine f(x)
  f0_best = cparams.f0(x, cparams.data);
  x_best.copy_from(x);

  // verify linear inequality constraints are kosher
  if (N_LININEQ != cparams.M.rows())
    throw std::runtime_error("Number of rows for linear inequality constraint matrix is mismatched with vector");

  // verify that f is non-positive for constraint functions and setup lambda
  // and fc
  // NOTE: lambda must be greater than zero
  if (lambda.size() != mm)
  {
    lambda.set_zero(mm);
    warm_start = false;
  }

  // compute fc and gc
  vector<Real> fc(mm), gc(nu_len);

  // evaluate constraints
  evaluate_inequality_constraints(x, cparams, fc);
  evaluate_equality_constraints(x, cparams, gc);

  // now determine lambda
  for (unsigned i=0; i< lambda.size(); i++)
  {
    if (fc[i] >= (Real) 0.0)
    {
      infeas_tol = std::max(fc[i], infeas_tol);
      fc[i] = -NEAR_ZERO;
    }

    if (!warm_start)
      lambda[i] = (Real) -1.0/fc[i];
    else if (lambda[i] <= (Real) 0.0)
      lambda[i] = (Real) 1.0;
  }

  // update infeas_tol
  if (!gc.empty())
  {
    typedef vector<Real>::const_iterator VIter;
    pair<VIter, VIter> mmax = boost::minmax_element(gc.begin(), gc.end());
    infeas_tol = std::max(infeas_tol, std::max(std::fabs(*mmax.first), *mmax.second));
  }

  // copy fc to fx
  fx.resize(fc.size());
  std::copy(fc.begin(), fc.end(), fx.begin());

  // copy gc to gx
  gx.resize(gc.size());
  std::copy(gc.begin(), gc.end(), gx.begin());

  // increase infeasibility tolerance, if necessary
  if (infeas_tol > (Real) 0.0)
    infeas_tol += NEAR_ZERO;

  // setup eta
  Real eta = -(fx.dot(lambda));

  // init t
  Real t = mum/eta;

  // init y and r (and r+)
  VectorN y;
  y.resize(n+mm+nu_len);
  r.resize(0);
  rplus.resize(n+mm+nu_len);

  // init vectors and matrices so we don't have to keep freeing and 
  // reallocating memory 
  VectorN x_plus, dx, lambda_plus, dlambda, nu_plus, dnu;
  dy.resize(n+mm+nu_len);
  x_plus.resize(n);
  dx.resize(n);
  lambda_plus.resize(mm);
  dlambda.resize(mm);
  nu_plus.resize(nu_len);
  dnu.resize(nu_len);

  FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_pd() entered" << endl;
  FILE_LOG(LOG_OPT) << "initial x: " << x << endl;
  FILE_LOG(LOG_OPT) << "infeasability tolerance: " << infeas_tol << endl;
  FILE_LOG(LOG_OPT) << "f1..m(x): " << fx << endl;
  FILE_LOG(LOG_OPT) << "initial t: " << t << endl;

  for (cparams.iterations = 0; cparams.iterations < cparams.max_iterations; cparams.iterations++)  
  {
    // reset t
    t = mum/eta;

    FILE_LOG(LOG_OPT) << "-------------------------------------------------" << endl;
    FILE_LOG(LOG_OPT) << "iteration: " << cparams.iterations << endl;
    FILE_LOG(LOG_OPT) << "t: " << t << endl;
    FILE_LOG(LOG_OPT) << " x: " << x << endl;
    FILE_LOG(LOG_OPT) << " lambda: " << lambda << endl;
    FILE_LOG(LOG_OPT) << " nu: " << nu << endl;
    FILE_LOG(LOG_OPT) << "  infeasability tolerance: " << infeas_tol << endl;
    FILE_LOG(LOG_OPT) << " f0: " << cparams.f0(x, cparams.data) << endl;
    FILE_LOG(LOG_OPT) << " f1..m: " << fx << endl;
    FILE_LOG(LOG_OPT) << " g1..r: " << gx << endl;

    // compute the constraint gradients
    calc_Df(x, cparams, Df);
    calc_Dg(x, cparams, Dg);

    // compute the residual
    calc_residual(x, lambda, nu, t, cparams, Df, Dg, r);

    // compute primal-dual search direction
    if (solve_KKT)
      solve_KKT(x, lambda, nu, r, Df, Dg, cparams, dy);
    else
      solve_KKT_pd(x, lambda, nu, r, Df, Dg, cparams, dy);
    dy.get_sub_vec(0, n, dx);
    dy.get_sub_vec(n, mm+n, dlambda);
    dy.get_sub_vec(mm+n, mm+n+nu_len, dnu);
    const Real dy_norm = dy.norm();

    FILE_LOG(LOG_OPT) << " dy: " << dy << endl;
    FILE_LOG(LOG_OPT) << "   dx: " << dx << endl;
    FILE_LOG(LOG_OPT) << "   dlambda: " << dlambda << endl;
    FILE_LOG(LOG_OPT) << " -- doing line search" << endl;
    FILE_LOG(LOG_OPT) << "   old residual norm: " << r.norm() << endl;

    // do line search
    Real smax = 1.0;
    for (unsigned i=0; i< mm; i++)
      if (dlambda[i] < 0)
      {
        smax = std::min(smax, -lambda[i]/dlambda[i]);
        if (smax < std::numeric_limits<Real>::epsilon())
        {
          lambda[i] += NEAR_ZERO;
          FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_pd() - smax effectively zero" << endl << "  increasing component of lambda slightly, and trying again..." << endl;
          continue;
        }
      }
    Real s = 0.99*smax;
    (x_plus.copy_from(dx) *= s) += x;
    (lambda_plus.copy_from(dlambda) *= s) += lambda;
    (nu_plus.copy_from(dnu) *= s) += nu;

    // do the backtracking line search -- first, satisfy inequality constraints
    FILE_LOG(LOG_OPT) << " -- backtracking to satisfy inequality constraints..." << endl;
    unsigned start = 0;
    while (!feasible(cparams, x_plus, infeas_tol, start))
    {
      FILE_LOG(LOG_OPT) << "s: " << s << " x+: " << x_plus << endl;
      s *= cparams.beta;
      if (s*dy_norm < std::numeric_limits<Real>::epsilon())
      {
        FILE_LOG(LOG_OPT) << " -- x iterate is not feasible; terminating optimization" << endl;
        x = x_best;
        return false;
      } 

      // update f 
      (x_plus.copy_from(dx) *= s) += x;
    }

    // update lambda and nu
    lambda_plus.copy_from(dlambda) *= s;
    lambda_plus += lambda;
    nu_plus.copy_from(dnu) *= s;
    nu_plus += nu;

    // re-evaluate the objective function
    Real f = cparams.f0(x_plus, cparams.data);
    if (f < f0_best)
    {
      f0_best = f;
      x_best.copy_from(x_plus); 
    }

    // determine norm of r (target)
    const Real rnorm = r.norm();
    FILE_LOG(LOG_OPT) << " target residual norm: " << rnorm << endl;
    FILE_LOG(LOG_OPT) << " target residual norm * (1 - alpha*s): " << (1-cparams.alpha*s)*rnorm << endl;

    // compute new residual
    calc_Df(x_plus, cparams, Df);
    calc_Dg(x_plus, cparams, Dg);
    calc_residual(x_plus, lambda_plus, nu_plus, t, cparams, Df, Dg, rplus);

    // setup best s and best residual
    Real best_s, best_residual;
    if (rplus.norm() < rnorm)
    {
      best_s = s;
      best_residual = rplus.norm();
    }
    else
    {
      best_s = (Real) 0.0;
      best_residual = rnorm;
    }

    // continue backtracking line search until residual decreased
    while (rplus.norm() > (1-cparams.alpha*s)*rnorm) 
    {
      FILE_LOG(LOG_OPT) << "s: " << s << " x+: " << x_plus << endl;
      FILE_LOG(LOG_OPT) << "residual norm: " << rplus.norm() << endl;
      s *= cparams.beta;

      // verify that s is not too small
      if (s*dy_norm < std::numeric_limits<Real>::epsilon())
      {
        FILE_LOG(LOG_OPT) << "s too small; using best_s instead" << endl;
        if (best_s == (Real) 0.0)
        {
          FILE_LOG(LOG_OPT) << "  -- best s also to small; terminating optimization..." << endl;
          x.copy_from(x_best);
          return true;
        }        

        // otherwise, update x+, lambda+, nu+, and break
        (x_plus.copy_from(dx) *= best_s) += x;
        (lambda_plus.copy_from(dlambda) *= best_s) += lambda;
        (nu_plus.copy_from(dnu) *= best_s) += nu;
        break;
      }

      // update x+, lambda+, nu+
      x_plus.copy_from(dx) *= s;
      x_plus += x;
      lambda_plus.copy_from(dlambda) *= s;
      lambda_plus += lambda;
      nu_plus.copy_from(dnu) *= s;
      nu_plus += nu;

      // re-evaluate the objective function
      Real f = cparams.f0(x_plus, cparams.data);
      if (f < f0_best)
      {
        f0_best = f;
        x_best.copy_from(x_plus); 
      }

      // update the residual
      calc_Df(x_plus, cparams, Df);
      calc_Dg(x_plus, cparams, Dg);
      calc_residual(x_plus, lambda_plus, nu_plus, t, cparams, Df, Dg, rplus);
    }

    // get x, lambda, and nu out
    x.swap(x_plus);
    lambda.swap(lambda_plus);
    nu.swap(nu_plus); 
    r.swap(rplus);

    // see whether we can quit
    if (cparams.tcheck && cparams.tcheck(x, cparams.data))
    {
      FILE_LOG(LOG_OPT) << "  -- user specified termination" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::optimize_convex_pd() exited" << endl;
      return true;
    }

    // compute feasibility of constraints 
    evaluate_inequality_constraints(x, cparams, fc);
    evaluate_equality_constraints(x, cparams, gc);

    // redetermine infeasibility tolerance
    infeas_tol = (Real) 0.0;
    for (unsigned i=0; i< mm; i++)
    {
      if (fc[i] > (Real) 0.0)
      {
        infeas_tol = std::max(infeas_tol, fc[i]);
        fc[i] = (Real) 0.0;
      }
    }

    // copy fc to fx and gc to gx
    std::copy(fc.begin(), fc.end(), fx.begin());
    std::copy(gc.begin(), gc.end(), gx.begin());

    // add a little numerical clearance to the infeasibility tolerance, if
    // necessary
    if (infeas_tol > (Real) 0.0)
      infeas_tol += NEAR_ZERO;

    // redetermine eta
    eta = -(fx.dot(lambda));

    // compute rpri and rdual norms
    rdual_nrm = CBLAS::nrm2(n, r.begin(), 1);
    rpri_nrm = CBLAS::nrm2(nu_len, r.begin()+n+mm, 1); 

    FILE_LOG(LOG_OPT) << "   new residual norm: " << r.norm() << endl;
    FILE_LOG(LOG_OPT) << "   smax: " << smax << endl;
    FILE_LOG(LOG_OPT) << "   s: " << s << endl;
    FILE_LOG(LOG_OPT) << "   rdual norm: " << rdual_nrm << endl;
    FILE_LOG(LOG_OPT) << "   eta: " << eta << endl;

    if (rpri_nrm <= EPS_FEAS && rdual_nrm <= EPS_FEAS && eta < EPS)
    {
      FILE_LOG(LOG_OPT) << " -- tolerances satisfied; optimize_convex_pd() exiting" << endl;
      x = x_best;
      return true;
    }
  }

  // if we're here, the maximum number of iterations has been exceeded
  FILE_LOG(LOG_OPT) << " -- maximum number of iterations exceeded; optimize_convex_pd() exiting" << endl;
  x = x_best;
  return false;
}

/// Determines the residual for the primal-dual interior point method
void Optimization::calc_residual(const VectorN& x, const VectorN& lambda, const VectorN& nu, Real t, const OptParams& cparams, const MatrixN& Df, const MatrixN& Dg, VectorN& r)
{
  SAFESTATIC VectorN g0;
  SAFESTATIC VectorN DfTlambda, DgTnu;
  SAFESTATIC vector<Real> fc, gc;

  // get box constraints and linear inequality constraints
  const unsigned N_LBOX = cparams.lb.size();
  const unsigned N_UBOX = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();

  // get m and n and the length of nu
  const unsigned m = cparams.m;
  const unsigned n = cparams.n;
  const unsigned nu_len = cparams.r + cparams.A.rows();
  const unsigned mm = m + N_LBOX + N_UBOX + N_LININEQ;

  // setup r
  r.resize(mm+n+nu_len);

  // determine the gradient of the objective function
  cparams.grad0(x, g0, cparams.data);
  FILE_LOG(LOG_OPT) << " gradient at x: " << g0 << endl;

  // setup the dual residual
  Df.transpose_mult(lambda, DfTlambda);
  Dg.transpose_mult(nu, DgTnu);
  g0 += DfTlambda;
  g0 += DgTnu;
  r.set_sub_vec(0, g0);

  // evaluate all constraints
  evaluate_inequality_constraints(x, cparams, fc);

  // setup central residual components
  fc.resize(mm);
  for (unsigned i=0; i< mm; i++)
  {
    if (fc[i] >= (Real) 0.0)
      fc[i] = (Real) 0.0;
    fc[i] *= -lambda[i];
  }

  // subtract 1/t from each component of lambda*f and use it to populate
  // rcent
  std::transform(fc.begin(), fc.end(), r.begin()+n, std::bind2nd(std::minus<Real>(), (Real) 1.0/t));

  // setup the primary residual
  evaluate_equality_constraints(x, cparams, gc);
  std::copy(gc.begin(), gc.end(), r.begin()+mm+n);

  FILE_LOG(LOG_OPT) << " lambda: " << lambda << endl;
  FILE_LOG(LOG_OPT) << " Df: " << endl << Df;
  FILE_LOG(LOG_OPT) << " Dg: " << endl << Dg;
  FILE_LOG(LOG_OPT) << " rdual: " << r.get_sub_vec(0, n) << " norm: " << r.get_sub_vec(0, n).norm() << endl;
  FILE_LOG(LOG_OPT) << " rcent: " << r.get_sub_vec(n, mm+n) << " norm: " << r.get_sub_vec(n, mm+n).norm() << endl;
  FILE_LOG(LOG_OPT) << " rpri: " << r.get_sub_vec(mm+n, r.size()) << " norm: " << r.get_sub_vec(mm+n, r.size()).norm() << endl;
}

/// Computes the matrix Df (gradients of equality constraints)
void Optimization::calc_Dg(const VectorN& x, const OptParams& cparams, MatrixN& Dg)
{
  SAFESTATIC MatrixN J;

  // get m and nu length
  const unsigned n = cparams.n;
  const unsigned r = cparams.r;
  const unsigned N_LINEQ = cparams.b.size();
  const unsigned rr = r + N_LINEQ;

  // determine Dg
  Dg.set_zero(rr, n);

  // get Jacobian constraints for nonlinear inequality constraints
  if (r > 0)
  {
    cparams.cJac_g(x, J, cparams.data);
    Dg.set_sub_mat(0, 0, J);
  }
  
  // setup nonlinear equalities
  Dg.set_sub_mat(r, 0, cparams.A);
}

/// Computes the matrix Df (gradients of inequality constraints)
void Optimization::calc_Df(const VectorN& x, const OptParams& cparams, MatrixN& Df)
{
  SAFESTATIC MatrixN J;

  // get m and nu length
  const unsigned m = cparams.m;
  const unsigned n = cparams.n;
  const unsigned N_LBOX = cparams.lb.size();
  const unsigned N_UBOX = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();
  const unsigned mm = m + N_LBOX + N_UBOX + N_LININEQ;

  // determine Df -- linear inequality constraints first
  Df.set_zero(mm, n);
  Df.set_sub_mat(m+N_LBOX+N_UBOX, 0, cparams.M);
  Df.negate();

  // get Jacobian constraints for nonlinear inequality constraints
  if (m > 0)
  {
    cparams.cJac_f(x, J, cparams.data);
    Df.set_sub_mat(0, 0, J);
  }

  // determine Df -- lower box constraints now
  for (unsigned i=0, j=m; i< N_LBOX; i++, j++)
    Df(j, i) = (Real) -1.0;

  // determine Df -- upper box constraints now
  for (unsigned i=0, j=m+N_LBOX; i< N_UBOX; i++, j++)
    Df(j, i) = (Real) 1.0;
}

/// Solves a KKT system with equality constraints for primal dual method
/*
 * This is the default solver.
 */
/*
void Optimization::solve_KKT_pd(const VectorN& x, const VectorN& lambda, const VectorN& nu, const VectorN& rvec, const MatrixN& Df, const MatrixN& Dg, const OptParams& cparams, VectorN& dy)
{
  SAFESTATIC MatrixNN H, workM, M;
  SAFESTATIC VectorN l_fc, dx, dnu, workv;
  SAFESTATIC VectorN rhs, rcent;
  SAFESTATIC VectorN f_m_lambda, dlambda, Dfdx, fx, gx;
  SAFESTATIC vector<Real> fc, gc;

  // get m and nu length
  const unsigned m = cparams.m;
  const unsigned n = cparams.n;
  const unsigned r = cparams.r;
  const unsigned N_LINEQ = cparams.b.size();
  const unsigned nu_len = N_LINEQ + r;
  const unsigned N_LBOX = cparams.lb.size();
  const unsigned N_UBOX = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();
  const unsigned mm = m + N_LBOX + N_UBOX + N_LININEQ;

  // evaluate constraints
  evaluate_inequality_constraints(x, cparams, fc); 
  evaluate_equality_constraints(x, cparams, gc);

  // correct fc as necessary
  for (unsigned i=0; i< mm; i++)
    if (fc[i] >= -NEAR_ZERO)
      fc[i] = -NEAR_ZERO;

  // compute inverse of fc
  std::transform(fc.begin(), fc.end(), fc.begin(), std::bind1st(std::divides<Real>(), (Real) 1.0));

  // compute Hessian terms
  cparams.hess(x, (Real) 1.0, lambda, nu, H, cparams.data);
  FILE_LOG(LOG_OPT) << "H" << endl << H;

  // setup Hessian from gradient terms for nonlinear inequalities 
  for (unsigned i=0; i< m; i++)
  {
    Df.get_row(i, workv);
    *VectorN::outer_prod(workv,workv,&workM) *= (lambda[i]*fc[i]);
    H -= workM;
  }

  // setup Hessian from gradient terms for nonlinear equalities
//  for (unsigned i=0; i< r; i++)
//  {
//    Dg.get_row(i, workv);
//    *VectorN::outer_prod(workv,workv,&workM) *= (nu[i]*1.0/gc[i]);
//    H += workM;
//  }

  // setup Hessian components from gradients for lower box constraints
  for (unsigned i=0, j=m; i< N_LBOX; i++, j++)
    H(i,i) -= (lambda[j]*fc[j]);

  // setup Hessian components from gradients for upper box constraints
  for (unsigned i=0, j=m+N_LBOX; i< N_UBOX; i++, j++)
    H(i,i) -= (lambda[j]*fc[j]);

  // setup Hessian for linear inequality constraints
  for (unsigned i=0, j=m+N_LBOX+N_UBOX; i< N_LININEQ; i++, j++)
  {
    cparams.M.get_row(i, workv);
    *VectorN::outer_prod(workv,workv,&workM) *= (lambda[j]*fc[j]);
    H -= workM;
  }

  // create the KKT matrix
  M.resize(n+nu_len);
  M.set_sub_mat(0, 0, H);
  M.set_sub_mat(0, n, Dg, true);
  M.set_sub_mat(n, 0, Dg, false);

  // zero the lower right corner of the KKT matrix
  BlockIterator block = M.block_start(n,n+nu_len,n,n+nu_len);
  std::fill_n(block, nu_len*nu_len, (Real) 0.0);

  // copy fc to fx
  fx.resize(fc.size());
  std::copy(fc.begin(), fc.end(), fx.begin());

  // copy gc to gx
  gx.resize(gc.size());
  std::copy(gc.begin(), gc.end(), gx.begin());

  // setup the right hand side
  rhs.resize(n+nu_len);
  rvec.get_sub_vec(0, n, workv);
  rhs.set_sub_vec(0, workv);
  rvec.get_sub_vec(n+mm, n+mm+nu_len, workv);
  rhs.set_sub_vec(n, workv);
  rhs.negate();

  FILE_LOG(LOG_OPT) << "** solving KKT equation **" << endl;
  FILE_LOG(LOG_OPT) << "M: " << endl << M;
  FILE_LOG(LOG_OPT) << "rhs: " << rhs << endl;

  // two cases: nu == 0 (no equality constraints), nu > 0 (equality constraints)
  if (nu_len == 0)
  {
    try
    {
      // Hessian *may* be positive-definite; try using Cholesky solver first
      LinAlg::solve_SPD_fast(M, rhs);
      rhs.get_sub_vec(0, n, dx);
    }
    catch (SingularException e)
    {
      // matrix is only PSD; retry using least squares solver
      // recreate the KKT matrix
      M.copy_from(H);

      // solve using robust least squares solver
      LinAlg::solve_LS_fast(M, rhs);
      rhs.get_sub_vec(0, n, dx);
    }
  }
  else
  {
    try
    {
      // try solving w/ symmetric LDL' solver ([Boyd 2004], p. 610, eq. 11.55)
      LinAlg::solve_symmetric_fast(M, rhs); 
      rhs.get_sub_vec(0, n, dx);
      rhs.get_sub_vec(n, n+nu_len, dnu);
    }
    catch (SingularException e)
    {    
      // re-create the KKT matrix
      M.set_sub_mat(0, 0, H);
      M.set_sub_mat(0, n, Dg, true);
      M.set_sub_mat(n, 0, Dg, false);

      // re-zero the lower right corner of the KKT matrix
      std::fill_n(block, nu_len*nu_len, (Real) 0.0);

      // solve using robust least squares solver
      LinAlg::solve_LS_fast(M, rhs);
      rhs.get_sub_vec(0, n, dx);
      rhs.get_sub_vec(n, n+nu_len, dnu);
    }
  }
  
  // solve for dlambda
  rvec.get_sub_vec(n, n+mm, rcent);
  f_m_lambda.resize(fc.size());
  std::transform(fc.begin(), fc.end(), lambda.begin(), f_m_lambda.begin(), std::multiplies<Real>());
  MatrixNN::diag_mult(f_m_lambda, Df.mult(dx,Dfdx), dlambda);
  dlambda.negate();
  dlambda += MatrixNN::diag_mult(fx, rcent, workv);

  FILE_LOG(LOG_OPT) << " Df: " << endl << Df;
  FILE_LOG(LOG_OPT) << " 1/f: " << VectorN(fc.begin(), fc.end()) << endl;
  FILE_LOG(LOG_OPT) << " lambda: " << lambda << endl;
  FILE_LOG(LOG_OPT) << " Hpd: " << endl << H;
  FILE_LOG(LOG_OPT) << " r: " << rvec << endl;
  FILE_LOG(LOG_OPT) << " dx: " << dx << endl;
  FILE_LOG(LOG_OPT) << " dnu: " << dnu << endl;
  FILE_LOG(LOG_OPT) << " dlambda: " << dlambda << endl;
  FILE_LOG(LOG_OPT) << "** KKT equation solved" << endl;

  // solve
  dy.resize(dx.size() + dlambda.size() + dnu.size());
  dy.set_sub_vec(0, dx);
  dy.set_sub_vec(n, dlambda);
  dy.set_sub_vec(mm+n, dnu);
} 
*/

/// Solves a KKT system with equality constraints for primal dual method
/*
 * This is the default solver.
 * This is based on the following KKT conditions:
 * 1. grad(f0(x)) + sum_{i=1}^{mm} lambda_i*f_i(x) + sum_{i=1}^{rr} nu_i*g_i(x) 
 * 2. -lambda_i * f_i(x) = 1/t (complementarity constraint), for i = 1 .. mm
 * 3. g_i(x) = 0, for i = 1 .. rr
 * 4. lambda >= 0
 * 5. f_i(x) < 0, for i = 1 .. mm
 * Condition #1 is the dual residual, #2 is the centrality residual, and #3
 * is the primary residual. Conditions #4 and #5 are handled by feasibility.
 * d#1/dx = Hess(f0(x)) + Df^T*lambda + Dg^T*nu
 * d#1/dlambda = Df^T 
 * d#1/dnu = Dg^T
 * d#2/dx = -diag(lambda) * Df 
 * d#2/dlambda = -diag(f)
 * d#2/dnu = 0
 * d#3/dx = Dg 
 * d#3/dlambda = 0 
 * d#3/dnu = 0 
 */
void Optimization::solve_KKT_pd(const VectorN& x, const VectorN& lambda, const VectorN& nu, const VectorN& rvec, const MatrixN& Df, const MatrixN& Dg, const OptParams& cparams, VectorN& dy)
{
  SAFESTATIC MatrixN H, M;
  SAFESTATIC VectorN nDfTlambda;
  SAFESTATIC vector<Real> fc;

  // get m and nu length
  const unsigned m = cparams.m;
  const unsigned n = cparams.n;
  const unsigned r = cparams.r;
  const unsigned N_LINEQ = cparams.b.size();
  const unsigned nu_len = N_LINEQ + r;
  const unsigned N_LBOX = cparams.lb.size();
  const unsigned N_UBOX = cparams.ub.size();
  const unsigned N_LININEQ = cparams.q.size();
  const unsigned mm = m + N_LBOX + N_UBOX + N_LININEQ;

  // compute Hessian terms
  cparams.hess(x, (Real) 1.0, lambda, nu, H, cparams.data);
  FILE_LOG(LOG_OPT) << "H" << endl << H;

  // compute some stuff we'll need:
  evaluate_inequality_constraints(x, cparams, fc);


  // 0. rdual
//  cparams.grad0(x, rdual, cparams.data);
//  rdual += Df.transpose_mult(fx, workv);
//  rdual += Dg.transpose_mult(gx, workv);

  // 1. -lambda_i * f_i - 1/t
//  rcent.resize(fc.size());
//  std::transform(fc.begin(), fc.end(), lambda.begin(), rcent.begin(), std::multiplies<Real>());
//  std::transform(rcent.begin(), rcent.end(), rcent.begin(), std::bind2nd(std::plus<Real>(), (Real) 1.0/t));
//  std::transform(rcent.begin(), rcent.end(), rcent.begin(), std::negate<Real>());

  // 2. g_i 
//  rpri.resize(gc.size());
//  std::transform(gc.begin(), gc.end(), rpri.begin(), std::bind2nd(std::minus<Real>(), (Real) 0.0));

  // compute -Df' * lambda
  Df.transpose_mult(lambda, nDfTlambda);
  nDfTlambda.negate();

  // create the KKT matrix
  M.resize(n+mm+nu_len,n+mm+nu_len);

  // setup the appropriate parts of the KKT matrix
  M.set_sub_mat(0, 0, H);
  M.set_sub_mat(0, n, Df, true);
  M.set_sub_mat(0, n+mm, Dg, true);
  M.set_sub_mat(n, 0, nDfTlambda, true);
  M.set_sub_mat(n+mm, 0, Dg);

  // zero out appropriate block of the KKT matrix
  BlockIterator bi = M.block_start(n, n+mm+nu_len, n, n+mm+nu_len); 
  std::fill_n(bi, (mm+nu_len)*(mm+nu_len), (Real) 0.0);

  // setup diagonal
  for (unsigned i=0, j=n; i< mm; i++, j++)
    M(j, j) = -fc[i]; 

  // create the right hand side
  dy.copy_from(rvec).negate();

  FILE_LOG(LOG_OPT) << "** solving KKT equation **" << endl;
  FILE_LOG(LOG_OPT) << "M: " << endl << M;
  FILE_LOG(LOG_OPT) << "rhs: " << dy << endl;

  // solve using stanard solver
  try
  {
    LinAlg::solve_fast(M, dy);
  }
  catch (SingularException e)
  {
    // reform the KKT matrix
    M.set_sub_mat(0, 0, H);
    M.set_sub_mat(0, n, Df, true);
    M.set_sub_mat(0, n+mm, Dg, true);
    M.set_sub_mat(n, 0, nDfTlambda, true);
    M.set_sub_mat(n+mm, 0, Dg);

    // zero out appropriate block of the KKT matrix
    BlockIterator bi = M.block_start(n, n+mm+nu_len, n, n+mm+nu_len); 
    std::fill_n(bi, (mm+nu_len)*(mm+nu_len), (Real) 0.0);

    // setup diagonal
    for (unsigned i=0, j=n; i< mm; i++, j++)
      M(j, j) = -fc[i]; 

    // use the least squares solver 
    try
    {
      LinAlg::solve_LS_fast1(M, dy);
    }
    catch (NumericalException e)
    {
      // reform the KKT matrix
      M.set_sub_mat(0, 0, H);
      M.set_sub_mat(0, n, Df, true);
      M.set_sub_mat(0, n+mm, Dg, true);
      M.set_sub_mat(n, 0, nDfTlambda, true);
      M.set_sub_mat(n+mm, 0, Dg);

      // zero out appropriate block of the KKT matrix
      BlockIterator bi = M.block_start(n, n+mm+nu_len, n, n+mm+nu_len); 
      std::fill_n(bi, (mm+nu_len)*(mm+nu_len), (Real) 0.0);

      // setup diagonal
      for (unsigned i=0, j=n; i< mm; i++, j++)
        M(j, j) = -fc[i]; 

      // last resort...
      LinAlg::solve_LS_fast2(M, dy);
    }
  }

  FILE_LOG(LOG_OPT) << " r: " << rvec << endl;
  FILE_LOG(LOG_OPT) << " dx: " << dy.get_sub_vec(0, n) << endl;
  FILE_LOG(LOG_OPT) << " dnu: " << dy.get_sub_vec(n+mm, dy.size()) << endl;
  FILE_LOG(LOG_OPT) << " dlambda: " << dy.get_sub_vec(n, n+mm) << endl;
  FILE_LOG(LOG_OPT) << "** KKT equation solved" << endl;
}

/// Solves a mixed linear complementarity problem
/**
 * Solves the MLCP: M11y + M12z + q1 = 0
 *                  M21y + M22z + q2 = w
 *                  w >= 0
 *                  z >= 0
 *                  z'w = 0
 * \param y the vector of unconstrained variables is returned here; the size
 *        of y must be set correctly before calling
 * \param z the vector of complementarity constrained variables is returned 
 *        here; the size of z must be set correctly before calling
 * \param M the MLCP matrix
 * \param q the MLCP vector
 * \return <b>false</b> if underlying LCP solver fails; <b>true</b> otherwise
 * \note you must check that the solutions are valid unless the MLCP is convex!
 * \note sizes of y and z on entry determines the number of 
 *       unconstrained/constrained variables 
 */
bool Optimization::mlcp(VectorN& y, VectorN& z, const MatrixN& M, const VectorN& q)
{
  MatrixN M11, M22, M12, M21, iM11, MM, workM;
  VectorN q1, q2, q1bar, qq;

  // setup data to go to the interior point method
  unsigned nfree = y.size();
  unsigned n = z.size();

  // verify that nfree + n = q.size()
  if (nfree + n != q.size())
    throw std::runtime_error("len(y) + len(z) != len(q)");

  // partition M
  M.get_sub_mat(0, nfree, 0, nfree, M11);
  M.get_sub_mat(0, nfree, nfree, n+nfree, M12);
  M.get_sub_mat(nfree, n+nfree, 0, nfree, M21);
  M.get_sub_mat(nfree, n+nfree, nfree, n+nfree, M22);

  // partition q
  q.get_sub_vec(0, nfree, q1);
  q.get_sub_vec(nfree, n+nfree, q2);

  // compute pseudo-inverse of M11 
  iM11.copy_from(M11);
  LinAlg::pseudo_inverse(iM11, LinAlg::svd1);

  // setup bar_M and b_q
  iM11.mult(M12, workM);
  M21.mult(workM, MM);
  MM.negate();
  MM += M22;
  iM11.mult(q1, q1bar);
  M21.mult(q1bar, qq);
  qq.negate();
  qq += q2;

  // init z to something sensible
  z.set_zero();

  FILE_LOG(LOG_OPT) << "Optimization::solve_MLCP() entered" << endl;
  FILE_LOG(LOG_OPT) << "MM: " << endl << MM;
  FILE_LOG(LOG_OPT) << "qq: " << qq << endl;
  FILE_LOG(LOG_OPT) << "z: " << z << endl;

  // call Lemke's algorithm
  bool result = lcp_lemke(MM, qq, z);
  if (!result)
  {
    FILE_LOG(LOG_OPT) << "-- lemke's algorithm failed!" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::solve_MLCP() exited" << endl;
    return false;
  }

  // determine y
  workM.mult(z, y) += q1bar;
  y.negate(); 

  FILE_LOG(LOG_OPT) << "-- lemke's algorithm successful" << endl;
  if (LOGGING(LOG_OPT))
    workM.copy_from(M11);
  FILE_LOG(LOG_OPT) << "rank(M11): " << LinAlg::calc_rank(workM) << " size M11: " << M11.rows() << endl;
  FILE_LOG(LOG_OPT) << "M11*y + M12*z + q1: " << (M11*y + M12*z + q1) << endl;
  FILE_LOG(LOG_OPT) << "z: " << z << endl;
  FILE_LOG(LOG_OPT) << "y: " << y << endl;
  FILE_LOG(LOG_OPT) << "equality constraint evaluation: " << (M11*y + M12*z + q1).norm() << endl;

  FILE_LOG(LOG_OPT) << "Optimization::solve_MLCP() exited" << endl;

  // indicate that LCP solver was successful 
  return true;
}

/// Determines whether the constraint functions are feasible
bool Optimization::feasible(const OptParams& oparams, const VectorN& x, Real infeas_tol, unsigned& start)
{
  SAFESTATIC VectorN work;

  const unsigned n = x.size();
  const unsigned N_LB = oparams.lb.size();
  const unsigned N_UB = oparams.ub.size();
  const unsigned N_LININEQ = oparams.q.size();

  // check lower box constraints first
  if (start < N_LB)
    for (unsigned i=start; i< N_LB; i++, start++)
      if (oparams.lb[i] - x[i] > infeas_tol)
        return false;

  // check upper box constraints next 
  if (start < N_LB + N_UB)
    for (unsigned i=start - N_LB; i< N_UB; i++, start++)
      if (x[i] - oparams.ub[i] > infeas_tol)
        return false;

  // check linear inequality constraints next
  if (start < N_LB + N_UB + N_LININEQ)
  {
    oparams.M.mult(x, work) -= oparams.q;
    work.negate();
    for (unsigned i=start - N_LB - N_UB; i< work.size(); i++, start++)
      if (work[i] > infeas_tol)
        return false;
  }

  // check nonlinear equality constraints next
  if (oparams.r > 0 && oparams.gx)
  {
    oparams.gx(x, work, oparams.data);
    pair<Real*, Real*> mmax = boost::minmax_element(work.begin(), work.end());
    if (mmax.first != work.end() && std::max(std::fabs(*mmax.first), *mmax.second) > oparams.eps_feas)
      return false;
  }

  // check nonlinear inequality constraints last
  if (oparams.m > 0 && oparams.fx)
  {
    oparams.fx(x, work, oparams.data);
    Real *maxele = std::max_element(work.begin(), work.end());
    if (maxele != work.end() && *maxele > infeas_tol)
      return false;
  }

  return true;
}

/// Computes the gradient of a function numerically
VectorN Optimization::ngrad(const VectorN& x, Real t, Real h, void* data, Real (*ofn)(const VectorN&, Real, void*))
{
  unsigned n = x.size();
  const Real inv_h2 = 0.5 / h;
  VectorN grad;

  // make a copy of x
  VectorN xx = x;

  grad.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    xx[i] += h;
    Real v1 = (*ofn)(xx, t, data);
    xx[i] -= 2*h;
    Real v2 = (*ofn)(xx, t, data);
    xx[i] += h;
    grad[i] = (v1 - v2)*inv_h2;
  }

  return grad;
}

/// Computes the hessian of a function numerically
MatrixN Optimization::nhess(const VectorN& x, Real t, Real h, void* data, Real (*ofn)(const VectorN&, Real, void*))
{
  unsigned n = x.size();
  const Real inv_h2 = 0.5 / h;
  MatrixN hess;

  // make a copy of x
  VectorN xx = x;

  hess.resize(n,n);
  for (unsigned i=0; i< n; i++)
  {
    xx[i] += h;
    VectorN v1 = ngrad(xx, t, h, data, ofn);
    xx[i] -= 2*h;
    VectorN v2 = ngrad(xx, t, h, data, ofn);
    xx[i] += h;
    hess.set_column(i, (v1 - v2)*inv_h2);
  }

  // average values of the Hessian
  for (unsigned i=0; i< n; i++)
    for (unsigned j=i+1; j< n; j++)
      hess(i,j) = hess(j,i) = 0.5*(hess(i,j) + hess(j,i));

  return hess;
}

/// Makes a (possibly) infinite value finite again
/**
 * Converts values of -inf to -DBL_MAX and inf to DBL_MAX
 * \note utility function for lp()
 */
Real Optimization::finitize(Real x)
{
  if (x == std::numeric_limits<Real>::infinity())
    return std::numeric_limits<Real>::max();
  else if (x == -std::numeric_limits<Real>::infinity())
    return -std::numeric_limits<Real>::max();
  else
    return x;
}

/// Solves a linear program using the method of Seidel
/**
 * This method exhibits complexity of O(d!n), where d is the dimension of the
 * variables and n is the number of constraints.
 * \param A the matrix for which Ax < b
 * \param b the vector for which Ax < b
 * \param c the optimization vector (maximizes c'x)
 * \param l the lower variable constraints on x
 * \param u the upper variable constraints on x
 * \param x the optimal solution (on successful return)
 * \return <b>true</b> if successful, <b>false</b> otherwise
 * \note using limits of +/- inf can result in overflow with this algorithm
 *       and is not recommended; use lower limits
 */
bool Optimization::lp(const MatrixN& A, const VectorN& b, const VectorN& c, const VectorN& l, const VectorN& u, VectorN& x)
{
  // get number of rows and columns in A
  unsigned n = A.rows();
  unsigned d = A.columns();

  FILE_LOG(LOG_OPT) << "Optimization::lp() entered" << endl;
  FILE_LOG(LOG_OPT) << "A: " << endl << A;
  FILE_LOG(LOG_OPT) << "b: " << b << endl;
  FILE_LOG(LOG_OPT) << "c: " << c << endl;
  FILE_LOG(LOG_OPT) << "l: " << l << endl;
  FILE_LOG(LOG_OPT) << "u: " << u << endl;

  // base case d = 1
  if (d == 1)
  {
    FILE_LOG(LOG_OPT) << "base case, d = 1" << endl;

    Real high = u[0]; 
    Real low = l[0];

    for (unsigned i=0; i< n; i++)
    {
      if (A(i,0) > std::numeric_limits<Real>::epsilon())
        high = std::min(high, b[i]/A(i,0));
      else if (A(i,0) < -std::numeric_limits<Real>::epsilon())
        low = std::max(low, b[i]/A(i,0));
      else if (b[i] < -std::numeric_limits<Real>::epsilon())
      {
        FILE_LOG(LOG_OPT) << "infeasible; b is negative and A is zero" << endl;

        return false; 
      }
    }

    // do a check for infeasibility
    if (high < low)
    {
      FILE_LOG(LOG_OPT) << "infeasible; high (" << high << ") < low (" << low << ")" << endl;

      return false; 
    }
  
    // set x
    x = VectorN(1);
    x[0] = (c[0] >= 0.0) ? high : low;

    FILE_LOG(LOG_OPT) << "optimal 1D x=" << x << endl;
    FILE_LOG(LOG_OPT) << "Optimization::lp() exited" << endl;

    // otherwise, good return
    return true; 
  }

  // pick a random shuffle for A and b
  vector<unsigned> permut(n);
  for (unsigned i=0; i< n; i++)
    permut[i] = i;
  std::random_shuffle(permut.begin(), permut.end());

  // setup aa and bb
  vector<VectorN> aa;
  VectorN bb;
  aa.resize(n);
  bb.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    aa[i] = A.get_row(permut[i]);
    bb[i] = b[permut[i]];
  }

  FILE_LOG(LOG_OPT) << "A (permuted): " << endl;
  for (unsigned i=0; i< n; i++)
    FILE_LOG(LOG_OPT) << aa[i] << endl;
  FILE_LOG(LOG_OPT) << "b (permuted): " << bb << endl;

  // setup optimum vector
  x.resize(d);
  for (unsigned i=0; i< d; i++)
    if (c[i] > 0.0)
        x[i] = u[i];
    else if (c[i] < 0.0) 
        x[i] = l[i];
    else
        x[i] = (std::fabs(l[i]) < std::fabs(u[i])) ? l[i] : u[i];

  FILE_LOG(LOG_OPT) << "initial x=" << x << endl;

  // process half-space constraints 
  for (unsigned i=0; i< n; i++)
  {
    FILE_LOG(LOG_OPT) << "-- processing halfspace constraint " << i << endl;

    // if x respects new halfspace constraint, nothing else to do..
    Real val = bb[i] - VectorN::dot(aa[i], x);
    if (val >= -NEAR_ZERO)
    {    
      FILE_LOG(LOG_OPT) << "  ** constraint already satisfied!" << endl;
      continue;
    }

    FILE_LOG(LOG_OPT) << "  -- constraint not satisfied (" << val << "); solving recursively" << endl;

    // search for maximum value in the a vector
    unsigned k = std::numeric_limits<unsigned>::max();
    Real maximal = -std::numeric_limits<Real>::max();
    for (unsigned j=0; j< d; j++)
      if (std::fabs(aa[i][j]) > maximal && aa[i][j] != 0.0)
      {
        maximal = std::fabs(aa[i][j]);
        k = j;
      }

    FILE_LOG(LOG_OPT) << "  -- determined k: " << k << endl;

    // look for infeasibility
    if (k == std::numeric_limits<unsigned>::max())
    {
      FILE_LOG(LOG_OPT) << "  -- k is infeasible; problem is infeasible" << endl;

      return false; 
    }

    // setup useful vector and constant
    VectorN aak = aa[i]/aa[i][k];
    Real bak = bb[i]/aa[i][k];

    FILE_LOG(LOG_OPT) << "  -- vector a/a(k): " << aak << " " << bak << endl;

    // copy vectors aa and bb
    vector<VectorN> aac;
    vector<Real> bbc;
    aac.resize(i);
    bbc.resize(i);
    for (unsigned j=0; j< i; j++)
    {
      aac[j].copy_from(aa[j]);
      bbc[j] = bb[j];
    }

    // modify copy of vector aa
    for (unsigned j=0; j< i; j++)
    {
      aac[j] -= aak*aa[j][k];
      assert(std::fabs(aac[j][k]) < NEAR_ZERO);
      aac[j] = remove_component(aac[j], k);
    }

    // modify copy of vector bb
    for (unsigned j=0; j< i; j++)
    {
      bbc[j] -= bak * aa[j][k];
      if (std::isinf(bbc[j]))
        bbc[j] = finitize(bbc[j]);
    }

    // modify copy of c
    assert(std::fabs((c - aak*c[k])[k]) < NEAR_ZERO);
    VectorN cprime = remove_component(c - aak*c[k], k);  

    // generate new lower and upper bounds for variables
    VectorN lprime = remove_component(l, k);
    VectorN uprime = remove_component(u, k);

    // setup two new constraints 
    VectorN f, g;
    f.set_zero(d);
    g.set_zero(d);
    f[k] = 1;
    g[k] = -1;
    assert(std::fabs((f - aak)[k]) < NEAR_ZERO);
    assert(std::fabs((g + aak)[k]) < NEAR_ZERO);
    f = remove_component(f - aak, k);
    g = remove_component(g + aak, k);
    Real bf = u[k] - bak;
    Real bg = -l[k] + bak;
    if (std::isinf(bf))
      bf = finitize(bf);
    if (std::isinf(bg))
      bg = finitize(bg);
    aac.push_back(f);
    bbc.push_back(bf);
    aac.push_back(g);
    bbc.push_back(bg);

    // create the Aprime matrix from aac
    MatrixN Aprime;
    Aprime.resize(aac.size(), d-1);  
    for (unsigned j=0; j< aac.size(); j++)
      Aprime.set_row(j, aac[j]);

    // create the bprime vector from bbc
    VectorN bprime;
    bprime.resize(bbc.size());
    for (unsigned j=0; j< bbc.size(); j++)
      bprime[j] = bbc[j];

    FILE_LOG(LOG_OPT) << "  -- A': " << endl << Aprime;
    FILE_LOG(LOG_OPT) << "  -- b': " << bprime << endl;
    FILE_LOG(LOG_OPT) << "  -- c': " << cprime << endl;
    FILE_LOG(LOG_OPT) << "  -- u': " << uprime << endl;
    FILE_LOG(LOG_OPT) << "  -- l': " << lprime << endl;
    FILE_LOG(LOG_OPT) << "  -- f: " << f << " " << bf << endl;
    FILE_LOG(LOG_OPT) << "  -- g: " << g << " " << bg << endl;
    FILE_LOG(LOG_OPT) << " + solving recursive subproblem" << endl;

    // solve the (d-1)-dimensional problem and ``lift'' the solution
    if (!lp(Aprime,bprime,cprime,lprime,uprime,x))
      return false;

    FILE_LOG(LOG_OPT) << "  -- recursively determined x: " << x << endl;
    FILE_LOG(LOG_OPT) << "  -- k: " << k << endl;

    // insert a zero into the k'th dimension of x
    x = insert_component(x, k);
    FILE_LOG(LOG_OPT) << "  -- x w/inserted component at k: " << x << endl;

    // solve for the proper k'th value of x
    x[k] = (bb[i] - VectorN::dot(aa[i], x))/aa[i][k];
    FILE_LOG(LOG_OPT) << "  -- optimal x (to this point): " << x << endl;

    // verify that half-plane constraints still met
    for (unsigned j=0; j<= i; j++)
      assert(VectorN::dot(aa[j], x) - bb[j] <= NEAR_ZERO);

    // verify that lower constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] >= l[j] - NEAR_ZERO);

    // verify that upper constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] <= u[j] + NEAR_ZERO);
  }

  // verify that half-plane constraints still met
  FILE_LOG(LOG_OPT) << "b - A*x: " << (b - A*x) << endl;
  for (unsigned j=0; j < n; j++)
    assert(VectorN::dot(aa[j], x) - bb[j] <= NEAR_ZERO);

  // verify that lower constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] >= l[i] - NEAR_ZERO);

  // verify that upper constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] <= u[i] + NEAR_ZERO);
  FILE_LOG(LOG_OPT) << "all halfspace constraints satisfied; optimum found!" << endl;
  FILE_LOG(LOG_OPT) << "optimum = " << x << endl;
  FILE_LOG(LOG_OPT) << "Optimization::lp() exited" << endl;

  return true; 
}

/// Inserts a component (value will be zero) at the k'th position in the given vector and returns a new vector
VectorN Optimization::insert_component(const VectorN& x, unsigned k)
{
  VectorN xn;
  xn.resize(x.size()+1);

  if (k == 0)
    xn.set_sub_vec(1, x);
  else if (k == x.size())
    xn.set_sub_vec(0, x);
  else
  {
    xn.set_sub_vec(0,x.get_sub_vec(0,k));
    xn.set_sub_vec(k+1,x.get_sub_vec(k,x.size()));
  }  
  xn[k] = 0;
  
  return xn;
}

/// Gets a subvector with the k'th component removed
/**
 * \note for use by lp()
 */
VectorN Optimization::remove_component(const VectorN& v, unsigned k)
{
  const unsigned d = v.size();

  if (k == 0)
    return v.get_sub_vec(1,d);
  else if (k == d-1)
    return v.get_sub_vec(0,d-1);
  else
  {
    // setup v w/component k removed
    VectorN vpm1, vpp1, vp;
    vp.resize(d-1);    
    v.get_sub_vec(0,k, vpm1);
    v.get_sub_vec(k+1,d, vpp1);
    vp.set_sub_vec(0,vpm1);
    vp.set_sub_vec(k,vpp1);
    return vp;
  }
}

/// Finds a Cauchy point for gradient projection QP method
VectorN& Optimization::find_cauchy_point(const MatrixN& G, const VectorN& c, const VectorN& l, const VectorN& u, const VectorN& gradient, VectorN& x)
{
  const Real INF = std::numeric_limits<Real>::max();

  if (G.rows() != G.columns())
    throw NonsquareMatrixException();

  // get size of gradient / x
  const unsigned n = x.size();
  assert(gradient.size() == n);
  assert(l.size() == n);
  assert(u.size() == n);

  // add the value 0 to t
  vector<Real> t;
  t.push_back((Real) 0.0);

  FILE_LOG(LOG_OPT) << "Optimization::find_cauchy_point() entered" << endl;
  FILE_LOG(LOG_OPT) << "  gradient: " << gradient << endl;
  FILE_LOG(LOG_OPT) << "  l: " << l << endl;
  FILE_LOG(LOG_OPT) << "  u: " << u << endl;
  FILE_LOG(LOG_OPT) << "  x: " << x << endl;

  // identify values of t for which each component reaches its bound along -g
  for (unsigned i=0; i< n; i++)
  {
    if (gradient[i] < (Real) 0.0 && u[i] < INF)
    {
      Real value = (x[i] - u[i])/gradient[i];
      if (value > (Real) 0.0)
        t.push_back(value);
    }
    else if (gradient[i] > (Real) 0.0 && l[i] > -INF)
    {
      Real value = (x[i] - l[i])/gradient[i];
      if (value > (Real) 0.0)
        t.push_back(value);
    }
    else
      t.push_back(INF);
  }

  // sort t and remove duplicate values
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());

  if (LOGGING(LOG_OPT))
  {
    std::ostringstream sstr;
    sstr << "  t: ";
    for (unsigned i=0; i< t.size(); i++)
      sstr << t[i] << " ";
    FILE_LOG(LOG_OPT) << sstr.str() << endl;
  } 
  FILE_LOG(LOG_OPT) << "  searching for local minimizer..." << endl;

  // determine first local minimizer
  SAFESTATIC VectorN p, Gp;
  p.resize(n);
  for (unsigned j=1; j< t.size(); j++)
  {
    unsigned jm1 = j-1;

    // compute p^{j-1}
    for (unsigned i=0; i< n; i++)
    {
      // (re)compute ti
      Real ti;
      if (gradient[i] < (Real) 0.0 && u[i] < INF)
        ti = (x[i] - u[i])/gradient[i];
      else if (gradient[i] > (Real) 0.0 && l[i] > -INF)
        ti = (x[i] - l[i])/gradient[i];
      else
        ti = INF;

      p[i] = (t[jm1] < ti) ? -gradient[i] : (Real) 0.0;
    }

    FILE_LOG(LOG_OPT) << "    search interval: [" << t[jm1] << ", " << t[j] << "]" << endl;

    // look for local minimizer
    Real fprime = c.dot(p) + x.dot(G.mult(p, Gp));

    FILE_LOG(LOG_OPT) << "    f': " << fprime << endl;
    // if f' > 0 there is a local minimizer at t[j-1]
    if (fprime > 0)
    {
      // there is a local minimizer of 0.5x'Gx + x'c at t = t[j-1]
      FILE_LOG(LOG_OPT) << "    local minimizer found!" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::find_cauchy_point() exiting" << endl;
      return x;
    }
    else
    {
      // compute dt*
      Real fpprime = p.dot(Gp);
      Real dt_star = -fprime / fpprime;
      FILE_LOG(LOG_OPT) << "    dt*: " << dt_star << endl;
      if (dt_star >= (Real) 0.0 && dt_star < t[j] - t[jm1])
      {
        FILE_LOG(LOG_OPT) << "    dt* in [0, t_j - t_{j-1}]" << endl;
        FILE_LOG(LOG_OPT) << "    local minimizer found!" << endl;
        FILE_LOG(LOG_OPT) << "Optimization::find_cauchy_point() exiting" << endl;

        // there is a local minimizer at t[jm1] + dt_star
        x += p*dt_star;

        // manually enforce constraints
        for (unsigned i=0; i< x.size(); i++)
          if (x[i] < l[i])
            x[i] = l[i];
          else if (x[i] > u[i])
            x[i] = u[i];

        return x;
      }
      else
      {
        // update x and keep looking
        Real dt = t[j] - t[jm1];
        if (!std::isinf(dt) && !std::isnan(dt))
          x += p*dt;

        FILE_LOG(LOG_OPT) << "    updating x: " << x << endl;
        FILE_LOG(LOG_OPT) << "    continuing search..." << endl;
      }
    }
  }

  FILE_LOG(LOG_OPT) << "  local minimizer must be last point examined" << endl;
  FILE_LOG(LOG_OPT) << "Optimization::find_cauchy_point() exiting" << endl;

  // manually enforce constraints
  for (unsigned i=0; i< x.size(); i++)
    if (x[i] < l[i])
      x[i] = l[i];
    else if (x[i] > u[i])
      x[i] = u[i];

  return x;
}

/// Minimizes 0.5x'Gx + x'c with box constraints
/**
 * Gradient projection method for quadratic programming with box constraints.
 * Convergence is guaranteed.
 * \param G a n x n symmetric, positive semi-definite matrix
 * \param c a n-dimensional vector
 * \param l a n-dimensional vector of lower limit constraints
 * \param u a n-dimensional vector of upper limit constraints
 * \param max_iter the maximum number of iterations
 * \param x the vector at which the minimum occurs
 * \return the number of iterations
 */
unsigned Optimization::qp_gradproj(const MatrixN& G, const VectorN& c, const VectorN& l, const VectorN& u, unsigned max_iter, VectorN& x, Real tol)
{
  if (G.rows() != G.columns())
    throw NonsquareMatrixException();

  const Real INF = std::numeric_limits<Real>::max();

  // minimum objective function value determined so far
  Real qmin = INF;

  // setup the number of inner iterations
  unsigned niter = 0;

  // determine n
  const unsigned n = c.size();

  // pick starting x halfway between l and u, if necessary
  for (unsigned i=0; i< n; i++)
    if (x[i] < l[i] || x[i] > u[i])
      x[i] = l[i]*0.5 + u[i]*0.5;

  // setup some work variables
  VectorN Gx_c, Gx, H, Wzz;
  MatrixN ZtGZ;
  VectorN xz, cz, dz, lz, uz, dxz;
  VectorN y, r1, r2, w, w1, w2, v, xstar;
  vector<bool> inactive_set;

  // inactive set
  inactive_set.resize(n);

  // the preconditioner
  H.resize(n);
  for (unsigned i=0; i< n; i++)
    H[i] = (G(i,i) > (Real) 0.0) ? (Real) 1.0/G(i,i) : (Real) 1.0;

  FILE_LOG(LOG_OPT) << "qp_gradproj() entered" << endl;
  FILE_LOG(LOG_OPT) << "  G: " << endl << G;
  FILE_LOG(LOG_OPT) << "  c: " << c << endl;
  FILE_LOG(LOG_OPT) << "  l: " << l << endl;
  FILE_LOG(LOG_OPT) << "  u: " << u << endl;
  FILE_LOG(LOG_OPT) << "  initial x: " << x << endl;
  FILE_LOG(LOG_OPT) << "  preconditioner (diagonal matrix): " << H << endl;

  // iterate..
  for (; niter < max_iter; niter++)
  {
    // compute the gradient -- if it is zero, quit now
    G.mult(x, Gx_c) += c;
    if (Gx_c.norm() < tol)
      break;
    
    // find Cauchy point
    find_cauchy_point(G, c, l, u, Gx_c, x);

    // evaluate Cauchy point
    Real q = x.dot(G.mult(x, Gx)) * (Real) 0.5 + x.dot(c);

    // determine which constraints are in the active set
    unsigned nactive = 0;
    for (unsigned i=0; i< n; i++)
      if (std::fabs(x[i] - l[i]) < std::numeric_limits<Real>::epsilon() || 
          std::fabs(x[i] - u[i]) < std::numeric_limits<Real>::epsilon())
      {
        inactive_set[i] = false;
        nactive++;
      }
      else
        inactive_set[i] = true;

    // form preconditioner inv(Wzz)
    H.select(inactive_set, Wzz);

    // form Z'GZ
    G.select_square(inactive_set, ZtGZ);

    FILE_LOG(LOG_OPT) << "  iteration: " << niter << endl;
    FILE_LOG(LOG_OPT) << "    preconditioner: " << Wzz << endl;
    FILE_LOG(LOG_OPT) << "    Z'GZ: " << endl << ZtGZ;
    FILE_LOG(LOG_OPT) << "    gradient: " << Gx_c << endl;
    FILE_LOG(LOG_OPT) << "    Cauchy point: " << x << endl;
    FILE_LOG(LOG_OPT) << "    q(x_c): " << q << endl;
    if (LOGGING(LOG_OPT))
    {
      std::ostringstream sstr;
      sstr << "   active set: ";
      for (unsigned i=0; i< n; i++)
        if (!inactive_set[i])
          sstr << i << " ";
      FILE_LOG(LOG_OPT) << sstr.str() << endl;
    }

    // prepare to do conjugate gradient iteration
    dxz.resize(n-nactive);
    xz.resize(n-nactive);
    lz.resize(n-nactive);
    uz.resize(n-nactive); 
    cz.resize(n-nactive);
    for (unsigned i=0, j=0; i< n; i++)
      if (inactive_set[i])
      {
        cz[j] = c[i];
        xz[j] = x[i];
        lz[j] = l[i];
        uz[j] = u[i];
        j++;
      }

    // *** Solve Z'GZ*xz = -cz (subspace minimization)
    if (cz.size() < 100)
    {
      dxz.copy_from(cz);
      dxz.negate();
      FILE_LOG(LOG_OPT) << " -- about to do positive definite subspace solve" << endl;
      FILE_LOG(LOG_OPT) << "ZtGZ: " << endl << ZtGZ;
      FILE_LOG(LOG_OPT) << "cz: " << dxz << endl;
      condition_and_factor_PD(ZtGZ);
      LinAlg::solve_chol_fast(ZtGZ, dxz);
      FILE_LOG(LOG_OPT) << "dxz: " << dxz << endl;

      // determine delta xz
      dxz -= xz;
      Real max_phi = 1.0;
      for (unsigned j=0; j< xz.size(); j++)
      {
        if (dxz[j] < (Real) 0.0 && lz[j] > -INF)
          max_phi = std::min(max_phi, (lz[j] - xz[j])/dxz[j]);
        else if (dxz[j] > (Real) 0.0 && uz[j] < INF)
          max_phi = std::min(max_phi, (uz[j] - xz[j])/dxz[j]);
        assert(max_phi >= -NEAR_ZERO);
      }
      
      if (max_phi < (Real) 0.0)
        max_phi = (Real) 0.0;
      FILE_LOG(LOG_OPT) << " -- maximum phi: " << max_phi << endl;

      // update xz
      dxz *= max_phi;
      xz += dxz;
    }     
    // *** do MINRES method (subspace minimization)
    else
    {
      y = -cz;
      r1 = -cz;
      MatrixN::diag_mult(Wzz, -cz, y);
      Real beta1 = -(cz.dot(y));

      // look for x = 0 solution
      if (beta1 == 0)
        xz.set_zero();
      else
      {
        // normalize y to get v1 later
        beta1 = std::sqrt(beta1);

        // init other quantities
        Real oldb = (Real) 0.0;
        Real beta = beta1;
        Real dbar = (Real) 0.0;
        Real epsln = (Real) 0.0;
        Real qrnorm = beta1;
        Real phibar = beta1;
        Real rhs1 = beta; 
        Real rhs2 = (Real) 0.0;
        Real tnorm2 = (Real) 0.0;
        Real ynorm2 = (Real) 0.0;
        Real cs = (Real) -1.0;
        Real sn = (Real) 0.0;
        w.set_zero(cz.size());
        w2.set_zero(cz.size());
        r2.copy_from(r1);

        // reset stop flag and setup gmax and gmin
        bool stop = false;
        Real gmax = std::numeric_limits<Real>::max();
        Real gmin = -gmax;

        // do the iterative method
//      for (unsigned i=0; i< n-nactive; i++)
        for (unsigned i=0;;i++)
        {
          Real s = (Real) 1.0/beta;
          v.copy_from(y) *= s;
          ZtGZ.mult(v, y);
          if (i > 0)
            y -= r1*(beta/oldb);

          Real alfa = v.dot(y);
          y += r2*(-alfa/beta);
          r1.copy_from(r2);
          r2.copy_from(y);
          MatrixN::diag_mult(Wzz, r2, y);
          oldb = beta;
          beta = r2.dot(y);
          if (beta < 0)
            break;
          beta = std::sqrt(beta);
          tnorm2 += alfa*alfa + oldb*oldb + beta*beta;

          // init a few things
          if (i == 0)
          {
            if (beta/beta1 <= 10*std::numeric_limits<Real>::epsilon())
              stop = true;
            gmax = std::fabs(alfa);
            gmin = gmax;
          }

          // apply previous rotation
          Real oldeps = epsln;
          Real delta = cs*dbar + sn*alfa;
          Real gbar = sn*dbar - cs*alfa;
          epsln = sn*beta;
          dbar = -cs*beta;
          Real root = std::sqrt(gbar*gbar + dbar*dbar);

          // compute the next plane rotation
          Real gamma = std::sqrt(gbar*gbar + beta*beta);
          gamma = std::max(gamma, std::numeric_limits<Real>::epsilon());
          cs = gbar/gamma;
          sn = beta/gamma;
          Real phi = cs*phibar;
          phibar = sn*phibar; 

          // prepare to update xz
          Real denom = (Real) 1.0/gamma;
          w1.copy_from(w2);
          w2.copy_from(w);
          w = (v - w1*oldeps - w2*delta)*denom;
 
          // determine the maximum phi 
          Real box_phi = std::numeric_limits<Real>::max();
          for (unsigned j=0; j< xz.size(); j++)
            if (w[j] < (Real) 0.0 && lz[j] > -INF)
              box_phi = std::min(box_phi, (lz[j] - xz[j])/w[j]);
            else if (w[j] > (Real) 0.0 && uz[j] < INF)
              box_phi = std::min(box_phi, (uz[j] - xz[j])/w[j]);

          // update xz
          xz += w*std::min(box_phi, phi);

          // look for stopping condition
          if (box_phi < phi)
            break;

          // prepare to loop again
          gmax = std::max(gmax, gamma);
          gmin = std::min(gmin, gamma);
          Real z = rhs1/gamma;
          ynorm2 = z*z + ynorm2;
          rhs1 = rhs2 - delta*z;
     
          // estimate various norms to see whether we can quit
          Real Anorm = std::sqrt(tnorm2);
          Real ynorm = std::sqrt(ynorm2);
          Real epsx = Anorm*ynorm*std::numeric_limits<Real>::epsilon();
          qrnorm = phibar;
          Real rnorm = qrnorm;
          Real test1 = rnorm/(Anorm*ynorm);
          Real test2 = root/Anorm;

          // see whether we can stop
          if (stop)
            break;

          // check to see whether we can quit
          Real t1 = (Real) 1.0 + test1;
          Real t2 = (Real) 1.0 + test2;
          if (t1 <= (Real) 1.0 || t2 <= (Real) 1.0 || epsx >= beta1 || 
              test2 <= NEAR_ZERO || test1 <= NEAR_ZERO)
          {
            FILE_LOG(LOG_OPT) << " -- MINRES terminating; " << i << " iterations used" << endl;
            break;
          }
        }
      }
    }

    // project xz back to x
    xstar.copy_from(x);
    for (unsigned i=0, j=0; i< n; i++)
      if (inactive_set[i])
        xstar[i] = xz[j++];

    // manually enforce constraints
    for (unsigned i=0; i< xstar.size(); i++)
      if (xstar[i] < l[i])
        xstar[i] = l[i];
      else if (xstar[i] > u[i])
        xstar[i] = u[i];

    // evaluate new x; if it is better than old x, keep it
    Real qstar = xstar.dot(G.mult(xstar, Gx)) * (Real) 0.5 + xstar.dot(c);
    if (qstar < q)
    {
      x.copy_from(xstar);
      q = qstar;
    }

    // if objective function did not decrease sufficiently, quit
    if (qmin - q < tol)
    {
      FILE_LOG(LOG_OPT) << " -- objective function did not decrease sufficiently; quitting!" << endl;
      break;
    }
    else // if (q <= qmin)
    {
      if (std::isnan(q))
        return niter;
      else
        assert(q <= qmin);
      FILE_LOG(LOG_OPT) << " -- objective function decreased by " << (qmin - q) << std::endl;
      qmin = q;
    }

    FILE_LOG(LOG_OPT) << "    x*: " << xstar << endl;
    FILE_LOG(LOG_OPT) << "    q(x*): " << qstar << endl;
  }

  return niter;
}

// auxiliary function for LCP solver: prints out variable type (z or w)
std::string var(unsigned v, unsigned n)
{
  if (v == n*2)
    return std::string("z0"); 
  std::ostringstream out;
  if (v < n)
    out << "z";
  else
  {
    out << "w";
    v -= n;
  }
  out << v;
  return out.str();
}

#ifdef _OPENMP
/// Converts an unsigned long long to a binary vector
void to_binary(unsigned long long i, vector<bool>& bin)
{
  // compute highest representable integer using boolean array
  unsigned n = bin.size();
  unsigned long long n2 = 1;
  for (unsigned j=n-1; j > 0; j--)
    n2 *= 2;

  // setup the binary number
  for (unsigned j=0; j< n; j++)
  {
    bin[j] = (i / n2 > 0);
    unsigned long long rem = i % n2;
    i = rem;
    n2 /= 2;
  }
} 

/// Enumerative algorithm for solving linear complementarity problems
/**
 */
void Optimization::lcp_enum(const MatrixN& M, const VectorN& q, vector<VectorN>& z)
{  
  const Real MIN_TOLERANCE = -NEAR_ZERO;
  const unsigned n = q.size();

  // get the number of threads
  const unsigned NTHREADS = omp_get_max_threads();

  // make NTHREADS copies of everything we need
  vector<VectorN> qq(NTHREADS);
  vector<MatrixN> Mij(NTHREADS);
  vector<MatrixN> Mnij(NTHREADS);
  vector<VectorN> qx(NTHREADS);
  vector<VectorN> qni(NTHREADS);
  vector<VectorN> result(NTHREADS);
  vector<vector<unsigned> > basicw(NTHREADS);
  vector<vector<unsigned> > nbasicz(NTHREADS);
  vector<vector<unsigned> > nbasicw(NTHREADS);

  // setup an array of bool's; false indicates variable wi non-basic/zi basic
  vector<vector<bool> > binary(NTHREADS);
  for (unsigned i=0; i< binary.size(); i++)
    binary[i] = vector<bool>(n);

  // start iterating through all 2^n possibilities
  const unsigned long long TWON = intpow(2,n);
  #pragma omp parallel for 
  for (unsigned long long i=0; i< TWON; i++)
  {
    // get the thread ID
    unsigned tID = omp_get_thread_num();

    // determine binary representation of i
    to_binary(i, binary[tID]);

    // setup indices for basic w variables
    unsigned num_basicw = 0;
    basicw[tID].clear();
    nbasicw[tID].clear();
    nbasicz[tID].clear();
    for (unsigned k=0; k< n; k++)
      if (!binary[tID][k])
        nbasicw[tID].push_back(k);
      else
      {
        basicw[tID].push_back(k);
        nbasicz[tID].push_back(k);
        num_basicw++;
      }

    // resize matrices
    Mij[tID].resize(num_basicw, num_basicw);
    Mnij[tID].resize(n-num_basicw, num_basicw);
    qx[tID].resize(num_basicw);

    // get qni
    qni[tID].resize(n-num_basicw);
    q.select(nbasicw[tID].begin(), nbasicw[tID].end(), qni[tID]);

    // select appropriate components of M
    M.select(basicw[tID].begin(), basicw[tID].end(), nbasicz[tID].begin(), nbasicz[tID].end(), Mij[tID]);        
    M.select(nbasicw[tID].begin(), nbasicw[tID].end(), nbasicz[tID].begin(), nbasicz[tID].end(), Mnij[tID]);

    // select components of q corresponding to basic indices of w
    q.select(basicw[tID].begin(), basicw[tID].end(), qx[tID]);
    qx[tID] *= -1.0;

    // determine inv(Mij)*q;
    bool singular = false;
    try
    {
      LinAlg::solve_fast(Mij[tID], qx[tID]);
    }
    catch (SingularException e)
    {
      singular = true;
    }

    // only proceed if not singular
    if (!singular)
    {
      // compute Mnij * qx + qni
      Mnij[tID].mult(qx[tID], result[tID]);
      result[tID] += qni[tID];        

      // form qq
      qq[tID].resize(n);
      qq[tID].set_sub_vec(0, qx[tID]);
      qq[tID].set_sub_vec(qx[tID].size(), result[tID]);

      // see whether there are any negative components of qq; if there are,
      // this is not a solution
      Real minqq = *std::min_element(qq[tID].begin(), qq[tID].end());
      if (minqq > MIN_TOLERANCE)
      {
        // this is a viable solution; get z
        #pragma omp critical
        {
          z.push_back(VectorN::zero(n));
          unsigned idx = 0;
          for (vector<unsigned>::const_iterator k = nbasicz[tID].begin(); k != nbasicz[tID].end(); k++)
            z.back()[*k] = qq[tID][idx++];
        }
      }
    }
  }
}
#else
/// Enumerative algorithm for solving linear complementarity problems
/**
 */
void Optimization::lcp_enum(const MatrixN& M, const VectorN& q, vector<VectorN>& z)
{
  const Real MIN_TOLERANCE = -NEAR_ZERO;
  vector<unsigned> basicw, nbasicz, nbasicw;
  const unsigned n = q.size();
  VectorN qq(n);
  MatrixN Mij;
  MatrixN Mnij;
  VectorN qx, qni, result;

  // setup an array of bool's; false indicates variable wi non-basic/zi basic
  vector<bool> counter(n, false);

  // start iterating through all 2^n possibilities
  const unsigned long long TWON = intpow(2,n);
  for (unsigned long i=0; i< TWON; i++)
  {
    // setup indices for basic w variables
    unsigned num_basicw = 0;
    basicw.clear();
    nbasicw.clear();
    nbasicz.clear();
    for (unsigned k=0; k< n; k++)
      if (!counter[k])
        nbasicw.push_back(k);
      else
      {
        basicw.push_back(k);
        nbasicz.push_back(k);
        num_basicw++;
      }

    // resize matrices
    Mij.resize(num_basicw, num_basicw);
    Mnij.resize(n-num_basicw, num_basicw);
    qx.resize(num_basicw);

    // get qni
    qni.resize(n-num_basicw);
    q.select(nbasicw.begin(), nbasicw.end(), qni);

    // select appropriate components of M
    M.select(basicw.begin(), basicw.end(), nbasicz.begin(), nbasicz.end(), Mij);        
    M.select(nbasicw.begin(), nbasicw.end(), nbasicz.begin(), nbasicz.end(), Mnij);

    // select components of q corresponding to basic indices of w
    q.select(basicw.begin(), basicw.end(), qx);
    qx *= -1.0;

    // determine inv(Mij)*q;
    bool singular = false;
    try
    {
      LinAlg::solve_fast(Mij, qx);
    }
    catch (SingularException e)
    {
      singular = true;
    }

    // only proceed if not singular
    if (!singular)
    {
      // compute Mnij * qx + qni
      Mnij.mult(qx, result);
      result += qni;        

      // form qq
      qq.set_sub_vec(0, qx);
      qq.set_sub_vec(qx.size(), result);

      // see whether there are any negative components of qq; if there are,
      // this is not a solution
      Real minqq = *std::min_element(qq.begin(), qq.end());
      if (minqq > MIN_TOLERANCE)
      {
        // this is a viable solution; get z
        z.push_back(VectorN::zero(n));
        unsigned idx = 0;
        for (vector<unsigned>::const_iterator k = nbasicz.begin(); k != nbasicz.end(); k++)
          z.back()[*k] = qq[idx++];
      }
    }

    // increment the counter
    for (int k=n-1; k >= 0; k--)
      if (counter[k])
        counter[k] = false;
      else
      {
        counter[k] = true;
        break;
      }
  }
}
#endif

/// Conjugate gradient-type algorithm for solving convex LCP problems
bool Optimization::lcp_gradproj(const MatrixN& M, const VectorN& q, VectorN& x, Real tol, unsigned max_iter)
{
  const unsigned n = q.size();
  const unsigned CHECK_FREQ = 100;

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // setup some variables we'll need
  VectorN c, d, w;

  // setup c and d
  c.set_zero(n);
  d.resize(n);
  for (unsigned i=0; i< n; i++)
    d[i] = std::numeric_limits<Real>::max();

  // setup x; try to use current value
  x.resize(n);
  for (unsigned i=0; i< n; i++)
    if (x[i] < (Real) 0.0)
      x[i] = (Real) 0.0;

  // run projected gradient algorithm
  for (unsigned iter=1; iter<= max_iter; iter+= CHECK_FREQ)
  {
    FILE_LOG(LOG_OPT) << "Optimization::lcp_projgrad() iteration: " << iter << std::endl;
    if (qp_gradproj(M, q, c, d, CHECK_FREQ, x, std::numeric_limits<Real>::epsilon()));
      return true;

    // determine w = Mx + q
    M.mult(x, w) += q;

    // verify w > 0 (to tolerance)
    if (*std::min_element(w.begin(), w.end()) < -tol)
      continue;

    // verify x'w = 0 (to tolerance)
    if (std::fabs(x.dot(w)) > tol)
      continue;
    else
      return true;
  }

  return false;
}

/// Converts a QP into the LCP (qq, MM) 
/**
 * Minimizes the QP 0.5*x'*G*x + x'*c subject to x >= 0 and M*x >= q
 */
void Optimization::qp_to_lcp1(const MatrixN& G, const VectorN& c, const MatrixN& M, const VectorN& q, MatrixN& MM, VectorN& qq)
{
  const unsigned N = G.rows();

  if (G.columns() != N)
    throw NonsquareMatrixException();

  // if n is zero, do nothing
  if (N == 0)
  {
    MM.resize(0,0);
    qq.resize(0);
    return;
  }

  // form the LCP matrix
  MM.set_zero(N + M.rows(), N + M.rows());
  MM.set_sub_mat(0, 0, G);
  MM.set_sub_mat(0, N, M, true);  // set submatrix of MM to M'
  MM.set_sub_mat(N, 0, M);

  // negate upper RHS of MM
  for (unsigned i=0; i< N; i++)
    for (unsigned j=N; j< MM.columns(); j++)
      MM(i,j) = -MM(i,j);

  // form the LCP vector
  qq.resize(MM.rows());
  qq.set_sub_vec(0, c);
  qq.set_sub_vec(N, -q);
}

/// Converts a QP into the LCP (qq, MM)
/**
 * Minimizes the QP 0.5*x'*G*x + x'*c subject to M*x >= q
 */
void Optimization::qp_to_lcp2(const MatrixN& G, const VectorN& c, const MatrixN& M, const VectorN& q, MatrixN& MM, VectorN& qq)
{
  SAFESTATIC MatrixN workM;

  // determine # of LCP variables
  const unsigned N = G.rows();
  const unsigned NLCP = (N + M.rows())*2;

  if (G.columns() != N)
    throw NonsquareMatrixException();

  // setup the LCP
  MM.set_zero(NLCP,NLCP);
  qq.set_zero(NLCP);

  // setup c components of qq
  qq.set_sub_vec(0, c);
  for (unsigned i=N, j=0; i< N+N; i++, j++)
    qq[i] = -c[j];

  // setup q components of qq
  for (unsigned i=N+N, j=0; i< qq.size(); i++, j++)
  {
    qq[i] = -q[j];
    qq[i+M.rows()] = -q[j];
  }

  // setup MM
  workM.copy_from(G).negate();
  MM.set_sub_mat(0, 0, G);
  MM.set_sub_mat(N, N, G);
  MM.set_sub_mat(0, N, workM);
  MM.set_sub_mat(N, 0, workM);
  workM.copy_from(M).negate();
  MM.set_sub_mat(0, N*2, workM, true);
  MM.set_sub_mat(0, N*2+M.rows(), M, true);
  MM.set_sub_mat(N*2, 0, M);
  MM.set_sub_mat(N*2+M.rows(), 0, workM);
}

/// Regularized wrapper around Lemke's algorithm
bool Optimization::lcp_lemke_regularized(const MatrixN& M, const VectorN& q, VectorN& z, int min_exp, unsigned step_exp, int max_exp, Real piv_tol, Real zero_tol)
{
  FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke_regularized() entered" << endl;
  SAFESTATIC FastThreadable<VectorN> w_x, qq_x, qe_x;
  SAFESTATIC FastThreadable<MatrixN> MM_x, Me_x;
  VectorN& w = w_x();
  VectorN& qq = qq_x();
  MatrixN& MM = MM_x();
  MatrixN& Me = Me_x();
  VectorN& qe = qe_x();

  // look for fast exit
  if (q.size() == 0)
  {
    z.resize(0);
    return true;
  }

  // copy M and q
  Me.copy_from(M);
  qe.copy_from(q);

  // equilibrate M and q
  equilibrate(Me, qe);

  // copy MM and qq  
  MM.copy_from(Me);
  qq.copy_from(qe);

  // assign value for zero tolerance, if necessary
  const Real ZERO_TOL = (zero_tol > (Real) 0.0) ? zero_tol : q.size() * std::numeric_limits<Real>::epsilon();

  // try non-regularized version first
  bool result = lcp_lemke(MM, qq, z, piv_tol, zero_tol);
  if (result)
  {
    // verify that solution truly is a solution -- check z
    if (*std::min_element(z.begin(), z.end()) >= -ZERO_TOL)
    {
      // check w
      Me.mult(z, w) += qe;
      if (*std::min_element(w.begin(), w.end()) >= -ZERO_TOL)
      {
        // check z'w
        std::transform(z.begin(), z.end(), w.begin(), w.begin(), std::multiplies<Real>());
        pair<Real*, Real*> mmax = boost::minmax_element(w.begin(), w.end());
        if (*mmax.first >= -ZERO_TOL && *mmax.second < ZERO_TOL)
        {
          FILE_LOG(LOG_OPT) << "  solved with no regularization necessary!" << endl;
          FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke_regularized() exited" << endl;
          return true;
        }
      }
    }
  }

  // start the regularization process
  int rf = min_exp;
  while (rf < max_exp)
  {
    // setup regularization factor
    Real lambda = std::pow((Real) 10.0, (Real) rf);

    // regularize M
    MM.copy_from(Me);
    for (unsigned i=0; i< MM.rows(); i++)
      MM(i,i) += lambda;

    // recopy q
    qq.copy_from(qe);

    // try to solve the LCP
    if ((result = lcp_lemke(MM, qq, z, piv_tol, zero_tol)))
    {
      // verify that solution truly is a solution -- check z
      if (*std::min_element(z.begin(), z.end()) > -ZERO_TOL)
      {
        // check w
        MM.copy_from(Me);
        for (unsigned i=0; i< MM.rows(); i++)
          MM(i,i) += lambda;
        qq.copy_from(qe);
        MM.mult(z, w) += qq;
        if (*std::min_element(w.begin(), w.end()) > -ZERO_TOL)
        {
          // check z'w
          std::transform(z.begin(), z.end(), w.begin(), w.begin(), std::multiplies<Real>());
          pair<Real*, Real*> mmax = boost::minmax_element(w.begin(), w.end());
          if (*mmax.first > -ZERO_TOL && *mmax.second < ZERO_TOL)
          {
            FILE_LOG(LOG_OPT) << "  solved with regularization factor: " << lambda << endl;
            FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke_regularized() exited" << endl;

            return true;
          }
        }
      }
    }

    // increase rf
    rf += step_exp;
  }

  FILE_LOG(LOG_OPT) << "  unable to solve given any regularization!" << endl;
  FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke_regularized() exited" << endl;

  // still here?  failure...
  return false;
}

/// Lemke's algorithm for solving linear complementarity problems
/**
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
bool Optimization::lcp_lemke(const MatrixN& M, const VectorN& q, VectorN& z, Real piv_tol, Real zero_tol)
{
  const unsigned n = q.size();
  const unsigned MAXITER = std::min((unsigned) 1000, 50*n);

  // look for immediate exit
  if (n == 0)
  {
    z.resize(0);
    return true;
  }

  // setup work variables
  SAFESTATIC FastThreadable<VectorN> Be_x, U_x, z0_x, x_x, d_x, xj_x, dj_x, w_x, result_x;
  SAFESTATIC FastThreadable<MatrixN> B_x, A_x, t1_x, t2_x;
  SAFESTATIC FastThreadable<vector<unsigned> > all_x, tlist_x, bas_x, nonbas_x, j_x; 

  // get references to all variables
  VectorN& Be = Be_x();
  VectorN& U = U_x();
  VectorN& z0 = z0_x();
  VectorN& x = x_x();
  VectorN& d = d_x();
  VectorN& xj = xj_x();
  VectorN& dj = dj_x();
  VectorN& w = w_x();
  VectorN& result = result_x();
  MatrixN& B = B_x();
  MatrixN& A = A_x();
  MatrixN& t1 = t1_x();
  MatrixN& t2 = t2_x();
  vector<unsigned>& all = all_x();
  vector<unsigned>& tlist = tlist_x();
  vector<unsigned>& bas = bas_x();
  vector<unsigned>& nonbas = nonbas_x();
  vector<unsigned>& j = j_x();

  // clear all vectors
  all.clear();
  tlist.clear();
  bas.clear();
  nonbas.clear();
  j.clear();

  // copy z to z0
  z0.copy_from(z);

  // come up with a sensible value for zero tolerance if none is given
  if (zero_tol <= (Real) 0.0)
    zero_tol = std::numeric_limits<Real>::epsilon() * M.norm_inf() * n;

  FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() entered" << endl;
  FILE_LOG(LOG_OPT) << "  M: " << endl << M;
  FILE_LOG(LOG_OPT) << "  q: " << q << endl;

  // see whether trivial solution exists
  if (*std::min_element(q.begin(), q.end()) > -zero_tol)
  {
    FILE_LOG(LOG_OPT) << " -- trivial solution found" << endl;
    FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << endl;
    z.set_zero(n);
    return true;
  }

  // initialize variables
  z.set_zero(n*2);
  unsigned t = 2*n;
  unsigned entering = t;
  unsigned leaving = 0;
  all.clear();
  for (unsigned i=0; i< n; i++)
    all.push_back(i);
  unsigned lvindex;
  unsigned idx;
  vector<unsigned>::iterator iiter;
  tlist.clear();

  // determine initial basis
  bas.clear();
  nonbas.clear();
  if (z0.size() != n)
    for (unsigned i=0; i< n; i++)
      nonbas.push_back(i);
  else
    for (unsigned i=0; i< n; i++)
      if (z0[i] > 0)
        bas.push_back(i);
      else
        nonbas.push_back(i);

  // B should ideally be a sparse matrix
  B.set_identity(n);
  B.negate();

  // determine initial values
  if (!bas.empty())
  {
    // select columns of M corresponding to z vars in the basis
    M.select(all.begin(), all.end(), bas.begin(), bas.end(), t1);

    // select columns of I corresponding to z vars not in the basis
    B.select(all.begin(), all.end(), nonbas.begin(), nonbas.end(), t2);

    // setup the basis matrix
    B.resize(n, t1.columns() + t2.columns());
    B.set_sub_mat(0,0,t1);
    B.set_sub_mat(0,t1.columns(),t2);
  }

  // solve B*x = -q
  try
  {
    A.copy_from(B);
    x.copy_from(q);
    LinAlg::solve_fast(A, x);
  }
  catch (SingularException e)
  {
    try
    {
      // use slower SVD pseudo-inverse
      A.copy_from(B);
      x.copy_from(q);
      LinAlg::solve_LS_fast1(A, x);
    }
    catch (NumericalException e)
    {
      A.copy_from(B);
      LinAlg::solve_LS_fast2(A, x);
    }
  }
  x.negate();

  // check whether initial basis provides a solution
  if (std::find_if(x.begin(), x.end(), std::bind2nd(std::less<Real>(), 0.0)) == x.end())
  {
    for (idx = 0, iiter = bas.begin(); iiter != bas.end(); iiter++, idx++)
      z[*iiter] = x[idx];
    z.resize(n, true);

    // check to see whether tolerances are satisfied
    FILE_LOG(LOG_OPT) << " -- initial basis provides a solution!" << std::endl;
    if (LOGGING(LOG_OPT))
    {
      M.mult(z, w) += q;
      Real minw = *std::min_element(w.begin(), w.end());
      Real w_dot_z = std::fabs(w.dot(z));
      FILE_LOG(LOG_OPT) << "  z: " << z << std::endl;
      FILE_LOG(LOG_OPT) << "  w: " << w << std::endl;
      FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
      FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
    }
    FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << endl;

    return true; 
  }

  // determine initial leaving variable
  Real* minx = std::min_element(x.begin(), x.begin() + n);
  Real tval = -*minx;
  BOOST_FOREACH(unsigned i, nonbas) // add w variables to basis
    bas.push_back(i+n);
  lvindex = std::distance(x.begin(), minx);
  iiter = bas.begin();
  std::advance(iiter, lvindex);
  leaving = *iiter;
  FILE_LOG(LOG_OPT) << " pivoting " << var(leaving,n) << " and " << var(entering, n) << endl;

  // pivot in the artificial variable
  *iiter = t;    // replace w var with z0 in basic indices
  U.resize(n);
  for (unsigned i=0; i< n; i++)
    U[i] = (x[i] < 0.0) ? 1.0 : 0.0;
  B.mult(U, Be);
  Be.negate();
  x += U*tval;
  x[lvindex] = tval;
  B.set_column(lvindex, Be);
  FILE_LOG(LOG_OPT) << "  new q: " << x << endl;

  // main iterations begin here
  for (unsigned iter=0; iter < MAXITER; iter++)
  {
    // check whether done; if not, get new entering variable
    if (leaving == t)
    {
      FILE_LOG(LOG_OPT) << "pivoting " << var(leaving,n) << " and " << var(entering,n) << endl;
      FILE_LOG(LOG_OPT) << "-- solved LCP successfully!" << endl;
      unsigned idx;
      for (idx = 0, iiter = bas.begin(); iiter != bas.end(); iiter++, idx++)
        z[*iiter] = x[idx];
      z.resize(n, true);

      // verify tolerances
      if (LOGGING(LOG_OPT))
      {
        M.mult(z, w) += q;
        Real minw = *std::min_element(w.begin(), w.end());
        Real w_dot_z = std::fabs(w.dot(z));
        FILE_LOG(LOG_OPT) << "  found solution!" << std::endl;
        FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
        FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
      }
      FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << endl;

      return true; 
    }
    else if (leaving < n)
    {
      entering = n + leaving;
      Be.set_zero(n);
      Be[leaving] = -1;
    }
    else
    {
      entering = leaving - n;
      M.get_column(entering, Be);
    }
    d.copy_from(Be);
    try
    {
      A.copy_from(B);
      LinAlg::solve_fast(A, d);
    }
    catch (SingularException e)
    {
      try
      {
        // use slower SVD pseudo-inverse
        A.copy_from(B);
        d.copy_from(Be);
        LinAlg::solve_LS_fast1(A, d);
      }
      catch (NumericalException e)
      {
        A.copy_from(B);
        LinAlg::solve_LS_fast2(A, d);
      }
    }

    // use a new pivot tolerance if necessary
    const Real PIV_TOL = (piv_tol > (Real) 0.0) ? piv_tol : std::numeric_limits<Real>::epsilon() * n * std::max((Real) 1.0, Be.norm_inf());

    // ** find new leaving variable
    j.clear();
    for (unsigned i=0; i< d.size(); i++)
      if (d[i] > PIV_TOL)
        j.push_back(i);
    // check for no new pivots; ray termination
    if (j.empty())
    {
      FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() - no new pivots (ray termination)" << endl;
      FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << endl;

      z.resize(n, true);
      return false;
    }

    FILE_LOG(LOG_OPT) << " -- column of M': " << d << endl;

    // select elements j from x and d
    xj.resize(j.size());
    dj.resize(xj.size());
    select(x.begin(), j.begin(), j.end(), xj.begin());
    select(d.begin(), j.begin(), j.end(), dj.begin());

    // compute minimal ratios x(j) + NEAR_ZERO ./ d(j), d > 0
    result.resize(xj.size());
    std::transform(xj.begin(), xj.end(), result.begin(), std::bind2nd(std::plus<Real>(), zero_tol));
    std::transform(result.begin(), result.end(), dj.begin(), result.begin(), std::divides<Real>());
    Real theta = *std::min_element(result.begin(), result.end());

    // NOTE: lexicographic ordering does not appear to be used here to prevent
    // cycling (see [Cottle 1992], pp. 340-342)
    // find indices of minimal ratios, d> 0
    //   divide x(j) ./ d(j) -- remove elements above the minimum ratio
    std::transform(xj.begin(), xj.end(), dj.begin(), result.begin(), std::divides<Real>());
    for (iiter = j.begin(), idx = 0; iiter != j.end(); )
      if (result[idx++] <= theta)
        iiter++;
      else
        iiter = j.erase(iiter);

    // if j is empty, then likely the zero tolerance is too low
    if (j.empty())
    {
      FILE_LOG(LOG_OPT) << "zero tolerance too low?" << std::endl;
      FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << std::endl;
      z.resize(n, true);
      return false;
    }

    // check whether artificial index among these
    tlist.clear();
    select(bas.begin(), j.begin(), j.end(), std::back_inserter(tlist));
    iiter = std::find(tlist.begin(), tlist.end(), t);
    if (iiter != tlist.end()) 
      lvindex = std::distance(tlist.begin(), iiter);
    else
    {
      // redetermine dj
      dj.resize(j.size());
      select(d.begin(), j.begin(), j.end(), dj.begin());

      // get the maximum
      Real* maxdj = std::max_element(dj.begin(), dj.end());
      lvindex = std::distance(dj.begin(), maxdj);
    }
    assert(lvindex < j.size());
    select(j.begin(), &lvindex, &lvindex+1, &lvindex);

    // set leaving = bas(lvindex)
    iiter = bas.begin();
    std::advance(iiter, lvindex);
    leaving = *iiter;

    // ** perform pivot
    Real ratio = x[lvindex]/d[lvindex];
    d*= ratio;
    x -= d;
    x[lvindex] = ratio;
    B.set_column(lvindex, Be);
    *iiter = entering;

    FILE_LOG(LOG_OPT) << "pivoting " << var(leaving,n) << " and " << var(entering,n) << endl;
  }

  FILE_LOG(LOG_OPT) << " -- maximum number of iterations exceeded" << endl;
  FILE_LOG(LOG_OPT) << "Optimization::lcp_lemke() exited" << std::endl;

  // max iterations exceeded
  z.resize(n, true);
  
  return false;
}

/// Interior point method for solving convex linear complementarity problems
bool Optimization::lcp_convex_ip(const MatrixN& M, const VectorN& q, VectorN& z, Real tol, Real eps, Real eps_feas, unsigned max_iterations)
{
  const unsigned n = q.size();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // setup convex optimization parameters
  OptParams cparams;
  cparams.n = q.size();
  cparams.m = cparams.r = 0;
  cparams.max_iterations = max_iterations;
  cparams.eps = eps;
  cparams.eps_feas = eps_feas;
  cparams.lb.set_zero(n);
  cparams.q.copy_from(q).negate();
  cparams.M.copy_from(M);

  // do the convex QP IP solver
  if (!qp_convex_ip(M, q, cparams, z))
  {
    FILE_LOG(LOG_OPT) << "Optimization::lcp_convex_ip() - convex optimization failed!" << endl;
    return false; 
  }

  if (LOGGING(LOG_OPT))
  {
    // check z >= 0, w>= 0, w'z = 0
    VectorN w = M*z + q;
    Real minz = *std::min_element(z.begin(), z.end());
    Real minw = *std::min_element(w.begin(), w.end());
    Real zTw = z.dot(w);
    FILE_LOG(LOG_OPT) << "Optimization::lcp_convex_ip() successful!" << endl;
    FILE_LOG(LOG_OPT) << "  -- minimum z: " << minz << endl;
    FILE_LOG(LOG_OPT) << "  -- minimum w: " << minw << endl;
    FILE_LOG(LOG_OPT) << "  -- z'w: " << zTw << endl;
  }

  return true;
}

/// Iterative method for solving linear complementarity problems w/symmetric M
/**
 * Solves problems of the form Ax + b = w, where x'w = 0, x >= 0 and w >= 0
 * \note this method is via [Murty 1988]
 * \note convergence is not guaranteed
 */
bool Optimization::lcp_iter_symm(const MatrixN& M, const VectorN& q, VectorN& z, Real tol, const unsigned iter)
{
  const unsigned n = q.size();

  if (M.rows() != n)
    throw MissizeException();

  if (M.columns() != n)
    throw NonsquareMatrixException();

  FILE_LOG(LOG_OPT) << "Optimization::solve_lcp_symm_iter() entered" << endl;
  FILE_LOG(LOG_OPT) << "M: " << std::endl << M;
  FILE_LOG(LOG_OPT) << "q: " << q << std::endl;

  // resize z, if necessary
  if (z.size() != n)
    z.set_zero(n);

  // NOTE: uses the projected SOR scheme (Murty 1988, p. 373)
  // NOTE: we use E as identity
  // NOTE: we use L as K 
 
  // setup work vectors and matrices
  VectorN row, znew;
  vector<VectorN> m;
  MatrixN L, G, U;
 
  // get rows of M
  m.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    M.get_row(i, row);
    m[i] = row;
  } 

  // compute L and G
  L.set_zero(n,n);
  G.set_zero(n,n);
  U.set_zero(n,n);
  for (unsigned i=0; i< n; i++)
    for (unsigned j=0; j< n; j++)
    {
      if (i > j)
        L(i,j) = M(i,j);
//      else if (i < j)
//        U(i,j) = M(i,j);
//      else
        else if (i == j)
          G(i,j) = M(i,j);
    }

  // compute info to determine lambda/omega
  Real max_Gjj = -1;
  for (unsigned i=0; i< n; i++)
    max_Gjj = std::max(max_Gjj, G(i,i));
  if (max_Gjj <= 0.0)
  {
    FILE_LOG(LOG_OPT) << " -- maximum diagonal element of G is not positive!" << endl;
    return false;
  }

  // determine lambda and omega; we'll set lambda arbitrarily
  const Real lambda = (Real) 0.5;
  const Real omega = (Real) 1.0/(lambda*max_Gjj);
  assert(lambda*omega < 2/max_Gjj);

  // setup new z
  znew.resize(n);

  // iterate..
  for (unsigned i=0; i< iter; i++)
  {
    // determine the new z (p. 367)
    znew[0] = std::max(lambda*(z[0] - omega*(m[0].dot(z)+q[0])),(Real) 0.0) + (1-lambda)*z[0];
    for (unsigned j=1; j< n; j++)
    {
      Real subsum = 0;
      for (unsigned l=0; l< j; l++)
        subsum += L(j,l)*(znew[l] - z[l]);
      znew[j] = std::max(lambda*(z[j] - omega*(m[j].dot(z)+q[j]+subsum)),(Real) 0.0) + (1-lambda)*z[j];
    }

    // swap the two
    std::swap(z, znew);

    // if there is little difference between the two, quit
    if ((znew - z).norm() < tol)
      break;

    FILE_LOG(LOG_OPT) << "iteration: " << i << endl;
    FILE_LOG(LOG_OPT) << "  new z: " << z << endl;
    FILE_LOG(LOG_OPT) << "  new w: " << (M*z + q) << endl;
    FILE_LOG(LOG_OPT) << "  new z'w: " << z.dot(M*z+q) << endl;
  }

  // check the solution
  if (LOGGING(LOG_OPT))
  {
    VectorN w = M * z + q;
    Real minz = *std::min_element(z.begin(), z.end());
    Real minw = *std::min_element(w.begin(), w.end());
    Real zw = z.dot(w);
    FILE_LOG(LOG_OPT) << "minimum z: " << minz << endl;
    FILE_LOG(LOG_OPT) << "minimum w: " << minw << endl;
    FILE_LOG(LOG_OPT) << "z'w: " << zw << endl;
  }

  FILE_LOG(LOG_OPT) << "Optimization::solve_lcp_symm_iter() exiting" << endl; 

  return true;
}

/// Iterative method for solving linear complementarity problems w/PD M
/**
 * Solves problems of the form Ax + b = w, where x'w = 0, x >= 0 and w >= 0
 * \note this method is via [Murty 1988]
 * \note convergence is guaranteed
 */
bool Optimization::lcp_iter_PD(const MatrixN& M, const VectorN& q, VectorN& x, Real tol, const unsigned iter)
{
  const unsigned n = q.size();

  FILE_LOG(LOG_OPT) << "Optimization::lcp_iter_PD() entered" << endl;

  MatrixN eye, eyepM, B;
  VectorN r, oldx;

  // resize x if necessary
  x.resize(n);

  // setup I + M and I - M (we call the latter 'B')
  eye.set_identity(q.size());
  eyepM.copy_from(M);
  B.copy_from(M);
  B.negate();
  for (unsigned i=0; i< n; i++)
  {
    eyepM(i,i) += (Real) 1.0;
    B(i,i) += (Real) 1.0;
  }

  // attempt to do Cholesky factorization of M + I 
  if (!LinAlg::factor_chol(eyepM))
  {
    FILE_LOG(LOG_OPT) << " -- M is not positive definite!" << endl;
    return false;
  }

  // solve (I + M) \ (I - M) and (I + M) \ -q
  LinAlg::solve_chol_fast(eyepM, B);
  r.copy_from(q);
  r.negate();
  LinAlg::solve_chol_fast(eyepM, r);

  // iterate the specified number of times
  for (unsigned k=0; k< iter; k++)
  {
    // store last x
    oldx.copy_from(x);

    // take the absolute value of x
    for (unsigned i=0; i< q.size(); i++)
      x[i] = std::fabs(x[i]);

    // update x
    B.mult(x, x) += r;

    if (LOGGING(LOG_OPT))
    {
      VectorN z = (const VectorN&) x;
      for (unsigned j=0; j< n; j++) z[j] += std::fabs(z[j]);
      VectorN w = M*z + q;
      Real minx = *std::min_element(z.begin(), z.end());
      Real minw = *std::min_element(w.begin(), w.end());
      Real xw = x.dot(w);
      FILE_LOG(LOG_OPT) << " iter " << k << ":  min x: " << minx << "  min w: " << minw << "  x'w: " << xw << endl;
    }

    // see whether we can stop
    if ((oldx - x).norm() < tol)
    {
      FILE_LOG(LOG_OPT) << " -- found a solution in " << k << " iterations!" << endl;
      break;
    }
  }

  // convert x to z
  for (unsigned i=0; i< q.size(); i++)
    x[i] += std::fabs(x[i]);

  if (LOGGING(LOG_OPT))
  {
    // check that x is >= 0
    Real minx = *std::min_element(x.begin(), x.end());
    FILE_LOG(LOG_OPT) << " -- minimum value of x: " << minx << endl;

    // check that r is >= 0
    M.mult(x, r) += q;
    Real minr = *std::min_element(r.begin(), r.end());
    FILE_LOG(LOG_OPT) << " -- minimum value of w: " << minr << endl;

    // check complementarity condition
    Real rdotx = r.dot(x);
    FILE_LOG(LOG_OPT) << " -- value of x'w: " << rdotx << endl;
  }

  return true;
}

/// Solves a linear program using the simplex method
/**
 *  Solves the linear program:  min c'x
 *                 subject to:  Ax = b
 *                              Mx >= q
 * \param the optimal point (if any) on return
 * \return <b>false</b> if problem infeasible; <b>true</b> otherwise
 */
bool Optimization::lp_simplex(const LPParams& lpparams, VectorN& x)
{
  #ifdef USE_GLPK
  // create the GLPK problem
  glp_prob* lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);

  // setup the objective function coefficients
  glp_add_cols(lp, lpparams.n);
  for (unsigned i=0; i< lpparams.n; i++)
    glp_set_obj_coef(lp, i+1, (double) lpparams.c[i]);

  // setup bounds on x
  const Real INF = std::numeric_limits<Real>::max();
  for (unsigned i=0; i< lpparams.n; i++)
  {
    // check for unbounded x
    if (lpparams.u.size() == 0 || lpparams.u[i] == INF)
    {
      // check for unbounded x
      if (lpparams.l.size() == 0 || lpparams.l[i] == -INF)
        glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
      else // lower bounded x
        glp_set_col_bnds(lp, i+1, GLP_LO, lpparams.l[i], 0.0);
    }
    else
    {
      // check for upper bounded x
      if (lpparams.l.size() == 0 || lpparams.l[i] == -INF)
        glp_set_col_bnds(lp, i+1, GLP_UP, 0.0, lpparams.u[i]);
      else // double-bounded x
        glp_set_col_bnds(lp, i+1, GLP_DB, lpparams.l[i], lpparams.u[i]);
    }
  }

  // setup number of constraints
  glp_add_rows(lp, lpparams.b.size() + lpparams.q.size());

  // setup equality constraints
  for (unsigned i=0; i< lpparams.b.size(); i++)
    glp_set_row_bnds(lp, i+1, GLP_FX, (double) lpparams.b[i], (double) lpparams.b[i]);

  // setup inequality constraints
  for (unsigned i=0; i< lpparams.q.size(); i++)
    glp_set_row_bnds(lp, lpparams.b.size()+i+1, GLP_LO, (double) lpparams.q[i], (double) lpparams.q[i]);

  // setup the coefficients
  unsigned sz = (lpparams.b.size() + lpparams.q.size()) * lpparams.n + 1;
  shared_array<int> ia(new int[sz]);
  shared_array<int> ja(new int[sz]);
  shared_array<double> ar(new double[sz]);

  // first setup the coefficients for A
  int m = 1;
  for (unsigned i=0; i< lpparams.b.size(); i++)
    for (unsigned j=0; j< lpparams.n; j++, m++)
    {
      ia[m] = (int) i+1;
      ja[m] = (int) j+1;
      ar[m] = (double) lpparams.A(i,j);
    }

  // now setup the coefficients for M
  for (unsigned i=0; i< lpparams.q.size(); i++)
    for (unsigned j=0; j< lpparams.n; j++, m++)
    {
      ia[m] = (int) i+1+lpparams.b.size();
      ja[m] = (int) j+1;
      ar[m] = (double) lpparams.M(i,j);
    }

  // load the matrix
  glp_load_matrix(lp, (lpparams.b.size()+lpparams.q.size())*lpparams.n, ia.get(), ja.get(), ar.get());

  // turn off terminal output
  glp_smcp sparams;
  glp_init_smcp(&sparams);
  sparams.msg_lev = GLP_MSG_OFF;

  // call the interior-point method
  bool result =  (glp_simplex(lp, &sparams) == 0 && glp_get_status(lp) == GLP_OPT);

  // get x out
  x.resize(lpparams.n);
  for (unsigned i=0; i< x.size(); i++)
    x[i] = (Real) glp_get_col_prim(lp, (int) i+1);

  // delete the GLPK problem
  glp_delete_prob(lp);

  return result;
  #else
  std::cerr << "Optimization::lp_simplex() - built without link to glpk!  Unable to solve LP!" << endl; 
  return false;
  #endif
}

/// Solves a convex quadratic program with an infeasible starting point
/**
 * \param A a matrix of equality constraints
 * \param b the RHS of the equality constraints
 * \param M a matrix of inequality constraints
 * \param q the RHS of the inequality constraints
 * \param x on input, an estimate of the solution; on output, a feasible point
 * \return <b>true</b> if successful, <b>false</b> otherwise
 */
bool Optimization::make_feasible_qp(const MatrixN& A, const VectorN& b, const MatrixN& M, const VectorN& q, VectorN& x, Real tol)
{
  // determine the dimension of x
  const unsigned xdim = A.columns();

  // determine the dimension of z
  const unsigned zdim = A.rows() + M.rows();

  // setup the vector c for linear programming
  VectorN c = VectorN::concat(VectorN::zero(xdim), VectorN::one(zdim));

  // compute Ax - b using the estimate
  VectorN Ax_m_b;
  if (x.size() == xdim)
    Ax_m_b = A*x - b;
  else
    Ax_m_b.set_zero(b.size());

  // compute q - Mx using the estimate
  VectorN q_m_Mx;
  if (x.size() == xdim)
    q_m_Mx = q - M*x;
  else
    q_m_Mx.set_zero(q.size());

  // setup the equality constraints
  MatrixN A_ = MatrixN::zero(A.rows(), xdim + zdim);
  A_.set_sub_mat(0, 0, A);
  for (unsigned i=0; i< A.rows(); i++)
    A_(i,xdim+i) = -sign(Ax_m_b[i]); 

  // setup the inequality constraints
  MatrixN M_ = MatrixN::zero(M.rows() + zdim, xdim + zdim);
  VectorN q_ = VectorN::zero(q.size() + zdim);
  M_.set_sub_mat(0, 0, M);
  q_.set_sub_vec(0, q);
  for (unsigned i=0; i< M.rows(); i++)
    M_(i,xdim+A.rows()+i) = M_(M.rows()+i,xdim+A.rows()+i) = (Real) 1.0;

  // setup feasible point for the linear program
  VectorN y(xdim + zdim);
  y.set_sub_vec(0, x);
  for (unsigned i=0; i< A.rows(); i++)
    y[xdim+i] = std::fabs(Ax_m_b[i]);
  for (unsigned i=0; i< M.rows(); i++)
    y[xdim+A.rows()+i] = std::max((Real) 0.0, q_m_Mx[i]);

  // setup the linear program parameters
  SAFESTATIC LPParams lpdata;
  lpdata.c.copy_from(c);
  lpdata.A.copy_from(A_);
  lpdata.b.copy_from(b);
  lpdata.M.copy_from(M_);
  lpdata.q.copy_from(q_);
  lpdata.l.resize(0);
  lpdata.u.resize(0);

  // solve the linear program
  bool result = lp_simplex(lpdata, y);
  if (!result)
    return false;

  // verify that solution is feasible
  if (c.dot(y) > tol)
    return false;

  // get the feasible point out
  y.get_sub_vec(0, xdim, x);

  FILE_LOG(LOG_OPT) << "Optimization::make_feasible_qp() entered" << endl;
  FILE_LOG(LOG_OPT) << " -- objective function value: " << c.dot(y) << endl;
  FILE_LOG(LOG_OPT) << " -- z: " << y.get_sub_vec(xdim, xdim+zdim) << endl;

  // indicate success
  return true;
}

/// Eliminates redundant equality constraints for a QP
void Optimization::eliminate_redundant_constraints(MatrixN& A)
{
  SAFESTATIC MatrixN Ar;

  // do Gaussian elimination
  Ar.copy_from(A);
  LinAlg::gauss_elim(Ar);

  // determine how many non-zero rows we have
  unsigned nz = 0;
  unsigned min_dim = std::min(Ar.rows(), Ar.columns());
  for (unsigned i=0; i< min_dim; i++)
    if (std::fabs(Ar(i,i)) > std::numeric_limits<Real>::epsilon())
      nz++;
    else
      break;

  // get the desired submatrix
  Ar.get_sub_mat(0,nz,0,Ar.columns(), A);
}

/// Eliminates redundant constraint from constraint equations for a QP or LP
void Optimization::eliminate_redundant_constraints(MatrixN& A, VectorN& b)
{
  SAFESTATIC MatrixN Ab;

  // augment A with b
  Ab.resize(A.rows(), A.columns()+1);
  Ab.set_sub_mat(0, 0, A);
  Ab.set_column(A.columns(), b);

  // do Gaussian elimination
  LinAlg::gauss_elim(Ab);

  // determine how many non-zero rows we have
  unsigned nz = 0;
  unsigned min_dim = std::min(Ab.rows(), Ab.columns());
  for (unsigned i=0; i< min_dim; i++)
    if (std::fabs(Ab(i,i)) > std::numeric_limits<Real>::epsilon())
      nz++;
    else
      break;

  // if we have sufficiently many non-zero rows, do nothing
  if (nz == A.rows())
    return;

  // put the Gaussian eliminated rows in A and b
  Ab.get_sub_mat(0, nz, 0, Ab.columns()-1, A);
  b.resize(nz);
  for (unsigned i=0; i< nz; i++)
    b[i] = Ab(i,Ab.columns()-1);
}

/// Hessian helper function for solving convex LCPs
void Optimization::qp_ip_hess(const VectorN& x, Real objscal, const VectorN& lambda, const VectorN& nu, MatrixN& H, void* data)
{
  // get the optimization data and QP data 
  const QPData& qpdata = *((QPData*) data);
  const OptParams& oparams = *((OptParams*) qpdata.data);

  // get G
  const MatrixN& G = qpdata.G;

  // copy precomputed Hessian
  H.copy_from(G) *= objscal;
}

/// Gradient helper function for solving convex LCPs
void Optimization::qp_ip_grad0(const VectorN& x, VectorN& g, void* data)
{
  // get the optimization data and QP data 
  const QPData& qpdata = *((QPData*) data);
  const OptParams& oparams = *((OptParams*) qpdata.data);

  // get QP data necessary for computing gradients 
  const MatrixN& G = qpdata.G;
  const VectorN& c = qpdata.c;
  const MatrixN& M = oparams.M;

  // set evaluation function gradient 
  G.mult(x, g) += c;
}

/// Function evaluation helper function for solving convex QPs
Real Optimization::qp_ip_f0(const VectorN& x, void* data)
{
  const Real S_BUFFER = std::numeric_limits<Real>::epsilon();

  // get the optimization data and QP data 
  const QPData& qpdata = *((QPData*) data);
  const OptParams& oparams = *((OptParams*) qpdata.data);

  // setup work vector
  SAFESTATIC VectorN workv;

  // get the QP data 
  const MatrixN& G = qpdata.G;
  const VectorN& c = qpdata.c;
  const MatrixN& M = oparams.M;
  const VectorN& q = oparams.q;

  // evaluate the function itself
  G.mult(x, workv) *= (Real) 0.5;
  workv += c;
  return x.dot(workv);
}

/// Solves a convex quadratic program using a primal-dual interior point method
/**
 * Minimzes: 0.5*x'*G*x + x'c subject to
 *           Ax = b
 *           Mx >= q
 * \param G a nxn symmetric, positive semi-definite matrix
 * \param c a nx1 vector
 * \param A a rxn matrix
 * \param b a rx1 vector
 * \param M a sxn matrix 
 * \param q a sx1 vector
 * \param x contains the optimal value on return; method can be "warm started"
 *        by setting x to a close to optimal value before calling.  x need not
 *        be a feasible point on entry.
 * \param eps the tolerance to solve the QP to using the interior-point solver
 * \param eps_feas the feasibility tolerance for the IP solver
 * \param max_iterations the maximum number of iterations for the IP solver
 * \return true if successful, false otherwise
 */
bool Optimization::qp_convex_ip(const MatrixN& G, const VectorN& c, OptParams& oparams, VectorN& x)
{
  // get number of variables
  const unsigned n = c.size();

  if (G.rows() != n)
    throw MissizeException();
  if (G.columns() != n)
    throw NonsquareMatrixException();

  // setup x 
  x.set_zero(n);

  // setup convex optimization parameters
  SAFESTATIC OptParams cparams;
  cparams.n = n;
  cparams.m = 0;
  cparams.r = 0;
  cparams.f0 = &qp_ip_f0;
  cparams.grad0 = &qp_ip_grad0;
  cparams.hess = &qp_ip_hess;
  cparams.max_iterations = oparams.max_iterations;
  cparams.eps = oparams.eps;
  cparams.eps_feas = oparams.eps_feas;
  cparams.zero_tol = oparams.zero_tol;
  cparams.alpha = oparams.alpha;
  cparams.beta = oparams.beta;
  cparams.mu = oparams.mu;
  cparams.M.copy_from(oparams.M);
  cparams.q.copy_from(oparams.q); 
  cparams.A.copy_from(oparams.A);
  cparams.b.copy_from(oparams.b); 
  cparams.lb.copy_from(oparams.lb);
  cparams.ub.copy_from(oparams.ub);

  // create QP data
  SAFESTATIC QPData qpd;
  qpd.G.copy_from(G);
  qpd.c.copy_from(c);
  qpd.data = (void*) &cparams;
  cparams.data = &qpd;

  // first, make feasible...
  if (!make_feasible_convex(cparams, x))
  {
    FILE_LOG(LOG_OPT) << "Optimization::qp_convex_ip() - could not find feasible point!" << endl;
    return false;
  }

  // optimize
  if (!optimize_convex_pd(cparams, x))
  {
    FILE_LOG(LOG_OPT) << "Optimization::qp_convex_ip() - convex optimization failed!" << endl;
    return false; 
  }

  if (LOGGING(LOG_OPT))
  {
    // check constraints
    VectorN w = oparams.M*x - oparams.q;
    Real minw = *std::min_element(w.begin(), w.end());
    Real eq = (oparams.A*x - oparams.b).norm();
    FILE_LOG(LOG_OPT) << "Optimization::qp_convex_ip() successful!" << endl;
    FILE_LOG(LOG_OPT) << "  -- maximum inequality constraint violation: " << -minw << endl;
    FILE_LOG(LOG_OPT) << "  -- equality constraint norm (lower is better): " << eq << endl;
  }

  return true;
}

/// Determines whether a convex quadratic program is feasible (for early termination)
bool Optimization::qp_convex_activeset_infeas_tcheck(const VectorN& x, void* data)
{
  SAFESTATIC VectorN y, feas;

  // get the QP parameters
  const OptParams& oparams = *((OptParams*) data);

  FILE_LOG(LOG_OPT) << "qp_convex_activeset_infeas_tcheck(): " << endl;
  FILE_LOG(LOG_OPT) << " x: " << x << endl; 
  FILE_LOG(LOG_OPT) << " lb: " << oparams.lb << endl; 
  FILE_LOG(LOG_OPT) << " ub: " << oparams.ub << endl; 
  FILE_LOG(LOG_OPT) << " ci: " << feas << std::endl;

  // first check lower and upper bounds
  x.get_sub_vec(0, oparams.M.columns(), y);

  // check inequality feasibility
  oparams.M.mult(y, feas) -= oparams.q;
  FILE_LOG(LOG_OPT) << " ci: " << feas << std::endl;
  Real min_ci = (feas.size() > 0) ? *std::min_element(feas.begin(), feas.end()): (Real) 1.0;

  // check equality feasibility
  Real max_ce = (Real) 0.0;
  if (oparams.b.size() > 0)
  {
    oparams.A.mult(y, feas) -= oparams.b;
    FILE_LOG(LOG_OPT) << "ce: " << feas << std::endl;
    max_ce = feas.norm_inf();
  }

  FILE_LOG(LOG_OPT) << "minimum ci feasibility: " << min_ci << std::endl; 
  FILE_LOG(LOG_OPT) << "maximum ce feasibility: " << max_ce << std::endl; 

  return (min_ci > -oparams.eps_feas && max_ce < oparams.eps_feas);
}

/// Solves a convex quadratic program for the case where initial point is infeasible
/**
 * This method is described in Nocedal and Wright, 1999, pg. 537.
 */
void Optimization::qp_convex_activeset_infeas(const MatrixN& G, const VectorN& c, Real upsilon, OptParams& qparams, VectorN& x, bool hot_start)
{
  SAFESTATIC VectorN c2, q2, y, infeas_eq, infeas_ineq;
  SAFESTATIC MatrixN M2, A2;
  SAFESTATIC MatrixN G2;
  SAFESTATIC vector<bool> eq_indices, ineq_indices;

  if (G.rows() != qparams.n)
    throw MissizeException();
  if (G.columns() != qparams.n)
    throw NonsquareMatrixException();

  // don't do anything if we don't need to 
  if (qparams.n == 0)
    return;

  // determine infeasibility
  qparams.M.mult(x, infeas_ineq) -= qparams.q;
  qparams.A.mult(x, infeas_eq) -= qparams.b;

  // determine number of new variables for equality constraints
  eq_indices.resize(infeas_eq.size());
  unsigned NEQ = 0;
  for (unsigned i=0; i< infeas_eq.size(); i++)
  {
    eq_indices[i] = (std::fabs(infeas_eq[i]) > std::numeric_limits<Real>::epsilon() * std::max(std::fabs(qparams.b[i]), (Real) 1.0));
    if (eq_indices[i])
      NEQ++;
  }

  // determine number of new variables for inequality constraints
  ineq_indices.resize(infeas_ineq.size());
  unsigned NINEQ = 0;
  for (unsigned i=0; i< infeas_ineq.size(); i++)
  {
    ineq_indices[i] = (infeas_ineq[i] < std::numeric_limits<Real>::epsilon() * -std::max(std::fabs(qparams.q[i]), (Real) 1.0));
    if (ineq_indices[i])
      NINEQ++;
  }

  // compute new variables necessary 
  const unsigned NAUG = (NINEQ > 0) ? NEQ*2+1 : NEQ*2;
  const unsigned NVARS = qparams.n + NAUG;

  FILE_LOG(LOG_OPT) << "-- executing infeasible active set method" << endl;
  FILE_LOG(LOG_OPT) << "-- upsilon: " << upsilon << endl;

  // verify that it is necessary to run this method
  if (NAUG == 0)
  {
    FILE_LOG(LOG_OPT) << "-- no infeasibilities detected!" << endl;
    return;
  }

  // setup a new QP optimization
  SAFESTATIC OptParams qparams2;
  qparams2.iterations = 0;
  qparams2.r = qparams2.m = 0;
  qparams2.max_iterations = qparams.max_iterations;
  qparams2.n = NVARS;
  qparams2.b.copy_from(qparams.b);
  qparams2.eps_feas = qparams.eps_feas;
  qparams2.eps = qparams.eps;
  qparams2.zero_tol = qparams.zero_tol;
  qparams2.tcheck = &qp_convex_activeset_infeas_tcheck;
  qparams2.data = (void*) &qparams;
  qparams2.ub.resize(0);

  // augment c
  const unsigned S_IDX = qparams.n;
  c2.resize(qparams2.n);
  c2.set_sub_vec(0, c);
  std::fill_n(c2.begin()+qparams.n, NAUG, upsilon * NAUG);

  // setup q
  qparams2.q.copy_from(qparams.q);

  // determine maximum inequality violation
  Real* min_ineq = std::min_element(infeas_ineq.begin(), infeas_ineq.end());

  // augment M -- this part just adds s variable to standard inequalities
  qparams2.M.resize(qparams.M.rows(), NVARS);
  qparams2.M.set_sub_mat(0,0,qparams.M);
  BlockIterator bi = qparams2.M.block_start(0, qparams.M.rows(), qparams.n, NVARS);
  std::fill_n(bi, NVARS*qparams.M.rows(), (Real) 0.0);
  for (unsigned i=0; i< qparams.q.size(); i++)
    if (ineq_indices[i])
      qparams2.M(i,qparams.n) = (Real) 1.0;

  // setup lb and ub - standard variables
  qparams2.lb.resize(NVARS);
  const Real INF = std::numeric_limits<Real>::max();
  if (qparams.lb.size() > 0)
    qparams2.lb.set_sub_vec(0, qparams.lb);
  else
    std::fill_n(qparams2.lb.begin(), qparams.n, -INF); 
  if (qparams.ub.size() > 0)
  {
    qparams2.ub.resize(NVARS);
    qparams2.ub.set_sub_vec(0, qparams.ub);
    std::fill_n(qparams2.ub.begin()+qparams.n, NVARS, INF); 
  }
  else
    qparams2.ub.resize(0);

  // setup non-negativity constraints on s, v, w
  std::fill(qparams2.lb.begin()+qparams.n, qparams2.lb.end(), (Real) 0.0);

  // prepare to determine initial variable values
  y.resize(NVARS);
  y.set_sub_vec(0, x);
  if (NINEQ > 0)
    // next line ensures that s is sufficiently large
    y[S_IDX] = std::max((Real) 1.0, nextafter(-*min_ineq, INF));

  // augment A (and determine initial values for v,w)
  const unsigned VW_IDX_START = (NINEQ == 0) ? qparams.n : qparams.n+1;
  qparams2.b.copy_from(qparams.b);
  qparams2.A.resize(qparams.b.size(), NVARS);
  qparams2.A.set_sub_mat(0,0,qparams.A);
  bi = qparams2.A.block_start(0, qparams.A.rows(), qparams.n, NVARS);
  std::fill_n(bi, qparams.A.rows()*NAUG, (Real) 0.0);
  for (unsigned i=0, j=VW_IDX_START; i< qparams.b.size(); i++)
    if (eq_indices[i])
    {
      qparams2.A(i, j+0) = (Real) -1.0;
      qparams2.A(i, j+1) = (Real) 1.0;
      if (infeas_eq[i] > (Real) 0.0)
      {
        y[j+0] = infeas_eq[i];
        y[j+1] = (Real) 0.0;
      }
      else
      {
        y[j+0] = (Real) 0.0;
        y[j+1] = -infeas_eq[i];
      }
    } 

  // augment G
  G2.resize(NVARS,NVARS);
  G2.set_sub_mat(0,0,G);
  BlockIterator bs = G2.block_start(qparams.n, qparams.n+NAUG, 0, qparams.n);
  BlockIterator be = G2.block_end(qparams.n, qparams.n+NAUG, 0, qparams.n);
  std::fill(bs, be, (Real) 0.0);
  bs = G2.block_start(0, qparams.n, qparams.n, qparams.n+NAUG);
  be = G2.block_end(0, qparams.n, qparams.n, qparams.n+NAUG);
  std::fill(bs, be, (Real) 0.0);
  bs = G2.block_start(qparams.n, qparams.n+NAUG, qparams.n, qparams.n+NAUG);
  be = G2.block_end(qparams.n, qparams.n+NAUG, qparams.n, qparams.n+NAUG);
  std::fill(bs, be, (Real) 0.0);

  // run the active set method
  qp_convex_activeset(G2, c2, qparams2, y, hot_start);

  // get the result out
  y.get_sub_vec(0, qparams.n, x);

  FILE_LOG(LOG_OPT) << "infeasible QP method()" << endl;
  FILE_LOG(LOG_OPT) << "final s: " << y.get_sub_vec(qparams.n, y.size()) << endl;
  FILE_LOG(LOG_OPT) << "M*x - q: " << (qparams.M * x - qparams.q) << endl;
}

/// Computes whether AMr (free) is full rank using the QR factorization
bool WorkingSet::full_rank() const
{
  unsigned k = std::min(_R.rows(), _R.columns());

  // first, determine the minimum and maximum coefficients 
  const Real INF = std::numeric_limits<Real>::max();
  Real min_rii = INF, max_rii = (Real) -INF;
  for (unsigned i=0; i< k; i++)
  {
    min_rii = std::min(min_rii, std::fabs(_R(i,i)));
    max_rii = std::max(max_rii, std::fabs(_R(i,i)));
  }

  // now check whether all values w/in tolerance
  return (min_rii > std::max(_R.rows(), _R.columns()) * max_rii * std::numeric_limits<Real>::epsilon());
}

/// Resets with a working set
void WorkingSet::reset(const MatrixN& A, const MatrixN& M, const std::vector<bool>& vworking, const std::vector<bool>& working)
{
  // initialize _updates
  _updates = 0;

  // get n and m
  _n = A.columns();
  _m = M.rows();
  assert(_n == M.columns());

  // clear vectors
  _working = working;
  _vworking = vworking;
  _free.clear();
  assert(_working.size() == _m);
  assert(_vworking.size() == _n);

  // compute number of working constraints
  _nworking = 0;
  for (unsigned i=0; i< _m; i++)
    if (_working[i])
      _nworking++;

  // make appropriate variables free
  _nvworking = 0;
  for (unsigned i=0; i< _n; i++)
  {
    if (_vworking[i])
      _nvworking++;
    else
      _free.push_back(i);
  }

  // initialize reduced A
  _Ar.resize(0, _n);

  // determine reduced A
  for (unsigned i=0; i< A.rows(); i++)
  {
    // get row of A
    A.get_row(i, _workv);

    // resize new Ar 
    _workM.resize(_Ar.rows()+1, _n);
    _workM.set_sub_mat(0,0,_Ar);
    _workM.set_row(_Ar.rows(), _workv);
    _workM2.copy_from(_workM);
    if (LinAlg::calc_rank(_workM2) < _workM.rows())
      continue;

    // setup Ar
    _Ar.copy_from(_workM);
  }

  // copy M and q
  _M = M;

  // verify that reduced A is not too large
  const unsigned R = _Ar.rows();
  if (_Ar.rows() >= _n)
    throw std::runtime_error("Too many rows in reduced A!");

  // reform AMr
  reform_AMr();

  // compute the QR factorization of AMr
  MatrixN::transpose(_AMr_free, _R);
  LinAlg::factor_QR(_R, _Q);

  // setup the range and nullspace
  if (_AMr_free.rows() == 0)
  {
    const unsigned NFREE = _AMr_free.columns();
    _Z.set_identity(NFREE);
  }
  else
  {
    unsigned ns = _AMr_free.columns() - _AMr_free.rows();
    _Q.get_sub_mat(0, _Q.rows(), _Q.columns()-ns, _Q.columns(), _Z);
  }
}

/// Sets up the working set 
void WorkingSet::reset(const MatrixN& A, const MatrixN& M)
{
  SAFESTATIC vector<bool> empty_var, empty_constraints;

  // setup empty working sets
  empty_var.resize(M.columns());
  empty_constraints.resize(M.rows());
  std::fill_n(empty_var.begin(), M.columns(), false);
  std::fill_n(empty_constraints.begin(), M.rows(), false);

  // call standard reset function
  reset(A, M, empty_var, empty_constraints);  
}

/// Computes whether the working set is full (no more constraints can be added)
bool WorkingSet::full() const
{
  // fast check
  assert(_nvworking <= _n - _nvworking);
  return (_nworking == _n - _nvworking - 1);
}

// attempts to add a variable to the working set; returns 'true' if successful
bool WorkingSet::add_var_to_working_set(unsigned var)
{
  // look for fast exit
  if (full())
    return false;

  // set vworking correctly
  assert(!_vworking[var]);
  _vworking[var] = true;
  _nvworking++;

  // update free
  vector<unsigned>::iterator free_iter = std::lower_bound(_free.begin(), _free.end(), var);
  assert(free_iter != _free.end());
  unsigned free_index = std::distance(_free.begin(), free_iter);
  _free.erase(free_iter);

  // reform AMr
  reform_AMr();

  // check whether we need a total recomputation of Q and R or just an update
  bool low_rank_update = true;
  if (_updates == max_updates)
  {
    // indicate no low-rank update was performed
    low_rank_update = false;

    // recompute QR factorization of AMr (free variables only)
    MatrixN::transpose(_AMr_free, _R);
    LinAlg::factor_QR(_R, _Q);
  }
  else
  {
    // save current Q and R
    _workM.copy_from(_Q);
    _workM2.copy_from(_R);

    // do the rank-1 update here
    LinAlg::update_QR_delete_rows(_Q, _R, free_index, 1);
  }

  // determine whether _AMr is full rank
  if (!full_rank())
  {
    // not full rank, revert old _Q and _R (if necessary)
    if (low_rank_update)
    {
      _Q.copy_from(_workM);
      _R.copy_from(_workM2);
    }

    // revert vworking
    _vworking[var] = false;
    _nvworking--;

    // revert free
    _free.insert(std::lower_bound(_free.begin(), _free.end(), var), var);

    // indicate failure
    return false;
  }

  // form nullspace 
  if (_AMr_free.rows() == 0)
  {
    // must form nullspace ourself
    const unsigned NFREE = _AMr_free.columns();
    _Z.set_identity(NFREE);
  }
  else
  {
    unsigned ns = _AMr_free.columns() - _AMr_free.rows();
    _Q.get_sub_mat(0, _Q.rows(), _Q.columns()-ns, _Q.columns(), _Z);
  }

  // update _updates
  if (low_rank_update)
    _updates++;
  else
    _updates = 0;

  // indicate success
  return true;
}

// attempts to add a constraint to the working set; returns 'true' if successful
bool WorkingSet::add_constraint_to_working_set(unsigned constraint)
{
  // look for fast exit
  if (full())
    return false;

  // set working correctly
  assert(!_working[constraint]);
  _working[constraint] = true;
  _nworking++;

  // reform AMr
  reform_AMr();

  // check whether we need a total recomputation of Q and R or just an update
  bool low_rank_update = true;
  if (_updates == max_updates)
  {
    // indicate no low-rank update was performed
    low_rank_update = false;

    // recompute QR factorization of AMr (free variables only)
    MatrixN::transpose(_AMr_free, _R);
    LinAlg::factor_QR(_R, _Q);
  }
  else
  {
    // save current Q and R
    _workM.copy_from(_Q);
    _workM2.copy_from(_R);

    // do the rank-1 update
    unsigned cidx = _Ar.rows();
    for (unsigned i=0; i< constraint; i++)
      if (_working[i])
        cidx++;

    // get the constraint
    _AMr_free.get_sub_mat(cidx, cidx+1, 0, _AMr_free.columns(), _workM);
    MatrixN::transpose(_workM, _workM2);
    LinAlg::update_QR_insert_cols(_Q, _R, _workM, cidx);
  }

  // determine whether _AMr is full rank
  if (!full_rank())
  {
    // not full rank, revert old _Q and _R (if necessary)
    if (low_rank_update)
    {
      _Q.copy_from(_workM);
      _R.copy_from(_workM2);
    }

    // revert vworking
    _working[constraint] = false;
    _nworking--;

    // indicate failure
    return false;
  }

  // form nullspace 
  if (_AMr_free.rows() == 0)
  {
    // must form nullspace ourself
    const unsigned NFREE = _AMr_free.columns();
    _Z.set_identity(NFREE);
  }
  else
  {
    unsigned ns = _AMr_free.columns() - _AMr_free.rows();
    _Q.get_sub_mat(0, _Q.rows(), _Q.columns()-ns, _Q.columns(), _Z);
  }

  // update _updates
  if (low_rank_update)
    _updates++;
  else
    _updates = 0;

  // indicate success
  return true;
}

// removes the variable from the working set
void WorkingSet::remove_var_from_working_set(unsigned var)
{
  // set vworking correctly
  assert(_vworking[var]);
  _vworking[var] = false;
  _nvworking--;

  // update free
  vector<unsigned>::iterator free_iter = std::lower_bound(_free.begin(), _free.end(), var);
  unsigned free_index = std::distance(_free.begin(), free_iter);
  _free.insert(free_iter, var);

  // reform AMr
  reform_AMr();

  // check whether we need a total recomputation of Q and R or just an update
  if (_updates++ == max_updates)
  {
    // reset _updates 
    _updates = 0;

    // recompute QR factorization of AMr (free variables only)
    MatrixN::transpose(_AMr_free, _R);
    LinAlg::factor_QR(_R, _Q);
  }
  else
  {
    // get the appropriate column of AMr
    _AMr_free.get_sub_mat(0, _AMr_free.rows(), free_index, free_index+1, _workM);
    MatrixN::transpose(_workM, _workM2);

    // do the rank-1 update
    LinAlg::update_QR_insert_rows(_Q, _R, _workM, free_index);
  }

  // form nullspace 
  if (_AMr_free.rows() == 0)
  {
    // must form nullspace ourself
    const unsigned NFREE = _AMr_free.columns();
    _Z.set_identity(NFREE);
  }
  else
  {
    unsigned ns = _AMr_free.columns() - _AMr_free.rows();
    _Q.get_sub_mat(0, _Q.rows(), _Q.columns()-ns, _Q.columns(), _Z);
  }
}

// removes the inequality constraint from the working set
void WorkingSet::remove_constraint_from_working_set(unsigned constraint)
{
  // set working correctly
  assert(_working[constraint]);
  _working[constraint] = false;
  _nworking--;

  // reform AMr_free and AMr
  reform_AMr();

  // check whether we need a total recomputation of Q and R or just an update
  if (_updates++ == max_updates)
  {
    // reset _updates 
    _updates = 0;

    // recompute QR factorization of AMr (free variables only)
    MatrixN::transpose(_AMr_free, _R);
    LinAlg::factor_QR(_R, _Q);
  }
  else
  {
    // determine the constraint index
    unsigned cidx = _Ar.rows();
    for (unsigned i=0; i< constraint; i++)
      if (_working[i])
        cidx++;

    // do the rank-1 update here
    LinAlg::update_QR_delete_cols(_Q, _R, cidx, 1);
  }

  // form nullspace 
  if (_AMr_free.rows() == 0)
  {
    // must form nullspace ourself
    const unsigned NFREE = _AMr_free.columns();
    _Z.set_identity(NFREE);
  }
  else
  {
    unsigned ns = _AMr_free.columns() - _AMr_free.rows();
    _Q.get_sub_mat(0, _Q.rows(), _Q.columns()-ns, _Q.columns(), _Z);
  }
} 

// Reforms the AMr matrix
void WorkingSet::reform_AMr()
{
  // determine R -- number of linear equality constraints
  const unsigned R = _Ar.rows();
  const unsigned NFREE = _n - _nvworking;

  // setup equality constraints portion of AMr
  _AMr.resize(R+_nworking, _n);
  _AMr_free.resize(R+_nworking, NFREE);
  _Ar.select_columns(_free.begin(), _free.end(), _workM);
  BlockIterator Ar_bs = _Ar.block_start(0, _Ar.rows(), 0, NFREE);
  BlockIterator Ar_be = _Ar.block_end(0, _Ar.rows(), 0, NFREE);
  BlockIterator AMr_bs = _AMr.block_start(0, _Ar.rows(), 0, NFREE); 
  std::copy(Ar_bs, Ar_be, AMr_bs); 
  _AMr_free.set_sub_mat(0,0, _workM);

  // set inequality constrained portions of AMr
  for (unsigned i=0,j=R; i< _m; i++)
  {
    if (!_working[i])
      continue;
    _M.get_row(i, _workv);
    _AMr.set_row(j, _workv);
    _workv.select(_free.begin(), _free.end(), _workv2);
    _AMr_free.set_row(j++, _workv2);
  }
}

/// Equilibrate a matrix
void Optimization::equilibrate(MatrixN& A)
{
  if (A.rows() == 0)
    return;

  // scale
  for (unsigned i=0; i< A.rows(); i++)
  {
    BlockIterator bi = A.block_start(i, i+1, 0, A.columns());
    pair<BlockIterator, BlockIterator> mmax = boost::minmax_element(bi, bi.end());
    Real scalar = std::max(-*mmax.first, *mmax.second);
    assert(scalar > (Real) 0.0);
    Real inv_scalar = (Real) 1.0/scalar;
    std::transform(bi, bi.end(), bi, std::bind2nd(std::multiplies<Real>(), inv_scalar));
  }
}

/// Equilibrates a linear system 
void Optimization::equilibrate(MatrixN& A, VectorN& b)
{
  if (A.rows() != b.size())
    throw MissizeException();

  if (A.rows() == 0)
    return;

  // scale
  for (unsigned i=0; i< A.rows(); i++)
  {
    BlockIterator bi = A.block_start(i, i+1, 0, A.columns());
    pair<BlockIterator, BlockIterator> mmax = boost::minmax_element(bi, bi.end());
    Real scalar = std::max(std::max(-*mmax.first, *mmax.second), std::fabs(b[i]));
    assert(scalar > (Real) 0.0);
    Real inv_scalar = (Real) 1.0/scalar;
    std::transform(bi, bi.end(), bi, std::bind2nd(std::multiplies<Real>(), inv_scalar));
    b[i] *= inv_scalar; 
  }
}

/// Solves a convex quadratic program using a primal active set method
/**
 * Minimzes: 0.5*x'*G*x + x'c subject to
 *           Ax = b
 *           Mx >= q
 * \param G a nxn symmetric, positive semi-definite or positive-definite matrix
 * \param c a nx1 vector
 * \param A a rxn matrix (r < n); if r >= n, then A must contain redundant
 *          rows (and must then be factorized to reduce its size), is
 *          infeasible, or contains a unique point. A should be of full row
 *          rank
 * \param b a rx1 vector
 * \param M a sxn matrix 
 * \param q a sx1 vector
 * \param x contains the optimal value on return; method can be "warm started"
 *        by setting x to a close to optimal value before calling.  x must be
 *        a feasible point on entry.
 * \note method comes via [Nocedal, 2006], pp. 468-480
 */
void Optimization::qp_convex_activeset(const MatrixN& G, const VectorN& c, OptParams& qparams, VectorN& x, bool hot_start)
{
  const Real ZERO_TOL = qparams.zero_tol * qparams.zero_tol; 
  const Real INF = std::numeric_limits<Real>::max();

  if (G.rows() != qparams.n)
    throw MissizeException();
  if (G.columns() != qparams.n)
    throw NonsquareMatrixException();

  // get A, b, M, q
  VectorN& lb = qparams.lb;
  VectorN& ub = qparams.ub;
  SAFESTATIC WorkingSet ws; 
  SAFESTATIC MatrixN M, A;  // scaled M and A
  SAFESTATIC VectorN q, b;  // scaled q and b
  SAFESTATIC vector<bool> working, vworking;
  const vector<unsigned>& free = ws.get_free_indices();
  const MatrixN& Z = ws.get_null();
  const MatrixN& AR_Q = ws.get_Q();
  const MatrixN& AR_R = ws.get_R();
  const MatrixN& AMr_free = ws.get_A_M_free();
  const MatrixN& AMr_full = ws.get_A_M();

  // determine n, r, s
  const unsigned N = G.rows();

  FILE_LOG(LOG_OPT) << "qp_convex_activeset() entered" << endl;
  FILE_LOG(LOG_OPT) << "G: " << endl << G;
  FILE_LOG(LOG_OPT) << "c: " << c << endl;
  FILE_LOG(LOG_OPT) << "lb: " << lb << endl;
  FILE_LOG(LOG_OPT) << "ub: " << ub << endl;
  FILE_LOG(LOG_OPT) << "M: " << endl << M;
  FILE_LOG(LOG_OPT) << "q: " << q << endl;
  FILE_LOG(LOG_OPT) << "A: " << endl << A;
  FILE_LOG(LOG_OPT) << "b: " << b << endl;

  // trivial exit
  if (N == 0)
    return;

  // check things...
  if (N != c.size())
    throw MissizeException();
  if (A.rows() != b.size())
    throw MissizeException();
  if (A.columns() != N && A.rows() != 0)
    throw MissizeException();
  if (qparams.M.rows() != qparams.q.size())
    throw MissizeException();
  if (lb.size() != N)
  {
    if (lb.size() == 0)
      lb.set_one(N) *= -INF;
    else
      throw MissizeException();
  }
  if (ub.size() != N)
  {
    if (ub.size() == 0)
      ub.set_one(N) *= INF;
    else
      throw MissizeException();
  }
  if (x.size() != N)
  {
    if (x.size() == 0)
      x.set_zero(N);
    else
      throw MissizeException();
  }

  // make x respect bounds constraints
  for (unsigned i=0; i< N; i++)
    if (x[i] < lb[i])
      x[i] = lb[i];

  // verify user understands what this method is capable of
  if (qparams.m > 0)
    throw std::runtime_error("qp_convex_activeset() incapable of handling nonlinear inequality constraints!");
  if (qparams.r > 0)
    throw std::runtime_error("qp_convex_activeset() incapable of handling nonlinear equality constraints!");

  // copy M,A,q
  M.copy_from(qparams.M);
  A.copy_from(qparams.A);
  q.copy_from(qparams.q);
  const unsigned S = M.rows();

  // equilibrate (normalize M,q and A) 
  equilibrate(M, q);
  equilibrate(A);

  // setup working variables
  SAFESTATIC MatrixN Gf, Lhs;
  SAFESTATIC MatrixN workM;
  SAFESTATIC VectorN workv, workv2, dx, dxf, lambda, vlambda, grad, gradf, pz;
  SAFESTATIC vector<bool> blocking;
  SAFESTATIC vector<Real> alphas;

  // setup blocking vector
  blocking.resize(S+N);

  // evaluate the objective function
  G.mult(x, workv) *= (Real) 0.5;
  workv += c;
  Real last_eval = x.dot(workv);

  // reset the working set
  if (hot_start)
  {
    unsigned total_working = A.rows();

    // determine the variables in the working set
    vworking.resize(N);
    for (unsigned i=0; i< N && total_working < N-1; i++)
    {
      const Real TOL = std::numeric_limits<Real>::epsilon() * std::max((Real) 1.0, std::fabs(x[i]));
      if (std::fabs(x[i] - lb[i]) < TOL ||
          std::fabs(ub[i] - x[i]) < TOL)
      {
        total_working++;
        vworking[i] = true;
      }
      else
        vworking[i] = false;
    }

    // determine working set (active lin. indep. inequality constraints)
    M.mult(x, workv) -= q;
    working.resize(S);
    for (unsigned i=0; i< S && total_working < N-1; i++)
    {
      const Real TOL = std::numeric_limits<Real>::epsilon() * std::max((Real) 1.0, std::fabs(q[i]));
      working[i] = (workv[i] < TOL);
      if (working[i])
        total_working++;
    }

    // now, reset the working set
    ws.reset(A, M, vworking, working);
  }
  else
    ws.reset(A, M); 

  // get number of equalities
  const unsigned R = ws.num_equalities();

  // do the active set method
  for (qparams.iterations = 0; qparams.iterations < qparams.max_iterations; qparams.iterations++)
  {
    // output working set
    if (LOGGING(LOG_OPT))
    {
      std::ostringstream str, str2;
      str << "working set: ";
      for (unsigned j=0; j< S; j++)
        if (ws.is_constraint_working(j))
          str << " " << j;
      str2 << "variable working set: ";
      for (unsigned j=0; j< N; j++)
        if (ws.is_var_working(j))
          str2 << " " << j;
      FILE_LOG(LOG_OPT) << str.str() << endl;
      FILE_LOG(LOG_OPT) << str2.str() << endl;
    }

    // compute gradient
    G.mult(x, grad) += c;

    // get reduced matrices and vectors
    grad.select(free.begin(), free.end(), gradf);
    G.select_square(free.begin(), free.end(), Gf);

    // compute the left hand side matrix Z'*Gf*Z
    Z.transpose_mult(Gf, workM);
    workM.mult(Z, Lhs);

    // compute the rhs vector -Z'*g
    Z.transpose_mult(gradf, pz);
    pz.negate();

    // condition and solve the system
    condition_and_factor_PD(Lhs);
    LinAlg::solve_chol_fast(Lhs, pz);
    Z.mult(pz, dxf);
    dx.set_zero(N);
    dx.set(free.begin(), free.end(), dxf);

    // G*Z*pz + A'*lambda = -g
    // A'*lambda = -g - Gf*Z*pz
    // compute lambda for linear inequality constraints by solving (using least
    // squares): Ar' * lambda = -gradf - Gf * Z * pz
    // equivalent to: Qr * Rr * lambda = gradf - Gf * Z * pz
    const unsigned NVARS = std::min(AR_Q.rows(), AR_R.columns());
    Z.mult(pz, workv);
    Gf.mult(workv, workv2);
    workv2 += gradf; 
    workv2.negate();   // this line concludes the right hand side
    lambda.set_zero(NVARS);
    if (NVARS > 0)
    {
      CBLAS::gemv(CblasTrans, AR_Q.rows(), NVARS, AR_Q, AR_Q.rows(), workv2, 1, (Real) 1.0, (Real) 0.0, lambda, 1);
      CBLAS::trsv(CblasUpper, CblasNoTrans, (int) NVARS, AR_R, AR_R.rows(), lambda, 1);  // triangular solve rhs 
    }

    // compute vlambda: G*dx + A'*lambda + V'*vlambda = -g
    // V'*vlambda = -g + -G*dx - A'*lambda
    AMr_full.transpose_mult(lambda, workv2);
    G.mult(dx, workv);
    workv += grad;
    workv += workv2;
    workv.negate();
    vlambda.resize(N - free.size());
    for (unsigned i=0, j=0; i< N; i++)
      if (!std::binary_search(free.begin(), free.end(), i))
        vlambda[j++] = workv[i];

    // negate lambda variables 
    lambda.negate();
    vlambda.negate();

    FILE_LOG(LOG_OPT) << "x: " << x << endl;
    FILE_LOG(LOG_OPT) << "gradient: " << grad << endl;
    FILE_LOG(LOG_OPT) << "gradient (reduced): " << gradf << endl;
    FILE_LOG(LOG_OPT) << "G (reduced): " << endl << Gf;
    FILE_LOG(LOG_OPT) << "homogeneous solution: " << pz << endl;
    FILE_LOG(LOG_OPT) << "dx (reduced): " << dxf << endl;
    FILE_LOG(LOG_OPT) << "constraint lambda: " << lambda << endl;
    FILE_LOG(LOG_OPT) << "variable lambda: " << vlambda << endl;

    // if dx is zero, examine the Lagrange multipliers
    if (dx.norm() < std::numeric_limits<Real>::epsilon() * std::max(x.norm(), (Real) 1.0))
    {
      FILE_LOG(LOG_OPT) << "dx = 0; examining Lagrange multipliers" << endl;

      // if the working set is empty, we're done...
      if (ws.get_total_working() == 0)
      {
        FILE_LOG(LOG_OPT) << "  -- working set is empty -> optimized!" << endl;
        return;
      }

      // find the minimum lambda and minimum vlambda
      Real* min_lambda = std::min_element(lambda.begin()+R, lambda.end());
      Real* min_vlambda = std::min_element(vlambda.begin(), vlambda.end());
      if (min_lambda == lambda.end())
        min_lambda = (Real*) &INF;
      if (min_vlambda == vlambda.end())
        min_vlambda = (Real*) &INF;
      Real min_min = std::min(*min_lambda, *min_vlambda);
      if (min_min > -std::numeric_limits<Real>::epsilon())
      {
        FILE_LOG(LOG_OPT) << "  -- minimum lambda = " << min_min << " -> optimized!" << endl;
        return;
      }

      FILE_LOG(LOG_OPT) << "  -- minimum lambda = " << min_min << endl;

      // remove the minimum constraint from the working set
      if (*min_lambda < *min_vlambda)
      {
        unsigned min_constraint = std::distance(lambda.begin()+R, min_lambda);
        for (unsigned j=0; j< S; j++)
          if (ws.is_constraint_working(j) && min_constraint-- == 0)
          {
            FILE_LOG(LOG_OPT) << "  -- removing constraint " << j << " from working set" << endl;

            // remove the constraint from the working set
            ws.remove_constraint_from_working_set(j);

            // quit looping
            break;
          }
      }
      else
      {
        unsigned min_constraint = std::distance(vlambda.begin(), min_vlambda);
        for (unsigned j=0; j< N; j++)
        {
          if (ws.is_var_working(j) && min_constraint-- == 0)
          {
            FILE_LOG(LOG_OPT) << "  -- removing variable constraint " << j << " from working set" << endl;

            // remove the bound from the working set
            ws.remove_var_from_working_set(j);

            // quit looping 
            break;
          }
        }
      }

      // continue iterating
      continue;
    }
    else
    {
      // compute all values of alpha
      alphas.resize(N+S);
      for (unsigned i=0; i< N+S; i++)
        alphas[i] = (Real) 1.0;

      // set all constraints to non-blocking
      for (unsigned j=0; j< N+S; j++)
        blocking[j] = false;

      // loop over x bounds first
      for (unsigned j=0, k=0; j< N; j++, k++)
      {
        // constraint cannot be blocking if it is in the working set
        if (!ws.is_var_working(j))
        {
          alphas[k] = (Real) 1.0;
          if (dx[j] > (Real) 0.0)
            alphas[k] = (ub[j] - x[j])/dx[j];
          else if (dx[j] < (Real) 0.0)
            alphas[k] = (lb[j] - x[j])/dx[j];
          if (alphas[k] < (Real) 0.0)
            alphas[k] = (Real) 1.0;
          FILE_LOG(LOG_OPT) << "  -- bounds alpha(" << j << ") = " << alphas[k] << endl;
        }
      }

      // loop over linear inequalities
      for (unsigned j=0, k=N; j< S; j++, k++)
      {
        // constraint cannot be blocking if it is in the working set
        if (!ws.is_constraint_working(j))
        {
          M.get_row(j, workv);
          Real dot = workv.dot(dx);
          if (dot >= (Real) 0.0)
            continue;
          alphas[k] = (q[j] - workv.dot(x))/dot;
          if (alphas[k] < (Real) 0.0)
            alphas[k] = (Real) 0.0;

          FILE_LOG(LOG_OPT) << "  -- alpha(" << j << ") = " << alphas[k] << endl;
          FILE_LOG(LOG_OPT) << "    m[" << j << "]: " << workv << "  q[" << j << "]: " << q[j] << "  m'x: " << workv.dot(x) << "  m'dx: " << workv.dot(dx) << endl;
        }
      }

      // find smallest value of alpha
      vector<Real>::const_iterator alphas_iter = std::min_element(alphas.begin(), alphas.end());
      Real alpha = (alphas_iter != alphas.end()) ? *alphas_iter : (Real) 1.0;

      // determine blocking constraints
      unsigned nblocking = 0;
      if (alpha + std::numeric_limits<Real>::epsilon() < (Real) 1.0)
      {
        for (unsigned i=0; i< alphas.size(); i++)
          if (alphas[i] < alpha + std::numeric_limits<Real>::epsilon())
          {
            blocking[i] = true;
            nblocking++;
          }
      }

      FILE_LOG(LOG_OPT) << "-- alpha: " << alpha << endl;
      FILE_LOG(LOG_OPT) << "  number of blocking constraints: " << nblocking << endl;

      // compute x
      dx *= alpha;
      dx += x;

      // verify that the objective function did not increase 
      G.mult(dx, workv) *= (Real) 0.5;
      workv += c;
      Real eval = dx.dot(workv);
      if (eval > last_eval + qparams.eps)
      {
        FILE_LOG(LOG_OPT) << " -- step increased objective function significantly!  terminating..." << endl;
        FILE_LOG(LOG_OPT) << "  new objective function value: " << eval << endl; 
        FILE_LOG(LOG_OPT) << "  old objective function value: " << last_eval << endl; 
        return;
      }
      last_eval = eval;

      // we can now get x
      x.copy_from(dx); 

      FILE_LOG(LOG_OPT) << "  qp iteration: " << qparams.iterations << endl;
      FILE_LOG(LOG_OPT) << "  new x: " << x << endl;
      FILE_LOG(LOG_OPT) << "  new optimization value: " << (0.5*x.dot(G*x) + x.dot(c)) << endl;

      // see whether we can quit
      if (qparams.tcheck && (*qparams.tcheck)(x, qparams.data))
      {
        FILE_LOG(LOG_OPT) << " -- user specified QP termination" << endl;
        FILE_LOG(LOG_OPT) << "qp_convex_activeset() exiting" << endl;
        return;
      }

      // if there are blocking constraints, add one of them to the working set
      while (nblocking > 0)
      {
        // pick the constraint to add to the working set
        unsigned chosen = rand() % nblocking;

        // get the constraint
        for (unsigned j=0; j< S+N; j++)
          if (blocking[j])
          {
            if (chosen-- == 0)
            {
              chosen = j;
              break;
            }
          }

         // see whether we are adding a bound constraint on a linear constraint
         if (chosen < N)
         {
           if (ws.add_var_to_working_set(chosen))
             break;
           else
           {
             FILE_LOG(LOG_OPT) << "Unable to add variable " << chosen << " to working set!" << endl;
             nblocking--;
             blocking[chosen] = false;
             if (nblocking == 0)
             {
               FILE_LOG(LOG_OPT) << "Unexpectedly unable to add constraint to working set!" << endl;
               return;
             }
           }
         }
         else
         {
           FILE_LOG(LOG_OPT) << "  attempt to add linear inequality constraint " << (chosen-N) << " to the working set" << endl;
           if (ws.add_constraint_to_working_set(chosen-N))
           {
             FILE_LOG(LOG_OPT) << "  constraint successfully added to working set!" << endl;
             break;
           }
           else
           {
             FILE_LOG(LOG_OPT) << "Unable to add constraint " << (chosen-N) << " to working set!" << endl;
             nblocking--;
             blocking[chosen-N] = false;
             if (nblocking == 0)
             {
               FILE_LOG(LOG_OPT) << "Unexpectedly unable to add constraint to working set!" << endl;
               return;
             } 
           }
         }
       }
     }
   
   // now continue looping...
   }
        
  FILE_LOG(LOG_OPT) << "-- qp terminated after maximum number of iterations" << endl; 
}

