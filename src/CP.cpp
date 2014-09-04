/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <limits>
#include <Moby/CompGeom.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/Log.h>
#include <Moby/LCP.h>
#include <Moby/CP.h>

using std::endl;
using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;


// static variables for strictly convex QP solver
Ravelin::MatrixNd CP::_R, CP::_J;
Ravelin::VectorNd CP::_d, CP::_z, CP::_s, CP::_r, CP::_np, CP::_u, CP::_xold, CP::_uold;
std::vector<unsigned> CP::_a, CP::_aold, CP::_iai;
std::vector<bool> CP::_iaexcl;


/// Finds closest point/common point for two shapes 
double CP::find_cpoint(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB, Point3d& closestA, Point3d& closestB)
{
  MatrixNd A, M, R, S; 
  VectorNd b, q, r, s;
  const double INF = 1e8;

  FILE_LOG(LOG_COLDET) << "CP::find_cpoint() entered" << std::endl;

  // the idea behind this is to solve the linear program:
  // minimize t (with variables x and t)
  // s.t. Mx <= q + t is satisfied
  // this will try to put a point furthest into the interior of one of the
  // polyhedra, where R*x <= r defines one polyhedron and S*x <= s defines
  // the other and:
  // M = | R | and q = | r |
  //     | S |         | s |
  // if t is positive, then polyhedra are separated

  // first, setup the facets
  pA->get_facets(poseA, R, r);
  pB->get_facets(poseB, S, s);
  M.resize(pA->num_facets() + pB->num_facets(), 4);
  q.resize(pA->num_facets() + pB->num_facets());
  M.block(0, pA->num_facets(), 0, 3) = R;
  M.block(pA->num_facets(), M.rows(), 0, 3) = S;
  q.segment(0, pA->num_facets()) = r;
  q.segment(pA->num_facets(), M.rows()) = s;

  // now setup t
  SharedVectorNd one_block = M.column(3);
  one_block.set_one();
  one_block.negate();

  // setup the objective function
  VectorNd c(4);
  c[0] = 0.0;  c[1] = 0.0;  c[2] = 0.0;
  c[3] = -1.0;

  // setup x
  VectorNd x(4);

  // setup lower and upper bounds
  const double Q_NORM = q.norm_inf();
  VectorNd l(4), u(4);
  l.set_one() *= -Q_NORM * 10.0;
  u.set_one() *= Q_NORM * 10.0;

  // now solve the LP
  lp_seidel(M, q, c, l, u, x);

  // determine the distance
  double d = x[3];
  FILE_LOG(LOG_COLDET) << " -- signed distance (via LP): " << d << std::endl;

  // if distance is positive, then the polyhedra are separated
  if (d > 0.0)
  {
    // we want to find the closest point for each polyhedron
    // setup the identity matrix (this is passed in to the solver)
    MatrixNd H(6, 6);
    H.set_identity();

    // there is no linear term
    c.set_zero(6);

    // make A and b empty
    A.resize(0, 6);
    b.resize(0);

    // setup new M and b
    M.resize(R.rows()+S.rows(), 6);
    q.resize(R.rows()+S.rows());
    M.block(0, pA->num_facets(), 0, 3) = R;
    M.block(0, pA->num_facets(), 3, 6).set_zero();
    M.block(pA->num_facets(), M.rows(), 3, 6) = S;
    M.block(pA->num_facets(), M.rows(), 0, 3).set_zero();
    q.segment(0, pA->num_facets()) = r;
    q.segment(pA->num_facets(), q.rows()) = s; 

    // strict convex solver expects things to be in the form M*x >= q 
    M.negate();
    q.negate();

    // y is the variable we are solving for
    VectorNd y(6);

    // output information about the QP solve
    FILE_LOG(LOG_COLDET) << "geometries are separated; solving with QP" << std::endl;
    FILE_LOG(LOG_COLDET) << "M: " << std::endl << M;
    FILE_LOG(LOG_COLDET) << "q: " << q << std::endl;

    // TODO: replace this with strict convex solver
/*
    static LCP lcp;
    MatrixNd MM(H.rows()+M.rows(), H.rows()+M.rows());
    VectorNd qq(H.rows()+M.rows());
    MM.block(0, H.rows(), 0, H.rows()) = H;
    SharedMatrixNd MT_block = MM.block(0, H.rows(), H.rows(), MM.columns());
    MatrixNd::transpose(M, MT_block);
    MM.block(0, H.rows(), H.rows(), MM.columns()).negate();
    MM.block(H.rows(), MM.rows(), 0, H.columns()) = M;
    MM.block(H.rows(), MM.rows(), H.columns(), MM.columns()).set_zero();
    qq.segment(0, H.rows()) = c;
    qq.segment(H.rows(), qq.rows()) = q;
    qq.segment(H.rows(), qq.rows()).negate();
    y.resize(12);
    lcp.lcp_lemke_regularized(MM, qq, y);
    FILE_LOG(LOG_COLDET) << "QP solution: " << y << std::endl;
*/
    if (!qp_strict_convex(H, c, A, b, M, q, y, true))
      throw std::runtime_error("Strictly convex solver failed!");

    // setup closest points in global frame
    Vector3d cpA_global(y[0], y[1], y[2], GLOBAL);
    Vector3d cpB_global(y[3], y[4], y[5], GLOBAL);

    // set closest points
    closestA = Pose3d::transform_point(poseA, cpA_global);
    closestB = Pose3d::transform_point(poseB, cpB_global);

    // compute distance between closest points 
    d = (cpA_global - cpB_global).norm(); 
    FILE_LOG(LOG_COLDET) << " -- closest point on A (via QP): " << cpA_global << std::endl;
    FILE_LOG(LOG_COLDET) << " -- closest point on B (via QP): " << cpB_global << std::endl;
    FILE_LOG(LOG_COLDET) << " -- signed distance (via QP): " << d << std::endl;
  }
  else
  { 
    // the closest point is shared
    Vector3d cp_global(x[0], x[1], x[2], GLOBAL);
    closestA = Pose3d::transform_point(poseA, cp_global);
    closestB = Pose3d::transform_point(poseB, cp_global);
  }
  
  FILE_LOG(LOG_COLDET) << "CP::find_cpoint() exited" << std::endl;

  return d;
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
bool CP::lp_seidel(const MatrixNd& A, const VectorNd& b, const VectorNd& c, const VectorNd& l, const VectorNd& u, VectorNd& x)
{
  // get number of rows and columns in A
  unsigned n = A.rows();
  unsigned d = A.columns();

  FILE_LOG(LOG_COLDET) << "CP::lp() entered" << endl;
  FILE_LOG(LOG_COLDET) << "A: " << endl << A;
  FILE_LOG(LOG_COLDET) << "b: " << b << endl;
  FILE_LOG(LOG_COLDET) << "c: " << c << endl;
  FILE_LOG(LOG_COLDET) << "l: " << l << endl;
  FILE_LOG(LOG_COLDET) << "u: " << u << endl;

  // base case d = 1
  if (d == 1)
  {
    FILE_LOG(LOG_COLDET) << "base case, d = 1" << endl;

    double high = u[0]; 
    double low = l[0];

    for (unsigned i=0; i< n; i++)
    {
      if (A(i,0) > std::numeric_limits<double>::epsilon())
        high = std::min(high, b[i]/A(i,0));
      else if (A(i,0) < -std::numeric_limits<double>::epsilon())
        low = std::max(low, b[i]/A(i,0));
      else if (b[i] < -std::numeric_limits<double>::epsilon())
      {
        FILE_LOG(LOG_COLDET) << "infeasible; b is negative and A is zero" << endl;

        return false; 
      }
    }

    // do a check for infeasibility
    if (high < low)
    {
      FILE_LOG(LOG_COLDET) << "infeasible; high (" << high << ") < low (" << low << ")" << endl;

      return false; 
    }
  
    // set x
    x.resize(1);
    x[0] = (c[0] >= 0.0) ? high : low;

    FILE_LOG(LOG_COLDET) << "optimal 1D x=" << x << endl;
    FILE_LOG(LOG_COLDET) << "CP::lp() exited" << endl;

    // otherwise, good return
    return true; 
  }

  // all work variables
  vector<VectorNd> aac;
  vector<VectorNd> aa;
  vector<unsigned> permut;
  vector<double> bbc;
  MatrixNd Aprime;
  VectorNd bb, aak, cprime, lprime, uprime, bprime, workv, f, g;

  // pick a random shuffle for A and b
  permut.resize(n);
  for (unsigned i=0; i< n; i++)
    permut[i] = i;
  std::random_shuffle(permut.begin(), permut.end());

  // setup aa and bb
  aa.resize(n);
  bb.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    A.get_row(permut[i], aa[i]);
    bb[i] = b[permut[i]];
  }

  FILE_LOG(LOG_COLDET) << "A (permuted): " << endl;
  for (unsigned i=0; i< n; i++)
    FILE_LOG(LOG_COLDET) << aa[i] << endl;
  FILE_LOG(LOG_COLDET) << "b (permuted): " << bb << endl;

  // setup optimum vector
  x.resize(d);
  for (unsigned i=0; i< d; i++)
    if (c[i] > 0.0)
        x[i] = u[i];
    else if (c[i] < 0.0) 
        x[i] = l[i];
    else
        x[i] = (std::fabs(l[i]) < std::fabs(u[i])) ? l[i] : u[i];

  FILE_LOG(LOG_COLDET) << "initial x=" << x << endl;

  // process half-space constraints 
  for (unsigned i=0; i< n; i++)
  {
    FILE_LOG(LOG_COLDET) << "-- processing halfspace constraint " << i << endl;

    // if x respects new halfspace constraint, nothing else to do..
    double val = bb[i] - VectorNd::dot(aa[i], x);
    if (val >= -EPS_DOUBLE)
    {    
      FILE_LOG(LOG_COLDET) << "  ** constraint already satisfied!" << endl;
      continue;
    }

    FILE_LOG(LOG_COLDET) << "  -- constraint not satisfied (" << val << "); solving recursively" << endl;

    // search for maximum value in the a vector
    unsigned k = std::numeric_limits<unsigned>::max();
    double maximal = -std::numeric_limits<double>::max();
    for (unsigned j=0; j< d; j++)
      if (std::fabs(aa[i][j]) > maximal && aa[i][j] != 0.0)
      {
        maximal = std::fabs(aa[i][j]);
        k = j;
      }

    FILE_LOG(LOG_COLDET) << "  -- determined k: " << k << endl;

    // look for infeasibility
    if (k == std::numeric_limits<unsigned>::max())
    {
      FILE_LOG(LOG_COLDET) << "  -- k is infeasible; problem is infeasible" << endl;

      return false; 
    }

    // setup useful vector and constant
    aak = aa[i];
    aak /= aa[i][k];
    double bak = bb[i]/aa[i][k];

    FILE_LOG(LOG_COLDET) << "  -- vector a/a(k): " << aak << " " << bak << endl;

    // copy vectors aa and bb
    aac.resize(i);
    bbc.resize(i);
    for (unsigned j=0; j< i; j++)
    {
      aac[j] = aa[j];
      bbc[j] = bb[j];
    }

    // modify copy of vector aa
    for (unsigned j=0; j< i; j++)
    {
      workv = aak;
      workv *= aa[j][k];
      aac[j] -= workv;
      assert(std::fabs(aac[j][k]) < EPS_DOUBLE);
      aac[j] = remove_component(aac[j], k, workv);
    }

    // modify copy of vector bb 
    for (unsigned j=0; j< i; j++)
    {
      bbc[j] -= bak * aa[j][k];
      if (std::isinf(bbc[j]))
        bbc[j] = finitize(bbc[j]);
    }

    // modify copy of c
    assert(std::fabs((c[k] - aak[k]*c[k])) < EPS_DOUBLE);
    workv = aak;
    workv *= c[k];
    workv -= c;
    workv.negate(); 
    remove_component(workv, k, cprime);  

    // generate new lower and upper bounds for variables
    remove_component(l, k, lprime);
    remove_component(u, k, uprime);

    // setup two new constraints 
    f.set_zero(d);
    g.set_zero(d);
    f[k] = 1;
    g[k] = -1;
    assert(std::fabs((f[k] - aak[k])) < EPS_DOUBLE);
    assert(std::fabs((g[k] + aak[k])) < EPS_DOUBLE);
    f = remove_component(f -= aak, k, workv);
    g = remove_component(g += aak, k, workv);
    double bf = u[k] - bak;
    double bg = -l[k] + bak;
    if (std::isinf(bf))
      bf = finitize(bf);
    if (std::isinf(bg))
      bg = finitize(bg);
    aac.push_back(f);
    bbc.push_back(bf);
    aac.push_back(g);
    bbc.push_back(bg);

    // create the Aprime matrix from aac
    Aprime.resize(aac.size(), d-1);  
    for (unsigned j=0; j< aac.size(); j++)
      Aprime.set_row(j, aac[j]);

    // create the bprime vector from bbc
    bprime.resize(bbc.size());
    for (unsigned j=0; j< bbc.size(); j++)
      bprime[j] = bbc[j];

    FILE_LOG(LOG_COLDET) << "  -- A': " << endl << Aprime;
    FILE_LOG(LOG_COLDET) << "  -- b': " << bprime << endl;
    FILE_LOG(LOG_COLDET) << "  -- c': " << cprime << endl;
    FILE_LOG(LOG_COLDET) << "  -- u': " << uprime << endl;
    FILE_LOG(LOG_COLDET) << "  -- l': " << lprime << endl;
    FILE_LOG(LOG_COLDET) << "  -- f: " << f << " " << bf << endl;
    FILE_LOG(LOG_COLDET) << "  -- g: " << g << " " << bg << endl;
    FILE_LOG(LOG_COLDET) << " + solving recursive subproblem" << endl;

    // solve the (d-1)-dimensional problem and ``lift'' the solution
    if (!lp_seidel(Aprime,bprime,cprime,lprime,uprime,x))
      return false;

    FILE_LOG(LOG_COLDET) << "  -- recursively determined x: " << x << endl;
    FILE_LOG(LOG_COLDET) << "  -- k: " << k << endl;

    // insert a zero into the k'th dimension of x
    x = insert_component(x, k, workv);
    FILE_LOG(LOG_COLDET) << "  -- x w/inserted component at k: " << x << endl;

    // solve for the proper k'th value of x
    x[k] = (bb[i] - VectorNd::dot(aa[i], x))/aa[i][k];
    FILE_LOG(LOG_COLDET) << "  -- optimal x (to this point): " << x << endl;
/*
    // verify that half-plane constraints still met
    for (unsigned j=0; j<= i; j++)
      assert(VectorNd::dot(aa[j], x) - bb[j] <= EPS_DOUBLE);

    // verify that lower constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] >= l[j] - EPS_DOUBLE);

    // verify that upper constraints still met
    for (unsigned j=0; j< l.size(); j++)
      assert(x[j] <= u[j] + EPS_DOUBLE);
*/
  }

  // verify that half-plane constraints still met
  if (LOGGING(LOG_COLDET))
  {
    for (unsigned i=0; i< b.rows(); i++)
      FILE_LOG(LOG_COLDET) << i << ": b - A*x = " << (b[i] - VectorNd::dot(A.row(i), x)) << endl;
  }
/*
  for (unsigned j=0; j < n; j++)
    assert(VectorNd::dot(aa[j], x) - bb[j] <= EPS_DOUBLE);
  // verify that lower constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] >= l[i] - EPS_DOUBLE);

  // verify that upper constraints still met
  for (unsigned i=0; i< l.size(); i++)
    assert(x[i] <= u[i] + EPS_DOUBLE);
*/

  FILE_LOG(LOG_COLDET) << "all halfspace constraints satisfied; optimum found!" << endl;
  FILE_LOG(LOG_COLDET) << "optimum = " << x << endl;
  FILE_LOG(LOG_COLDET) << "CP::lp_seidel() exited" << endl;

  return true; 
}

/// Inserts a component (value will be zero) at the k'th position in the given vector and returns a new vector
VectorNd& CP::insert_component(const VectorNd& x, unsigned k, VectorNd& xn)
{
  xn.resize(x.size()+1);

  if (k == 0)
    xn.set_sub_vec(1, x);
  else if (k == x.size())
    xn.set_sub_vec(0, x);
  else
  {
    xn.set_sub_vec(0,x.segment(0,k));
    xn.set_sub_vec(k+1,x.segment(k,x.size()));
  }  
  xn[k] = 0;
  
  return xn;
}

/// Gets a subvector with the k'th component removed
/**
 * \note for use by lp()
 */
VectorNd& CP::remove_component(const VectorNd& v, unsigned k, VectorNd& vp)
{
  const unsigned d = v.size();

  if (k == 0)
    return v.get_sub_vec(1,d,vp);
  else if (k == d-1)
    return v.get_sub_vec(0,d-1,vp);
  else
  {
    // setup v w/component k removed
    vp.resize(d-1);    
    vp.set_sub_vec(0, v.segment(0, k));
    vp.set_sub_vec(k, v.segment(k+1, d)); 
    return vp;
  }
}

/// Makes a (possibly) infinite value finite again
/**
 * Converts values of -inf to -DBL_MAX and inf to DBL_MAX
 * \note utility function for lp()
 */
double CP::finitize(double x)
{
  if (x == std::numeric_limits<double>::infinity())
    return std::numeric_limits<double>::max();
  else if (x == -std::numeric_limits<double>::infinity())
    return -std::numeric_limits<double>::max();
  else
    return x;
}

// The Solving function, implementing the Goldfarb-Idnani method
bool CP::qp_strict_convex(MatrixNd& G, const VectorNd& c, 
                      const MatrixNd& A, const VectorNd& b,  
                      const MatrixNd& M, const VectorNd& q, 
                      VectorNd& x, bool G_factored)
{
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());
  const double INF = std::numeric_limits<double>::max();
  unsigned n = G.columns();
  unsigned p = A.rows();
  unsigned m = M.rows();

  // check matrix sizes
  if (G.rows() != n)
    throw NonsquareMatrixException();
  if (A.columns() != n)
    throw MissizeException();
  if (b.rows() != p)
    throw MissizeException();
  if (M.columns() != n)
    throw MissizeException();
  if (q.size() != m)
    throw MissizeException();

  // make x the correct size
  x.resize(n);
  double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1 
    * and the full step length t2 */

  // initialize the iteration variable
  unsigned iter = 0;

  // size vectors
  _a.resize(m+p);
  _aold.resize(m+p);
  _iai.resize(m+p);
  _s.resize(m+p);
  _r.resize(m+p);
  _np.resize(n);
  _u.resize(m+p);
  _xold.resize(n);
  _uold.resize(m+p);
  _iaexcl.resize(m + p);
  
  /*
   * Preprocessing phase
   */
  
  /* compute the trace of G */
  double c1 = 0.0;
  RowIteratord Giter = G.row_iterator_begin(); 
  for (unsigned i = 0; i < n; i++, Giter+= n)
    c1 += *Giter;

  // if G has been factored, square the trace
  if (G_factored)
    c1 *= c1;

  // do Cholesky factorization if necessary
  if (!G_factored)
  {
    if (!LinAlgd::factor_chol(G))
      return false;
  }
 
  /* initialize the matrix R */
  _R.set_zero(n, n);
  _d.set_zero(n);
  double R_norm = 1.0; /* this variable will hold the norm of the matrix R */
  
  /* compute the inverse of the factorized matrix G^-1, this is the initial value for _J */
  _J = G;
  LinAlgd::inverse_chol(_J);
  _z = _J.column(_J.columns()-1);

  // compute the trace of _J
  double c2 = 0.0;
  RowIteratord Jiter = _J.row_iterator_begin(); 
  for (unsigned i = 0; i < n; i++, Jiter+= n)
    c2 += *Jiter;
  
  /* c1 * c2 is an estimate for cond(G) */
  
  /* 
    * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + c x 
   * this is a feasible point in the dual space
   * x = G^-1 * c
   */
/*
  x = c;
  LinAlgd::solve_chol_fast(G, x);
  x.negate();
*/
  /* and compute the current solution value */ 
  double f_value = 0.5 * c.dot(x);
  
  /* Add equality constraints to the working set A */
  unsigned iq = 0;
  for (unsigned i = 0; i < p; i++)
  {
    _np = A.row(i);
    _J.transpose_mult(_np, _d);
    update_z(_J, _d, iq, _z);
    update_r(_R, _d, iq, _r);
    
    /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
      becomes feasible */
    double t2 = 0.0;
    if (std::fabs(_z.dot(_z)) > NEAR_ZERO) // i.e. z != 0
      t2 = (-_np.dot(x) + b[i]) / _z.dot(_np);
    
    /* set x = x + t2 * z */
    for (unsigned k = 0; k < n; k++)
      x[k] += t2 * _z[k];
    
    /* set u = u+ */
    _u[iq] = t2;
    for (unsigned k = 0; k < iq; k++)
      _u[k] -= t2 * _r[k];
    
    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * _z.dot(_np);
    _a[i] = -i - 1;
    
    if (!add_constraint(_R, _J, _d, iq, R_norm))
    {    
      // Equality constraints are linearly dependent
      throw std::runtime_error("Constraints are linearly dependent");
      return false;
    }
  }
  
  /* set _iai = K \ A */
  for (unsigned i = 0; i < m; i++)
    _iai[i] = i;
  
l1:  iter++;
  /* step 1: choose a violated constraint */
  unsigned ip = std::numeric_limits<unsigned>::max();
  for (unsigned i = p; i < iq; i++)
  {
    ip = _a[i];
    _iai[ip] = -1;
  }
  
  /* compute s[x] = M * x - q for all elements of K \ A */
  double ss = 0.0;
  double psi = 0.0; /* this value will contain the sum of all infeasibilities */
  ip = 0; /* ip will be the index of the chosen violated constraint */
  for (unsigned i = 0; i < m; i++)
  {
    ColumnIteratord_const Miter = M.row(i).column_iterator_begin();
    _iaexcl[i] = true;
    double sum = 0.0;
    for (unsigned j = 0; j < n; j++)
      sum += Miter[j] * x[j];
    sum -= q[i];
    _s[i] = sum;
    psi += std::min(0.0, sum);
  }
  
  if (std::fabs(psi) <= m * std::numeric_limits<double>::epsilon() * c1 * c2* 100.0)
  {
    /* numerically there are not infeasibilities anymore */
    return true;
  }
  
  /* save old values for u and A */
  for (unsigned i = 0; i < iq; i++)
  {
    _uold[i] = _u[i];
    _aold[i] = _a[i];
  }
  /* and for x */
  _xold = x;
  
l2: /* Step 2: check for feasibility and determine a new S-pair */
    for (unsigned i = 0; i < m; i++)
    {
      if (_s[i] < ss && _iai[i] < std::numeric_limits<unsigned>::max() && _iaexcl[i])
      {
        ss = _s[i];
        ip = i;
      }
    }
  if (ss >= 0.0)
    return true;
  
  /* set _np = n[ip] */
  _np = M.row(ip);
  /* set u = [u 0]^T */
  _u[iq] = 0.0;
  /* add ip to the active set A */
  _a[iq] = ip;
  
 
l2a:/* Step 2a: determine step direction */
    /* compute z = H _np: the step direction in the primal space (through J, see the paper) */
   _J.transpose_mult(_np, _d);
  update_z(_J, _d, iq, _z);
  /* compute N* _np (if q > 0): the negative of the step direction in the dual space */
  update_r(_R, _d, iq, _r);
  
  /* Step 2b: compute step length */
  unsigned l = 0;
  /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
  t1 = INF; /* +INF */
  /* find the index l s.t. it reaches the minimum of u+[x] / r */
  for (unsigned k = p; k < iq; k++)
  {
    if (_r[k] > 0.0)
    {
      if (_u[k] / _r[k] < t1)
      {
        t1 = _u[k] / _r[k];
        l = _a[k];
      }
    }
  }
  /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
  if (std::fabs(_z.dot(_z))  > NEAR_ZERO) // i.e. z != 0
    t2 = -_s[ip] / _z.dot(_np);
  else
    t2 = INF; /* +INF */
  
  /* the step is chosen as the minimum of t1 and t2 */
  t = std::min(t1, t2);
  
  /* Step 2c: determine new S-pair and take step: */
  
  /* case (i): no step in primal or dual space */
  if (t >= INF)
    return false;
  /* case (ii): step in dual space */
  if (t2 >= INF)
  {
    /* set u = u +  t * [-r 1] and drop constraint l from the active set _a */
    for (unsigned k = 0; k < iq; k++)
      _u[k] -= t * _r[k];
    _u[iq] += t;
    _iai[l] = l;
    delete_constraint(_R, _J, _a, _u, n, p, iq, l);
    goto l2a;
  }
  
  /* case (iii): step in primal and dual space */
  
  /* set x = x + t * z */
  for (unsigned k = 0; k < n; k++)
    x[k] += t * _z[k];
  /* update the solution value */
  f_value += t * _z.dot(_np) * (0.5 * t + _u[iq]);
  /* u = u + t * [-r 1] */
  for (unsigned k = 0; k < iq; k++)
    _u[k] -= t * _r[k];
  _u[iq] += t;

  if (std::fabs(t - t2) < NEAR_ZERO)
  {
    /* full step has taken */
    /* add constraint ip to the active set*/
    if (!add_constraint(_R, _J, _d, iq, R_norm))
    {
      _iaexcl[ip] = false;
      delete_constraint(_R, _J, _a, _u, n, p, iq, ip);
      for (unsigned i = 0; i < m; i++)
        _iai[i] = i;
      for (unsigned i = p; i < iq; i++)
      {
        _a[i] = _aold[i];
        _u[i] = _uold[i];
        _iai[_a[i]] = -1;
      }
      for (unsigned i = 0; i < n; i++)
        x[i] = _xold[i];
      goto l2; /* go to step 2 */
    }    
    else
      _iai[ip] = std::numeric_limits<unsigned>::max();
    goto l1;
  }
  
  /* a patial step has taken */
  /* drop constraint l */
  _iai[l] = l;
  delete_constraint(_R, _J, _a, _u, n, p, iq, l);

  /* update s[ip] = M * x - q */
  double sum = 0.0;
  for (unsigned k = 0; k < n; k++)
    sum += M.row(ip).dot(x);
  _s[ip] = sum - q[ip];
  
  goto l2a;
}

void CP::update_z(const MatrixNd& J, const VectorNd& d, unsigned iq, VectorNd& z)
{
  const unsigned N = z.size();
  
  /* setting of z = H * d */
  for (unsigned i = 0; i < N; i++)
  {
    z[i] = 0.0;
    RowIteratord_const Jiter = J.row(i).row_iterator_begin()+iq;
    for (unsigned j = iq; j < N; j++, Jiter++)
      z[i] += *Jiter * d[j];
  }
}

void CP::update_r(const MatrixNd& R, const VectorNd& d, unsigned iq, VectorNd& r)
{
  const unsigned N = d.size();

  /* setting of r = R^-1 d */
  for (unsigned i = iq - 1; i < N; i--)
  {
    double sum = 0.0;
    RowIteratord_const Riter = R.row(i).row_iterator_begin()+i+1;
    for (unsigned j = i + 1; j < iq; j++, Riter++)
      sum += *Riter * r[j];
    r[i] = (d[i] - sum) / R(i,i);
  }
}

bool CP::add_constraint(MatrixNd& R, MatrixNd& J, VectorNd& d, unsigned& iq, double& R_norm)
{
  const unsigned N = d.size();
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());  
  
  /* we have to find the Givens rotation which will reduce the element
    d[j] to zero.
    if it is already zero we don't have to do anything, except of
    decreasing j */  
  for (unsigned j = N - 1; j >= iq + 1 && j < N; j--)
  {
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
    If cc is one, then element (j) of d is zero compared with element
    (j - 1). Hence we don't have to do anything. 
    If cc is zero, then we just have to switch column (j) and column (j - 1) 
    of J. Since we only switch columns in J, we have to be careful how we
    update d depending on the sign of gs.
    Otherwise we have to apply the Givens rotation to these columns.
    The i - 1 element of d has to be updated to h. */
    double cc = d[j - 1];
    double ss = d[j];
    double h = distance(cc, ss);
    if (std::fabs(h) < NEAR_ZERO) // h == 0
      continue;
    d[j] = 0.0;
    ss = ss / h;
    cc = cc / h;
    if (cc < 0.0)
    {
      cc = -cc;
      ss = -ss;
      d[j - 1] = -h;
    }
    else
      d[j - 1] = h;
    double xny = ss / (1.0 + cc);
    RowIteratord Jjm1iter = J.column(j-1).row_iterator_begin();
    RowIteratord Jjiter = J.column(j).row_iterator_begin();
    for (unsigned k = 0; k < N; k++)
    {
      double t1 = Jjm1iter[k];
      double t2 = Jjiter[k];
      Jjm1iter[k] = t1 * cc + t2 * ss;
      Jjiter[k] = xny * (t1 + Jjm1iter[k]) - t2;
    }
  }
  /* update the number of constraints added*/
  iq++;
  /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
  RowIteratord Riter = R.column(iq-1).row_iterator_begin();
  for (unsigned i = 0; i < iq; i++)
    Riter[i] = d[i];
  
  if (std::fabs(d[iq - 1]) <= std::numeric_limits<double>::epsilon() * R_norm) 
  {
    // problem degenerate
    return false;
  }
  R_norm = std::max(R_norm, std::fabs(d[iq - 1]));
  return true;
}

void CP::delete_constraint(MatrixNd& R, MatrixNd& J, std::vector<unsigned>& A, VectorNd& u, unsigned n, unsigned p, unsigned& iq, unsigned l)
{
  unsigned qq = std::numeric_limits<unsigned>::max(); 
  double cc, ss, h, xny, t1, t2;
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());  

  /* Find the index qq for active constraint l to be removed */
  for (unsigned i = p; i < iq; i++)
    if (A[i] == l)
    {
      qq = i;
      break;
    }
      
  /* remove the constraint from the active set and the duals */
  for (unsigned i = qq; i < iq - 1; i++)
  {
      A[i] = A[i + 1];
      u[i] = u[i + 1];
      R.column(i) = R.column(i+1);
    }
      
  A[iq - 1] = A[iq];
  u[iq - 1] = u[iq];
  A[iq] = 0; 
  u[iq] = 0.0;
  R.column(iq-1).set_zero();
  /* constraint has been fully removed */
  iq--;
  
  if (iq == 0)
    return;
  
  for (unsigned j = qq; j < iq; j++)
  {
    cc = R(j,j);
    ss = R(j + 1,j);
    h = distance(cc, ss);
    if (std::fabs(h) < NEAR_ZERO) // h == 0
      continue;
    cc = cc / h;
    ss = ss / h;
    R(j+1,j) = 0.0;
    if (cc < 0.0)
    {
      R(j,j) = -h;
      cc = -cc;
      ss = -ss;
    }
    else
      R(j,j) = h;
    
    xny = ss / (1.0 + cc);
    ColumnIteratord Riterj = R.row(j).column_iterator_begin();
    ColumnIteratord Riterjp1 = R.row(j+1).column_iterator_begin();
    for (unsigned k = j + 1; k < iq; k++)
    {
      t1 = Riterj[k];
      t2 = Riterjp1[k];
      Riterj[k] = t1 * cc + t2 * ss;
      Riterjp1[k] = xny * (t1 + R(j,k)) - t2;
    }

    RowIteratord Jiterj = J.column(j).row_iterator_begin();
    RowIteratord Jiterjp1 = J.column(j+1).row_iterator_begin();
    for (unsigned k = 0; k < n; k++)
    {
      t1 = Jiterj[k];
      t2 = Jiterjp1[k];
      Jiterj[k] = t1 * cc + t2 * ss;
      Jiterjp1[k] = xny * (Jiterj[k] + t1) - t2;
    }
  }
}

double CP::distance(double a, double b)
{
  double a1 = std::fabs(a);
  double b1 = std::fabs(b);
  if (a1 > b1) 
  {
    double t = (b1 / a1);
    return a1 * std::sqrt(1.0 + t * t);
  }
  else
    if (b1 > a1)
    {
      double t = (a1 / b1);
      return b1 * std::sqrt(1.0 + t * t);
    }
  return a1 * std::sqrt(2.0);
}


