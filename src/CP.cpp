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


/// Finds closest point/common point for two shapes 
double CP::find_cpoint(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB, Point3d& closestA, Point3d& closestB)
{
  MatrixNd A, M, R, S; 
  VectorNd b, q, r, s;
  const double INF = 1e8;
  static VectorNd l, u, x;

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

  // setup lower and upper bounds
  const double Q_NORM = q.norm_inf();
  l.resize(4);
  u.resize(4);
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
    // setup linear term 
    c.set_zero(8);
    c[6] = c[7] = -1.0;

    // setup new M and b
    M.resize(R.rows()+S.rows(), 8);
    q.resize(R.rows()+S.rows());
    M.block(0, pA->num_facets(), 0, 3) = R;
    M.block(0, pA->num_facets(), 3, 6).set_zero();
    M.block(pA->num_facets(), M.rows(), 3, 6) = S;
    M.block(pA->num_facets(), M.rows(), 0, 3).set_zero();
    q.segment(0, pA->num_facets()) = r;
    q.segment(pA->num_facets(), q.rows()) = s; 

    // setup s and t components of M
    M.column(6).segment(0, pA->num_facets()).set_one();
    M.column(6).segment(pA->num_facets(), M.rows()).set_zero();
    M.column(7).segment(0, pA->num_facets()).set_zero();
    M.column(7).segment(pA->num_facets(), M.rows()).set_one();

    // setup l and u again
    l.resize(8);
    u.resize(8);
    l.set_one() *= -Q_NORM * 10.0;
    u.set_one() *= Q_NORM * 10.0;

    // resize x 
    x.resize(8);

    // output information about the LP solve
    FILE_LOG(LOG_COLDET) << "geometries are separated; solving with LP" << std::endl;
    FILE_LOG(LOG_COLDET) << "M: " << std::endl << M;
    FILE_LOG(LOG_COLDET) << "q: " << q << std::endl;

    // do the LP solve
    lp_seidel(M, q, c, l, u, x);
    FILE_LOG(LOG_COLDET) << "LP solution: " << x << std::endl;

    // setup closest points in global frame
    Vector3d cpA_global(x[0], x[1], x[2], GLOBAL);
    Vector3d cpB_global(x[3], x[4], x[5], GLOBAL);

    // set closest points
    closestA = Pose3d::transform_point(poseA, cpA_global);
    closestB = Pose3d::transform_point(poseB, cpB_global);

    // compute distance between closest points 
    d = (cpA_global - cpB_global).norm(); 
    FILE_LOG(LOG_COLDET) << " -- closest point on A (via LP): " << cpA_global << std::endl;
    FILE_LOG(LOG_COLDET) << " -- closest point on B (via LP): " << cpB_global << std::endl;
    FILE_LOG(LOG_COLDET) << " -- signed distance (via LP): " << d << std::endl;
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


