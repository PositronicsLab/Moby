/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <limits>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/LP.h>

using std::endl;
using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

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
bool LP::lp_seidel(const MatrixNd& A, const VectorNd& b, const VectorNd& c, const VectorNd& l, const VectorNd& u, VectorNd& x)
{
  // get number of rows and columns in A
  unsigned n = A.rows();
  unsigned d = A.columns();

  FILE_LOG(LOG_OPT) << "LP::lp() entered" << endl;
  FILE_LOG(LOG_OPT) << "A: " << endl << A;
  FILE_LOG(LOG_OPT) << "b: " << b << endl;
  FILE_LOG(LOG_OPT) << "c: " << c << endl;
  FILE_LOG(LOG_OPT) << "l: " << l << endl;
  FILE_LOG(LOG_OPT) << "u: " << u << endl;

  // base case d = 1
  if (d == 1)
  {
    FILE_LOG(LOG_OPT) << "base case, d = 1" << endl;

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
    x.resize(1);
    x[0] = (c[0] >= 0.0) ? high : low;

    FILE_LOG(LOG_OPT) << "optimal 1D x=" << x << endl;
    FILE_LOG(LOG_OPT) << "LP::lp() exited" << endl;

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
    double val = bb[i] - VectorNd::dot(aa[i], x);
    if (val >= -EPS_DOUBLE)
    {    
      FILE_LOG(LOG_OPT) << "  ** constraint already satisfied!" << endl;
      continue;
    }

    FILE_LOG(LOG_OPT) << "  -- constraint not satisfied (" << val << "); solving recursively" << endl;

    // search for maximum value in the a vector
    unsigned k = std::numeric_limits<unsigned>::max();
    double maximal = -std::numeric_limits<double>::max();
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
    aak = aa[i];
    aak /= aa[i][k];
    double bak = bb[i]/aa[i][k];

    FILE_LOG(LOG_OPT) << "  -- vector a/a(k): " << aak << " " << bak << endl;

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

    FILE_LOG(LOG_OPT) << "  -- A': " << endl << Aprime;
    FILE_LOG(LOG_OPT) << "  -- b': " << bprime << endl;
    FILE_LOG(LOG_OPT) << "  -- c': " << cprime << endl;
    FILE_LOG(LOG_OPT) << "  -- u': " << uprime << endl;
    FILE_LOG(LOG_OPT) << "  -- l': " << lprime << endl;
    FILE_LOG(LOG_OPT) << "  -- f: " << f << " " << bf << endl;
    FILE_LOG(LOG_OPT) << "  -- g: " << g << " " << bg << endl;
    FILE_LOG(LOG_OPT) << " + solving recursive subproblem" << endl;

    // solve the (d-1)-dimensional problem and ``lift'' the solution
    if (!lp_seidel(Aprime,bprime,cprime,lprime,uprime,x))
      return false;

    FILE_LOG(LOG_OPT) << "  -- recursively determined x: " << x << endl;
    FILE_LOG(LOG_OPT) << "  -- k: " << k << endl;

    // insert a zero into the k'th dimension of x
    x = insert_component(x, k, workv);
    FILE_LOG(LOG_OPT) << "  -- x w/inserted component at k: " << x << endl;

    // solve for the proper k'th value of x
    x[k] = (bb[i] - VectorNd::dot(aa[i], x))/aa[i][k];
    FILE_LOG(LOG_OPT) << "  -- optimal x (to this point): " << x << endl;
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
  if (LOGGING(LOG_OPT))
  {
    for (unsigned i=0; i< b.rows(); i++)
      FILE_LOG(LOG_OPT) << i << ": b - A*x = " << (b[i] - VectorNd::dot(A.row(i), x)) << endl;
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

  FILE_LOG(LOG_OPT) << "all halfspace constraints satisfied; optimum found!" << endl;
  FILE_LOG(LOG_OPT) << "optimum = " << x << endl;
  FILE_LOG(LOG_OPT) << "LP::lp_seidel() exited" << endl;

  return true; 
}

/// Inserts a component (value will be zero) at the k'th position in the given vector and returns a new vector
VectorNd& LP::insert_component(const VectorNd& x, unsigned k, VectorNd& xn)
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
VectorNd& LP::remove_component(const VectorNd& v, unsigned k, VectorNd& vp)
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
double LP::finitize(double x)
{
  if (x == std::numeric_limits<double>::infinity())
    return std::numeric_limits<double>::max();
  else if (x == -std::numeric_limits<double>::infinity())
    return -std::numeric_limits<double>::max();
  else
    return x;
}


