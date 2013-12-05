/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <numeric>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <iostream>
#include <boost/algorithm/minmax.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/foreach.hpp>
#include <Ravelin/select>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/SingularException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/NumericalException.h>
#include <Moby/Log.h>
#include <Moby/Constants.h>
#include <Moby/LCP.h>

using namespace Ravelin;
using namespace Moby;
using std::endl;
using std::vector;
using std::pair;
using std::map;
using std::make_pair;
using boost::shared_ptr;

// Sole constructor
LCP::LCP()
{
}

/// Fast pivoting algorithm for denerate, monotone LCPs with few nonzero, nonbasic variables 
bool LCP::lcp_fast(const MatrixNd& M, const VectorNd& q, VectorNd& z, double zero_tol)
{
  const unsigned N = q.rows();

  // look for trivial solution
  if (N == 0)
  {
    z.set_zero(0);
    return true;
  }

  // set zero tolerance if necessary
  if (zero_tol < 0.0)
    zero_tol = M.rows() * M.norm_inf() * std::numeric_limits<double>::epsilon();

  // get minimum element of q (really w)
  unsigned minw = std::min_element(q.begin(), q.end()) - q.begin();
  if (q[minw] > -zero_tol)
  {
    z.set_zero(N);
    return true;
  }

  // setup basic and nonbasic variable indices
  _nonbas.clear();
  _nonbas.push_back(minw); 
  _bas.clear();
  _bas.resize(N-1);
  for (unsigned i=0, j=0; i< N; i++)
    if (i != minw)
      _bas[j++] = i;

  // loop for maximum number of pivots
  const unsigned MAX_PIV = std::max(N*N, (unsigned) 1000);
  for (unsigned piv=0; piv < MAX_PIV; piv++)
  {
    // select nonbasic indices
    M.select_square(_nonbas.begin(), _nonbas.end(), _Msub);
    M.select_columns(_nonbas.begin(), _nonbas.end(), _Mcolsub);
    q.select(_nonbas.begin(), _nonbas.end(), _z);
    _z.negate();

    // solve for nonbasic z
    _LA.solve_fast(_Msub, _z);

    // compute w and find minimum value
    _Mcolsub.mult(_z, _w) += q;
    minw = rand_min(_w, zero_tol);

    // if w >= 0, check whether any component of z < 0
    if (minw > -zero_tol)
    {
      // find the (a) minimum of z
      unsigned minz = rand_min(_z, zero_tol); 
      if (_z[minz] < -zero_tol)
      {
        // get the original index and remove it from the nonbasic set
        unsigned idx = _nonbas[minz];
        _nonbas[minz] = _nonbas.back();
        _nonbas.pop_back();

        // move index to basic set and continue looping
        _bas.push_back(idx);
      }
      else
      {
        // found the solution
        z.set_zero(N);

        // set values of z corresponding to _z
        for (unsigned i=0, j=0; j < _bas.size(); i++, j++)
          z[_bas[j]] = _z[i];

        return true;
      }
    }
    else
    {
      // one or more components of w violating w >= 0
      // move component of w from basic set to nonbasic set
      unsigned idx = _bas[minw];
      _bas[minw] = _bas.back();
      _bas.pop_back();
      _nonbas.push_back(idx);

      // look whether any component of z needs to move to basic set
      unsigned minz = rand_min(_z, zero_tol); 
      if (_z[minz] < -zero_tol)
      {
        // move index to basic set and continue looping
        unsigned idx = _nonbas[minz];
        _nonbas[minz] = _nonbas.back();
        _nonbas.pop_back();
        _bas.push_back(idx);
      }
    }
  }

  // if we're here, then the maximum number of pivots has been exceeded
  return false;
}

/// Get the minimum index of vector v; if there are multiple minima (within zero_tol), returns one randomly 
unsigned LCP::rand_min(const VectorNd& v, double zero_tol)
{
  static vector<unsigned> minima;
  minima.clear();
  unsigned minv = std::min_element(v.begin(), v.end()) - v.begin();
  for (unsigned i=0; i< v.rows(); i++)
    if (v[i] < v[minv] + zero_tol)
      minima.push_back(i);
  return minima[rand() % minima.size()];
}

/// Regularized wrapper around Lemke's algorithm
bool LCP::lcp_lemke_regularized(const MatrixNd& M, const VectorNd& q, VectorNd& z, int min_exp, unsigned step_exp, int max_exp, double piv_tol, double zero_tol)
{
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() entered" << endl;
  pair<ColumnIteratord, ColumnIteratord> mmax;
 
  // look for fast exit
  if (q.size() == 0)
  {
    z.resize(0);
    return true;
  }

  // copy MM
  _MM = M;

  // assign value for zero tolerance, if necessary
  const double ZERO_TOL = (zero_tol > (double) 0.0) ? zero_tol : q.size() * M.norm_inf() * NEAR_ZERO;

  FILE_LOG(LOG_OPT) << " zero tolerance: " << ZERO_TOL << endl;

  // try non-regularized version first
  bool result = lcp_lemke(_MM, q, z, piv_tol, zero_tol);
  if (result)
  {
    // verify that solution truly is a solution -- check z
    if (*std::min_element(z.begin(), z.end()) >= -ZERO_TOL)
    {
      // check w
      M.mult(z, _wx) += q;
      if (*std::min_element(_wx.begin(), _wx.end()) >= -ZERO_TOL)
      {
        // check z'w
        std::transform(z.begin(), z.end(), _wx.begin(), _wx.begin(), std::multiplies<double>());
        mmax = boost::minmax_element(_wx.begin(), _wx.end());
        if (*mmax.first >= -ZERO_TOL && *mmax.second < ZERO_TOL)
        {
          FILE_LOG(LOG_OPT) << "  solved with no regularization necessary!" << endl;
          FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;
          return true;
        }
        else
        {
          FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() - '<w, z> not within tolerance(min value: " << *mmax.first << " max value: " << *mmax.second << ")" << std::endl; 
        }
      }
      else
      {
        FILE_LOG(LOG_OPT) << "  LCP::lcp_lemke() - 'w' not solved to desired tolerance" << std::endl;
        FILE_LOG(LOG_OPT) << "  minimum w: " << *std::min_element(_wx.column_iterator_begin(), _wx.column_iterator_end()) << std::endl;
      }
    }
    else
    {
      FILE_LOG(LOG_OPT) << "  LCP::lcp_lemke() - 'z' not solved to desired tolerance" << std::endl;
      FILE_LOG(LOG_OPT) << "  minimum z: " << *std::min_element(z.column_iterator_begin(), z.column_iterator_end()) << std::endl;
    }
  }

  // start the regularization process
  int rf = min_exp;
  while (rf < max_exp)
  {
    // setup regularization factor
    double lambda = std::pow((double) 10.0, (double) rf);

    // regularize M
    _MM = M;
    for (unsigned i=0; i< M.rows(); i++)
      _MM(i,i) += lambda;

    // try to solve the LCP
    if ((result = lcp_lemke(_MM, q, z, piv_tol, zero_tol)))
    {
      // verify that solution truly is a solution -- check z
      if (*std::min_element(z.begin(), z.end()) > -ZERO_TOL)
      {
        // check w
        _MM.mult(z, _wx) += q;
        if (*std::min_element(_wx.begin(), _wx.end()) > -ZERO_TOL)
        {
          // check z'w
          std::transform(z.begin(), z.end(), _wx.begin(), _wx.begin(), std::multiplies<double>());
          mmax = boost::minmax_element(_wx.begin(), _wx.end());
          if (*mmax.first > -ZERO_TOL && *mmax.second < ZERO_TOL)
          {
            FILE_LOG(LOG_OPT) << "  solved with regularization factor: " << lambda << endl;
            FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;

            return true;
          }
          else
          {
            FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() - '<w, z> not within tolerance(min value: " << *mmax.first << " max value: " << *mmax.second << ")" << std::endl; 
          }
        }
        else
        {
          FILE_LOG(LOG_OPT) << "  LCP::lcp_lemke() - 'w' not solved to desired tolerance" << std::endl;
          FILE_LOG(LOG_OPT) << "  minimum w: " << *std::min_element(_wx.column_iterator_begin(), _wx.column_iterator_end()) << std::endl;
        }
      }
      else
      {
        FILE_LOG(LOG_OPT) << "  LCP::lcp_lemke() - 'z' not solved to desired tolerance" << std::endl;
        FILE_LOG(LOG_OPT) << "  minimum z: " << *std::min_element(z.column_iterator_begin(), z.column_iterator_end()) << std::endl;
      }
    }

    // increase rf
    rf += step_exp;
  }

  FILE_LOG(LOG_OPT) << "  unable to solve given any regularization!" << endl;
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;

  // still here?  failure...
  return false;
}

/// Sets a basis
void LCP::set_basis(unsigned n, unsigned count, vector<unsigned>& bas, vector<unsigned>& nbas)
{
  // clear bas and nbas
  bas.clear();
  nbas.clear();

  unsigned long long n2 = 1 << (n-1);
  for (unsigned i=0; i< n; i++)
  {
    if (count / n2 > 0)
      bas.push_back(i);
    else
      nbas.push_back(i);
    count = count % n2;
    n2 = n2 >> 1;
  } 
}

/// Lemke's algorithm for solving linear complementarity problems
/**
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
bool LCP::lcp_lemke(const MatrixNd& M, const VectorNd& q, VectorNd& z, double piv_tol, double zero_tol)
{
  const unsigned n = q.size();
  const unsigned MAXITER = std::min((unsigned) 1000, 50*n);

  // look for immediate exit
  if (n == 0)
  {
    z.resize(0);
    return true;
  }

  // clear all vectors
  _all.clear();
  _tlist.clear();
  _bas.clear();
  _nonbas.clear();
  _j.clear();

  // copy z to z0
z.set_zero(n);
  _z0 = z;

  // come up with a sensible value for zero tolerance if none is given
  if (zero_tol <= (double) 0.0)
    zero_tol = std::numeric_limits<double>::epsilon() * M.norm_inf() * n;

  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() entered" << endl;
  FILE_LOG(LOG_OPT) << "  M: " << endl << M;
  FILE_LOG(LOG_OPT) << "  q: " << q << endl;

  // see whether trivial solution exists
  if (*std::min_element(q.begin(), q.end()) > -zero_tol)
  {
    FILE_LOG(LOG_OPT) << " -- trivial solution found" << endl;
    FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;
    z.set_zero(n);
    return true;
  }

  // initialize variables
  z.set_zero(n*2);
  unsigned t = 2*n;
  unsigned entering = t;
  unsigned leaving = 0;
  _all.clear();
  for (unsigned i=0; i< n; i++)
    _all.push_back(i);
  unsigned lvindex;
  unsigned idx;
  vector<unsigned>::iterator iiter;
  _tlist.clear();

  // determine initial basis
  _bas.clear();
  _nonbas.clear();
  if (_z0.size() != n)
    for (unsigned i=0; i< n; i++)
      _nonbas.push_back(i);
  else
    for (unsigned i=0; i< n; i++)
      if (_z0[i] > 0)
        _bas.push_back(i);
      else
        _nonbas.push_back(i);

  // B should ideally be a sparse matrix
  _Bl.set_identity(n);
  _Bl.negate();

  // determine initial values
  if (!_bas.empty())
  {
    // select columns of M corresponding to z vars in the basis
    M.select(_all.begin(), _all.end(), _bas.begin(), _bas.end(), _t1);

    // select columns of I corresponding to z vars not in the basis
    _Bl.select(_all.begin(), _all.end(), _nonbas.begin(), _nonbas.end(), _t2);

    // setup the basis matrix
    _Bl.resize(n, _t1.columns() + _t2.columns());
    _Bl.set_sub_mat(0,0,_t1);
    _Bl.set_sub_mat(0,_t1.columns(),_t2);
  }

  // solve B*x = -q
  unsigned basis_count = std::numeric_limits<unsigned>::max();
  while (true)
  {
    try
    {
      _Al = _Bl;
      _x = q;
      _LA.solve_fast(_Al, _x);
      break;
    }
    catch (SingularException e)
    {
      // if initial basis didn't work, prepare to iterate through all 2^n bases
      if (basis_count == std::numeric_limits<unsigned>::max())
        basis_count = 0;
        
      // cycle through all bases until we find one that works
      set_basis(n, basis_count++, _bas, _nonbas);

      // select columns of M corresponding to z vars in the basis
      M.select(_all.begin(), _all.end(), _bas.begin(), _bas.end(), _t1);

      // select columns of I corresponding to z vars not in the basis
      _Bl.select(_all.begin(), _all.end(), _nonbas.begin(), _nonbas.end(), _t2);

      // setup the basis matrix
      _Bl.resize(n, _t1.columns() + _t2.columns());
      _Bl.set_sub_mat(0,0,_t1);
      _Bl.set_sub_mat(0,_t1.columns(),_t2);
    }
  }
  _x.negate();

  // check whether initial basis provides a solution
  if (std::find_if(_x.begin(), _x.end(), std::bind2nd(std::less<double>(), 0.0)) == _x.end())
  {
    for (idx = 0, iiter = _bas.begin(); iiter != _bas.end(); iiter++, idx++)
      z[*iiter] = _x[idx];
    z.resize(n, true);

    // check to see whether tolerances are satisfied
    FILE_LOG(LOG_OPT) << " -- initial basis provides a solution!" << std::endl;
    if (LOGGING(LOG_OPT))
    {
      M.mult(z, _wl) += q;
      double minw = *std::min_element(_wl.begin(), _wl.end());
      double w_dot_z = std::fabs(_wl.dot(z));
      FILE_LOG(LOG_OPT) << "  z: " << z << std::endl;
      FILE_LOG(LOG_OPT) << "  _w: " << _wl << std::endl;
      FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
      FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
    }
    FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

    return true; 
  }

  // use a new pivot tolerance if necessary
  const double PIV_TOL = (piv_tol > (double) 0.0) ? piv_tol : std::numeric_limits<double>::epsilon() * n * std::max((double) 1.0, M.norm_inf());

  // determine initial leaving variable
  ColumnIteratord min_x = std::min_element(_x.begin(), _x.begin() + n);
  double tval = -*min_x;
  BOOST_FOREACH(unsigned i, _nonbas) // add w variables to basis
    _bas.push_back(i+n);
  lvindex = std::distance(_x.begin(), min_x);
  iiter = _bas.begin();
  std::advance(iiter, lvindex);
  leaving = *iiter;
  FILE_LOG(LOG_OPT) << " -- x: " << _x << endl;
  FILE_LOG(LOG_OPT) << " -- first pivot: leaving index=" << lvindex << "  entering index=" << entering << " minimum value: " << tval << endl;

  // pivot in the artificial variable
  *iiter = t;    // replace w var with _z0 in basic indices
  _u.resize(n);
  for (unsigned i=0; i< n; i++)
    _u[i] = (_x[i] < 0.0) ? 1.0 : 0.0;
  _Bl.mult(_u, _Be);
  _Be.negate();
  _u *= tval;
  _x += _u;
  _x[lvindex] = tval;
  _Bl.set_column(lvindex, _Be);
  FILE_LOG(LOG_OPT) << "  new q: " << _x << endl;

  // main iterations begin here
  for (unsigned iter=0; iter < MAXITER; iter++)
  {
    if (LOGGING(LOG_OPT))
    {
      std::ostringstream basic;
      for (unsigned i=0; i< _bas.size(); i++)
        basic << " " << _bas[i];
      FILE_LOG(LOG_OPT) << "basic variables:" << basic.str() << std::endl;
    }

    // check whether done; if not, get new entering variable
    if (leaving == t)
    {
      FILE_LOG(LOG_OPT) << "-- solved LCP successfully!" << endl;
      unsigned idx;
      for (idx = 0, iiter = _bas.begin(); iiter != _bas.end(); iiter++, idx++)
        z[*iiter] = _x[idx];
      z.resize(n, true);

      // verify tolerances
      if (LOGGING(LOG_OPT))
      {
        M.mult(z, _wl) += q;
        double minw = *std::min_element(_wl.begin(), _wl.end());
        double w_dot_z = std::fabs(_wl.dot(z));
        FILE_LOG(LOG_OPT) << "  found solution!" << std::endl;
        FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
        FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
      }
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

      return true; 
    }
    else if (leaving < n)
    {
      entering = n + leaving;
      _Be.set_zero(n);
      _Be[leaving] = -1;
    }
    else
    {
      entering = leaving - n;
      M.get_column(entering, _Be);
    }
    _dl = _Be;
    try
    {
      _Al = _Bl;
      _LA.solve_fast(_Al, _dl);
    }
    catch (SingularException e)
    {
      FILE_LOG(LOG_OPT) << " -- warning: regular linear system solver failed; trying LS solver..." << std::endl;
      try
      {
        // use slower SVD pseudo-inverse
        _Al = _Bl;
        _dl = _Be;
        _LA.solve_LS_fast1(_Al, _dl);
      }
      catch (NumericalException e)
      {
        _Al = _Bl;
        _LA.solve_LS_fast2(_Al, _dl);
      }
    }

    // ** find new leaving variable
    _j.clear();
    for (unsigned i=0; i< _dl.size(); i++)
      if (_dl[i] > PIV_TOL)
        _j.push_back(i);
    // check for no new pivots; ray termination
    if (_j.empty())
    {
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() - no new pivots (ray termination)" << endl;
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

      z.resize(n, true);
      return false;
    }

    if (LOGGING(LOG_OPT))
    {
      std::ostringstream j;
      for (unsigned i=0; i< _j.size(); i++)
        j << " " << _j[i];
      FILE_LOG(LOG_OPT) << "d: " << _dl << std::endl;
      FILE_LOG(LOG_OPT) << "j (before min ratio):" << j.str() << std::endl;
    }

    // select elements j from x and d
    _xj.resize(_j.size());
    _dj.resize(_xj.size());
    select(_x.begin(), _j.begin(), _j.end(), _xj.begin());
    select(_dl.begin(), _j.begin(), _j.end(), _dj.begin());

    // compute minimal ratios x(j) + EPS_DOUBLE ./ d(j), d > 0
    _result.resize(_xj.size());
    std::transform(_xj.begin(), _xj.end(), _result.begin(), std::bind2nd(std::plus<double>(), zero_tol));
    std::transform(_result.begin(), _result.end(), _dj.begin(), _result.begin(), std::divides<double>());
    double theta = *std::min_element(_result.begin(), _result.end());

    // NOTE: lexicographic ordering does not appear to be used here to prevent
    // cycling (see [Cottle 1992], pp. 340-342)
    // find indices of minimal ratios, d> 0
    //   divide _x(j) ./ d(j) -- remove elements above the minimum ratio
    std::transform(_xj.begin(), _xj.end(), _dj.begin(), _result.begin(), std::divides<double>());
    for (iiter = _j.begin(), idx = 0; iiter != _j.end(); )
      if (_result[idx++] <= theta)
        iiter++;
      else
        iiter = _j.erase(iiter);

    if (LOGGING(LOG_OPT))
    {
      std::ostringstream j;
      for (unsigned i=0; i< _j.size(); i++)
        j << " " << _j[i];
      FILE_LOG(LOG_OPT) << "j (after min ratio):" << j.str() << std::endl;
    }

    // if j is empty, then likely the zero tolerance is too low
    if (_j.empty())
    {
      FILE_LOG(LOG_OPT) << "zero tolerance too low?" << std::endl;
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << std::endl;
      z.resize(n, true);
      return false;
    }

    // check whether artificial index among these
    _tlist.clear();
    select(_bas.begin(), _j.begin(), _j.end(), std::back_inserter(_tlist));
    if (std::find(_tlist.begin(), _tlist.end(), t) != _tlist.end())
    {
      iiter = std::find(_bas.begin(), _bas.end(), t); 
      lvindex = iiter - _bas.begin();
    }
    else
    {
      // several indices pass the minimum ratio test, pick one randomly
//      lvindex = _j[rand() % _j.size()];
// NOTE: solver seems *much* more capable of solving when we pick the first
// element rather than picking a random one
      lvindex = _j[0];
    }

    // set leaving = bas(lvindex)
    iiter = _bas.begin();
    std::advance(iiter, lvindex);
    leaving = *iiter;

    // ** perform pivot
    double ratio = _x[lvindex]/_dl[lvindex];
    _dl *= ratio;
    _x -= _dl;
    _x[lvindex] = ratio;
    _Bl.set_column(lvindex, _Be);
    *iiter = entering;
    FILE_LOG(LOG_OPT) << " -- pivoting: leaving index=" << lvindex << "  entering index=" << entering << endl;
  }

  FILE_LOG(LOG_OPT) << " -- maximum number of iterations exceeded" << endl;
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << std::endl;

  // max iterations exceeded
  z.resize(n, true);
  
  return false;
}

/// Regularized wrapper around Lemke's algorithm for srpase matrices
bool LCP::lcp_lemke_regularized(const SparseMatrixNd& M, const VectorNd& q, VectorNd& z, int min_exp, unsigned step_exp, int max_exp, double piv_tol, double zero_tol)
{
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() entered" << endl;

  // look for fast exit
  if (q.size() == 0)
  {
    z.resize(0);
    return true;
  }

  // copy MM
  _MMs = M;

  // assign value for zero tolerance, if necessary
  const double ZERO_TOL = (zero_tol > (double) 0.0) ? zero_tol : q.size() * std::numeric_limits<double>::epsilon();

  // try non-regularized version first
  bool result = lcp_lemke(_MMs, q, z, piv_tol, zero_tol);
  if (result)
  {
    // verify that solution truly is a solution -- check z
    if (*std::min_element(z.begin(), z.end()) >= -ZERO_TOL)
    {
      // check w
      M.mult(z, _wx) += q;
      if (*std::min_element(_wx.begin(), _wx.end()) >= -ZERO_TOL)
      {
        // check z'w
        std::transform(z.begin(), z.end(), _wx.begin(), _wx.begin(), std::multiplies<double>());
        pair<ColumnIteratord, ColumnIteratord> mmax = boost::minmax_element(_wx.begin(), _wx.end());
        if (*mmax.first >= -ZERO_TOL && *mmax.second < ZERO_TOL)
        {
          FILE_LOG(LOG_OPT) << "  solved with no regularization necessary!" << endl;
          FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;
          return true;
        }
      }
    }
  }

  // add a zero sparse diagonal matrix to _MMs
  _eye = SparseMatrixNd::identity(q.size());
  (_zero = _eye) *= 0.0;
  _MMs += _zero; 

  // start the regularization process
  int rf = min_exp;
  while (rf < max_exp)
  {
    // setup regularization factor
    double lambda = std::pow((double) 10.0, (double) rf);
    (_diag_lambda = _eye) *= lambda;

    // regularize M
    (_MMx = _MMs) += _diag_lambda;

    // try to solve the LCP
    if ((result = lcp_lemke(_MMx, q, z, piv_tol, zero_tol)))
    {
      // verify that solution truly is a solution -- check z
      if (*std::min_element(z.begin(), z.end()) > -ZERO_TOL)
      {
        // check w
        _MMx.mult(z, _wx) += q;
        if (*std::min_element(_wx.begin(), _wx.end()) > -ZERO_TOL)
        {
          // check z'w
          std::transform(z.begin(), z.end(), _wx.begin(), _wx.begin(), std::multiplies<double>());
          pair<ColumnIteratord, ColumnIteratord> mmax = boost::minmax_element(_wx.begin(), _wx.end());
          if (*mmax.first > -ZERO_TOL && *mmax.second < ZERO_TOL)
          {
            FILE_LOG(LOG_OPT) << "  solved with regularization factor: " << lambda << endl;
            FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;

            return true;
          }
        }
      }
    }

    // increase rf
    rf += step_exp;
  }

  FILE_LOG(LOG_OPT) << "  unable to solve given any regularization!" << endl;
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke_regularized() exited" << endl;

  // still here?  failure...
  return false;
}

/// Lemke's algorithm for solving linear complementarity problems using sparse matrices
/**
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
bool LCP::lcp_lemke(const SparseMatrixNd& M, const VectorNd& q, VectorNd& z, double piv_tol, double zero_tol)
{
  const unsigned n = q.size();
  const unsigned MAXITER = std::min((unsigned) 1000, 50*n);

  // look for immediate exit
  if (n == 0)
  {
    z.resize(0);
    return true;
  }

  // clear all vectors
  _all.clear();
  _tlist.clear();
  _bas.clear();
  _nonbas.clear();
  _j.clear();

  // copy z to z0
  _z0 = z;

  // come up with a sensible value for zero tolerance if none is given
  if (zero_tol <= (double) 0.0)
    zero_tol = std::numeric_limits<double>::epsilon() * M.norm_inf() * n;

  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() entered" << endl;
  FILE_LOG(LOG_OPT) << "  M: " << endl << M;
  FILE_LOG(LOG_OPT) << "  q: " << q << endl;

  // see whether trivial solution exists
  if (*std::min_element(q.begin(), q.end()) > -zero_tol)
  {
    FILE_LOG(LOG_OPT) << " -- trivial solution found" << endl;
    FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;
    z.set_zero(n);
    return true;
  }

  // initialize variables
  z.set_zero(n*2);
  unsigned t = 2*n;
  unsigned entering = t;
  unsigned leaving = 0;
  _all.clear();
  for (unsigned i=0; i< n; i++)
    _all.push_back(i);
  unsigned lvindex;
  unsigned idx;
  vector<unsigned>::iterator iiter;
  _tlist.clear();

  // determine initial basis
  _bas.clear();
  _nonbas.clear();
  if (_z0.size() != n)
    for (unsigned i=0; i< n; i++)
      _nonbas.push_back(i);
  else
    for (unsigned i=0; i< n; i++)
      if (_z0[i] > 0)
        _bas.push_back(i);
      else
        _nonbas.push_back(i);

  // determine initial values
  if (!_bas.empty())
  {
    typedef map<pair<unsigned, unsigned>, double> ValueMap;
    ValueMap values, newvalues;

    // select columns of M corresponding to z vars in the basis
    M.get_values(values);
    for (ValueMap::const_iterator i = values.begin(); i != values.end(); i++)
    {
      vector<unsigned>::const_iterator j = std::find(_bas.begin(), _bas.end(), i->first.second);
      if (j == _bas.end())
        continue;
      else
        newvalues[make_pair(i->first.first, j - _bas.begin())] = i->second;
    }

    // "select" columns of eye corresponding to z vars not in the basis
    // select_columns(_nonbas.begin(), _nonbas.end(), _st2);
    for (unsigned i=0, j=_bas.size(); i< _nonbas.size(); i++, j++)
      newvalues[make_pair(_nonbas[i],j)] = 1.0;

    // setup the basis matrix
    _sBl = SparseMatrixNd(SparseMatrixNd::eCSC, n, n, newvalues);
  }
  else
  {
    _sBl = SparseMatrixNd::identity(SparseMatrixNd::eCSC, n);
    _sBl.negate();
  }

  // solve B*x = -q
  _LA.solve_sparse_direct(_sBl, q, Ravelin::eNoTranspose, _x);
  _x.negate();

  // check whether initial basis provides a solution
  if (std::find_if(_x.begin(), _x.end(), std::bind2nd(std::less<double>(), 0.0)) == _x.end())
  {
    for (idx = 0, iiter = _bas.begin(); iiter != _bas.end(); iiter++, idx++)
      z[*iiter] = _x[idx];
    z.resize(n, true);

    // check to see whether tolerances are satisfied
    FILE_LOG(LOG_OPT) << " -- initial basis provides a solution!" << std::endl;
    if (LOGGING(LOG_OPT))
    {
      M.mult(z, _wl) += q;
      double minw = *std::min_element(_wl.begin(), _wl.end());
      double w_dot_z = std::fabs(_wl.dot(z));
      FILE_LOG(LOG_OPT) << "  z: " << z << std::endl;
      FILE_LOG(LOG_OPT) << "  w: " << _wl << std::endl;
      FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
      FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
    }
    FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

    return true; 
  }

  // determine initial leaving variable
  ColumnIteratord min_x = std::min_element(_x.begin(), _x.begin() + n);
  double tval = -*min_x;
  BOOST_FOREACH(unsigned i, _nonbas) // add w variables to basis
    _bas.push_back(i+n);
  lvindex = std::distance(_x.begin(), min_x);
  iiter = _bas.begin();
  std::advance(iiter, lvindex);
  leaving = *iiter;

  // pivot in the artificial variable
  *iiter = t;    // replace w var with _z0 in basic indices
  _u.resize(n);
  for (unsigned i=0; i< n; i++)
    _u[i] = (_x[i] < 0.0) ? 1.0 : 0.0;
  _sBl.mult(_u, _Be);
  _Be.negate();
  _u *= tval;
  _x += _u;
  _x[lvindex] = tval;
  _sBl.set_column(lvindex, _Be);
  FILE_LOG(LOG_OPT) << "  new q: " << _x << endl;

  // main iterations begin here
  for (unsigned iter=0; iter < MAXITER; iter++)
  {
    // check whether done; if not, get new entering variable
    if (leaving == t)
    {
      FILE_LOG(LOG_OPT) << "-- solved LCP successfully!" << endl;
      unsigned idx;
      for (idx = 0, iiter = _bas.begin(); iiter != _bas.end(); iiter++, idx++)
        z[*iiter] = _x[idx];
      z.resize(n, true);

      // verify tolerances
      if (LOGGING(LOG_OPT))
      {
        M.mult(z, _wl) += q;
        double minw = *std::min_element(_wl.begin(), _wl.end());
        double w_dot_z = std::fabs(_wl.dot(z));
        FILE_LOG(LOG_OPT) << "  found solution!" << std::endl;
        FILE_LOG(LOG_OPT) << "  minimum w: " << minw << std::endl;
        FILE_LOG(LOG_OPT) << "  w'z: " << w_dot_z << std::endl;
      }
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

      return true; 
    }
    else if (leaving < n)
    {
      entering = n + leaving;
      _Be.set_zero(n);
      _Be[leaving] = -1;
    }
    else
    {
      entering = leaving - n;
      M.get_column(entering, _Be);
    }
    _LA.solve_sparse_direct(_sBl, _Be, Ravelin::eNoTranspose, _dl);

    // use a new pivot tolerance if necessary
    const double PIV_TOL = (piv_tol > (double) 0.0) ? piv_tol : std::numeric_limits<double>::epsilon() * n * std::max((double) 1.0, _Be.norm_inf());

    // ** find new leaving variable
    _j.clear();
    for (unsigned i=0; i< _dl.size(); i++)
      if (_dl[i] > PIV_TOL)
        _j.push_back(i);
    // check for no new pivots; ray termination
    if (_j.empty())
    {
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() - no new pivots (ray termination)" << endl;
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << endl;

      z.resize(n, true);
      return false;
    }

    FILE_LOG(LOG_OPT) << " -- column of M': " << _dl << endl;

    // select elements j from x and d
    _xj.resize(_j.size());
    _dj.resize(_xj.size());
    select(_x.begin(), _j.begin(), _j.end(), _xj.begin());
    select(_d.begin(), _j.begin(), _j.end(), _dj.begin());

    // compute minimal ratios x(j) + EPS_DOUBLE ./ d(j), d > 0
    _result.resize(_xj.size());
    std::transform(_xj.begin(), _xj.end(), _result.begin(), std::bind2nd(std::plus<double>(), zero_tol));
    std::transform(_result.begin(), _result.end(), _dj.begin(), _result.begin(), std::divides<double>());
    double theta = *std::min_element(_result.begin(), _result.end());

    // NOTE: lexicographic ordering does not appear to be used here to prevent
    // cycling (see [Cottle 1992], pp. 340-342)
    // find indices of minimal ratios, d> 0
    //   divide _x(j) ./ d(j) -- remove elements above the minimum ratio
    std::transform(_xj.begin(), _xj.end(), _dj.begin(), _result.begin(), std::divides<double>());
    for (iiter = _j.begin(), idx = 0; iiter != _j.end(); )
      if (_result[idx++] <= theta)
        iiter++;
      else
        iiter = _j.erase(iiter);

    // if j is empty, then likely the zero tolerance is too low
    if (_j.empty())
    {
      FILE_LOG(LOG_OPT) << "zero tolerance too low?" << std::endl;
      FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << std::endl;
      z.resize(n, true);
      return false;
    }

    // check whether artificial index among these
    _tlist.clear();
    select(_bas.begin(), _j.begin(), _j.end(), std::back_inserter(_tlist));
    if (std::find(_tlist.begin(), _tlist.end(), t) != _tlist.end())
    {
      iiter = std::find(_bas.begin(), _bas.end(), t); 
      lvindex = iiter - _bas.begin();
    }
    else
    {
      // several indices pass the minimum ratio test, pick one randomly
//      lvindex = _j[rand() % _j.size()];

      // NOTE: solver seems *much* more capable of solving when we pick the first
      // element rather than picking a random one
      lvindex = _j[0];
    }

    // set leaving = bas(lvindex)
    iiter = _bas.begin();
    std::advance(iiter, lvindex);
    leaving = *iiter;

    // ** perform pivot
    double ratio = _x[lvindex]/_dl[lvindex];
    _dl *= ratio;
    _x -= _dl;
    _x[lvindex] = ratio;
    _sBl.set_column(lvindex, _Be);
    *iiter = entering;
    FILE_LOG(LOG_OPT) << " -- pivoting: leaving index=" << lvindex << "  entering index=" << entering << endl;
  }

  FILE_LOG(LOG_OPT) << " -- maximum number of iterations exceeded" << endl;
  FILE_LOG(LOG_OPT) << "LCP::lcp_lemke() exited" << std::endl;

  // max iterations exceeded
  z.resize(n, true);
  
  return false;
}


