/// Generic += function
template <class U, class V>
void add_to(const U& source, V& target)
{
  assert(source.size() == target.size());
  Ravelin::ColumnIteratord_const sourcei = source.begin();
  Ravelin::ColumnIteratord targeti = target.begin();
  while (sourcei != source.end())
  {
    *targeti += *sourcei;
    targeti++;
    sourcei++; 
  }
}

/// Generic copy function
template <class U, class V>
void copy_to(const U& source, V& target)
{
  target.resize(source.size());
  Ravelin::ColumnIteratord_const sourcei = source.begin();
  Ravelin::ColumnIteratord targeti = target.begin();
  while (sourcei != source.end())
  {
    *targeti = *sourcei;
    targeti++;
    sourcei++; 
  }
}

/// Generic dot product function
template <class U, class V>
double dot(const U& u, const V& v)
{
  assert(u.size() == v.size());
  double sum = 0.0;
  Ravelin::ColumnIteratord_const ui = u.begin();
  Ravelin::ColumnIteratord_const vi = v.begin();
  while (ui != u.end())
  {
    sum += *vi * *ui;
    vi++;
    ui++; 
  }
  return sum;
}

/// Finds a Cauchy point for gradient projection QP method
template <class M, class U, class V, class W>
V& find_cauchy_point(const M& G, const U& c, const V& l, const V& u, const Ravelin::VectorNd& gradient, W& x)
{
  const double INF = std::numeric_limits<double>::max();

  if (G.rows() != G.columns())
    throw Ravelin::NonsquareMatrixException();

  // get size of gradient / x
  const unsigned n = x.size();
  assert(gradient.size() == n);
  assert(l.size() == n);
  assert(u.size() == n);

  // add the value 0 to t
  _alphas.clear();
  _alphas.push_back((double) 0.0);

  // identify values of t for which each component reaches its bound along -g
  for (unsigned i=0; i< n; i++)
  {
    if (gradient[i] < (double) 0.0 && u[i] < INF)
    {
      double value = (x[i] - u[i])/gradient[i];
      if (value > (double) 0.0)
        _alphas.push_back(value);
    }
    else if (gradient[i] > (double) 0.0 && l[i] > -INF)
    {
      double value = (x[i] - l[i])/gradient[i];
      if (value > (double) 0.0)
        _alphas.push_back(value);
    }
    else
      _alphas.push_back(INF);
  }

  // sort alphas and remove duplicate values
  std::sort(_alphas.begin(), _alphas.end());
  _alphas.erase(std::unique(_alphas.begin(), _alphas.end()), _alphas.end());

  // determine first local minimizer
  _p.resize(n);
  for (unsigned j=1; j< _alphas.size(); j++)
  {
    unsigned jm1 = j-1;

    // compute p^{j-1}
    for (unsigned i=0; i< n; i++)
    {
      // (re)compute alphasi
      double alphasi;
      if (gradient[i] < (double) 0.0 && u[i] < INF)
        alphasi = (x[i] - u[i])/gradient[i];
      else if (gradient[i] > (double) 0.0 && l[i] > -INF)
        alphasi = (x[i] - l[i])/gradient[i];
      else
        alphasi = INF;

      _p[i] = (_alphas[jm1] < alphasi) ? -gradient[i] : (double) 0.0;
    }

    // look for local minimizer
    double fprime = dot(c, _p) + dot(x, G.mult(_p, _Gp));

    // if f' > 0 there is a local minimizer at t[j-1]
    if (fprime > 0)
      return x;
    else
    {
      // compute dt*
      double fpprime = dot(_p, _Gp);
      double dt_star = -fprime / fpprime;
      if (dt_star >= (double) 0.0 && dt_star < _alphas[j] - _alphas[jm1])
      {
        // there is a local minimizer at t[jm1] + dt_star
        _workv = _p;
        _workv *= dt_star;
        add_to(_workv, x);

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
        double dt = _alphas[j] - _alphas[jm1];
        if (!std::isinf(dt) && !std::isnan(dt))
        {
          _workv = _p;
          _workv *= dt;
          add_to(_workv, x);
        }
      }
    }
  }

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
template <class M, class V>
unsigned qp_gradproj(const M& G, const V& c, const V& l, const V& u, unsigned max_iter, V& x, double tol)
{
  if (G.rows() != G.columns())
    throw Ravelin::NonsquareMatrixException();

  const double INF = std::numeric_limits<double>::max();

  // minimum objective function value determined so far
  double qmin = INF;

  // setup the number of inner iterations
  unsigned niter = 0;

  // determine n
  const unsigned n = c.size();

  // pick starting x halfway between l and u, if necessary
  for (unsigned i=0; i< n; i++)
    if (x[i] < l[i] || x[i] > u[i])
      x[i] = l[i]*0.5 + u[i]*0.5;

  // the preconditioner
  _H.resize(n);
  for (unsigned i=0; i< n; i++)
    _H[i] = (G(i,i) > (double) 0.0) ? (double) 1.0/G(i,i) : (double) 1.0;

  // iterate..
  for (; niter < max_iter; niter++)
  {
    // compute the gradient -- if it is zero, quit now
    G.mult(x, _Gx_c) += c;
    if (_Gx_c.norm() < tol)
      break;
    
    // find Cauchy point
    find_cauchy_point(G, c, l, u, _Gx_c, x);

    // evaluate Cauchy point
    double q = dot(x, G.mult(x, _Gx)) * (double) 0.5 + dot(x, c);

    // determine which constraints are in the active set
    _inactive_set.clear();
    unsigned nactive = 0;
    for (unsigned i=0; i< n; i++)
      if (std::fabs(x[i] - l[i]) < std::numeric_limits<double>::epsilon() || 
          std::fabs(x[i] - u[i]) < std::numeric_limits<double>::epsilon())
        nactive++;
      else
        _inactive_set.push_back(i);

    // form preconditioner inv(_Wzz)
    _H.select(_inactive_set.begin(), _inactive_set.end(), _Wzz);

    // form Z'GZ
    G.select_square(_inactive_set.begin(), _inactive_set.end(), _ZtGZ);

    // prepare to do conjugate gradient iteration
    _dxz.resize(n-nactive);
    _xz.resize(n-nactive);
    _lz.resize(n-nactive);
    _uz.resize(n-nactive); 
    _cz.resize(n-nactive);
    for (unsigned i=0; i< _inactive_set.size(); i++)
    {
      unsigned j = _inactive_set[i];
      _cz[i] = c[j];
      _xz[i] = x[j];
      _lz[i] = l[j];
      _uz[i] = u[j];
    }

    // *** Solve Z'GZ*xz = -cz (subspace minimization)
    if (_cz.size() < 100)
    {
      _dxz = _cz;
      _dxz.negate();
      _LA.factor_chol(_ZtGZ);
      _LA.solve_chol_fast(_ZtGZ, _dxz);

      // determine delta xz
      _dxz -= _xz;
      double max_phi = 1.0;
      for (unsigned j=0; j< _xz.size(); j++)
      {
        if (_dxz[j] < (double) 0.0 && _lz[j] > -INF)
          max_phi = std::min(max_phi, (_lz[j] - _xz[j])/_dxz[j]);
        else if (_dxz[j] > (double) 0.0 && _uz[j] < INF)
          max_phi = std::min(max_phi, (_uz[j] - _xz[j])/_dxz[j]);
        assert(max_phi >= -NEAR_ZERO);
      }
      
      if (max_phi < (double) 0.0)
        max_phi = (double) 0.0;

      // update xz
      _dxz *= max_phi;
      _xz += _dxz;
    }     
    // *** do MINRES method (subspace minimization)
    else
    {
      (_workv = _cz).negate();
      _y = _workv;
      _r1 =_workv;
      M::diag_mult(_Wzz, _workv, _y);
      double beta1 = dot(_workv, _y);

      // look for x = 0 solution
      if (beta1 == 0)
        _xz.set_zero();
      else
      {
        // normalize y to get v1 later
        beta1 = std::sqrt(beta1);

        // init other quantities
        double oldb = (double) 0.0;
        double beta = beta1;
        double dbar = (double) 0.0;
        double epsln = (double) 0.0;
        double qrnorm = beta1;
        double phibar = beta1;
        double rhs1 = beta; 
        double rhs2 = (double) 0.0;
        double tnorm2 = (double) 0.0;
        double ynorm2 = (double) 0.0;
        double cs = (double) -1.0;
        double sn = (double) 0.0;
        _w.set_zero(_cz.size());
        _w2.set_zero(_cz.size());
        _r2 = _r1;

        // reset stop flag and setup gmax and gmin
        bool stop = false;
        double gmax = std::numeric_limits<double>::max();
        double gmin = -gmax;

        // do the iterative method
//      for (unsigned i=0; i< n-nactive; i++)
        for (unsigned i=0;;i++)
        {
          double s = (double) 1.0/beta;
          (_v = _y) *= s;
          _ZtGZ.mult(_v, _y);
          if (i > 0)
          {
            _workv = _r1;
            _workv *= (beta/oldb);
            _y -= _workv;
          }
          double alfa = dot(_v, _y);
          _workv = _r2;
          _workv *= (-alfa/beta);
          _y += _workv;
          _r1 = _r2;
          _r2 = _y;
          M::diag_mult(_Wzz, _r2, _y);
          oldb = beta;
          beta = dot(_r2, _y);
          if (beta < 0)
            break;
          beta = std::sqrt(beta);
          tnorm2 += alfa*alfa + oldb*oldb + beta*beta;

          // init a few things
          if (i == 0)
          {
            if (beta/beta1 <= 10*std::numeric_limits<double>::epsilon())
              stop = true;
            gmax = std::fabs(alfa);
            gmin = gmax;
          }

          // apply previous rotation
          double oldeps = epsln;
          double delta = cs*dbar + sn*alfa;
          double gbar = sn*dbar - cs*alfa;
          epsln = sn*beta;
          dbar = -cs*beta;
          double root = std::sqrt(gbar*gbar + dbar*dbar);

          // compute the next plane rotation
          double gamma = std::sqrt(gbar*gbar + beta*beta);
          gamma = std::max(gamma, std::numeric_limits<double>::epsilon());
          cs = gbar/gamma;
          sn = beta/gamma;
          double phi = cs*phibar;
          phibar = sn*phibar; 

          // prepare to update xz
          double denom = (double) 1.0/gamma;
          _w1 = _w2;
          _w2 = _w;
          (_workv = _w1) *= oldeps;
          (_w = _w2) *= delta; 
          _w += _v;
          _w += _workv;
          _w *= denom;
 
          // determine the maximum phi 
          double box_phi = std::numeric_limits<double>::max();
          for (unsigned j=0; j< _xz.size(); j++)
            if (_w[j] < (double) 0.0 && _lz[j] > -INF)
              box_phi = std::min(box_phi, (_lz[j] - _xz[j])/_w[j]);
            else if (_w[j] > (double) 0.0 && _uz[j] < INF)
              box_phi = std::min(box_phi, (_uz[j] - _xz[j])/_w[j]);

          // update xz
          (_workv = _w) *= std::min(box_phi, phi);
          _xz += _workv;

          // look for stopping condition
          if (box_phi < phi)
            break;

          // prepare to loop again
          gmax = std::max(gmax, gamma);
          gmin = std::min(gmin, gamma);
          double z = rhs1/gamma;
          ynorm2 = z*z + ynorm2;
          rhs1 = rhs2 - delta*z;
     
          // estimate various norms to see whether we can quit
          double Anorm = std::sqrt(tnorm2);
          double ynorm = std::sqrt(ynorm2);
          double epsx = Anorm*ynorm*std::numeric_limits<double>::epsilon();
          qrnorm = phibar;
          double rnorm = qrnorm;
          double test1 = rnorm/(Anorm*ynorm);
          double test2 = root/Anorm;

          // see whether we can stop
          if (stop)
            break;

          // check to see whether we can quit
          double t1 = (double) 1.0 + test1;
          double t2 = (double) 1.0 + test2;
          if (t1 <= (double) 1.0 || t2 <= (double) 1.0 || epsx >= beta1 || 
              test2 <= NEAR_ZERO || test1 <= NEAR_ZERO)
            break;
        }
      }
    }

    // project xz back to x
    copy_to(x,_xstar);
    for (unsigned i=0; i< _inactive_set.size(); i++)
      _xstar[_inactive_set[i]] = _xz[i];

    // manually enforce constraints
    for (unsigned i=0; i< _xstar.size(); i++)
      if (_xstar[i] < l[i])
        _xstar[i] = l[i];
      else if (_xstar[i] > u[i])
        _xstar[i] = u[i];

    // evaluate new x; if it is better than old x, keep it
    double qstar = dot(_xstar, G.mult(_xstar, _Gx)) * (double) 0.5 + dot(_xstar, c);
    if (qstar < q)
    {
      copy_to(_xstar, x);
      q = qstar;
    }

    // if objective function did not decrease sufficiently, quit
    if (qmin - q < tol)
      break;
    else // if (q <= qmin)
    {
      if (std::isnan(q))
        return niter;
      else
        assert(q <= qmin);
      qmin = q;
    }
  }

  return niter;
}


