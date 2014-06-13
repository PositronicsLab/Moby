/// Active-set QP algorithm for solving QPs 
/**
 * \param H the quadratic term
 * \param c the linear term
 * \param lb the lower bound
 * \param ub the lower bound
 * \param M the linear inequality constraints (M*z >= q)
 * \param q the linear inequality bound (M*z >= q)
 * \param A the linear equality constraint (A*z = b)
 * \param b the linear equality constraint (A*z = b)
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
template <class Mat1, class Vec1, class Vec2, class Vec3, class Mat2, class Vec4, class Mat3, class Vec5, class Vec6>
bool QLCPD::qp_activeset(const Mat1& H, const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat2& M, const Vec4& q, const Mat3& A, const Vec5& b, Vec6& z)
{
  int n = c.rows();
  const double INF = 1e+29;
  
  // resize z, if necessary
  z.resize((unsigned) n);

  // setup the constraint matrix
  assert(M.columns() == (unsigned) n);
  assert(A.columns() == (unsigned) n);
  _X.resize(z.rows(), 1+M.rows()+A.rows());
  _X.set_column(0, c);
  _X.set_sub_mat(0, 1, M, Ravelin::eTranspose);
  _X.set_sub_mat(0, 1+M.rows(), A, Ravelin::eTranspose);  

  // setup storage
  int la = (int) _X.leading_dim();
  int m = M.rows() + A.rows();

  // set k to dummy value (should not be used since mode < 2)
  int k = 0;
  int kmax = n;
  int maxg = std::min(6,kmax+1);

  // setup 'common' variables (generally tolerance parameters)
  set_values(n, m);

  // setup lower bounds and upper bounds
  assert(lb.size() == ub.size());
  assert(lb.size() == (unsigned) n);
  _lb.set_zero(lb.rows() + _X.columns());
  _ub.set_zero(ub.rows() + _X.columns());
  _lb.set_sub_vec(0, lb);
  _ub.set_sub_vec(0, ub);
  _lb.set_sub_vec((int) lb.rows(), q);
  _lb.set_sub_vec((int) (lb.rows()+q.rows()), b);
  std::fill_n(_ub.begin()+ub.rows(), q.size(), INF);
  _ub.set_sub_vec((int) (lb.rows()+q.rows()), b);

  // bound lower/upper bound constraints properly
  for (unsigned i=0; i< (unsigned) n; i++)
  {
    if (_ub[i] > INF)
      _ub[i] = INF;
    if (_lb[i] < -INF)
      _lb[i] = -INF;
  }

  // set the mode to warm start
  int mode = 1;

  // set maximum number of gradient calls
  int mxgr = 100000000;

  // set maximum level of recursion
  int mlp = std::min(20,std::max(m,3));

  // setup rgtol
  double rgtol = 1e-8;

  // setup lower bound on objective function
  double fmin = -INF;

  // indicate this problem is not a LP
  int linear = 0;

  // setup output options
  int iprint = 0;    // diagnostic output off
  int nout = 0;      // "channel" for output (stdout?)
  int ifail = 0;

  // setup workspaces
  _r.resize(n+m);
  _w.resize(n+m);
  _e.resize(n+m);
  _alp.resize(n+m);

  // setup linear workspaces
  _ls.resize(n+m);
  _lp.resize(mlp);

  // setup other miscellaneous variables
  int peq = 0;  // pointer to end of equality constraint indices in ls

  // initial estimates of eigenvalues of Hessian
  int nv = 1;
  _v.resize(maxg);
  _v[0] = 1.0;

  // setup workspace
  _ws.resize(_mxws);
  _lws.resize(_mxlws);
  std::copy(c.column_iterator_begin(), c.column_iterator_end(), _ws.begin());
  std::copy(H.column_iterator_begin(), H.column_iterator_end(), _ws.begin()+n);

  // resize the gradient
  _g.resize(n);

  // the function value at the end
  double f;

  try
  {
    // call the function
    qlcpd(&n, &m, &k, &kmax, &maxg, 
          _X.data(), &la, z.data(), _lb.data(), _ub.data(), 
          &f, &fmin, _g.data(), _r.data(), _w.data(), _e.data(), &_ls[0],
          _alp.data(), &_lp[0], &mlp, &peq, &_ws[0], &_lws[0], 
          _v.data(), &nv, &linear, &rgtol, &mode, &ifail, &mxgr, &iprint, &nout);
  }
  catch (std::runtime_error e)
  {
    return false;
  }

  // look whether failure is indicated
  if (ifail != 0)
    return false;

  // verify feasibility (M*z >= q)
  M.mult(z, _r) -= q;
  if (std::find_if(_r.begin(), _r.end(), std::bind2nd(std::less<double>(), -NEAR_ZERO)) != _r.end())
    return false;

  // verify feasibility (A*z = b)
  A.mult(z, _r) -= b;
  for (Ravelin::RowIteratord_const i = _r.row_iterator_begin(); i != _r.row_iterator_end(); i++)
    if (std::fabs(*i) > NEAR_ZERO)
      return false;

  // verify z >= lb and z <= ub
  for (unsigned i=0; i< n; i++)
    if ((z[i] < lb[i] && !rel_equal(z[i], lb[i])) ||
        (z[i] > ub[i] && !rel_equal(z[i], ub[i])))
      return false; 

  // all checks passed 
  return true; 
}

/// Active-set LP algorithm for finding closest point to feasibility 
/**
 * \param lb the lower bound
 * \param ub the lower bound
 * \param M the linear inequality constraints (M*z >= q)
 * \param q the linear inequality bound (M*z >= q)
 * \param A the linear equality constraint (A*z = b)
 * \param b the linear equality constraint (A*z = b)
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
template <class Vec1, class Vec2, class Mat1, class Vec3, class Mat2, class Vec4, class Vec5>
bool QLCPD::find_closest_feasible(const Vec1& lb, const Vec2& ub, const Mat1& M, const Vec3& q, const Mat2& A, const Vec4& b, Vec5& z)
{
  int n = lb.rows();
  
  const double INF = 1e+29;
  
  // resize z, if necessary
  z.resize((unsigned) n+1);

  // setup the constraint matrix
  assert(M.columns() == (unsigned) n);
  assert(A.columns() == (unsigned) n);
  _X.resize(z.rows(), 2+M.rows()+A.rows());
  _X.set_sub_mat(0, 1, M, Ravelin::eTranspose);
  _X.set_sub_mat(0, 1+M.rows(), A, Ravelin::eTranspose);  

  // setup c
  for (unsigned i=0; i< n; i++)
    _X(i,0) = 0.0;
  _X(n,0) = 1.0;

  // setup components of M with s (M affected by s)
  for (unsigned i=0, j=1; i< M.rows(); i++, j++)
    _X(n,j) = 1.0;

  // setup components of A with s (A not affected by s)
  for (unsigned i=0, j=1+M.rows(); i< A.rows(); i++, j++)
    _X(n,j) = 0.0; 

  // setup storage
  int la = (int) _X.leading_dim();
  int m = M.rows() + A.rows();

  // set k to dummy value (should not be used since mode < 2)
  int k = 0;
  int kmax = n+1;
  int maxg = std::min(6,kmax+1);

  // setup 'common' variables (generally tolerance parameters)
  set_values(n+1, m);

  // setup lower bounds and upper bounds
  assert(lb.size() == ub.size());
  assert(lb.size() == (unsigned) n);
  _lb.set_zero(lb.rows() + _X.columns() + 1);
  _ub.set_zero(ub.rows() + _X.columns() + 1);
  _lb.set_sub_vec(0, lb);
  _ub.set_sub_vec(0, ub);
  _lb.set_sub_vec((int) lb.rows()+1, q);
  _lb.set_sub_vec((int) (lb.rows()+q.rows()+1), b);
  std::fill_n(_ub.begin()+ub.rows()+1, q.size(), INF);
  _ub.set_sub_vec((int) (lb.rows()+q.rows()+1), b);

  // bound lower/upper bound constraints properly
  for (unsigned i=0; i< (unsigned) n; i++)
  {
    if (_ub[i] > INF)
      _ub[i] = INF;
    if (_lb[i] < -INF)
      _lb[i] = -INF;
  }

  // setup bound for s
  _lb[n] = 0.0;
  _ub[n] = INF;

  // set the mode to warm start
  int mode = 1;

  // set maximum number of gradient calls
  int mxgr = 100000000;

  // set maximum level of recursion
  int mlp = std::min(20,std::max(m,3));

  // setup rgtol
  double rgtol = 1e-8;

  // setup lower bound on objective function
  double fmin = -INF;

  // indicate this problem is a LP
  int linear = 1;

  // setup output options
  int iprint = 0;    // diagnostic output off
  int nout = 0;      // "channel" for output (stdout?)
  int ifail = 0;

  // setup workspaces
  _r.resize(n+m+1);
  _w.resize(n+m+1);
  _e.resize(n+m+1);
  _alp.resize(n+m+1);

  // setup linear workspaces
  _ls.resize(n+m+1);
  _lp.resize(mlp);

  // setup other miscellaneous variables
  int peq = 0;  // pointer to end of equality constraint indices in ls

  // initial estimates of eigenvalues of Hessian
  int nv = 1;
  _v.resize(maxg);
  _v[0] = 1.0;

  // setup workspace
  _ws.resize(_mxws);
  _lws.resize(_mxlws);

  // resize the gradient
  _g.resize(n+1);

  // the function value at the end
  double f;

  // call the function
  int nplus = n+1;
  try
  {
    qlcpd(&nplus, &m, &k, &kmax, &maxg, 
          _X.data(), &la, z.data(), _lb.data(), _ub.data(), 
          &f, &fmin, _g.data(), _r.data(), _w.data(), _e.data(), &_ls[0],
          _alp.data(), &_lp[0], &mlp, &peq, &_ws[0], &_lws[0], 
          _v.data(), &nv, &linear, &rgtol, &mode, &ifail, &mxgr, &iprint, &nout);
  }
  catch (std::runtime_error e)
  {
    return false;
  }

  // get the 's' value
  const double S = z[n];

  // resize z
  z.resize(n, true);

  // look whether failure is indicated
  if (ifail != 0)
    return false;

  // verify feasibility (M*z >= q)
  M.mult(z, _r) -= q;
  if (std::find_if(_r.begin(), _r.end(), std::bind2nd(std::less<double>(), -(S+NEAR_ZERO))) != _r.end())
    return false;

  // verify feasibility (A*z = b)
  A.mult(z, _r) -= b;
  for (Ravelin::RowIteratord_const i = _r.row_iterator_begin(); i != _r.row_iterator_end(); i++)
    if (std::fabs(*i) > NEAR_ZERO)
      return false;

  // verify z >= lb and z <= ub
  for (unsigned i=0; i< n; i++)
    if ((z[i] < lb[i] && !rel_equal(z[i], lb[i])) ||
        (z[i] > ub[i] && !rel_equal(z[i], ub[i])))
      return false; 

  // all checks passed 
  return true; 
}

/// Active-set QP algorithm for solving LPs 
/**
 * \param c the linear term
 * \param lb the lower bound
 * \param ub the lower bound
 * \param M the linear inequality constraints (M*z >= q)
 * \param q the linear inequality bound (M*z >= q)
 * \param A the linear equality constraint (A*z = b)
 * \param b the linear equality constraint (A*z = b)
 * \param z a vector "close" to the solution on input (optional); contains
 *        the solution on output
 */
template <class Vec1, class Vec2, class Vec3, class Mat1, class Vec4, class Mat2, class Vec5, class Vec6>
bool QLCPD::lp_activeset(const Vec1& c, const Vec2& lb, const Vec3& ub, const Mat1& M, const Vec4& q, const Mat2& A, const Vec5& b, Vec6& z)
{
  int n = c.rows();
  const double INF = 1e+29;
  
  // resize z, if necessary
  z.resize((unsigned) n);

  // setup the constraint matrix
  assert(M.columns() == (unsigned) n);
  assert(A.columns() == (unsigned) n);
  _X.resize(z.rows(), 1+M.rows()+A.rows());
  _X.set_column(0, c);
  _X.set_sub_mat(0, 1, M, Ravelin::eTranspose);
  _X.set_sub_mat(0, 1+M.rows(), A, Ravelin::eTranspose);  

  // setup storage
  int la = (int) _X.leading_dim();
  int m = M.rows() + A.rows();

  // set k to dummy value (should not be used since mode < 2)
  int k = 0;
  int kmax = n;
  int maxg = std::min(6,kmax+1);

  // setup 'common' variables (generally tolerance parameters)
  set_values(n, m);

  // setup lower bounds and upper bounds
  assert(lb.size() == ub.size());
  assert(lb.size() == (unsigned) n);
  _lb.set_zero(lb.rows() + _X.columns());
  _ub.set_zero(ub.rows() + _X.columns());
  _lb.set_sub_vec(0, lb);
  _ub.set_sub_vec(0, ub);
  _lb.set_sub_vec((int) lb.rows(), q);
  _lb.set_sub_vec((int) (lb.rows()+q.rows()), b);
  std::fill_n(_ub.begin()+ub.rows(), q.size(), INF);
  _ub.set_sub_vec((int) (lb.rows()+q.rows()), b);

  // bound lower/upper bound constraints properly
  for (unsigned i=0; i< (unsigned) n; i++)
  {
    if (_ub[i] > INF)
      _ub[i] = INF;
    if (_lb[i] < -INF)
      _lb[i] = -INF;
  }

  // set the mode to warm start
  int mode = 1;

  // set maximum number of gradient calls
  int mxgr = 100000000;

  // set maximum level of recursion
  int mlp = std::min(20,std::max(m,3));

  // setup rgtol
  double rgtol = 1e-8;

  // setup lower bound on objective function
  double fmin = -INF;

  // indicate this problem is a LP
  int linear = 1;

  // setup output options
  int iprint = 0;    // diagnostic output off
  int nout = 0;      // "channel" for output (stdout?)
  int ifail = 0;

  // setup workspaces
  _r.resize(n+m);
  _w.resize(n+m);
  _e.resize(n+m);
  _alp.resize(n+m);

  // setup linear workspaces
  _ls.resize(n+m);
  _lp.resize(mlp);

  // setup other miscellaneous variables
  int peq = 0;  // pointer to end of equality constraint indices in ls

  // initial estimates of eigenvalues of Hessian
  int nv = 1;
  _v.resize(maxg);
  _v[0] = 1.0;

  // setup workspace
  _ws.resize(_mxws);
  _lws.resize(_mxlws);
  std::copy(c.column_iterator_begin(), c.column_iterator_end(), _ws.begin());

  // resize the gradient
  _g.resize(n);

  // the function value at the end
  double f;

  // call the function
  try 
  {
    qlcpd(&n, &m, &k, &kmax, &maxg, 
          _X.data(), &la, z.data(), _lb.data(), _ub.data(), 
          &f, &fmin, _g.data(), _r.data(), _w.data(), _e.data(), &_ls[0],
          _alp.data(), &_lp[0], &mlp, &peq, &_ws[0], &_lws[0], 
          _v.data(), &nv, &linear, &rgtol, &mode, &ifail, &mxgr, &iprint, &nout);
  }
  catch (std::runtime_error e)
  {
    return false;
  }

  // look whether failure is indicated
  if (ifail != 0)
    return false;

  // verify feasibility (M*z >= q)
  M.mult(z, _r) -= q;
  if (std::find_if(_r.begin(), _r.end(), std::bind2nd(std::less<double>(), -NEAR_ZERO)) != _r.end())
    return false;

  // verify feasibility (A*z = b)
  A.mult(z, _r) -= b;
  for (Ravelin::RowIteratord_const i = _r.row_iterator_begin(); i != _r.row_iterator_end(); i++)
    if (std::fabs(*i) > NEAR_ZERO)
      return false;

  // verify z >= lb and z <= ub
  for (unsigned i=0; i< n; i++)
    if ((z[i] < lb[i] && !rel_equal(z[i], lb[i])) ||
        (z[i] > ub[i] && !rel_equal(z[i], ub[i])))
      return false; 

  // all checks passed 
  return true; 
}


