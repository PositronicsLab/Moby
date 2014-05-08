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

  // call the function
  qlcpd(&n, &m, &k, &kmax, &maxg, 
        _X.data(), &la, z.data(), _lb.data(), _ub.data(), 
        &f, &fmin, _g.data(), _r.data(), _w.data(), _e.data(), &_ls[0],
        _alp.data(), &_lp[0], &mlp, &peq, &_ws[0], &_lws[0], 
        _v.data(), &nv, &linear, &rgtol, &mode, &ifail, &mxgr, &iprint, &nout);

  // check result
  return (ifail == 0);
}


