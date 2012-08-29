/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

using std::endl;
using std::vector;

/// Method for 4th-order implicit Runge-Kutta integration
template <class T>
void Rosenbrock4Integrator<T>::integrate(T& x, T (*f)(const T&, Real, Real, void*), Real& time, Real step_size, void* data)
{
  // determine desired time
  Real tdes = time + step_size;

  while (time < tdes)
  {
    step(x, f, time, step_size, step_size, data);
    if (time + step_size > tdes)
      step_size = tdes - time;
  }

  return;
}

/// Does all the work
template <class T>
void Rosenbrock4Integrator<T>::step(T& x, T (*f)(const T&, Real, Real, void*), Real& time, Real step_size, Real& tnext, void* data)
{
  SAFESTATIC MatrixNN A, J;
  SAFESTATIC VectorN xscal, g1, g2, g3, g4, xsav, dxsav;
  SAFESTATIC vector<int> ipiv;
  const unsigned n = x.size();

  // general constants
  const Real SAFETY = (Real) 0.9;
  const Real GROW = (Real) 1.5;
  const Real PGROW = (Real) -.25;
  const Real SHRINK = (Real) 0.5;
  const Real PSHRINK = (Real) -1.0/3.0;
  const Real ERRCON = (Real) 0.1296;
  const unsigned MAXTRY = 40;

  // Shampine constants
  const Real GAM = (Real) 0.5;
  const Real A21 = (Real) 2.0;
  const Real A31 = (Real) 48.0/25.0;
  const Real A32 = (Real) 6.0/25.0;
  const Real C21 = (Real) -8.0;
  const Real C31 = (Real) 372.0/25.0;
  const Real C32 = (Real) 12.0/5.0;
  const Real C41 = (Real) -112.0/125.0;
  const Real C42 = (Real) -54.0/125.0;
  const Real C43 = (Real) -2.0/5.0;
  const Real B1 = (Real) 19.0/9.0;
  const Real B2 = (Real) 0.5;
  const Real B3 = (Real) 25.0/108.0;
  const Real B4 = (Real) 125.0/108.0;
  const Real E1 = (Real) 17.0/54.0;
  const Real E2 = (Real) 7.0/36.0;
  const Real E3 = (Real) 0.0;
  const Real E4 = (Real) 125.0/108.0;
  const Real C1X = (Real) 0.5;
  const Real C2X = (Real) -3.0/2.0;
  const Real C3X = (Real) 121.0/50.0;
  const Real C4X = (Real) 29.0/250.0;
  const Real A2X = (Real) 1.0;
  const Real A3X = (Real) 3.0/5.0;

  // setup g1, g2, g3, g4
  g1.resize(n);
  g2.resize(n);
  g3.resize(n);
  g4.resize(n);

  // setup xscal
  xscal.resize(n);
  for (unsigned i=0; i< n; i++)
    xscal[i] = std::max(rel_err_tol, std::fabs(x[i]));

  // get the derivative at time
  T dxdt = f(x, time, step_size, data);

  // save initial values
  Real tsav = time;
  xsav.copy_from(x);
  dxsav.copy_from(dxdt);

  // determine the vector df/dt
  Real sqrt_eps = std::sqrt(std::numeric_limits<Real>::epsilon());
  time += sqrt_eps;
  T dfdt = f(x, time, step_size, data);
  time -= sqrt_eps;

  // determine the Jacobian using forward differencing
  J.resize(n);
  for (unsigned i=0; i< n; i++)
  {
    x[i] += sqrt_eps;
    T column_i = f(x, time, step_size, data);
    x[i] -= sqrt_eps;
    for (unsigned j=0; j< n; j++)
      J(j, i) = column_i[j];
  }

  // set stepsize to initial trial value
  Real h = step_size;

  // init A matrix and pivots for LU decomposition
  A.resize(n);
  ipiv.resize(n);

  // iterate up to the maximal number of tries
  for (unsigned j=0; j< MAXTRY; j++)
  {
    // setup the matrix 1 - gamma*h*f'
    for (unsigned i=0; i< n; i++)
    {
      for (unsigned k=0; k< n; k++)
        A(i,k) = -J(i,k);
      A(i,i) += (Real) 1.0/(GAM*h);
    }
    
    // do LU decomposition of the matrix A
    LinAlg::factor_LU(A, ipiv);

    // setup r.h.s. for g1
    for (unsigned i=0; i< n; i++)
      g1[i] = dxsav[i] + h*C1X*dfdt[i];

    // solve for g1
    LinAlg::solve_LU_fast(A, false, ipiv, g1);

    // compute intermediate values of x and t
    for (unsigned i=0; i< n; i++)
      x[i] = g1[i]*A21 + xsav[i];
    time = tsav + A2X*h;

    // compute dxdt at the intermediate values
    dxdt = f(x, time, step_size, data);

    // setup r.h.s. for g2
    for (unsigned i=0; i< n; i++)
      g2[i] = dxdt[i] + h*C2X*dfdt[i]+C21*g1[i]/h;

    // solve for g2
    LinAlg::solve_LU_fast(A, false, ipiv, g2);

    // compute intermediate values of x and t
    for (unsigned i=0; i< n; i++)
      x[i] = xsav[i] + A31*g1[i]+A32*g2[i];
    time = tsav + A3X*h;

    // compute dxdt at the intermediate values
    dxdt = f(x, time, step_size, data);

    // setup r.h.s. for g3
    for (unsigned i=0; i< n; i++)
      g3[i] = dxdt[i] + h*C3X*dfdt[i] + (C31*g1[i]+C32*g2[i])/h;

    // solve for g3
    LinAlg::solve_LU_fast(A, false, ipiv, g3);

    // setup r.h.s. for g4
    for (unsigned i=0; i< n; i++)
      g4[i] = dxdt[i] + h*C4X*dfdt[i] + (C41*g1[i] + C42*g2[i] + C43*g3[i])/h;

    // solve for g4
    LinAlg::solve_LU_fast(A, false, ipiv, g4);

    // get fourth-order estimate of x and error estimate
    VectorN err(n);
    for (unsigned i=0; i< n; i++)
    {
      x[i] = xsav[i] + B1*g1[i] + B2*g2[i] + B3*g3[i] + B4*g4[i];
      err[i] = E1*g1[i] + E2*g2[i] + E3*g3[i] + E4*g4[i];
    }

    // update t
    time = tsav + h;

    // evaluate accuracy
    Real errmax = (Real) 0.0;
    for (unsigned i=0; i< n; i++)
      errmax = std::max(errmax, std::fabs(err[i]/xscal[i]));
    errmax /= rel_err_tol; 
    if (errmax <= (Real) 1.0)
    {
      // step succeeded; compute size of next step and return
      tnext = (errmax > ERRCON ? SAFETY*h*std::pow(errmax,PGROW) : GROW*h);
      return;
    }
    else
    {
      // truncation error too large, reduce step size
      tnext = SAFETY*h*std::pow(errmax, PSHRINK);
      h = (h >= (Real) 0.0 ? std::max(tnext, SHRINK*h) : std::min(tnext, SHRINK*h));
    }
  }

  std::cerr << "Rosenbrock4Integrator::step() - maximum number of iterations exceeded (40)" << endl;
}

template <class T>
void Rosenbrock4Integrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "Rosenbrock4Integrator") == 0);
  Integrator<T>::load_from_xml(node, id_map); 
  
  // get the relative error tolerance
  const XMLAttrib* rerr_attrib = node->get_attrib("rel-err-tol");
  if (rerr_attrib)
    rel_err_tol = rerr_attrib->get_real_value();
}

template <class T>
void Rosenbrock4Integrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator<T>::save_to_xml(node, shared_objects); 
  node->name = "Rosenbrock4Integrator";
  node->attribs.insert(XMLAttrib("rel-err-tol", rel_err_tol));
}


