/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/Rosenbrock4Integrator.h>

using std::endl;
using std::vector;
using boost::shared_ptr;
using Ravelin::VectorNd;
using Ravelin::MatrixNd;
using namespace Moby;

/// Method for 4th-order implicit Runge-Kutta integration
void Rosenbrock4Integrator::integrate(VectorNd& x, VectorNd (*f)(const VectorNd&, double, double, void*), double time, double step_size, void* data)
{
  // determine desired time
  double tdes = time + step_size;

  while (time < tdes)
  {
    step(x, f, time, step_size, step_size, data);
    if (time + step_size > tdes)
      step_size = tdes - time;
  }

  return;
}

/// Does all the work
void Rosenbrock4Integrator::step(VectorNd& x, VectorNd (*f)(const VectorNd&, double, double, void*), double& time, double step_size, double& tnext, void* data)
{
  SAFESTATIC MatrixNd A, J;
  SAFESTATIC VectorNd xscal, g1, g2, g3, g4, xsav, dxsav;
  SAFESTATIC vector<int> ipiv;
  const unsigned n = x.size();

  // general constants
  const double SAFETY = (double) 0.9;
  const double GROW = (double) 1.5;
  const double PGROW = (double) -.25;
  const double SHRINK = (double) 0.5;
  const double PSHRINK = (double) -1.0/3.0;
  const double ERRCON = (double) 0.1296;
  const unsigned MAXTRY = 40;

  // Shampine constants
  const double GAM = (double) 0.5;
  const double A21 = (double) 2.0;
  const double A31 = (double) 48.0/25.0;
  const double A32 = (double) 6.0/25.0;
  const double C21 = (double) -8.0;
  const double C31 = (double) 372.0/25.0;
  const double C32 = (double) 12.0/5.0;
  const double C41 = (double) -112.0/125.0;
  const double C42 = (double) -54.0/125.0;
  const double C43 = (double) -2.0/5.0;
  const double B1 = (double) 19.0/9.0;
  const double B2 = (double) 0.5;
  const double B3 = (double) 25.0/108.0;
  const double B4 = (double) 125.0/108.0;
  const double E1 = (double) 17.0/54.0;
  const double E2 = (double) 7.0/36.0;
  const double E3 = (double) 0.0;
  const double E4 = (double) 125.0/108.0;
  const double C1X = (double) 0.5;
  const double C2X = (double) -3.0/2.0;
  const double C3X = (double) 121.0/50.0;
  const double C4X = (double) 29.0/250.0;
  const double A2X = (double) 1.0;
  const double A3X = (double) 3.0/5.0;

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
  VectorNd dxdt = f(x, time, step_size, data);

  // save initial values
  double tsav = time;
  xsav = x;
  dxsav = dxdt;

  // determine the vector df/dt
  double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  time += sqrt_eps;
  VectorNd dfdt = f(x, time, step_size, data);
  time -= sqrt_eps;

  // determine the Jacobian using forward differencing
  J.resize(n,n);
  for (unsigned i=0; i< n; i++)
  {
    x[i] += sqrt_eps;
    VectorNd column_i = f(x, time, step_size, data);
    x[i] -= sqrt_eps;
    for (unsigned j=0; j< n; j++)
      J(j, i) = column_i[j];
  }

  // set stepsize to initial trial value
  double h = step_size;

  // init A matrix and pivots for LU decomposition
  A.resize(n,n);
  ipiv.resize(n);

  // iterate up to the maximal number of tries
  for (unsigned j=0; j< MAXTRY; j++)
  {
    // setup the matrix 1 - gamma*h*f'
    for (unsigned i=0; i< n; i++)
    {
      for (unsigned k=0; k< n; k++)
        A(i,k) = -J(i,k);
      A(i,i) += (double) 1.0/(GAM*h);
    }
    
    // do LU decomposition of the matrix A
    _LA.factor_LU(A, ipiv);

    // setup r.h.s. for g1
    for (unsigned i=0; i< n; i++)
      g1[i] = dxsav[i] + h*C1X*dfdt[i];

    // solve for g1
    _LA.solve_LU_fast(A, false, ipiv, g1);

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
    _LA.solve_LU_fast(A, false, ipiv, g2);

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
    _LA.solve_LU_fast(A, false, ipiv, g3);

    // setup r.h.s. for g4
    for (unsigned i=0; i< n; i++)
      g4[i] = dxdt[i] + h*C4X*dfdt[i] + (C41*g1[i] + C42*g2[i] + C43*g3[i])/h;

    // solve for g4
    _LA.solve_LU_fast(A, false, ipiv, g4);

    // get fourth-order estimate of x and error estimate
    VectorNd err(n);
    for (unsigned i=0; i< n; i++)
    {
      x[i] = xsav[i] + B1*g1[i] + B2*g2[i] + B3*g3[i] + B4*g4[i];
      err[i] = E1*g1[i] + E2*g2[i] + E3*g3[i] + E4*g4[i];
    }

    // update t
    time = tsav + h;

    // evaluate accuracy
    double errmax = (double) 0.0;
    for (unsigned i=0; i< n; i++)
      errmax = std::max(errmax, std::fabs(err[i]/xscal[i]));
    errmax /= rel_err_tol; 
    if (errmax <= (double) 1.0)
    {
      // step succeeded; compute size of next step and return
      tnext = (errmax > ERRCON ? SAFETY*h*std::pow(errmax,PGROW) : GROW*h);
      return;
    }
    else
    {
      // truncation error too large, reduce step size
      tnext = SAFETY*h*std::pow(errmax, PSHRINK);
      h = (h >= (double) 0.0 ? std::max(tnext, SHRINK*h) : std::min(tnext, SHRINK*h));
    }
  }

  std::cerr << "Rosenbrock4Integrator::step() - maximum number of iterations exceeded (40)" << endl;
}

void Rosenbrock4Integrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "Rosenbrock4Integrator") == 0);
  Integrator::load_from_xml(node, id_map); 
  
  // get the relative error tolerance
  XMLAttrib* rerr_attrib = node->get_attrib("rel-err-tol");
  if (rerr_attrib)
    rel_err_tol = rerr_attrib->get_real_value();
}

void Rosenbrock4Integrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects); 
  node->name = "Rosenbrock4Integrator";
  node->attribs.insert(XMLAttrib("rel-err-tol", rel_err_tol));
}


