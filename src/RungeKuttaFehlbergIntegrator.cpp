/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/RungeKuttaFehlbergIntegrator.h>

using boost::shared_ptr;
using Ravelin::VectorNd;
using namespace Moby;

RungeKuttaFehlbergIntegrator::RungeKuttaFehlbergIntegrator()
{
}

void RungeKuttaFehlbergIntegrator::fehl ( void f ( double t, double* y, double* yp, void* data), int neqn, double* y, double t, double h, double* yp, double* f1, double* f2, double* f3, double* f4, double* f5, double* s)
//****************************************************************************80
//
//  Purpose:
//
//    R8_FEHL takes one Fehlberg fourth-fifth order step.
//
//  Discussion:
//
//    This version of the routine uses DOUBLE real arithemtic.
//
//    This routine integrates a system of NEQN first order ordinary differential
//    equations of the form
//      dY(i)/dT = F(T,Y(1:NEQN))
//    where the initial values Y and the initial derivatives
//    YP are specified at the starting point T.
//
//    The routine advances the solution over the fixed step H and returns
//    the fifth order (sixth order accurate locally) solution
//    approximation at T+H in array S.
//
//    The formulas have been grouped to control loss of significance.
//    The routine should be called with an H not smaller than 13 units of
//    roundoff in T so that the various independent arguments can be
//    distinguished.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Erwin Fehlberg,
//    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
//    NASA Technical Report R-315, 1969.
//
//    Lawrence Shampine, Herman Watts, S Davenport,
//    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
//    SIAM Review,
//    Volume 18, pages 376-411, 1976.
//
//  Parameters:
//
//    Input, external F, a user-supplied subroutine to evaluate the
//    derivatives Y'(T), of the form:
//
//      void f ( double t, double y[], double yp[] )
//
//    Input, int NEQN, the number of equations to be integrated.
//
//    Input, double Y[NEQN], the current value of the dependent variable.
//
//    Input, double T, the current value of the independent variable.
//
//    Input, double H, the step size to take.
//
//    Input, double YP[NEQN], the current value of the derivative of the
//    dependent variable.
//
//    Output, double F1[NEQN], F2[NEQN], F3[NEQN], F4[NEQN], F5[NEQN], derivative
//    values needed for the computation.
//
//    Output, double S[NEQN], the estimate of the solution at T+H.
//
{
  double ch;

  ch = h / (double) 4.0;

  for ( int i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * yp[i];
  }

  f ( t + ch, f5, f1, this );

  ch = (double) 3.0 * h / (double) 32.0;

  for ( int i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * ( yp[i] + (double) 3.0 * f1[i] );
  }

  f ( t + (double) 3.0 * h / (double) 8.0, f5, f2, this );

  ch = h / (double) 2197.0;

  for (int i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( (double) 1932.0 * yp[i] 
    + ( (double) 7296.0 * f2[i] - (double) 7200.0 * f1[i] ) 
    );
  }

  f ( t + (double) 12.0 * h / (double) 13.0, f5, f3, this );

  ch = h / (double) 4104.0;

  for (int i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( 
      ( (double) 8341.0 * yp[i] - (double) 845.0 * f3[i] ) 
    + ( (double) 29440.0 * f2[i] - (double) 32832.0 * f1[i] ) 
    );
  }

  f ( t + h, f5, f4, this );

  ch = h / (double) 20520.0;

  for (int i = 0; i < neqn; i++ )
  {
    f1[i] = y[i] + ch * 
    ( 
      ( (double) -6080.0 * yp[i] 
      + ( (double) 9295.0 * f3[i] - (double) 5643.0 * f4[i] ) 
      ) 
    + ( (double) 41040.0 * f1[i] - (double) 28352.0 * f2[i] ) 
    );
  }

  f ( t + h / (double) 2.0, f1, f5, this );
//
//  Ready to compute the approximate solution at T+H.
//
  ch = h / (double) 7618050.0;

  for (int i = 0; i < neqn; i++ )
  {
    s[i] = y[i] + ch * 
    ( 
      ( (double) 902880.0 * yp[i] 
      + ( (double) 3855735.0 * f3[i] - (double) 1371249.0 * f4[i] ) ) 
    + ( (double) 3953664.0 * f2[i] + (double) 277020.0 * f5[i] ) 
    );
  }

  return;
}

void RungeKuttaFehlbergIntegrator::rkf45 ( void f ( double t, double* y, double* yp, void* data), int neqn, double* y, double* yp, double& t, double tout, bool init)

//****************************************************************************80
//
//  Purpose:
//
//    RKF45 carries out the Runge-Kutta-Fehlberg method.
//
//  Discussion:
//
//    This routine is primarily designed to solve non-stiff and mildly stiff
//    differential equations when derivative evaluations are inexpensive.
//    It should generally not be used when the user is demanding
//    high accuracy.
//
//    This routine integrates a system of NEQN first-order ordinary differential
//    equations of the form:
//
//      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
//
//    where the Y(1:NEQN) are given at T.
//
//    Typically the subroutine is used to integrate from T to TOUT but it
//    can be used as a one-step integrator to advance the solution a
//    single step in the direction of TOUT.  On return, the parameters in
//    the call list are set for continuing the integration.  The user has
//    only to call again (and perhaps define a new value for TOUT).
//
//    Before the first call, the user must 
//
//    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
//      and declare F in an EXTERNAL statement;
//
//    * initialize the parameters:
//      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
//      In particular, T should initially be the starting point for integration,
//      Y should be the value of the initial conditions, and FLAG should 
//      normally be +1.
//
//    Normally, the user only sets the value of FLAG before the first call, and
//    thereafter, the program manages the value.  On the first call, FLAG should
//    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
//    have been reset by the program to the value of 2 (or -2 in single 
//    step mode), and the user can continue to call the routine with that 
//    value of FLAG.
//
//    (When the input magnitude of FLAG is 1, this indicates to the program 
//    that it is necessary to do some initialization work.  An input magnitude
//    of 2 lets the program know that that initialization can be skipped, 
//    and that useful information was computed earlier.)
//
//    The routine returns with all the information needed to continue
//    the integration.  If the integration reached TOUT, the user need only
//    define a new TOUT and call again.  In the one-step integrator
//    mode, returning with FLAG = -2, the user must keep in mind that 
//    each step taken is in the direction of the current TOUT.  Upon 
//    reaching TOUT, indicated by the output value of FLAG switching to 2,
//    the user must define a new TOUT and reset FLAG to -2 to continue 
//    in the one-step integrator mode.
//
//    In some cases, an error or difficulty occurs during a call.  In that case,
//    the output value of FLAG is used to indicate that there is a problem
//    that the user must address.  These values include:
//
//    * 3, integration was not completed because the input value of RELERR, the 
//      relative error tolerance, was too small.  RELERR has been increased 
//      appropriately for continuing.  If the user accepts the output value of
//      RELERR, then simply reset FLAG to 2 and continue.
//      (NOTE: this is now disabled)
//
//    * 4, integration was not completed because more than MAXNFE derivative 
//      evaluations were needed.  This is approximately (MAXNFE/6) steps.
//      The user may continue by simply calling again.  The function counter 
//      will be reset to 0, and another MAXNFE function evaluations are allowed.
//      (NOTE: this is now disabled)
//
//    * 5, integration was not completed because the solution vanished, 
//      making a pure relative error test impossible.  The user must use 
//      a non-zero ABSERR to continue.  Using the one-step integration mode 
//      for one step is a good way to proceed.
//
//    * 6, integration was not completed because the requested accuracy 
//      could not be achieved, even using the smallest allowable stepsize. 
//      The user must increase the error tolerances ABSERR or RELERR before
//      continuing.  It is also necessary to reset FLAG to 2 (or -2 when 
//      the one-step integration mode is being used).  The occurrence of 
//      FLAG = 6 indicates a trouble spot.  The solution is changing 
//      rapidly, or a singularity may be present.  It often is inadvisable 
//      to continue.
//
//    * 7, it is likely that this routine is inefficient for solving
//      this problem.  Too much output is restricting the natural stepsize
//      choice.  The user should use the one-step integration mode with 
//      the stepsize determined by the code.  If the user insists upon 
//      continuing the integration, reset FLAG to 2 before calling 
//      again.  Otherwise, execution will be terminated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Erwin Fehlberg,
//    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
//    NASA Technical Report R-315, 1969.
//
//    Lawrence Shampine, Herman Watts, S Davenport,
//    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
//    SIAM Review,
//    Volume 18, pages 376-411, 1976.
//
//  Parameters:
//
//    Input, external F, a user-supplied subroutine to evaluate the
//    derivatives Y'(T), of the form:
//
//      void f ( double t, double y[], double yp[] )
//
//    Input, int NEQN, the number of equations to be integrated.
//
//    Input/output, double Y[NEQN], the current solution vector at T.
//
//    Input/output, double YP[NEQN], the derivative of the current solution 
//    vector at T.  The user should not set or alter this information!
//
//    Input/output, double T, the current value of the independent variable.
//
//    Input, double TOUT, the output point at which solution is desired.  
//    TOUT = T is allowed on the first call only, in which case the routine
//    returns with FLAG = 2 if continuation is possible.
//
//    Input, double RELERR, ABSERR, the relative and absolute error tolerances
//    for the local error test.  At each step the code requires:
//      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
//    for each component of the local error and the solution vector Y.
//    RELERR cannot be "too small".  If the routine believes RELERR has been
//    set too small, it will reset RELERR to an acceptable value and return
//    immediately for user action.
//
//    Input, int FLAG, indicator for status of integration. On the first call, 
//    set FLAG to +1 for normal use, or to -1 for single step mode.  On 
//    subsequent continuation steps, FLAG should be +2, or -2 for single 
//    step mode.
//
//    Output, int RKF45_D, indicator for status of integration.  A value of 2 
//    or -2 indicates normal progress, while any other value indicates a 
//    problem that should be addressed.
//
{
  double esttol = (double) -1.0;
  double et;
  double s;
  double toln;
//
//  Check the input parameters.
//
  const double EPS = std::numeric_limits<double>::epsilon();

  if ( neqn < 1 )
    throw std::runtime_error("At least one equation required for Runge-Kutta-Fehlberg integration");

  if (this->rerr_tolerance < (double) 0.0 )
    throw std::runtime_error("Relative error tolerance must be non-negative");

  if (this->aerr_tolerance < (double) 0.0 )
    throw std::runtime_error("Absolute error tolerance must be non-negative");

  // resize f1..f5, if necessary
  if (f1.size() != neqn)
  {
    f1.resize(neqn);
    f2.resize(neqn);
    f3.resize(neqn);
    f4.resize(neqn);
    f5.resize(neqn);
  }

//  Is this a continuation call?
  if (!init && t == tout)
      return;

//
//  Restrict the relative error tolerance to be at least 
//
//    2*EPS+REMIN 
//
//  to avoid limiting precision difficulties arising from impossible 
//  accuracy requests.
//
  const double REMIN = 2.0 * std::numeric_limits<double>::epsilon() + 1e-12;
//
//  Is the relative error tolerance too small?
//
  if (this->rerr_tolerance < REMIN)
    std::cerr << "RungeKuttaFehlbergIntegrator warning- relative tolerance too small; " << std::endl << "using " << REMIN << " instead" << std::endl;

  // setup relerr
  const double RELERR = std::max(this->rerr_tolerance, REMIN); 

  // setup the step size
  double dt = tout - t;
//
//  Initialization:
//
//  Set the initialization completion indicator, INIT;
//  set the indicator for too many output points, KOP;
//  evaluate the initial derivatives
//  set the counter for function evaluations, NFE;
//  estimate the starting stepsize.
//

  if (init)
  {
    _kop = 0;
    f ( t, y, yp, this );

    if (t == tout )
      return;

    _h = std::fabs( dt );
    toln = (double) 0.0;

    for (int k = 0; k < neqn; k++ )
    {
      double tol = RELERR * std::fabs ( y[k] ) + this->aerr_tolerance;
      if ( (double) 0.0 < tol )
      {
        toln = tol;
        double ypk = std::fabs ( yp[k] );
        if ( tol < ypk * std::pow ( _h, (double) 5 ) )
        {
          _h = std::pow ( ( tol / ypk ), (double) 0.2 );
        }
      }
    }

    if ( toln <= (double) 0.0 )
    {
      _h = (double) 0.0;
    }

    _h = std::max ( _h, (double) 26.0 * EPS * std::max ( std::fabs ( t ), std::fabs ( dt ) ) );
  }
//
//  Set stepsize for integration in the direction from T to TOUT.
//
  _h = sign ( dt ) * std::fabs ( _h );
//
//  Test to see if too may output points are being requested.
//
  if ( (double) 2.0 * std::fabs ( dt ) <= std::fabs ( _h ) )
  {
    _kop = _kop + 1;
  }
//
//  Unnecessary frequency of output.
//
  if ( _kop == 100 )
  {
    _kop = 0;
    std::cerr << "RungeKuttaFehlbergIntegrator::integrate()- function troublesome for integrator!" << std::endl;
    return;
  }
//
//  If we are too close to the output point, then simply extrapolate and return.
//
  if ( std::fabs ( dt ) <= (double) 26.0 * EPS * std::fabs ( t ) )
  {
    t = tout;
    for (int i = 0; i < neqn; i++ )
    {
      y[i] = y[i] + dt * yp[i];
    }
    f ( t, y, yp, this );

    // normal return
    return;
  }
//
//  Initialize the output point indicator.
//
  bool output = false;
//
//  To avoid premature underflow in the error tolerance function,
//  scale the error tolerances.
//
  const double SCALE = (double) 2.0 / RELERR;
  const double AE = SCALE * this->aerr_tolerance;
//
//  Step by step integration.
//
  for ( ; ; )
  {
    bool hfaild = false;
//
//  Set the smallest allowable stepsize.
//
    double hmin = (double) 26.0 * EPS * std::fabs ( t );
//
//  Adjust the stepsize if necessary to hit the output point.
//
//  Look ahead two steps to avoid drastic changes in the stepsize and
//  thus lessen the impact of output points on the code.
//
    dt = tout - t;

    if ( (double) 2.0 * std::fabs ( _h ) <= std::fabs ( dt ) )
    {
    }
    else
//
//  Will the next successful step complete the integration to the output point?
//
    {
      if ( std::fabs ( dt ) <= std::fabs ( _h ) )
      {
        output = true;
        _h = dt;
      }
      else
      {
        _h = (double) 0.5 * dt;
      }

    }
//
//  Here begins the core integrator for taking a single step.
//
//  The tolerances have been scaled to avoid premature underflow in
//  computing the error tolerance function ET.
//  To avoid problems with zero crossings, relative error is measured
//  using the average of the magnitudes of the solution at the
//  beginning and end of a step.
//  The error estimate formula has been grouped to control loss of
//  significance.
//
//  To distinguish the various arguments, H is not permitted
//  to become smaller than 26 units of roundoff in T.
//  Practical limits on the change in the stepsize are enforced to
//  smooth the stepsize selection process and to avoid excessive
//  chattering on problems having discontinuities.
//  To prevent unnecessary failures, the code uses 9/10 the stepsize
//  it estimates will succeed.
//
//  After a step failure, the stepsize is not allowed to increase for
//  the next attempted step.  This makes the code more efficient on
//  problems having discontinuities and more effective in general
//  since local extrapolation is being used and extra caution seems
//  warranted.
//
//  Test the number of derivative function evaluations.
//  If okay, try to advance the integration from T to T+H.
//
    for ( ; ; )
    {
//
//  Advance an approximate solution over one step of length H.
//
      fehl ( f, neqn, y, t, _h, yp, &f1[0], &f2[0], &f3[0], &f4[0], &f5[0], &f1[0] );
//
//  Compute and test allowable tolerances versus local error estimates
//  and remove scaling of tolerances.  The relative error is
//  measured with respect to the average of the magnitudes of the
//  solution at the beginning and end of the step.
//
      double eeoet = (double) 0.0;
 
      for (int k = 0; k < neqn; k++ )
      {
        const double ET = std::fabs ( y[k] ) + std::fabs ( f1[k] ) + AE;

        // solution vanished: pure relative error test impossible
        if ( ET <= (double) 0.0 )
        {
          std::cerr << "RungeKuttaFehlbergIntegrator::integrate() - solution vanished;" << std::endl << "pure relative error test impossible" << std::endl;
          return;
        }

        const double EE = std::fabs 
        ( ((double)  -2090.0 * yp[k] 
          + ((double)  21970.0 * f3[k] - (double) 15048.0 * f4[k] ) 
          ) 
        + ( (double) 22528.0 * f2[k] - (double) 27360.0 * f5[k] ) 
        );

        eeoet = std::max ( eeoet, EE / ET );

      }

      esttol = std::fabs ( _h ) * eeoet * SCALE / (double) 752400.0;

      if ( esttol <= (double) 1.0 )
      {
        break;
      }
//
//  Unsuccessful step.  Reduce the stepsize, try again.
//  The decrease is limited to a factor of 1/10.
//
      hfaild = true;
      output = false;

      if ( esttol < (double) 59049.0 )
      {
        s = (double) 0.9 / std::pow ( esttol, (double) 0.2 );
      }
      else
      {
        s = (double) 0.1;
      }

      _h = s * _h;

      if ( std::fabs (_h ) < hmin )
      {
        // step size too small
        std::cerr << "RungeKuttaFehlbergIntegrator::integrate()- step size " << _h << " smaller than minimum (" << hmin << ")" << std::endl;
        return;
      }

    }
//
//  We exited the loop because we took a successful step.  
//  Store the solution for T+H, and evaluate the derivative there.
//
    t = t + _h;
    for (int i = 0; i < neqn; i++ )
    {
      y[i] = f1[i];
    }
    f ( t, y, yp, this );
//
//  Choose the next stepsize.  The increase is limited to a factor of 5.
//  If the step failed, the next stepsize is not allowed to increase.
//
    if ( (double) 0.0001889568 < esttol )
    {
      s = (double) 0.9 / std::pow ( esttol, (double) 0.2 );
    }
    else
    {
      s = (double) 5.0;
    }

    if ( hfaild )
    {
      s = std::min ( s, (double) 1.0 );
    }

    _h = sign ( _h ) * std::max ( s * std::fabs ( _h ), hmin );
//
//  End of core integrator
//
//  Should we take another step?
//
    if ( output )
    {
      // normal return
      t = tout;
      return;
    }
  }
}

double RungeKuttaFehlbergIntegrator::sign ( double x )
{
  if ( x < (double) 0.0 )
  {
    return ( (double) -1.0 );
  } 
  else
  {
    return ( (double) +1.0 );
  }
}

void RungeKuttaFehlbergIntegrator::fprime(double t, double* y, double* yprime, void* data)
{
  // get the object
  RungeKuttaFehlbergIntegrator* rkf = (RungeKuttaFehlbergIntegrator*) data;

  // get the derivative
  std::copy(y, y+rkf->_N, rkf->_x.data());
  std::copy(yprime, yprime+rkf->_N, rkf->_xprime.data());
  rkf->_f(rkf->_x, t, rkf->_dt, rkf->_data, rkf->_xprime);
}

void RungeKuttaFehlbergIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  const double desired_time = time + step_size;

  // setup data
  _N = x.size();
  _x.resize(_N);
  _xprime.resize(_N);
  _dt = step_size;
  _f = f;
  _data = data;

  // get derivative at t
  _xprime.resize(x.size());
  f(x, time, step_size, data, _xprime);

  // initialize first
  rkf45(fprime, x.size(), x.data(), _xprime.data(), time, desired_time, true);

  // continue running
  while (time < desired_time)
  {
    _dt = desired_time - time;
    f(x, time, _dt, data, _xprime);
    rkf45(fprime, x.size(), x.data(), _xprime.data(), time, desired_time, false);
  }
}

void RungeKuttaFehlbergIntegrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "RungeKuttaFehlbergIntegrator") == 0);

  // call the parent method
  VariableStepIntegrator::load_from_xml(node, id_map);
}

void RungeKuttaFehlbergIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent method 
  VariableStepIntegrator::save_to_xml(node, shared_objects); 

  // rename the node
  node->name = "RungeKuttaFehlbergIntegrator";
}


