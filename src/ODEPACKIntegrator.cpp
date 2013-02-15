/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/ODEPACKIntegrator.h>

#ifdef THREADSAFE
pthread_mutex_t Moby::ODEPACKIntegratorMutex::_odepack_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

namespace Moby {

/// Method for integration
template <>
void ODEPACKIntegrator<Vector3>::integrate(Vector3& x, Vector3& (*f)(const Vector3&, Real, Real, void*, Vector3&), Real& time, Real step_size, void* data)
{
  #ifndef USE_ODEPACK
  std::cerr << "Moby built without ODEPACK support!  Using explicit Euler integration instead..." << std::endl;

  // save the old time
  const Real old_time = time;

  // update the time
  time += step_size;

  // compute new state
  Vector3 dx;
  x += (f(x, old_time, step_size, data, dx) * step_size);
  #else

  // get the maximum number of steps
  const unsigned MAX_STEPS = 100000;

  // setup non-precision specific inputs
  int neq = (int) x.size();
  int itol = 1;   // indicate atol is a scalar
  int itask = 1;  // normal computation of output
  int istate = 1; // recompute everything
  int iopt = 1;   // use optional arguments
  int lrw = (use_stiff_integrator) ? (22+9*neq+neq*neq) : (20+16*neq);
  int liw = (use_stiff_integrator) ? (20+neq) : 20;
  int mf = (use_stiff_integrator) ? 22 : 10;

  // setup integer work array
  boost::shared_array<int> iwork(new int[liw]);

  #ifdef BUILD_SINGLE

  // setup float parameters
  float tout = time + step_size;
  float atol = std::numeric_limits<float>::epsilon();

  // setup float work array
  boost::shared_array<float> rwork(new float[lrw]);

  // setup optional inputs (0 indicates default)
  rwork[4] = 0f;            // step size attempted on the first step
  rwork[5] = 0f;            // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Vector3 = f;
  fn_data = data;  

  // call ODEPACK - note that it updates t
  slsode_(fsingledouble_Vector3, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          rwork.get(), &lrw, iwork.get(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  #else // BUILD_DOUBLE or BUILD_ARBITRARY_PRECISION

  // setup double parameters
  double tout = time + step_size;

  // setup double work array
  boost::shared_array<double> rwork(new double[lrw]);

  // setup optional inputs (0 indicates default)
  rwork[4] = (Real) 0.0;    // step size attempted on the first step
  rwork[5] = (Real) 0.0;    // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  #ifdef BUILD_DOUBLE

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Vector3 = f;
  fn_data = data;  

  // call ODEPACK -- note that it updates t
  dlsode_(fsingledouble_Vector3, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          rwork.get(), &lrw, iwork.get(), &liw, NULL, &mf);

  #else // BUILD_ARBITRARY_PRECISION

  // setup t
  double t = (double) time;

  // setup double version of x
  boost::shared_array<double> xdouble(new double[neq]);
  for (int i=0; i< neq; i++)
    xdouble[i] = x[i];

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Vector3 = f;
  fn_data = data;  

  // call ODEPACK -- note that it updates t
  dlsode_(fsingledouble_Vector3, &neq, xdouble, &t, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          rwork.get(), &lrw, iwork.get(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  // get new time
  time = (Real) t;

  // get x
  for (int i=0; i< neq; i++)
    x[i] = xdouble[i];

  #endif // BUILD DOUBLE / BUILD_ARBITRARY_PRECISION
  #endif // BUILD_SINGLE

  // output appropriate error messages, based on istate
  switch (istate)
  {
    case 1:
      std::cerr << "ODEPACKIntegrator::integrate(): nothing was done b/c t=tout" << std::endl;
      break;

    case 2:
      FILE_LOG(LOG_SIMULATOR) << "ODEPACKIntegrator::integrate(): ODEPACK reports success" << std::endl;
      break;

    case -1:
      std::cerr << "ODEPACKIntegrator::integrate(): excess work done on this call" << std::endl;
      break;

    case -2:
      std::cerr << "ODEPACKIntegrator::integrate(): excess accuracy requested (tolerances too small)" << std::endl;
      break;

    case -3:
      std::cerr << "ODEPACKIntegrator::integrate(): illegal inputs detected" << std::endl;
      break;

    case -4:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated error test failures (check all inputs)" << std::endl;
      break;

    case -5:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated convergence failures (perhaps wrong " << std::endl << "choice of integration method or tolerances)." << std::endl;
      break;

    case -6:
      std::cerr << "ODEPACKIntegrator::integrate(): error weight became zero and pure relative " << std::endl << "error control (absolute error tolerance = 0.0) used" << std::endl;
      break;
  }

  #endif // USE_ODEPACK
}

/// Method for integration
template <>
void ODEPACKIntegrator<VectorN>::integrate(VectorN& x, VectorN& (*f)(const VectorN&, Real, Real, void*, VectorN&), Real& time, Real step_size, void* data)
{
  #ifndef USE_ODEPACK
  std::cerr << "Moby built without ODEPACK support!  Using explicit Euler integration instead..." << std::endl;

  // save the old time
  const Real old_time = time;

  // update the time
  time += step_size;

  // compute new state
  VectorN dx;
  x += (f(x, old_time, step_size, data, dx) * step_size);
  #else

  // get the maximum number of steps
  const unsigned MAX_STEPS = 100000;

  // setup non-precision specific inputs
  int neq = (int) x.size();
  int itol = 1;  // indicate atol is a scalar
  int itask = 1; // normal computation of output
  int istate = 1; // recompute everything
  int iopt = 1; // use options
  int lrw = (use_stiff_integrator) ? (22+9*neq+neq*neq) : (20+16*neq);
  int liw = (use_stiff_integrator) ? (20+neq) : 20;
  int mf = (use_stiff_integrator) ? 22 : 10;

  // setup integer work array
  SAFESTATIC std::vector<int> iwork;
  iwork.resize(liw);

  #ifdef BUILD_SINGLE

  // setup float parameters
  float tout = time + step_size;

  // setup float work array
  SAFESTATIC std::vector<float> rwork;
  rwork.resize(lrw);

  // setup optional inputs (0 indicates default)
  rwork[4] = 0f;            // step size attempted on the first step
  rwork[5] = 0f;            // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_VectorN = f;
  fn_data = data;

  // call ODEPACK - note that it updates t
  slsode_(fsingledouble_VectorN, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          &rwork.front(), &lrw, &iwork.front(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  #else // BUILD_DOUBLE or BUILD_ARBITRARY_PRECISION

  // setup double parameters
  double tout = time + step_size;

  // setup double work array
  std::vector<double> rwork;
  rwork.resize(lrw);

  // setup optional inputs (0 indicates default)
  rwork[4] = (Real) 0.0;    // step size attempted on the first step
  rwork[5] = (Real) 0.0;    // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  #ifdef BUILD_DOUBLE

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_VectorN = f;
  fn_data = data;

  // call ODEPACK -- note that it updates t
  dlsode_(fsingledouble_VectorN, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          &rwork.front(), &lrw, &iwork.front(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  #else // BUILD_ARBITRARY_PRECISION

  // setup t
  double t = (double) time;

  // setup double version of x
  boost::shared_array<double> xdouble(new double[neq]);
  for (int i=0; i< neq; i++)
    xdouble[i] = x[i];

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_VectorN = f;
  fn_data = data;

  // call ODEPACK -- note that it updates t
  dlsode_(farb_VectorN, &neq, xdouble, &t, &tout, &itol, &rerr_tolerance,
          &aerr_tolerance, &itask, &istate, &iopt, &rwork.front(), &lrw, 
          &iwork.front(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  // get new time
  time = (Real) t;

  // get x
  for (int i=0; i< neq; i++)
    x[i] = xdouble[i];

  #endif // BUILD DOUBLE / BUILD_ARBITRARY_PRECISION
  #endif // BUILD_SINGLE

  // output appropriate error messages, based on istate
  switch (istate)
  {
    case 1:
      std::cerr << "ODEPACKIntegrator::integrate(): nothing was done b/c t=tout" << std::endl;
      break;

    case 2:
      FILE_LOG(LOG_SIMULATOR) << "ODEPACKIntegrator::integrate(): ODEPACK reports success" << std::endl;
      break;

    case -1:
      std::cerr << "ODEPACKIntegrator::integrate(): excess work done on this call" << std::endl;
      break;

    case -2:
      std::cerr << "ODEPACKIntegrator::integrate(): excess accuracy requested (tolerances too small)" << std::endl;
      break;

    case -3:
      std::cerr << "ODEPACKIntegrator::integrate(): illegal inputs detected" << std::endl;
      break;

    case -4:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated error test failures (check all inputs)" << std::endl;
      break;

    case -5:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated convergence failures (perhaps wrong " << std::endl << "choice of integration method or tolerances)." << std::endl;
      break;

    case -6:
      std::cerr << "ODEPACKIntegrator::integrate(): error weight became zero and pure relative " << std::endl << "error control (absolute error tolerance = 0.0) used" << std::endl;
      break;
  }

  #endif // USE_ODEPACK
}

/// Method for integration
template <>
void ODEPACKIntegrator<Quat>::integrate(Quat& x, Quat& (*f)(const Quat&, Real, Real, void*, Quat&), Real& time, Real step_size, void* data)
{
  #ifndef USE_ODEPACK
  std::cerr << "Moby built without ODEPACK support!  Using explicit Euler integration instead..." << std::endl;

  // save the old time
  const Real old_time = time;

  // update the time
  time += step_size;

  // compute new state
  Quat dx;
  x += (f(x, old_time, step_size, data, dx) * step_size);
  #else

  // get the maximum number of steps
  const unsigned MAX_STEPS = 100000;

  // setup non-precision specific inputs
  int neq = (int) x.size();
  int itol = 1;   // indicate atol is a scalar
  int itask = 1;  // normal computation of output
  int istate = 1; // recompute everything
  int iopt = 1;   // use options
  int lrw = (use_stiff_integrator) ? (22+9*neq+neq*neq) : (20+16*neq);
  int liw = (use_stiff_integrator) ? (20+neq) : 20;
  int mf = (use_stiff_integrator) ? 22 : 10;

  // setup integer work array
  boost::shared_array<int> iwork(new int[liw]);

  #ifdef BUILD_SINGLE

  // setup float parameters
  float tout = time + step_size;

  // setup float work array
  boost::shared_array<float> rwork(new float[lrw]);

  // setup optional inputs (0 indicates default)
  rwork[4] = 0f;            // step size attempted on the first step
  rwork[5] = 0f;            // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Quat = f;
  fn_data = data;  

  // call ODEPACK - note that it updates t
  slsode_(fsingledouble_Quat, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          rwork.get(), &lrw, iwork.get(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  #else // BUILD_DOUBLE or BUILD_ARBITRARY_PRECISION

  // setup double parameters
  double tout = time + step_size;

  // setup double work array
  boost::shared_array<double> rwork(new double[lrw]);

  // setup optional inputs (0 indicates default)
  rwork[4] = (Real) 0.0;    // step size attempted on the first step
  rwork[5] = (Real) 0.0;    // maximum absolute step size allowed
  rwork[6] = min_step_size; // minimum absolute step size allowed 
  iwork[4] = 0;             // maximum order allowed
  iwork[5] = MAX_STEPS;     // maximum number of internal steps
  iwork[6] = 1;             // maximum number of messages printed

  #ifdef BUILD_DOUBLE

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Quat = f;
  fn_data = data;  

  // call ODEPACK -- note that it updates t
  dlsode_(fsingledouble_Quat, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          rwork.get(), &lrw, iwork.get(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  #else // BUILD_ARBITRARY_PRECISION

  // setup t
  double t = (double) time;

  // setup double version of x
  boost::shared_array<double> xdouble(new double[neq]);
  for (int i=0; i< neq; i++)
    xdouble[i] = x[i];

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_Quat = f;
  fn_data = data;  

  // call ODEPACK -- note that it updates t
  dlsode_(farb_Quat, &neq, xdouble, &t, &tout, &itol, &rerr_tolerance,
          &aerr_tolerance, &itask, &istate, &iopt, rwork.get(), &lrw, 
          iwork.get(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

  // get new time
  time = (Real) t;

  // get x
  for (int i=0; i< neq; i++)
    x[i] = xdouble[i];

  #endif // BUILD DOUBLE / BUILD_ARBITRARY_PRECISION
  #endif // BUILD_SINGLE

  // output appropriate error messages, based on istate
  switch (istate)
  {
    case 1:
      std::cerr << "ODEPACKIntegrator::integrate(): nothing was done b/c t=tout" << std::endl;
      break;

    case 2:
      FILE_LOG(LOG_SIMULATOR) << "ODEPACKIntegrator::integrate(): ODEPACK reports success" << std::endl;
      break;

    case -1:
      std::cerr << "ODEPACKIntegrator::integrate(): excess work done on this call" << std::endl;
      break;

    case -2:
      std::cerr << "ODEPACKIntegrator::integrate(): excess accuracy requested (tolerances too small)" << std::endl;
      break;

    case -3:
      std::cerr << "ODEPACKIntegrator::integrate(): illegal inputs detected" << std::endl;
      break;

    case -4:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated error test failures (check all inputs)" << std::endl;
      break;

    case -5:
      std::cerr << "ODEPACKIntegrator::integrate(): repeated convergence failures (perhaps wrong " << std::endl << "choice of integration method or tolerances)." << std::endl;
      break;

    case -6:
      std::cerr << "ODEPACKIntegrator::integrate(): error weight became zero and pure relative " << std::endl << "error control (absolute error tolerance = 0.0) used" << std::endl;
      break;
  }

  #endif // USE_ODEPACK
}

} // end namespace

