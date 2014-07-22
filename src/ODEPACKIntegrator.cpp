/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <algorithm>
#include <strings.h>
#include <Moby/ODEPACKIntegrator.h>

using boost::shared_ptr;
using namespace Moby;
using Ravelin::VectorNd;

#ifdef THREADSAFE
pthread_mutex_t Moby::ODEPACKIntegratorMutex::_odepack_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef USE_ODEPACK
extern "C"  {

// prototype for ODEPACK (double version)
void dlsode_(void (*f)(int*, double*, double*, double*), int* neq, double* y,
            double* t, double* tout, int* itol, double* rtol, double* atol,
            int* itask, int* istate, int* iopt, double* rwork, int* lrw,
            int* iwork, int* liw, 
            void (*jac)(int*, double*, double*, int*, int*, double*, int*),
            int* mf);

} // extern "C" 

static double tnew;
static void* fn_data;
static VectorNd& (*fn_VectorNd)(const VectorNd&, double, double, void*, VectorNd&);

static void fdouble_VectorNd(int* neq, double* t, double* y, double* ydot)
{
  static VectorNd x, xdot;
  x.resize(*neq);
  std::copy(y, y+*neq, x.data());
  (*fn_VectorNd)(x, *t, tnew - *t, fn_data, xdot);
  std::copy(xdot.data(), xdot.data()+*neq, ydot);
}
#endif // #ifdef USE_ODEPACK

/// Method for integration
void ODEPACKIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  #ifndef USE_ODEPACK
  std::cerr << "Moby built without ODEPACK support!  Using explicit Euler integration instead..." << std::endl;

  // save the old time
  const double old_time = time;

  // update the time
  time += step_size;

  // compute new state
  VectorNd dx;
  x += (f(x, old_time, step_size, data, dx) *= step_size);
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
  _iwork.resize(liw);

  // setup double parameters
  double tout = time + step_size;

  // setup double work array
  _rwork.resize(lrw);

  // setup optional inputs (0 indicates default)
  _rwork[4] = (double) 0.0;    // step size attempted on the first step
  _rwork[5] = (double) 0.0;    // maximum absolute step size allowed
  _rwork[6] = min_step_size; // minimum absolute step size allowed 
  _iwork[4] = 0;             // maximum order allowed
  _iwork[5] = MAX_STEPS;     // maximum number of internal steps
  _iwork[6] = 1;             // maximum number of messages printed

  // save the derivative function and data
  #ifdef THREADSAFE
  pthread_mutex_lock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif
  tnew = time + step_size;
  fn_VectorNd = f;
  fn_data = data;

  // call ODEPACK -- note that it updates t
  dlsode_(fdouble_VectorNd, &neq, x.data(), &time, &tout, &itol, 
          &rerr_tolerance, &aerr_tolerance, &itask, &istate, &iopt, 
          &_rwork.front(), &lrw, &_iwork.front(), &liw, NULL, &mf);

  // release the mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&ODEPACKIntegratorMutex::_odepack_mutex);
  #endif

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


ODEPACKIntegrator::ODEPACKIntegrator()
{
  use_stiff_integrator = false;
}

/// Implements Base::load_from_xml()
void ODEPACKIntegrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "ODEPACKIntegrator") == 0);
  VariableStepIntegrator::load_from_xml(node, id_map); 

  // get whether the stiff integrator should be used
  XMLAttrib* stiff_attrib = node->get_attrib("stiff");
  if (stiff_attrib)
    use_stiff_integrator = stiff_attrib->get_bool_value();
}

/// Implements Base::save_to_xml()
void ODEPACKIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{ 
  VariableStepIntegrator::save_to_xml(node, shared_objects); 
  node->name = "ODEPACKIntegrator";

  // save the method data
  node->attribs.insert(XMLAttrib("stiff", use_stiff_integrator));
}

