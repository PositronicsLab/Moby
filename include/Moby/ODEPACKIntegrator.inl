/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef USE_ODEPACK
extern "C"  {

// prototype for ODEPACK (double version)
void dlsode_(void (*f)(int*, double*, double*, double*), int* neq, double* y,
            double* t, double* tout, int* itol, double* rtol, double* atol,
            int* itask, int* istate, int* iopt, double* rwork, int* lrw,
            int* iwork, int* liw, 
            void (*jac)(int*, double*, double*, int*, int*, double*, int*),
            int* mf);

// prototype for ODEPACK (float version)
void slsode_(void (*f)(int*, float*, float*, float*), int* neq, float* y,
            float* t, float* tout, int* itol, float* rtol, float* atol,
            int* itask, int* istate, int* iopt, float* rwork, int* lrw,
            int* iwork, int* liw, 
            void (*jac)(int*, float*, float*, int*, int*, float*, int*),
            int* mf);

} // extern "C" 
#endif

#ifdef USE_ODEPACK

static Real tnew;
static void* fn_data;
static Quat& (*fn_Quat)(const Quat&, Real, Real, void*, Quat&);
static Vector3& (*fn_Vector3)(const Vector3&, Real, Real, void*, Vector3&);
static VectorN& (*fn_VectorN)(const VectorN&, Real, Real, void*, VectorN&);

#ifndef BUILD_ARBITRARY_PRECISION
static void fsingledouble_Vector3(int* neq, Real* t, Real* y, Real* ydot)
{
  Vector3 x(y), xdot;
  (*fn_Vector3)(x, *t, tnew - *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = xdot[i];
}

static void fsingledouble_VectorN(int* neq, Real* t, Real* y, Real* ydot)
{
  VectorN x((unsigned) *neq, y), xdot;
  (*fn_VectorN)(x, *t, tnew - *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = xdot[i];
}

static void fsingledouble_Quat(int* neq, Real* t, Real* y, Real* ydot)
{
  Quat x, xdot;
  for (int i=0; i< *neq; i++)
    x[i] = y[i];
  (*fn_Quat)(x, *t, tnew - *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = xdot[i];
}
#else
static void farb_Vector3(int* neq, double* t, double* y, double* ydot)
{
  Vector3 x, xdot;
  for (int i=0; i< *neq; i++)
    x[i] = (Real) y[i];
  (*fn_Vector3)(x, (Real) *t, tnew - (Real) *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = (double) xdot[i];
}

static void farb_VectorN(int* neq, double* t, double* y, double* ydot)
{
  VectorN x(*neq), xdot;
  for (int i=0; i< *neq; i++)
    x[i] = (Real) y[i];
  (*fn_VectorN)(x, (Real) *t, tnew - (Real) *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = (double) xdot[i];
}

static void farb_Quat(int* neq, double* t, double* y, double* ydot)
{
  Quat x, xdot;
  for (int i=0; i< *neq; i++)
    x[i] = (Real) y[i];
  (*fn_Quat)(x, (Real) *t, tnew - (Real) *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = (double) xdot[i];
}
#endif // #else
#endif // #ifdef USE_ODEPACK

template <class T>
ODEPACKIntegrator<T>::ODEPACKIntegrator()
{
  use_stiff_integrator = false;
}

/// Implements Base::load_from_xml()
template <class T>
void ODEPACKIntegrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "ODEPACKIntegrator") == 0);
  VariableStepIntegrator<T>::load_from_xml(node, id_map); 

  // get whether the stiff integrator should be used
  const XMLAttrib* stiff_attrib = node->get_attrib("stiff");
  if (stiff_attrib)
    use_stiff_integrator = stiff_attrib->get_bool_value();
}

/// Implements Base::save_to_xml()
template <class T>
void ODEPACKIntegrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  VariableStepIntegrator<T>::save_to_xml(node, shared_objects); 
  node->name = "ODEPACKIntegrator";

  // save the method data
  node->attribs.insert(XMLAttrib("stiff", use_stiff_integrator));
}

