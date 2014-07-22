/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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

} // extern "C" 
#endif

#ifdef USE_ODEPACK

static double tnew;
static void* fn_data;
static Ravelin::VectorNd& (*fn_Ravelin::VectorNd)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&);

static void fdouble_Ravelin::VectorNd(int* neq, double* t, double* y, double* ydot)
{
  Ravelin::VectorNd x((unsigned) *neq, y), xdot;
  (*fn_Ravelin::VectorNd)(x, *t, tnew - *t, fn_data, xdot);
  for (int i=0; i< *neq; i++)
    ydot[i] = xdot[i];
}

#endif // #ifdef USE_ODEPACK

ODEPACKIntegrator::ODEPACKIntegrator()
{
  use_stiff_integrator = false;
}

/// Implements Base::load_from_xml()
void ODEPACKIntegrator::load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "ODEPACKIntegrator") == 0);
  VariableStepIntegrator::load_from_xml(node, id_map); 

  // get whether the stiff integrator should be used
  const XMLAttrib* stiff_attrib = node->get_attrib("stiff");
  if (stiff_attrib)
    use_stiff_integrator = stiff_attrib->get_bool_value();
}

/// Implements Base::save_to_xml()
void ODEPACKIntegrator::save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const
{ 
  VariableStepIntegrator::save_to_xml(node, shared_objects); 
  node->name = "ODEPACKIntegrator";

  // save the method data
  node->attribs.insert(XMLAttrib("stiff", use_stiff_integrator));
}

