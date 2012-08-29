/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Variable step Euler is explicit by default 
template <class T>
VariableEulerIntegrator<T>::VariableEulerIntegrator()
{
  semi_implicit = false;
}

/// Method for 1st-order VariableEuler integration
template <class T>
void VariableEulerIntegrator<T>::integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data)
{
  const Real desired_time = time + step_size;

  while (time < desired_time)
    integrate_variable(x, f, time, desired_time - time, data);
}

/// Does variable step integration
template <class T>
void VariableEulerIntegrator<T>::integrate_variable(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data)
{
  Real H2 = step_size * (Real) 0.5;

  // compute a step at step_size and two steps at step_size/2
  x0.copy_from(x);
  x0 += (f(x, time, step_size, data, dx) * step_size); 
  x1.copy_from(x);
  x1 += (f(x,time, H2, data, dx) * H2);
  x1 += (f(x1,time+H2, H2, data,dx) * H2);

  // compute error estimate
  Real aerr = VectorN::calc_abs_err(x0, x1);
  Real rerr = VectorN::calc_rel_err(x0, x1);

  // if the error estimate is within tolerance, quit
  // NOTE: second conditional is due to relative error possibly being NaN
  if (aerr < this->aerr_tolerance && !(rerr > this->rerr_tolerance))
  {
    // update x
    x.copy_from(x1) *= (Real) 2.0;
    x -= x0;

    // update the time
    time += step_size;

    return;
  }

  // otherwise, compute the new step size
  step_size *= 0.9 * this->aerr_tolerance/aerr;

  // if the new step size is too small, update and quit
  if (step_size < this->min_step_size)
  {
    // output error message
    std::cerr << "VariableEulerIntegrator::integrate() warning- minimum step size exceeded!" << std::endl;
    std::cerr << "  absolute error estimate: " << aerr << std::endl;
    std::cerr << "  relative error estimate: " << rerr << std::endl;
    std::cerr << "  step size: " << step_size << std::endl;
  
    // update x and time
    x.copy_from(x1) *= (Real) 2.0;
    x -= x0;
    time += step_size;
    
    // quit
    return;
  }

  // call integrator again, with lower step size
  integrate_variable(x, f, time, step_size, data); 
}

/// Implements Base::load_from_xml()
template <class T>
void VariableEulerIntegrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "VariableEulerIntegrator") == 0);

  // call the parent method
  VariableStepIntegrator<T>::load_from_xml(node, id_map); 

  // determine whether the integrator is semi-implicit, if given
  const XMLAttrib* symp_attr = node->get_attrib("semi-implicit");
  if (symp_attr)
    semi_implicit = symp_attr->get_bool_value();
}

/// Implements Base::save_to_xml()
template <class T>
void VariableEulerIntegrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // call the parent method 
  VariableStepIntegrator<T>::save_to_xml(node, shared_objects); 

  // rename the node
  node->name = "VariableEulerIntegrator";

  // save the semi-implicit variable
  node->attribs.insert(XMLAttrib("semi-implicit", semi_implicit));
}

