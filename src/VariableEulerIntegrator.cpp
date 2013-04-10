/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/VariableEulerIntegrator.h>

using namespace Moby;
using namespace Ravelin;

/// Variable step Euler is explicit by default 
VariableEulerIntegrator::VariableEulerIntegrator()
{
  semi_implicit = false;
}

/// Method for 1st-order VariableEuler integration
void VariableEulerIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double& time, double step_size, void* data)
{
  const double desired_time = time + step_size;

  while (time < desired_time)
    integrate_variable(x, f, time, desired_time - time, data);
}

/// Does variable step integration
void VariableEulerIntegrator::integrate_variable(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double& time, double step_size, void* data)
{
  double H2 = step_size * (double) 0.5;

  // compute a step at step_size and two steps at step_size/2
  x0 = x;
  x0 += (f(x, time, step_size, data, dx) *= step_size); 
  x1 = x;
  x1 += (f(x,time, H2, data, dx) *= H2);
  x1 += (f(x1,time+H2, H2, data,dx) *= H2);

  // compute error estimate
  double aerr = VectorNd::calc_abs_err(x0, x1);
  double rerr = VectorNd::calc_rel_err(x0, x1);

  // if the error estimate is within tolerance, quit
  // NOTE: second conditional is due to relative error possibly being NaN
  if (aerr < this->aerr_tolerance && !(rerr > this->rerr_tolerance))
  {
    // update x
    (x = x1) *= (double) 2.0;
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
    (x = x1) *= (double) 2.0;
    x -= x0;
    time += step_size;
    
    // quit
    return;
  }

  // call integrator again, with lower step size
  integrate_variable(x, f, time, step_size, data); 
}

/// Implements Base::load_from_xml()
void VariableEulerIntegrator::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "VariableEulerIntegrator") == 0);

  // call the parent method
  VariableStepIntegrator::load_from_xml(node, id_map); 

  // determine whether the integrator is semi-implicit, if given
  const XMLAttrib* symp_attr = node->get_attrib("semi-implicit");
  if (symp_attr)
    semi_implicit = symp_attr->get_bool_value();
}

/// Implements Base::save_to_xml()
void VariableEulerIntegrator::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // call the parent method 
  VariableStepIntegrator::save_to_xml(node, shared_objects); 

  // rename the node
  node->name = "VariableEulerIntegrator";

  // save the semi-implicit variable
  node->attribs.insert(XMLAttrib("semi-implicit", semi_implicit));
}

