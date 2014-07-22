/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/RungeKuttaIntegrator.h>

using boost::shared_ptr;
using Ravelin::VectorNd;
using namespace Moby;

/// Method for 4th-order Runge-Kutta integration
void RungeKuttaIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  const double ONE_SIXTH = 1.0/6.0;
  const double ONE_THIRD = 1.0/3.0;
  SAFESTATIC VectorNd k1, k2, k3, k4, tmp;
  
  // compute k1
  f(x, time, step_size, data, k1) *= step_size;

  // compute k2
  const double HALF_STEP = step_size * (double) 0.5;
  ((tmp = k1) *= 0.5) += x;
  f(tmp, time + HALF_STEP, HALF_STEP, data, k2) *= step_size;

  // compute k3
  ((tmp = k2) *= 0.5) += x;
  f(tmp, time + HALF_STEP, HALF_STEP, data, k3) *= step_size;

  // compute k4
  ((tmp = k3) *= 0.5) += x;
  f(tmp, time + step_size, step_size, data, k4) *= step_size;

  // update the time
  time += step_size;

  // compute new state
  k1 *= ONE_SIXTH;
  k2 *= ONE_THIRD;
  k3 *= ONE_THIRD;
  k4 *= ONE_SIXTH;
  x += k1;
  x += k2;
  x += k3;
  x += k4;

  return;
}

void RungeKuttaIntegrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "RungeKuttaIntegrator") == 0);
  Integrator::load_from_xml(node, id_map); 
}

void RungeKuttaIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects);
  node->name = "RungeKuttaIntegrator";
}


