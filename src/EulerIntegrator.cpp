/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/EulerIntegrator.h>

using boost::shared_ptr;
using Ravelin::VectorNd;
using namespace Moby;

/// Method for 1st-order Euler integration
void EulerIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  // save the old time
  const double old_time = time;

  // update the time
  time += step_size;

  // compute new state
  x += (f(x, old_time, step_size, data, _dx) *= step_size);
}

/// Implements Base::save_to_xml()
void EulerIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects); 
  node->name = "EulerIntegrator";
}

