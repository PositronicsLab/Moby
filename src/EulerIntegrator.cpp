/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/EulerIntegrator.h>

using namespace Moby;

/// Method for 1st-order Euler integration
void EulerIntegrator::integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double& time, double step_size, void* data)
{
  // save the old time
  const double old_time = time;

  // update the time
  time += step_size;

  // compute new state
  x += (f(x, old_time, step_size, data, _dx) *= step_size);
}

/// Implements Base::load_from_xml()
void EulerIntegrator::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "EulerIntegrator") == 0);
  Integrator::load_from_xml(node, id_map); 

  // determine whether the integrator is semi-implicit, if given
  const XMLAttrib* symp_attr = node->get_attrib("semi-implicit");
  if (symp_attr)
    semi_implicit = symp_attr->get_bool_value();
}

/// Implements Base::save_to_xml()
void EulerIntegrator::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects); 
  node->name = "EulerIntegrator";

  // save the semi-implicit variable
  node->attribs.insert(XMLAttrib("semi-implicit", semi_implicit));
}

