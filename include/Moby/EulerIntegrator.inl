/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Method for 1st-order Euler integration
template <class T>
void EulerIntegrator<T>::integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data)
{
  // save the old time
  const Real old_time = time;

  // update the time
  time += step_size;

  // compute new state
  T dx;
  x += (f(x, old_time, step_size, data, dx) * step_size);
}

/// Implements Base::load_from_xml()
template <class T>
void EulerIntegrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "EulerIntegrator") == 0);
  Integrator<T>::load_from_xml(node, id_map); 

  // determine whether the integrator is semi-implicit, if given
  const XMLAttrib* symp_attr = node->get_attrib("semi-implicit");
  if (symp_attr)
    semi_implicit = symp_attr->get_bool_value();
}

/// Implements Base::save_to_xml()
template <class T>
void EulerIntegrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator<T>::save_to_xml(node, shared_objects); 
  node->name = "EulerIntegrator";

  // save the semi-implicit variable
  node->attribs.insert(XMLAttrib("semi-implicit", semi_implicit));
}

