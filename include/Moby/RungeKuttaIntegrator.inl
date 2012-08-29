/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Method for 4th-order Runge-Kutta integration
template <class T>
void RungeKuttaIntegrator<T>::integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data)
{
  const Real ONE_SIXTH = 1.0/6.0;
  const Real ONE_THIRD = 1.0/3.0;
  SAFESTATIC T k1, k2, k3, k4, tmp;
  
  // compute k1
  f(x, time, step_size, data, k1) *= step_size;

  // compute k2
  const Real HALF_STEP = step_size * (Real) 0.5;
  (tmp.copy_from(k1) *= 0.5) += x;
  f(tmp, time + HALF_STEP, HALF_STEP, data, k2) *= step_size;

  // compute k3
  (tmp.copy_from(k2) *= 0.5) += x;
  f(tmp, time + HALF_STEP, HALF_STEP, data, k3) *= step_size;

  // compute k4
  (tmp.copy_from(k3) *= 0.5) += x;
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

template <class T>
void RungeKuttaIntegrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "RungeKuttaIntegrator") == 0);
  Integrator<T>::load_from_xml(node, id_map); 
}

template <class T>
void RungeKuttaIntegrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator<T>::save_to_xml(node, shared_objects);
  node->name = "RungeKuttaIntegrator";
}


