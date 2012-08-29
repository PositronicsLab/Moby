/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Method for 4th-order implicit Runge-Kutta integration
template <class T>
void RungeKuttaImplicitIntegrator<T>::integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data)
{
  // create dx/dt1 and dx/dt2
  T dxdt1;
  f(x, time, step_size, data, dxdt1);
  T dxdt2 = (const T&) dxdt1;

  // now traverse two steps with half step length
  const Real HALF_STEP = step_size * (Real) 0.5;
  step(x, f, time, HALF_STEP, data, dxdt1, dxdt2);
  f(x, time + step_size/2.0, HALF_STEP, data, dxdt1);
  dxdt2 = (const T&) dxdt1;
  step(x, f, time+HALF_STEP, HALF_STEP, data, dxdt1, dxdt2);
  
  // update the time
  time += step_size;

  return;
}

/// Takes a step
template <class T>
void RungeKuttaImplicitIntegrator<T>::step(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real time, Real step_size, void* data, T& dxdt1, T& dxdt2)
{
  const Real SQRT3 = std::sqrt(3.0);
  const Real IR3 = 1.0 / SQRT3;
  const unsigned ITER_STEPS = 3;

  SAFESTATIC T vtmp, xtmp1, xtmp2;
  
  // iterative solution of Y1 and Y2
  for (unsigned nu=0; nu < ITER_STEPS; nu++)
  {
    vtmp.copy_from(dxdt2) *= (0.5 * (0.5 - IR3));
    xtmp1.copy_from(dxdt1) *= (Real) 0.25;
    xtmp1 += vtmp;
    xtmp1 *= step_size; 
    xtmp1 += x;
    vtmp.copy_from(dxdt1) *= (0.5 * (0.5 + IR3));
    xtmp2.copy_from(dxdt2) *= (Real) 0.25;
    xtmp2 *= step_size;
    xtmp2 += x;

    // evaluate the ode
    Real tx1 = time + 0.5 * step_size * (1.0 - IR3);
    Real tx2 = time + 0.5 * step_size * (1.0 + IR3);
    f(xtmp1, tx1, step_size - tx1 - time, data, dxdt1);
    f(xtmp2, tx2, step_size - tx2 - time, data, dxdt2);
  }
    
  // update x copy
  vtmp.copy_from(dxdt1);
  vtmp += dxdt2;
  vtmp *= ((Real) 0.5 * step_size);
  x += vtmp;
}

template <class T>
void RungeKuttaImplicitIntegrator<T>::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "RungeKuttaImplicitIntegrator") == 0);
  Integrator<T>::load_from_xml(node, id_map); 
}

template <class T>
void RungeKuttaImplicitIntegrator<T>::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator<T>::save_to_xml(node, shared_objects); 
  node->name = "RungeKuttaImplicitIntegrator";
}


