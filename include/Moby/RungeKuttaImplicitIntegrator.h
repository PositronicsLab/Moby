/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_IMPLICIT_INTEGRATOR_H
#define _RUNGE_KUTTA_IMPLICIT_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order implicit Runge-Kutta integration
template <class T>
class RungeKuttaImplicitIntegrator : public Integrator<T>
{
  public:
    virtual void integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
    static void step(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real time, Real step_size, void* data, T& dxdt1, T& dxdt2);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
}; // end class def

// include inline functions
#include "RungeKuttaImplicitIntegrator.inl"

} // end namespace

#endif

