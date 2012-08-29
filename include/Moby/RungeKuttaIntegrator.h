/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_INTEGRATOR_H
#define _RUNGE_KUTTA_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order Runge-Kutta integration
template <class T>
class RungeKuttaIntegrator : public Integrator<T>
{
  public:
    virtual void integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
}; // end class

// include inline functions
#include "RungeKuttaIntegrator.inl"

}  // end namespace

#endif

