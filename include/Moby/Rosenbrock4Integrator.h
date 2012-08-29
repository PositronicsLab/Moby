/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _ROSENBROCK_IMPLICIT_INTEGRATOR_H
#define _ROSENBROCK_IMPLICIT_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order implicit Rosenbrock integration with adaptive step size control
template <class T>
class Rosenbrock4Integrator : public Integrator<T>
{
  public:
    Rosenbrock4Integrator() { rel_err_tol = NEAR_ZERO; }
    virtual void integrate(T& x, T (*f)(const T&, Real, Real, void*), Real& time, Real step_size, void* data);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// The relative error tolerance to be achieved
    Real rel_err_tol;

  private:
    void step(T& x, T (*f)(const T&, Real, Real, void*), Real& time, Real step_size, Real& tnext, void* data);

}; // end class def

// include inline functions
#include "Rosenbrock4Integrator.inl"

} // end namespace

#endif

