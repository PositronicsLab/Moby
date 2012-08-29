/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VARIABLE_EULER_INTEGRATOR_H
#define _VARIABLE_EULER_INTEGRATOR_H

#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing 4th-order variable-step Runge-Kutta integration
template <class T>
class VariableEulerIntegrator : public VariableStepIntegrator<T>
{
  public:
    VariableEulerIntegrator();
    virtual void integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// If set to <b>true</b>, semi-implicit integration is used
    /**
     * Semi-implicit integration integrates the velocities first, then uses the
     * updated velocities to integrate the position; standard integration uses
     * the velocities computed from the last time step.
     */
    bool semi_implicit;

  private:
    T x0, x1, dx;
    void integrate_variable(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
}; // end class def

// include inline functions
#include "VariableEulerIntegrator.inl"

} // end namespace

#endif

