/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _ODEPACK_INTEGRATOR_H
#define _ODEPACK_INTEGRATOR_H

#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing integration using the ODEPACK library (optional)
template <class T>
class ODEPACKIntegrator : public VariableStepIntegrator<T>
{
  public:
    ODEPACKIntegrator();
    virtual void integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// If 'true', implicit integration is used (for "stiff" differential eqns)
    bool use_stiff_integrator;
}; // end class

// include inline functions
#include "ODEPACKIntegrator.inl"

} // end namespace

#endif

