/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VARIABLE_STEP_INTEGRATOR_H
#define _VARIABLE_STEP_INTEGRATOR_H

#include <Moby/Integrator.h>
#include <Moby/Constants.h>

namespace Moby {

/// An abstract class for an ODE variable-step integration mechanism
template <class T>
class VariableStepIntegrator : public Integrator<T>
{
  public:
    VariableStepIntegrator()
    {
      min_step_size = (Real) 0.0;
      rerr_tolerance = aerr_tolerance = NEAR_ZERO;
    }

    virtual ~VariableStepIntegrator() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// Determines whether this is a variable-stepping integrator
    virtual bool is_variable() const { return true; }

    /// The minimum step to take (default is 0.0)
    Real min_step_size;

    /// The relative error tolerance (default is NEAR_ZERO)
    Real rerr_tolerance;

    /// The absolute error tolerance (default is NEAR_ZERO)
    Real aerr_tolerance;
}; // end class

// include inline functions
#include "VariableStepIntegrator.inl"

} // end namespace

#endif

