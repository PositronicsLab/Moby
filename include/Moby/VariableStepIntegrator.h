/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _VARIABLE_STEP_INTEGRATOR_H
#define _VARIABLE_STEP_INTEGRATOR_H

#include <Moby/Integrator.h>
#include <Moby/Constants.h>

namespace Moby {

/// An abstract class for an ODE variable-step integration mechanism
class VariableStepIntegrator : public Integrator
{
  public:
    VariableStepIntegrator()
    {
      min_step_size = (double) 0.0;
      rerr_tolerance = aerr_tolerance = NEAR_ZERO;
    }

    virtual ~VariableStepIntegrator() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// Determines whether this is a variable-stepping integrator
    virtual bool is_variable() const { return true; }

    /// The minimum step to take (default is 0.0)
    double min_step_size;

    /// The relative error tolerance (default is NEAR_ZERO)
    double rerr_tolerance;

    /// The absolute error tolerance (default is NEAR_ZERO)
    double aerr_tolerance;
}; // end class

} // end namespace

#endif

