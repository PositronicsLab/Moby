/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _EULER_INTEGRATOR_H
#define _EULER_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 1st-order Euler integration
template <class T>
class EulerIntegrator : public Integrator<T>
{
  public:
    /// Euler integrator is explicit by default
    EulerIntegrator() { semi_implicit = false; }
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
}; // end class

// include inline functions
#include "EulerIntegrator.inl"

} // end namespace

#endif

