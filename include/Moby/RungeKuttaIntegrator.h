/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_INTEGRATOR_H
#define _RUNGE_KUTTA_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order Runge-Kutta integration
class RungeKuttaIntegrator : public Integrator
{
  public:
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
}; // end class

}  // end namespace

#endif

