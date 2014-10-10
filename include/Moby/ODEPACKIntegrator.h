/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _ODEPACK_INTEGRATOR_H
#define _ODEPACK_INTEGRATOR_H

#include <cstring>
#include <Moby/Log.h>
#include <Moby/VariableStepIntegrator.h>
#include <cstring>
#ifdef THREADSAFE
#include <pthread.h>
#endif

namespace Moby {

class ODEPACKIntegratorMutex
{
  public:
    #ifdef THREADSAFE
    static pthread_mutex_t _odepack_mutex;
    #endif
};

/// A class for performing integration using the ODEPACK library (optional)
class ODEPACKIntegrator : public VariableStepIntegrator
{
  public:
    ODEPACKIntegrator();
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// If 'true', implicit integration is used (for "stiff" differential eqns)
    bool use_stiff_integrator;

  private:
    // temporary variables
    std::vector<int> _iwork;
    std::vector<double> _rwork;
}; // end class

} // end namespace

#endif

