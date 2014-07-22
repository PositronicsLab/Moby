/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_EULER_INTEGRATOR_H
#define _MOBY_EULER_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 1st-order Euler integration
class EulerIntegrator : public Integrator
{
  public:
    EulerIntegrator() {  }
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

  private:

    /// temporary variable
    Ravelin::VectorNd _dx;    
}; // end class

} // end namespace

#endif

