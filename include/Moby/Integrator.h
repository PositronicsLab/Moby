/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * Ravelin::VectorNdhis library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include <Moby/Base.h>
#include <Moby/XMLTree.h>

namespace Moby {

/// An abstract class for an ODE integration mechanism
class Integrator : public virtual Base
{
  public:
    Integrator() { }
    virtual ~Integrator() {}
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// Determines whether this is a variable-stepping integrator (false, by default, for inherited classes)
    virtual bool is_variable() const { return false; }

    /// Methods that an ODE integrator must implement
    /**
     * \param x the current state variable; on return, the new state variable
     * \param f a pointer to the function of state and time to be integrated
     * \param time the current time, contains the new time on return
     * \param step_size the step size for integration
     */
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data) = 0;
}; // end class

} // end namespace

#endif

