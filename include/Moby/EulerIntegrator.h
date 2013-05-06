/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_EULER_INTEGRATOR_H
#define _MOBY_EULER_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 1st-order Euler integration
class EulerIntegrator : public Integrator
{
  public:
    /// Euler integrator is explicit by default
    EulerIntegrator() { semi_implicit = false; }
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double& time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;


    /// If set to <b>true</b>, semi-implicit integration is used
    /**
     * Semi-implicit integration integrates the velocities first, then uses the
     * updated velocities to integrate the position; standard integration uses
     * the velocities computed from the last time step.
     */
    bool semi_implicit;

  private:

    /// temporary variable
    Ravelin::VectorNd _dx;    
}; // end class

} // end namespace

#endif

