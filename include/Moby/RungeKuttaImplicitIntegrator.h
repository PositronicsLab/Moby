/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_IMPLICIT_INTEGRATOR_H
#define _RUNGE_KUTTA_IMPLICIT_INTEGRATOR_H

#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order implicit Runge-Kutta integration
class RungeKuttaImplicitIntegrator : public Integrator
{
  public:
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double& time, double step_size, void* data);
    static void step(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data, Ravelin::VectorNd& dxdt1, Ravelin::VectorNd& dxdt2);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
}; // end class def

} // end namespace

#endif

