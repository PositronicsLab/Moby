/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _VARIABLE_EULER_INTEGRATOR_H
#define _VARIABLE_EULER_INTEGRATOR_H

#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing 4th-order variable-step Runge-Kutta integration
class VariableEulerIntegrator : public VariableStepIntegrator
{
  public:
    VariableEulerIntegrator();
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
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
    Ravelin::VectorNd x0, x1, dx;
    double calc_abs_err(const Ravelin::VectorNd& x0, const Ravelin::VectorNd& x1);
    double calc_rel_err(const Ravelin::VectorNd& x0, const Ravelin::VectorNd& x1);
    void integrate_variable(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
}; // end class def

} // end namespace

#endif

