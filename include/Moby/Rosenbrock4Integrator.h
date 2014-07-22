/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _ROSENBROCK_IMPLICIT_INTEGRATOR_H
#define _ROSENBROCK_IMPLICIT_INTEGRATOR_H

#include <Ravelin/LinAlgd.h>
#include <Moby/Constants.h>
#include <Moby/Integrator.h>

namespace Moby {

/// A class for performing 4th-order implicit Rosenbrock integration with adaptive step size control
class Rosenbrock4Integrator : public Integrator
{
  public:
    Rosenbrock4Integrator() { rel_err_tol = NEAR_ZERO; }
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd (*f)(const Ravelin::VectorNd&, double, double, void*), double time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The relative error tolerance to be achieved
    double rel_err_tol;

  private:
    void step(Ravelin::VectorNd& x, Ravelin::VectorNd (*f)(const Ravelin::VectorNd&, double, double, void*), double& time, double step_size, double& tnext, void* data);
    Ravelin::LinAlgd _LA;
}; // end class def

} // end namespace

#endif

