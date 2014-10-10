/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_FEHLBERG_INTEGRATOR_H
#define _RUNGE_KUTTA_FEHLBERG_INTEGRATOR_H

#include <vector>
#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing 4th-order variable-step Runge-Kutta integration
class RungeKuttaFehlbergIntegrator : public VariableStepIntegrator
{
  public:
    RungeKuttaFehlbergIntegrator();
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

  private:
    void rkf45 ( void f ( double t, double* y, double* yp, void* data), int neqn, double* y, double* yp, double& t, double tout, bool init);
    void fehl (void f ( double t, double* y, double* yp, void* data), int neqn, double* y, double t, double h, double* yp, double* f1, double* f2, double* f3, double* f4, double* f5, double* s);
    static void fprime(double t, double* y, double* yprime, void* data);
    static double sign(double x);
    Ravelin::VectorNd _x, _xprime;
    Ravelin::VectorNd& (*_f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&);
    void* _data;
    unsigned _N;
    double _dt;
    std::vector<double> f1, f2, f3, f4, f5;
    double _h;
    int _kop;
}; // end class def

} // end namespace

#endif

