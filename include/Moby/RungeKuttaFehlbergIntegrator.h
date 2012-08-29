/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RUNGE_KUTTA_FEHLBERG_INTEGRATOR_H
#define _RUNGE_KUTTA_FEHLBERG_INTEGRATOR_H

#include <vector>
#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing 4th-order variable-step Runge-Kutta integration
template <class T>
class RungeKuttaFehlbergIntegrator : public VariableStepIntegrator<T>
{
  public:
    RungeKuttaFehlbergIntegrator();
    virtual void integrate(T& x, T& (*f)(const T&, Real, Real, void*, T&), Real& time, Real step_size, void* data);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

  private:
    void rkf45 ( void f ( Real t, Real* y, Real* yp, void* data), int neqn, Real* y, Real* yp, Real& t, Real tout, bool init);
    void fehl (void f ( Real t, Real* y, Real* yp, void* data), int neqn, Real* y, Real t, Real h, Real* yp, Real* f1, Real* f2, Real* f3, Real* f4, Real* f5, Real* s);
    static void fprime(Real t, Real* y, Real* yprime, void* data);
    static Real sign(Real x);
    T _x, _xprime;
    T& (*_f)(const T&, Real, Real, void*, T&);
    void* _data;
    unsigned _N;
    Real _dt;
    std::vector<Real> f1, f2, f3, f4, f5;
    Real _h;
    int _kop;
}; // end class def

// include inline functions
#include "RungeKuttaFehlbergIntegrator.inl"

} // end namespace

#endif

