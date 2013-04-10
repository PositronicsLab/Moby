/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/RungeKuttaImplicitIntegrator.h>

using namespace Moby;

/// Method for 4th-order implicit Runge-Kutta integration
void RungeKuttaImplicitIntegrator::integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double& time, double step_size, void* data)
{
  // create dx/dt1 and dx/dt2
  Ravelin::VectorNd dxdt1;
  f(x, time, step_size, data, dxdt1);
  Ravelin::VectorNd dxdt2 = (const Ravelin::VectorNd&) dxdt1;

  // now traverse two steps with half step length
  const double HALF_STEP = step_size * (double) 0.5;
  step(x, f, time, HALF_STEP, data, dxdt1, dxdt2);
  f(x, time + step_size/2.0, HALF_STEP, data, dxdt1);
  dxdt2 = (const Ravelin::VectorNd&) dxdt1;
  step(x, f, time+HALF_STEP, HALF_STEP, data, dxdt1, dxdt2);
  
  // update the time
  time += step_size;

  return;
}

/// Takes a step
void RungeKuttaImplicitIntegrator::step(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data, Ravelin::VectorNd& dxdt1, Ravelin::VectorNd& dxdt2)
{
  const double SQRT3 = std::sqrt(3.0);
  const double IR3 = 1.0 / SQRT3;
  const unsigned ITER_STEPS = 3;

  SAFESTATIC Ravelin::VectorNd vtmp, xtmp1, xtmp2;
  
  // iterative solution of Y1 and Y2
  for (unsigned nu=0; nu < ITER_STEPS; nu++)
  {
    (vtmp = dxdt2) *= (0.5 * (0.5 - IR3));
    (xtmp1 = dxdt1) *= (double) 0.25;
    xtmp1 += vtmp;
    xtmp1 *= step_size; 
    xtmp1 += x;
    (vtmp = dxdt1) *= (0.5 * (0.5 + IR3));
    (xtmp2 = dxdt2) *= (double) 0.25;
    xtmp2 *= step_size;
    xtmp2 += x;

    // evaluate the ode
    double tx1 = time + 0.5 * step_size * (1.0 - IR3);
    double tx2 = time + 0.5 * step_size * (1.0 + IR3);
    f(xtmp1, tx1, step_size - tx1 - time, data, dxdt1);
    f(xtmp2, tx2, step_size - tx2 - time, data, dxdt2);
  }
    
  // update x copy
  vtmp = dxdt1;
  vtmp += dxdt2;
  vtmp *= ((double) 0.5 * step_size);
  x += vtmp;
}

void RungeKuttaImplicitIntegrator::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map) 
{ 
  assert(strcasecmp(node->name.c_str(), "RungeKuttaImplicitIntegrator") == 0);
  Integrator::load_from_xml(node, id_map); 
}

void RungeKuttaImplicitIntegrator::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects); 
  node->name = "RungeKuttaImplicitIntegrator";
}


