/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <Moby/BulirschStoerIntegrator.h>

using std::vector;
using boost::shared_ptr;
using namespace boost::numeric::odeint;
using Ravelin::VectorNd;
using namespace Moby;

// static definitions
VectorNd& (*BulirschStoerIntegrator::_f)(const VectorNd&, double, double, void*, VectorNd&);
double BulirschStoerIntegrator::_dt;
vector<double> BulirschStoerIntegrator::_y;
VectorNd BulirschStoerIntegrator::_x;
VectorNd BulirschStoerIntegrator::_dxdt;
void* BulirschStoerIntegrator::_data;

BulirschStoerIntegrator::BulirschStoerIntegrator()
{
}

void BulirschStoerIntegrator::f(const vector<double>& y, vector<double>& dydt, const double t)
{
  // copy
  std::copy(y.begin(), y.end(), _x.begin());

  // get the derivative
  _f(_x, t, _dt, _data, _dxdt);

  // copy back to dydt
  dydt.resize(_dxdt.size());
  std::copy(_dxdt.begin(), _dxdt.end(), dydt.begin());
}

void BulirschStoerIntegrator::integrate(VectorNd& x, VectorNd& (*ode)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  // setup the integrator
  bulirsch_stoer<std::vector<double> > integrator(aerr_tolerance, rerr_tolerance, 1.0, 1.0);
  
  // setup y
  _y.resize(x.size());
  std::copy(x.begin(), x.end(), _y.begin());

  // setup data
  _x.resize(x.size());
  _dxdt.resize(x.size());
  _dt = step_size;
  _f = ode;
  _data = data;

  // integrate adaptively
  integrate_adaptive(integrator, f, _y, time, time+step_size, step_size); 

  // copy _y back to x
  std::copy(_y.begin(), _y.end(), x.begin());
}

void BulirschStoerIntegrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "BulirschStoerIntegrator") == 0);

  // call the parent method
  VariableStepIntegrator::load_from_xml(node, id_map);
}

void BulirschStoerIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent method 
  VariableStepIntegrator::save_to_xml(node, shared_objects); 

  // rename the node
  node->name = "BulirschStoerIntegrator";
}


