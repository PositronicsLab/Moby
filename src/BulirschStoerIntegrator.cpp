/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <strings.h>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <Moby/BulirschStoerIntegrator.h>

using namespace boost::numeric::odeint;
using boost::shared_ptr;
using Ravelin::VectorNd;
using namespace Moby;

// static definitions
VectorNd& (*BulirschStoerIntegrator::_f)(const VectorNd&, double, double, void*, VectorNd&);
double BulirschStoerIntegrator::_dt;
void* BulirschStoerIntegrator::_data;

BulirschStoerIntegrator::BulirschStoerIntegrator()
{
}

void BulirschStoerIntegrator::f(const VectorNd& y, VectorNd& dydt, const double t)
{
  // get the derivative
  _f(y, t, _dt, _data, dydt);
}

void BulirschStoerIntegrator::integrate(VectorNd& x, VectorNd& (*f)(const VectorNd&, double, double, void*, VectorNd&), double time, double step_size, void* data)
{
  // integrate adaptively
  bulirsch_stoer<VectorNd> integrator(aerr_tolerance, rerr_tolerance, 1.0, 1.0);

  // setup data
  _dt = step_size;
  _f = f;
  _data = data;
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


