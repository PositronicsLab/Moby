/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _BULIRSCH_STOER_INTEGRATOR_H
#define _BULIRSCH_STOER_INTEGRATOR_H

#include <vector>
#include <Moby/VariableStepIntegrator.h>

namespace Moby {

/// A class for performing adaptive Bulirsch Stoer integration 
class BulirschStoerIntegrator : public VariableStepIntegrator
{
  public:
    BulirschStoerIntegrator();
    virtual void integrate(Ravelin::VectorNd& x, Ravelin::VectorNd& (*f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&), double time, double step_size, void* data);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

  private:
    static void f(const Ravelin::VectorNd&, Ravelin::VectorNd&, const double);
    static Ravelin::VectorNd& (*_f)(const Ravelin::VectorNd&, double, double, void*, Ravelin::VectorNd&);
    static double _dt;
    static void* _data;
}; // end class def

} // end namespace

#endif

