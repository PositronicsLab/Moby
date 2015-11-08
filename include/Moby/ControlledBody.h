/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CONTROLLED_BODY_H
#define _CONTROLLED_BODY_H

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <Moby/Base.h>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/Visualizable.h>
#include <Moby/RecurrentForce.h>

namespace Moby {

/// Superclass for controlled bodies 
class ControlledBody : public virtual Visualizable
{
  public:

    ControlledBody() 
    { 
      controller = NULL; 
    }

    virtual ~ControlledBody() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The controller callback, if any, for this body
    Ravelin::VectorNd& (*controller)(Ravelin::VectorNd&, double, void*);

    /// Clone of this body to be passed to the controller
    boost::shared_ptr<ControlledBody> clone;

    /// Argument to be passed to the controller
    void* controller_arg;

    /// Gets the set of recurrent forces applied to this body
    const std::list<RecurrentForcePtr>& get_recurrent_forces() const { return _rfs; }

    /// Gets the set of recurrent forces applied to this body
    std::list<RecurrentForcePtr>& get_recurrent_forces() { return _rfs; }

    /// Prepares to compute the derivative of the body (sustained constraints) 
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) = 0;

    /// Prepares to compute the derivative of the body
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) = 0;

    /// Computes the derivative of the body
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx) = 0;

  private:

    /// Set of recurrent forces applied to this body
    std::list<RecurrentForcePtr> _rfs;

  protected:

    /// Pointer to the simulator (necessary for applying impulses w/constraints)
    boost::weak_ptr<Simulator> simulator;
}; // end class

} // end namespace

#endif
