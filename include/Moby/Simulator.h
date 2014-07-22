/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include <sys/times.h>
#include <list>
#include <map>
#include <set>
#include <boost/shared_ptr.hpp>

#include <Moby/Base.h>
#include <Moby/Log.h>
#include <Moby/Integrator.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>

namespace osg { 
  class Node;
  class Group; 
}

namespace Moby {

class RigidBody;
class ArticulatedBody;
class VisualizationData;
class DynamicBody;

/// Simulator for both unarticulated and articulated rigid bodies without contact
/**
 * Class used for performing dynamics simulation of rigid bodies without contact.
 * Rigid body simulation of articulated bodies is supported using both
 * maximal and reduced coordinate approaches.
 */
class Simulator : public virtual Base
{
  public:
    Simulator();
    virtual ~Simulator(); 
    virtual double step(double step_size);
    DynamicBodyPtr find_dynamic_body(const std::string& name) const;
    void add_dynamic_body(DynamicBodyPtr body);
    void remove_dynamic_body(DynamicBodyPtr body);
    void update_visualization();
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);  

    /// The current simulation time
    double current_time;

    /// The integrator used to step the simulation
    boost::shared_ptr<Integrator> integrator;

    /// Gets the list of dynamic bodies in the simulator
    /**
     * \note if a dynamic body is articulated, only the articulated body is
     *       returned, not the links
     */
    const std::vector<DynamicBodyPtr>& get_dynamic_bodies() const { return _bodies; }

    void add_transient_vdata(osg::Node* vdata);

    /// Gets the persistent visualization data
    osg::Node* get_persistent_vdata() const { return (osg::Node*) _persistent_vdata; } 

    /// Gets the transient (one-step) visualization data
    osg::Node* get_transient_vdata() const { return (osg::Node*) _transient_vdata; } 

    /// Callback function after a step is completed
    void (*post_step_callback_fn)(Simulator* s);

    /// User time spent by dynamics on the last step
    double dynamics_time;

  protected:
    virtual double check_pairwise_constraint_violations(double t) { return 0.0; }
    osg::Group* _persistent_vdata;
    osg::Group* _transient_vdata;
    void update_bounds() const;

    /// The set of bodies in the simulation
    std::vector<DynamicBodyPtr> _bodies;
  
    /// The derivative at the current time
    Ravelin::VectorNd _current_dx;

    template <class ForwardIterator>
    double integrate(double step_size, ForwardIterator begin, ForwardIterator end);

    /// Integrates all dynamic bodies
    double integrate(double step_size) { return integrate(step_size, _bodies.begin(), _bodies.end()); }

  private:
    static Ravelin::VectorNd& ode(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
}; // end class

// include inline functions
#include "Simulator.inl"

} // end namespace 

#endif

