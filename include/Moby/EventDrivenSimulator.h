/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _EVENT_SIMULATOR_H
#define _EVENT_SIMULATOR_H

#include <map>
#include <Moby/sorted_pair>
#include <Moby/Simulator.h>
#include <Moby/ImpactEventHandler.h>
#include <Moby/Event.h>

namespace Moby {

class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// An event-driven simulator
class EventDrivenSimulator : public Simulator
{
  friend class CollisionDetection;
  friend class DeformableCCD;

  public:
    EventDrivenSimulator();
    virtual ~EventDrivenSimulator() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void output_object_state(std::ostream& out) const;
    virtual Real step(Real dt);

    /// Gets the shared pointer for this
    boost::shared_ptr<EventDrivenSimulator> get_this() { return boost::dynamic_pointer_cast<EventDrivenSimulator>(shared_from_this()); }
    
    /// The amount of constraint violation allowable per simulated second 
    /**
     * Reducing this scalar will will generally increase the time
     * taken to solve the event formulations but will decrease
     * constraint violation- to a point. If this scalar is too small,
     * the simulation may be unable to progress (i.e., a virtual Zeno point
     * is encountered). If this scalar is too large, (obviously) undesirable
     * constraint violation may result. The default value is machine eps^(1/4).
     */
    Real constraint_violation_tolerance;

    /// The tolerance over the first time of impact with which additional points of contact are treated as impacting
    /**
     * If this value is set too low, then contact points might be missed; this
     * is particularly a problem for objects in resting contact.  If this value
     * is set too high, then parts of objects not immediately contacting 
     * (though impacting soon) will be considered.  Default value is NEAR_ZERO. 
     */
    Real toi_tolerance;

    /// The maximum step size taken when handling a Zeno point (default is INF)
    Real max_Zeno_step;

    /// The tolerance of whether a velocity change indicates a Zeno point
    /**
     * If this value is set too low, then the simulator may get "stuck" on an
     * unidentified Zeno point. If this value is set too high, the simulator
     * will not be quite as accurate.  Default value is NEAR_ZERO.
     */
    Real Zeno_vel_tolerance;

    /// The collision detection mechanisms
    std::list<boost::shared_ptr<CollisionDetection> > collision_detectors;

    /// Callback function after a mini-step is completed
    void (*post_mini_step_callback_fn)(EventDrivenSimulator* s);

    /// The callback function (called when events have been determined)
    /**
     * The callback function can remove events from the list, which will disable
     * their processing (however, doing so may prevent the simulation from
     * making progress, as the simulator attempts to disallow violations.
     */
    void (*event_callback_fn)(std::vector<Event>&, boost::shared_ptr<void>);

    /// The callback function (called after event impulses are applied)
    void (*event_post_impulse_callback_fn)(const std::vector<Event>&, boost::shared_ptr<void>);

    /// Data passed to event callback
    boost::shared_ptr<void> event_callback_data;
    
    /// Data passed to event impulse callback
    boost::shared_ptr<void> event_post_impulse_callback_data;
 
    /// Gets the (sorted) event data
    std::vector<Event>& get_events() { return _events; }

    /// Mapping from objects to contact parameters
    std::map<sorted_pair<BasePtr>, boost::shared_ptr<ContactParameters> > contact_params;

    bool render_contact_points;

  private:
    void handle_Zeno_point(Real dt, const std::vector<std::pair<VectorN, VectorN> >& q0, std::vector<std::pair<VectorN, VectorN> >& q1);
    static void copy(const std::vector<std::pair<VectorN, VectorN> >& source, std::vector<std::pair<VectorN, VectorN> >& dest);
    static void determine_treated_bodies(std::list<std::list<Event*> >& groups, std::vector<DynamicBodyPtr>& bodies);
    Real find_and_handle_events(Real dt, const std::vector<std::pair<VectorN, VectorN> >& q0, const std::vector<std::pair<VectorN, VectorN> >& q1, bool& Zeno);
    bool will_impact(Event& e, const std::vector<std::pair<VectorN, VectorN> >& q0, const std::vector<std::pair<VectorN, VectorN> >& q1, Real dt) const;
    void get_coords_and_velocities(std::vector<std::pair<VectorN, VectorN> >& q) const;
    void set_coords_and_velocities(const std::vector<std::pair<VectorN, VectorN> >& q0, const std::vector<std::pair<VectorN, VectorN> >& q1, Real t) const;
    void set_coords_and_velocities(const std::vector<std::pair<VectorN, VectorN> >& q) const;
    void preprocess_event(Event& e);
    void check_violation();
    void find_events(Real dt);
    void find_limit_events(const std::vector<std::pair<VectorN, VectorN> >& q0, const std::vector<std::pair<VectorN, VectorN> >& q1, Real dt, std::vector<Event>& limit_events);
    Real find_TOI(Real dt, const std::vector<std::pair<VectorN, VectorN> >& q0, const std::vector<std::pair<VectorN, VectorN> >& q1); 
    void handle_events();
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;

    // Visualization functions
    void visualize_contact( Event& event );

    /// Determines whether the simulation constraints have been violated
    bool _simulation_violated;

    /// The vector of events
    std::vector<Event> _events;

    /// Object for handling impact events
    ImpactEventHandler _impact_event_handler;
}; // end class

} // end namespace

#endif


