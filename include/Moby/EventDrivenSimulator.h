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
#include <Moby/CCD.h>
#include <Moby/Event.h>

namespace Moby {

class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// An event-driven simulator
class EventDrivenSimulator : public Simulator
{
  friend class CollisionDetection;

  public:
    EventDrivenSimulator();
    virtual ~EventDrivenSimulator() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual double step(double dt);

    /// Determines whether two geometries are not checked
    std::set<sorted_pair<CollisionGeometryPtr> > unchecked_pairs;

    /// The coordinates vector before and after the step
    std::vector<Ravelin::VectorNd> _q0, _qf, _qsave;

    /// The velocities vector before and after the step
    std::vector<Ravelin::VectorNd> _qd0, _qdf, _qdsave;

    /// Vectors set and passed to collision detection
    std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> > _x0, _x1;

    /// Gets the shared pointer for this
    boost::shared_ptr<EventDrivenSimulator> get_this() { return boost::dynamic_pointer_cast<EventDrivenSimulator>(shared_from_this()); }
    
    /// Callback function for getting contact parameters
    boost::shared_ptr<ContactParameters> (*get_contact_parameters_callback_fn)(CollisionGeometryPtr g1, CollisionGeometryPtr g2);

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

    /// If set to 'true' event driven simulator will process contact points for rendering
    bool render_contact_points;

    /// User time spent by collision detection on the last step
    double coldet_time;

    /// User time spent by event handling on the last step
    double event_time;

    /// The relative error tolerance for adaptive Euler stepping (default=1e-8)
    double rel_err_tol;

    /// The absolute error tolerance for adaptive Euler stepping (default=1e-8)
    double abs_err_tol;

    /// The minimum step size (default=1e-5)
    double minimum_step;

    /// The maximum time allowed for processing events
    double max_event_time;

  protected:
    void check_pairwise_constraint_violations();

  private:
    double compute_next_event_time() const;
    void integrate_velocities_Euler(double dt);
    void integrate_positions_Euler(double dt);
    void save_state();
    void restore_state();
    void calc_fwd_dyn() const;
    void step_si_Euler(double dt);
    void integrate_DAE(double dt);
    static void determine_treated_bodies(std::list<std::list<Event*> >& groups, std::vector<DynamicBodyPtr>& bodies);
    void find_events();
    void preprocess_event(Event& e);
    void handle_events();
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;
    double calc_CA_step();
    void update_constraint_violations();
    void determine_geometries();
    void calculate_bounds() const;

    // Visualization functions
    void visualize_contact( Event& event );

    /// Determines whether the simulation constraints have been violated
    bool _simulation_violated;

    /// The continuous collision detection mechanism
    mutable CCD _ccd;

    /// Work vector
    Ravelin::VectorNd _workV;

    /// The vector of events
    std::vector<Event> _events;

    /// Interpenetration constraint violation tolerances
    std::map<sorted_pair<CollisionGeometryPtr>, double> _ip_tolerances;

    /// Object for handling impact events
    ImpactEventHandler _impact_event_handler;

    /// The Euler step size
    double euler_step;

    /// The geometries in the simulator
    std::list<CollisionGeometryPtr> _geometries;

    /// Geometric pairs that should be checked for events (according to broad phase collision detection)
    std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> > _pairs_to_check;
}; // end class

} // end namespace

#endif


