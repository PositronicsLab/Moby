/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _EVENT_SIMULATOR_H
#define _EVENT_SIMULATOR_H

#include <map>
#include <Moby/sorted_pair>
#include <Moby/Simulator.h>
#include <Moby/ImpactEventHandler.h>
#include <Moby/AccelerationEventHandler.h>
#include <Moby/PairwiseDistInfo.h>
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

    /// stepping timings
    double step_times[6];

    /// stepping statistics
    unsigned step_stats[6];

    /// the minimum integration step over a single step(.) call
    double int_min_step_stat;

    /// the maximum integration step over a single step(.) call
    double int_max_step_stat;

    /// the mean integration step over a single step(.) call
    double int_mean_step_stat;

    /// The relative error tolerance for adaptive stepping (default=1e-8)
    double rel_err_tol;

    /// The absolute error tolerance for adaptive stepping (default=1e-8)
    double abs_err_tol;

    /// The minimum step size (default=1e-5)
    double minimum_step;

    /// The maximum time allowed for processing events
    double max_event_time;

  protected:
    virtual double check_pairwise_constraint_violations(double t);
    void validate_limit_estimates();
    void find_events(double min_contact_dist);
    void integrate_velocities_Euler(double dt);
    void integrate_positions_Euler(double dt);
    void calc_fwd_dyn() const;
    void preprocess_event(Event& e);
    void determine_geometries();
    void broad_phase(double dt);
    void calc_pairwise_distances();
    void visualize_contact( Event& event );

    /// Object for handling impact events
    ImpactEventHandler _impact_event_handler;

    /// The vector of events
    std::vector<Event> _events;

  private:
    enum IntegrationResult { eIntegrationSuccessful, eVelocityLimitExceeded, eMinStepReached };
    IntegrationResult integrate_generic(double dt, clock_t& start);

    struct EventCmp
    {
      bool operator()(const Event& e1, const Event& e2) const;
    };

    template <class ForwardIterator>
    double integrate_with_accel_events(double step_size, ForwardIterator begin, ForwardIterator end);

    /// Integrates all dynamic bodies
    double integrate_with_accel_events(double step_size) { return integrate_with_accel_events(step_size, _bodies.begin(), _bodies.end()); }

    void handle_acceleration_events();
    void check_constraint_velocity_violations(double t);
    static Ravelin::VectorNd& ode_accel_events(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
    void save_state();
    void restore_state();
    void step_si_Euler(double dt);
    static void determine_treated_bodies(std::list<std::list<Event*> >& groups, std::vector<DynamicBodyPtr>& bodies);
    void handle_events();
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;
    double calc_CA_step();
    double calc_next_CA_step(double contact_dist_thresh) const;
    double calc_next_CA_Euler_step(double contact_dist_thresh) const;
    void update_constraint_violations(const std::vector<PairwiseDistInfo>& pairwise_distances);
    void reset_limit_estimates() const;

    /// The continuous collision detection mechanism
    mutable CCD _ccd;

    /// Pairwise distances at bodies' current configurations
    std::vector<PairwiseDistInfo> _pairwise_distances;

    /// The derivative at the current time
    Ravelin::VectorNd _current_accel_dx;

    /// Work vector
    Ravelin::VectorNd _workV;

    /// Interpenetration constraint violation tolerances
    std::map<sorted_pair<CollisionGeometryPtr>, double> _ip_tolerances;

    /// Velocity tolerances
    std::map<Event, double, EventCmp> _zero_velocity_tolerances;

    /// Object for handling acceleration events
    AccelerationEventHandler _accel_event_handler;

    /// The Euler step size
    double euler_step;

    /// The minimum step for advancement
    double min_advance;

    /// The distance threshold for a contact to be handled as an impact
    double impacting_contact_dist_thresh; 

    /// The distance threshold for a contact to be handled as a sustained contact 
    double sustained_contact_dist_thresh;

    /// The geometries in the simulator
    std::vector<CollisionGeometryPtr> _geometries;

    /// Geometric pairs that should be checked for events (according to broad phase collision detection)
    std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> > _pairs_to_check;
}; // end class

#include "EventDrivenSimulator.inl"

} // end namespace

#endif


