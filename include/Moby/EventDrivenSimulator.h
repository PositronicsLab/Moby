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
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/PenaltyConstraintHandler.h>
#include <Moby/SustainedUnilateralConstraintHandler.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/CCD.h>
#include <Moby/UnilateralConstraint.h>

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

    /// The callback function (called when constraints have been determined)
    /**
     * The callback function can remove constraints from the list, which will disable
     * their processing (however, doing so may prevent the simulation from
     * making progress, as the simulator attempts to disallow violations.
     */
    void (*constraint_callback_fn)(std::vector<UnilateralConstraint>&, boost::shared_ptr<void>);

    /// The callback function (called after forces/impulses are applied)
    void (*constraint_post_callback_fn)(const std::vector<UnilateralConstraint>&, boost::shared_ptr<void>);

    /// Data passed to unilateral constraint callback
    boost::shared_ptr<void> constraint_callback_data;
    
    /// Data passed to post-constraint callback
    boost::shared_ptr<void> constraint_post_callback_data;
 
    /// Gets the (sorted) compliant constraint data
    std::vector<UnilateralConstraint>& get_compliant_constraints() { return _compliant_constraints; }

    /// Gets the (sorted) rigid constraint data
    std::vector<UnilateralConstraint>& get_rigid_constraints() { return _rigid_constraints; }

    /// Mapping from objects to contact parameters
    std::map<sorted_pair<BasePtr>, boost::shared_ptr<ContactParameters> > contact_params;

    /// If set to 'true' event driven simulator will process contact points for rendering
    bool render_contact_points;

    /// User time spent by collision detection on the last step
    double coldet_time;

    /// User time spent by constraint handling on the last step
    double constraint_time;

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

    /// The maximum time allowed for processing constraints
    double max_constraint_time;

  protected:
    void calc_impacting_unilateral_constraint_forces(double dt);
    virtual double check_pairwise_constraint_violations(double t);
    void validate_limit_estimates();
    void find_unilateral_constraints(double min_contact_dist);
    void integrate_velocities_Euler(double dt);
    void integrate_positions_Euler(double dt);
    void calc_fwd_dyn();
    void calc_compliant_unilateral_constraint_forces();
    void preprocess_constraint(UnilateralConstraint& e);
    void determine_geometries();
    void broad_phase(double dt);
    void calc_pairwise_distances();
    void visualize_contact( UnilateralConstraint& constraint );

    /// Object for handling impact constraints
    ImpactConstraintHandler _impact_constraint_handler;

    /// Object for handling penalty constraints
    PenaltyConstraintHandler _penalty_constraint_handler;
    
    /// The vector of rigid constraints
    std::vector<UnilateralConstraint> _rigid_constraints;

    /// The vector of compliant constraints
    std::vector<UnilateralConstraint> _compliant_constraints;

  private:
    enum IntegrationResult { eIntegrationSuccessful, eVelocityLimitExceeded, eMinStepReached };
    IntegrationResult integrate_generic(double dt, clock_t& start);

    struct UnilateralConstraintCmp
    {
      bool operator()(const UnilateralConstraint& e1, const UnilateralConstraint& e2) const;
    };

    template <class ForwardIterator>
    double integrate_with_sustained_constraints(double step_size, ForwardIterator begin, ForwardIterator end);

    /// Integrates all dynamic bodies
    double integrate_with_sustained_constraints(double step_size) { return integrate_with_sustained_constraints(step_size, _bodies.begin(), _bodies.end()); }

    void calc_rigid_sustained_unilateral_constraint_forces();
    void check_constraint_velocity_violations(double t);
    static Ravelin::VectorNd& ode_sustained_constraints(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
    void save_state();
    void restore_state();
    void step_si_Euler(double dt);
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
    std::map<UnilateralConstraint, double, UnilateralConstraintCmp> _zero_velocity_tolerances;

    /// Object for handling sustained rigid unilateral constraints 
    SustainedUnilateralConstraintHandler _rigid_unilateral_constraint_handler;

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

    /// Geometric pairs that should be checked for unilateral constraints (according to broad phase collision detection)
    std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> > _pairs_to_check;
}; // end class

#include "EventDrivenSimulator.inl"

} // end namespace

#endif


