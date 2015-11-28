/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _EVENT_SIMULATOR_H
#define _EVENT_SIMULATOR_H

#include <map>
#include <Ravelin/sorted_pair>
#include <Moby/Simulator.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/PenaltyConstraintHandler.h>
#include <Moby/SustainedUnilateralConstraintHandler.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/CCD.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/TimeSteppingSimulator.h>

namespace Moby {

class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// An event-driven simulator
class EventDrivenSimulator : public TimeSteppingSimulator
{
  friend class CollisionDetection;
  friend class ConstraintStabilization;

  public:
    EventDrivenSimulator();
    virtual ~EventDrivenSimulator() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual double step(double dt);
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;

    /// The coordinates vector before and after the step
    std::vector<Ravelin::VectorNd> _q0, _qf, _qsave;

    /// The velocities vector before and after the step
    std::vector<Ravelin::VectorNd> _qd0, _qdf, _qdsave;

    /// Vectors set and passed to collision detection
    std::vector<std::pair<ControlledBodyPtr, Ravelin::VectorNd> > _x0, _x1;

    /// Gets the shared pointer for this
    boost::shared_ptr<EventDrivenSimulator> get_this() { return boost::dynamic_pointer_cast<EventDrivenSimulator>(shared_from_this()); }
    
    /// The Euler step size
    /**
     * Below this step size, a time-stepping method is triggered. The
     * stepping method might be necessary if inconsistent contact 
     * configurations (Painleve' Paradox type configurations) are encountered.
     * The time-stepping method is first order.
     */ 
    double euler_step;

  protected:
    virtual double check_pairwise_constraint_violations(double t);
    void validate_limit_estimates();

  private:
    enum IntegrationResult { eIntegrationSuccessful, eVelocityLimitExceeded, eMinStepReached };
    IntegrationResult integrate_generic(double& dt);

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
    double calc_CA_step();
    double calc_next_CA_step(double contact_dist_thresh) const;
    void reset_limit_estimates() const;
    void update_constraint_violations(const std::vector<PairwiseDistInfo>& pairwise_distances);

    /// Interpenetration constraint violation tolerances
    std::map<sorted_pair<CollisionGeometryPtr>, double> _ip_tolerances;

    /// Interpenetration constraint violation tolerances
    std::map<Ravelin::sorted_pair<CollisionGeometryPtr>, double> _ip_tolerances;

    /// Velocity tolerances
    /**
     * The constraint solver can only guarantee constraints are satisfied
     * to some tolerance near zero. This map allows the zero velocity tolerance
     * to be specified for any constraint. 
     */
    std::map<UnilateralConstraint, double, UnilateralConstraintCmp> _zero_velocity_tolerances;

    /// Object for handling sustained rigid unilateral constraints 
    SustainedUnilateralConstraintHandler _rigid_unilateral_constraint_handler;
}; // end class

#include "EventDrivenSimulator.inl"

} // end namespace

#endif


