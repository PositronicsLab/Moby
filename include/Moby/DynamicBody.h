/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _DYNAMIC_BODY_H
#define _DYNAMIC_BODY_H

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <Moby/Base.h>
#include <Moby/Constants.h>
#include <Moby/Integrator.h>
#include <Moby/Log.h>
#include <Moby/Visualizable.h>
#include <Moby/EventProblemData.h>
#include <Moby/RecurrentForce.h>

namespace Moby {

/// Superclass for deformable bodies and single and multi-rigid bodies  
class DynamicBody : public Visualizable
{
  public:
    enum GeneralizedCoordinateType { eRodrigues, eAxisAngle };

    DynamicBody() 
    { 
      controller = NULL; 
    }

    virtual ~DynamicBody() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);

    /// Forces a recalculation of forward dynamics
    virtual void calc_fwd_dyn(double dt) = 0;

    /// Updates the event problem data matrices and vectors
    virtual void update_event_data(EventProblemData& epd) = 0;

    /// Updates the body velocity using event problem data
    virtual void update_velocity(const EventProblemData& epd) = 0;

    /// Resets the force and torque accumulators on the dynamic body
    virtual void reset_accumulators() = 0;

    /// Transforms the dynamic body by the given transform
    virtual void transform(const Ravelin::Pose3d& T) = 0;

    /// Calculates the kinetic energy of the body
    virtual double calc_kinetic_energy() const = 0;

    /// Gets the number of generalized coordinates
    virtual unsigned num_generalized_coordinates(GeneralizedCoordinateType gctype) const = 0;

    /// Adds a generalized force to the body
    virtual void add_generalized_force(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gf) = 0;

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gj) = 0;

    /// Gets the generalized coordinates of this body
    virtual Ravelin::VectorNd& get_generalized_coordinates(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc) = 0;

    /// Gets the generalized velocity of this body
    virtual Ravelin::VectorNd& get_generalized_velocity(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv) = 0;

    /// Gets the generalized acceleration of this body
    virtual Ravelin::VectorNd& get_generalized_acceleration(GeneralizedCoordinateType gctype, Ravelin::VectorNd& ga) = 0;

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc) = 0;

    /// Sets the generalized velocity of this body
    /**
      * \param gv the generalized velocity
      * \note uses the current generalized coordinates
      */
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv) = 0;

    /// Gets the generalized inertia of this body
    virtual Ravelin::MatrixNd& get_generalized_inertia(GeneralizedCoordinateType gctype, Ravelin::MatrixNd& M) = 0;

    /// Solves using the inverse generalized inertia
    virtual Ravelin::MatrixNd& solve_generalized_inertia(GeneralizedCoordinateType gctype, const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X) = 0;

    /// Solves using the inverse generalized inertia
    virtual Ravelin::VectorNd& solve_generalized_inertia(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& b, Ravelin::VectorNd& x) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& get_generalized_forces(GeneralizedCoordinateType gctype, Ravelin::VectorNd& f) = 0;

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the wrench 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& convert_to_generalized_force(GeneralizedCoordinateType gctype, SingleBodyPtr body, const Ravelin::Wrenchd& w, Ravelin::VectorNd& gf) = 0;

    /// The controller callback, if any, for this body
    void (*controller)(boost::shared_ptr<DynamicBody>, double, void*);

    /// Argument to be passed to the controller
    void* controller_arg;

    /// Gets the set of recurrent forces applied to this body
    const std::list<RecurrentForcePtr>& get_recurrent_forces() const { return _rfs; }

    /// Gets the set of recurrent forces applied to this body
    std::list<RecurrentForcePtr>& get_recurrent_forces() { return _rfs; }

    /// Gets the angular speed of this body (or maximum angular speed of the links, if this body is articulated)
    virtual double get_aspeed() const = 0; 

  private:

    /// Set of recurrent forces applied to this body
    std::list<RecurrentForcePtr> _rfs;

    /// Temporaries for use with integration
    Ravelin::VectorNd gc, gv, gcgv, xp, xv, xa;

    static Ravelin::VectorNd& ode_both(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
}; // end class

} // end namespace

#endif
