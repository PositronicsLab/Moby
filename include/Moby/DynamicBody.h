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
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Optimization.h>
#include <Moby/EventProblemData.h>
#include <Moby/RecurrentForce.h>
#include <Moby/LinAlg.h>

namespace Moby {

/// Superclass for deformable bodies and single and multi-rigid bodies  
class DynamicBody : public Visualizable
{
  friend class EventDrivenSimulator;

  public:
    enum GeneralizedCoordinateType { eRodrigues, eAxisAngle };

    DynamicBody() 
    { 
      controller = NULL; 
      _enabled = true;
    }

    virtual ~DynamicBody() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void integrate(Real t, Real h, boost::shared_ptr<Integrator<VectorN> > integrator);

    /// Forces a recalculation of forward dynamics
    virtual void calc_fwd_dyn(Real dt) = 0;

    /// Updates the event problem data matrices and vectors
    virtual void update_event_data(EventProblemData& epd) = 0;

    /// Updates the body velocity using event problem data
    virtual void update_velocity(const EventProblemData& epd) = 0;

    /// Resets the force and torque accumulators on the dynamic body
    virtual void reset_accumulators() = 0;

    /// Transforms the dynamic body by the given transform
    virtual void transform(const Matrix4& T) = 0;

    /// Calculates the kinetic energy of the body
    virtual Real calc_kinetic_energy() const = 0;

    /// Gets the number of generalized coordinates
    virtual unsigned num_generalized_coordinates(GeneralizedCoordinateType gctype) const = 0;

    /// Adds a generalized force to the body
    virtual void add_generalized_force(GeneralizedCoordinateType gctype, const VectorN& gf) = 0;

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(GeneralizedCoordinateType gctype, const VectorN& gj) = 0;

    /// Gets the generalized coordinates of this body
    virtual VectorN& get_generalized_coordinates(GeneralizedCoordinateType gctype, VectorN& gc) = 0;

    /// Gets the generalized velocity of this body
    virtual VectorN& get_generalized_velocity(GeneralizedCoordinateType gctype, VectorN& gv) = 0;

    /// Gets the generalized acceleration of this body
    virtual VectorN& get_generalized_acceleration(GeneralizedCoordinateType gctype, VectorN& ga) = 0;

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const VectorN& gc) = 0;

    /// Sets the generalized velocity of this body
    /**
      * \param gv the generalized velocity
      * \note uses the current generalized coordinates
      */
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const VectorN& gv) = 0;

    /// Gets the generalized inertia of this body
    virtual MatrixN& get_generalized_inertia(GeneralizedCoordinateType gctype, MatrixN& M) = 0;

    /// Solves using the inverse generalized inertia
    virtual MatrixN& solve_generalized_inertia(GeneralizedCoordinateType gctype, const MatrixN& B, MatrixN& X) = 0;

    /// Solves using the inverse generalized inertia
    virtual VectorN& solve_generalized_inertia(GeneralizedCoordinateType gctype, const VectorN& b, VectorN& x) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual VectorN& get_generalized_forces(GeneralizedCoordinateType gctype, VectorN& f) = 0;

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param p the point at which to apply the force and torque
     * \param f the force component
     * \param t the torque component
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual VectorN& convert_to_generalized_force(GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& p, const Vector3& f, const Vector3& t, VectorN& gf) = 0;

    /// The controller callback, if any, for this body
    void (*controller)(boost::shared_ptr<DynamicBody>, Real, void*);

    /// Argument to be passed to the controller
    void* controller_arg;

    /// Gets the set of recurrent forces applied to this body
    const std::list<RecurrentForcePtr>& get_recurrent_forces() const { return _rfs; }

    /// Gets the set of recurrent forces applied to this body
    std::list<RecurrentForcePtr>& get_recurrent_forces() { return _rfs; }

    /// Gets the angular speed of this body (or maximum angular speed of the links, if this body is articulated)
    virtual Real get_aspeed() const = 0; 

    /// Sets the body to be enabled or disabled
    virtual void set_enabled(bool flag) { _enabled = flag; }

    /// Gets whether the body is enabled
    bool is_enabled() const { return _enabled; }

  private:

    /// Whether the body is enabled or disabled
    bool _enabled;

    /// Set of recurrent forces applied to this body
    std::list<RecurrentForcePtr> _rfs;

    /// For use with integration
    VectorN gc, gv, gcgv;

    static VectorN& ode_both(const VectorN& x, Real t, Real dt, void* data, VectorN& dx);
}; // end class

} // end namespace

#endif
