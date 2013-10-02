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
#include <Moby/Event.h>
#include <Moby/RecurrentForce.h>

namespace Moby {

/// Superclass for deformable bodies and single and multi-rigid bodies  
class DynamicBody : public Visualizable
{
  friend class EventDrivenSimulator;

  public:
    enum GeneralizedCoordinateType { eEuler, eSpatial };

    DynamicBody() 
    { 
      controller = NULL; 
    }

    virtual ~DynamicBody() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);

    /// The Jacobian transforms from the generalized coordinate from to the given frame
    virtual Ravelin::MatrixNd& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J) = 0;
    virtual Ravelin::MatrixNd& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J) = 0;

    /// Sets the computation frame type for this body
    virtual void set_computation_frame_type(ReferenceFrameType rftype) = 0;

    /// Gets the computation frame type for this body
    ReferenceFrameType get_computation_frame_type() const { return _rftype; }

    /// Forces a recalculation of forward dynamics
    virtual void calc_fwd_dyn() = 0;

    /// Resets the force and torque accumulators on the dynamic body
    virtual void reset_accumulators() = 0;

    /// Rotates the dynamic body by the given orientation 
    virtual void rotate(const Ravelin::Quatd& q) = 0;

    /// Translates the dynamic body by the given translation 
    virtual void translate(const Ravelin::Origin3d& o) = 0;

    /// Calculates the kinetic energy of the body
    virtual double calc_kinetic_energy() = 0;

    /// Gets the frame for generalized coordinates
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_gc_pose() const = 0;

    /// Gets the number of generalized coordinates
    virtual unsigned num_generalized_coordinates(GeneralizedCoordinateType gctype) const = 0;

    /// Sets the generalized forces on the body
    virtual void set_generalized_forces(const Ravelin::VectorNd& gf) = 0;

    /// Adds a generalized force to the body
    virtual void add_generalized_force(const Ravelin::VectorNd& gf) = 0;

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(const Ravelin::VectorNd& gj) = 0;

    /// Gets the generalized coordinates of this body
    virtual Ravelin::VectorNd& get_generalized_coordinates(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc) = 0;

    /// Gets the generalized velocity of this body
    virtual Ravelin::VectorNd& get_generalized_velocity(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv) = 0;

    /// Gets the generalized acceleration of this body
    virtual Ravelin::VectorNd& get_generalized_acceleration(Ravelin::VectorNd& ga) = 0;

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc) = 0;

    /// Sets the generalized velocity of this body
    /**
      * \param gv the generalized velocity
      * \note uses the current generalized coordinates
      */
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv) = 0;

    /// Gets the generalized inertia of this body
    virtual Ravelin::MatrixNd& get_generalized_inertia(Ravelin::MatrixNd& M) = 0;

    /// Solves using the inverse generalized inertia
    virtual Ravelin::MatrixNd& solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X) = 0;

    /// Solves using the inverse generalized inertia
    Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);

    /// Solves using the inverse generalized inertia
    Ravelin::MatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::MatrixNd& X);

    /// Solves using the inverse generalized inertia
    Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::SharedMatrixNd& X);

    /// Solves using the inverse generalized inertia
    virtual Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::VectorNd& x) = 0;

    /// Solves using the inverse generalized inertia
    Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::VectorNd& x);

    /// Solves using the inverse generalized inertia
    Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::SharedVectorNd& x);

    /// Solves using the inverse generalized inertia
    Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::SharedVectorNd& x);

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual Ravelin::MatrixNd& transpose_solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X) = 0;

    /// Solves the transpose matrix using the inverse generalized inertia
    Ravelin::MatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::MatrixNd& X);

    /// Solves the transpose matrix using the inverse generalized inertia
    Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::SharedMatrixNd& X);

    /// Solves the transpose matrix using the inverse generalized inertia
    Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X);

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& get_generalized_forces(Ravelin::VectorNd& f) = 0;

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the force 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, const Point3d& p, Ravelin::VectorNd& gf) = 0;

    /// The controller callback, if any, for this body
    void (*controller)(boost::shared_ptr<DynamicBody>, double, void*);

    /// Argument to be passed to the controller
    void* controller_arg;

    /// Gets the set of recurrent forces applied to this body
    const std::list<RecurrentForcePtr>& get_recurrent_forces() const { return _rfs; }

    /// Gets the set of recurrent forces applied to this body
    std::list<RecurrentForcePtr>& get_recurrent_forces() { return _rfs; }

    /// Gets the angular speed of this body (or maximum angular speed of the links, if this body is articulated)
    virtual double get_aspeed() = 0; 

  protected:

    /// The computation frame type
    ReferenceFrameType _rftype;

  private:

    /// Set of recurrent forces applied to this body
    std::list<RecurrentForcePtr> _rfs;

    /// Temporaries for use with integration
    Ravelin::VectorNd gc, gv, gcgv, xp, xv, xa;

    static Ravelin::VectorNd& ode_both(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
}; // end class

} // end namespace

#endif
