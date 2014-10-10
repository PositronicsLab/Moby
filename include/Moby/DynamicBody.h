/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
      _kinematic_update = false;
    }

    virtual ~DynamicBody() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The Jacobian transforms from the generalized coordinate from to the given frame
    virtual Ravelin::MatrixNd& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J) = 0;
    virtual Ravelin::MatrixNd& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J) = 0;

    /// Returns true if one or more of the limit estimates has been exceeded
    virtual bool limit_estimates_exceeded() const = 0;

    /// Validates the limit estimates
    virtual void validate_limit_estimates() = 0;

    /// Validates position-based variables (potentially dangerous for a user to call)
    virtual void validate_position_variables() { };

    /// Validates velocity-based variables (potentially dangerous for a user to call)
    virtual void validate_velocity_variables() { };

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
    virtual void set_generalized_forces(const Ravelin::SharedVectorNd& gf) = 0;

    /// Sets the generalized forces on the body
    virtual void set_generalized_forces(const Ravelin::VectorNd& gf)
    {
      const Ravelin::SharedVectorNd gf_shared = gf.segment(0, gf.size()).get();
      set_generalized_forces(gf_shared);
    }

    /// Adds a generalized force to the body
    virtual void add_generalized_force(const Ravelin::SharedVectorNd& gf) = 0;

    /// Adds a generalized force to the body
    virtual void add_generalized_force(const Ravelin::VectorNd& gf)
    {
      const Ravelin::SharedVectorNd gf_shared = gf.segment(0, gf.size()).get();
      add_generalized_force(gf_shared);
    }

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(const Ravelin::SharedVectorNd& gj) = 0;

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(const Ravelin::VectorNd& gj)
    {
      const Ravelin::SharedVectorNd gj_shared = gj.segment(0, gj.size()).get();
      apply_generalized_impulse(gj_shared);
    }

    /// Gets the generalized coordinates of this body
    virtual Ravelin::SharedVectorNd& get_generalized_coordinates(GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gc) = 0;

    /// Gets the generalized coordinates of this body
    virtual Ravelin::VectorNd& get_generalized_coordinates(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc)
    {
      const unsigned NGC = num_generalized_coordinates(gctype);
      gc.resize(NGC);
      Ravelin::SharedVectorNd gc_shared = gc.segment(0, gc.size());
      get_generalized_coordinates(gctype, gc_shared);
      return gc;
    }

    /// Gets the generalized velocity of this body
    virtual Ravelin::SharedVectorNd& get_generalized_velocity(GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gv) = 0;

    /// Gets the generalized velocity of this body
    virtual Ravelin::VectorNd& get_generalized_velocity(GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv)
    {
      const unsigned NGC = num_generalized_coordinates(gctype);
      gv.resize(NGC);
      Ravelin::SharedVectorNd gv_shared = gv.segment(0, gv.size());
      get_generalized_velocity(gctype, gv_shared);
      return gv;
    }

    /// Gets the generalized velocity of this body
    virtual Ravelin::SharedVectorNd& get_generalized_acceleration(Ravelin::SharedVectorNd& ga) = 0;

    /// Gets the generalized acceleration of this body
    virtual Ravelin::VectorNd& get_generalized_acceleration(Ravelin::VectorNd& ga) 
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      ga.resize(NGC);
      Ravelin::SharedVectorNd ga_shared = ga.segment(0, ga.size());
      get_generalized_acceleration(ga_shared);
      return ga; 
    }

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gc) = 0;

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc) 
    {
      const Ravelin::SharedVectorNd gc_shared = gc.segment(0, gc.size()).get();
      set_generalized_coordinates(gctype, gc_shared);
    }

    /// Sets the generalized velocity of this body
    /**
      * \param gv the generalized velocity
      * \note uses the current generalized coordinates
      */
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gv) = 0;

    /// Sets the generalized velocity of this body
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv) 
    {
      const Ravelin::SharedVectorNd gv_shared = gv.segment(0, gv.size()).get();
      set_generalized_velocity(gctype, gv_shared);
    }

    /// Gets the generalized inertia of this body
    Ravelin::MatrixNd& get_generalized_inertia(Ravelin::MatrixNd& M)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      M.resize(NGC, NGC);
      Ravelin::SharedMatrixNd X = M.block(0, NGC, 0, NGC);
      get_generalized_inertia(X);
      return M;
    }

    /// Gets the generalized inertia of this body
    virtual Ravelin::SharedMatrixNd& get_generalized_inertia(Ravelin::SharedMatrixNd& M) = 0;

    /// Solves using the inverse generalized inertia
    virtual Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X) = 0;

    /// Solves using the inverse generalized inertia
    virtual Ravelin::MatrixNd& solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedMatrixNd Bx = B.block(0, B.rows(), 0, B.columns()).get();
      X.resize(NGC, B.columns());
      Ravelin::SharedMatrixNd Xx = X.block(0, X.rows(), 0, X.columns());
      solve_generalized_inertia(Bx, Xx);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::MatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::MatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      X.resize(NGC, B.columns());
      Ravelin::SharedMatrixNd Xx = X.block(0, X.rows(), 0, X.columns());
      solve_generalized_inertia(B, Xx);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::SharedMatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedMatrixNd Bx = B.block(0, B.rows(), 0, B.columns()).get();
      solve_generalized_inertia(Bx, X);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::VectorNd& x)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedVectorNd bx = b.segment(0, b.rows()).get();
      x.resize(NGC);
      Ravelin::SharedVectorNd xx = x.segment(0, x.rows());
      solve_generalized_inertia(bx, xx);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::VectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::VectorNd& x)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      x.resize(NGC);
      Ravelin::SharedVectorNd xx = x.segment(0, x.rows());
      solve_generalized_inertia(b, xx);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::VectorNd& b, Ravelin::SharedVectorNd& x)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedVectorNd bx = b.segment(0, b.rows()).get();
      solve_generalized_inertia(bx, x);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::SharedVectorNd& x) = 0;

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual Ravelin::MatrixNd& transpose_solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedMatrixNd Bx = B.block(0, B.rows(), 0, B.columns()).get();
      X.resize(NGC, B.rows());
      Ravelin::SharedMatrixNd Xx = X.block(0, X.rows(), 0, X.columns());
      transpose_solve_generalized_inertia(Bx, Xx);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual Ravelin::MatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::MatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      X.resize(NGC, B.rows());
      Ravelin::SharedMatrixNd Xx = X.block(0, X.rows(), 0, X.columns());
      transpose_solve_generalized_inertia(B, Xx);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::MatrixNd& B, Ravelin::SharedMatrixNd& X)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      const Ravelin::SharedMatrixNd Bx = B.block(0, B.rows(), 0, B.columns()).get();
      transpose_solve_generalized_inertia(Bx, X);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::SharedVectorNd& get_generalized_forces(Ravelin::SharedVectorNd& f) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& get_generalized_forces(Ravelin::VectorNd& f)
    {
      f.resize(num_generalized_coordinates(DynamicBody::eSpatial));
      Ravelin::SharedVectorNd f_shared = f.segment(0, f.size());
      get_generalized_forces(f_shared);
      return f;
    }

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the force 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::SharedVectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::SharedVectorNd& gf) = 0;

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the force 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual Ravelin::VectorNd& convert_to_generalized_force(SingleBodyPtr body, const Ravelin::SForced& w, Ravelin::VectorNd& gf)
    {
      const unsigned NGC = num_generalized_coordinates(DynamicBody::eSpatial);
      gf.resize(NGC);
      Ravelin::SharedVectorNd gf_shared = gf.segment(0, gf.size());
      convert_to_generalized_force(body, w, gf_shared);
      return gf;
    }

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

    /// Gets whether this body is kinematically updated (rather than having its dynamics integrated); default is false
    virtual bool get_kinematic() const { return _kinematic_update; }

    /// Sets whether this body is kinematically updated (rather than having its dynamics integrated); default is false
    virtual void set_kinematic(bool flag) { _kinematic_update = flag; }

    /// Prepares to compute the derivative of the body (sustained constraints) 
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) = 0;

    /// Prepares to compute the derivative of the body
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) = 0;

    /// Computes the derivative of the body
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx) = 0;

    /// Computes the derivative of the body without throwing any exceptions
    virtual void ode_noexcept(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data, Ravelin::SharedVectorNd& dx) = 0;

    /// Resets limit estimates for continuous constraint detection
    virtual void reset_limit_estimates() = 0;

  protected:

    /// The computation frame type
    ReferenceFrameType _rftype;

    /// Kinematic update flag
    bool _kinematic_update;

    /// Temporaries for use with integration
    Ravelin::VectorNd gc, gv, gcgv, xp, xv, xa;

  private:

    /// Set of recurrent forces applied to this body
    std::list<RecurrentForcePtr> _rfs;
}; // end class

} // end namespace

#endif
