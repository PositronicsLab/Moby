/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _JOINT_H
#define _JOINT_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Moby/Base.h>
#include <Moby/Vector3.h>
#include <Moby/SMatrix6N.h>
#include <Moby/DynamicBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Visualizable.h>

namespace Moby {

class VisualizationData;

/// Defines a joint used in articulated rigid body dynamics simulation
/**
 * \todo implement a rest position for q?
 */
class Joint : public Visualizable
{
  public:
    enum ConstraintType { eUnknown, eImplicit, eExplicit };
    enum DOFs { DOF_1=0, DOF_2=1, DOF_3=2, DOF_4=3, DOF_5=4, DOF_6=5 };

    Joint();
    Joint(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual ~Joint() {}
    void add_force(const VectorN& force);
    void reset_force();    
    virtual const SMatrix6N& get_spatial_axes(ReferenceFrameType type);
    virtual const SMatrix6N& get_spatial_axes_complement();
    void reset_spatial_axis();
    Vector3 get_position_global(bool use_outboard = false) const;
    VectorN get_scaled_force() { VectorN f; return get_scaled_force(f); }
    VectorN& get_scaled_force(VectorN& f);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_link(RigidBodyPtr link);
    virtual void set_outboard_link(RigidBodyPtr link);
    virtual void update_spatial_axes();
    SMatrix6N& get_spatial_constraints(ReferenceFrameType rftype, SMatrix6N& s);
    ConstraintType get_constraint_type() const { return _constraint_type; }
    void evaluate_constraints_dot(Real C[6]);
    virtual void determine_Q_dot();

    /// Sets whether this constraint is implicit or explicit (or unknown)
    void set_constraint_type(ConstraintType type) { _constraint_type = type; }

    /// Gets whether this joint is in a singular configuration
    /**
     * \note only used by reduced-coordinate articulated bodies
     */
    virtual bool is_singular_config() const = 0;

    /// Gets the shared pointer to this joint
    JointPtr get_this() { return boost::dynamic_pointer_cast<Joint>(shared_from_this()); }

    /// Gets the constant shared pointer to this joint
    JointConstPtr get_this() const { return boost::dynamic_pointer_cast<const Joint>(shared_from_this()); }

    /// Gets the number of constraint equations for this joint
    virtual unsigned num_constraint_eqns() const { return 6 - num_dof(); }

    /// Evaluates the constraint equations for this joint
    /**
     * When the joint constraints are exactly satisfied, the result will be
     * the zero vector.
     * \param C a vector of size num_constraint_eqns(); contains the evaluation
     *        (on return)
     * \note only used by maximal-coordinate articulated bodies
     */
    virtual void evaluate_constraints(Real C[]) = 0;

    /// Computes the constraint Jacobian for this joint with respect to the given body
    /**
     * \param gctype the type of generalized coordinates 
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    void calc_constraint_jacobian(DynamicBody::GeneralizedCoordinateType gctype, RigidBodyPtr body, unsigned index, Real Cq[]) 
    {
      if (gctype == DynamicBody::eAxisAngle)
      {
        // compute the constraint Jacobian using rodrigues parameters
        Real Cq2[7];
        calc_constraint_jacobian_rodrigues(body, index, Cq2);

        // convert to axis-angle representation
        Cq[0] = Cq2[0];
        Cq[1] = Cq2[1];
        Cq[2] = Cq2[2];
        Vector3 ang = body->get_orientation().G_mult(Cq2[3], Cq2[4], Cq2[5], Cq2[6]) * (Real) 0.5;
        Cq[3] = ang[0];
        Cq[4] = ang[1];
        Cq[5] = ang[2];
      }
      else
        calc_constraint_jacobian_rodrigues(body, index, Cq);
    }

    /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body
    /**
     * \param gctype the type of generalized coordinates 
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    void calc_constraint_jacobian_dot(DynamicBody::GeneralizedCoordinateType gctype, RigidBodyPtr body, unsigned index, Real Cq[]) 
    {
      if (gctype == DynamicBody::eAxisAngle)
      {
        // compute the constraint Jacobian using rodrigues parameters
        Real Cq2[7];
        calc_constraint_jacobian_dot_rodrigues(body, index, Cq2);

        // convert to axis-angle representation
        Cq[0] = Cq2[0];
        Cq[1] = Cq2[1];
        Cq[2] = Cq2[2];
        Vector3 ang = body->get_orientation().G_mult(Cq2[3], Cq2[4], Cq2[5], Cq2[6]);
        Cq[3] = ang[0];
        Cq[4] = ang[1];
        Cq[5] = ang[2];
      }
      else
        calc_constraint_jacobian_dot_rodrigues(body, index, Cq);
    }

    /// Abstract method to get the spatial axes derivatives for this joint
    /**
     * Only applicable for reduced-coordinate articulated bodies
     */
    virtual const SMatrix6N& get_spatial_axes_dot(ReferenceFrameType type) = 0;

    /// Abstract method to get the local transform for this joint
    /**
     * The local transform for the joint transforms the coordinate frame
     * attached to the joint center and aligned with the inner link frame.
     */
    virtual const Matrix4& get_transform() = 0;

    /// Abstract method to determine the value of Q (joint position) from current transforms
    virtual void determine_Q() = 0;

    /// Gets the number of degrees-of-freedom for this joint
    virtual unsigned num_dof() const = 0;
  
    /// Gets the inboard link for this joint
    RigidBodyPtr get_inboard_link() const { return (_inboard_link.expired()) ? RigidBodyPtr() : RigidBodyPtr(_inboard_link); }

    /// Gets the outboard link for this joint
    RigidBodyPtr get_outboard_link() const { return (_outboard_link.expired()) ? RigidBodyPtr() : RigidBodyPtr(_outboard_link); }

    /// Gets the articulated body corresponding to this body
    /**
     * \return a pointer to the articulated body, or NULL if this body is not 
     *         a link an articulated body
     */
    ArticulatedBodyPtr get_articulated_body() { return (_abody.expired()) ? ArticulatedBodyPtr() : ArticulatedBodyPtr(_abody); }

    /// Sets the articulated body corresponding to this body
    /**
     * \param body a pointer to the articulated body or NULL if this body is
     *        not a link in an articulated body
     */
    void set_articulated_body(ArticulatedBodyPtr abody) { _abody = abody; }

    /// The lower joint limit (currently unenforced)
    VectorN lolimit;

    /// The upper joint limit (currently unenforced)
    VectorN hilimit;

    /// The maximum force that can be applied to each DOF of the joint
    VectorN maxforce;

    /// The position of this joint
    VectorN q;

    /// The velocity of this joint
    VectorN qd;

    /// The acceleration of this joint
    VectorN qdd;

    /// The coefficient of restitution applied when this joint reaches a limit
    Real limit_restitution;

    /// The coulomb friction coefficient for this joint
    Real mu_fc;

    /// The viscous friction coefficient for this joint
    Real mu_fv;

    /// The total friction applied to this joint (computed by forward dynamics / impulse application)
    VectorN ff;

    /// The actuator force (user/controller sets this)
    VectorN force;

    /// Constraint forces calculated by forward dynamics
    VectorN lambda;

    /// Gets the joint index (returns UINT_MAX if not set)
    unsigned get_index() const { return _joint_idx; }

    /// Sets the joint index
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_index(unsigned index) { _joint_idx = index; }

    /// Gets the constraint index for this joint
    unsigned get_constraint_index() const { return _constraint_idx; }

    /// Sets the constraint index for this joint
    void set_constraint_index(unsigned idx) { _constraint_idx = idx; } 

    /// Sets the coordinate index for this joint
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_coord_index(unsigned index) { _coord_idx = index; }

    /// Gets the starting coordinate index for this joint
    unsigned get_coord_index() const { return _coord_idx; }

  protected:
    void calc_s_bar_from_si();

    /// Computes the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    virtual void calc_constraint_jacobian_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[]) = 0;
 
     /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    virtual void calc_constraint_jacobian_dot_rodrigues(RigidBodyPtr body, unsigned index, Real Cq[]) = 0;

    virtual const Matrix4* get_visualization_transform();

    /// Method for initializing all variables in the joint
    /**
     * This method should be called at the beginning of all constructors of
     * all derived classes.
     */
    void init_data();

    /// The spatial axes (in inner link frame) for the joint
    /**
     * Spatial axes are used in the dynamics equations for reduced-coordinate
     * articulated bodies only.
     */
    SMatrix6N _si;

    /// The complement of the spatial axes (in joint position frame) for the joint
    /**
     * Spatial axes are used in the dynamic equations for reduced-coordinate
     * articulated bodies only.
     */
    SMatrix6N _s_bar;

    /// The spatial axes in global frame for the joint
    /**
     * Spatial axes are used in the dynamics equations for reduced-coordinate
     * articulated bodies only.
     */
    SMatrix6N _s0;

    /// The stored tare value for the initial joint configuration
    /**
     * Spatial axes are used in the dynamics equations for reduced-coordinate
     * articulated bodies only.
     */
    VectorN _q_tare;
   
  private:
    ConstraintType _constraint_type;
    unsigned _joint_idx;
    unsigned _coord_idx;
    unsigned _constraint_idx;
    bool _s0_valid;
    Matrix4 _vtransform;
    boost::weak_ptr<RigidBody> _inboard_link;
    boost::weak_ptr<RigidBody> _outboard_link;
    boost::weak_ptr<ArticulatedBody> _abody;
}; // end class
} // end namespace

#endif

