/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _JOINT_H
#define _JOINT_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
#include <Moby/Base.h>
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
    void add_force(const Ravelin::VectorNd& force);
    void reset_force();    
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes();
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes_complement();
    Ravelin::Point3d get_location(bool use_outboard = false) const;
    Ravelin::VectorNd& get_scaled_force(Ravelin::VectorNd& f);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_link(RigidBodyPtr link);
    virtual void set_outboard_link(RigidBodyPtr link);
    virtual void update_spatial_axes();
    std::vector<Ravelin::Twistd>& get_spatial_constraints(ReferenceFrameType rftype, std::vector<Ravelin::Twistd>& s);
    ConstraintType get_constraint_type() const { return _constraint_type; }
    void evaluate_constraints_dot(double C[6]);
    virtual void determine_q_dot();
    void determine_q_tare();
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _F; };

    /// Gets the pose of this joint relative to the outboard link (rather than the inboard link as is standard)
    boost::shared_ptr<const Ravelin::Pose3d> get_pose_outboard() const { return _Fb; };

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
    boost::shared_ptr<const Joint> get_this() const { return boost::dynamic_pointer_cast<const Joint>(shared_from_this()); }

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
    virtual void evaluate_constraints(double C[]) = 0;

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
    void calc_constraint_jacobian(DynamicBody::GeneralizedCoordinateType gctype, RigidBodyPtr body, unsigned index, double Cq[]) 
    {
      if (gctype == DynamicBody::eSpatial)
      {
        // compute the constraint Jacobian using euler parameters
        double Cq2[7];
        calc_constraint_jacobian_euler(body, index, Cq2);

        // convert to axis-angle representation
        Cq[0] = Cq2[0];
        Cq[1] = Cq2[1];
        Cq[2] = Cq2[2];
        Ravelin::Vector3d ang = body->get_pose()->q.G_mult(Cq2[3], Cq2[4], Cq2[5], Cq2[6]) * (double) 0.5;
        Cq[3] = ang[0];
        Cq[4] = ang[1];
        Cq[5] = ang[2];
      }
      else
        calc_constraint_jacobian_euler(body, index, Cq);
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
    void calc_constraint_jacobian_dot(DynamicBody::GeneralizedCoordinateType gctype, RigidBodyPtr body, unsigned index, double Cq[]) 
    {
      if (gctype == DynamicBody::eSpatial)
      {
        // compute the constraint Jacobian using euler parameters
        double Cq2[7];
        calc_constraint_jacobian_dot_euler(body, index, Cq2);

        // convert to axis-angle representation
        Cq[0] = Cq2[0];
        Cq[1] = Cq2[1];
        Cq[2] = Cq2[2];
        Ravelin::Vector3d ang = body->get_pose()->q.G_mult(Cq2[3], Cq2[4], Cq2[5], Cq2[6]);
        Cq[3] = ang[0];
        Cq[4] = ang[1];
        Cq[5] = ang[2];
      }
      else
        calc_constraint_jacobian_dot_euler(body, index, Cq);
    }

    /// Abstract method to get the spatial axes derivatives for this joint
    /**
     * Only applicable for reduced-coordinate articulated bodies
     */
    virtual const std::vector<Ravelin::Twistd>& get_spatial_axes_dot() = 0;

    /// Abstract method to get the local transform for this joint
    /**
     * The local transform for the joint transforms the coordinate frame
     * attached to the joint center and aligned with the inner link frame.
     */
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose() = 0;

    /// Abstract method to determine the value of Q (joint position) from current transforms
    virtual void determine_q(Ravelin::VectorNd& q) = 0;

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

    /// The lower joint limit 
    Ravelin::VectorNd lolimit;

    /// The upper joint limit 
    Ravelin::VectorNd hilimit;

    /// The maximum force that can be applied to each DOF of the joint
    Ravelin::VectorNd maxforce;

    /// The position of this joint
    Ravelin::VectorNd q;

    /// The velocity of this joint
    Ravelin::VectorNd qd;

    /// The acceleration of this joint
    Ravelin::VectorNd qdd;

    /// The coefficient of restitution applied when this joint reaches a limit
    double limit_restitution;

    /// The coulomb friction coefficient for this joint
    double mu_fc;

    /// The viscous friction coefficient for this joint
    double mu_fv;

    /// The total friction applied to this joint (computed by forward dynamics / impulse application)
    Ravelin::VectorNd ff;

    /// The actuator force (user/controller sets this)
    Ravelin::VectorNd force;

    /// Constraint forces calculated by forward dynamics
    Ravelin::VectorNd lambda;

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
    void calc_s_bar_from_s();

    /// The frame induced by the joint 
    boost::shared_ptr<Ravelin::Pose3d> _Fprime;

    /// The frame of this joint
    boost::shared_ptr<Ravelin::Pose3d> _F;

    /// The frame of this joint _backward_ from the outboard link
    boost::shared_ptr<Ravelin::Pose3d> _Fb;

    /// Computes the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    virtual void calc_constraint_jacobian_euler(RigidBodyPtr body, unsigned index, double Cq[]) = 0;
 
     /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
    virtual void calc_constraint_jacobian_dot_euler(RigidBodyPtr body, unsigned index, double Cq[]) = 0;

    virtual boost::shared_ptr<const Ravelin::Pose3d> get_visualization_pose();

    /// Method for initializing all variables in the joint
    /**
     * This method should be called at the beginning of all constructors of
     * all derived classes.
     */
    void init_data();

    /// The spatial axes (in joint frame) for the joint
    /**
     * Spatial axes are used in the dynamics equations for reduced-coordinate
     * articulated bodies only.
     */
    std::vector<Ravelin::Twistd> _s;

    /// The complement of the spatial axes (in joint position frame) for the joint
    /**
     * Spatial axes are used in the dynamic equations for reduced-coordinate
     * articulated bodies only.
     */
    std::vector<Ravelin::Twistd> _s_bar;

    /// The stored "tare" value for the initial joint configuration
    /**
     * The tare value is the value that the joint assumes in the an
     * articulated body's initial configuration. This value is necessary
     * so that- when the body's joints are set to the zero vector- the body
     * re-enters the initial configuration.
     */
    Ravelin::VectorNd _q_tare;

    /// Set whether _q_tare needs to be determined
    bool _determine_q_tare;

  private:
    // working variables for calc_s_bar_from_s()
    Ravelin::MatrixNd _ns;

    ConstraintType _constraint_type;
    unsigned _joint_idx;
    unsigned _coord_idx;
    unsigned _constraint_idx;
    boost::shared_ptr<Ravelin::Pose3d> _vtransform;
    boost::weak_ptr<RigidBody> _inboard_link;
    boost::weak_ptr<RigidBody> _outboard_link;
    boost::weak_ptr<ArticulatedBody> _abody;
}; // end class
} // end namespace

#endif

