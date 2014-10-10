/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _ARTICULATED_BODY_H
#define _ARTICULATED_BODY_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/SMomentumd.h>
#include <Moby/sorted_pair>
#include <Moby/UnilateralConstraint.h>
#include <Moby/DynamicBody.h>
#include <Moby/Joint.h>

namespace Moby {

class RigidBody;

/// Abstract class for articulated bodies
class ArticulatedBody : public DynamicBody
{
  public:
    ArticulatedBody();
    virtual ~ArticulatedBody() {}
    virtual bool is_floating_base() const = 0;
    virtual RigidBodyPtr get_base_link() const = 0;
    unsigned num_constraint_eqns_explicit() const;
    unsigned num_constraint_eqns_implicit() const;
    virtual void rotate(const Ravelin::Quatd& q);
    virtual void translate(const Ravelin::Origin3d& o);
    virtual double calc_kinetic_energy();
    virtual void update_visualization();
    RigidBodyPtr find_link(const std::string& id) const; 
    JointPtr find_joint(const std::string& id) const; 
    void get_adjacent_links(std::list<sorted_pair<RigidBodyPtr> >& links) const;
    virtual void set_links_and_joints(const std::vector<RigidBodyPtr>& links, const std::vector<JointPtr>& joints);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_joint_dof() const;
    void find_loops(std::vector<unsigned>& loop_indices, std::vector<std::vector<unsigned> >& loop_links) const;
    virtual Ravelin::MatrixNd& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J);
    virtual Ravelin::MatrixNd& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, Ravelin::MatrixNd& J);
    void update_joint_constraint_violations();
    bool is_joint_constraint_violated() const;
    double calc_CA_time_for_joints() const;
    virtual void ode_noexcept(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void reset_limit_estimates();
    virtual bool limit_estimates_exceeded() const;
    double find_next_joint_limit_time() const;
    void update_joint_vel_limits();
    virtual void validate_limit_estimates();

    /// Gets the number of degrees-of-freedom permitted by explicit constraints
    virtual unsigned num_joint_dof_explicit() const = 0;

    /// Gets the number of degrees-of-freedom permitted by implicit constraints
    virtual unsigned num_joint_dof_implicit() const = 0;

    /// Finds (joint) limit constraints 
    template <class OutputIterator>
    OutputIterator find_limit_constraints(OutputIterator begin) const;

    /// Gets the set of links
    virtual const std::vector<RigidBodyPtr>& get_links() const { return _links; }

    /// Gets the set of joints
    virtual const std::vector<JointPtr>& get_joints() const { return _joints; }

    /// Gets shared pointer to this object as type ArticulatedBody
    ArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<ArticulatedBody>(shared_from_this()); }

    /// Gets shared pointer to this object as type const ArticulateBody
    boost::shared_ptr<const ArticulatedBody> get_this() const { return boost::dynamic_pointer_cast<const ArticulatedBody>(shared_from_this()); }

    /// Abstract method for applying an impulse to this articulated body
    /**
     * \param w the impulsive force 
     * \param link link in the articulated body where the impulse is applied
     */
    virtual void apply_impulse(const Ravelin::SMomentumd& w, RigidBodyPtr link) = 0;
      
    /// Method for resetting the force and torque accumulators on all links
    virtual void reset_accumulators() = 0;

    /// Multiplies Jc' for this body by the given vector 
    virtual Ravelin::VectorNd& transpose_Jc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) = 0; 

    /// Multiplies Jc' for this body by the given matrix
    virtual Ravelin::MatrixNd& transpose_Jc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) = 0;

    /// Multiplies Dc' for this body by the given vector 
    virtual Ravelin::VectorNd& transpose_Dc_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) = 0;

    /// Multiplies Dc' for this body by the given matrix 
    virtual Ravelin::MatrixNd& transpose_Dc_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) = 0;

    /// Multiplies Jl' for this body by the given vector 
    virtual Ravelin::VectorNd& transpose_Jl_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) = 0;

    /// Multiplies Jl' for this body by the given matrix
    virtual Ravelin::MatrixNd& transpose_Jl_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) = 0;

    /// Multiplies Dx' for this body by the given vector 
    virtual Ravelin::VectorNd& transpose_Dx_mult(const Ravelin::VectorNd& v, Ravelin::VectorNd& result) = 0;

    /// Multiplies Dx' for this body by the given matrix 
    virtual Ravelin::MatrixNd& transpose_Dx_mult(const Ravelin::MatrixNd& m, Ravelin::MatrixNd& result) = 0; 

  protected:
    /// Vector for processing links
    std::vector<unsigned> _processed;

    /// Method for "compiling" the body
    /**
     * Compilation is necessary anytime the structure or topolgy of the body
     * changes.  Generally, this will only be necessary once - after a call
     * to set_links() or set_joints().
     */
    virtual void compile();

    /// The set of links for this articulated body
    std::vector<RigidBodyPtr> _links;

    /// The set of joints for this articulated body
    std::vector<JointPtr> _joints;

    // the limit bound expansion for updating joint velocity limit estimates
    double limit_bound_expansion;

  private:
    // joint constraint violation
    std::vector<double> _cvio;

    // joint velocity tolerances (for joints at constraints)
    std::vector<double> _cvel_vio;

    // lower velocity limits for this body
    std::vector<double> _vel_limits_lo;

    // upper acceleration limits for this body
    std::vector<double> _vel_limits_hi;

    // temporary variables
    Ravelin::VectorNd _dq;

    // indicates whether velocity bounds exceeded limits since being reset
    bool _vel_limit_exceeded;

    void check_joint_vel_limit_exceeded_and_update();
    virtual double get_aspeed();
}; // end class

#include "ArticulatedBody.inl"

} // end namespace Moby

#endif
