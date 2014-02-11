/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _ARTICULATED_BODY_H
#define _ARTICULATED_BODY_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/SMomentumd.h>
#include <Moby/sorted_pair>
#include <Moby/Event.h>
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
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);

    /// Determines whether joint acceleration limits were exceeded on last integration
    bool joint_accel_limit_exceeded() const { return _acc_limits_exceeded; }

    /// Gets the number of degrees-of-freedom permitted by explicit constraints
    virtual unsigned num_joint_dof_explicit() const = 0;

    /// Gets the number of degrees-of-freedom permitted by implicit constraints
    virtual unsigned num_joint_dof_implicit() const = 0;

    template <class OutputIterator>
    OutputIterator find_limit_events(double dt, OutputIterator begin) const;

    /// Finds (joint) limit events
    template <class OutputIterator>
    OutputIterator find_limit_events(const Ravelin::VectorNd& q0, const Ravelin::VectorNd& q1, double dt, OutputIterator begin);

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
    virtual void compile() = 0;

    /// The set of links for this articulated body
    std::vector<RigidBodyPtr> _links;

    /// The set of joints for this articulated body
    std::vector<JointPtr> _joints;

  private:
    // joint constraint violation
    std::vector<double> _cvio;

    // lower acceleration limits for this body
    std::vector<double> _acc_limits_lo;

    // upper acceleration limits for this body
    std::vector<double> _acc_limits_hi;

    // temporary variables
    Ravelin::VectorNd _dq;

    // indicates whether acceleration bounds exceeded limits since being reset
    bool _acc_limits_exceeded;

    void check_joint_accel_limit_exceeded();
    virtual double get_aspeed();
    static Ravelin::VectorNd& ode_both(const Ravelin::VectorNd& x, double t, double dt, void* data, Ravelin::VectorNd& dx);
}; // end class

#include "ArticulatedBody.inl"

} // end namespace Moby

#endif
