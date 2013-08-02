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

/// Convex optimization data for computing forward dynamics
struct ABFwdDynOptData
{
  unsigned N_IMPLICIT_DOF; // # degrees-of-freedom for implicit joints
  unsigned N_EXPLICIT_CONSTRAINT_EQNS;
  unsigned N_JOINT_DOF;    // total # of joint degrees-of-freedom
  unsigned N_LOOPS;        // # of kinematic loops

  Ravelin::VectorNd z;               // the homogeneous solution
  Ravelin::MatrixNd R;               // the nullspace
  Ravelin::MatrixNd G;              // the quadratic objective function
  Ravelin::MatrixNd Dx;              // the movement Jacobian for explicit joints
  Ravelin::VectorNd c;               // the linear objective function
  Ravelin::VectorNd fext;            // external applied force
  std::vector<unsigned> true_indices; // maps from [implicit; explicit] indices
                                      // to true (body) joint indices
  std::vector<unsigned> loop_indices; // the loop that each joint belongs to
  std::vector<Ravelin::MatrixNd> Z;   // the Z matrices (motion for joints after loops)
  std::vector<Ravelin::MatrixNd> Zd;  // Z matrices to be multiplied by delta
  std::vector<Ravelin::MatrixNd> Z1d; // Z matrices to be multiplied by (1-delta)
  std::vector<double> mu_c;       // squared Coulomb joint friction coefficients
  std::vector<double> visc;       // the viscous joint friction forces

  // for speeding up computations
  Ravelin::MatrixNd Rff;             // components of R corresponding to implicit frict
  Ravelin::MatrixNd DxTRbetax;       // Dx' * components of R corresponding to expl. fric
  Ravelin::VectorNd zff;             // components of z corresponding to implicit frict
  Ravelin::VectorNd zbetax;          // components of z corresponding to explicit frict
};
 
/// Abstract class for articulated bodies
class ArticulatedBody : public DynamicBody
{
  public:
    ArticulatedBody();
    virtual ~ArticulatedBody() {}
    virtual bool is_floating_base() const = 0;
    virtual RigidBodyPtr get_base_link() const = 0;
    unsigned num_constraint_eqns_implicit() const;
    unsigned num_constraint_eqns_explicit() const;
    virtual void rotate(const Ravelin::Quatd& q);
    virtual void translate(const Ravelin::Origin3d& o);
    virtual double calc_kinetic_energy();
    virtual void update_visualization();
    RigidBodyPtr find_link(const std::string& id) const; 
    JointPtr find_joint(const std::string& id) const; 
    void get_adjacent_links(std::list<sorted_pair<RigidBodyPtr> >& links) const;
    virtual void set_links(const std::vector<RigidBodyPtr>& links);
    virtual void set_joints(const std::vector<JointPtr>& joints);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual unsigned num_joint_dof() const;
    void find_loops(std::vector<unsigned>& loop_indices, std::vector<std::vector<unsigned> >& loop_links) const;
    void compute_Z_matrices(const std::vector<unsigned>& loop_indices, const std::vector<std::vector<unsigned> >& loop_links, std::vector<Ravelin::MatrixNd>& Zd, std::vector<Ravelin::MatrixNd>& Z1d, std::vector<Ravelin::MatrixNd>& Z) const;
    virtual std::vector<Ravelin::SVelocityd>& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> frame, DynamicBodyPtr body, std::vector<Ravelin::SVelocityd>& J);

    /// Gets the number of degrees-of-freedom permitted by implicit constraints
    virtual unsigned num_joint_dof_implicit() const = 0;

    /// Gets the number of degrees-of-freedom permitted by explicit constraints
    virtual unsigned num_joint_dof_explicit() const = 0;

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

    /// Use the advanced (but relatively slow) full friction model?
    bool use_advanced_friction_model;

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

    Ravelin::MatrixNd& determine_F(unsigned link_idx, boost::shared_ptr<const Ravelin::Pose3d> Tf, const std::vector<unsigned>& loop_indices, Ravelin::MatrixNd& F) const;
    static double calc_fwd_dyn_f0(const Ravelin::VectorNd& x, void* data);
    static void calc_fwd_dyn_fx(const Ravelin::VectorNd& x, Ravelin::VectorNd& fc, void* data);
    static void calc_fwd_dyn_grad0(const Ravelin::VectorNd& x, Ravelin::VectorNd& grad, void* data);
    static void calc_fwd_dyn_cJac(const Ravelin::VectorNd& x, Ravelin::MatrixNd& J, void* data);
    static void calc_fwd_dyn_hess(const Ravelin::VectorNd& x, double objscal, const Ravelin::VectorNd& lambda, const Ravelin::VectorNd& nu, Ravelin::MatrixNd& H, void* data);
    void calc_joint_constraint_forces(const std::vector<unsigned>& loop_indices, const Ravelin::VectorNd& delta, const std::vector<Ravelin::MatrixNd>& Zd, const std::vector<Ravelin::MatrixNd>& Z1d, const std::vector<Ravelin::MatrixNd>& Z, const Ravelin::VectorNd& ff) const;

    /// The set of links for this articulated body
    std::vector<RigidBodyPtr> _links;

    /// The set of joints for this articulated body
    std::vector<JointPtr> _joints;

  private:
    // temporary variables
    Ravelin::VectorNd _dq;

    virtual double get_aspeed();
/*    Ravelin::SForced transform_force(RigidBodyPtr link, const Ravelin::Vector3& x) const;
*/
    static void objective_grad(const Ravelin::VectorNd& x, void* data, Ravelin::VectorNd& g);
}; // end class

#include "ArticulatedBody.inl"

} // end namespace Moby

#endif
