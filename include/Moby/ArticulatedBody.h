/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _ARTICULATED_BODY_H
#define _ARTICULATED_BODY_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/sorted_pair>
#include <Moby/Event.h>
#include <Moby/DynamicBody.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/SVector6.h>
#include <Moby/Optimization.h>
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

  VectorN z;               // the homogeneous solution
  MatrixN R;               // the nullspace
  MatrixNN G;              // the quadratic objective function
  MatrixN Dx;              // the movement Jacobian for explicit joints
  VectorN c;               // the linear objective function
  VectorN fext;            // external applied force
  std::vector<unsigned> true_indices; // maps from [implicit; explicit] indices
                                      // to true (body) joint indices
  std::vector<unsigned> loop_indices; // the loop that each joint belongs to
  std::vector<MatrixN> Z;   // the Z matrices (motion for joints after loops)
  std::vector<MatrixN> Zd;  // Z matrices to be multiplied by delta
  std::vector<MatrixN> Z1d; // Z matrices to be multiplied by (1-delta)
  std::vector<Real> mu_c;       // squared Coulomb joint friction coefficients
  std::vector<Real> visc;       // the viscous joint friction forces

  // for speeding up computations
  MatrixN Rff;             // components of R corresponding to implicit frict
  MatrixN DxTRbetax;       // Dx' * components of R corresponding to expl. fric
  VectorN zff;             // components of z corresponding to implicit frict
  VectorN zbetax;          // components of z corresponding to explicit frict
};
 
/// Abstract class for articulated bodies
class ArticulatedBody : public DynamicBody
{
  public:
    ArticulatedBody();
    virtual ~ArticulatedBody() {}
    unsigned num_constraint_eqns_implicit() const;
    unsigned num_constraint_eqns_explicit() const;
    virtual void transform(const Matrix4& T);
    virtual Real calc_kinetic_energy() const;
    virtual void update_visualization();
    RigidBodyPtr find_link(const std::string& id) const; 
    JointPtr find_joint(const std::string& id) const; 
    void get_adjacent_links(std::list<sorted_pair<RigidBodyPtr> >& links) const;
    virtual void set_links(const std::vector<RigidBodyPtr>& links);
    virtual void set_joints(const std::vector<JointPtr>& joints);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual unsigned num_joint_dof() const;
    void get_constraint_events(std::vector<Event>& cevents) const;
    void find_loops(std::vector<unsigned>& loop_indices, std::vector<std::vector<unsigned> >& loop_links) const;
    void compute_Z_matrices(const std::vector<unsigned>& loop_indices, const std::vector<std::vector<unsigned> >& loop_links, std::vector<MatrixN>& Zd, std::vector<MatrixN>& Z1d, std::vector<MatrixN>& Z) const;

    /// Gets the number of degrees-of-freedom permitted by implicit constraints
    virtual unsigned num_joint_dof_implicit() const = 0;

    /// Gets the number of degrees-of-freedom permitted by explicit constraints
    virtual unsigned num_joint_dof_explicit() const = 0;

    /// Finds (joint) limit events
    template <class OutputIterator>
    OutputIterator find_limit_events(const VectorN& q0, const VectorN& q1, Real dt, OutputIterator begin);

    /// Gets the set of links
    virtual const std::vector<RigidBodyPtr>& get_links() const { return _links; }

    /// Gets the set of joints
    virtual const std::vector<JointPtr>& get_joints() const { return _joints; }

    /// Returns true if any of the link positions have changed
    bool positions_invalidated() const { return !_positions_valid; }

    /// Returns true if any of the link velocities have changed
    bool velocities_invalidated() const { return !_velocities_valid; }

    /// Invalidates the link positions (manually)
    virtual void invalidate_positions() { _positions_valid = false; }

    /// Invalidates the link velocities (manually)
    virtual void invalidate_velocities() { _velocities_valid = false; }

    /// Validates the link positions (manually)
    void validate_positions() { _positions_valid = true; }

    /// Validates the link velocities (manually)
    void validate_velocities() { _velocities_valid = true; }

    /// Gets shared pointer to this object as type ArticulatedBody
    ArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<ArticulatedBody>(shared_from_this()); }

    /// Gets shared pointer to this object as type const ArticulateBody
    ArticulatedBodyConstPtr get_this() const { return boost::dynamic_pointer_cast<const ArticulatedBody>(shared_from_this()); }

    /// Abstract method for applying an impulse to this articulated body
    /**
     * \param j the linear component of an impulse
     * \param k the angular component of an impulse
     * \param p the point at which to apply the impulse
     * \param link link in the articulated body where the impulse is applied
     */
    virtual void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& p, RigidBodyPtr link) = 0;
      
    /// Method for resetting the force and torque accumulators on all links
    virtual void reset_accumulators() = 0;

    /// Use the advanced (but relatively slow) full friction model?
    bool use_advanced_friction_model;

    /// Multiplies Jc' for this body by the given vector 
    virtual VectorN& transpose_Jc_mult(const VectorN& v, VectorN& result) = 0; 

    /// Multiplies Jc' for this body by the given matrix
    virtual MatrixN& transpose_Jc_mult(const MatrixN& m, MatrixN& result) = 0;

    /// Multiplies Dc' for this body by the given vector 
    virtual VectorN& transpose_Dc_mult(const VectorN& v, VectorN& result) = 0;

    /// Multiplies Dc' for this body by the given matrix 
    virtual MatrixN& transpose_Dc_mult(const MatrixN& m, MatrixN& result) = 0;

    /// Multiplies Jl' for this body by the given vector 
    virtual VectorN& transpose_Jl_mult(const VectorN& v, VectorN& result) = 0;

    /// Multiplies Jl' for this body by the given matrix
    virtual MatrixN& transpose_Jl_mult(const MatrixN& m, MatrixN& result) = 0;

    /// Multiplies Dx' for this body by the given vector 
    virtual VectorN& transpose_Dx_mult(const VectorN& v, VectorN& result) = 0;

    /// Multiplies Dx' for this body by the given matrix 
    virtual MatrixN& transpose_Dx_mult(const MatrixN& m, MatrixN& result) = 0; 

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

    MatrixN& determine_F(unsigned link_idx, const Matrix4& Tf, const std::vector<unsigned>& loop_indices, MatrixN& F) const;
    static Real calc_fwd_dyn_f0(const VectorN& x, void* data);
    static void calc_fwd_dyn_fx(const VectorN& x, VectorN& fc, void* data);
    static void calc_fwd_dyn_grad0(const VectorN& x, VectorN& grad, void* data);
    static void calc_fwd_dyn_cJac(const VectorN& x, MatrixN& J, void* data);
    static void calc_fwd_dyn_hess(const VectorN& x, Real objscal, const VectorN& lambda, const VectorN& nu, MatrixNN& H, void* data);
    void calc_joint_constraint_forces(const std::vector<unsigned>& loop_indices, const VectorN& delta, const std::vector<MatrixN>& Zd, const std::vector<MatrixN>& Z1d, const std::vector<MatrixN>& Z, const VectorN& ff) const;

    /// The set of links for this articulated body
    std::vector<RigidBodyPtr> _links;

    /// The set of joints for this articulated body
    std::vector<JointPtr> _joints;

  private:
    virtual Real get_aspeed() const;
    SVector6 transform_force(RigidBodyPtr link, const Vector3& x) const;
    static void objective_grad(const VectorN& x, void* data, VectorN& g);

    bool _positions_valid;
    bool _velocities_valid;
}; // end class

#include "ArticulatedBody.inl"

} // end namespace Moby

#endif
