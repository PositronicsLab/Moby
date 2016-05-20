/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_CONSTRAINT_H
#define _MOBY_CONSTRAINT_H

#include <iostream>
#include <list>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SForced.h>
#include <Ravelin/DynamicBodyd.h>
#include <Ravelin/SingleBodyd.h>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/Joint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/ContactParameters.h>

namespace osg { class Node; }

namespace Moby {

class CollisionGeometry;

/// Class for describing unilateral constraints and implicit bilateral constraints in the simulation        
class Constraint 
{
  public:
    enum ConstraintType { eNone, eContact, eLimit, eInverseDynamics, eSpringDamper, eImplicitJoint };
    enum ConstraintEquationType { eInequality, eEquality, eComplementarity };
    enum ConstraintCoeffType { eImpulsiveConstraint, eVelocityConstraint };
    enum ConstraintVelocityClass { ePositive, eZero, eNegative };
    Constraint();
    Constraint(const Constraint& e) { *this = e; }
    static void determine_connected_constraints(const std::vector<Constraint>& constraints, std::list<std::pair<std::list<Constraint*>, std::list<boost::shared_ptr<Ravelin::SingleBodyd> > > >& groups, std::list<std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > >& remaining_islands);
    static void remove_inactive_groups(std::list<std::pair<std::list<Constraint*>, std::list<boost::shared_ptr<Ravelin::SingleBodyd> > > >& groups);
    Constraint& operator=(const Constraint& e);
    ConstraintCoeffType get_constraint_coeff_type(unsigned constraint_eqn_index) const;
    unsigned get_velocity_constraint_index(unsigned constraint_eqn_index) const;
    void get_impulse_constraint_coeffs(unsigned constraint_eqn_index, Ravelin::SharedVectorNd& a);
    bool slack_allowed() const;
    unsigned num_variables() const;
    unsigned num_impulsive_variables() const;
    unsigned num_constraint_equations() const;
    void apply_restitution(Ravelin::VectorNd& x);
    double get_damping(unsigned constraint_eqn_index) const;
    double get_lower_variable_limit(unsigned var_index) const;
    double get_upper_variable_limit(unsigned var_index) const;
    void apply_test_impulse(unsigned var_index);
    void apply_impulses(const Ravelin::VectorNd& x, std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, Ravelin::VectorNd>& gj);
    ConstraintEquationType get_constraint_equation_type(unsigned constraint_eqn_index) const; 
    std::vector<double>& get_constraint_coeffs(unsigned constraint_eqn_index, std::vector<double>& coeffs); 
    double get_constraint_rhs(unsigned constraint_index, double inv_dt);

    double calc_contact_vel(const Ravelin::Vector3d& v) const;
    double calc_projected_vel(unsigned var_index) const;
    double calc_constraint_vel(unsigned constraint_eqn_index) const;
    double calc_constraint_tol() const;
    ConstraintVelocityClass determine_constraint_velocity_class() const;
    bool is_impacting() const { return determine_constraint_velocity_class() == eNegative; }
    bool is_resting() const { return determine_constraint_velocity_class() == eZero; }
    bool is_separating() const { return determine_constraint_velocity_class() == ePositive; }
    void set_contact_parameters(const ContactParameters& cparams);
    void determine_contact_tangents();

    template <class OutputIterator>
    OutputIterator get_super_bodies(OutputIterator begin) const;

    /// The coefficient for when a contact is to be treated without slip
    static double NO_SLIP_COEFF;

    /// The type of constraint
    ConstraintType constraint_type;

    /// The joint for an implicit joint constraint
    JointPtr implicit_joint;

    /// The joint at which the limit is reached (for limit constraints)
    JointPtr limit_joint;

    /// The joint for inverse dynamics
    JointPtr inv_dyn_joint;

    /// The desired joint velocity for an inverse dynamics constraint 
    Ravelin::VectorNd qdot_des;

    /// Maximum allowable force for inverse dynamics application
    Ravelin::VectorNd inv_dyn_fmax;

    /// The joint for spring/damper joints
    JointPtr spring_damper_joint;

    // the contact stiffness (applicable only when contact is in compliant layer)
    double contact_stiffness;

    // the contact damping  (applicable only when contact is in compliant layer)
    double contact_damping;

    // the limit stiffness  (applicable only when limit is in compliant layer)
    double limit_stiffness;

    // the limit damping  (applicable only when limit is in compliant layer)
    double limit_damping;

    /// Signed violation for this constraint
    double signed_violation;

    /// The coefficient of restitution for this limit
    double limit_epsilon;

    /// The DOF at which the limit is reached (for limit constraints)
    unsigned limit_dof;

    /// Whether the upper/lower limit is reached (for limit constraints)
    bool limit_upper;

    /// Limit impulse magnitude (for limit constraints)
    double limit_impulse;

    /// The point contact (for contact constraints)
    Point3d contact_point;

    /// The vector pointing outward from the contact on the first body, in world coordinates (for contact constraints)
    Ravelin::Vector3d contact_normal;  

    /// The first tangent direction to the contact normal
    Ravelin::Vector3d contact_tan1;

    /// The second tangent direction to the contact normal
    Ravelin::Vector3d contact_tan2;

    /// Impulse that has been applied (for contact constraints)
    /**
     * Impulse applied to the body of the first geometry; 
     * the reverse of this force / impulse is applied to the body of the second 
     * geometry.
     */
    Ravelin::SMomentumd contact_impulse;

    /// Impulse to keep implicit joints together
    Ravelin::VectorNd implicit_joint_impulse;

    /// Impulse for spring and damper
    Ravelin::VectorNd spring_damper_impulse;

    /// Impulse for inverse dynamics constraints
    Ravelin::VectorNd inv_dyn_impulse;

    /// The collision geometries involved (for contact constraints)
    CollisionGeometryPtr contact_geom1, contact_geom2;

    /// The coefficient of Coulomb friction (for contact constraints)
    double contact_mu_coulomb;

    /// The coefficient of viscous friction (for contact constraints)
    double contact_mu_viscous;

    /// The coefficient of restitution (for contact constraints)
    double contact_epsilon;

    /// The number of friction directions >= 4 (for contact constraints)
    unsigned contact_NK;

    osg::Node* to_visualization_data() const;

    /// Tolerance for the constraint (users never need to modify this)
    double tol;

    void write_vrml(const std::string& filename, double sphere_radius = 0.1, double normal_length = 1.0) const;

  private:
    double eval_spring() const;
    double calc_spring_damper_fmax() const; 

    // structure for comparing pairs of doubles
    struct DblComp
    {
      bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) const
      {
        return (a.first < b.first - NEAR_ZERO && a.second < b.second - NEAR_ZERO);
      }
    };

    static bool is_linked(const Constraint& e1, const Constraint& e2);

    template <class BidirectionalIterator>
    static void insertion_sort(BidirectionalIterator begin, BidirectionalIterator end);
}; // end class

std::ostream& operator<<(std::ostream& out, const Constraint& e);

#include "Constraint.inl"

} // end namespace Moby

#endif

