/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_EVENT_H
#define _MOBY_EVENT_H

#include <iostream>
#include <list>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SForced.h>
#include <Ravelin/DynamicBodyd.h>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/ContactParameters.h>

namespace osg { class Node; }

namespace Moby {

class CollisionGeometry;

/// Container class for describing a unilateral constraint in the simulation        
class UnilateralConstraint 
{
  public:
    enum UnilateralConstraintType { eNone, eContact, eLimit };
    enum UnilateralConstraintClass { ePositive, eZero, eNegative };
    enum Compliance { eRigid, eCompliant};
    UnilateralConstraint();
    UnilateralConstraint(const UnilateralConstraint& e) { _contact_frame = boost::shared_ptr<Ravelin::Pose3d>(new Ravelin::Pose3d); *this = e; }
    static void determine_connected_constraints(const std::vector<UnilateralConstraint>& constraints, const std::vector<JointPtr>& implicit_joints, std::list<std::list<UnilateralConstraint*> >& groups, std::list<std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > >& remaining_islands);
    static void remove_inactive_groups(std::list<std::list<UnilateralConstraint*> >& groups);
    UnilateralConstraint& operator=(const UnilateralConstraint& e);
    double calc_contact_vel(const Ravelin::Vector3d& v) const;
    double calc_constraint_vel() const;
    double calc_constraint_tol() const;
    UnilateralConstraintClass determine_constraint_class() const;
    bool is_impacting() const { return determine_constraint_class() == eNegative; }
    bool is_resting() const { return determine_constraint_class() == eZero; }
    bool is_separating() const { return determine_constraint_class() == ePositive; }
    void set_contact_parameters(const ContactParameters& cparams);
    void determine_contact_tangents();
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return GLOBAL; }
    void compute_constraint_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void compute_cross_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M) const;

    template <class OutputIterator>
    OutputIterator get_super_bodies(OutputIterator begin) const;

    /// The type of constraint
    UnilateralConstraintType constraint_type;

    /// Compliance of the constraint
    Compliance compliance;

    /// The joint at which the limit is reached (for limit constraints)
    JointPtr limit_joint;

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

    /// The collision geometries involved (for contact constraints)
    CollisionGeometryPtr contact_geom1, contact_geom2;

    /// The coefficient of Coulomb friction (for contact constraints)
    double contact_mu_coulomb;

    /// The coefficient of viscous friction (for contact constraints)
    double contact_mu_viscous;

    /// Penalty Method Depth Penalty
    double contact_penalty_Kp;
    
    /// Penalty Method Interpenetration Speed
    double contact_penalty_Kv;
    
    /// The coefficient of restitution (for contact constraints)
    double contact_epsilon;

    /// The number of friction directions >= 4 (for contact constraints)
    unsigned contact_NK;

    osg::Node* to_visualization_data() const;

    /// Tolerance for the constraint (users never need to modify this)
    double tol;

    void write_vrml(const std::string& filename, double sphere_radius = 0.1, double normal_length = 1.0) const;

  private:
    // structure for comparing pairs of doubles
    struct DblComp
    {
      bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) const
      {
        return (a.first < b.first - NEAR_ZERO && a.second < b.second - NEAR_ZERO);
      }
    };

    // static variables
    boost::shared_ptr<Ravelin::Pose3d> _contact_frame;
    static Ravelin::MatrixNd JJ, J, Jx, Jy, J1, J2, dJ1, dJ2, workM1, workM2;
    static Ravelin::VectorNd v, workv, workv2;

    void compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M, boost::shared_ptr<Ravelin::DynamicBodyd> su) const;
    void compute_cross_contact_contact_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M, boost::shared_ptr<Ravelin::DynamicBodyd> su, const Ravelin::MatrixNd& J) const;
    void compute_cross_contact_limit_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_contact_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_limit_constraint_data(const UnilateralConstraint& e, Ravelin::MatrixNd& M) const;
    static bool is_linked(const UnilateralConstraint& e1, const UnilateralConstraint& e2);
    unsigned get_super_bodies(boost::shared_ptr<Ravelin::DynamicBodyd>& sb1, boost::shared_ptr<Ravelin::DynamicBodyd>& sb2) const;

    template <class BidirectionalIterator>
    static void insertion_sort(BidirectionalIterator begin, BidirectionalIterator end);
}; // end class

std::ostream& operator<<(std::ostream& out, const UnilateralConstraint& e);

#include "UnilateralConstraint.inl"

} // end namespace Moby

#endif

