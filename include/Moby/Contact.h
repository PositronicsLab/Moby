/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_CONTACT_H
#define _MOBY_CONTACT_H

#include <iostream>
#include <list>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SForced.h>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/ContactParameters.h>

namespace osg { class Node; }

namespace Moby {

class CollisionGeometry;

/// Container class for describing a contact in the simulation        
class Contact 
{
  public:
    enum ContactType { eNone, eContact };
    enum CoulombFrictionType { eUndetermined, eSticking, eSlipping };
    enum DistanceClass { eNegative, eZero, ePositive };
    Contact();
    Contact(const Contact& e) { _contact_frame = boost::shared_ptr<Ravelin::Pose3d>(new Ravelin::Pose3d); *this = e; }
    static void determine_connected_contacts(const std::vector<Contact>& contacts, std::list<std::list<Contact*> >& groups);
    static void remove_inactive_groups(std::list<std::list<Contact*> >& groups);
    Contact& operator=(const Contact& e);
    double calc_contact_vel() const;
    double calc_contact_accel() const;
    double calc_contact_tol() const;
    DistanceClass determine_distance_class() const;
    bool distance_is_negative() const { return determine_distance_class() == eNegative; }
    bool distance_is_zero() const { return determine_distance_class() == eZero; }
    bool distance_is_positive() const { return determine_distance_class() == ePositive; }
    void set_contact_parameters(const ContactParameters& cparams);
    void determine_contact_tangents();
    void update_contact_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void update_cross_contact_data(const Contact& e, Ravelin::MatrixNd& M) const;
    CoulombFrictionType get_friction_type() const { return _ftype; }
    static void determine_minimal_set(std::list<Contact*>& group);
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return GLOBAL; }

    template <class OutputIterator>
    OutputIterator get_super_bodies(OutputIterator begin) const;

    /// The type of contact
    ContactType contact_type;

    /// Constraint [normal] impulse magnitude (for constraint contacts)
    Ravelin::VectorNd constraint_nforce;

    /// Constraint [friction] impulse magnitude (for constraint contacts)
    Ravelin::VectorNd constraint_fforce;

    /// The joint (for constraint contacts)
    JointPtr constraint_joint;

    /// The point contact 
    Point3d contact_point;
    
    /// The vector pointing outward from the contact on the first body, in world coordinates
    Ravelin::Vector3d contact_normal;  

    /// The first tangent direction to the contact normal 
    Ravelin::Vector3d contact_tan1;

    /// The second tangent direction to the contact normal
    Ravelin::Vector3d contact_tan2;

    /// The time derivative of the contact normal 
    Ravelin::Vector3d contact_normal_dot;  

    /// The time derivative of the first tangent direction 
    Ravelin::Vector3d contact_tan1_dot;

    /// The time derivative of the second tangent direction 
    Ravelin::Vector3d contact_tan2_dot;

    /// Force that has been applied 
    /**
     * Force applied to the body of the first geometry; 
     * the reverse of this force is applied to the body of the second 
     * geometry.
     */
    Ravelin::SMomentumd contact_force;

    /// The collision geometries involved 
    CollisionGeometryPtr contact_geom1, contact_geom2;

    /// The coefficient of Coulomb friction 
    double contact_mu_coulomb;

    /// The number of friction directions >= 4 (for contact events)
    unsigned contact_NK;

    osg::Node* to_visualization_data() const;

    /// Tolerance for the contact (users never need to modify this)
    double tol;

    /// The sticking tolerance (users never need to modify this)
    double stick_tol;

    void write_vrml(const std::string& filename, double sphere_radius = 0.1, double normal_length = 1.0) const;
    void compute_contact_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void compute_cross_contact_data(const Contact& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_data(const Contact& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_data(const Contact& e, Ravelin::MatrixNd& M, DynamicBodyPtr su) const;
    void compute_cross_contact_contact_data(const Contact& e, Ravelin::MatrixNd& M, DynamicBodyPtr su, const Ravelin::MatrixNd& J) const;

  private:
    // static variables
    boost::shared_ptr<Ravelin::Pose3d> _contact_frame;
    static Ravelin::MatrixNd J1, J2, workM1, workM2, dJ1, dJ2, Jx, Jy;
    static std::vector<Ravelin::SVelocityd> vel;
    static Ravelin::VectorNd v, workv, workv2;

    // the type of friction contact (sticking or slipping)
    CoulombFrictionType _ftype;

    void compute_dotv_data(Ravelin::VectorNd& q) const;
    static bool is_linked(const Contact& e1, const Contact& e2);
    unsigned get_super_bodies(DynamicBodyPtr& sb1, DynamicBodyPtr& sb2) const;
    static void determine_convex_set(std::list<Contact*>& group);
    static bool is_contact_manifold_2D(const std::list<Contact*>& group);
}; // end class

std::ostream& operator<<(std::ostream& out, const Contact& e);

#include "Contact.inl"

} // end namespace Moby

#endif

