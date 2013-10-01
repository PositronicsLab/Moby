/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_EVENT_H
#define _MOBY_EVENT_H

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

/// Container class for describing an event in the simulation        
class Event 
{
  public:
    enum EventType { eNone, eContact, eLimit };
    enum EventClass { ePositive, eZero, eNegative };
    enum DerivType { eVel, eAccel };
    enum CoulombFrictionType { eUndetermined, eSlipping, eSticking }; 
    Event();
    Event(const Event& e) { _event_frame = boost::shared_ptr<Ravelin::Pose3d>(new Ravelin::Pose3d); *this = e; }
    static void determine_connected_events(const std::vector<Event>& events, std::list<std::list<Event*> >& groups);
    static void remove_inactive_groups(std::list<std::list<Event*> >& groups);
    Event& operator=(const Event& e);
    double calc_event_vel() const;
    double calc_event_accel() const;
    double calc_vevent_tol() const;
    double calc_aevent_tol() const;
    EventClass determine_event_class() const;
    CoulombFrictionType get_friction_type() const {return _ftype; }
    bool is_impacting() const { return determine_event_class() == eNegative; }
    bool is_resting() const { return determine_event_class() == eZero; }
    bool is_separating() const { return determine_event_class() == ePositive; }
    void set_contact_parameters(const ContactParameters& cparams);
    void determine_contact_tangents();
    static void determine_minimal_set(std::list<Event*>& group);
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return GLOBAL; }
    void compute_event_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void compute_cross_event_data(const Event& e, Ravelin::MatrixNd& M) const;

    template <class OutputIterator>
    OutputIterator get_super_bodies(OutputIterator begin) const;

    /// The type of event
    EventType event_type;

    /// The derivative type of the event
    DerivType deriv_type;

    /// The time that the event occurs [0,1]
    double t;

    /// The "real" time that the event occurs [0, infinity]
    double t_true;

    /// The joint at which the limit is reached (for limit events)
    JointPtr limit_joint;

    /// The coefficient of restitution for this limit
    double limit_epsilon;

    /// The DOF at which the limit is reached (for limit events)
    unsigned limit_dof;

    /// Whether the upper/lower limit is reached (for limit events)
    bool limit_upper;

    /// Limit impulse magnitude (for limit events)
    double limit_impulse;

    /// Constraint [normal] impulse magnitude (for constraint events)
    Ravelin::VectorNd constraint_nimpulse;

    /// Constraint [friction] impulse magnitude (for constraint events)
    Ravelin::VectorNd constraint_fimpulse;

    /// The joint (for constraint events)
    JointPtr constraint_joint;

    /// The point contact (for contact events)
    Point3d contact_point;

    /// The vector pointing outward from the contact on the first body, in world coordinates (for contact events)
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

    /// Impulse that has been applied (for contact events)
    /**
     * Impulse applied to the body of the first geometry; 
     * the reverse of this force / impulse is applied to the body of the second 
     * geometry.
     */
    Ravelin::SMomentumd contact_impulse;

    /// The collision geometries involved (for contact events)
    CollisionGeometryPtr contact_geom1, contact_geom2;

    /// The coefficient of Coulomb friction (for contact events)
    double contact_mu_coulomb;

    /// The coefficient of viscous friction (for contact events)
    double contact_mu_viscous;

    /// The coefficient of restitution (for contact events)
    double contact_epsilon;

    /// The number of friction directions >= 4 (for contact events)
    unsigned contact_NK;

    osg::Node* to_visualization_data() const;

    /// Tolerance for the event (users never need to modify this)
    double tol;

    /// Sticking tolerance for acceleration-level events (users never need to modify this)
    double stick_tol;

    void write_vrml(const std::string& filename, double sphere_radius = 0.1, double normal_length = 1.0) const;
    bool operator<(const Event& e) const { return t < e.t; }

  private:
    // structure for comparing pairs of doubles
    struct DblComp
    {
      bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b)
      {
        return (a.first < b.first - NEAR_ZERO && a.second < b.second - NEAR_ZERO);
      }
    };


    // static variables
    boost::shared_ptr<Ravelin::Pose3d> _event_frame;
    static Ravelin::MatrixNd JJ, J, Jx, Jy, J1, J2, dJ1, dJ2, workM1, workM2;
    static Ravelin::VectorNd v, workv, workv2;

    // the type of friction contact (acceleration-level only)
    CoulombFrictionType _ftype;
    void compute_vevent_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void compute_cross_vevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_aevent_data(Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    void compute_cross_aevent_data(const Event& e, Ravelin::MatrixNd& M) const;	
    void compute_cross_contact_contact_vevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_vevent_data(const Event& e, Ravelin::MatrixNd& M, DynamicBodyPtr su) const;
    void compute_cross_contact_contact_vevent_data(const Event& e, Ravelin::MatrixNd& M, DynamicBodyPtr su, const Ravelin::MatrixNd& J) const;
    void compute_cross_contact_limit_vevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_contact_vevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_limit_vevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_aevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_contact_contact_aevent_data(const Event& e, Ravelin::MatrixNd& M, DynamicBodyPtr su) const;
    void compute_cross_contact_contact_aevent_data(const Event& e, Ravelin::MatrixNd& M, DynamicBodyPtr su, const Ravelin::MatrixNd& J) const;
    void compute_cross_contact_limit_aevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_contact_aevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_cross_limit_limit_aevent_data(const Event& e, Ravelin::MatrixNd& M) const;
    void compute_dotv_data(Ravelin::VectorNd& q) const;
    static bool is_linked(const Event& e1, const Event& e2);
    unsigned get_super_bodies(DynamicBodyPtr& sb1, DynamicBodyPtr& sb2) const;
    static void determine_convex_set(std::list<Event*>& group);
    static void process_convex_set_group(std::list<Event*>& group);
    static bool is_contact_manifold_2D(const std::list<Event*>& group);

    template <class BidirectionalIterator>
    static void insertion_sort(BidirectionalIterator begin, BidirectionalIterator end);
}; // end class

std::ostream& operator<<(std::ostream& out, const Event& e);

#include "Event.inl"

} // end namespace Moby

#endif

