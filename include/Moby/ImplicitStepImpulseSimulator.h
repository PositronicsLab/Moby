/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _IMPLICIT_STEP_IMPULSE_SIMULATOR_H
#define _IMPLICIT_STEP_IMPULSE_SIMULATOR_H

#include <set>
#include <map>
#include <Moby/sorted_pair>
#include <Moby/RecurrentForce.h>
#include <Moby/ContactSimulator.h>
#include <Moby/CollisionMethod.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactData.h>
#include <Moby/ArticulatedBody.h>

namespace Moby {

class Contact;

/// A contact handling simulator based on the method of Anitescu and Potra 
class ImplicitStepImpulseSimulator : public ContactSimulator
{
  public:
    ImplicitStepImpulseSimulator();
    virtual ~ImplicitStepImpulseSimulator() {}
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void output_object_state(std::ostream& out) const;

    /// Steps the simulation forward
    /**
     * The simulator will take less than the step_size requested if
     * collision impulses are applied.  If the full step size is requested, call step_full().
     * \return the step size taken (<= step_size)
     */
    virtual Real step(Real step_size);

    /// The set of contacts determined on the last call to step()
    std::list<Contact> contacts;

    /// The tolerance over the first time of impact with which additional points of contact are treated as impacting
    /**
     * If this value is set too low, then contact points might be missed; this
     * is particularly a problem for objects in resting contact.  If this value
     * is set too high, then parts of objects not contacting will be 
     * considered.  Default value is NEAR_ZERO. 
     */
    Real toi_tolerance;

  protected:
    virtual void prep_contacts(std::list<Contact>& contacts);
    bool active_contact(Real dt, const std::list<Contact>& contacts);
    static void* find_contact(void* arg);

  private:
    static void load_states(const DynamicStates& state_t);
    static void save_states(DynamicStates& state_t);
    Real find_TOC(Real dt, std::multimap<Real, Contact>& contact_map, std::list<Contact>& contacts) const;
    bool is_active(const Contact& c) const;
    Real min_contact_vel(const std::list<Contact>& contacts) const;
    static void determine_contacting_bodies(const std::list<Contact>& contacts, std::set<DynamicBodyPtr>& contacting_bodies);

    void integrate_contacting(const std::list<Contact>& contacts, const std::set<DynamicBodyPtr>& contacting_bodies, Real dt);
    void integrate_non_contacting(const std::list<Contact>& contacts, const std::set<DynamicBodyPtr>& contacting_bodies, Real dt);
    void determine_contacts(std::multimap<Real, Contact>& contact_map, Real TOI, std::list<Contact>& contacts);
    static void get_treated_bodies(const std::list<Contact>& contacts, std::list<DynamicBodyPtr>& treated_dbs, std::list<RigidBodyPtr>& treated_rbs);
    void treat_impacts(std::list<Contact>& contacts);
}; // end class

} // end namespace Moby

#endif
