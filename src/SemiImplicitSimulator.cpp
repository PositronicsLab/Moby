/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <cmath>
#include <limits>
#include <list>
#include <Moby/Constants.h>
#include <Moby/Contact.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RecurrentForce.h>
#include <Moby/CollisionDetection.h>
#include <Moby/Integrator.h>
#include <Moby/XMLTree.h>
#include <Moby/SemiImplicitSimulator.h>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::list;
using std::vector;
using std::map;
using std::set;
using std::multimap;
using std::string;
using std::endl;

#define DO_TIMING

/// Constructs a SemiImplicitSimulator with no parameters
SemiImplicitSimulator::AnitescuPotraSimulator() 
{
  toi_tolerance = NEAR_ZERO;
}

/// This is the naive method
/*
Real SemiImplicitSimulator::step(Real step_size)
{
  Real TOC;
  typedef std::pair<DynamicBodyPtr, shared_ptr<void> > StatePair;
  multimap<Real, Contact> contact_map;

  // get reference to the list of contacts
  list<Contact>& contacts = this->contacts; 

  // setup the amount remaining to step
  Real dt = step_size;

  // determine external forces
  determine_external_forces();

  // if there are no collision detectors, don't do any contact handling...
  if (collision_detectors.empty())
    return ContactSimulator::step(step_size);

  // clear one-step visualization data
  #ifdef USE_INVENTOR
  _transient_vdata->removeAllChildren();
  #endif

  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << std::endl;

  // look for contact
  ContactSimulator::is_contact(dt, contact_map);
  TOC = find_TOC(dt, contact_map, this->contacts);

  // remove contacts that do not reference bodies
  _bodies.sort();
  for (list<Contact>::iterator i = contacts.begin(); i != contacts.end(); )
  {
    RigidBodyPtr rb1 = i->geom1->get_rigid_body();
    RigidBodyPtr rb2 = i->geom2->get_rigid_body();

    if (std::find(_bodies.begin(), _bodies.end(), rb1) == _bodies.end() ||
        std::find(_bodies.begin(), _bodies.end(), rb2) == _bodies.end())
      i = contacts.erase(i);
    else
      i++;
  }

  // if there is no contact, just call the underlying simulator method
  if (TOC > dt || contacts.empty())
  {
    FILE_LOG(LOG_CONTACT) << " no contacts -- stepping simulator regularly" << std::endl;
    return Simulator::step(step_size);
  }

  // integrate all bodies forward by the time-of-contact
  integrate(TOC);

  FILE_LOG(LOG_SIMULATOR) << "  -- integrating bodies forward by " << TOC << std::endl;

  // prep contacts and treat impacts
  prep_contacts(contacts);
  treat_impacts(contacts);

  FILE_LOG(LOG_SIMULATOR) << "  -- treated impacts" << std::endl;

  // update dt to indicate how much we stepped
  dt = TOC;

  // clear set of contacts 
  contacts.clear();

  // reset force and torque accumulators for all bodies
  BOOST_FOREACH(DynamicBodyPtr body, _bodies)
    body->reset_accumulators();

  // indicate that external forces have not been computed for the next step
  _external_forces_determined = false;

  // update the current time
  current_time += dt;

  FILE_LOG(LOG_SIMULATOR) << "+integration complete" << std::endl;
  
  return dt; 
}
*/

/// Implements virtual method Simulator::step() for stepping the simulation
/**
 * The integration process proceeds as follows:
 * 1. Velocities are integrated forward by step_size
 * 2. Time of impact (TOI) is determined by collision detection (using new
 *    velocities)
 * 3. Positions are stepped to TOI using the new velocities
 * 4. Impacts are treated
 * 5. (repeat #2-#4 until step_size has been taken)
 */
Real SemiImplicitSimulator::step(Real step_size)
{
  Real TOC;
  typedef std::pair<DynamicBodyPtr, shared_ptr<void> > StatePair;
  multimap<Real, Contact> contact_map;

  // get reference to the list of contacts
  list<Contact>& contacts = this->contacts; 

  // determine external forces
  determine_external_forces();

  // if there are no collision detectors, don't do any contact handling...
  if (collision_detectors.empty())
    return ContactSimulator::step(step_size);

  // clear one-step visualization data
  #ifdef USE_INVENTOR
  _transient_vdata->removeAllChildren();
  #endif

  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << std::endl;

  // integrate the velocities forward by the step size
  integrate_velocities(step_size);

  // look for contact
  ContactSimulator::is_contact(step_size, contact_map);
  Real TOI = find_TOI(step_size, contact_map, this->contacts);

  // if there is no impact, just call the underlying simulator method
  if (TOI > step_size || contacts.empty())
  {
    FILE_LOG(LOG_CONTACT) << " no contacts -- stepping simulator regularly" << std::endl;
    return Simulator::step(step_size);
  }

  // setup amount remaining to step
  Real dt = step_size;

  // loop until no time remaining
  while (dt > 0.0)
  { 
    // step positions to TOI
    integrate_positions(TOI);
    dt -= TOI;

    // treat impacts
    prep_contacts(contacts)
    treat_impacts(contacts);

    // determine treated rigid bodies
    std::list<RigidBodyPtr> treated_rbs;
    determine_treated_bodies(contacts, treated_rbs);

    // subtract TOI from all keys in the contact map
    subtract_TOI(contact_map, TOI);

    // update TOIs for treated bodies
    ContactSimulator::update_contact_map(dt, treated_rbs, contact_map);
    TOI = find_TOI(dt, contact_map, this->contacts);

    FILE_LOG(LOG_SIMULATOR) << " -- new contact map: " << std::endl;
    for (multimap<Real, Contact>::const_iterator i = contact_map.begin(); i != contact_map.end(); i++)
      FILE_LOG(LOG_SIMULATOR) << "    TOI: " << i->first << std::endl << i->second;

    // if there is no impact, just call the underlying simulator method
    if (TOI > dt || contacts.empty())
    {
      FILE_LOG(LOG_CONTACT) << " no contacts -- stepping simulator regularly" << std::endl;
      return Simulator::step(dt);
    }
  }

  // clear set of contacts 
  contacts.clear();

  // reset force and torque accumulators for all bodies
  BOOST_FOREACH(DynamicBodyPtr body, _bodies)
    body->reset_accumulators();

  // indicate that external forces have not been computed for the next step
  _external_forces_determined = false;

  // update the current time
  current_time += dt;

  FILE_LOG(LOG_SIMULATOR) << "+integration complete" << std::endl;
  
  return dt; 
}

/// Subtracts a TOI from all keys in a contact map
void SemiImplicitSimulator::subtract_TOI(multimap<Real, Contact>& contact_map, Real TOI)
{
  FILE_LOG(LOG_SIMULATOR) << "subtract_TOI() entered" << std::endl;

  multimap<Real, Contact> updated_contact_map;
  for (multimap<Real, Contact>::const_iterator i = contact_map.begin(); i != contact_map.end(); i++)
  {
    FILE_LOG(LOG_SIMULATOR) << " -- changing TOI for contact from " << i->first << " to " << (i->first - TOI) << std::endl;
    updated_contact_map.insert(std::make_pair(i->first - TOI, i->second));
  }
  contact_map = updated_contact_map;

  FILE_LOG(LOG_SIMULATOR) << "subtract_TOI() exited" << std::endl;
}

/// Determines the set of bodies from a set of contacts
void SemiImplicitSimulator::get_treated_bodies(const list<Contact>& contacts, list<RigidBodyPtr>& bodies)
{
  // loop through all contacts
  BOOST_FOREACH(const Contact& c, contacts)
  {
    // get the rigid bodies of the contacts
    RigidBodyPtr rb1 = c.geom1->get_rigid_body();
    RigidBodyPtr rb2 = c.geom2->get_rigid_body();

    // get the articulated bodies, if any
    ArticulatedBodyPtr ab1 = rb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = rb2->get_articulated_body();

    // add the bodies to the list
    if (ab1)
      bodies.insert(bodies.end(), ab1->get_links().begin(), ab1->get_links().end());
    else if (rb1->is_enabled())
      bodies.push_back(rb1);

    if (ab2)
      bodies.insert(bodies.end(), ab2->get_links().begin(), ab2->get_links().end());
    else if (rb2->is_enabled())
      bodies.push_back(rb2);
  }

  // sort the list and make it unique
  bodies.sort();
  bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());
}

/// Finds the time-of-impact of first impacting contact and sets up the contact map
Real SemiImplicitSimulator::find_TOI(Real dt, multimap<Real, Contact>& contact_map, std::list<Contact>& contacts) const
{
  const Real INF = std::numeric_limits<Real>::max();

  // clear the list of contacts
  contacts.clear();

  // get the iterator start
  multimap<Real, Contact>::iterator citer = contact_map.begin();

  // loop while the iterator does not point to the end 
  while (citer != contact_map.end())
  {
    // set tmin
    Real tmin = citer->first;

    // check for early exit
    if (tmin > dt)
      return INF;

    // check for impacting contact
    bool impacting = is_impacting(citer->second);

    // find the first contacts
    for (citer++; citer != contact_map.end(); citer++)
    {
      // see whether we are done
      if (citer->first > tmin + toi_tolerance)
        break;

      // see whether this contact is impacting (if we don't yet have an
      // impacting contact)
      if (!impacting)
        impacting = is_impacting(citer->second);
    }

    // see whether we are done
    if (impacting)
    {
      // make contact set
      contacts.clear();
      for (multimap<Real, Contact>::const_iterator j = contact_map.begin(); j != citer; j++)
      {
        FILE_LOG(LOG_SIMULATOR) << "Detected impacting contact at TOI " << j->first << ": " << endl << j->second;
        contacts.push_back(j->second);
      }

      // modify contact map
      contact_map.erase(contact_map.begin(), citer);

      return tmin;
    }
  }

  // contact map is empty, no contacts
  return INF;
}

/// Determines whether a contact is impacting 
bool SemiImplicitSimulator::is_impacting(const Contact& c) const 
{
  // get the two bodies
  RigidBodyPtr b1 = c.geom1->get_rigid_body();
  RigidBodyPtr b2 = c.geom2->get_rigid_body();

  // get the relative velocity at the contact point
  Vector3 lpv = b1->get_lvel() - b2->get_lvel();
  Vector3 av1 = Vector3::cross(b1->get_avel(), c.point - b1->get_position());
  Vector3 av2 = Vector3::cross(b2->get_avel(), c.point - b2->get_position());
  Real rvel = c.normal.dot(lpv + av1 - av2);

  FILE_LOG(LOG_SIMULATOR) << "is_impacting(): bodies " << b1->id << " and " << b2->id << " rvel: " << rvel << endl;

  // determine whether it is impacting 
  return rvel <= -collision_tolerance; 
}

/// Loads the dynamic state of all bodies
void SemiImplicitSimulator::load_states(const DynamicStates& state_t)
{
  typedef std::pair<DynamicBodyPtr, shared_ptr<void> > StatePair;

  BOOST_FOREACH(const StatePair& sp, state_t)
    sp.first->load_dynamic_state(sp.second);
}

/// Saves the dynamic state of all bodies
void SemiImplicitSimulator::save_states(DynamicStates& state_t)
{
  typedef std::pair<DynamicBodyPtr, shared_ptr<void> > StatePair;

  BOOST_FOREACH(StatePair& sp, state_t)
    sp.second = sp.first->save_dynamic_state();
}

/// Determines the set of contacts
/**
 * Updates contacts and modifies the contact map appropriately.
 */
void SemiImplicitSimulator::determine_contacts(multimap<Real, Contact>& contact_map, Real TOI, list<Contact>& contacts)
{
  // find all contacts within the desired tolerance of the TOI
  multimap<Real, Contact>::iterator end = contact_map.begin();
  for (end++; end != contact_map.end(); end++)
    if (end->first > TOI + toi_tolerance)
      break;

  FILE_LOG(LOG_SIMULATOR) << "  -- TOI: " << TOI << std::endl;
  FILE_LOG(LOG_SIMULATOR) << "  -- # of contact points: " << std::distance(contact_map.begin(), end) << std::endl;

  // update the set of contacts
  contacts.clear();
  for (multimap<Real, Contact>::const_iterator i = contact_map.begin(); i != end; i++)
    contacts.push_back(i->second);

  // remove those contacts that we will treat from the multimap
  contact_map.erase(contact_map.begin(), end);
}

/// Performs all contact handling
void SemiImplicitSimulator::treat_impacts(list<Contact>& contacts)
{
  FILE_LOG(LOG_SIMULATOR) << "SemiImplicitSimulator::treat_impacts() entered" << std::endl;

  // ** first, we want to remove some contact points (randomly) as necessary
  // determine maximum points for each pair of rigid bodies and contact
  // points for each pair of rigid bodies
  map<sorted_pair<RigidBodyPtr>, unsigned> max_points;
  map<sorted_pair<RigidBodyPtr>, vector<Contact> > rb_contacts;
  BOOST_FOREACH(const Contact& c, contacts)
  {
    // get the rigid bodies
    RigidBodyPtr rb1 = c.geom1->get_rigid_body();
    RigidBodyPtr rb2 = c.geom2->get_rigid_body();
    sorted_pair<RigidBodyPtr> sp(rb1, rb2);

    // determine the maximum # of contact points
    unsigned maxp = c.get_contact_data()->max_points;
    map<sorted_pair<RigidBodyPtr>, unsigned>::iterator max_points_iter;
    if ((max_points_iter = max_points.find(sp)) == max_points.end())
      max_points[sp] = maxp;
    else
      max_points_iter->second = std::min(max_points_iter->second, maxp);

    // add the contact to the list of contacts for the pair
    rb_contacts[sp].push_back(c);
  }  

  // clear the set of contacts -- it will be rebuilt with new points
  contacts.clear();

  // randomly remove some contact points as necessary
  for (map<sorted_pair<RigidBodyPtr>, vector<Contact> >::iterator i = rb_contacts.begin(); i != rb_contacts.end(); i++)
  {
    // determine how many contact points for this pair
    unsigned nc = i->second.size();

    // determine the maximum number of points for this pair
    unsigned maxp = max_points[i->first];

    // remove some if necessary
    if (nc > maxp)
    {
      std::random_shuffle(i->second.begin(), i->second.end());
      while (nc > maxp)
      {
        i->second.pop_back();
        nc--;
      }
    }

    // copy the points to the new list of contacts
    std::copy(i->second.begin(), i->second.end(), std::back_inserter(contacts));
  }

  // call the pre-impulse application callback
  if (contact_impulse_callback_fn)
    (*contact_impulse_callback_fn)(contact_impulse_callback_data);

  // **** process the impacts
  // add the visualization for contacts to the transient 
  // visualization data and output impacting pairs
  FILE_LOG(LOG_SIMULATOR) << "  -- processing contacts: " << std::endl;
  BOOST_FOREACH(const Contact& c, contacts)
    FILE_LOG(LOG_SIMULATOR) << c << std::endl;

  // get the map of all collision methods (and contacts to be processed
  // by them)
  FILE_LOG(LOG_SIMULATOR) << " -- associating contacts with collision methods..." << std::endl;

  map<shared_ptr<CollisionMethod>, list<Contact> > cms;
  BOOST_FOREACH(const Contact& c, contacts)
  {
    // verify that contact data exists
    assert(c.get_contact_data());

    // call the collision method for this pair
    cms[c.get_contact_data()->collision_method].push_back(c);

    FILE_LOG(LOG_SIMULATOR) << "  -- added contact " << c << " to collision method " << c.get_contact_data()->collision_method << "(" << c.get_contact_data()->collision_method->id << ")" << std::endl;
  }

  #ifdef DO_TIMING
  std::ofstream out("contact-timing", std::ofstream::app);
  clock_t start = clock();
  #endif

  // process all contacts, one contact method at a time, with local
  // methods going first
  for (map<shared_ptr<CollisionMethod>, list<Contact> >::const_iterator i = cms.begin(); i != cms.end(); i++)
    if (!i->first->is_global_method())
      i->first->process_impacts(i->second);

  // now do global methods
  for (map<shared_ptr<CollisionMethod>, list<Contact> >::const_iterator i = cms.begin(); i != cms.end(); i++)
    if (i->first->is_global_method())
      i->first->process_impacts(i->second);

  #ifdef DO_TIMING
  clock_t stop = clock();
  out << current_time << " " << (((double) (stop-start))/CLOCKS_PER_SEC) << std::endl;
  out.close();
  #endif

  FILE_LOG(LOG_SIMULATOR) << "SemiImplicitSimulator::treat_impacts() exited" << std::endl;
}

/// Processes set of contacts from current simulation state
void SemiImplicitSimulator::prep_contacts(list<Contact>& contacts)
{
  // call the contact callback (if any)
  if (contact_callback_fn)
    (*contact_callback_fn)(contacts, contact_callback_data);

  // store contacts with points and contact data
  for (list<Contact>::iterator i = contacts.begin(); i != contacts.end(); )
  {
    // skip contacts with no points
    if (i->contact_geometry == Contact::eNone)
    {
      i = contacts.erase(i);
      continue;
    }

    // preprocess the contact (get contact data)
    preprocess_contact(*i);

    // if the contact has no contact data or no collision method, ignore it
    if (!i->get_contact_data())
    {
      i = contacts.erase(i);
      continue;
    }

    // all tests here..  keep this contact
    i++;
  }
}

/// Implements Base::load_from_xml()
void SemiImplicitSimulator::load_from_xml(XMLTreeConstPtr node, map<string, BasePtr>& id_map)
{
  list<XMLTreeConstPtr> child_nodes;
  map<string, BasePtr>::const_iterator id_iter;

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "SemiImplicitSimulator") == 0);

  // first, load all data specified to the Simulator object
  ContactSimulator::load_from_xml(node, id_map);

  // load load all data specified to the SemiImplicitRestitutionModel object
  SemiImplicitRestitutionModel::load_from_xml(node, id_map);

  // load TOI tolerance, if specified
  const XMLAttrib* toi_tol_attr = node->get_attrib("TOI-tolerance");
  if (toi_tol_attr)
    this->toi_tolerance = toi_tol_attr->get_real_value();
}

/// Implements Base::save_to_xml()
void SemiImplicitSimulator::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // call Simulator's save method first
  ContactSimulator::save_to_xml(node, shared_objects);

  // call restitution model's save method
  SemiImplicitRestitutionModel::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "SemiImplicitSimulator";

  // save TOI tolerance
  node->attribs.insert(XMLAttrib("TOI-tolerance", toi_tolerance));
}

/// Outputs this class data to the stream
/**
 * This method outputs all of the low-level details to the stream; if
 * serialization is desired, use save_to_xml() instead.
 * \sa save_to_xml()
 */
void SemiImplicitSimulator::output_object_state(std::ostream& out) const
{
  // call parent method
  ContactSimulator::output_object_state(out);

  // indicate the object type
  out << "SemiImplicitSimulator object" << std::endl; 

  // output resting contacts
  out << "  contacts: " << std::endl;
  for (list<Contact>::const_iterator i = this->contacts.begin(); i != this->contacts.end(); i++)
    out << "    contact: " << &*i << std::endl;
}


