/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <boost/tuple/tuple.hpp>
#include <Moby/XMLTree.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/DynamicBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CollisionDetection.h>
#include <Moby/ContactParameters.h>
#include <Moby/VariableStepIntegrator.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/InvalidStateException.h>
#include <Moby/InvalidVelocityException.h>
#include <Moby/TimeSteppingSimulator.h>

#ifdef USE_OSG
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/PositionAttitudeTransform>
#include <osg/Quat>
#endif // USE_OSG

using std::endl;
using std::set;
using std::list;
using std::vector;
using std::map;
using std::make_pair;
using std::multimap;
using std::pair;
using boost::tuple;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Default constructor
TimeSteppingSimulator::TimeSteppingSimulator()
{
  constraint_callback_fn = NULL;
  constraint_post_callback_fn = NULL;
  post_mini_step_callback_fn = NULL;
  get_contact_parameters_callback_fn = NULL;
  render_contact_points = false;
  contact_dist_thresh = 1e-4;
}

/// Steps the simulator forward by the given step size
double TimeSteppingSimulator::step(double step_size)
{
  const double INF = std::numeric_limits<double>::max();
  current_time += step_size;

  // clear timings
  dynamics_time = (double) 0.0;
  constraint_time = (double) 0.0;
  coldet_time = (double) 0.0;

  // setup timer
  clock_t start = clock();

  // determine the set of collision geometries
  determine_geometries();

  // clear one-step visualization data
  #ifdef USE_OSG
  _transient_vdata->removeChildren(0, _transient_vdata->getNumChildren());
  #endif

  FILE_LOG(LOG_SIMULATOR) << "+stepping simulation from time: " << this->current_time << " by " << step_size << std::endl;
  if (LOGGING(LOG_SIMULATOR))
  {
    VectorNd q;
    BOOST_FOREACH(DynamicBodyPtr db, _bodies)
    {
      db->get_generalized_coordinates(DynamicBody::eEuler, q);
      FILE_LOG(LOG_SIMULATOR) << " body " << db->id << " coordinates (before): " << q << std::endl;
    }
  }

  // do broad phase collision detection
  broad_phase(step_size);

  // compute pairwise distances at the current configuration
  calc_pairwise_distances();

  // find unilateral constraints 
  find_unilateral_constraints(contact_dist_thresh);

  // integrate accelerations forward by dt to get new velocities
  integrate_velocities_Euler(step_size);

  // do the impact constraint handler to compute new velocities
  calc_impacting_unilateral_constraint_forces(step_size);

  // integrate positions forward using new velocities
  integrate_positions_Euler(step_size);

  // call the callback 
  if (post_step_callback_fn)
    post_step_callback_fn(this);

  return step_size;
}

/// Implements Base::load_from_xml()
void TimeSteppingSimulator::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  list<shared_ptr<const XMLTree> > child_nodes;
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify node name b/c this is abstract class
  assert(strcasecmp(node->name.c_str(), "TimeSteppingSimulator") == 0);

  // first, load all data specified to the Simulator object
  Simulator::load_from_xml(node, id_map);

  // get the contact distance threshold
  XMLAttrib* contact_dist_thresh_attr = node->get_attrib("contact-dist-thresh");
  if (contact_dist_thresh_attr)
    contact_dist_thresh = contact_dist_thresh_attr->get_real_value(); 

  // read in any ContactParameters
  child_nodes = node->find_child_nodes("ContactParameters");
  if (!child_nodes.empty())
    contact_params.clear();
  for (list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    boost::shared_ptr<ContactParameters> cd(new ContactParameters);
    cd->load_from_xml(*i, id_map);
    contact_params[cd->objects] = cd;
  }

  // read all disabled pairs
  child_nodes = node->find_child_nodes("DisabledPair");
  for (std::list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    // get the two ID attributes
    XMLAttrib* id1_attrib = (*i)->get_attrib("object1-id");
    XMLAttrib* id2_attrib = (*i)->get_attrib("object2-id");

    // make sure that they were read
    if (!id1_attrib || !id2_attrib)
    {
      std::cerr << "TimeSteppingSimulator::load_from_xml() - did not find ";
      std::cerr << "object1-id and/or object2-id" << std::endl;
      std::cerr << "  in offending node: " << std::endl << *node;
      continue;
    }

    // get the two IDs
    const std::string& ID1 = id1_attrib->get_string_value();
    const std::string& ID2 = id2_attrib->get_string_value();

    // setup pairs of geometries to disable
    std::list<CollisionGeometryPtr> disabled1, disabled2;

    // find the first object
    if ((id_iter = id_map.find(ID1)) == id_map.end())
    {
      std::cerr << "TimeSteppingSimulator::load_from_xml() - could not find ";
      std::cerr << "object with object1-id" << std::endl;
      std::cerr << "  '" << ID1 << "' in offending node: " << std::endl << *node;
      continue;
    }
    BasePtr o1 = id_iter->second;
    CollisionGeometryPtr g1 = dynamic_pointer_cast<CollisionGeometry>(o1);
    if (g1)
      disabled1.push_back(g1);
    else
    {
      RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(o1);
      if (rb1)
        disabled1 = rb1->geometries;
      else
      {
        ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(o1);
        if (ab1)
        {
          BOOST_FOREACH(RigidBodyPtr rb, ab1->get_links())
            disabled1.insert(disabled1.end(), rb->geometries.begin(), rb->geometries.end());
        }
        else
        {
          std::cerr << "TimeSteppingSimulator::load_from_xml() - object with object1-id is not a usable type!" << std::endl;
          continue;
        }
      }
    }

    // find the second object
    if ((id_iter = id_map.find(ID2)) == id_map.end())
    {
      std::cerr << "TimeSteppingSimulator::load_from_xml() - could not find ";
      std::cerr << "object with object2-id" << std::endl;
      std::cerr << "  '" << ID2 << "' in offending node: " << std::endl << *node;
      continue;
    }
    BasePtr o2 = id_iter->second;
    CollisionGeometryPtr g2 = dynamic_pointer_cast<CollisionGeometry>(o2);
    if (g2)
      disabled2.push_back(g2);
    else
    {
      RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(o2);
      if (rb2)
        disabled2 = rb2->geometries;
      else
      {
        ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(o2);
        if (ab2)
        {
          BOOST_FOREACH(RigidBodyPtr rb, ab2->get_links())
            disabled2.insert(disabled2.end(), rb->geometries.begin(), rb->geometries.end());
        }
        else
        {
          std::cerr << "TimeSteppingSimulator::load_from_xml() - object with object2-id is not a usable type!" << std::endl;
          continue;
        }
      }
    }

 
   // add the pairs to the unchecked pairs list
   BOOST_FOREACH(CollisionGeometryPtr cg1, disabled1)
     BOOST_FOREACH(CollisionGeometryPtr cg2, disabled2)
       if (cg1 != cg2 && cg1->get_single_body() != cg2->get_single_body())
         unchecked_pairs.insert(make_sorted_pair(cg1, cg2));
  }
}

/// Implements Base::save_to_xml()
void TimeSteppingSimulator::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call Simulator's save method first
  Simulator::save_to_xml(node, shared_objects);

  // reset the node's name
  node->name = "TimeSteppingSimulator";

  // save the contact distance threshold
  node->attribs.insert(XMLAttrib("contact-dist-thresh", contact_dist_thresh));

  // save all ContactParameters
  for (map<sorted_pair<BasePtr>, shared_ptr<ContactParameters> >::const_iterator i = contact_params.begin(); i != contact_params.end(); i++)
  {
    XMLTreePtr new_node(new XMLTree("ContactParameters"));
    node->add_child(new_node);
    i->second->save_to_xml(new_node, shared_objects);
  }

  // save all disabled pairs
  for (std::set<sorted_pair<CollisionGeometryPtr> >::const_iterator i = unchecked_pairs.begin(); i != unchecked_pairs.end(); i++)
  {
    XMLTreePtr child_node(new XMLTree("DisabledPair"));
    child_node->attribs.insert(XMLAttrib("object1-id", i->first->id));
    child_node->attribs.insert(XMLAttrib("object2-id", i->second->id));
    node->add_child(child_node);
  }
}


