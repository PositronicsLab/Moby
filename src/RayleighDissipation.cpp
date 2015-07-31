/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/RCArticulatedBody.h>
#include <Moby/DynamicBody.h>
#include <Moby/RayleighDissipation.h>
#include <Moby/XMLTree.h>

using std::endl;
using std::set;
using std::list;
using std::vector;
using std::map;
using std::make_pair;
using std::multimap;
using std::pair;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

RayleighDissipation::RayleighDissipation()
{
}

void RayleighDissipation::apply(const std::vector<DynamicBodyPtr>& bodies)
{
  const double LAMBDA = 0.99;
  VectorNd v, Mv;
  MatrixNd M;

  // loop through all bodies
  for (unsigned i=0; i< bodies.size(); i++)
  {
    // get lambda 
    double lambda = LAMBDA;
    std::map<DynamicBodyPtr, double>::const_iterator body_iter;
    if ((body_iter = _coeffs.find(bodies[i])) != _coeffs.end())
      lambda = body_iter->second;

    // get generalized velocity
    bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, v);

    // compute the amount of force necessary to bring the body to rest
    double vnrm_sq = v.norm_sq();

    if(vnrm_sq < NEAR_ZERO)
      continue;

    // get the generalized inertia
    bodies[i]->get_generalized_inertia(M);

    // get the momentum
    M.mult(v, Mv);

    double f = Mv.norm() / std::sqrt(vnrm_sq);

    // see whether the force is sufficient
    double diff = 0.5 * vnrm_sq * lambda - f;
    if (diff >= 0.0)
    {
      v.set_zero();
      bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, v);
    }
    else
    {
      f += diff;
      v *= -f;
      bodies[i]->apply_generalized_impulse(v);
    }
  }
}

/// Implements Base::load_from_xml()
void RayleighDissipation::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // call child method
  Base::load_from_xml(node, id_map);

  // verify that the name of this node is correct
  assert(strcasecmp(node->name.c_str(), "RayleighDissipation") == 0);

  // look for coefficients
  std::list<shared_ptr<const XMLTree> > body_nodes = node->find_child_nodes("Body");
  for (std::list<shared_ptr<const XMLTree> >::const_iterator i = body_nodes.begin(); i != body_nodes.end(); i++)
  {
    XMLAttrib* id_attr = (*i)->get_attrib("id");
    if (!id_attr)
    {
      std::cerr << "RayleighDissipation::load_from_xml() - no 'id' specified in 'Body' node" << std::endl;
      continue;
    }

    // look for the ID 
    std::map<std::string, BasePtr>::const_iterator id_map_iter = id_map.find(id_attr->get_string_value());
    if (id_map_iter == id_map.end())
    {
      std::cerr << "RayleighDissipation::load_from_xml() - id '" << id_attr->get_string_value() << "' not found in id map" << std::endl;
      continue;
    }

    // attempt to cast the ID as a body 
    DynamicBodyPtr db = dynamic_pointer_cast<DynamicBody>(id_map_iter->second);
    if (!db)
    {
      std::cerr << "RayleighDissipation::load_from_xml() - id '" << id_attr->get_string_value() << "' not castable to a dynamic body" << std::endl;
      continue;
    }

    // get lambda
    XMLAttrib* coeff_attr = (*i)->get_attrib("lambda");
    if (!coeff_attr)
    {
      std::cerr << "RayleighDissipation::load_from_xml() - no 'lambda' specified in 'Body' node" << std::endl;
      continue;
    }

    // all checks passed; add it to the map 
    _coeffs[db] = coeff_attr->get_real_value();
  }
}

/// Implements Base::save_to_xml()
void RayleighDissipation::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const 
{
  // call parent method
  Base::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "RayleighDissipation";

  // loop through all bodies
  for (std::map<DynamicBodyPtr, double>::const_iterator i = _coeffs.begin(); i != _coeffs.end(); i++)
  {
    XMLTreePtr child_node(new XMLTree("Body"));
    child_node->attribs.insert(XMLAttrib(i->first->id, i->second));
    node->add_child(child_node);
  }
}
 
