/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <strings.h>
#include <Moby/VariableStepIntegrator.h>

using boost::shared_ptr;
using namespace Moby;

/// Implements Base::load_from_xml()
void VariableStepIntegrator::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) 
{ 
  Integrator::load_from_xml(node, id_map); 
 
  // get the absolute error tolerance
  XMLAttrib* aerr_attrib = node->get_attrib("abs-err-tol");
  if (aerr_attrib)
    aerr_tolerance = aerr_attrib->get_real_value();

  // get the relative error tolerance
  XMLAttrib* rerr_attrib = node->get_attrib("rel-err-tol");
  if (rerr_attrib)
    rerr_tolerance = rerr_attrib->get_real_value();

  // get the minimum step size
  XMLAttrib* mss_attrib = node->get_attrib("min-step-size");
  if (mss_attrib)
    min_step_size = mss_attrib->get_real_value();
}

/// Implements Base::save_to_xml()
void VariableStepIntegrator::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{ 
  Integrator::save_to_xml(node, shared_objects); 
  node->name = "VariableStepIntegrator";

  // save the method data
  node->attribs.insert(XMLAttrib("rel-err-tol", rerr_tolerance));
  node->attribs.insert(XMLAttrib("abs-err-tol", aerr_tolerance));
  node->attribs.insert(XMLAttrib("min-step-size", min_step_size));
}

