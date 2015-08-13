/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _DISSIPATION_H
#define _DISSIPATION_H

#include <Moby/RecurrentForce.h>

namespace Moby {
class Dissipation : public Base 
{
  public:
    Dissipation();
    virtual ~Dissipation() {}
    void apply(const std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> >& bodies);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The mapping from bodies to decay coefficients 
    std::map<boost::shared_ptr<Ravelin::DynamicBodyd>, double> _coeffs; 
}; // end class
} // end namespace

#endif

