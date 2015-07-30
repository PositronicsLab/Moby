/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _EXP_DISSIPATION_H
#define _EXP_DISSIPATION_H

#include <Moby/Dissipation.h>

namespace Moby {

class ExponentialDissipation : public Dissipation 
{
  public:
    ExponentialDissipation();
    virtual ~ExponentialDissipation() {}
    void apply(const std::vector<DynamicBodyPtr>& bodies);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// The mapping from bodies to decay coefficients 
    std::map<DynamicBodyPtr, double> _coeffs; 
}; // end class
} // end namespace

#endif

