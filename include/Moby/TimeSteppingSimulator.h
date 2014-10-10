/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _TIME_STEPPING_SIMULATOR_H
#define _TIME_STEPPING_SIMULATOR_H

#include <Moby/EventDrivenSimulator.h>

namespace Moby {

class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// A time-stepping simulator
class TimeSteppingSimulator : public EventDrivenSimulator
{
  friend class CollisionDetection;

  public:
    TimeSteppingSimulator();
    virtual ~TimeSteppingSimulator() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual double step(double dt);

    /// the contact distance threshold
    double contact_dist_thresh;

    /// Gets the shared pointer for this
    boost::shared_ptr<TimeSteppingSimulator> get_this() { return boost::dynamic_pointer_cast<TimeSteppingSimulator>(shared_from_this()); }
}; // end class

} // end namespace

#endif


