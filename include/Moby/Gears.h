/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GEARS_H
#define _GEARS_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/Geard.h>
#include <Moby/Base.h>
#include <Moby/DynamicBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Visualizable.h>

namespace Moby {

class VisualizationData;

/// Defines a gear ("joint") constraint 
/**
 * \todo implement a rest position for q?
 */
class Gears : public Gear 
{
  public:

    Gear();
    Gear(boost::weak_ptr<RigidBody> inboard, boost::weak_ptr<RigidBody> outboard);
    virtual ~Gear() {}
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_link(RigidBodyPtr link, bool update_pose);
    virtual void set_outboard_link(RigidBodyPtr link, bool update_pose);

  private:
    double _ratio; // the input:output ratio for the gears
}; // end class
} // end namespace

#endif

