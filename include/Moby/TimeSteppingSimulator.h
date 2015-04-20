/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _TS_SIMULATOR_H
#define _TS_SIMULATOR_H

#include <map>
#include <Moby/sorted_pair>
#include <Moby/ConstraintSimulator.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/PenaltyConstraintHandler.h>
#include <Moby/SustainedUnilateralConstraintHandler.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/CCD.h>
#include <Moby/UnilateralConstraint.h>

namespace Moby {

class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// An event-driven simulator
class TimeSteppingSimulator : public ConstraintSimulator
{
  friend class CollisionDetection;

  public:
    TimeSteppingSimulator();
    virtual ~TimeSteppingSimulator() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual double step(double dt);
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;

    /// Determines whether two geometries are not checked
    std::set<sorted_pair<CollisionGeometryPtr> > unchecked_pairs;

    /// Vectors set and passed to collision detection
    std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> > _x0, _x1;

    /// Gets the shared pointer for this
    boost::shared_ptr<TimeSteppingSimulator> get_this() { return boost::dynamic_pointer_cast<TimeSteppingSimulator>(shared_from_this()); }
    
  protected:
    std::set<sorted_pair<CollisionGeometryPtr> > get_current_contact_geoms() const;
    double integrate_forward(double dt);
    void step_si_Euler(double dt);
    double calc_next_CA_Euler_step(double contact_dist_thresh) const;
    void step_si_Euler0(double dt);
    void integrate_velocities_Euler0(double dt);
    void integrate_positions_Euler0(double dt);
    double calc_next_CA_Euler_step0(double contact_dist_thresh) const;
}; // end class

} // end namespace

#endif


