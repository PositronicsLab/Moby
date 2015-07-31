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

    // the minimum step that the simulator should take (default = 1e-8)
    double min_step_size;

    /// Determines whether two geometries are not checked
    std::set<sorted_pair<CollisionGeometryPtr> > unchecked_pairs;

    /// Vectors set and passed to collision detection
    std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> > _x0, _x1;

    /// Gets the shared pointer for this
    boost::shared_ptr<TimeSteppingSimulator> get_this() { return boost::dynamic_pointer_cast<TimeSteppingSimulator>(shared_from_this()); }
    
  protected:
    bool constraints_met(const std::vector<PairwiseDistInfo>& current_pairwise_distances);
    std::set<sorted_pair<CollisionGeometryPtr> > get_current_contact_geoms() const;
    double do_mini_step(double dt);
    void step_si_Euler(double dt);
    double calc_next_CA_Euler_step(double contact_dist_thresh) const;
    void calc_impacting_unilateral_constraint_forces2(double dt);

  private:
    double record_error(RigidBodyPtr rb);
    double record_error(DynamicBodyPtr db);
    void store_state(RigidBodyPtr rb);
    void store_state(DynamicBodyPtr db);
    void update_state(RigidBodyPtr rb, double h);
    void update_state(DynamicBodyPtr db, double h);
    void integrate_bodies(const std::vector<DynamicBodyPtr>& bodies, double h);
    bool calc_integration_error(const std::vector<DynamicBodyPtr>& bodies, double h, const std::vector<Ravelin::VectorNd>& qe_large, const std::vector<Ravelin::VectorNd>& v_large, const std::vector<Ravelin::VectorNd>& qe_small, const std::vector<Ravelin::VectorNd>& v_small, double& min_k, std::string logging_flag = std::string());
}; // end class

} // end namespace

#endif


