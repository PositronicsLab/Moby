/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CONSTRAINT_SIMULATOR_H
#define _CONSTRAINT_SIMULATOR_H

#include <map>
#include <Moby/sorted_pair>
#include <Moby/Simulator.h>
#include <Moby/ImpactConstraintHandler.h>
#include <Moby/PenaltyConstraintHandler.h>
#include <Moby/SustainedUnilateralConstraintHandler.h>
#include <Moby/PairwiseDistInfo.h>
#include <Moby/CCD.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/ConstraintStabilization.h>

namespace Moby {

class Dissipation;
class ContactParameters;
class CollisionDetection;
class CollisionGeometry;

/// An virtual class for simulation with constraints 
class ConstraintSimulator : public Simulator
{
  friend class CollisionDetection;
  friend class ConstraintStabilization;

  public:
    ConstraintSimulator();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    boost::shared_ptr<ContactParameters> get_contact_parameters(CollisionGeometryPtr geom1, CollisionGeometryPtr geom2) const;

    /// Determines whether two geometries are not checked
    std::set<sorted_pair<CollisionGeometryPtr> > unchecked_pairs;

    /// Vectors set and passed to collision detection
    std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> > _x0, _x1;

    /// Callback function for getting contact parameters
    boost::shared_ptr<ContactParameters> (*get_contact_parameters_callback_fn)(CollisionGeometryPtr g1, CollisionGeometryPtr g2);

    /// Callback function after a mini-step is completed
    void (*post_mini_step_callback_fn)(ConstraintSimulator* s);

    /// The callback function (called when constraints have been determined)
    /**
     * The callback function can remove constraints from the list, which will disable
     * their processing (however, doing so may prevent the simulation from
     * making progress, as the simulator attempts to disallow violations.
     */
    void (*constraint_callback_fn)(std::vector<UnilateralConstraint>&, boost::shared_ptr<void>);

    /// The callback function (called after forces/impulses are applied)
    void (*constraint_post_callback_fn)(const std::vector<UnilateralConstraint>&, boost::shared_ptr<void>);

    /// Data passed to unilateral constraint callback
    boost::shared_ptr<void> constraint_callback_data;
    
    /// Data passed to post-constraint callback
    boost::shared_ptr<void> constraint_post_callback_data;
 
    /// Gets the (sorted) compliant constraint data
    std::vector<UnilateralConstraint>& get_compliant_constraints() { return _compliant_constraints; }

    /// Gets the (sorted) rigid constraint data
    std::vector<UnilateralConstraint>& get_rigid_constraints() { return _rigid_constraints; }

    /// Mapping from objects to contact parameters
    std::map<sorted_pair<BasePtr>, boost::shared_ptr<ContactParameters> > contact_params;

    /// If set to 'true' simulator will process contact points for rendering
    bool render_contact_points;

    /// Gets the collision detection mechanism
    boost::shared_ptr<CollisionDetection> get_collision_detection() const { return _coldet; }

    /// The distance threshold for a contact to be handled
    /**
     * Bodies are only considered to be in contact 
     * if their distance is less than this
     * threshold.
     */
    double contact_dist_thresh;

  protected:
    void calc_impacting_unilateral_constraint_forces(double dt);
    void find_unilateral_constraints(double min_contact_dist);
    void calc_fwd_dyn();
    void calc_compliant_unilateral_constraint_forces();
    void preprocess_constraint(UnilateralConstraint& e);
    void determine_geometries();
    void broad_phase(double dt);
    void calc_pairwise_distances();
    void visualize_contact( UnilateralConstraint& constraint );

    /// Object for handling impact constraints
    ImpactConstraintHandler _impact_constraint_handler;

    /// Object for handling penalty constraints
    PenaltyConstraintHandler _penalty_constraint_handler;
    
    /// The vector of rigid constraints
    std::vector<UnilateralConstraint> _rigid_constraints;

    /// The vector of compliant constraints
    std::vector<UnilateralConstraint> _compliant_constraints;

  protected:

    double calc_CA_step();
    double calc_next_CA_Euler_step(double contact_dist_thresh) const;
    void update_constraint_violations(const std::vector<PairwiseDistInfo>& pairwise_distances);

    /// The constraint stabilization mechanism
    ConstraintStabilization _cstab;

    /// The dissipation mechanism, if any
    boost::shared_ptr<Dissipation> _dissipator;

    /// Pairwise distances at bodies' current configurations
    std::vector<PairwiseDistInfo> _pairwise_distances;

    /// Work vector
    Ravelin::VectorNd _workV;

    /// Interpenetration constraint violation tolerances
    std::map<sorted_pair<CollisionGeometryPtr>, double> _ip_tolerances;

    /// The collision detection mechanism
    boost::shared_ptr<CollisionDetection> _coldet;

    /// The geometries in the simulator
    std::vector<CollisionGeometryPtr> _geometries;

    /// Geometric pairs that should be checked for unilateral constraints (according to broad phase collision detection)
    std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> > _pairs_to_check;
}; // end class

} // end namespace

#endif


