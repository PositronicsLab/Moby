/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _COLLISION_DETECTION_H
#define _COLLISION_DETECTION_H

#include <map>
#include <set>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Pose3d.h>
#include <Moby/sorted_pair>
#include <Moby/Base.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/RigidBody.h>

namespace Moby {

class EventDrivenSimulator;
class ArticulatedBody;
class CollisionGeometry;
class Triangle;
class IndexedTriArray;

/// Structure for holding a pair of colliding triangles
struct CollidingTriPair
{
  CollisionGeometryPtr geom1;    // geometry from which first triangle is taken
  CollisionGeometryPtr geom2;    // geometry from which second triangle is taken
  unsigned tri1;                 // index of first triangle in mesh
  unsigned tri2;                 // index of second triangle in mesh
  const IndexedTriArray* mesh1;  // first triangle mesh
  const IndexedTriArray* mesh2;  // second triangle mesh
};

/// Defines an abstract collision detection mechanism
/**
 * Contact finding and collision detection are two separate, but related, 
 * procedures. Collision detection quickly determines whether two bodies are 
 * geometrically intersecting at their current configurations. Contact finding
 * determines whether two bodies will come into contact (or are already in 
 * contact) in a given time interval; if the bodies do/will contact, contact
 * finding computes contact points and normals. 
 */
class CollisionDetection : public virtual Base
{
  public:
    enum DetectionMode { eFirstContact, eAllContacts };
    CollisionDetection();
    virtual ~CollisionDetection() {}
    void operator=(const CollisionDetection* source);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    void add_dynamic_body(DynamicBodyPtr body);
    void remove_dynamic_body(DynamicBodyPtr body);
    virtual void remove_collision_geometry(CollisionGeometryPtr geom);
    virtual void remove_all_collision_geometries();
    virtual void add_rigid_body(RigidBodyPtr body); 
    virtual void remove_rigid_body(RigidBodyPtr body); 
    virtual void add_articulated_body(ArticulatedBodyPtr abody, bool disable_adjacent);
    virtual void remove_articulated_body(ArticulatedBodyPtr abody);
    virtual void output_object_state(std::ostream& out) const;
    virtual void set_enabled(BasePtr b1, BasePtr b2, bool enabled);
    virtual void set_enabled(BasePtr b, bool enabled);
    bool is_checked(CollisionGeometryPtr cg1, CollisionGeometryPtr cg2) const;

    /// Get the shared pointer for this
    boost::shared_ptr<CollisionDetection> get_this() { return boost::dynamic_pointer_cast<CollisionDetection>(shared_from_this()); }

    /// Determines whether there is a contact between the given pairs of states 
    /**
     * Generalized coordinates are set to q1 on entry; generalized coordinates
     * should be set to q0 on return.
     * \param events the set of determined contacts, on return
     * \return <b>true</b> if there is contact in the time interval, 
     *           <b>false</b> otherwise
     */
    virtual bool is_contact(double dt, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q0, const std::vector<std::pair<DynamicBodyPtr, Ravelin::VectorNd> >& q1, std::vector<UnilateralConstraint>& contacts) = 0;

    /// Adds the specified geometry to the set of geometries checked for collision
    /**
     * \note derived classes will generally need to provide an implementation
     *       of this method (and must call this method 
     *       [CollisionGeometry::add_collision_geometry()] explicitly as well!)
     */
    virtual void add_collision_geometry(CollisionGeometryPtr geom) { _geoms.insert(geom); } 

    /// Gets the set of geometries checked by this collision detector
    const std::set<CollisionGeometryPtr>& get_collision_geometries() const { return _geoms; }
    
    /// Determines whether collision checking for the specified pair is enabled
    bool is_enabled(CollisionGeometryPtr g1, CollisionGeometryPtr g2) const { return disabled_pairs.find(make_sorted_pair(g1, g2)) == disabled_pairs.end(); }
      
    /// Determines whether collision checking for the specified geometry is enabled
    bool is_enabled(CollisionGeometryPtr g) const { return disabled.find(g) == disabled.end(); }

    /// Determines whether there is a collision at the current simulation state with the given tolerance
    /**
     * \param epsilon the tolerance to check for collision; two bodies are 
     *        considered to be in collision if the distance between the 
     *        bodies is less than epsilon
     * \note implementing classes must store the colliding pairs in 
     *         _colliding_pairs
     */
    virtual bool is_collision(double epsilon = 0.0) = 0;

    /// Adds an articulated body to the bodies checked for collision
    void add_articulated_body(ArticulatedBodyPtr abody) { add_articulated_body(abody, disable_adjacent_default); }

    /// Pointer to the event-driven simulator
    boost::weak_ptr<EventDrivenSimulator> simulator;

    /// The default setting for disabling adjacent links in added articulated bodies
    bool disable_adjacent_default;

    /// The pairs of triangles in collision after a call to is_collision()
    std::list<CollidingTriPair> colliding_tris;

    /// The collision detection mode
    /**
     * The mode indicates whether the 
     * collision detector should terminate the detection process as soon as its
     * is determined that two geometries are intersecting (i.e., eFirstContact) 
     * or if all regions of contact (e.g., triangles) should be processed.
     */
    DetectionMode mode;

    /// Determines whether all contacts or first contacts only are reported
    /**
     * The collision detection scheme will generally be much faster if only
     * first contacts are reported, but some simulation models (event driven)
     * may require all contacts to be reported for correct behavior.
     */
    bool return_all_contacts;

    /// The set of disabled pairs of CollisionGeometry objects
    std::set<sorted_pair<CollisionGeometryPtr> > disabled_pairs;

    /// The set of disabled CollisionGeometry objects 
    std::set<CollisionGeometryPtr> disabled;

    /// The set of geometries in collision (from last call to is_collision())
    std::set<sorted_pair<CollisionGeometryPtr> > colliding_pairs;

  protected:

    template <class OutputIterator>
    OutputIterator get_single_bodies(OutputIterator output_begin) const;

    template <class OutputIterator>
    OutputIterator get_rigid_bodies(OutputIterator output_begin) const;

    template <class OutputIterator>
    OutputIterator get_dynamic_bodies(OutputIterator output_begin) const;

    static double calc_distance(CollisionGeometryPtr a, CollisionGeometryPtr b, const Ravelin::Transform3d& aTb, Point3d& cpa, Point3d& cpb); 

    /// The set of geometries checked by the collision detector
    std::set<CollisionGeometryPtr> _geoms;
}; // end class

#include "CollisionDetection.inl"

} // end namespace Moby

#endif

