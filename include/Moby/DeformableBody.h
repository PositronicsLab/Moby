/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _DEFORMABLE_BODY_H
#define _DEFORMABLE_BODY_H

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Moby/Node.h>
#include <Moby/AABB.h>
#include <Moby/Visualizable.h>
#include <Moby/SingleBody.h>
#include <Moby/IndexedTetraArray.h>
#include <Moby/IndexedTriArray.h>

namespace Moby {

class TriangleMeshPrimitive;

/// Class for deformable bodies simulated using systems of nodes (currently mass-spring-damper systems) 
/**
 * Each DeformableBody object has a set of nodes, a tetrahedral mesh, a
 * triangular mesh, and a vertex map associated with it. The meshes and vertex
 * mapping correspond to the body in its undeformed state. Each vertex in the
 * triangle mesh maps to a set of barycentric coordinates for a single 
 * tetrahedron in the tetrahedral mesh. The nodes of the body correspond to 
 * the vertices of the tetrahedron. When the nodes move, the tetrahedra deform; 
 * the mapping in the undeformed state informs us how the triangular geometry 
 * deforms as well.
 */
class DeformableBody : public SingleBody
{
  public:
    DeformableBody();
    std::ostream& to_vrml(std::ostream& o, double scale) const;
    virtual void reset_accumulators();
    virtual void transform(const Ravelin::Pose3d& T);
    virtual double calc_kinetic_energy() const;
    virtual unsigned num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const;
    virtual void add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gf);
    virtual void apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gj);
    virtual Ravelin::VectorNd& get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gc);
    virtual Ravelin::VectorNd& get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& gv);
    virtual Ravelin::VectorNd& get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& ga);
    virtual void set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gc);
    virtual void set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& gv);
    virtual Ravelin::MatrixNd& get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::MatrixNd& M);
    virtual Ravelin::VectorNd& get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, Ravelin::VectorNd& f);
    virtual Ravelin::VectorNd& convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Ravelin::Wrenchd& w, Ravelin::VectorNd& gf);
    virtual void set_mesh(boost::shared_ptr<const IndexedTetraArray> tetra_mesh, boost::shared_ptr<Primitive> tri_mesh);
    virtual Ravelin::Vector3d calc_point_vel(const Ravelin::Point3d& p) const;
    virtual void add_wrench(const Ravelin::Wrenchd& w);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual Ravelin::MatrixNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::MatrixNd& B, Ravelin::MatrixNd& X);
    virtual Ravelin::VectorNd& solve_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, const Ravelin::VectorNd& b, Ravelin::VectorNd& x);
    virtual double calc_potential_energy() const = 0;

    /// Gets the pose of the deformable body 
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _F; }

    /// Gets the linear velocity of the center-of-mass of the body
    virtual const Ravelin::Twistd& velocity() const { return _xd; }

    /// Deformable bodies are always enabled
    virtual bool is_enabled() const { return true; } 

    /// Deformable bodies are not part of an articulated body
    virtual ArticulatedBodyPtr get_articulated_body() const { return ArticulatedBodyPtr(); }

    template <class OutputIterator>
    OutputIterator get_all_collision_geometries(OutputIterator output_begin) { if (_geometry) *output_begin++ = _geometry; return output_begin; }

    /// Get the nodes composing the deformable body
    const std::vector<boost::shared_ptr<Node> >& get_nodes() const { return _nodes; }

    /// Gets the collision geometry
    CollisionGeometryPtr get_collision_geometry() const { return _geometry; }

    /// Gets the shared pointer for <b>this</b>
    DeformableBodyPtr get_this() { return boost::dynamic_pointer_cast<DeformableBody>(shared_from_this()); }

  protected:
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_visualization_pose() { return _identity_pose; }
    void calc_com_and_vels();
    void update_geometries();
    unsigned find_closest_tetrahedron(const Ravelin::Point3d& p) const;
    Tetrahedron get_tetrahedron(unsigned i) const;

    template <class OutputIterator>
    OutputIterator get_tetrahedra(const Ravelin::Point3d& p, OutputIterator output_begin) const;

    /// The nodes comprising the deformable body
    std::vector<boost::shared_ptr<Node> > _nodes;

    /// Structure for mapping tetrahedral geometry to triangle geometry
    struct VertexMap
    {
      /// Index of the tetrahedron
      unsigned tetra;   

      /// Barycentric coordinates of the vertex
      double uvw[3];
    };

    /// The tetrahedral mesh for the body's geometry 
    /**
     * \note these are only the indices of the vertices that make up the
     *       tetrahedra
     */
    std::vector<IndexedTetra> _tetrahedra;

    /// The vertex map
    std::vector<VertexMap> _vertex_map;

    /// The velocity of the body (computation frame)
    Ravelin::Twistd _xd;

    /// The inertia matrix (computation frame)
    Ravelin::SpatialRBInertiad _J;

    /// Pose of the body 
    boost::shared_ptr<Ravelin::Pose3d> _F;

  private:
    /// The pose for visualization (always identity)
    boost::shared_ptr<Ravelin::Pose3d> _identity_pose;

    AABBPtr build_AABB_tree(std::map<BVPtr, std::list<unsigned> >& aabb_tetra_map);
    void split_tetra(const Ravelin::Point3d& point, unsigned axis, const std::list<unsigned>& otetra, std::list<unsigned>& ptetra, std::list<unsigned>& ntetra);
    bool split(AABBPtr source, AABBPtr& tgt1, AABBPtr& tgt2, unsigned axis, const std::list<unsigned>& tetra, std::list<unsigned>& ptetra, std::list<unsigned>& ntetra);

    template <class InputIterator, class OutputIterator>
    OutputIterator get_vertices(InputIterator begin, InputIterator end, OutputIterator output_begin) const;

    /// The tetrahedral mesh for the undeformed state of the body
    boost::shared_ptr<const IndexedTetraArray> _tetra_mesh;

    /// The triangle mesh for the undeformed state of the body
    boost::shared_ptr<Primitive> _tri_mesh;

    /// The AABB hierarchy around the deformed tetrahedral mesh
    AABBPtr _hroot;

    /// The mapping from AABB leafs to tetra contained within
    std::map<BVPtr, std::list<unsigned> > _aabb_tetra_map;

    /// The triangle mesh primitive used for the collision geometry
    boost::shared_ptr<TriangleMeshPrimitive> _cgeom_primitive;

    /// The collision geometry
    CollisionGeometryPtr _geometry;

    /// Determines whether configuration updates are enabled/disabled
    bool _config_updates_enabled;

  protected:
    void enable_config_updates();
    void disable_config_updates();
}; // end class

#include <Moby/DeformableBody.inl>

} // end namespace

#endif
