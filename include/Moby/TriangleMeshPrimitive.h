/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _TRIANGLE_MESH_PRIMITIVE_H
#define _TRIANGLE_MESH_PRIMITIVE_H

#include <queue>
#include <fstream>
#include <iostream>
#include <cmath>
#include <list>
#include <string>
#include <Moby/Types.h>
#include <Moby/Primitive.h>

namespace Moby {

/// Represents a triangle mesh "primitive" for inertia properties, collision detection, and visualization
class TriangleMeshPrimitive : public Primitive
{
  public: 
    TriangleMeshPrimitive();
    TriangleMeshPrimitive(const std::string& filename, bool center = true);
    TriangleMeshPrimitive(const std::string& filename, const Matrix4& T, bool center = true);
    void set_edge_sample_length(Real len);
    virtual osg::Node* create_visualization();

    /// Gets the length of an edge in the mesh above which point sub-samples are created
    Real get_edge_sample_length() const { return _edge_sample_length; }

    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);  
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual boost::shared_ptr<void> save_state() const;
    virtual void load_state(boost::shared_ptr<void> state);
    virtual BVPtr get_BVH_root();
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices);
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual bool intersect_seg(const Vector3* u, BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh() { return _mesh; }
    virtual void set_deformable(bool flag);
    virtual void set_intersection_tolerance(Real tol);
    void set_mesh(boost::shared_ptr<const IndexedTriArray> mesh);
    virtual void set_transform(const Matrix4& T);

  private:
    void center();
    virtual void calc_mass_properties();

    /// Determines whether we convexify the mesh for inertial calculations
    bool _convexify_inertia;

    /// The root bounding volume around the primitive; can differ based on whether the geometry is deformable
    BVPtr _root;
    
    /// The underlying mesh
    /**
     * \note the mesh changes when the primitive's transform changes
     */
    boost::shared_ptr<const IndexedTriArray> _mesh;

    /// Edge sample length above which pseudo-vertices are added
    Real _edge_sample_length;

    /// Thick triangle structure augmented with mesh data
    class AThickTri : public ThickTriangle
    {
      public:
        AThickTri(const Triangle& tri, Real tol) : ThickTriangle(tri, tol) {}
        boost::shared_ptr<const IndexedTriArray> mesh;  // the mesh that this triangle came from
        unsigned tri_idx;             // the index of this triangle
    };

    struct TriangleMeshPrimitiveState
    {
      boost::shared_ptr<void> pstate;  // state information for the primitive object
      BVPtr rootBVH;     // root of the bounding volume hierarchy
      std::map<BVPtr, std::list<unsigned> > mesh_tris;
      std::map<BVPtr, std::list<unsigned> > mesh_vertices;
      std::map<BVPtr, std::list<boost::shared_ptr<AThickTri> > > tris;
    };

    void construct_mesh_vertices(boost::shared_ptr<const IndexedTriArray> mesh);
    void build_BB_tree();
    void split_tris(const Vector3& point, const Vector3& normal, const IndexedTriArray& orig_mesh, const std::list<unsigned>& ofacets, std::list<unsigned>& pfacets, std::list<unsigned>& nfacets);
    bool split(boost::shared_ptr<const IndexedTriArray> mesh, BVPtr source, BVPtr& tgt1, BVPtr& tgt2, const Vector3& axis);
    static bool is_degen_point_on_tri(boost::shared_ptr<AThickTri> tri, const Vector3& p);

    template <class InputIterator, class OutputIterator>
    static OutputIterator get_vertices(const IndexedTriArray& tris, InputIterator fselect_begin, InputIterator fselect_end, OutputIterator output);

    /// Mapping from BVs to triangles contained within
    std::map<BVPtr, std::list<unsigned> > _mesh_tris;

    /// Vertices used by get_vertices() [and referenced by _mesh_vertices]
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// Mapping from BVs to vertex indices contained within
    std::map<BVPtr, std::list<unsigned> > _mesh_vertices;

    /// Mapping from BV leafs to thick triangles
    std::map<BVPtr, std::list<boost::shared_ptr<AThickTri> > > _tris;

    /// List of triangles covered by a bounding volume
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;

    /// Mapping from Vector3 pointers to mesh vertex indices
    std::map<const Vector3*, unsigned> _mesh_vertex_map;
}; // end class

#include "TriangleMeshPrimitive.inl"

} // end namespace 

#endif

