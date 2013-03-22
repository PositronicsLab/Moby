/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPHERE_PRIMITIVE_H
#define _SPHERE_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class BoundingSphere;

/// Represents a sphere primitive for inertia properties, collision detection, and visualization
class SpherePrimitive : public Primitive
{
  public: 
    SpherePrimitive();
    SpherePrimitive(Real radius);
    SpherePrimitive(const Matrix4& T);
    SpherePrimitive(Real radius, unsigned n);
    SpherePrimitive(Real radius, const Matrix4& T);
    SpherePrimitive(Real radius, unsigned n, const Matrix4& T);
    void set_radius(Real radius);
    void set_num_points(unsigned n);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void set_transform(const Matrix4& T);
    virtual BVPtr get_BVH_root();
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void set_intersection_tolerance(Real tol);
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices);
    virtual osg::Node* create_visualization();

    /// Gets the radius for this sphere
    Real get_radius() const { return _radius; }

    /// Gets the number of points used to create the sphere for visualization / collision checking
    unsigned get_num_points() const { return _npoints; }
    
  private:
    virtual void calc_mass_properties();

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// The bounding volume for the sphere
    boost::shared_ptr<BoundingSphere> _bsph; 

    /// Radius of the sphere primitive
    Real _radius;

    /// Number of points used to create collision geometry
    unsigned _npoints;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
}; // end class

} // end namespace

#endif
