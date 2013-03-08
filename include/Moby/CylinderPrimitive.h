/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CYLINDER_PRIMITIVE_H
#define _CYLINDER_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class OBB;

/// Defines a cylinder primitive
class CylinderPrimitive : public Primitive
{
  public:
    CylinderPrimitive();
    CylinderPrimitive(Real radius, Real height);
    CylinderPrimitive(Real radius, Real height, unsigned n, unsigned nrings, const Matrix4& T);
    CylinderPrimitive(Real radius, Real height, const Matrix4& T);
    void set_radius(Real radius);
    void set_height(Real height);
    void set_num_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    virtual void set_transform(const Matrix4& T);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual BVPtr get_BVH_root();
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices); 
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv); 
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void set_intersection_tolerance(Real tol);
    virtual osg::Node* create_visualization();

    /// Gets the radius of this cylinder
    Real get_radius() const { return _radius; }

    /// Gets the height of this cylinder
    Real get_height() const { return _height; }

    /// Gets the number of points in a circle on the cylinder 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    Real calc_penetration_depth(const Vector3& p) const;
    unsigned intersect_line(const Vector3& origin, const Vector3& dir, Real& t0, Real& t1) const;
    virtual void calc_mass_properties(); 
 
    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// Radius of the cylinder
    Real _radius;

    /// Height of the cylinder
    Real _height;

    /// Number of points used to create the cylinder
    unsigned _npoints;

    /// Number of rings used to create virtual vertices
    unsigned _nrings;

    /// The bounding volume around the cylinder
    boost::shared_ptr<OBB> _obb;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
};

} // end namespace

#endif
