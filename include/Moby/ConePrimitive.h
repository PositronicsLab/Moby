/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CONE_PRIMITIVE_H
#define _CONE_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class OBB;

/// Defines a cone primitive
/**
 * The axis of the untransformed cone lies along the global Y-axis, and the
 * center of the cone is coincident with the global origin.
 */
class ConePrimitive : public Primitive
{
  public:
    ConePrimitive();
    ConePrimitive(Real radius, Real height);
    ConePrimitive(Real radius, Real height, unsigned npoints, unsigned nrings, const Matrix4& T);
    ConePrimitive(Real radius, Real height, const Matrix4& T);
    void set_radius(Real radius);
    void set_height(Real height);
    void set_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual BVPtr get_BVH_root();
    virtual void set_transform(const Matrix4& T);
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void set_intersection_tolerance(Real tol);
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices);

    #ifdef USE_OSG
    virtual osg::Node* create_visualization();
    #endif

    /// Gets the number of rings on the cone
    unsigned get_num_rings() const { return _nrings; }

    /// Gets the radius of this cone
    Real get_radius() const { return _radius; }

    /// Gets the height of this cone
    Real get_height() const { return _height; }

    /// Gets the number of points in a circle on the cone 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    static Real sqr(Real x) { return x*x; }
    virtual void calc_mass_properties(); 
    Real calc_penetration_depth(const Vector3& p) const; 
    Vector3 determine_normal(const Vector3& query) const;

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// The bounding volume for the primitive
    boost::shared_ptr<OBB> _obb; 

    /// Radius of the bottom ring of the cone
    Real _radius;

    /// Height of the cone
    Real _height;

    /// Number of points per ring on the cone
    unsigned _npoints;

    /// Number of rings on the cone (default=1)
    unsigned _nrings;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
};

} // end namespace

#endif

