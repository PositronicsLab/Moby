/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _BOX_PRIMITIVE_H
#define _BOX_PRIMITIVE_H

#include <Moby/Triangle.h>
#include <Moby/Matrix3.h>
#include <Moby/Primitive.h>
#include <Moby/DummyBV.h>

namespace Moby {

/// Represents a solid box centered at the origin (by default)
class BoxPrimitive : public Primitive
{
  public:
    BoxPrimitive();
    BoxPrimitive(Real xlen, Real ylen, Real zlen);
    BoxPrimitive(Real xlen, Real ylen, Real zlen, const Matrix4& T);
    BoxPrimitive(const Matrix4& T);
    void set_size(Real xlen, Real ylen, Real zlen);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual BVPtr get_BVH_root();
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual void set_intersection_tolerance(Real tol);
    virtual void set_transform(const Matrix4& T);
    void set_edge_sample_length(Real len);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void get_vertices(BVPtr, std::vector<const Vector3*>& vertices);
    virtual osg::Node* create_visualization();

    /// Get the x-length of this box
    Real get_x_len() const { return _xlen; }

    /// Get the y-length of this box
    Real get_y_len() const { return _ylen; }

    /// Ge the z-length of this box
    Real get_z_len() const { return _zlen; }

  private:
    enum FaceID { ePOSX, eNEGX, ePOSY, eNEGY, ePOSZ, eNEGZ };
    static void determine_normal(const Vector3& lengths, const Vector3& p, Vector3& normal);
    static bool determine_normal_abs(const Vector3& lengths, const Vector3& p, Vector3& normal);

    virtual void calc_mass_properties();

    /// The maximum edge length for the box
    Real _edge_sample_length;

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// The bounding volume (surprise, it's just a box)
    OBBPtr _obb; 

    /// The box lengths
    Real _xlen, _ylen, _zlen;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;

}; // end class
} // end namespace

#endif
