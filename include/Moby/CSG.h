/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CSG_H
#define _CSG_H

#include <Moby/AABB.h>
#include <Moby/Primitive.h>

namespace Moby {

/// Represents a CSG (constructive solid geometry) for collision, visualization, and inertial properties 
class CSG : public Primitive
{
  public: 
    enum BooleanOperation { eUnion, eIntersection, eDifference };
    CSG();
    CSG(const Matrix4& T);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void set_intersection_tolerance(Real tol);
    virtual BVPtr get_BVH_root();
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices); 
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    #ifdef USE_OSG
    virtual osg::Node* create_visualization();
    #endif
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(); 
    void set_operator(BooleanOperation op);
    void set_operand1(PrimitivePtr op1);
    void set_operand2(PrimitivePtr op2);
    void set_mesh(boost::shared_ptr<const IndexedTriArray> mesh);
    PrimitivePtr get_operand1() const { return _op1; }
    PrimitivePtr get_operand2() const { return _op2; }
    BooleanOperation get_operator() const { return _op; }
    virtual void set_transform(const Matrix4& T);

  private:
    bool intersect_seg_union(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    bool intersect_seg_intersect(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    bool intersect_seg_diff(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual void calc_mass_properties();
    void center_mesh();

    /// The axis aligned bounding box for this primitive
    boost::shared_ptr<AABB> _aabb;

    /// The Boolean operation to be carried out upon the two primitives
    BooleanOperation _op;

    /// The two operands
    PrimitivePtr _op1, _op2;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;

    /// The set of vertices
    boost::shared_ptr<std::vector<Vector3> > _vertices;

    /// A pointer to the computed triangle mesh
    boost::shared_ptr<IndexedTriArray> _mesh;
}; // end class

} // end namespace

#endif
