/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
    CSG(const Ravelin::Pose3d& T);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices);
//    virtual void get_vertices(BVPtr bv, std::vector<const Point3d*>& vertices); 
    virtual double calc_dist_and_normal(const Point3d& p, Ravelin::Vector3d& normal) const;
//    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    virtual osg::Node* create_visualization();
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P); 
    void set_operator(BooleanOperation op);
    void set_operand1(PrimitivePtr op1);
    void set_operand2(PrimitivePtr op2);
    void set_mesh(boost::shared_ptr<const IndexedTriArray> mesh);
    PrimitivePtr get_operand1() const { return _op1; }
    PrimitivePtr get_operand2() const { return _op2; }
    BooleanOperation get_operator() const { return _op; }
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, boost::shared_ptr<const Ravelin::Pose3d> pose_this, boost::shared_ptr<const Ravelin::Pose3d> pose_p, Point3d& pthis, Point3d& pb) const;

  private:
//    bool intersect_seg_union(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
//    bool intersect_seg_intersect(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
//    bool intersect_seg_diff(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    virtual void calc_mass_properties();
    void center_mesh();

    /// The axis aligned bounding box for this primitive
    boost::shared_ptr<AABB> _aabb;

    /// The Boolean operation to be carried out upon the two primitives
    BooleanOperation _op;

    /// The two operands
    PrimitivePtr _op1, _op2;

    /// Map from the geometry to the vector of vertices (w/transform and intersection tolerance applied), if any
    std::map<CollisionGeometryPtr, std::vector<Point3d> > _vertices;

    /// A pointer to the computed triangle mesh
    boost::shared_ptr<IndexedTriArray> _mesh;
}; // end class

} // end namespace

#endif
