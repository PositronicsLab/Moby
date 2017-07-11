/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _PLANE_PRIMITIVE_H
#define _PLANE_PRIMITIVE_H

#include <Moby/Primitive.h>
#include <Moby/OBB.h>

namespace Moby {

class CylinderPrimitive;
class SpherePrimitive;
class BoxPrimitive;

/// Represents a heightmap with height zero on the xz plane (primitive can be transformed)
class PlanePrimitive : public Primitive
{
  friend class CCD;

  public: 
    PlanePrimitive();
    PlanePrimitive(const Ravelin::Pose3d& T);
    Point3d get_supporting_point(const Ravelin::Vector3d& dir) const override;
    void set_pose(const Ravelin::Pose3d& T) override;
    void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const override;
    void get_vertices(BVPtr bv, boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const;
    double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const override;
    double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const override;
    double calc_signed_dist(boost::shared_ptr<const CylinderPrimitive> s, Point3d& pthis, Point3d& psph) const;
    double calc_signed_dist(boost::shared_ptr<const SpherePrimitive> s, Point3d& pthis, Point3d& psph) const;
    double calc_signed_dist(const Point3d& p) const override;
    osg::Node* create_visualization() override;
    boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P) override;
    void calc_mass_properties() override { _density.reset(); _J.set_zero(); }
    BVPtr get_BVH_root(CollisionGeometryPtr geom) override;
    void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) override;
    void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const override;
    bool is_convex() const override { return true; }
    double get_bounding_radius() const override { return 0.0; }
    double get_maximum_compliant_layer_depth() const override { return get_compliant_layer_depth(); }

  protected:
    double calc_height(const Point3d& p) const;

    /// The bounding volumes for the heightmap 
    std::map<CollisionGeometryPtr, boost::shared_ptr<OBB> > _obbs; 

}; // end class

} // end namespace

#endif
