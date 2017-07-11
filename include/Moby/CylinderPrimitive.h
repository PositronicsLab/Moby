/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CYLINDER_PRIMITIVE_H
#define _CYLINDER_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class OBB;
class SpherePrimitive;

/// Defines a cylinder primitive
class CylinderPrimitive : public Primitive
{
  public:
    CylinderPrimitive();
    CylinderPrimitive(double radius, double height);
    CylinderPrimitive(double radius, double height, unsigned n, unsigned nrings, const Ravelin::Pose3d& T);
    CylinderPrimitive(double radius, double height, const Ravelin::Pose3d& T);
    bool is_convex() const override { return true; }
    void set_radius(double radius);
    void set_height(double height);
    void set_num_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    void set_pose(const Ravelin::Pose3d& T) override;
    void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map) override;
    void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const override;
    BVPtr get_BVH_root(CollisionGeometryPtr geom) override;
    void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const override;
    double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const override;
    double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const override;
    boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P) override;
    osg::Node* create_visualization() override;
    Point3d get_supporting_point(const Ravelin::Vector3d& d) const override;
    double calc_signed_dist(const Point3d& p) const override;
    double get_bounding_radius() const override { return _radius + _height/2.0; }
    double get_maximum_compliant_layer_depth() const override { throw std::runtime_error("CylinderPrimitive::get_compliant_layer_depth() not implemented."); }

    /// Gets the radius of this cylinder
    double get_radius() const { return _radius; }

    /// Gets the height of this cylinder
    double get_height() const { return _height; }

    /// Gets the number of points in a circle on the cylinder 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    bool intersect_seg(const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    bool point_inside(const Point3d& p, Ravelin::Vector3d& normal) const;
    double calc_dist(const SpherePrimitive* s, Point3d& pcyl, Point3d& psph) const;
    double calc_penetration_depth(const Point3d& p) const;
    unsigned intersect_line(const Point3d& origin, const Ravelin::Vector3d& dir, double& t0, double& t1) const;
    void calc_mass_properties() override;
 
    /// Radius of the cylinder
    double _radius;

    /// Height of the cylinder
    double _height;

    /// Number of points used to create the cylinder
    unsigned _npoints;

    /// Number of rings used to create virtual vertices
    unsigned _nrings;

    /// The bounding volumes around the cylinder, indexed by geometry
    std::map<CollisionGeometryPtr, OBBPtr> _obbs;
};

} // end namespace

#endif
