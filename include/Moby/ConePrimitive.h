/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CONE_PRIMITIVE_H
#define _CONE_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class OBB;
class SpherePrimitive;

/// Defines a cone primitive
/**
 * The axis of the untransformed cone lies along the global Y-axis, and the
 * center of the cone is coincident with the global origin.
 */
class ConePrimitive : public Primitive
{
  public:
    ConePrimitive();
    ConePrimitive(double radius, double height);
    ConePrimitive(double radius, double height, unsigned npoints, unsigned nrings, const Ravelin::Pose3d& T);
    ConePrimitive(double radius, double height, const Ravelin::Pose3d& T);
    virtual bool is_convex() const { return true; }
    void set_radius(double radius);
    void set_height(double height);
    void set_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const;
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const;
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P);
    virtual osg::Node* create_visualization();
    virtual Point3d get_supporting_point(const Ravelin::Vector3d& d) const;
    virtual double calc_signed_dist(const Point3d& p) const;

    /// Gets the number of rings on the cone
    unsigned get_num_rings() const { return _nrings; }

    /// Gets the radius of this cone
    double get_radius() const { return _radius; }

    /// Gets the height of this cone
    double get_height() const { return _height; }

    /// Gets the number of points in a circle on the cone 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    bool intersect_seg(const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    double calc_dist(const SpherePrimitive* s, Point3d& pcone, Point3d& psph) const;
    static double sqr(double x) { return x*x; }
    virtual void calc_mass_properties(); 
    double calc_penetration_depth(const Point3d& p) const; 
    Ravelin::Vector3d determine_normal(const Point3d& query) const;
    bool point_inside(const Point3d& p, Ravelin::Vector3d& normal) const;

    /// The bounding volumes for the primitive, indexed by geometry
    std::map<CollisionGeometryPtr, OBBPtr> _obbs; 

    /// Radius of the bottom ring of the cone
    double _radius;

    /// Height of the cone
    double _height;

    /// Number of points per ring on the cone
    unsigned _npoints;

    /// Number of rings on the cone (default=1)
    unsigned _nrings;
};

} // end namespace

#endif

