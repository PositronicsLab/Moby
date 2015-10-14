/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _TORUS_PRIMITIVE_H
#define _TORUS_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class PlanePrimitive;

/// Represents a solid box centered at the origin (by default)
class TorusPrimitive : public Primitive
{
  public:
    TorusPrimitive();
    TorusPrimitive(double major_radius, double minor_radius);
    TorusPrimitive(const Ravelin::Pose3d& T);
    void set_radii(double major_radius, double minor_radius);
    virtual bool is_convex() const { return false; }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const;
    double calc_closest_point(const Point3d& point, Point3d& closest) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P);
    virtual osg::Node* create_visualization();
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const;
    virtual double calc_signed_dist(boost::shared_ptr<const PlanePrimitive> p, Point3d& pthis, Point3d& pp) const;
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& p) const;
    virtual double calc_signed_dist(const Point3d& p) const;
    virtual double get_bounding_radius() const { return _major_radius + _minor_radius; }

    /// Gets the major radius
    double get_major_radius() const { return _major_radius; }

    /// Get the minor radius 
    double get_minor_radius() const { return _minor_radius; }

  private:
    virtual void calc_mass_properties();
    static double urand(double a, double b);

    /// Map from the geometry to the vector of vertices (w/transform and intersection tolerance applied), if any
    std::map<CollisionGeometryPtr, std::vector<Point3d> > _vertices;

    /// The bounding volume (it's an OBB)
    std::map<CollisionGeometryPtr, OBBPtr> _obbs; 

    /// The major radius of the torus 
    double _major_radius;

    /// The minor radius of the torus
    double _minor_radius;
}; // end class
} // end namespace

#endif
