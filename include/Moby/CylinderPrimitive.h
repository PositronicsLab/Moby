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
    CylinderPrimitive(double radius, double height);
    CylinderPrimitive(double radius, double height, unsigned n, unsigned nrings, const Ravelin::Pose3d& T);
    CylinderPrimitive(double radius, double height, const Ravelin::Pose3d& T);
    void set_radius(double radius);
    void set_height(double height);
    void set_num_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void get_vertices(BVPtr bv, std::vector<const Point3d*>& vertices); 
    virtual bool point_inside(BVPtr bv, const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv); 
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void set_intersection_tolerance(double tol);
    virtual osg::Node* create_visualization();

    /// Gets the radius of this cylinder
    double get_radius() const { return _radius; }

    /// Gets the height of this cylinder
    double get_height() const { return _height; }

    /// Gets the number of points in a circle on the cylinder 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    double calc_penetration_depth(boost::shared_ptr<const Ravelin::Pose3d> P, const Point3d& p) const;
    unsigned intersect_line(boost::shared_ptr<const Ravelin::Pose3d> P, const Point3d& origin, const Ravelin::Vector3d& dir, double& t0, double& t1) const;
    virtual void calc_mass_properties(); 
 
    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Point3d> > _vertices;

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

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
};

} // end namespace

#endif
