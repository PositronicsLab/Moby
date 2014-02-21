/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPHERE_PRIMITIVE_H
#define _SPHERE_PRIMITIVE_H

#include <Moby/Primitive.h>

namespace Moby {

class BoxPrimitive;
class BoundingSphere;

/// Represents a sphere primitive for inertia properties, collision detection, and visualization
class SpherePrimitive : public Primitive
{
  public: 
    SpherePrimitive();
    SpherePrimitive(double radius);
    SpherePrimitive(const Ravelin::Pose3d& T);
    SpherePrimitive(double radius, unsigned n);
    SpherePrimitive(double radius, const Ravelin::Pose3d& T);
    SpherePrimitive(double radius, unsigned n, const Ravelin::Pose3d& T);
    void set_radius(double radius);
    void set_num_points(unsigned n);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual void get_vertices(std::vector<Point3d>& vertices);
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual double calc_dist_and_normal(const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, boost::shared_ptr<const Ravelin::Pose3d> pose_this, boost::shared_ptr<const Ravelin::Pose3d> pose_p, Point3d& pthis, Point3d& pp) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual osg::Node* create_visualization();
    double calc_signed_dist(boost::shared_ptr<const SpherePrimitive> s, boost::shared_ptr<const Ravelin::Pose3d> pose_this, boost::shared_ptr<const Ravelin::Pose3d> pose_s, Point3d& pthis, Point3d& psph) const;

    /// Gets the radius for this sphere
    double get_radius() const { return _radius; }

    /// Gets the number of points used to create the sphere for visualization / collision checking
    unsigned get_num_points() const { return _npoints; }
    
  private:
    virtual void calc_mass_properties();

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Map from the geometry to the vector of vertices (w/transform and intersection tolerance applied), if any
    std::map<CollisionGeometryPtr, std::vector<Point3d> > _vertices;

    /// The bounding volumes for the sphere
    std::map<CollisionGeometryPtr, boost::shared_ptr<BoundingSphere> > _bsphs; 

    /// Radius of the sphere primitive
    double _radius;

    /// Number of points used to create collision geometry
    unsigned _npoints;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
}; // end class

} // end namespace

#endif
