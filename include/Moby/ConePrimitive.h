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
    void set_radius(double radius);
    void set_height(double height);
    void set_circle_points(unsigned n);
    void set_num_rings(unsigned n);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void get_vertices(std::vector<Point3d>& vertices);
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual double calc_dist_and_normal(const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, const Ravelin::Transform3d& thisTb, Point3d& pthis, Point3d& pp) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual osg::Node* create_visualization();
    double calc_dist(const SpherePrimitive* s, Point3d& pcone, Point3d& psph) const;

    /// Gets the number of rings on the cone
    unsigned get_num_rings() const { return _nrings; }

    /// Gets the radius of this cone
    double get_radius() const { return _radius; }

    /// Gets the height of this cone
    double get_height() const { return _height; }

    /// Gets the number of points in a circle on the cone 
    unsigned get_circle_points() const { return _npoints; }
    
  private:
    static double sqr(double x) { return x*x; }
    virtual void calc_mass_properties(); 
    double calc_penetration_depth(boost::shared_ptr<const Ravelin::Pose3d> P, const Point3d& p) const; 
    Ravelin::Vector3d determine_normal(boost::shared_ptr<const Ravelin::Pose3d> P, const Point3d& query) const;

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Map from the geometry to the vector of vertices (w/transform and intersection tolerance applied), if any
    std::map<CollisionGeometryPtr, std::vector<Point3d> > _vertices;

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

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;
};

} // end namespace

#endif

