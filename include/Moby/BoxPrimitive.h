/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _BOX_PRIMITIVE_H
#define _BOX_PRIMITIVE_H

#include <Moby/Triangle.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/DummyBV.h>

namespace Moby {

class SpherePrimitive;

/// Represents a solid box centered at the origin (by default)
class BoxPrimitive : public PolyhedralPrimitive
{
  public:
    BoxPrimitive();
    BoxPrimitive(double xlen, double ylen, double zlen);
    BoxPrimitive(double xlen, double ylen, double zlen, const Ravelin::Pose3d& T);
    BoxPrimitive(const Ravelin::Pose3d& T);
    void set_size(double xlen, double ylen, double zlen);
    virtual unsigned num_facets() const { return 6; }
    virtual void get_facets(boost::shared_ptr<const Ravelin::Pose3d> P, Ravelin::MatrixNd& M, Ravelin::VectorNd& q) const;
    virtual bool is_convex() const { return true; }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const;
    double calc_closest_point(const Point3d& point, Point3d& closest) const;
    virtual void set_pose(const Ravelin::Pose3d& T);
    void set_edge_sample_length(double len);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P);
    virtual osg::Node* create_visualization();
    double calc_signed_dist(boost::shared_ptr<const SpherePrimitive> s, Point3d& pthis, Point3d& psph) const;
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const;
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& p) const;
    virtual double calc_signed_dist(const Point3d& p) const;
    double calc_closest_points(boost::shared_ptr<const SpherePrimitive> s, Point3d& pbox, Point3d& psph) const;

    /// Get the x-length of this box
    double get_x_len() const { return _xlen; }

    /// Get the y-length of this box
    double get_y_len() const { return _ylen; }

    /// Ge the z-length of this box
    double get_z_len() const { return _zlen; }

  private:
    enum FaceID { ePOSX, eNEGX, ePOSY, eNEGY, ePOSZ, eNEGZ };
    static double sqr(double x) { return x*x; }

    virtual void calc_mass_properties();

    /// The maximum edge length for the box
    double _edge_sample_length;

    /// Map from the geometry to the vector of vertices (w/transform and intersection tolerance applied), if any
    std::map<CollisionGeometryPtr, std::vector<Point3d> > _vertices;

    /// The bounding volume (no surprise, it's just a box)
    std::map<CollisionGeometryPtr, OBBPtr> _obbs; 

    /// The box lengths
    double _xlen, _ylen, _zlen;
}; // end class
} // end namespace

#endif
