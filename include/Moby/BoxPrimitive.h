/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _BOX_PRIMITIVE_H
#define _BOX_PRIMITIVE_H

#include <Moby/Triangle.h>
#include <Moby/Primitive.h>
#include <Moby/DummyBV.h>

namespace Moby {

/// Represents a solid box centered at the origin (by default)
class BoxPrimitive : public Primitive
{
  public:
    BoxPrimitive();
    BoxPrimitive(double xlen, double ylen, double zlen);
    BoxPrimitive(double xlen, double ylen, double zlen, const Ravelin::Pose3d& T);
    BoxPrimitive(const Ravelin::Pose3d& T);
    void set_size(double xlen, double ylen, double zlen);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual bool point_inside(BVPtr bv, const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);
    virtual void set_intersection_tolerance(double tol);
    virtual void set_pose(const Ravelin::Pose3d& T);
    void set_edge_sample_length(double len);
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh();
    virtual void get_vertices(BVPtr, std::vector<const Point3d*>& vertices);
    virtual osg::Node* create_visualization();

    /// Get the x-length of this box
    double get_x_len() const { return _xlen; }

    /// Get the y-length of this box
    double get_y_len() const { return _ylen; }

    /// Ge the z-length of this box
    double get_z_len() const { return _zlen; }

  private:
    enum FaceID { ePOSX, eNEGX, ePOSY, eNEGY, ePOSZ, eNEGZ };
    void determine_normal(const Ravelin::Vector3d& lengths, const Point3d& p, boost::shared_ptr<const Ravelin::Pose3d> P, Ravelin::Vector3d& normal) const;
    bool determine_normal_abs(const Ravelin::Vector3d& lengths, const Point3d& p, boost::shared_ptr<const Ravelin::Pose3d> P, Ravelin::Vector3d& normal) const;

    virtual void calc_mass_properties();

    /// The maximum edge length for the box
    double _edge_sample_length;

    /// Pointer to the determined mesh (w/transform applied), if any
    boost::shared_ptr<IndexedTriArray> _mesh;

    /// Pointer to the vector of vertices (w/transform and intersection tolerance applied), if any
    boost::shared_ptr<std::vector<Point3d> > _vertices;

    /// The bounding volume (no surprise, it's just a box)
    std::map<CollisionGeometryPtr, OBBPtr> _obbs; 

    /// The box lengths
    double _xlen, _ylen, _zlen;

    /// The "sub" mesh 
    std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _smesh;

}; // end class
} // end namespace

#endif
