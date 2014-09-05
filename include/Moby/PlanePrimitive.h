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

class SpherePrimitive;
class BoxPrimitive;

/// Represents a heightmap with height zero on the xz plane (primitive can be transformed)
class PlanePrimitive : public Primitive
{
  friend class CCD;

  public: 
    PlanePrimitive();
    PlanePrimitive(const Ravelin::Pose3d& T);
    virtual Point3d get_supporting_point(const Ravelin::Vector3d& dir) const;
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const;
    void get_vertices(BVPtr bv, boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const;
    virtual double calc_dist_and_normal(const Point3d& point, std::vector<Ravelin::Vector3d>& normals) const;
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const;
    double calc_signed_dist(boost::shared_ptr<const BoxPrimitive> b, Point3d& pthis, Point3d& pb) const;
    double calc_signed_dist(boost::shared_ptr<const SpherePrimitive> s, Point3d& pthis, Point3d& psph) const;
    virtual double calc_signed_dist(const Point3d& p) const;
    virtual osg::Node* create_visualization();
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P);
    virtual void calc_mass_properties() { _density.reset(); _J.set_zero(); }
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual bool is_convex() const { return true; }

  protected:
    virtual double calc_height(const Point3d& p) const;

    /// The bounding volumes for the heightmap 
    std::map<CollisionGeometryPtr, boost::shared_ptr<OBB> > _obbs; 

}; // end class

} // end namespace

#endif
