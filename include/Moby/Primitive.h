/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H

#include <set>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Base.h>
#include <Moby/Triangle.h>
#include <Moby/ThickTriangle.h>
#include <Moby/Constants.h>
#include <Moby/IndexedTriArray.h>

namespace osg {
  class MatrixTransform;
  class Material;
  class Node;
  class Matrixd;
}

namespace Moby {

class CollisionGeometry;

/// Defines a triangle-mesh-based primitive type used for inertial property calculation and geometry provisions
/**
 * The center-of-mass of derived types may be at the origin of the world,
 * or not.  Additionally, Primitive can take a transformation matrix in its
 * constructor, with which the primitive data (com, inertia matrix, and geometry)
 * can be transformed.
 */
class Primitive : public virtual Base
{
  friend class CSG;

  public:
    Primitive();
    Primitive(const Ravelin::Pose3d& T);
    virtual ~Primitive();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    void update_visualization();
    void set_mass(double mass);
    void set_density(double density);
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual Point3d get_supporting_point(const Ravelin::Vector3d& d) const;
    virtual double calc_signed_dist(const Point3d& p) const;
    void add_collision_geometry(CollisionGeometryPtr cg);
    void remove_collision_geometry(CollisionGeometryPtr cg);
    boost::shared_ptr<const Ravelin::Pose3d> get_pose(CollisionGeometryPtr g) const;

    /// Determines whether this primitive is convex
    virtual bool is_convex() const { return false; }

    /// Computes the distance between a point and this primitive
    virtual double calc_dist_and_normal(const Point3d& p, std::vector<Ravelin::Vector3d>& normals) const = 0;

    /// Computes the signed distance between this and another primitive
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const = 0;

     /// Gets the visualization for this primitive
    virtual osg::Node* get_visualization();
    virtual osg::Node* create_visualization() = 0;

    /// Gets the root bounding volume for this primitive
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom) = 0; 

    /// Get vertices corresponding to this primitive
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices) const = 0;

    /// Gets the inertial frame of this primitive
    boost::shared_ptr<const Ravelin::Pose3d> get_inertial_pose() const { return _jF; }

    /// Gets the pose of this primitive 
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _F; } 

    /// Gets the underlying triangle mesh for this primitive 
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P) = 0;

    /// Gets the inertia for this primitive 
    const Ravelin::SpatialRBInertiad& get_inertia() const { return _J; }

  protected:
    virtual void calc_mass_properties() = 0;

    /// The pose of this primitive (relative to the global frame)
    boost::shared_ptr<Ravelin::Pose3d> _F;

    /// The inertial pose of this primitive (relative to the global frame)
    boost::shared_ptr<Ravelin::Pose3d> _jF;

    /// The density of this primitive
    boost::shared_ptr<double> _density;

    /// The inertia of the primitive
    Ravelin::SpatialRBInertiad _J;

  protected:

    /// The poses of this primitive, relative to a collision geometry
    std::map<boost::weak_ptr<CollisionGeometry>, boost::shared_ptr<Ravelin::Pose3d> > _cg_poses;

    /// The poses, relative to a particular collision geometry
    std::set<boost::shared_ptr<Ravelin::Pose3d> > _poses;

  private:

    /// The visualization transform (possibly NULL)
    osg::MatrixTransform* _vtransform;

    /// The material for the visualization (possibly NULL)
    osg::Material* _mat;
}; // end class

} // end namespace

#endif
