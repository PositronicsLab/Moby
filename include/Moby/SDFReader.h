/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _SDF_READER_H
#define _SDF_READER_H

#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>
#include <map>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/Types.h>
#include <Moby/EventDrivenSimulator.h>

namespace Moby {

class Simulator;
class RigidBody;
class RCArticulatedBody;
class MCArticulatedBody;
class Primitive;

/// Used to read the simulator state from XML
class SDFReader
{
  public:
    static boost::shared_ptr<EventDrivenSimulator> read(const std::string& fname);
    
  private:
    enum TupleType { eNone, eVectorN, eVector3, eQuat };

    struct SurfaceData
    {
      double epsilon;  // coefficient of restitution
      double mu_c;     // Coulomb friction coefficient
      double mu_v;     // viscous friction coefficient
      unsigned NK;     // edges in friction cone approximation
    };

    static void read_surface(boost::shared_ptr<const XMLTree> node, boost::shared_ptr<SurfaceData>& sd);
    static boost::shared_ptr<EventDrivenSimulator> read_world(boost::shared_ptr<const XMLTree> node);
    static std::vector<DynamicBodyPtr> read_models(boost::shared_ptr<const XMLTree> node, boost::shared_ptr<EventDrivenSimulator> sim);
    static unsigned read_uint(boost::shared_ptr<const XMLTree> node);
    static double read_double(boost::shared_ptr<const XMLTree> node);
    static bool read_bool(boost::shared_ptr<const XMLTree> node);
    static Ravelin::Vector3d read_Vector3(boost::shared_ptr<const XMLTree> node);
    static boost::shared_ptr<const XMLTree> find_subtree(boost::shared_ptr<const XMLTree> root, const std::string& name);
    static boost::shared_ptr<const XMLTree> find_one_tag(const std::string& tag, boost::shared_ptr<const XMLTree> root);
    static std::list<boost::shared_ptr<const XMLTree> > find_tag(const std::string& tag, boost::shared_ptr<const XMLTree> root);
    static void find_tag(const std::string& tag, boost::shared_ptr<const XMLTree> root, std::list<boost::shared_ptr<const XMLTree> >& l);
    static PrimitivePtr read_heightmap(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_plane(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_box(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_sphere(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_cylinder(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_cone(boost::shared_ptr<const XMLTree> node);
    static PrimitivePtr read_trimesh(boost::shared_ptr<const XMLTree> node);
    static DynamicBodyPtr read_model(boost::shared_ptr<const XMLTree> node, std::map<RigidBodyPtr, boost::shared_ptr<SurfaceData> >& sdata);
    static RigidBodyPtr read_link(boost::shared_ptr<const XMLTree> node, boost::shared_ptr<SurfaceData>& sdata);
    static void read_collision_node(boost::shared_ptr<const XMLTree> node, RigidBodyPtr rb, boost::shared_ptr<SurfaceData>& sd);
    static void read_visual_node(boost::shared_ptr<const XMLTree> node, RigidBodyPtr rb);
    static PrimitivePtr read_geometry(boost::shared_ptr<const XMLTree> node);
    static Ravelin::Pose3d read_pose(boost::shared_ptr<const XMLTree> node);
    static Ravelin::SpatialRBInertiad read_inertial(boost::shared_ptr<const XMLTree> node, RigidBodyPtr rb);

    static JointPtr read_joint(boost::shared_ptr<const XMLTree> node, const std::map<std::string, RigidBodyPtr>& link_map, RigidBodyPtr& base_link);
    static TupleType get_tuple(boost::shared_ptr<const XMLTree> node);
}; // end class
} // end namespace

#endif

