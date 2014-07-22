/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _URDF_READER_H
#define _URDF_READER_H

#ifdef USE_OSG
#include <osg/Node>
#endif
#include <map>
#include <string>
#include <queue>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>

namespace Moby {

class BoxPrimitive;
class SpherePrimitive;
class CylinderPrimitive;
class TriangleMeshPrimitive;
class Simulator;
class RigidBody;
class RCArticulatedBody;
class MCArticulatedBody;
class Primitive;

/// Used to read the simulator state from URDF
class URDFReader
{
  public:
    static bool read(const std::string& fname, std::string& name, std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints);
    
  private:
    class URDFData
    {
      public:
        ~URDFData()
        {
          #ifdef USE_OSG
          for (std::map<RigidBodyPtr, void*>::const_iterator i = visual_transform_nodes.begin(); i != visual_transform_nodes.end(); i++)
            ((osg::Node*) i->second)->unref();
          #endif
        }

        std::map<RigidBodyPtr, void*> visual_transform_nodes;
        std::map<JointPtr, RigidBodyPtr> joint_parent, joint_child;
        std::map<std::string, std::pair<Ravelin::VectorNd, std::string> > materials;
    };

    static void find_outboards(const URDFData& data, RigidBodyPtr link, std::vector<std::pair<JointPtr, RigidBodyPtr> >& outboards, std::map<RigidBodyPtr, RigidBodyPtr>& parents);
    static void output_data(const URDFData& data, RigidBodyPtr link);
    static JointPtr find_joint(const URDFData& data, RigidBodyPtr outboard_link);
    static void find_children(const URDFData& data, RigidBodyPtr link, std::queue<RigidBodyPtr>& q, std::map<RigidBodyPtr, RigidBodyPtr>& parents);
    static bool read_texture(boost::shared_ptr<const XMLTree> node, URDFData& data, std::string& fname);
    static bool read_color(boost::shared_ptr<const XMLTree> node, URDFData& data, Ravelin::VectorNd& color);
    static void read_material(boost::shared_ptr<const XMLTree> node, URDFData& data, void* osg_node);
    static PrimitivePtr read_primitive(boost::shared_ptr<const XMLTree> node, URDFData& data);
//    static boost::shared_ptr<TriangleMeshPrimitive> read_trimesh(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static boost::shared_ptr<SpherePrimitive> read_sphere(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static boost::shared_ptr<BoxPrimitive> read_box(boost::shared_ptr<const XMLTree> node, URDFData& data);
//    static boost::shared_ptr<CylinderPrimitive> read_cylinder(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static Ravelin::Matrix3d read_inertia(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static double read_mass(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static Ravelin::Pose3d read_origin(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static void read_collision(boost::shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link);
    static void read_visual(boost::shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link);
    static void read_inertial(boost::shared_ptr<const XMLTree> node, URDFData& data, RigidBodyPtr link);
    static void read_safety_controller(boost::shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint); 
    static void read_calibration(boost::shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint); 
    static void read_limits(boost::shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint); 
    static void read_dynamics(boost::shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint); 
    static void read_axis(boost::shared_ptr<const XMLTree> node, URDFData& data, JointPtr joint); 
    static RigidBodyPtr read_parent(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<RigidBodyPtr>& links); 
    static RigidBodyPtr read_child(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<RigidBodyPtr>& links); 
    static void read_joint(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
    static void read_joints(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
    static void read_links(boost::shared_ptr<const XMLTree> node, URDFData& data, std::vector<RigidBodyPtr>& links); 
    static void read_link(boost::shared_ptr<const XMLTree> node, URDFData& data, std::vector<RigidBodyPtr>& links); 
    static bool read_robot(boost::shared_ptr<const XMLTree> node, URDFData& data, std::string& name, std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
}; // end class
} // end namespace

#endif

