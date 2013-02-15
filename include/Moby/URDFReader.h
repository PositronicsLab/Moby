/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
#include <Moby/VectorN.h>
#include <Moby/MatrixN.h>

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

        std::map<RigidBodyPtr, Matrix4> inertia_transforms;
        std::map<RigidBodyPtr, Matrix4> visual_transforms;
        std::map<RigidBodyPtr, Matrix4> collision_transforms;
        std::map<JointPtr, Matrix4> joint_transforms;
        std::map<RigidBodyPtr, void*> visual_transform_nodes; 
        std::map<JointPtr, Vector3> joint_axes;
        std::map<JointPtr, RigidBodyPtr> joint_parent, joint_child;
        std::map<std::string, std::pair<VectorN, std::string> > materials;
    };

    static void find_outboards(const URDFData& data, RigidBodyPtr link, std::vector<std::pair<JointPtr, RigidBodyPtr> >& outboards, std::map<RigidBodyPtr, RigidBodyPtr>& parents);
    static void output_data(const URDFData& data, RigidBodyPtr link);
    static bool valid_transformation(const URDFData& data, RigidBodyPtr parent, JointPtr joint, RigidBodyPtr link);
    static void fix_Moby(URDFData& data, const std::vector<RigidBodyPtr>& links, const std::vector<JointPtr>& joints);
    static JointPtr find_joint(const URDFData& data, RigidBodyPtr outboard_link);
    static void find_children(const URDFData& data, RigidBodyPtr link, std::queue<RigidBodyPtr>& q, std::map<RigidBodyPtr, RigidBodyPtr>& parents);
    static bool read_texture(XMLTreeConstPtr node, URDFData& data, std::string& fname);
    static bool read_color(XMLTreeConstPtr node, URDFData& data, VectorN& color);
    static void read_material(XMLTreeConstPtr node, URDFData& data, void* osg_node);
    static PrimitivePtr read_primitive(XMLTreeConstPtr node, URDFData& data);
    static boost::shared_ptr<TriangleMeshPrimitive> read_trimesh(XMLTreeConstPtr node, URDFData& data);
    static boost::shared_ptr<SpherePrimitive> read_sphere(XMLTreeConstPtr node, URDFData& data);
    static boost::shared_ptr<BoxPrimitive> read_box(XMLTreeConstPtr node, URDFData& data);
    static boost::shared_ptr<CylinderPrimitive> read_cylinder(XMLTreeConstPtr node, URDFData& data);
    static bool transform_frames(URDFData& data, const std::vector<RigidBodyPtr>& links, const std::vector<JointPtr>& joints);
    static Matrix3 read_inertia(XMLTreeConstPtr node, URDFData& data);
    static Real read_mass(XMLTreeConstPtr node, URDFData& data);
    static Matrix4 read_origin(XMLTreeConstPtr node, URDFData& data);
    static void read_collision(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link);
    static void read_visual(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link);
    static void read_inertial(XMLTreeConstPtr node, URDFData& data, RigidBodyPtr link);
    static void read_safety_controller(XMLTreeConstPtr node, URDFData& data, JointPtr joint); 
    static void read_calibration(XMLTreeConstPtr node, URDFData& data, JointPtr joint); 
    static void read_limits(XMLTreeConstPtr node, URDFData& data, JointPtr joint); 
    static void read_dynamics(XMLTreeConstPtr node, URDFData& data, JointPtr joint); 
    static void read_axis(XMLTreeConstPtr node, URDFData& data, JointPtr joint); 
    static RigidBodyPtr read_parent(XMLTreeConstPtr node, URDFData& data, const std::vector<RigidBodyPtr>& links); 
    static RigidBodyPtr read_child(XMLTreeConstPtr node, URDFData& data, const std::vector<RigidBodyPtr>& links); 
    static void read_joint(XMLTreeConstPtr node, URDFData& data, const std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
    static void read_joints(XMLTreeConstPtr node, URDFData& data, const std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
    static void read_links(XMLTreeConstPtr node, URDFData& data, std::vector<RigidBodyPtr>& links); 
    static void read_link(XMLTreeConstPtr node, URDFData& data, std::vector<RigidBodyPtr>& links); 
    static bool read_robot(XMLTreeConstPtr node, URDFData& data, std::string& name, std::vector<RigidBodyPtr>& links, std::vector<JointPtr>& joints); 
}; // end class
} // end namespace

#endif

