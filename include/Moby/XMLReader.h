/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _XML_READER_H
#define _XML_READER_H

#include <map>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/VectorN.h>
#include <Moby/MatrixN.h>

namespace Moby {

class Simulator;
class RigidBody;
class RCArticulatedBody;
class MCArticulatedBody;
class Primitive;

/// Used to read the simulator state from XML
class XMLReader
{
  public:
    static std::map<std::string, BasePtr> read(const std::string& fname);
    
  private:
    enum TupleType { eNone, eVectorN, eVector3, eQuat };
    static XMLTreeConstPtr find_subtree(XMLTreeConstPtr root, const std::string& name);
    static void process_tag(const std::string& tag, XMLTreeConstPtr root, void (*fn)(XMLTreeConstPtr, std::map<std::string, BasePtr>&), std::map<std::string, BasePtr>& id_map);
    static void read_deformable_ccd(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_generalized_ccd(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_c2a_ccd(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_mesh_dcd(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_coldet_plugin(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_box(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_sphere(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_cylinder(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_cone(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_trimesh(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_tetramesh(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_CSG(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_primitive_plugin(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_odepack_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_euler_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_variable_euler_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rk4_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rkf4_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rk4i_integrator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_event_driven_simulator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_simulator(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_ps_deformable_body(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rigid_body(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_mc_abody(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rc_abody(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_rc_abody_symbolic(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_osg_group(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_collision_geometry(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_joint_plugin(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_prismatic_joint(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_fixed_joint(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_revolute_joint(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_spherical_joint(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_universal_joint(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_gravity_force(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_damping_force(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static void read_stokes_drag_force(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    static TupleType get_tuple(XMLTreeConstPtr node);
}; // end class
} // end namespace

#endif

