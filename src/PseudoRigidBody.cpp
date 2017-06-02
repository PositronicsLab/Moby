#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/LinAlgd.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/XMLTree.h>
#include <Moby/PseudoRigidBody.h>

using std::map;
using std::list;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

PseudoRigidBody::PseudoRigidBody()
{
  // Must setup visualization and collision geometry here; inertia will be
  // read in with the rigid body. 

  // Setup initial point mass locations.
  const double XC = 0.5;
  const double YC = 0.5;
  const double ZC = 0.5;
  const double DD = 1;
  std::vector<PointMass>& pms = _point_mass_states;

/*
  // This is the top of the rigid core (view is from y = positive infinity).
  // .-----.
  // |     |
  // |     |
  // .-----.
  //
  // ^ -z
  // |
  // |
  // -----> +x

  // The compliant layer will be extended, and the rigid core will no longer
  // be visible when we add vertices- without subdivisions:
  // * *     * *
  // * *     * *
  //   
  //   
  // * *     * *
  // * *     * *
  // Springs will connect each adjacent pair of vertices.

  // Pistons will be added here:
  // * *     * *
  // *-|     |-*
  //   
  //   
  // *-|     |-*
  // * *     * *
  
  // Vertices will be labeled as:
  // A B     E F
  // C D     H G
  //   
  //   
  // M N     I J
  // P O     L K
  // There will be four such sets.

  // Default dimensions of the rigid core are a 1m x 1m x 1m cube.
  // Default distance between A/B, A/C, etc. is 1cm.
  
  // Plane YC + DD
  pms.push_back(PointMass(Point3d(-XC - DD, YC + DD, -ZC - DD))); // A
  pms.push_back(PointMass(Point3d(-XC,      YC + DD, -ZC - DD))); // B
  pms.push_back(PointMass(Point3d(-XC - DD, YC + DD, -ZC     ))); // C
  pms.push_back(PointMass(Point3d(-XC     , YC + DD, -ZC     ))); // D 
   
  pms.push_back(PointMass(Point3d(+XC     , YC + DD, -ZC - DD))); // E
  pms.push_back(PointMass(Point3d(+XC + DD, YC + DD, -ZC - DD))); // F
  pms.push_back(PointMass(Point3d(+XC + DD, YC + DD, -ZC     ))); // G
  pms.push_back(PointMass(Point3d(+XC     , YC + DD, -ZC     ))); // H 
  
  pms.push_back(PointMass(Point3d(+XC     , YC + DD, +ZC     ))); // I 
  pms.push_back(PointMass(Point3d(+XC + DD, YC + DD, +ZC     ))); // J
  pms.push_back(PointMass(Point3d(+XC + DD, YC + DD, +ZC + DD))); // K
  pms.push_back(PointMass(Point3d(+XC     , YC + DD, +ZC + DD))); // L 
   
  pms.push_back(PointMass(Point3d(-XC - DD, YC + DD, +ZC     ))); // M 
  pms.push_back(PointMass(Point3d(-XC     , YC + DD, +ZC     ))); // N
  pms.push_back(PointMass(Point3d(-XC     , YC + DD, +ZC + DD))); // O
  pms.push_back(PointMass(Point3d(-XC - DD, YC + DD, +ZC + DD))); // P 

  // Plane YC
  pms.push_back(PointMass(Point3d(-XC - DD, YC, -ZC - DD))); // A
  pms.push_back(PointMass(Point3d(-XC,      YC, -ZC - DD))); // B
  pms.push_back(PointMass(Point3d(-XC - DD, YC, -ZC     ))); // C
  pms.push_back(PointMass(Point3d(-XC     , YC, -ZC     ))); // D 
   
  pms.push_back(PointMass(Point3d(+XC     , YC, -ZC - DD))); // E
  pms.push_back(PointMass(Point3d(+XC + DD, YC, -ZC - DD))); // F
  pms.push_back(PointMass(Point3d(+XC + DD, YC, -ZC     ))); // G
  pms.push_back(PointMass(Point3d(+XC     , YC, -ZC     ))); // H 
  
  pms.push_back(PointMass(Point3d(+XC     , YC, +ZC     ))); // I 
  pms.push_back(PointMass(Point3d(+XC + DD, YC, +ZC     ))); // J
  pms.push_back(PointMass(Point3d(+XC + DD, YC, +ZC + DD))); // K
  pms.push_back(PointMass(Point3d(+XC     , YC, +ZC + DD))); // L 
   
  pms.push_back(PointMass(Point3d(-XC - DD, YC, +ZC     ))); // M 
  pms.push_back(PointMass(Point3d(-XC     , YC, +ZC     ))); // N
  pms.push_back(PointMass(Point3d(-XC     , YC, +ZC + DD))); // O
  pms.push_back(PointMass(Point3d(-XC - DD, YC, +ZC + DD))); // P 

  // Plane -YC
  pms.push_back(PointMass(Point3d(-XC - DD, -YC, -ZC - DD))); // A
  pms.push_back(PointMass(Point3d(-XC,      -YC, -ZC - DD))); // B
  pms.push_back(PointMass(Point3d(-XC - DD, -YC, -ZC     ))); // C
  pms.push_back(PointMass(Point3d(-XC     , -YC, -ZC     ))); // D 
   
  pms.push_back(PointMass(Point3d(+XC     , -YC, -ZC - DD))); // E
  pms.push_back(PointMass(Point3d(+XC + DD, -YC, -ZC - DD))); // F
  pms.push_back(PointMass(Point3d(+XC + DD, -YC, -ZC     ))); // G
  pms.push_back(PointMass(Point3d(+XC     , -YC, -ZC     ))); // H 
  
  pms.push_back(PointMass(Point3d(+XC     , -YC, +ZC     ))); // I 
  pms.push_back(PointMass(Point3d(+XC + DD, -YC, +ZC     ))); // J
  pms.push_back(PointMass(Point3d(+XC + DD, -YC, +ZC + DD))); // K
  pms.push_back(PointMass(Point3d(+XC     , -YC, +ZC + DD))); // L 
   
  pms.push_back(PointMass(Point3d(-XC - DD, -YC, +ZC     ))); // M 
  pms.push_back(PointMass(Point3d(-XC     , -YC, +ZC     ))); // N
  pms.push_back(PointMass(Point3d(-XC     , -YC, +ZC + DD))); // O
  pms.push_back(PointMass(Point3d(-XC - DD, -YC, +ZC + DD))); // P 

  // Plane -YC - DD
  pms.push_back(PointMass(Point3d(-XC - DD, -YC - DD, -ZC - DD))); // A
  pms.push_back(PointMass(Point3d(-XC,      -YC - DD, -ZC - DD))); // B
  pms.push_back(PointMass(Point3d(-XC - DD, -YC - DD, -ZC     ))); // C
  pms.push_back(PointMass(Point3d(-XC     , -YC - DD, -ZC     ))); // D 
   
  pms.push_back(PointMass(Point3d(+XC     , -YC - DD, -ZC - DD))); // E
  pms.push_back(PointMass(Point3d(+XC + DD, -YC - DD, -ZC - DD))); // F
  pms.push_back(PointMass(Point3d(+XC + DD, -YC - DD, -ZC     ))); // G
  pms.push_back(PointMass(Point3d(+XC     , -YC - DD, -ZC     ))); // H 
  
  pms.push_back(PointMass(Point3d(+XC     , -YC - DD, +ZC     ))); // I 
  pms.push_back(PointMass(Point3d(+XC + DD, -YC - DD, +ZC     ))); // J
  pms.push_back(PointMass(Point3d(+XC + DD, -YC - DD, +ZC + DD))); // K
  pms.push_back(PointMass(Point3d(+XC     , -YC - DD, +ZC + DD))); // L 
   
  pms.push_back(PointMass(Point3d(-XC - DD, -YC - DD, +ZC     ))); // M 
  pms.push_back(PointMass(Point3d(-XC     , -YC - DD, +ZC     ))); // N
  pms.push_back(PointMass(Point3d(-XC     , -YC - DD, +ZC + DD))); // O
  pms.push_back(PointMass(Point3d(-XC - DD, -YC - DD, +ZC + DD))); // P 

  // Prepare to set up pistons
  Vector3d u_x(1,0,0), v1_x, v2_x;
  Vector3d u_z(0,0,1), v1_z, v2_z;
  Vector3d::determine_orthonormal_basis(u_x, v1_x, v2_x);
  Vector3d::determine_orthonormal_basis(u_z, v1_z, v2_z);

  // -- Set up pistons in the two y planes y = +/- YC.
  // Plane YC
  _constraints.push_back(PistonConstraint(-u_x, v1_x, -v2_x, 
     Point3d(-XC, YC, -ZC), 18));  // D -> C
  _constraints.push_back(PistonConstraint(-u_z, v1_z, -v2_z, 
     Point3d(-XC, YC, -ZC), 17));  // D -> B 
  _constraints.push_back(PistonConstraint(u_x, v1_x, v2_x, 
     Point3d(XC, YC, -ZC), 22));  // H -> G
  _constraints.push_back(PistonConstraint(-u_z, v1_z, -v2_z, 
     Point3d(XC, YC, -ZC), 20));  // H -> E 
  _constraints.push_back(PistonConstraint(-u_x, v1_x, -v2_x, 
     Point3d(-XC, YC, ZC), 28));  // N -> M
  _constraints.push_back(PistonConstraint(u_z, v1_z, v2_z, 
     Point3d(-XC, YC, ZC), 30));  // N -> O 
  _constraints.push_back(PistonConstraint(u_x, v1_x, v2_x, 
     Point3d(XC, YC, ZC), 25));  // I -> J
  _constraints.push_back(PistonConstraint(u_z, v1_z, v2_z, 
     Point3d(XC, YC, ZC), 27));  // I -> L 

  // Plane -YC
  _constraints.push_back(PistonConstraint(-u_x, v1_x, -v2_x, 
     Point3d(-XC, -YC, -ZC), 34));  // D -> C
  _constraints.push_back(PistonConstraint(-u_z, v1_z, -v2_z, 
     Point3d(-XC, -YC, -ZC), 33));  // D -> B 
  _constraints.push_back(PistonConstraint(u_x, v1_x, v2_x, 
     Point3d(XC, -YC, -ZC), 38));  // H -> G
  _constraints.push_back(PistonConstraint(-u_z, v1_z, -v2_z, 
     Point3d(XC, -YC, -ZC), 36));  // H -> E 
  _constraints.push_back(PistonConstraint(-u_x, v1_x, -v2_x, 
     Point3d(-XC, -YC, ZC), 44));  // N -> M
  _constraints.push_back(PistonConstraint(u_z, v1_z, v2_z, 
     Point3d(-XC, -YC, ZC), 46));  // N -> O 
  _constraints.push_back(PistonConstraint(u_x, v1_x, v2_x, 
     Point3d(XC, -YC, ZC), 41));  // I -> J
  _constraints.push_back(PistonConstraint(u_z, v1_z, v2_z, 
     Point3d(XC, -YC, ZC), 43));  // I -> L 

  // Set up vertical pistons (emanating from y = YC)
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(XC, YC, ZC), 3));  // D -> (YC+DD:D)
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(XC, YC, ZC), 7));  // H -> (YC+DD:H)
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(XC, YC, ZC), 8));  // I -> (YC+DD:I)
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(XC, YC, ZC), 13))  // N -> (YC+DD:N)

  // Set up vertical pistons (emanating from y = -YC)
  _constraints.push_back(PistonConstraint(-u_y, v1_y, -v2_y, 
     Point3d(XC, -YC, ZC), 195));  // D -> (YC-DD:D)
  _constraints.push_back(PistonConstraint(-u_y, v1_y, -v2_y, 
     Point3d(XC, -YC, ZC), 199));  // H -> (YC-DD:H)
  _constraints.push_back(PistonConstraint(-u_y, v1_y, -v2_y, 
     Point3d(XC, -YC, ZC), 200));  // I -> (YC-DD:I)
  _constraints.push_back(PistonConstraint(-u_y, v1_y, -v2_y, 
     Point3d(XC, -YC, ZC), 205))  // N -> (YC-DD:N)

  // TODO: set spring connections 
      
  // TODO: set spring resting lengths

  // TODO: set all visualization geometry. There will be eight boxes (one at
  //       each corner) plus twelve boxes (one at each edge) = 20 boxes.

  // TODO: get the rigid body core visualization geometry.
*/

  // Simplified box: only top face is compliant. Note that these points are
  // in the global frame, so we have to anticipate that the pseudo-rigid body
  // box has been pushed upward.
  pms.push_back(PointMass(Point3d(-XC, YC*2 + DD, -ZC)));
  pms.push_back(PointMass(Point3d(-XC, YC*2 + DD, +ZC)));
  pms.push_back(PointMass(Point3d(+XC, YC*2 + DD, +ZC)));
  pms.push_back(PointMass(Point3d(+XC, YC*2 + DD, -ZC)));
 
  // Prepare to set up pistons.
  Vector3d u_y(0,1,0), v1_y, v2_y;
  Vector3d::determine_orthonormal_basis(u_y, v1_y, v2_y);

  // Set up vertical pistons (emanating from y = YC)
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(-XC, YC, -ZC), 0));
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(-XC, YC, +ZC), 1));
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(+XC, YC, +ZC), 2));
  _constraints.push_back(PistonConstraint(u_y, v1_y, v2_y, 
     Point3d(+XC, YC, -ZC), 3));

  // Set spring stiffness and damping.
  // TODO: check and adjust units
  _spring_k = 1;  
  _spring_c = .01; 
  
  // Set up visualization.
  osg::Group* group = get_visualization_data(); 

  // 1. Create the visualization for the rigid core.
  osg::Box* box = new osg::Box;
  osg::Vec3 half_lens(XC, YC, ZC);
  box->setHalfLengths(half_lens);
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(box));
  group->addChild(geode);

  // 2. Create the visualization for the deformable vertices.
  // 2a. Create the geode nonsense first.
  geode = new osg::Geode;
  osg::Geometry* geom = new osg::Geometry;
  geode->addDrawable(geom);
  group->addChild(geode);

  // 2b. Create the vertex array.
// TODO: replace once vertices on the box no longer necessary (see below).
//  _vert_array = new osg::Vec3Array(_point_mass_states.size()); 
  const unsigned X = 0, Y = 1, Z = 2;
 _vert_array = new osg::Vec3Array(_point_mass_states.size() + 4); 
  update_visualization();

  // TODO: remove this once vertices on the top corners of the box no longer
  //       necessary.
  // 2b*. Create more vertices corresponding to the corners of the box.
  //      Presumably these will not be necessary as the compliant layer fully
  //      surrounds the box.
  (*_vert_array)[4] = osg::Vec3((float) -XC, (float) YC, (float) -ZC);
  (*_vert_array)[5] = osg::Vec3((float) -XC, (float) YC, (float) +ZC);
  (*_vert_array)[6] = osg::Vec3((float) +XC, (float) YC, (float) +ZC);
  (*_vert_array)[7] = osg::Vec3((float) +XC, (float) YC, (float) -ZC);

  // Point the geometry to the vertex array.
  geom->setVertexArray(_vert_array);

  // 2c. Create the faces.
  auto add_face = [this, geom](int a, int b, int c, int d)
  {
    osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS, 0);
    face->push_back(a); 
    face->push_back(b); 
    face->push_back(c);
    face->push_back(d);
    _osg_faces.push_back(face);
    geom->addPrimitiveSet(face);
  };
  add_face(0, 1, 2, 3);  // bottom face
  add_face(0, 3, 7, 4);  // vertical face at z = -ZC
  add_face(1, 2, 6, 5);  // vertical face at z = +ZC
  add_face(2, 3, 7, 6);  // vertical face at x = +XC
  add_face(0, 1, 5, 4);  // vertical face at x = -XC

  // 3. Update the poses (if possible).
  update_poses();

  // 4. Construct piston springs.
  _springs.emplace_back(0, DD);
  _springs.emplace_back(1, DD/2);
  _springs.emplace_back(2, DD);
  _springs.emplace_back(3, DD/2);

  // 5. Construct springs between point masses.
  _springs.emplace_back(0, 1, DD);
  _springs.emplace_back(1, 2, DD);
  _springs.emplace_back(2, 3, DD);
  _springs.emplace_back(3, 0, DD);
}

// TODO: remove this function once rigid body and constraints are constructed
// simultaneously.
void PseudoRigidBody::update_poses()
{
  if (!_core)
    return;

  for (unsigned i=0; i< _constraints.size(); i++)
  {
    _constraints[i].u.pose = _core->get_pose();
    _constraints[i].v1.pose = _core->get_pose();
    _constraints[i].v2.pose = _core->get_pose();
    _constraints[i].p.pose = _core->get_pose();
  }
}

void PseudoRigidBody::load_from_xml(shared_ptr<const XMLTree> node, map<std::string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;

  // Load parent data.
  ControlledBody::load_from_xml(node, id_map);

  // Get the name of the rigid body core ID.
  XMLAttrib* core_id_attr = node->get_attrib("rigid-core-id");
  if (core_id_attr)
  {
    // get the ID
    const std::string& ID = core_id_attr->get_string_value();

    // attempt to find the ID
    if ((id_iter = id_map.find(ID)) == id_map.end())
    {
      std::cerr << "PseudoRigidBody::load_from_xml() - RigidBody id: ";
      std::cerr << ID << " not found!" << std::endl << "  offending node: ";
      std::cerr << std::endl << *node;
      return;
    }

    // get the rigid body
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(id_iter->second);
    if (!rb)
    {
      std::cerr << "PseudoRigidBody::load_from_xml() - RigidBody id: ";
      std::cerr << ID << " is not a rigid body!" << std::endl;
      std::cerr << "  offending node: " << std::endl << *node;
      return;
    }

    // save the core
    _core = rb;
  }

  update_poses();
}

// TODO: Implement me
void PseudoRigidBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  assert(false);
}

// Transforms the force applied at a piston to a generalized force.
// A positive force applied at the piston should push the bodies apart (this
// behavior is identical to how increasing distance between bodies maps to
// values of piston positions closer to positive infinity).
VectorNd& PseudoRigidBody::transform_piston_force_to_generalized_force(unsigned piston_index, double f, VectorNd& gf) const
{
  const PistonConstraint& constraint = _constraints[piston_index];

  // Get the index of the point mass.
  const unsigned pm_idx = constraint.point_mass_index;

  // Clear the generalized force
  gf.set_zero();

  // Get the vector describing the axis direction in the rigid core frame and
  // transform it to the point mass force frame.
  Vector3d f_global = Pose3d::transform_vector(GLOBAL, constraint.u) * f;
  const unsigned SPATIAL_D = 6, THREE_D = 3;
  const unsigned index_start = SPATIAL_D + THREE_D * pm_idx;
  gf.segment(index_start, index_start + THREE_D) =  f_global;

  // Apply the reaction force to the rigid body.
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = Origin3d(Pose3d::transform_point(GLOBAL, constraint.p));
  SForced sforce(P);
  sforce.set_force(-f_global);
  gf.segment(0, THREE_D) = sforce.get_force();
  gf.segment(THREE_D, SPATIAL_D) = sforce.get_torque();
  return gf;
}

/*
// TODO: We need to compute ALL piston positions and velocities.
// Computes the connecting piston positions.
void PseudoRigidBody::compute_connecting_piston_positions_and_velocities()
{
  for (unsigned i=0; i< _connecting_piston_indices.size(); i++)
  {
    // Get the index of the connecting piston in the set of piston indices.
    unsigned j = _connecting_piston_indices[i];

    // Get the indices of the key pistons connected by each connecting piston.
    unsigned k1 = _connecting_piston_key_piston_indices[i].first;
    unsigned k2 = _connecting_piston_key_piston_indices[i].second;

    // Compute the position (the compression/extension) of the connecting piston
    // TODO: fix this to account for the rigid core being one of the bodies.
    const Point3d& a = _mass_locations[k1];
    const Point3d& b = _mass_locations[k2];
    const Vector3d& adot = _mass_velocities[k1];
    const Vector3d& bdot = _mass_velocities[k2];
    double len = (a - b).norm(); 

    // NOTE: the distance accounts for the Euclidean distance between the
    //       two attachment points. This distance will only be equal to the
    //       piston length if C(q) = 0. Similarly, the time derivative of that
    //       quantity will only equal the piston velocity if C(q, \dot{q}) = 0. 

    // Compute the velocity of the connecting piston.
    // d/dt ((a - b) dot (a - b)) ^ 1/2 = 
    //       ((a - b) dot (a - b)) ^ -1/2 * d/dt ((a - b) dot (a - b))
    //                                    = 2(a - b) dot d/dt (a - b)
    double vel = 1/len * 2*(a - b).dot(adot - bdot);

    // Set the piston position and velocity.
    _piston_positions[j] = len;
    _piston_velocities[j] = vel;
  }
}
*/

/// TODO: implement me.
MatrixNd& PseudoRigidBody::calc_jacobian(shared_ptr<const Pose3d> source_pose, shared_ptr<const Pose3d> target_pose, shared_ptr<DynamicBodyd> body, MatrixNd& J)
{
  assert(false);
  return J;
}

/// TODO: implement me.
MatrixNd& PseudoRigidBody::calc_jacobian_dot(shared_ptr<const Pose3d> source_pose, shared_ptr<const Pose3d> target_pose, shared_ptr<DynamicBodyd> body, MatrixNd& Jdot)
{
  assert(false);
  return Jdot;
}

/// Resets all force accumulators on the body.
void PseudoRigidBody::reset_accumulators()
{
  _core->reset_accumulators();
  for (unsigned i=0; i< _point_mass_states.size(); i++)
    _point_mass_states[i].f.set_zero();
}

/// Translates the entire body by vector @p o.
void PseudoRigidBody::translate(const Origin3d& o)
{
  _core->translate(o);

  // o is in the global frame.
  Point3d o_global(o, GLOBAL);

  for (unsigned i=0; i< _point_mass_states.size(); i++) {
    _point_mass_states[i].x +=
        Pose3d::transform_point(_point_mass_states[i].x.pose,
                                o_global);
  }
}

/// TODO: implement me.
void PseudoRigidBody::rotate(const Quatd& q)
{
  assert(false);
}

/// TODO: implement me.
double PseudoRigidBody::calc_kinetic_energy(boost::shared_ptr<const Pose3d> P)
{
  assert(false);
}

// The generalized coordinates pose for this body is the generalized coordinates
// pose for the rigid core.
boost::shared_ptr<const Pose3d> PseudoRigidBody::get_gc_pose() const
{
  return _core->get_gc_pose();
}

/// The computation frame for this body is the computation frame for the rigid
/// core. 
void PseudoRigidBody::set_computation_frame_type(ReferenceFrameType rftype)
{
  _core->set_computation_frame_type(rftype);
}

/// Gets the generalized inertia for this body.
SharedMatrixNd& PseudoRigidBody::get_generalized_inertia(SharedMatrixNd& M)
{
  const unsigned SPATIAL_D = 6;

  // Zero the matrix (it's mostly sparse).
  M.set_zero();

  // Set the rigid core.
  SharedMatrixNd coreJ = M.block(0, SPATIAL_D, 0, SPATIAL_D);
  _core->get_inertia().to_matrix(coreJ);

  // Set the components corresponding to all point masses.
  for (unsigned i=SPATIAL_D; i< M.rows(); i++)
    M(i,i) = get_compliant_layer_mass() / num_masses();

  return M;
}

// Solves the system M*X = B for unknown matrix X, where M is the generalized
// inertia and B is an arbitrary matrix.
SharedMatrixNd& PseudoRigidBody::solve_generalized_inertia(const SharedMatrixNd& B,
                                                           SharedMatrixNd& X)
{
  const unsigned THREE_D = 3, SPATIAL_D = 6;

  // Solve the rigid core first.
  shared_ptr<const Pose3d> P = get_gc_pose();
  for (unsigned i=0; i< B.columns(); i++)
  {
    SForced fcore(P);
    fcore.set_force(Vector3d(B.column(i).segment(0, THREE_D), P));
    fcore.set_torque(Vector3d(B.column(i).segment(THREE_D, SPATIAL_D), P));
    SAcceld acore = _core->get_inertia().inverse_mult(fcore);
    X.column(i).segment(0, THREE_D) = acore.get_linear();
    X.column(i).segment(THREE_D, SPATIAL_D) = acore.get_angular();

    // Solve all other masses.
    (X.column(i).segment(SPATIAL_D, X.rows()) =
        B.column(i).segment(SPATIAL_D, B.rows())) /=
            (get_compliant_layer_mass() / num_masses());
  }

  return X;
}

// Solves the system M*X = Bᵀ for unknown matrix X, where M is the generalized
// inertia and B is an arbitrary matrix.
SharedMatrixNd& PseudoRigidBody::transpose_solve_generalized_inertia(const SharedMatrixNd& B, SharedMatrixNd& X)
{
  const unsigned THREE_D = 3, SPATIAL_D = 6;

  shared_ptr<const Pose3d> P = get_gc_pose();
  for (unsigned i=0; i< B.rows(); i++)
  {
    // Solve the rigid core first.
    SForced fcore(P);
    fcore.set_force(Vector3d(B.row(i).segment(0, THREE_D), P));
    fcore.set_torque(Vector3d(B.row(i).segment(THREE_D, SPATIAL_D), P));
    SAcceld acore = _core->get_inertia().inverse_mult(fcore);
    X.column(i).segment(0, THREE_D) = acore.get_linear();
    X.column(i).segment(THREE_D, SPATIAL_D) = acore.get_angular();

    // Solve all other masses.
    (X.column(i).segment(SPATIAL_D, X.rows()) =
        B.row(i).segment(SPATIAL_D, B.columns())) /=
            (get_compliant_layer_mass() / num_masses());
  }

  return X;
}

// Solves the system M*x = b for unknown vector x, where M is the generalized
// inertia and b is an arbitrary vector.
SharedVectorNd& PseudoRigidBody::solve_generalized_inertia(const SharedVectorNd& b,
                                                           SharedVectorNd& x)
{
  const unsigned THREE_D = 3, SPATIAL_D = 6;

  // Solve the rigid core first.
  SForced fcore(get_gc_pose());
  fcore.set_force(Vector3d(b.segment(0, THREE_D), get_gc_pose()));
  fcore.set_torque(Vector3d(b.segment(THREE_D, SPATIAL_D), get_gc_pose()));
  SAcceld acore = _core->get_inertia().inverse_mult(fcore);
  x.segment(0, THREE_D) = acore.get_linear();
  x.segment(THREE_D, SPATIAL_D) = acore.get_angular();

  // Solve all other masses.
  (x.segment(SPATIAL_D, x.size()) = b.segment(SPATIAL_D, x.size())) /=
      (get_compliant_layer_mass() / num_masses());

  return x;
}

/// TODO: implement me.
SharedVectorNd& PseudoRigidBody::convert_to_generalized_force(shared_ptr<SingleBodyd> body, const SForced& w, SharedVectorNd& gf)
{
  assert(false);
  return gf;
}

// Gets the generalized coordinates for this body.
SharedVectorNd& PseudoRigidBody::get_generalized_coordinates_euler(SharedVectorNd& gc)
{
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eEuler);
  const unsigned THREE_D = 3;

  assert(gc.size() == num_masses()*3 + N_CORE_GC);
  SharedVectorNd gc_sub = gc.segment(0, N_CORE_GC);
  _core->get_generalized_coordinates_euler(gc_sub);
  for (unsigned i=0; i< num_masses(); i++)
    gc.segment(N_CORE_GC + i*THREE_D, N_CORE_GC + (i+1)*THREE_D) = _point_mass_states[i].x;

  return gc;
}

/// Gets the number of generalized coordinates for this body.
unsigned PseudoRigidBody::num_generalized_coordinates(GeneralizedCoordinateType gctype) const
{
  return _core->num_generalized_coordinates(gctype) + num_masses()*3;
}

void PseudoRigidBody::prepare_to_calc_ode(SharedConstVectorNd& x, double t, double dt, void* data)
{
  assert(false);
}

void PseudoRigidBody::prepare_to_calc_ode_sustained_constraints(SharedConstVectorNd& x, double t, double dt, void* data)
{
  assert(false); 
}

// TODO: implement me
void PseudoRigidBody::ode(double t, double dt, void* data, SharedVectorNd& dx)
{
  assert(false);
}

/// Sets the generalized coordinates for this body.
void PseudoRigidBody::set_generalized_coordinates_euler(const SharedVectorNd& gc)
{
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eEuler);
  const unsigned THREE_D = 3;

  assert(gc.size() == num_masses()*3 + N_CORE_GC);
  _core->set_generalized_coordinates_euler(gc.segment(0, N_CORE_GC));
  for (unsigned i=0; i< num_masses(); i++)
    _point_mass_states[i].x = gc.segment(N_CORE_GC + i*THREE_D, N_CORE_GC + (i+1)*THREE_D);
}

// Gets the generalized velocity for this body.
SharedVectorNd& PseudoRigidBody::get_generalized_velocity(GeneralizedCoordinateType gctype, SharedVectorNd& gv)
{
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(gctype);
  const unsigned THREE_D = 3;

  assert(gv.size() == num_masses()*3 + N_CORE_GC);
  SharedVectorNd gv_sub = gv.segment(0, N_CORE_GC);
  _core->get_generalized_velocity(gctype, gv_sub);
  for (unsigned i=0; i< num_masses(); i++)
    gv.segment(N_CORE_GC + i*THREE_D, N_CORE_GC + (i+1)*THREE_D) = _point_mass_states[i].xdot;

  return gv;
}

// Sets the generalized velocity for this body.
void PseudoRigidBody::set_generalized_velocity(GeneralizedCoordinateType gctype,
                                               const SharedVectorNd& gv)
{
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(gctype);
  const unsigned THREE_D = 3;

  assert(gv.size() == num_masses()*3 + N_CORE_GC);
  _core->set_generalized_velocity(gctype, gv.segment(0, N_CORE_GC));
  for (unsigned i=0; i< num_masses(); i++)
    _point_mass_states[i].xdot = gv.segment(N_CORE_GC + i*THREE_D, N_CORE_GC + (i+1)*THREE_D);
}

/// Gets the generalized forces on the body.
Ravelin::SharedVectorNd& PseudoRigidBody::get_generalized_forces(SharedVectorNd& gf)
{
  const unsigned THREE_D = 3;
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eSpatial);
  SharedVectorNd gf_sub = gf.segment(0, N_CORE_GC);
  _core->get_generalized_forces(gf_sub);
  for (unsigned i=0, j=N_CORE_GC; i< num_masses(); i++, j+= THREE_D)
    gf.segment(j, j+THREE_D) = _point_mass_states[i].f;

  return gf;
}

/// Adds the generalized force to the body.
void PseudoRigidBody::add_generalized_force(const SharedVectorNd& gf)
{
  const unsigned THREE_D = 3;
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eSpatial);
  _core->add_generalized_force(gf.segment(0, N_CORE_GC));
  for (unsigned i=0, j=N_CORE_GC; i< num_masses(); i++, j+= THREE_D)
    _point_mass_states[i].f += gf.segment(j, j+THREE_D);
}

/// Sets the generalized forces on the body.
void PseudoRigidBody::set_generalized_forces(const SharedVectorNd& gf)
{
  const unsigned THREE_D = 3;
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eSpatial);
  _core->set_generalized_forces(gf.segment(0, N_CORE_GC));
  for (unsigned i=0, j=N_CORE_GC; i< num_masses(); i++, j+= THREE_D)
    _point_mass_states[i].f = gf.segment(j, j+THREE_D);
}

/// Applies the generalized impulse to the body.
void PseudoRigidBody::apply_generalized_impulse(const SharedVectorNd& gj)
{
  const double point_mass_mass = get_compliant_layer_mass() / num_masses();
  const unsigned THREE_D = 3;
  const unsigned N_CORE_GC = _core->num_generalized_coordinates(DynamicBodyd::eSpatial);
  _core->apply_generalized_impulse(gj.segment(0, N_CORE_GC));
  for (unsigned i=0, j=N_CORE_GC; i< num_masses(); i++, j+= THREE_D)
  {
    Origin3d gj_sub = gj.segment(j, j+THREE_D);
    gj_sub /= point_mass_mass;
    _point_mass_states[i].xdot += Vector3d(gj_sub, GLOBAL);
  } 
}

/// Gets the generalized acceleration.
SharedVectorNd& PseudoRigidBody::get_generalized_acceleration(SharedVectorNd& ga)
{
  const unsigned N_SPATIAL = 6;
  const unsigned THREE_D = 3;

  assert(ga.size() == num_masses()*THREE_D + N_SPATIAL);
  SharedVectorNd ga_sub = ga.segment(0, N_SPATIAL);
  _core->get_generalized_acceleration(ga_sub);
  for (unsigned i=0; i< num_masses(); i++)
    ga.segment(N_SPATIAL + i*THREE_D, N_SPATIAL + (i+1)*THREE_D) = _point_mass_states[i].xddot;

  return ga;
}

/// Sets the generalized acceleration. 
void PseudoRigidBody::set_generalized_acceleration(const SharedVectorNd& ga)
{
  const unsigned N_SPATIAL = 6;
  const unsigned THREE_D = 3;

  assert(ga.size() == num_masses()*THREE_D + N_SPATIAL);
  SharedVectorNd* ga_nonconst = (SharedVectorNd*) &ga;
  SharedVectorNd ga_sub = ga_nonconst->segment(0, N_SPATIAL);
  _core->set_generalized_acceleration(ga_sub);
  for (unsigned i=0; i< num_masses(); i++)
    _point_mass_states[i].xddot = ga.segment(N_SPATIAL + i*THREE_D, N_SPATIAL + (i+1)*THREE_D); 
}


// Evaluates all constraints.
void PseudoRigidBody::evaluate_constraints(VectorNd& C)
{
  assert(C.size() == _constraints.size()*2);

  // Get the frame for the rigid body (pistons will be in this frame too).
  shared_ptr<const Pose3d> P = _core->get_pose();
 
  for (unsigned i=0, j=0; i< _constraints.size(); i++)
  {
    // Get the point mass position in the proper frame
    Point3d x = Pose3d::transform_point(P, _point_mass_states[i].x);

    // Compute r = x - (core_x + R*u)
    Vector3d r = x - _constraints[i].p;

    // Evaluate.
    C[j++] = r.dot(_constraints[i].v1);
    C[j++] = r.dot(_constraints[i].v2);
  }
}

// Evaluates the time derivative of all constraints.
void PseudoRigidBody::evaluate_constraints_dot(VectorNd& Cdot)
{
  assert(Cdot.size() == _constraints.size()*2);

  // Get the frame for the rigid body (pistons will be in this frame too).
  shared_ptr<const Pose3d> P = _core->get_pose();

  for (unsigned i=0, j=0; i< _constraints.size(); i++)
  {
    // Compute dr/dt = xdot - (core_xd + w x r)
    const Vector3d r = Pose3d::transform_point(GLOBAL, _constraints[i].p) -
                       Pose3d::calc_relative_pose(P, GLOBAL).x;
    const Vector3d xd0 = Pose3d::transform_vector(GLOBAL, 
                         _core->get_velocity().get_linear());
    const Vector3d w0 = Pose3d::transform_vector(GLOBAL,
                        _core->get_velocity().get_angular());
    const Vector3d rdot = Pose3d::transform_vector(P, 
        _point_mass_states[i].xdot - xd0 + Vector3d::cross(w0, r));

    // Evaluate.
    Cdot[j++] = rdot.dot(_constraints[i].v1);
    Cdot[j++] = rdot.dot(_constraints[i].v2);
  }
}

// Evaluates the second time derivative of all constraints.
void PseudoRigidBody::evaluate_constraints_ddot(VectorNd& Cddot)
{
  assert(Cddot.size() == _constraints.size()*2);

  // Get the frame for the rigid body (pistons will be in this frame too).
  shared_ptr<const Pose3d> P = _core->get_pose();

  for (unsigned i=0, j=0; i< _constraints.size(); i++)
  {
    // Compute d^2r/dt^2 = xddot - (core_xdd + a x r + w x w x r)
    const Vector3d r = Pose3d::transform_point(GLOBAL, _constraints[i].p) -
                       Pose3d::calc_relative_pose(P, GLOBAL).x;
    const Vector3d xdd0 = Pose3d::transform_vector(GLOBAL, 
                          _core->get_accel().get_linear());
    const Vector3d w0 = Pose3d::transform_vector(GLOBAL,
                        _core->get_velocity().get_angular());
    const Vector3d a0 = Pose3d::transform_vector(GLOBAL,
                        _core->get_accel().get_angular());
    const Vector3d rddot = Pose3d::transform_vector(P, 
        _point_mass_states[i].xddot - xdd0 + Vector3d::cross(a0, r) +
        Vector3d::cross(w0, Vector3d::cross(w0, r)));

    // Evaluate.
    Cddot[j++] = rddot.dot(_constraints[i].v1);
    Cddot[j++] = rddot.dot(_constraints[i].v2);
  }
}

// Computes the Jacobian matrix for the constraints.
// Constraint equation is C(q) = 0.
// dC(q)/dt = xC/xq \cdot \dot{q}
// So, the Jacobian matrix can be computed by setting the generalized velocity
// each of the unit vectors in R^n, evaluating dC(q)/dt (which is usually
// simple to evaluate), and populating an output matrix CC, which happens to
// be equal to the Jacobian matrix.
void PseudoRigidBody::calc_jacobian(MatrixNd& J)
{
  // Save the current generalized velocities.
  VectorNd gv_save;
  DynamicBodyd::get_generalized_velocity(DynamicBodyd::eSpatial, gv_save);
  
  // Resize the matrix J appropriately.
  J.set_zero(num_constraint_eqns(), gv_save.size());

  // Probe every velocity.
  VectorNd dq_probe(gv_save.size());
  for (unsigned i=0; i < gv_save.size(); i++)
  {
    // Set the i'th basis vector.
    dq_probe.set_zero(gv_save.size());
    dq_probe[i] = 1;

    // Set the generalized velocity.
    DynamicBodyd::set_generalized_velocity(DynamicBodyd::eSpatial, dq_probe);

    // Evaluate the constraints.
    dq_probe.resize(num_constraint_eqns());
    evaluate_constraints_dot(dq_probe);
    J.column(i) = dq_probe;
  }

  // Restore the generalized velocities.
  DynamicBodyd::set_generalized_velocity(DynamicBodyd::eSpatial, gv_save);
}

// Computes the time derivative of the Jacobian matrix of the constraint
// equations times the generalized velocity.
// Constraint equation is C(q) = 0.
// d^2C(q)/dt^2 = xC/xq \cdot \ddot{q} + d/dt(xC/xq) \cdot \dot{q} 
// So, the vector can be computed by setting the generalized
// acceleration to zero and evaluating d^2C(q)/dt^2 (which is usually
// simple to evaluate).
void PseudoRigidBody::calc_jacobian_dot_v(VectorNd& Jdot_v)
{
  // Save the current generalized acceleration
  VectorNd ga_save;
  DynamicBodyd::get_generalized_acceleration(ga_save);

  // Make the generalized acceleration zero.
  VectorNd zero(ga_save.size());
  zero.set_zero();
  DynamicBodyd::set_generalized_acceleration(zero);
  
  // Resize Jdot_v appropriately.
  Jdot_v.resize(num_constraint_eqns());

  // Evaluate the constraints.
  evaluate_constraints_ddot(Jdot_v);

  // Restore the generalized acceleration.
  DynamicBodyd::set_generalized_acceleration(ga_save);
}

// Computes forward dynamics for this body.
void PseudoRigidBody::calc_fwd_dyn()
{
  MatrixNd J, iM_JT, J_iM_JT;
  VectorNd f, iM_f, J_iM_f, Jdot_v;
  VectorNd v, lambda, vdot, C, Cdot;
  const double ALPHA = 1000, BETA = 100;

  // Evaluate the constraints and their time derivatives.
  C.resize(num_constraint_eqns());
  Cdot.resize(num_constraint_eqns());
  evaluate_constraints(C);
  evaluate_constraints_dot(Cdot);

  // Get the generalized velocity.
  DynamicBodyd::get_generalized_velocity(DynamicBodyd::eSpatial, v);

  // Get the generalized force vector.
  DynamicBodyd::get_generalized_forces(f);

  // Update f with spring and damper forces
  update_generalized_forces_with_spring_damper_forces(f);

  // Compute the constraint Jacobian and its time derivative.
  calc_jacobian(J);
  calc_jacobian_dot_v(Jdot_v);

  // Compute the constraint forces using:
  // | M  -J' | | dv/dt | = | f        |
  // | J   0  | | λ     |   | -dJ/dt⋅v |
  // dv/dt = M⁻¹(f + J'λ)
  // JM⁻¹J'λ = -JM⁻¹f - dJ/dt⋅xd
  DynamicBodyd::transpose_solve_generalized_inertia(J, iM_JT);
  J.mult(iM_JT, J_iM_JT);
  DynamicBodyd::solve_generalized_inertia(f, iM_f);
  J.mult(iM_f, J_iM_f);
std::cout << "C: " << C << std::endl;
std::cout << "dC/dt: " << Cdot << std::endl;
std::cout << "J: " << std::endl << J;
std::cout << "v: " << v << std::endl;
std::cout << "f: " << f << std::endl;
std::cout << "inv(M)*f: " << iM_f << std::endl;
std::cout << "J * inv(M)*f: " << J_iM_f << std::endl;
std::cout << "dot(J) * v: " << Jdot_v << std::endl;
  C *= ALPHA;
  Cdot *= BETA;
  ((lambda = J_iM_f) += Jdot_v += C += Cdot).negate();
std::cout << "J*inv(M)*J': " << std::endl << J_iM_JT;
std::cout << "rhs: " << lambda << std::endl;
  _LA.solve_LS_fast(J_iM_JT, lambda, LinAlgd::eSVD1, 1e-12);
std::cout << "lambda: " << lambda << std::endl;

  // Solve for the joint and floating base accelerations using:
  // dv/dt = M⁻¹(J'λ + fext)
  iM_JT.mult(lambda, vdot) += iM_f;
std::cout << "vdot: " << vdot << std::endl;

  // Set the acceleration
  SAcceld a(get_gc_pose());
  Vector3d xdd(vdot.segment(0,3), get_gc_pose());
  Vector3d alpha(vdot.segment(3,6), get_gc_pose());
  a.set_linear(xdd);
  a.set_angular(alpha);
  _core->set_accel(a);
  for (unsigned i=0, j=6; i< _point_mass_states.size(); i++, j+= 3)
    _point_mass_states[i].xddot = vdot.segment(j,j+3);

  // Evaluate the constraint forces
VectorNd Cddot(num_constraint_eqns());
evaluate_constraints_ddot(Cddot);
std::cout << "d^2C/dt^2: " << Cddot << std::endl;
}

// Determines the extension of a spring.
double PseudoRigidBody::determine_spring_extension(unsigned i) const
{
  // See whether the spring is attached to a piston.
  if (_springs[i].constraint_type == SpringDamperData::ConstraintType::ePiston)
  {
    // Get the piston index.
    const unsigned j = _springs[i].piston_index;

    // Get the vector emanating along the piston direction on the rigid body,
    // which will define a plane (using the piston emanation point on the rigid
    // body as well). By projecting the point mass onto this plane, the signed
    // distance, meaning the "extension" can be computed. Plane will correspond
    // to <u,x> = <u, px> = d 
    const Vector3d u0 = Pose3d::transform_vector(GLOBAL, _constraints[j].u);
    const Point3d p0 = Pose3d::transform_point(GLOBAL, _constraints[j].p);
    const double d = u0.dot(p0);

    // Get x
    const Point3d& x = _point_mass_states[_constraints[j].point_mass_index].x;

    // Compute the projection.
    return u0.dot(x) - d;
  }
  else
  {
    // Get the locations of the two point masses and the distance between
    // them.
    const unsigned idx1 = _springs[i].pm1_index;
    const unsigned idx2 = _springs[i].pm2_index;
    const Point3d& loc1 = _point_mass_states[idx1].x;
    const Point3d& loc2 = _point_mass_states[idx2].x;
    return (loc1 - loc2).norm();
  }
}

// Determines the time derivative of the extension of a spring.
double PseudoRigidBody::determine_spring_extension_dot(unsigned i) const
{
  // See whether the spring is attached to a piston.
  if (_springs[i].constraint_type == SpringDamperData::ConstraintType::ePiston)
  {
    // This is the time derivative of the routine in
    // determine_spring_extension() 
    // 
    // Get the piston index.
    const unsigned j = _springs[i].piston_index;

    // u should be defined in the rigid core frame. u0 (u defined in the global
    // frame) is R*u. du0/dt = w0 x u0.
    assert(_constraints[j].u.pose == _core->get_pose());
    const Vector3d u0 = Pose3d::transform_vector(GLOBAL, _constraints[j].u);
    const Vector3d w0 = Pose3d::transform(GLOBAL, _core->get_velocity()).get_angular(); 
    const Vector3d du0dt = Vector3d::cross(w0, u0);

    // Get x and xdot
    const Vector3d& x = _point_mass_states[_constraints[j].point_mass_index].x;
    const Vector3d& xdot = _point_mass_states[_constraints[j].point_mass_index].xdot;

    // Original operation was u0.dot(x) = u0.dot(p0). Time derivative of this 
    // is du0/dt.dot(x) + u0.dot(xdot) = du0/dt.dot(p0) + 
    //                                   u0.dot(dx'/dt + w0 x r),
    // where x' is the center of mass of the rigid core, and r is the vector
    // from the center of mass of the rigid core to p in the global frame.
    const Point3d p0 = Pose3d::transform_point(GLOBAL, _constraints[j].p);
    const Vector3d r = p0 - Pose3d::calc_relative_pose(_core->get_pose(), GLOBAL).x;
    const Vector3d xprime_dot = Pose3d::transform(GLOBAL, _core->get_velocity()).get_linear(); 
    return du0dt.dot(x) + u0.dot(xdot) - du0dt.dot(p0) - u0.dot(xprime_dot + Vector3d::cross(w0, r));
  }
  else
  {
    // Get the locations of the two point masses and the distance between
    // them.
    const unsigned idx1 = _springs[i].pm1_index;
    const unsigned idx2 = _springs[i].pm2_index;
    const Point3d& loc1 = _point_mass_states[idx1].x;
    const Point3d& loc2 = _point_mass_states[idx2].x;

    // Get the velocities of the two point masses and the time derivative of
    // the distance between them.
    const Vector3d& v1 = _point_mass_states[idx1].xdot;
    const Vector3d& v2 = _point_mass_states[idx2].xdot;
    return 1.0/(loc1 - loc2).norm() * 2 * (loc1 - loc2).dot(v1 - v2);
  }
}

void PseudoRigidBody::convert_spring_force_to_generalized_force(unsigned i, double f, VectorNd& gf) const
{
  // See whether the spring is attached to a piston.
  if (_springs[i].constraint_type == SpringDamperData::ConstraintType::ePiston)
  {
    transform_piston_force_to_generalized_force(i, f, gf);
    return;
  }
  else
  {
    const unsigned THREE_D = 3, SPATIAL_D = 6;

    // Get the vector between the two point masses, which will end up becoming
    // the force vector.
    unsigned idx1 = _springs[i].pm1_index;
    unsigned idx2 = _springs[i].pm2_index;
    const Point3d& loc1 = _point_mass_states[idx1].x;
    const Point3d& loc2 = _point_mass_states[idx2].x;
    const Vector3d v = loc2 - loc1;
    double v_nrm = v.norm();
    assert(v_nrm > NEAR_ZERO);
    const Vector3d vf = v * (f / v_nrm);

    // Apply the force vector to the two point masses.
    gf.segment(SPATIAL_D + idx1*THREE_D, SPATIAL_D + (idx1+1)*THREE_D) = Origin3d(vf);
    gf.segment(SPATIAL_D + idx2*THREE_D, SPATIAL_D + (idx2+1)*THREE_D) = Origin3d(-vf);
  } 
}

// Updates the generalized forces on the rigid body and the mass elements due
// to the springs and dampers.
void PseudoRigidBody::update_generalized_forces_with_spring_damper_forces(VectorNd& gf)
{
  // Use a temporary generalized force vector.
  VectorNd tmp_gf = gf;

  // Loop over all springs 
  for (unsigned i=0; i< _springs.size(); i++)
  {
    // determine the extension/compression of the spring.
    double x = determine_spring_extension(i);

    // determine the time derivative of the spring extension/compression.
    double xdot = determine_spring_extension_dot(i);

    // compute the spring / damping force.
    const double f = _spring_k * (_springs[i].rest_len - x) +
                     _spring_c * -xdot;

    // convert this to a generalized force.
    convert_spring_force_to_generalized_force(i, f, tmp_gf);

    // Sum with the existing generalized forces.
    gf += tmp_gf;
  }

// TODO: delete!
/*
  // Loop over all point mass connections
  for (unsigned i=0; i< _pairwise_point_mass_springs.size(); i++)
  {
    // Get the locations of the two point masses and the distance between
    // them.
    unsigned idx1 = _pairwise_point_mass_springs[i].pm1_index;
    unsigned idx2 = _pairwise_point_mass_springs[i].pm2_index;
    const Point3d& loc1 = _point_mass_states[idx1].x;
    const Point3d& loc2 = _point_mass_states[idx2].x;
    const double dist = (loc1 - loc2).norm();

    // Get the velocities of the two point masses and the time derivative of
    // the distance between them.
    const Vector3d& v1 = _point_mass_states[idx1].xdot;
    const Vector3d& v2 = _point_mass_states[idx2].xdot;
    const double dot_dist = 1.0/(loc1 - loc2).norm() * (v1 - v2);

    // Get the spring index.
    const unsigned sidx = _pairwise_point_mass_springs[i].spring_index;

    // Compute the force to be applied to the bodies.
    const double f = _spring_k * (_spring_resting_lengths[sidx])
  }

  // loop over all springs
  for (unsigned i=0; i< _spring_resting_lengths.size(); i++) {
    // compute the spring and damping force.
    const double
        f = _spring_k * (_spring_resting_lengths[i] - _joint_positions[i]) +
        _spring_c * (-_joint_velocities[i]);

    // TODO: apply the force to the joint.
  }
*/
}

void PseudoRigidBody::update_visualization()
{
  const unsigned X = 0, Y = 1, Z = 2;
  assert(_vert_array);

  // Update the vertex array points using the location of vertices at their
  // deformed locations. 
  for (unsigned i=0; i< _point_mass_states.size(); i++)
  {
    (*_vert_array)[i] = osg::Vec3((float) _point_mass_states[i].x[X],
                                  (float) _point_mass_states[i].x[Y], 
                                  (float) _point_mass_states[i].x[Z]);
  }
}

