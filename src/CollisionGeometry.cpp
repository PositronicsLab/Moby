/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <stdexcept>
#include <iostream>
#include <stack>
#include <fstream>
#include <Moby/Polyhedron.h>
#include <Moby/CompGeom.h>
#include <Moby/RigidBody.h>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/CollisionGeometry.h>

using std::vector;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Constructs a CollisionGeometry with no triangle mesh, identity transformation and relative transformation
CollisionGeometry::CollisionGeometry()
{
  _F = shared_ptr<Pose3d>(new Pose3d);
}

/// Gets a supporting point for this geometry in a particular direction
Point3d CollisionGeometry::get_supporting_point(const Vector3d& d) const
{
  // get the primitive from this
  PrimitivePtr primitive = get_geometry();

  // get this
  shared_ptr<const CollisionGeometry> cg_const = dynamic_pointer_cast<const CollisionGeometry>(shared_from_this());
  CollisionGeometryPtr cg = const_pointer_cast<CollisionGeometry>(cg_const);

  // get the pose for this
  shared_ptr<const Pose3d> P = primitive->get_pose(cg);

  // transform the vector
  Vector3d dir = Pose3d::transform_vector(P, d);

  // get the supporting point from the primitive
  return primitive->get_supporting_point(dir);
}

/// Gets the farthest point from this geometry
double CollisionGeometry::get_farthest_point_distance() const
{
  // get the primitive from this
  PrimitivePtr primitive = get_geometry();

  // get the vertices
  vector<Point3d> verts;
  get_vertices(verts);
  if (verts.empty())
    return 0.0;

  // get the rigid body pose in P's frame
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(get_single_body());
  Point3d rbX = Pose3d::transform_point(verts.front().pose, Point3d(0,0,0,rb->get_pose()));
 
  // find which point is closest
  double max_dist = 0.0;
  for (unsigned i=0; i< verts.size(); i++)
    max_dist = std::max(max_dist, (verts[i] - rbX).norm()); 

  return max_dist; 
}

/// Sets the single body associated with this CollisionGeometry
void CollisionGeometry::set_single_body(SingleBodyPtr s)
{
  _single_body = s;
  _F->rpose = s->get_pose();
}

/// Sets the collision geometry via a primitive
/**
 * The primitive is not cloned, nor is it unaltered; <b>this</b> points to
 * the pointer returned by this method (typically <b>primitive</b>).
 * \return primitive <i>unless the geometry of the underlying primitive is 
 * inconsistent, degenerate, or non-convex<i>; in that case,  
 * a corrected primitive will be returned.
 */
PrimitivePtr CollisionGeometry::set_geometry(PrimitivePtr primitive)
{
  Quatd EYE;

  if (_single_body.expired())
    throw std::runtime_error("CollisionGeometry::set_geometry() called before single body set!");

  SingleBodyPtr sb(_single_body);
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(sb);
  if (rb && !Quatd::rel_equal(rb->get_pose()->q, EYE))
  {
    std::cerr << "CollisionGeometry::set_primitive() warning - rigid body's orientation is not identity." << std::endl;
    std::cerr << "  At the rigid body's current orientation (" << AAngled(rb->get_pose()->q) << ")" << std::endl;
    std::cerr << "  the primitive wll have the orientation (" << AAngled(primitive->get_pose()->q) << ")" << std::endl;
  }

  // save the primitive
  _geometry = primitive;

  // add this to the primitive
  CollisionGeometryPtr cg = dynamic_pointer_cast<CollisionGeometry>(shared_from_this());
  primitive->add_collision_geometry(cg);

  return primitive;
}

/// Gets vertices for a primitive
void CollisionGeometry::get_vertices(std::vector<Point3d>& vertices) const
{
  // get the primitive from this
  PrimitivePtr primitive = get_geometry();

  // get a pointer to this
  shared_ptr<const CollisionGeometry> cg_const = dynamic_pointer_cast<const CollisionGeometry>(shared_from_this());
  CollisionGeometryPtr cg = const_pointer_cast<CollisionGeometry>(cg_const);

  // get the pose for this
  shared_ptr<const Pose3d> P = primitive->get_pose(cg);
 
  // get the vertices from the primitive
  vertices.clear();
  primitive->get_vertices(P, vertices);
}

/// Writes the collision geometry mesh to the specified VRML file
/**
 * \note the mesh is transformed using the current transformation
 */
void CollisionGeometry::write_vrml(const std::string& filename) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // open the file for writing
  std::ofstream out(filename.c_str());
  assert(!out.fail());

  // write the VRML header
  out << "#VRML V2.0 utf8" << std::endl << std::endl;

  // write the transform
  Origin3d P = _F->x;
  AAngled aa = _F->q;
  out << "Transform { " << std::endl;
  out << "  rotation " << aa.x << " " << aa.y << " " << aa.z << " " << aa.angle << std::endl;
  out << "  translation " << P[X] << " " << P[Y] << " " << P[Z] << std::endl;
  out << "  children [" << std::endl;
 
  // get the vertices according to this pose
  shared_ptr<const CollisionGeometry> cg_const = dynamic_pointer_cast<const CollisionGeometry>(shared_from_this());
  CollisionGeometryPtr cg = const_pointer_cast<CollisionGeometry>(cg_const);
  shared_ptr<const Pose3d> pose = _geometry->get_pose(cg);
 
  // get the geometry
  const std::vector<Origin3d>& vertices = _geometry->get_mesh(pose)->get_vertices();
  const std::vector<IndexedTri>& facets = _geometry->get_mesh(pose)->get_facets();

  // write the mesh
  out << "   Shape {" << std::endl;
  out << "      appearance Appearance { material Material {" << std::endl;
  out << "        transparency 0" << std::endl;
  out << "        shininess 0.2" << std::endl;
  out << "        ambientIntensity 0.2" << std::endl;
  out << "        emissiveColor 0 0 0" << std::endl;
  out << "        specularColor 0 0 0" << std::endl;
  out << "        diffuseColor .8 .8 .8" << std::endl;
  out << "        }}" << std::endl;

  // write the geometry
  out << "      geometry IndexedFaceSet {" << std::endl; 
  out << "        coord Coordinate { point [ ";
  for (unsigned i=0; i< vertices.size(); i++)
     out << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << ", "; 
  out << " ] }" << std::endl;
  out << "        coordIndex [ ";
  for (unsigned i=0; i< facets.size(); i++)
    out << facets[i].a << " " << facets[i].b << " " << facets[i].c << " -1, ";

  out << " ] } } ] }" << std::endl;
  
  // close the file
  out.close();
}

/// Sets the relative pose of this geometry.
/**
 * \param P the relative pose (P.pose must be set relative to the single body pose)
 */
void CollisionGeometry::set_relative_pose(const Pose3d& P)
{
  // verify the P's relative pose is correct 
  if (P.rpose != _F->rpose)
    throw std::runtime_error("CollisionGeometry::set_relative_pose() - pose not relative to the underlying single body");

  // set the pose
  *_F = P;
}

/// Calculates the (unsigned) distance of a point from this collision geometry
double CollisionGeometry::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  // get the primitive from this
  PrimitivePtr primitive = get_geometry();

  // get the collision geometry
  shared_ptr<const CollisionGeometry> cg_const = dynamic_pointer_cast<const CollisionGeometry>(shared_from_this());
  CollisionGeometryPtr cg = const_pointer_cast<CollisionGeometry>(cg_const);

  // setup a new point with a new pose
  Point3d px = Pose3d::transform_point(primitive->get_pose(cg), p);

  // call the primitive function
  return primitive->calc_dist_and_normal(px, normals);
}

/// Calculates the signed distance for a primitive
double CollisionGeometry::calc_signed_dist(const Point3d& p)
{
  // get the primitive from this
  PrimitivePtr primitive = get_geometry();

  // get the collision geometry
  CollisionGeometryPtr cg = dynamic_pointer_cast<CollisionGeometry>(shared_from_this());

  // setup a new point with a new pose
  Point3d px = Pose3d::transform_point(primitive->get_pose(cg), p);

  // call the primitive function
  return primitive->calc_signed_dist(px);
}

/// Calculates the signed distances between two geometries and returns closest points if geometries are not interpenetrating
double CollisionGeometry::calc_signed_dist(CollisionGeometryPtr gA, CollisionGeometryPtr gB, Point3d& pA, Point3d& pB) 
{
  // get the two primitives
  PrimitivePtr primA = gA->get_geometry();
  PrimitivePtr primB = gB->get_geometry();

  // setup poses for the points
  pA.pose = primA->get_pose(gA);
  pB.pose = primB->get_pose(gB);

  FILE_LOG(LOG_COLDET) << "CollisionGeometry::calc_signed_dist() - computing signed distance between " << gA->get_single_body()->id << " and " << gB->get_single_body()->id << std::endl;

  // now compute the signed distance
  return primA->calc_signed_dist(primB, pA, pB);
}

/*
/// Calculates the closest points (and squared distance) between geometries a and b
 * \param a the first collision geometry
 * \param b the second collision geometry
 * \param aTb the transform from b's frame to a's frame
 * \param cpa the closest point to b on a (in a's frame)
 * \param cpb the closest point to a on b (in b's frame)
 * \return the squared distance between cpa and cpb
double CollisionGeometry::calc_dist(CollisionGeometryPtr a, CollisionGeometryPtr b, const Transform3d& aTb, Point3d& cpa, Point3d& cpb)
{
  // get the poses of a and b
  shared_ptr<const Pose3d> Pa = a->get_pose();
  shared_ptr<const Pose3d> Pb = b->get_pose();

  // setup the minimum distance
  double min_dist = std::numeric_limits<double>::max();

  // determine the transform from b's frame to a's frame
  Transform3d bTa = Transform3d::invert(aTb);

  // get the primitives 
  PrimitivePtr a_primitive = a->get_geometry(); 
  PrimitivePtr b_primitive = b->get_geometry();

  // get the mesh data from the plugins
  const IndexedTriArray& a_mesh = *a_primitive->get_mesh();
  const IndexedTriArray& b_mesh = *b_primitive->get_mesh(); 

  // do pairwise distance checks
  for (unsigned i=0; i< a_mesh.num_tris(); i++)
  {
    // get the triangle
    Triangle ta = a_mesh.get_triangle(i, Pa);

    // loop over all triangles in b
    for (unsigned j=0; j< b_mesh.num_tris(); j++)
    {
      // get the untransformed second triangle
      Triangle utb = b_mesh.get_triangle(j, Pb);

      // transform the second triangle
      Triangle tb = Triangle::transform(utb, aTb);
 
      // get distance between the two triangles
      Point3d cpa_tmp, cpb_tmp;
      double dist = Triangle::calc_sq_dist(ta, tb, cpa_tmp, cpb_tmp);

      // if it's the minimum distance, save the closest points
      if (dist < min_dist)
      {
        min_dist = dist;
        cpa = cpa_tmp;
        cpb = bTa.transform_point(cpb_tmp);
      }
    }
  }

  return min_dist;
}
*/

/// Implements Base::load_from_xml()
void CollisionGeometry::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // load Base specific data
  Base::load_from_xml(node, id_map);

  // verify that this node is of type CollisionGeometry
  assert (strcasecmp(node->name.c_str(), "CollisionGeometry") == 0);

  // read relative pose, if specified
  Pose3d TR;
  TR.rpose = _F->rpose;
  XMLAttrib* rel_origin_attr = node->get_attrib("relative-origin");
  XMLAttrib* rel_rpy_attr = node->get_attrib("relative-rpy");
  XMLAttrib* rel_quat_attr = node->get_attrib("relative-quat");
  if (rel_origin_attr)
    TR.x = rel_origin_attr->get_origin_value();
  if (rel_quat_attr)
    TR.q = rel_quat_attr->get_quat_value();
  else if (rel_rpy_attr)
    TR.q = rel_rpy_attr->get_rpy_value();
  set_relative_pose(TR);

  // read the primitive ID, if any
  XMLAttrib* primitive_id_attrib = node->get_attrib("primitive-id");
  if (primitive_id_attrib)
  {
    // get the primitive ID
    const std::string& id = primitive_id_attrib->get_string_value();

    // search for it
    std::map<std::string, BasePtr>::const_iterator id_iter = id_map.find(id);
    if (id_iter == id_map.end())
    {
      std::cerr << "CollisionGeometry::load_from_xml() - primitive with ID '";
      std::cerr << id << "'" << std::endl << "  not found in offending node: ";
      std::cerr << std::endl << *node;
    }
    else
    {
      PrimitivePtr pc = boost::dynamic_pointer_cast<Primitive>(id_iter->second);
      PrimitivePtr newpc = set_geometry(pc);

      // now, we're going to reset the geometry in the ID map if it has changed
      if (newpc != pc)
        id_map[id] = newpc;
    }  
  }
}

/// Implements Base::save_to_xml()
void CollisionGeometry::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save Base data
  Base::save_to_xml(node, shared_objects);

  // set the node name
  node->name = "CollisionGeometry";

  // add the relative pose 
  node->attribs.insert(XMLAttrib("relative-origin", _F->x));
  node->attribs.insert(XMLAttrib("relative-quat", _F->q));

  // save the ID of the primitive and add the primitive to the shared list
  if (_geometry)
  {
    node->attribs.insert(XMLAttrib("primitive-id", _geometry->id));
    shared_objects.push_back(_geometry);
  }
}

