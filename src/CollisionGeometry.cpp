/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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

using boost::dynamic_pointer_cast;
using Ravelin::Pose3d;
using Ravelin::Origin3d;
using Ravelin::AAngled;
using Ravelin::Point3d;
using boost::shared_ptr;
using namespace Moby;

/// Constructs a CollisionGeometry with no triangle mesh, identity transformation and relative transformation
CollisionGeometry::CollisionGeometry()
{
  _F = shared_ptr<Pose3d>(new Pose3d);
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
  if (_single_body.expired())
    throw std::runtime_error("CollisionGeometry::set_geometry() called before single body set!");

  SingleBodyPtr sb(_single_body);
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(sb);
  if (rb && !Pose3d::rel_equal(*rb->get_pose(), *shared_ptr<const Pose3d>()))    
    std::cerr << "CollisionGeometry::set_primitive() warning - rigid body's transform is not identity!" << std::endl;

  // save the primitive
  _geometry = primitive;

  return primitive;
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
  
  // get the geometry
  const std::vector<Point3d>& vertices = _geometry->get_mesh()->get_vertices();
  const std::vector<IndexedTri>& facets = _geometry->get_mesh()->get_facets();

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

/// Sets the transform of this CollisionGeometry
/**
 * This method sets the base transform of this geometry. This geometry's true 
 * transform is the result of the relative transform applied to this transform.
 * This method recursively calls itself on each of its children, so that the 
 * user need only make one call to set_transform() at the root of a tree of 
 * transforms.  Note that the calling method must consider the relative 
 * transform applied to this geometry.  If the transform of the parent / 
 * related body is T1 and the relative transform is T2 then set_transform() 
 * should be called with T1 * T2.
 * \param transform the transformation
 * \param relF_accounted determines whether the relative transform is accounted for in transform; if 
 *        <b>false</b>, transform will be transformed by the relative transform before storing and propagating
 * \sa get_relF()
 * \sa set_relF() 
 */
void CollisionGeometry::set_relative_pose(shared_ptr<const Pose3d> P)
{
  // update the transform
  *_F = *P;
}

/// Implements Base::load_from_xml()
void CollisionGeometry::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // load Base specific data
  Base::load_from_xml(node, id_map);

  // verify that this node is of type CollisionGeometry
  assert (strcasecmp(node->name.c_str(), "CollisionGeometry") == 0);

  // read relative pose, if specified
  shared_ptr<Pose3d> TR(new Pose3d);
  const XMLAttrib* rel_origin_attr = node->get_attrib("relative-origin");
  const XMLAttrib* rel_rpy_attr = node->get_attrib("relative-rpy");
  const XMLAttrib* rel_quat_attr = node->get_attrib("relative-quat");
  if (rel_origin_attr)
    TR->x = rel_origin_attr->get_origin_value();
  if (rel_quat_attr)
    TR->q = rel_quat_attr->get_quat_value();
  else if (rel_rpy_attr)
    TR->q = rel_rpy_attr->get_rpy_value();
  set_relative_pose(TR);

  // read the primitive ID, if any
  const XMLAttrib* primitive_id_attrib = node->get_attrib("primitive-id");
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

