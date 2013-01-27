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
#include <Moby/AAngle.h>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/CollisionGeometry.h>

using namespace Moby;
using boost::dynamic_pointer_cast;

/// Constructs a CollisionGeometry with no triangle mesh, identity transformation and relative transformation
CollisionGeometry::CollisionGeometry()
{
  _transform = IDENTITY_4x4;
  _rel_transform = IDENTITY_4x4;
  _rel_transform_identity = true;
}

/// Sets the relative transform (from its parent CollisionGeometry or dynamic body) for this CollisionGeometry
void CollisionGeometry::set_rel_transform(const Matrix4& transform, bool update_global_transform)
{  
  // see whether to update the global transform
  if (update_global_transform)
  {
    // determine how the transform will change
    Matrix4 update = Matrix4::inverse_transform(_rel_transform) * transform;

    // update the global transform
    _transform = _transform * update;
  
    // update transforms for all of the children
    BOOST_FOREACH(CollisionGeometryPtr cg, _children)
      cg->set_transform(_transform, false);
  }

  // set the relative transform
  _rel_transform = transform;

  // check whether the relative transform is equal (to within floating point tolerance) of identity
  _rel_transform_identity = Matrix4::epsilon_equals(_rel_transform, IDENTITY_4x4, NEAR_ZERO);
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
  if (rb && !Matrix4::epsilon_equals(IDENTITY_4x4, rb->get_transform(), NEAR_ZERO))
    std::cerr << "CollisionGeometry::set_primitive() - rigid body's transform is not identity!" << std::endl;

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
  Vector3 P = _transform.get_translation();
  AAngle aa(&_transform);
  out << "Transform { " << std::endl;
  out << "  rotation " << aa.x << " " << aa.y << " " << aa.z << " " << aa.angle << std::endl;
  out << "  translation " << P[X] << " " << P[Y] << " " << P[Z] << std::endl;
  out << "  children [" << std::endl;
  
  // get the geometry
  const std::vector<Vector3>& vertices = _geometry->get_mesh()->get_vertices();
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
 * \param rel_transform_accounted determines whether the relative transform is accounted for in transform; if 
 *        <b>false</b>, transform will be transformed by the relative transform before storing and propagating
 * \sa get_rel_transform()
 * \sa set_rel_transform() 
 */
void CollisionGeometry::set_transform(const Matrix4& transform, bool rel_transform_accounted)
{
  // determine whether to account for the relative transform
  if (rel_transform_accounted || _rel_transform_identity)
    _transform = transform;
  else
    _transform = transform * _rel_transform;
  
  // set the transform for all of the children
  BOOST_FOREACH(CollisionGeometryPtr cg, _children)
    cg->set_transform(transform, false);
}

/// Implements Base::load_from_xml()
void CollisionGeometry::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // load Base specific data
  Base::load_from_xml(node, id_map);

  // verify that this node is of type CollisionGeometry
  assert (strcasecmp(node->name.c_str(), "CollisionGeometry") == 0);

  // read transform, if specified
  const XMLAttrib* transform_attrib = node->get_attrib("transform");
  if (transform_attrib)
  {
    Matrix4 T;
    transform_attrib->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
    {
      std::cerr << "CollisionGeometry::load_from_xml() warning: bad transform? ";
      std::cerr << std::endl << T << " when reading node " << std::endl;
      std::cerr << *node << std::endl;
      std::cerr << "  --> possibly a floating-point error..." << std::endl;
    }
    set_transform(T, true);
  }

  // read relative transform, if specified
  const XMLAttrib* rel_transform_attrib = node->get_attrib("rel-transform");
  if (rel_transform_attrib)
  {
    Matrix4 T;
    rel_transform_attrib->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
    {
      std::cerr << "CollisionGeometry::load_from_xml() warning: bad transform? ";
      std::cerr << std::endl << T << " when reading node " << std::endl;
      std::cerr << *node << std::endl;
      std::cerr << "  --> possibly a floating-point error..." << std::endl;
    }
    set_rel_transform(T, true);
  }

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

  // read any sub-collision geometry nodes
  std::list<XMLTreeConstPtr> subcg_nodes = node->find_child_nodes("CollisionGeometry");
  if (!subcg_nodes.empty())
  {
    // clear all existing child collision geometry nodes
    _children.clear();

    // create child nodes and read from XML
    for (std::list<XMLTreeConstPtr>::const_iterator i = subcg_nodes.begin(); i != subcg_nodes.end(); i++)
    {
      // create the geometry
      CollisionGeometryPtr cg(new CollisionGeometry);

      // populate it from XML
      cg->load_from_xml(*i, id_map);

      // set the single body of the child
      cg->set_single_body(get_single_body());

      // add the child to this
      children.push_back(cg);
    }
  }
}

/// Implements Base::save_to_xml()
void CollisionGeometry::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // save Base data
  Base::save_to_xml(node, shared_objects);

  // set the node name
  node->name = "CollisionGeometry";

  // add the transform attribute
  node->attribs.insert(XMLAttrib("transform", _transform));

  // add the rel-transform attribute
  node->attribs.insert(XMLAttrib("rel-transform", _rel_transform));

  // save the ID of the primitive and add the primitive to the shared list
  if (_geometry)
  {
    node->attribs.insert(XMLAttrib("primitive-id", _geometry->id));
    shared_objects.push_back(_geometry);
  }

  // create nodes for the children, and save them
  for (unsigned i=0; i< _children.size(); i++)
  {
    XMLTreePtr child_node(new XMLTree("CollisionGeometry"));
    node->add_child(child_node);
    _children[i]->save_to_xml(child_node, shared_objects);
  }
}

