/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Array>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osg/Geometry>
#endif
#include <stack>
#include <cctype>
#include <string>
#include <queue>
#include <iostream>
#include <fstream>
#include <Moby/Log.h>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/OBB.h>
#include <Moby/BoundingSphere.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/GJK.h>
#include <Moby/TriangleMeshPrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using std::cerr;
using std::endl;
using std::list;
using std::map;
using std::string;
using std::queue;
using std::vector;
using std::pair;
using std::make_pair;
using std::stack;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;

// TODO: determine whether the mesh makes a polyhedron

/// Creates the triangle mesh primitive
TriangleMeshPrimitive::TriangleMeshPrimitive()
{
  _convexify_inertia = false;
  _edge_sample_length = std::numeric_limits<double>::max();
}

/// Creates the triangle mesh from a geometry file and optionally centers it
TriangleMeshPrimitive::TriangleMeshPrimitive(const string& filename, bool center) 
{
  // do not convexify inertia by default
  _convexify_inertia = false;

  // do not sample edges by default
  _edge_sample_length = std::numeric_limits<double>::max();

  // construct a new triangle mesh from the filename
  if (filename.find(".obj") == filename.size() - 4)
    set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(IndexedTriArray::read_from_obj(filename))));
  else
    throw std::runtime_error("TriangleMeshPrimitive (constructor): unknown mesh file type!");

  // center the mesh
  if (center)
    this->center();

  // update the visualization, if necessary
  update_visualization();
}

/// Creates the triangle mesh from a geometry file and optionally centers it
TriangleMeshPrimitive::TriangleMeshPrimitive(const string& filename, const Pose3d& T, bool center) : Primitive(T) 
{ 
  // do not convexify inertia by default
  _convexify_inertia = false;

  // do not sample edges by default
  _edge_sample_length = std::numeric_limits<double>::max();

  // construct a new triangle mesh from the filename
  if (filename.find("obj") == filename.size() - 4)
    set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(IndexedTriArray::read_from_obj(filename))));
  else
    throw std::runtime_error("TriangleMeshPrimitive (constructor): unknown mesh file type!");

  // center the mesh
  if (center)
    this->center();

  // update the visualization, if necessary
  update_visualization();
}

/// Sets the edge sample length for this triangle mesh 
void TriangleMeshPrimitive::set_edge_sample_length(double len)
{
  _edge_sample_length = len;

  // vertices are no longer valid
  _vertices.clear();
  _mesh_vertices.clear();
}

/// Creates the visualization for this primitive
osg::Node* TriangleMeshPrimitive::create_visualization()
{
  #ifdef USE_OSG
  const unsigned X = 0, Y = 1, Z = 2;

  // create a new group to hold the geometry
  osg::Group* group = new osg::Group;

  // only create stuff if necessary
  if (_mesh)
  {
    // create necessary OSG elements for visualization
    osg::Geode* geode = new osg::Geode;
    osg::Geometry* geom = new osg::Geometry;
    geode->addDrawable(geom);
    group->addChild(geode);

    // NOTE: we no longer do this b/c transforms are *not* applied directly
    // to the primitive

    // get the inverse of the current transformation
    Transform3d T_inv = Transform3d::identity();

    // back the transformation out of the mesh; NOTE: we have to do this
    // b/c the base Primitive class uses the transform in the visualization
    IndexedTriArray mesh = _mesh->transform(T_inv);

    // get the vertices and facets
    const std::vector<Origin3d>& verts = mesh.get_vertices();
    const std::vector<IndexedTri>& facets = mesh.get_facets();

    // create the vertex array
    osg::Vec3Array* varray = new osg::Vec3Array(verts.size()); 
    for (unsigned i=0; i< verts.size(); i++)
      (*varray)[i] = osg::Vec3((float) verts[i][X], (float) verts[i][Y], (float) verts[i][Z]);
    geom->setVertexArray(varray);

    // create the faces
    for (unsigned i=0; i< facets.size(); i++)
    {
      osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
      face->push_back(facets[i].a);
      face->push_back(facets[i].b);
      face->push_back(facets[i].c);
      geom->addPrimitiveSet(face);
    }
  }

  return group;
  #else
  return NULL;
  #endif
}

/// Centers the triangle mesh
void TriangleMeshPrimitive::center()
{
  // if there is no mesh, quit now
  if (!_mesh)
    return;

  // we want to transform the mesh so that the inertial pose is coincident
  // with the primitive pose; first, prepare the transformation
  Transform3d T = Pose3d::calc_relative_pose(_jF, _F);

  // do the transformation
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(T)));

  // re-calculate mass properties 
  calc_mass_properties();

  // verify that the inertial frame is near zero
  assert(Point3d(_jF->x, GLOBAL).norm() < NEAR_ZERO);

  // update the visualization
  update_visualization();
}

/// Implements Base::load_from_xml()
/**
 * \note if centering is done, it is done <i<before</i> any transform is applied
 */
void TriangleMeshPrimitive::load_from_xml(shared_ptr<const XMLTree> node, map<string, BasePtr>& id_map)
{
  // load data from the Primitive 
  Primitive::load_from_xml(node, id_map);

  // check that this node name is correct
  assert(strcasecmp(node->name.c_str(), "TriangleMesh") == 0);

  // determine whether to convexify for inertial calculation
  XMLAttrib* cvx_mesh_attr = node->get_attrib("convexify-inertia");
  if (cvx_mesh_attr)
    _convexify_inertia = cvx_mesh_attr->get_bool_value();

  // read in the edge sample length
  XMLAttrib* esl_attr = node->get_attrib("edge-sample-length");
  if (esl_attr)
    set_edge_sample_length(esl_attr->get_real_value());

  // make sure that this Triangle array has a filename specified
  XMLAttrib* fname_attr = node->get_attrib("filename");
  if (!fname_attr)
  {
    cerr << "TriangleMeshPrimitive::load_from_xml() - trying to load a ";
    cerr << " triangle mesh w/o a filename!" << endl;
    cerr << "  offending node: " << endl << *node << endl;
    return;
  }

  // setup the file extensions
  const char* OBJ_EXT = ".obj";

  // get the filename
  string fname(fname_attr->get_string_value());

  // get the lowercase version of the filename
  string fname_lower = fname;
  std::transform(fname_lower.begin(), fname_lower.end(), fname_lower.begin(), (int(*)(int)) std::tolower);

  // get the type of file and construct the triangle mesh appropriately
  if (fname_lower.find(string(OBJ_EXT)) == fname_lower.size() - strlen(OBJ_EXT))
    set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(IndexedTriArray::read_from_obj(fname))));
  else
  {
    cerr << "TriangleMeshPrimitive::load_from_xml() - unrecognized filename extension" << endl;
    cerr << "  for attribute 'filename'.  Valid extensions are '.obj' (Wavefront OBJ)" << endl;
  }
  
  // see whether to center the mesh
  XMLAttrib* center_attr = node->get_attrib("center");
  if (center_attr && center_attr->get_bool_value())
    this->center();

  // recompute mass properties
  calc_mass_properties();

  // update the visualization, if necessary
  update_visualization();
}

/// Implements Base::save_to_xml()
void TriangleMeshPrimitive::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call this parent's save_to_xml() method
  Primitive::save_to_xml(node, shared_objects);

  // set the name for this node
  node->name = "TriangleMesh";

  // add the edge sample length as an attribute
  node->attribs.insert(XMLAttrib("edge-sample-length", _edge_sample_length));

  // save convexification for inertial calculation
  node->attribs.insert(XMLAttrib("convexify-inertia", _convexify_inertia));

  // make a filename using "this"
  const unsigned MAX_DIGITS = 28;
  char buffer[MAX_DIGITS+1];
  sprintf(buffer, "%p", this);
  string filename = "triarray" + string(buffer) + ".obj";

  // add the filename as an attribute
  node->attribs.insert(XMLAttrib("filename", filename));

  // do not save the array to the OBJ file if it already exists (which we
  // crudely check for using std::ifstream to avoid OS-specific calls -- note
  // that it is possible that opening a file may fails for other reasons than
  // the file does not exist)
  std::ifstream in(filename.c_str());
  if (in.fail())
  {
    // make sure there is a mesh to write
    if (!_mesh)
      return;

    // get the transform from the global pose to the primitive pose
    shared_ptr<const Pose3d> P = get_pose();
    Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);
    IndexedTriArray mesh_xform = _mesh->transform(T);

    // write the mesh
    mesh_xform.write_to_obj(filename);
  }
  else
    in.close();
}

/// Sets the mesh
void TriangleMeshPrimitive::set_mesh(boost::shared_ptr<const IndexedTriArray> mesh)
{
  // TODO: remove this calculation (replaces mesh with convex hull)
  const vector<Origin3d>& verts = mesh->get_vertices();
  PolyhedronPtr poly = CompGeom::calc_convex_hull(verts.begin(), verts.end());
  _mesh = shared_ptr<const IndexedTriArray>(new IndexedTriArray(poly->get_mesh()));

  // TODO: restore this
  // set the mesh
//  _mesh = mesh;

  // vertices and bounding volumes are no longer valid
  _vertices.clear();
  _mesh_vertices.clear();
  _roots.clear();;

  // update visualization
  update_visualization();
}

/// Calculates mass properties of this primitive
/**
 * Computes the mass, center-of-mass, and inertia of this primitive.
 */
void TriangleMeshPrimitive::calc_mass_properties()
{
  const unsigned X = 0, Y = 1, Z = 2;
  double volume_ints[10];

  // if there is no mesh, set things to some defaults
  if (!_mesh)
  {
    _J.m = 0.0;
    _J.h.set_zero();
    _J.J.set_zero();
    return;
  }

  // determine which mesh to use
  PolyhedronPtr poly;
  const IndexedTriArray* mesh = NULL;
  if (_convexify_inertia)
  {
    const vector<Origin3d>& verts = _mesh->get_vertices();
    poly = CompGeom::calc_convex_hull(verts.begin(), verts.end());
    mesh = &poly->get_mesh();
  }
  else
    mesh = _mesh.get();

  // get triangles
  std::list<Triangle> tris;
  mesh->get_tris(std::back_inserter(tris), get_pose());

  // compute the centroid of the triangle mesh
  _jF->x = CompGeom::calc_centroid_3D(tris.begin(), tris.end());

  // calculate volume integrals
  mesh->calc_volume_ints(volume_ints);

  // we'll need the volume
  const double volume = volume_ints[0];

  // compute the mass if density is given
  if (_density)
    _J.m = *_density * volume;

// NOTE: we no longer transform the inertial components, since we recompute
  // compute the center-of-mass
  // NOTE: we use this instead of the COM calculated from calc_volume_ints()
  // b/c the latter seems to be less accurate; NOTE: we need to check to see
  // why that is the case
  // Point3d com(volume_ints[1], volume_ints[2], volume_ints[3]);
  // com *= (1.0/volume);
  //  Point3d com = calc_com();

  // compute the inertia tensor relative to world origin
  _J.J(X,X) = (_J.m/volume) * (volume_ints[5] + volume_ints[6]);
  _J.J(Y,Y) = (_J.m/volume) * (volume_ints[4] + volume_ints[6]);
  _J.J(Z,Z) = (_J.m/volume) * (volume_ints[4] + volume_ints[5]);
  _J.J(X,Y) = (-_J.m/volume) * volume_ints[7];
  _J.J(Y,Z) = (-_J.m/volume) * volume_ints[8];
  _J.J(X,Z) = (-_J.m/volume) * volume_ints[9];

  // set the symmetric values for J
  _J.J(Y,X) = _J.J(X,Y);
  _J.J(Z,Y) = _J.J(Y,Z);
  _J.J(Z,X) = _J.J(X,Z);

  // if one or more values of J is NaN, don't verify anything
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=i; j<= Z; j++)
      if (std::isnan(_J.J(i,j)))
      {
        _J.J.set_zero();
        return;
      }
}

/// Gets the pointer to the root bounding box
BVPtr TriangleMeshPrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  // build the bounding box if necessary
  BVPtr& root = _roots[geom]; 
  if (!root)
    build_BB_tree(geom);

  return root; 
}

/// Returns whether the mesh is convex (currently mesh must be convex)
bool TriangleMeshPrimitive::is_convex() const
{
  return true;
}

/// Computes the signed distance to a point from the mesh 
double TriangleMeshPrimitive::calc_signed_dist(const Point3d& p) const
{
  // verify that the point is defined with respect to one of the poses
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  // if the primitive is convex and the point is outside, return the closest
  // distance
  if (is_convex())
  {
    // loop through all triangles of the mesh
    bool inside = false;
    double min_dist = std::numeric_limits<double>::max();
    double max_dist = -min_dist;
    for (unsigned i=0; i< _mesh->num_tris(); i++)
    {
      double dist = _mesh->get_triangle(i, p.pose).calc_signed_dist(p);
      if (dist < 0.0)
      {
        max_dist = std::max(max_dist, dist);
        inside = true;
      }
      else
        min_dist = std::min(min_dist, dist);
    }

    // if the point is outside, return the minimum distance
    if (!inside)
      return min_dist;
    else
      return max_dist;
  } 
  else
  {
    throw std::runtime_error("TriangleMeshPrimitive::calc_signed_dist() - non-convex meshes not currently supported!");
  }
}

/// Computes the distance and normal from a point on the mesh 
double TriangleMeshPrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  // setup the normal
  normals.push_back(Vector3d());
  Vector3d& normal = normals.back();

  // verify that the point is defined with respect to one of the poses
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  // if the primitive is convex and the point is outside, return the closest
  // distance
  if (is_convex())
  {
    // loop through all triangles of the mesh
    bool inside = false;
    double min_dist = std::numeric_limits<double>::max();
    double max_dist = -min_dist;
    for (unsigned i=0; i< _mesh->num_tris(); i++)
    {
      double dist = _mesh->get_triangle(i, p.pose).calc_signed_dist(p);
      if (dist < 0.0)
      {
        if (dist > max_dist)
        {
          normal = _mesh->get_triangle(i, p.pose).calc_normal();
          max_dist = dist;
        }
        inside = true;
      }
      else if (!inside)
      {
        if (dist < min_dist)
        {
          normal = _mesh->get_triangle(i, p.pose).calc_normal();
          min_dist = dist;
        }
      }
    }

    // if the point is outside, return the minimum distance
    if (!inside)
      return min_dist;
    else
      return max_dist;
  } 
  else
  {
    throw std::runtime_error("TriangleMeshPrimitive::calc_signed_dist() - non-convex meshes not currently supported!");
  }
}

/// Computes the distance and normal from a point on the mesh 
double TriangleMeshPrimitive::calc_signed_dist(shared_ptr<const Primitive> primitive, Point3d& pthis, Point3d& pprimitive) const
{
  if (is_convex() && primitive->is_convex())
  {
    shared_ptr<const Pose3d> Ptri = pthis.pose;
    shared_ptr<const Pose3d> Pgeneric = pprimitive.pose;
    shared_ptr<const Primitive> tthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return GJK::do_gjk(tthis, primitive, Ptri, Pgeneric, pthis, pprimitive);
  }
  else
  {
    throw std::runtime_error("TriangleMeshPrimitive::calc_signed_dist() - non-convex meshes not currently supported!");
  }
}

/// Determines whether a point is inside / on one of the thick triangles
/*
bool TriangleMeshPrimitive::point_inside(BVPtr bv, const Point3d& p, Vector3d& normal) const
{
  const double EXPANSION_CONST = 0.01;

  // expand the BV 
  BVPtr ebv;
  if (!is_deformable())
  {
    assert(dynamic_pointer_cast<OBB>(bv));
    OBB eobb = *dynamic_pointer_cast<OBB>(bv);
    eobb.l *= ((double) 1.0 + EXPANSION_CONST);
    ebv = BVPtr(new OBB(eobb));
  }
  else
  {
    assert(dynamic_pointer_cast<BoundingSphere>(bv));
    BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
    bs.radius *= ((double) 1.0 + EXPANSION_CONST);
    ebv = BVPtr(new BoundingSphere(bs));
  }

  // handle the two cases: leaf and internal node
  if (bv->is_leaf())
  {
    // get the triangles of the BV
    assert(_tris.find(bv) != _tris.end());
    const list<shared_ptr<AThickTri> >& tris = _tris.find(bv)->second;
    
    // see whether the point is inside/on one of the thick triangles
    BOOST_FOREACH(shared_ptr<AThickTri> tri, tris)
      if (tri->point_inside(p) && !is_degen_point_on_tri(tri, p))
      {
        normal = tri->determine_normal(p);
        return true;
      }

    // still here?  not inside/on any of the triangles
    return false;
  }
  else
  {
    BOOST_FOREACH(BVPtr child, bv->children)
      if (!ebv->outside(p) && point_inside(child, p, normal))
        return true;
  }  

  // still here?  not inside/on...
  return false;
}
*/

/// Intersects a line segment against the triangle mesh and returns first point of intersection
/*
bool TriangleMeshPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  const unsigned LEAF_TRIS_CUTOFF = 5;
  const double EXPANSION_CONST = 0.01;

  // setup statistics variables
  unsigned n_bv_tests = 0;
  unsigned n_tri_tests = 0;

  // setup first point of intersection
  double tfirst = std::numeric_limits<double>::max();

  // create a stack and add the BV to it
  stack<BVPtr> S;
  S.push(bv);

  // process until stack is empty
  while (!S.empty())
  {
    // get the BV from the top of the stack
    bv = S.top();
    S.pop();

    // if the BV is a leaf, add all thick triangles to output if there are
    // not many thick triangles; otherwise, see whether the line segment
    // intersects the OBB
    if (bv->is_leaf())
    {
      // get the list of thick triangles 
      assert(_tris.find(bv) != _tris.end());
      const list<shared_ptr<AThickTri> >& tris = _tris.find(bv)->second;
    
      // expand the BV 
      BVPtr ebv;
      if (!is_deformable())
      {
        assert(dynamic_pointer_cast<OBB>(bv));
        OBB eobb = *dynamic_pointer_cast<OBB>(bv);
        eobb.l *= ((double) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((double) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      // if the list is sufficiently small or there is an intersection with
      // the BV, check all triangles
      double tmin = 0, tmax = 1;
      Point3d pt;
      if (tris.size() <= LEAF_TRIS_CUTOFF || (++n_bv_tests && ebv->intersects(seg, tmin, tmax, pt)))
      {
        BOOST_FOREACH(shared_ptr<AThickTri> tri, tris)
        {
          n_tri_tests++;
          FILE_LOG(LOG_COLDET) << "   -- intersecting thick tri: " << endl << tri->tri <<"  and segment " << seg.first << " / " << seg.second << endl;
          double tx;
          if (tri->intersect_seg(seg, tx, pt) && tx < tfirst)
          {
            // project the point onto the actual triangle
            Vector3d tri_normal = tri->tri.calc_normal();
            double tri_offset = tri->tri.calc_offset(tri_normal);

            // compute P = I - n*n'
            Matrix3d P;
            Opsd::outer_prod(tri_normal, -tri_normal, P);
            P += Matrix3d::identity();
            Point3d proj_point(P * Origin3d(pt), tri_normal.pose);
            double remainder = tri_offset - proj_point.dot(tri_normal);
            proj_point += tri_normal * remainder;

            // ensure that projected point does not occur on vertex / edge
            if (!is_degen_point_on_tri(tri, proj_point))
            {
              tfirst = tx;
              isect = pt;
              normal = tri->determine_normal(isect);
              FILE_LOG(LOG_COLDET) << "  intersection detected!  time of impact: " << t << endl;
            }
            else
            {
              FILE_LOG(LOG_COLDET) << "  intersection detected, but point of intersection on vertex/edge" << endl;
              FILE_LOG(LOG_COLDET) << "    intersection point: " << isect << endl;
            }
          }
        }
      }
    }
    else
    {
      // copy the BV to an expanded one
      BVPtr ebv;
      if (!is_deformable())
      {
        assert(dynamic_pointer_cast<OBB>(bv));
        OBB eobb = *dynamic_pointer_cast<OBB>(bv);
        eobb.l *= ((double) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((double) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      double tmin = 0, tmax = 1;
      Point3d pt;
      n_bv_tests++;
      if (ebv->intersects(seg, tmin, tmax, pt))
      {
        FILE_LOG(LOG_COLDET) << " -- internal BV " << bv << " and line segment intersect; testing " << bv->children.size() << " child OBB nodes" << std::endl;
        BOOST_FOREACH(BVPtr child, bv->children)
        {
          FILE_LOG(LOG_COLDET) << "    -- child: " << child << std::endl;
          S.push(child);
        }
      }
    }
  }

  // look for valid intersection
  t = tfirst;

  FILE_LOG(LOG_COLDET) << "# of BV / line segment tests: " << n_bv_tests << endl;
  FILE_LOG(LOG_COLDET) << "# of tri / line segment tests: " << n_tri_tests << endl;

  return (tfirst > -NEAR_ZERO && tfirst < 1.0 + NEAR_ZERO);
}

/// Gets vertices corresponding to the bounding volume
void TriangleMeshPrimitive::get_vertices(BVPtr bv, vector<const Point3d*>& vertices) 
{
  // get the vertices for this geometry
  vector<Point3d>& verts = _vertices[bv->geom];

  // get the mesh covered by the BV
  map<BVPtr, list<unsigned> >::const_iterator v_iter = _mesh_vertices.find(bv);
  const list<unsigned>& vlist = v_iter->second;

  // get the vertex indices
  for (list<unsigned>::const_iterator i = vlist.begin(); i != vlist.end(); i++)
    vertices.push_back(&verts[*i]);
}
*/

/// Gets vertices corresponding to a particular pose 
void TriangleMeshPrimitive::get_vertices(shared_ptr<const Pose3d> P, vector<Point3d>& vertices) const 
{
  // verify that the primitive knows about this mesh
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // get the mesh vertices
  vertices = _vertices; 

  // set the pose for each vertex
  for (unsigned i=0; i< vertices.size(); i++)
    vertices[i].pose = P;
}

/// Transforms this primitive
void TriangleMeshPrimitive::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the global frame to the old pose
  Transform3d cTg = Pose3d::calc_relative_pose(GLOBAL, _F);

  // determine the transformation from the old to the new pose
  Transform3d xTc = Pose3d::calc_relative_pose(_F, x);

  // determine the transformation from the new pose to the global frame 
  Transform3d gTx = Pose3d::calc_relative_pose(x, GLOBAL);

  // compute the transformation
  Transform3d T = gTx * xTc * cTg;

  // go ahead and set the new transform
  Primitive::set_pose(p);

  // reset mesh, vertices, and bounding volumes 
  _mesh.reset();
  _vertices.clear();
  _mesh_vertices.clear();
  _roots.clear();

  // recalculate the mass properties
  calc_mass_properties();
}

/****************************************************************************
 Methods for building bounding box trees begin 
****************************************************************************/

/// Builds an bounding volume tree (OBB or BoundingSphere) from an indexed triangle mesh using a top-down approach 
/*
 * \return the root of the bounding box tree
 */
void TriangleMeshPrimitive::build_BB_tree(CollisionGeometryPtr geom)
{
  const unsigned THREE_D = 3;
  BVPtr child1, child2;

  FILE_LOG(LOG_BV) << "TriangleMeshPrimitive::build_BB_tree() entered" << endl;

  // clear any existing data for BVs for this geometry
  BVPtr& geom_root = _roots[geom];
  if (geom_root)
  {
    std::queue<BVPtr> q;
    q.push(geom_root);
    while (!q.empty())
    {
      // get the bv off of the front of the queue
      BVPtr bv = q.front();
      q.pop();

      // add children of the bv to the queue
      BOOST_FOREACH(BVPtr child, bv->children)
        q.push(child);

      // clear data associated with this bv
      _mesh_tris.erase(bv);
      _mesh_vertices.erase(bv);
      _tris.erase(bv);
    }
  }

  // get the vertices from the mesh
  const vector<Origin3d>& verts = _mesh->get_vertices();

  // transform the vertices into Point3d objects
  vector<Point3d> vertices(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = Point3d(verts[i], GLOBAL);

  // build an BV around all vertices 
  BVPtr root = BVPtr(new OBB(vertices.begin(), vertices.end()));

  // point the root to the geometry
  root->geom = geom;

  // set root to point to all facet indices
  list<unsigned> tris_idx;
  for (unsigned i=0; i< _mesh->get_facets().size(); i++)
    tris_idx.push_back(i);

  // setup mapping from BV to mesh
  _mesh_tris[root] = tris_idx;

  FILE_LOG(LOG_BV) << "  -- created root: " << root << endl;

  // add the root to a splitting queue
  queue<BVPtr> Q;
  Q.push(root);

  // split until we can't split further
  while (!Q.empty())
  {
    // get the bounding box off of the top of the queue
    BVPtr bb = Q.front();
    Q.pop();

    FILE_LOG(LOG_BV) << "  -- splitting BV: " << bb << endl;

    // split the bounding box across each of the three axes
    for (unsigned i=0; i< 3; i++)
    {
      Vector3d axis(GLOBAL);

      // get the i'th column of R
      OBBPtr obb = dynamic_pointer_cast<OBB>(bb);
      assert(obb);
      obb->R.get_column(i, axis);

      // split the bounding box across the axis
      if (split(_mesh, bb, child1, child2, axis))
        break;
    }

    // make sure that this BV was divisible
    if (!child1)
      continue;

    // child was divisible; remove thick triangles
    _tris.erase(bb);

    // setup child pointers
    bb->children.push_back(child1);
    bb->children.push_back(child2);

    // get lists of triangles for children
    assert(_mesh_tris.find(child1) != _mesh_tris.end());
    assert(_mesh_tris.find(child2) != _mesh_tris.end());
    const std::list<unsigned>& c1tris = _mesh_tris.find(child1)->second;
    const std::list<unsigned>& c2tris = _mesh_tris.find(child2)->second;

    // create thick triangles for child1
    list<shared_ptr<AThickTri> >& ttris1 = _tris[child1];
    BOOST_FOREACH(unsigned idx, c1tris)
    {
      try
      {
        ttris1.push_back(shared_ptr<AThickTri>(new AThickTri(_mesh->get_triangle(idx, get_pose()), 0.0)));
        ttris1.back()->mesh = _mesh;
        ttris1.back()->tri_idx = idx;
      }
      catch (NumericalException e)
      {
        // we won't do anything...  we just won't add the triangle
      }
    }
    
    // create thick triangles for child2
    list<shared_ptr<AThickTri> >& ttris2 = _tris[child2];
    BOOST_FOREACH(unsigned idx, c2tris)
    {
      try
      {
        ttris2.push_back(shared_ptr<AThickTri>(new AThickTri(_mesh->get_triangle(idx, get_pose()), 0.0)));
        ttris2.back()->mesh = _mesh;
        ttris2.back()->tri_idx = idx;
      }
      catch (NumericalException e)
      {
        // we won't do anything...  we just won't add the triangle
      }
    }

    // add children to the queue for processing if they have more than one tri
    if (c1tris.size() > 1)
      Q.push(child1);
    if (c2tris.size() > 1)
      Q.push(child2);

    FILE_LOG(LOG_BV) << "  -- built children " << child1 << ", " << child2 << " from " << bb << endl;
    FILE_LOG(LOG_BV) << "  -- ID: " << child1 << endl;
    FILE_LOG(LOG_BV) << "  -- ID: " << child2 << endl;
  }

  // now, collapse the tree
  Q.push(root);
  while (!Q.empty())
  {
    // for any children with a greater volume than the obb in question,
    // remove the grandchildren and add them as children
    BVPtr bb = Q.front();
    double vol = bb->calc_volume();
    bool erased_one = false;
    for (list<BVPtr>::iterator i = bb->children.begin(); i != bb->children.end(); )
    {
      // get the volume of this child
      double voli = (*i)->calc_volume();
      if (!(*i)->is_leaf() && voli > vol + NEAR_ZERO)
      {
        erased_one = true;
        BOOST_FOREACH(BVPtr gchild, (*i)->children)
          bb->children.push_back(gchild);
        i = (*i)->children.erase(i);
      }
      else
        i++;
    }

    if (!erased_one)
    {
      Q.pop();
      BOOST_FOREACH(BVPtr child, bb->children)
        if (!child->is_leaf())
          Q.push(child);
    }
  }

  Q.push(root);
  while (!Q.empty())
  {
    // get the element off of the front of the queue
    BVPtr bb = Q.front();
    Q.pop();

    // fatten the bounding volume
    // cast it to an OBB
    OBBPtr obb = dynamic_pointer_cast<OBB>(bb);
    assert(obb);

    // add all children to the queue
    if (!bb->is_leaf())
    {
      BOOST_FOREACH(BVPtr child, bb->children)
        Q.push(child);
    }

    // wipe out userdata
    bb->userdata = shared_ptr<void>();
  }

  // save the root
  geom_root = root;

  // output how many triangles are in each bounding box
  if (LOGGING(LOG_BV))
  {
    stack<pair<BVPtr, unsigned> > S;
    S.push(make_pair(root, 0));
    while (!S.empty())
    {
      // get the node off of the top of the stack
      BVPtr node = S.top().first;
      unsigned depth = S.top().second;
      S.pop();
      
      // get the triangles in this BV
      const list<unsigned>& tris = _mesh_tris.find(node)->second;
      std::ostringstream out;
      for (unsigned i=0; i< depth; i++)
        out << " ";
      out << "triangles covered by BV: " << tris.size() << endl;
      FILE_LOG(LOG_BV) << out.str();

      // put all children onto the stack
      BOOST_FOREACH(BVPtr child, node->children)
        S.push(make_pair(child, depth+1));
    }
  }

  // build set of mesh vertices
  construct_mesh_vertices(_mesh, geom);

  FILE_LOG(LOG_BV) << "Primitive::build_BB_tree() exited" << endl;
}

/// Creates the set of mesh vertices
void TriangleMeshPrimitive::construct_mesh_vertices(shared_ptr<const IndexedTriArray> mesh, CollisionGeometryPtr geom)
{
  const unsigned EDGES_PER_TRI = 3;

  // get the sets of vertices and facets from the mesh
  const vector<Origin3d>& mesh_vertices = mesh->get_vertices();
  const vector<IndexedTri>& mesh_facets = mesh->get_facets();

  // also, we'll need the vertex-to-facet map
  vector<list<unsigned> > vf_map = mesh->determine_vertex_facet_map();

  // now, modify the vertices based on the intersection tolerance
  for (unsigned i=0; i< mesh_vertices.size(); i++)
  {
    // get list of facets indicent to vertex i
    const list<unsigned>& ifacets = vf_map[i];

    // setup the current normal
    Vector3d normal = Vector3d::zero(GLOBAL);

    // add all coincident normals together, then normalize
    BOOST_FOREACH(unsigned j, ifacets)
      normal += mesh->get_triangle(j, GLOBAL).calc_normal();

    // if we can't normalize, skip this vertex 
    if (normal.norm() < NEAR_ZERO)
      continue;

    // otherwise, normalize the normal and add intersection tolerance (in dir
    // of normal) to vertex i
    normal.normalize();
  }

  // iterate over all mesh triangles
  for (map<BVPtr, list<unsigned> >::const_iterator i = _mesh_tris.begin(); i != _mesh_tris.end(); i++)
  {
    // get the list of facets
    const list<unsigned>& covered_facets = i->second;

    // create the list of vertices for this BV
    list<unsigned>& vlist = _mesh_vertices[i->first];

    // get the edges referenced by each facet
    BOOST_FOREACH(unsigned j, covered_facets)
    {
      // add vertices a, b, c of the covered facet to the list
      vlist.push_back(mesh_facets[j].a);
      vlist.push_back(mesh_facets[j].b);
      vlist.push_back(mesh_facets[j].c);

      // add all vertex samples to the list
      sorted_pair<unsigned> e[EDGES_PER_TRI];
      e[0] = make_sorted_pair(mesh_facets[j].a, mesh_facets[j].b);
      e[1] = make_sorted_pair(mesh_facets[j].b, mesh_facets[j].c);
      e[2] = make_sorted_pair(mesh_facets[j].c, mesh_facets[j].a);
    }
  }
}

/// Splits a collection of triangles along a splitting plane into 2 new meshes 
void TriangleMeshPrimitive::split_tris(const Point3d& point, const Vector3d& normal, const IndexedTriArray& orig_mesh, const list<unsigned>& ofacets, list<unsigned>& pfacets, list<unsigned>& nfacets) 
{
  // get original vertices and facets
  const vector<Origin3d>& vertices = orig_mesh.get_vertices();
  const vector<IndexedTri>& facets = orig_mesh.get_facets();

  // determine the splitting plane: ax + by + cz = d
  double offset = Vector3d::dot(point, normal);

  // determine the side of the splitting plane of the triangles
  Plane plane(normal, offset);
  BOOST_FOREACH(unsigned i, ofacets)
  {
    // setup the vertices in the same pose as the normal
    Point3d pa(vertices[facets[i].a], normal.pose);
    Point3d pb(vertices[facets[i].b], normal.pose);
    Point3d pc(vertices[facets[i].c], normal.pose);

    // get the three signed distances
    double sa = plane.calc_signed_distance(pa);
    double sb = plane.calc_signed_distance(pb);
    double sc = plane.calc_signed_distance(pc);
    double min_s = std::min(sa, std::min(sb, sc));
    double max_s = std::max(sa, std::max(sb, sc));    

    // see whether we can cleanly put the triangle into one side
    if (min_s > 0)
      pfacets.push_back(i);
    else if (max_s < 0)
      nfacets.push_back(i);
    else
    {
      // triangle is split down the middle; get its centroid
      Triangle tri(pa, pb, pc);
      Point3d tri_centroid = tri.calc_centroid();
      double scent = plane.calc_signed_distance(tri_centroid);
      if (scent > 0)
        pfacets.push_back(i);
      else
        nfacets.push_back(i);
    }
  }
}

/// Splits a bounding box  along a given axis into two new bounding boxes; returns true if split successful
bool TriangleMeshPrimitive::split(shared_ptr<const IndexedTriArray> mesh, shared_ptr<BV> source, shared_ptr<BV>& tgt1, shared_ptr<BV>& tgt2, const Vector3d& axis) 
{
  // setup two lists of triangles
  list<unsigned> ptris, ntris;

  // clear both targets
  tgt1 = shared_ptr<BV>();
  tgt2 = shared_ptr<BV>();

  // get the mesh and the list of triangles
  assert(_mesh_tris.find(source) != _mesh_tris.end());
  const list<unsigned>& tris = _mesh_tris.find(source)->second;

  // make sure that not trying to split a single triangle
  assert(tris.size() > 1); 

  // determine the centroid of this set of triangles
  list<Triangle> t_tris;
  BOOST_FOREACH(unsigned idx, tris)
    t_tris.push_back(mesh->get_triangle(idx, source->get_relative_pose())); 
  Point3d centroid = CompGeom::calc_centroid_3D(t_tris.begin(), t_tris.end());

  // get the side of the splitting plane of the triangles
  split_tris(centroid, axis, *mesh, tris, ptris, ntris);
  if (ptris.empty() || ntris.empty())
    return false;

  // get vertices from both meshes
  vector<Point3d> pverts, nverts;
  get_vertices(*mesh, ptris.begin(), ptris.end(), std::back_inserter(pverts), centroid.pose);
  get_vertices(*mesh, ntris.begin(), ntris.end(), std::back_inserter(nverts), centroid.pose);

  // create two new BVs 
  tgt1 = OBBPtr(new OBB(pverts.begin(), pverts.end()));
  tgt2 = OBBPtr(new OBB(nverts.begin(), nverts.end()));

  // setup geometry pointers
  tgt1->geom = source->geom;
  tgt2->geom = source->geom;

  // setup mesh data for the BVs
  _mesh_tris[tgt1] = ptris;
  _mesh_tris[tgt2] = ntris;

  return true;
}

/****************************************************************************
 Methods for building bounding box trees end 
****************************************************************************/

