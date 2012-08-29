/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
#include <Moby/LinAlg.h>
#include <Moby/OBB.h>
#include <Moby/BoundingSphere.h>
#include <Moby/TriangleMeshPrimitive.h>

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

/// Creates the triangle mesh primitive
TriangleMeshPrimitive::TriangleMeshPrimitive()
{
  _convexify_inertia = false;
  _edge_sample_length = std::numeric_limits<Real>::max();
}

/// Creates the triangle mesh from a geometry file and optionally centers it
TriangleMeshPrimitive::TriangleMeshPrimitive(const string& filename, bool center) 
{
  // do not convexify inertia by default
  _convexify_inertia = false;

  // do not sample edges by default
  _edge_sample_length = std::numeric_limits<Real>::max();

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

/// Creates the triangle mesh from a geometry file and optionally centers it
TriangleMeshPrimitive::TriangleMeshPrimitive(const string& filename, const Matrix4& T, bool center) : Primitive(T) 
{ 
  // do not convexify inertia by default
  _convexify_inertia = false;

  // do not sample edges by default
  _edge_sample_length = std::numeric_limits<Real>::max();

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

/// Sets whether this primitive is used for a deformable body
void TriangleMeshPrimitive::set_deformable(bool flag)
{
  Primitive::set_deformable(flag);

  // vertices, mesh, and BVH are no longer valid 
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Vector3> >();
  _mesh_vertices.clear();
  _root = BVPtr();
  _invalidated = true;
}

/// Sets the edge sample length for this box
void TriangleMeshPrimitive::set_edge_sample_length(Real len)
{
  _edge_sample_length = len;

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
  _mesh_vertices.clear();
  _invalidated = true;
}

/// Creates the visualization for this primitive
#ifdef USE_OSG
osg::Node* TriangleMeshPrimitive::create_visualization()
{
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
    Matrix4 T_inv = IDENTITY_4x4;

    // back the transformation out of the mesh; NOTE: we have to do this
    // b/c the base Primitive class uses the transform in the visualization
    IndexedTriArray mesh = _mesh->transform(T_inv);

    // get the vertices and facets
    const std::vector<Vector3>& verts = mesh.get_vertices();
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
}
#endif

/// Centers the triangle mesh
void TriangleMeshPrimitive::center()
{
  // if there is no mesh, quit now
  if (!_mesh)
    return;

  // NOTE: we no longer do this b/c transforms are *not* applied directly
  // to the primitive
  Matrix4 Tinv = IDENTITY_4x4;

  // back the transform out of the mesh
  IndexedTriArray mesh = _mesh->transform(Tinv);

  // get the c.o.m. of this new mesh
  std::list<Triangle> tris;
  mesh.get_tris(std::back_inserter(tris));
  Vector3 centroid = CompGeom::calc_centroid_3D(tris.begin(), tris.end());

  // translate this mesh so that its origin is at its c.o.m.
  mesh = mesh.translate(-centroid);

  // re-transform the mesh
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(mesh.transform(get_transform())));

  // re-calculate mass properties 
  calc_mass_properties();

  // center-of-mass should be approximately zero
  assert(_com.norm() < NEAR_ZERO);

  // update the visualization
  update_visualization();
}

/// Implements Base::load_from_xml()
/**
 * \note if centering is done, it is done <i<before</i> any transform is applied
 */
void TriangleMeshPrimitive::load_from_xml(XMLTreeConstPtr node, map<string, BasePtr>& id_map)
{
  // load data from the Primitive 
  Primitive::load_from_xml(node, id_map);

  // check that this node name is correct
  assert(strcasecmp(node->name.c_str(), "TriangleMesh") == 0);

  // determine whether to convexify for inertial calculation
  const XMLAttrib* cvx_mesh_attr = node->get_attrib("convexify-inertia");
  if (cvx_mesh_attr)
    _convexify_inertia = cvx_mesh_attr->get_bool_value();

  // read in the edge sample length
  const XMLAttrib* esl_attr = node->get_attrib("edge-sample-length");
  if (esl_attr)
    set_edge_sample_length(esl_attr->get_real_value());

  // make sure that this Triangle array has a filename specified
  const XMLAttrib* fname_attr = node->get_attrib("filename");
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
  const XMLAttrib* center_attr = node->get_attrib("center");
  if (center_attr && center_attr->get_bool_value())
    this->center();

  // update the visualization, if necessary
  update_visualization();
}

/// Implements Base::save_to_xml()
void TriangleMeshPrimitive::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
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

    // transform the mesh w/transform backed out
    Matrix4 iT = Matrix4::inverse_transform(_T);
    IndexedTriArray mesh_xform = _mesh->transform(iT);

    // write the mesh
    mesh_xform.write_to_obj(filename);
  }
  else
    in.close();
}

/// Sets the mesh
void TriangleMeshPrimitive::set_mesh(boost::shared_ptr<const IndexedTriArray> mesh)
{
  // set the mesh
  _mesh = mesh;

  // vertices and bounding volumes are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
  _mesh_vertices.clear();
  _root = BVPtr();
  _invalidated = true;

  // map pointers to vertices
  _mesh_vertex_map.clear();
  const vector<Vector3>& verts = _mesh->get_vertices();
  for (unsigned i=0; i< verts.size(); i++)
    _mesh_vertex_map[&verts[i]] = i;

  // recalculate the mass properties
  if (!is_deformable())
    calc_mass_properties();

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
  Real volume_ints[10];

  // if there is no mesh, set things to some defaults
  if (!_mesh)
  {
    _com = ZEROS_3;
    _J = ZEROS_3x3;
    return;
  }

  // determine which mesh to use
  PolyhedronPtr poly;
  const IndexedTriArray* mesh = NULL;
  if (_convexify_inertia)
  {
    const vector<Vector3>& verts = _mesh->get_vertices();
    poly = CompGeom::calc_convex_hull_3D(verts.begin(), verts.end());
    mesh = &poly->get_mesh();
  }
  else
    mesh = _mesh.get();

  // get triangles
  std::list<Triangle> tris;
  mesh->get_tris(std::back_inserter(tris));

  // compute the centroid of the triangle mesh
  _com = CompGeom::calc_centroid_3D(tris.begin(), tris.end());

  // calculate volume integrals
  mesh->calc_volume_ints(volume_ints);

  // we'll need the volume
  const Real volume = volume_ints[0];

  // compute the mass if density is given
  if (_density)
    _mass = *_density * volume;

// NOTE: we no longer transform the inertial components, since we recompute
  // compute the center-of-mass
  // NOTE: we use this instead of the COM calculated from calc_volume_ints()
  // b/c the latter seems to be less accurate; NOTE: we need to check to see
  // why that is the case
  // Vector3 com(volume_ints[1], volume_ints[2], volume_ints[3]);
  // com *= (1.0/volume);
  //  Vector3 com = calc_com();

  // compute the inertia tensor relative to world origin
  _J(X,X) = (_mass/volume) * (volume_ints[5] + volume_ints[6]);
  _J(Y,Y) = (_mass/volume) * (volume_ints[4] + volume_ints[6]);
  _J(Z,Z) = (_mass/volume) * (volume_ints[4] + volume_ints[5]);
  _J(X,Y) = (-_mass/volume) * volume_ints[7];
  _J(Y,Z) = (-_mass/volume) * volume_ints[8];
  _J(X,Z) = (-_mass/volume) * volume_ints[9];

  // set the symmetric values for J
  _J(Y,X) = _J(X,Y);
  _J(Z,Y) = _J(Y,Z);
  _J(Z,X) = _J(X,Z);

  // rotate/scale the inertia tensor
  transform_inertia(_mass, _J, _com, _T, _J, _com);

  // if one or more values of J is NaN, don't verify anything
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=i; j<= Z; j++)
      if (std::isnan(_J(i,j)))
      {
        _J = ZEROS_3x3;
        return;
      }
}

/// Gets the pointer to the root bounding box
BVPtr TriangleMeshPrimitive::get_BVH_root()
{
  // build the bounding box if necessary
  if (!_root)
    build_BB_tree();

  return _root; 
}

/// Determines whether the point on a thick triangle is degenerate
bool TriangleMeshPrimitive::is_degen_point_on_tri(shared_ptr<AThickTri> tri, const Vector3& p)
{
  // get the barycentric coordinates for the point
  Real s, t;
  tri->tri.determine_barycentric_coords(p, s, t);

  // correct s and t, if necessary
  if (s < (Real) 0.0)
  {
    if (s < -NEAR_ZERO)
      FILE_LOG(LOG_COLDET) << "Primitive::is_degen_point_on_tri() warning- s=" << s << endl;
    s = (Real) 0.0;
  }
  if (t < (Real) 0.0)
  {
    if (t < -NEAR_ZERO)
      FILE_LOG(LOG_COLDET) << "Primitive::is_degen_point_on_tri() warning- t=" << t << endl;
    t = (Real) 0.0;
  }
  if (s + t > (Real) 1.0)
  {
    if (s + t > (Real) 1.0 + NEAR_ZERO)
      FILE_LOG(LOG_COLDET) << "Primitive::is_degen_point_on_tri() warning- s+t=" << (s+t) << endl;
    Real isum = (Real) 1.0/(s + t);
    s *= isum;
    t *= isum;
  }

  // get the feature of the determined point
  Triangle::FeatureType feat = tri->tri.determine_feature(s, t);

  // if the feature is a face, we can exit now
  if (feat == Triangle::eFace)
    return false;

  // for the other features, we'll need the indexed triangle from the original 
  // triangle mesh 
  const IndexedTriArray& mesh = *tri->mesh;
  const IndexedTri& f = mesh.get_facets()[tri->tri_idx];

  // now process the feature type
  switch (feat)
  {
    // see whether vertex A is coincident to all coplanar facets
    case Triangle::eVertexA:
      return !mesh.is_coplanar(f.a);

    // see whether vertex B is coincident to all coplanar facets
    case Triangle::eVertexB:
      return !mesh.is_coplanar(f.b);
    
    // see whether vertex C is coincident to all coplanar facets
    case Triangle::eVertexC:
      return !mesh.is_coplanar(f.c);

    // see whether edge AB is coincident to all coplanar facets
    case Triangle::eEdgeAB:
      return !mesh.is_coplanar(f.a, f.b);
      
    // see whether edge BC is coincident to all coplanar facets
    case Triangle::eEdgeBC:
      return !mesh.is_coplanar(f.b, f.c);
     
    // see whether edge AC is coincident to all coplanar facets
    case Triangle::eEdgeAC:
      return !mesh.is_coplanar(f.a, f.c);
 
    // this should already be handled
    case Triangle::eFace:
      assert(false);
      return true;

    // this should not happen
    case Triangle::eNone:
      assert(false);
      return true;
  }

  // should never get here...
  assert(false);
  return true;
}

/// Gets mesh data for the geometry with the specified bounding volume
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& TriangleMeshPrimitive::get_sub_mesh(BVPtr bv)
{
  assert(_mesh_tris.find(bv) != _mesh_tris.end());
  _smesh = make_pair(_mesh, _mesh_tris.find(bv)->second);
  return _smesh;
}

/// Determines whether a point is inside / on one of the thick triangles
bool TriangleMeshPrimitive::point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const
{
  const Real EXPANSION_CONST = 0.01;

  // expand the BV 
  BVPtr ebv;
  if (!is_deformable())
  {
    assert(dynamic_pointer_cast<OBB>(bv));
    OBB eobb = *dynamic_pointer_cast<OBB>(bv);
    eobb.l *= ((Real) 1.0 + EXPANSION_CONST);
    ebv = BVPtr(new OBB(eobb));
  }
  else
  {
    assert(dynamic_pointer_cast<BoundingSphere>(bv));
    BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
    bs.radius *= ((Real) 1.0 + EXPANSION_CONST);
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

/// Intersects a line segment against the triangle mesh and returns first point of intersection (special method for self collision checks)
bool TriangleMeshPrimitive::intersect_seg(const Vector3* u, BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const
{
  const unsigned LEAF_TRIS_CUTOFF = 5;
  const Real EXPANSION_CONST = 0.01;

  // setup statistics variables
  unsigned n_bv_tests = 0;
  unsigned n_tri_tests = 0;

  // setup first point of intersection
  Real tfirst = std::numeric_limits<Real>::max();

  // create a stack and add the BV to it
  stack<BVPtr> S;
  S.push(bv);

  // get the vertex corresponding to u and get all facets incident to u
  assert(_mesh_vertex_map.find(u) != _mesh_vertex_map.end());
  unsigned vidx = _mesh_vertex_map.find(u)->second;
  const list<unsigned>& ifacets = _mesh->get_incident_facets(vidx);
  vector<unsigned> incident_facets(ifacets.begin(), ifacets.end());
  std::sort(incident_facets.begin(), incident_facets.end());

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
        eobb.l *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      // if the list is sufficiently small or there is an intersection with
      // the BV, check all triangles
      Real tmin = 0, tmax = 1;
      Vector3 pt;
      if (tris.size() <= LEAF_TRIS_CUTOFF || (++n_bv_tests && ebv->intersects(seg, tmin, tmax, pt)))
      {
        BOOST_FOREACH(shared_ptr<AThickTri> tri, tris)
        {
          // don't do intersection test if u comes from tri
          if (std::binary_search(incident_facets.begin(), incident_facets.end(), tri->tri_idx))
            continue;

          // otherwise, check...
          n_tri_tests++;
          FILE_LOG(LOG_COLDET) << "   -- intersecting thick tri: " << endl << tri->tri <<"  and segment " << seg.first << " / " << seg.second << endl;
          Real tx;
          if (tri->intersect_seg(seg, tx, pt) && tx < tfirst)
          {
            // project the point onto the actual triangle
            Vector3 tri_normal = tri->tri.calc_normal();
            Real tri_offset = tri->tri.calc_offset(tri_normal);

            // compute P = I - n*n'
            Matrix3 P;
            Vector3::outer_prod(tri_normal, -tri_normal, &P);
            P += Matrix3::identity();
            Vector3 proj_point = P * pt;
            Real remainder = tri_offset - proj_point.dot(tri_normal);
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
        eobb.l *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      Real tmin = 0, tmax = 1;
      Vector3 pt;
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


/// Intersects a line segment against the triangle mesh and returns first point of intersection
bool TriangleMeshPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const
{
  const unsigned LEAF_TRIS_CUTOFF = 5;
  const Real EXPANSION_CONST = 0.01;

  // setup statistics variables
  unsigned n_bv_tests = 0;
  unsigned n_tri_tests = 0;

  // setup first point of intersection
  Real tfirst = std::numeric_limits<Real>::max();

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
        eobb.l *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      // if the list is sufficiently small or there is an intersection with
      // the BV, check all triangles
      Real tmin = 0, tmax = 1;
      Vector3 pt;
      if (tris.size() <= LEAF_TRIS_CUTOFF || (++n_bv_tests && ebv->intersects(seg, tmin, tmax, pt)))
      {
        BOOST_FOREACH(shared_ptr<AThickTri> tri, tris)
        {
          n_tri_tests++;
          FILE_LOG(LOG_COLDET) << "   -- intersecting thick tri: " << endl << tri->tri <<"  and segment " << seg.first << " / " << seg.second << endl;
          Real tx;
          if (tri->intersect_seg(seg, tx, pt) && tx < tfirst)
          {
            // project the point onto the actual triangle
            Vector3 tri_normal = tri->tri.calc_normal();
            Real tri_offset = tri->tri.calc_offset(tri_normal);

            // compute P = I - n*n'
            Matrix3 P;
            Vector3::outer_prod(tri_normal, -tri_normal, &P);
            P += Matrix3::identity();
            Vector3 proj_point = P * pt;
            Real remainder = tri_offset - proj_point.dot(tri_normal);
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
        eobb.l *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new OBB(eobb));
      }
      else
      {
        assert(dynamic_pointer_cast<BoundingSphere>(bv));
        BoundingSphere bs = *dynamic_pointer_cast<BoundingSphere>(bv);
        bs.radius *= ((Real) 1.0 + EXPANSION_CONST);
        ebv = BVPtr(new BoundingSphere(bs));
      }

      Real tmin = 0, tmax = 1;
      Vector3 pt;
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
void TriangleMeshPrimitive::get_vertices(BVPtr bv, vector<const Vector3*>& vertices) 
{
  // if there are no vertices, we need to build them
  if (!_vertices)
    build_BB_tree();

  // get the mesh covered by the BV
  map<BVPtr, list<unsigned> >::const_iterator v_iter = _mesh_vertices.find(bv);
  const list<unsigned>& vlist = v_iter->second;

  // get the vertex indices
  for (list<unsigned>::const_iterator i = vlist.begin(); i != vlist.end(); i++)
    vertices.push_back(&(*_vertices)[*i]);
}

/// Transforms this primitive
void TriangleMeshPrimitive::set_transform(const Matrix4& T)
{
  // determine the transformation from the old to the new transform 
  Matrix4 Trel = T * Matrix4::inverse_transform(_T);

  // go ahead and set the new transform
  Primitive::set_transform(T);

  // transform mesh
  if (_mesh)
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(Trel)));

  // vertices and bounding volumes are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
  _mesh_vertices.clear();
  _root = BVPtr();
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();
}

/// Loads the state of this primitive
void TriangleMeshPrimitive::load_state(shared_ptr<void> state)
{
  shared_ptr<TriangleMeshPrimitiveState> tps = boost::static_pointer_cast<TriangleMeshPrimitiveState>(state);
  Primitive::load_state(tps->pstate);
  _root = tps->rootBVH;
  _mesh_tris = tps->mesh_tris;
  _mesh_vertices = tps->mesh_vertices;
  _tris = tps->tris;
}

/// Saves the state of this primitive
shared_ptr<void> TriangleMeshPrimitive::save_state() const
{
  shared_ptr<TriangleMeshPrimitiveState> tps(new TriangleMeshPrimitiveState);
  tps->pstate = Primitive::save_state();
  tps->rootBVH = _root;
  tps->mesh_tris = _mesh_tris;
  tps->mesh_vertices = _mesh_vertices;
  tps->tris = _tris;

  return tps;
}

/****************************************************************************
 Methods for building bounding box trees begin 
****************************************************************************/

/// Builds an bounding volume tree (OBB or BoundingSphere) from an indexed triangle mesh using a top-down approach 
/*
 * \return the root of the bounding box tree
 */
void TriangleMeshPrimitive::build_BB_tree()
{
  const unsigned THREE_D = 3;
  BVPtr child1, child2;

  FILE_LOG(LOG_BV) << "TriangleMeshPrimitive::build_BB_tree() entered" << endl;

  // clear any existing data
  _mesh_tris.clear();
  _mesh_vertices.clear();
  _tris.clear();

  // get the vertices from the mesh
  const vector<Vector3>& vertices = _mesh->get_vertices();

  // build an BV around all vertices 
  BVPtr root;
  if (!is_deformable())
    root = BVPtr(new OBB(vertices.begin(), vertices.end()));
  else
    root = BVPtr(new BoundingSphere(vertices.begin(), vertices.end()));

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
      Vector3 axis;

      // get the i'th column of R if an OBB
      if (!is_deformable())
      {
        OBBPtr obb = dynamic_pointer_cast<OBB>(bb);
        assert(obb);
        obb->R.get_column(i, axis.begin());
      }
      else
      {
        if (i == 0) axis = Vector3(1,0,0);
        else if (i == 1) axis = Vector3(0,1,0);
        else axis = Vector3(0,0,1); 
      }

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
        ttris1.push_back(shared_ptr<AThickTri>(new AThickTri(_mesh->get_triangle(idx), _intersection_tolerance)));
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
        ttris2.push_back(shared_ptr<AThickTri>(new AThickTri(_mesh->get_triangle(idx), _intersection_tolerance)));
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
    Real vol = bb->calc_volume();
    bool erased_one = false;
    for (list<BVPtr>::iterator i = bb->children.begin(); i != bb->children.end(); )
    {
      // get the volume of this child
      Real voli = (*i)->calc_volume();
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
    if (!is_deformable())
    {
      // cast it to an OBB
      OBBPtr obb = dynamic_pointer_cast<OBB>(bb);
      assert(obb);

      // fatten the OBB
      obb->l[0] += _intersection_tolerance;    
      obb->l[1] += _intersection_tolerance;    
      obb->l[2] += _intersection_tolerance;    
    }
    else
    {
      // cast it to a bounding sphere
      shared_ptr<BoundingSphere> bs = dynamic_pointer_cast<BoundingSphere>(bb);
      assert(bs);

      // fatten the bounding sphere
      for (unsigned i=0; i< THREE_D; i++)
        bs->radius += _intersection_tolerance;
    }

    // add all children to the queue
    if (!bb->is_leaf())
      BOOST_FOREACH(BVPtr child, bb->children)
        Q.push(child);

    // wipe out userdata
    bb->userdata = shared_ptr<void>();
  }

  // save the root
  _root = root;

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
  construct_mesh_vertices(_mesh);

  FILE_LOG(LOG_BV) << "Primitive::build_BB_tree() exited" << endl;
}

/// Sets the intersection tolerance
void TriangleMeshPrimitive::set_intersection_tolerance(Real tol)
{
  Primitive::set_intersection_tolerance(tol);

  // mesh, vertices, and BVH are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Vector3> >();
  _mesh_vertices.clear();
  _root = BVPtr();
}

/// Creates the set of mesh vertices
void TriangleMeshPrimitive::construct_mesh_vertices(shared_ptr<const IndexedTriArray> mesh)
{
  const unsigned EDGES_PER_TRI = 3;

  // clear the map of mesh vertices
  _mesh_vertices.clear();

  // get the sets of vertices and facets from the mesh
  const vector<Vector3>& mesh_vertices = mesh->get_vertices();
  const vector<IndexedTri>& mesh_facets = mesh->get_facets();

  // also, we'll need the vertex-to-facet map
  vector<list<unsigned> > vf_map = mesh->determine_vertex_facet_map();

  // create a new vector of vertices
  _vertices = shared_ptr<vector<Vector3> >(new vector<Vector3>(mesh_vertices));

  // now, modify the vertices based on the intersection tolerance
  for (unsigned i=0; i< mesh_vertices.size(); i++)
  {
    // get list of facets indicent to vertex i
    const list<unsigned>& ifacets = vf_map[i];

    // setup the current normal
    Vector3 normal = ZEROS_3;

    // add all coincident normals together, then normalize
    BOOST_FOREACH(unsigned j, ifacets)
      normal += mesh->get_triangle(j).calc_normal();

    // if we can't normalize, skip this vertex 
    if (normal.norm() < NEAR_ZERO)
      continue;

    // otherwise, normalize the normal and add intersection tolerance (in dir
    // of normal) to vertex i
    normal.normalize();
    (*_vertices)[i] += normal*_intersection_tolerance;
  }

  // now, add additional samples based on edges in the mesh
  map<sorted_pair<unsigned>, list<unsigned> > edge_subsamples;
  for (unsigned i=0; i< mesh_facets.size(); i++)
  {
    // setup sorted pairs for the three edges
    sorted_pair<unsigned> e[EDGES_PER_TRI];
    e[0] = make_sorted_pair(mesh_facets[i].a, mesh_facets[i].b);
    e[1] = make_sorted_pair(mesh_facets[i].b, mesh_facets[i].c);
    e[2] = make_sorted_pair(mesh_facets[i].c, mesh_facets[i].a);

    // check the three edges
    for (unsigned j=0; j< EDGES_PER_TRI; j++)
    {
      if (edge_subsamples.find(e[j]) == edge_subsamples.end())
      {
        // edge does not already exist..  add vertices as necessary
        list<unsigned>& ess = edge_subsamples[e[j]];

        // subdivide edge to create new vertices as necessary
        queue<sorted_pair<unsigned> > q;
        q.push(e[j]);
        while (!q.empty())
        {
          // get the two vertices of the "edge"
          unsigned vi = q.front().first;
          unsigned vj = q.front().second;
          q.pop();
          const Vector3& v1 = (*_vertices)[vi];
          const Vector3& v2 = (*_vertices)[vj];

          // subdivide, adding a vertex as necessary
          if ((v1-v2).norm() > _edge_sample_length)
          {
            unsigned vk = _vertices->size();
            _vertices->push_back((v1+v2) * (Real) 0.5);
            ess.push_back(vk);
            q.push(make_sorted_pair(vi,vk));
            q.push(make_sorted_pair(vk,vj));
          }
        }
      }
    }
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
      for (unsigned k=0; k< EDGES_PER_TRI; k++)
      {
        assert(edge_subsamples.find(e[k]) != edge_subsamples.end());
        const list<unsigned>& ess = edge_subsamples[e[k]];
        vlist.insert(vlist.end(), ess.begin(), ess.end());
      }
    }
  }
}

/// Splits a collection of triangles along a splitting plane into 2 new meshes 
void TriangleMeshPrimitive::split_tris(const Vector3& point, const Vector3& normal, const IndexedTriArray& orig_mesh, const list<unsigned>& ofacets, list<unsigned>& pfacets, list<unsigned>& nfacets) 
{
  // get original vertices and facets
  const vector<Vector3>& vertices = orig_mesh.get_vertices();
  const vector<IndexedTri>& facets = orig_mesh.get_facets();

  // determine the splitting plane: ax + by + cz = d
  Real offset = Vector3::dot(point, normal);

  // determine the side of the splitting plane of the triangles
  Plane plane(normal, offset);
  BOOST_FOREACH(unsigned i, ofacets)
  {
    // get the three signed distances
    Real sa = plane.calc_signed_distance(vertices[facets[i].a]);
    Real sb = plane.calc_signed_distance(vertices[facets[i].b]);
    Real sc = plane.calc_signed_distance(vertices[facets[i].c]);
    Real min_s = std::min(sa, std::min(sb, sc));
    Real max_s = std::max(sa, std::max(sb, sc));    

    // see whether we can cleanly put the triangle into one side
    if (min_s > 0)
      pfacets.push_back(i);
    else if (max_s < 0)
      nfacets.push_back(i);
    else
    {
      // triangle is split down the middle; get its centroid
      Triangle tri(vertices[facets[i].a], vertices[facets[i].b], vertices[facets[i].c]);
      Vector3 tri_centroid = tri.calc_centroid();
      Real scent = plane.calc_signed_distance(tri_centroid);
      if (scent > 0)
        pfacets.push_back(i);
      else
        nfacets.push_back(i);
    }
  }
}

/// Splits a bounding box  along a given axis into two new bounding boxes; returns true if split successful
bool TriangleMeshPrimitive::split(shared_ptr<const IndexedTriArray> mesh, shared_ptr<BV> source, shared_ptr<BV>& tgt1, shared_ptr<BV>& tgt2, const Vector3& axis) 
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
    t_tris.push_back(mesh->get_triangle(idx)); 
  Vector3 centroid = CompGeom::calc_centroid_3D(t_tris.begin(), t_tris.end());

  // get the side of the splitting plane of the triangles
  split_tris(centroid, axis, *mesh, tris, ptris, ntris);
  if (ptris.empty() || ntris.empty())
    return false;

  // get vertices from both meshes
  vector<Vector3> pverts, nverts;
  get_vertices(*mesh, ptris.begin(), ptris.end(), std::back_inserter(pverts));
  get_vertices(*mesh, ntris.begin(), ntris.end(), std::back_inserter(nverts));

  // create two new BVs 
  if (!is_deformable())
  {
    tgt1 = OBBPtr(new OBB(pverts.begin(), pverts.end()));
    tgt2 = OBBPtr(new OBB(nverts.begin(), nverts.end()));
  }
  else
  {
    tgt1 = shared_ptr<BoundingSphere>(new BoundingSphere(pverts.begin(), pverts.end()));
    tgt2 = shared_ptr<BoundingSphere>(new BoundingSphere(nverts.begin(), nverts.end()));
  }

  // setup mesh data for the BVs
  _mesh_tris[tgt1] = ptris;
  _mesh_tris[tgt2] = ntris;

  return true;
}

/****************************************************************************
 Methods for building bounding box trees end 
****************************************************************************/

