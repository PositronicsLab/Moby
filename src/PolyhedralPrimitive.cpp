/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/MatrixTransform>
#include <osg/Geode>
#include <osg/Geometry>
#endif
#include <Moby/SpherePrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/GJK.h>
#include <Moby/XMLTree.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/AABB.h>
#include <Moby/PolyhedralPrimitive.h>

using std::cerr;
using std::string;
using std::endl;
using std::list;
using std::vector;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;
using boost::shared_ptr;
using boost::weak_ptr;
using namespace Ravelin;
using namespace Moby;

/// The pointer to the test vectors
vector<Vector3d>* PolyhedralPrimitive::test_vecs;

/// Special comparison routine for separating axis test 
bool PolyhedralPrimitive::compare_vecs(unsigned i, unsigned j)
{
  return (std::fabs((*test_vecs)[i][0]) < std::fabs((*test_vecs)[j][0]));
}

/// Gets face normals for the polyhedron
void PolyhedralPrimitive::add_to_face_vector(const Transform3d& wTe, vector<Vector3d>& normals) const
{
  // get the polyhedron
  const Polyhedron& poly = get_polyhedron();

  // get the vector of faces 
  const vector<shared_ptr<Polyhedron::Face> >& faces = poly.get_faces();

  // setup normals
  const unsigned NORMALS_ST = normals.size();
  normals.resize(NORMALS_ST + faces.size());
  for (unsigned i=0; i< faces.size(); i++) 
    normals[i+NORMALS_ST] = wTe.transform_vector(Vector3d(faces[i]->get_plane().get_normal().data(), wTe.source));  
}

/// Gets edge direction vectors for the polyhedron
void PolyhedralPrimitive::create_edge_vector(const Transform3d& wTe, vector<Vector3d>& edge_vectors) const
{
  // get the polyhedron
  const Polyhedron& poly = get_polyhedron();

  // get the vector of edges 
  const vector<shared_ptr<Polyhedron::Edge> >& edges = poly.get_edges();

  // clear the edge vector
  edge_vectors.resize(edges.size()); 

  // setup a loop index
  unsigned i = 0;
  BOOST_FOREACH(boost::shared_ptr<Polyhedron::Edge> edge, edges)
  {
    Ravelin::Vector3d v1_e(edge->v1->o,wTe.source);
    Ravelin::Vector3d v2_e(edge->v2->o,wTe.source);

    Ravelin::Vector3d v1_w = wTe.transform_point(v1_e);
    Ravelin::Vector3d v2_w = wTe.transform_point(v2_e);
    
    Ravelin::Vector3d edge_v = (v2_w - v1_w);
    edge_v.normalize();
    edge_vectors[i++] = edge_v;
  }
}

/// Determines whether the primitive is convex
bool PolyhedralPrimitive::is_convex() const
{
  // abuse the const keyword
  Polyhedron* pnc = (Polyhedron*) &_poly;
  return pnc->is_convex();
}

double PolyhedralPrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end()); 

  // see whether the point is inside or outside the primitive
  unsigned closest_facet;
  double dist = _poly.calc_signed_distance(Origin3d(p), closest_facet);

  // get the closest feature to the point
  shared_ptr<Polyhedron::Feature> closest_feature = _poly.find_closest_feature(Origin3d(p), closest_facet); 

  // try to cast it as a vertex first
  shared_ptr<Polyhedron::Vertex> v = dynamic_pointer_cast<Polyhedron::Vertex>(closest_feature);
  if (v)
  {
    // get all faces incident to the vertex
    for (list<weak_ptr<Polyhedron::Edge> >::const_iterator ei = v->e.begin(); ei != v->e.end(); ei++)
    {
      shared_ptr<Polyhedron::Edge> e(*ei);

      // get both faces incident to the edge
      shared_ptr<Polyhedron::Face> f1(e->face1);
      shared_ptr<Polyhedron::Face> f2(e->face2);

      // get the planes
      Plane plane1 = f1->get_plane();
      Plane plane2 = f2->get_plane();
      Vector3d n1 = plane1.get_normal();
      Vector3d n2 = plane2.get_normal();
      n1.pose = p.pose;
      n2.pose = p.pose;
      normals.push_back(n1);      
      normals.push_back(n2);      
    } 
  }
  else
  {
    shared_ptr<Polyhedron::Edge> e = dynamic_pointer_cast<Polyhedron::Edge>(closest_feature);
    if (e)
    {
      // get both faces incident to the edge
      shared_ptr<Polyhedron::Face> f1(e->face1);
      shared_ptr<Polyhedron::Face> f2(e->face2);

      // get the planes
      Plane plane1 = f1->get_plane();
      Plane plane2 = f2->get_plane();
      Vector3d n1 = plane1.get_normal();
      Vector3d n2 = plane2.get_normal();
      n1.pose = p.pose;
      n2.pose = p.pose;
      normals.push_back(n1);      
      normals.push_back(n2);      
    }
    else
    {
      shared_ptr<Polyhedron::Face> f = dynamic_pointer_cast<Polyhedron::Face>(closest_feature);
      assert(f);

      // get the plane
      Plane plane = f->get_plane();
      Vector3d n = plane.get_normal();
      n.pose = p.pose;
      normals.push_back(n);      
    }
  }

  return dist;
}

/// Dirties the osg vertices.
void PolyhedralPrimitive::dirty_vertex_visualization_data()
{
  _vert_array->dirty();
}

/// Gets the visualization vertex.
osg::Vec3* PolyhedralPrimitive::get_visualization_vertex(boost::shared_ptr<Polyhedron::Vertex> p) const
{
#ifdef USE_OSG
  assert(_vertex_mapping.find(p) != _vertex_mapping.end());
  return &(*_vert_array)[_vertex_mapping.find(p)->second];
#else
  return nullptr;
#endif
}

/// creates the visualization for the primitive
osg::Node* PolyhedralPrimitive::create_visualization()
{
  #ifdef USE_OSG
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // create the inverse transformation (to backout the applied transform)
  osg::MatrixTransform* transform = new osg::MatrixTransform;
  osg::Matrixd tgt;

  // get the rotation matrix
  Matrix3d M = Quatd::invert(_F->q);

  // setup the rotation components of tgt
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = M(i,j);

  // invert the translation
  Origin3d ix = -M * _F->x;

  // setup the translation components of tgt
  for (unsigned i=X; i<= Z; i++)
    tgt(W,i) = ix[i];

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (double) 0.0;
  tgt(W,W) = (double) 1.0;

  // set the transform
  transform->setMatrix(tgt);

  // do the necessary OSG stuff
  osg::Geode* geode = new osg::Geode;
  osg::Geometry* geometry = new osg::Geometry;
  geometry->setDataVariance(osg::Object::DYNAMIC);
  geometry->setUseDisplayList(false);
  geometry->setUseVertexBufferObjects(false);
  geode->addDrawable(geometry);

  // create vertices and setup a mapping
  _vertex_mapping.clear();
  _vert_array = new osg::Vec3Array(_poly.get_vertices().size());
  const vector<shared_ptr<Polyhedron::Vertex> >& v = _poly.get_vertices();
  for (unsigned i=0; i< v.size(); i++)
  {
    (*_vert_array)[i] = osg::Vec3(v[i]->o[X], v[i]->o[Y], v[i]->o[Z]);
    _vertex_mapping[v[i]] = i;
  }

  // associate the vertex array with the geometry
  geometry->setVertexArray(_vert_array);

  // iterate over all faces
  const vector<shared_ptr<Polyhedron::Face> >& f = _poly.get_faces();
  for (unsigned i=0; i< f.size(); i++)
  {
    // prepare to setup the indices
    osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON, 0);

    // iterate over the vertices in the face
    Polyhedron::VertexFaceIterator vfi(f[i], true);
    face->push_back(_vertex_mapping[*vfi]);
    while (vfi.has_next())
    {
      vfi.advance();
      face->push_back(_vertex_mapping[*vfi]);
    }  

    // add the polygon to the face    
    geometry->addPrimitiveSet(face);
  }

  transform->addChild(geode);
  return transform;
  #else
  return NULL;
  #endif
}

BVPtr PolyhedralPrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  // first get the vertices for the polyhedron
  shared_ptr<const Pose3d> geom_pose = get_pose(geom);
  vector<Point3d> vertices;
  get_vertices(geom_pose, vertices);

  // create the bounding box
  shared_ptr<AABB> bb(new AABB(vertices.begin(), vertices.end()));
  return bb;
}

void PolyhedralPrimitive::get_vertices(shared_ptr<const Pose3d> P, std::vector<Point3d>& vertices) const
{
  // clear the set of vertices
  vertices.clear();

  // verify that the pose is found
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // iterate through all vertices of the polyhedron
  const std::vector<shared_ptr<Polyhedron::Vertex> >& v = _poly.get_vertices();
  for (unsigned i=0; i< v.size(); i++)
    vertices.push_back(Point3d(v[i]->o, P));
}

shared_ptr<const IndexedTriArray> PolyhedralPrimitive::get_mesh(boost::shared_ptr<const Pose3d> P)
{
  return _mesh; 
}

/// Sets the pose of this primitive
void PolyhedralPrimitive::set_pose(const Pose3d& P)
{
  // setup P as a pointer so that we can compute a relative pose
  shared_ptr<Pose3d> pP(new Pose3d(P));

  // get the transformation between the last and this pose
  Transform3d T = Pose3d::calc_relative_pose(_F, pP);

  // call the primitive function
  Primitive::set_pose(P);

  // transform the polyhedron
  _poly = _poly.transform(T);
} 

/// Sets the polyhedron corresponding to this primitive
/**
 * Should only be done when the primitive hasn't been transformed.
 */
void PolyhedralPrimitive::set_polyhedron(const Polyhedron& p)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // verify no transformation
  Matrix3d R = _F->q;
  double trace = R(X,X) + R(Y,Y) + R(Z,Z);
  if (_F->x.norm() > NEAR_ZERO || std::fabs(trace - 3.0) > NEAR_ZERO)
    throw std::runtime_error("set_polyhedron() should only be called with identity pose");

  // set the polyhedron
  _poly = p;

  // get the mesh
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_poly.get_mesh()));

  // calculate mass properties
  calc_mass_properties();
}

/// Calculates mass properties for the polyhedron
void PolyhedralPrimitive::calc_mass_properties()
{
  std::cerr << "Polyhedral::calc_mass_properties() - implement me!" << std::endl;
}

/// Computes the distance between two polyhedra
double PolyhedralPrimitive::calc_dist(shared_ptr<const PolyhedralPrimitive> p, Point3d& pthis, Point3d& pp) const
{
  const double INF = std::numeric_limits<double>::max();
  shared_ptr<TessellatedPolyhedron> tpoly;

  // verify that p is convex
  if (!is_convex() && !p->is_convex())
    throw std::runtime_error("Polyhedron is not convex!");

  // if the primitive is polyhedral and convex, can use vclip 
  shared_ptr<const PolyhedralPrimitive> bthis = dynamic_pointer_cast<const PolyhedralPrimitive>(shared_from_this());
  shared_ptr<const Pose3d> poseA = pthis.pose;
  shared_ptr<const Pose3d> poseB = pp.pose;
  shared_ptr<const Polyhedron::Feature> closestA, closestB;

  // attempt to use vclip
  double dist = Polyhedron::vclip(bthis, p, poseA, poseB, closestA, closestB); 
  return std::max(dist, 0.0);
}

/// Computes the signed distance between two polyhedra
double PolyhedralPrimitive::calc_signed_dist(shared_ptr<const PolyhedralPrimitive> p, Point3d& pthis, Point3d& pp) const
{
  const double INF = std::numeric_limits<double>::max();
  shared_ptr<TessellatedPolyhedron> tpoly;

  // verify that p is convex
  if (!is_convex() && !p->is_convex())
    throw std::runtime_error("Polyhedron is not convex!");

  // if the primitive is polyhedral and convex, can use vclip 
  shared_ptr<const PolyhedralPrimitive> bthis = dynamic_pointer_cast<const PolyhedralPrimitive>(shared_from_this());
  shared_ptr<const Pose3d> poseA = pthis.pose;
  shared_ptr<const Pose3d> poseB = pp.pose;
  shared_ptr<const Polyhedron::Feature> closestA, closestB;

  // attempt to use vclip
  double dist = Polyhedron::vclip(bthis, p, poseA, poseB, closestA, closestB); 
  if (dist >= 0.0)
    return dist;

  // get the two polyhedra
  const Polyhedron& polyA = bthis->get_polyhedron();
  const Polyhedron& polyB = p->get_polyhedron();

  // get transforms to global frame
  Ravelin::Transform3d wTa = Ravelin::Pose3d::calc_relative_pose(poseA, GLOBAL);
  Ravelin::Transform3d wTb = Ravelin::Pose3d::calc_relative_pose(poseB, GLOBAL);

  // Compute first set of testing vectors from face normals
  std::vector<Ravelin::Vector3d> test_vectors;
  bthis->add_to_face_vector(wTa, test_vectors);
  p->add_to_face_vector(wTb, test_vectors);

  // create testing axes (cross-products of edges from A and B)
  std::vector <Ravelin::Vector3d> evA, evB;
  bthis->create_edge_vector(wTa, evA);
  p->create_edge_vector(wTb, evB);
  for (std::vector<Ravelin::Vector3d>::iterator evAi = evA.begin();
       evAi != evA.end(); ++evAi) {
    for (std::vector<Ravelin::Vector3d>::iterator evBi = evB.begin();
         evBi != evB.end(); ++evBi) {
      Ravelin::Vector3d xv = Ravelin::Vector3d::cross(*evAi, *evBi);
      double nrm = xv.norm();
      if (nrm > NEAR_ZERO) {
        xv /= nrm;
        test_vectors.push_back(xv);
      }
    }
  }

  // We want to find parallel and anti-parallel normalized vectors. The
  // algorithm works like this: 
  // Sort based on the magnitude of the first dimension of the vector.
  // Check the distance of all vectors within distance TOL of the given
  // dimension.
  const unsigned X = 0; 
  const double TOL = 1e-2;  // large tolerance
  std::vector<unsigned> vecs(test_vectors.size());
  for (unsigned i=0; i< test_vectors.size(); i++)
    vecs[i] = i;
  test_vecs = &test_vectors;
  std::sort(vecs.begin(), vecs.end(), compare_vecs);

  unsigned end = vecs.size()-1;
  for (int i=end-1; i>= 0; i--)
  {
    // look for the first vector that is farther than than tolerance from
    // end's element
    if (std::fabs(test_vectors[end][X]) - std::fabs(test_vectors[i][X]) > TOL)
    {
      // compare all of the vectors up to, but not including, i 
      for (int j=end-1; j > i; j--)
      {
        if (std::fabs(std::fabs(test_vectors[end].dot(test_vectors[j]) - 1.0)) < 1e-4)
        {
          test_vectors.erase(test_vectors.begin()+j);
          end--;  // we must now decrease end and i 
          i--;  // 
        }
      }

      // now must update end
      end--; 
    }
  }

  // create Vector3d for all vertices
  std::vector <boost::shared_ptr<Polyhedron::Vertex> > vAa = polyA.get_vertices();
  std::vector <boost::shared_ptr<Polyhedron::Vertex> > vBb = polyB.get_vertices();

  std::vector <Ravelin::Vector3d> vector_a;
  BOOST_FOREACH(boost::shared_ptr < Polyhedron::Vertex > vertex, vAa)
  {
    Ravelin::Vector3d v(vertex->o, wTa.source);
    Ravelin::Vector3d vw = wTa.transform_point(v);
    vector_a.push_back(vw);
  }
  std::vector <Ravelin::Vector3d> vector_b;
  BOOST_FOREACH(boost::shared_ptr < Polyhedron::Vertex > vertex, vBb)
  {
    Ravelin::Vector3d v(vertex->o, wTb.source);
    Ravelin::Vector3d vw = wTb.transform_point(v);
    vector_b.push_back(vw);
  }

  // ***********************************************************************
  // find the minimum overlap
  // ***********************************************************************
  double min_overlap = std::numeric_limits<double>::max();
  Ravelin::Vector3d min_axis;
  boost::shared_ptr <Polyhedron::Vertex> a_vertex, b_vertex;
  int direction = 1;
  for (std::vector<Ravelin::Vector3d>::iterator test_i = test_vectors.begin(); test_i != test_vectors.end(); ++test_i) {
    double min_a, max_a, min_b, max_b;
    boost::shared_ptr <Polyhedron::Vertex> minV_a, minV_b, maxV_a, maxV_b;
    unsigned min_index_a, min_index_b, max_index_a, max_index_b;

    // project vertices onto the candidate axis
    // NOTE: this could be done in O(lg n) time of the number of vertices rather
    // than the current O(n) time operation.
    project(vector_a, *test_i, min_a, max_a, min_index_a, max_index_a);
    project(vector_b, *test_i, min_b, max_b, min_index_b, max_index_b);

    // Compute the amount of overlap.
    double o1 = max_a - min_b;
    double o2 = max_b - min_a;

    if (o1 > 0.0 && o2 > 0.0) {
      // there is an overlap
      double overlap = std::min(o1, o2);
      boost::shared_ptr <Polyhedron::Vertex> v1, v2;

      if (min_overlap - overlap > NEAR_ZERO) {
        min_overlap = overlap;
        min_axis = *test_i;
        if (fabs(overlap - o1) < NEAR_ZERO) {
          direction = -1;
          a_vertex = vAa[max_index_a];
          assert(a_vertex);
          b_vertex = vBb[min_index_b];
          assert(b_vertex);
        } else {
          direction = 1;
          b_vertex = vBb[max_index_b];
          assert(b_vertex);
          a_vertex = vAa[min_index_a];
          assert(a_vertex);
        }
      }
    }
  }

  // ensure that the distance is negated
  assert(min_overlap >= 0.0);
  dist = -min_overlap;

  return dist;
/*
  // compute transforms
  Transform3d wTa = Pose3d::calc_relative_pose(poseA, GLOBAL);
  Transform3d wTb = Pose3d::calc_relative_pose(poseB, GLOBAL);

  // compute half-spaces
  vector<std::pair<Vector3d, double> > hs;
  PolyhedralPrimitive::get_halfspaces(bthis->get_polyhedron(), poseA, wTa, std::back_inserter(hs));
  PolyhedralPrimitive::get_halfspaces(p->get_polyhedron(), poseB, wTb, std::back_inserter(hs));

  // find a point interior to both
  Ravelin::Origin3d ip;
  double dist2 = -CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);
  if (dist2 >= NEAR_ZERO)
    return dist2;
  assert(dist2 < NEAR_ZERO);

  // attempt to calculate the half-space intersection
  try
  {
    tpoly = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
  }
  catch (NumericalException e)
  {
    // if we're here, then the volume of intersection is an area, find
    // closest points and quit
    Vector3d ip0(ip, GLOBAL);
    pthis = wTa.inverse_transform_point(ip0);
    pp = wTb.inverse_transform_point(ip0);

    return 0.0;
  } 

 // setup minimum distance and the contact plane offset
  double min_dist = INF;

  // get the vertices from the polyhedron
  const std::vector<Ravelin::Origin3d>& vertices = tpoly->get_vertices();

  // determine the signed distance 
  const std::vector<IndexedTri>& facets = tpoly->get_facets();
  for (unsigned i=0; i< facets.size(); i++)
  {
    // setup a triangle
    Triangle tri(Point3d(vertices[facets[i].a], GLOBAL), 
                 Point3d(vertices[facets[i].b], GLOBAL),
                 Point3d(vertices[facets[i].c], GLOBAL));

    // get the reverse of the normal and the offset
    Ravelin::Vector3d ncand = -tri.calc_normal();
    double offset = tri.calc_offset(ncand);

    // NOTE: this is currently an O(n) operation, but it could be turned into
    //       an O(lg N) one
    Point3d p(tpoly->find_extreme_vertex(Ravelin::Origin3d(ncand)), GLOBAL);

    // compute the distance of the extremal point from the face
    double dist = ncand.dot(p) - offset;
    assert(dist > -NEAR_ZERO);
    if (dist < min_dist)
    {
      // if the absolute distance is less than the minimum, we've found our
      // normal
      min_dist = dist;
    }
  }

  // return the negative minimum distance; don't worry about closest points
  // during interpenetration
  return -min_dist;
  */
}

/// Computes the distance between this primitive and another
double PolyhedralPrimitive::calc_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // now try polyhedron/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
  {
    throw std::runtime_error("Sphere/polyhedron distance not yet implemented");
    return 0.0;
  } 

  // now try polyhedron/plane
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return planep->calc_signed_dist(bthis, pp, pthis);
  }

  // now try polyhedron/heightmap
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
  {
    throw std::runtime_error("Heightmap/polyhedron distance not yet implemented");
    return 0.0;
  }

  // now try convex polyhedron/convex polyhedron
  shared_ptr<const PolyhedralPrimitive> polyp = dynamic_pointer_cast<const PolyhedralPrimitive>(p);
  if (polyp)
    return calc_dist(polyp, pthis, pp);

  throw std::runtime_error("Unanticipated signed distance types!");
  return 0.0;
}

/// Computes the signed distance between this primitive and another
double PolyhedralPrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // now try polyhedron/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
  {
    throw std::runtime_error("Sphere/polyhedron distance not yet implemented");
    return 0.0;
  } 

  // now try polyhedron/plane
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return planep->calc_signed_dist(bthis, pp, pthis);
  }

  // now try polyhedron/heightmap
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
  {
    throw std::runtime_error("Heightmap/polyhedron distance not yet implemented");
    return 0.0;
  }

  // now try convex polyhedron/convex polyhedron
  shared_ptr<const PolyhedralPrimitive> polyp = dynamic_pointer_cast<const PolyhedralPrimitive>(p);
  if (polyp)
    return calc_signed_dist(polyp, pthis, pp);

  throw std::runtime_error("Unanticipated signed distance types!");
  return 0.0;
}

/// Implements Base::load_from_xml() for serialization
void PolyhedralPrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is PolyhedralPrimitive
  // DO NOT DO THIS B/C OF CHILDREN
  // assert(strcasecmp(node->name.c_str(), "Polyhedron") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // make sure that a filename is specified
  XMLAttrib* fname_attr = node->get_attrib("filename");
  if (!fname_attr)
  {
    cerr << "PolyhedralPrimitive::load_from_xml() - trying to load a ";
    cerr << " polyhedron w/o a filename!" << endl;
    cerr << "  offending node: " << endl << *node << endl;
    return;
  }

  // setup the file extensions
  const char* OBJ_EXT = ".obj";

  // get the filename
  string fname(fname_attr->get_string_value());

  // get the lowercase version of the filename
  string fname_lower = fname;
  std::transform(fname_lower.begin(), fname_lower.end(), fname_lower.begin(), ::tolower);

  // get the type of file and construct the triangle mesh appropriately
  if (fname_lower.find(string(OBJ_EXT)) == fname_lower.size() - strlen(OBJ_EXT))
  {
    // setup a transform
    Transform3d T = Pose3d::calc_relative_pose(_F, GLOBAL);

    // read in the file using an indexed triangle array
    IndexedTriArray ita = IndexedTriArray::read_from_obj(fname).transform(T);

    // get all of the vertices and compute the convex hull (yielding a
    // tessellated polyhedron)
    const std::vector<Origin3d>& vertices = ita.get_vertices();
    TessellatedPolyhedronPtr tessellated_poly = CompGeom::calc_convex_hull(vertices.begin(), vertices.end());   

    // convert the tessellated polyhedron to a standard polyhedron and set it
    // NOTE: we avoid the set function b/c a transform may have been applied
    tessellated_poly->to_polyhedron(_poly);
  }
  else
  {
    cerr << "PolyhedralPrimitive::load_from_xml() - unrecognized filename extension" << endl;
    cerr << "  for attribute 'filename'.  Valid extensions are '.obj' (Wavefront OBJ)" << endl;
  }

  // update the visualization
  update_visualization();
 
  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void PolyhedralPrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Polyhedron";

  // make a filename using "this"
  const unsigned MAX_DIGITS = 28;
  char buffer[MAX_DIGITS+1];
  sprintf(buffer, "%p", this);
  string filename = "polyhedron" + string(buffer) + ".obj";

  // add the filename as an attribute
  node->attribs.insert(XMLAttrib("filename", filename));

  // do not save the array to the OBJ file if it already exists (which we
  // crudely check for using std::ifstream to avoid OS-specific calls -- note
  // that it is possible that opening a file may fails for other reasons than
  // the file does not exist)
  std::ifstream in(filename.c_str());
  if (in.fail())
  {
    // get the transform from the global pose to the primitive pose
    shared_ptr<const Pose3d> P = get_pose();
    Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);
    Polyhedron poly_xform = _poly.transform(T);

    // write the mesh
    poly_xform.write_to_obj(filename);
  }
  else
    in.close();
}

void PolyhedralPrimitive::project(const std::vector<Ravelin::Vector3d>& vectors, const Ravelin::Vector3d& axis, double& min_dot, double& max_dot, unsigned& min_index, unsigned& max_index)
{
  // initialize outputs to indicator and safe values
  min_index = max_index = std::numeric_limits<unsigned>::max();
  min_dot = std::numeric_limits<double>::max();
  max_dot = -std::numeric_limits<double>::max();

  for(unsigned i = 0 ; i < vectors.size(); ++i)
  {
    double value = axis.dot(vectors[i]);
    if (value < min_dot)
    {
      min_dot = value;
      min_index = i;
    }

    if (value > max_dot)
    {
      max_dot = value;
      max_index = i;
    }
  }
}

