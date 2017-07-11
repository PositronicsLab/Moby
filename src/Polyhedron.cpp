/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <map>
#include <Moby/Types.h>
#include <Moby/Log.h>
#include <Moby/CompGeom.h>
#include <Moby/PolyhedralPrimitive.h>
#include <Moby/Polyhedron.h>

using namespace Ravelin;
using namespace Moby;
using boost::weak_ptr;
using boost::shared_ptr;
using boost::shared_array;
using boost::const_pointer_cast; 
using boost::dynamic_pointer_cast; 
using std::map;
using std::vector;
using std::make_pair;
using std::endl;
using std::list;

#define CREATE_LOOKUP(vA, vB, eAB) { \
if ((vei = v_edges.find(std::make_pair(vA, vB))) != v_edges.end()) \
{ \
  eAB = vei->second; \
  if (cw) \
  { \
    assert(!eAB->face1); \
    eAB->face1 = f; \
  } \
  else \
  { \
    assert(!eAB->face2); \
    eAB->face2 = f; \
  } \
} \
else if ((vei = v_edges.find(std::make_pair(vB, vA))) != v_edges.end()) \
{ \
  eAB = vei->second; \
  if (cw) \
  { \
    assert(!eAB->face2); \
    eAB->face2 = f; \
  } \
  else \
  { \
    assert(!eAB->face1); \
    eAB->face1 = f; \
  } \
} \
else \
{ \
  eAB = boost::shared_ptr<Polyhedron::Edge>(new Polyhedron::Edge); \
  v_edges[std::make_pair(vA, vB)] = eAB; \
  if (cw) \
    eAB->face1 = f; \
  else \
    eAB->face2 = f; \
  eAB->v1 = vertex_map[vA->point]; \
  eAB->v2 = vertex_map[vB->point]; \
  poly._edges.push_back(eAB); \
  eAB->v1->e.push_back(eAB); \
  eAB->v2->e.push_back(eAB); \
} \
}

/// Sums coincident normals to the vertex on a polyhedron.
Origin3d Polyhedron::Vertex::sum_coincident_normals() const
{
  // Get all coincident faces.
  std::set<shared_ptr<Polyhedron::Face>> coincident_faces;
  for (auto weak_edge : this->e)
  {
    shared_ptr<Polyhedron::Edge> edge(weak_edge);
    coincident_faces.insert(edge->face1);
    coincident_faces.insert(edge->face2);
  }

  // Sum normals.
  Vector3d sum(0, 0, 0, GLOBAL);
  for (auto f : coincident_faces)
    sum += f->get_plane().get_normal();
  return Origin3d(sum);
}

/// Gets the polyhedron as a triangle mesh
IndexedTriArray Polyhedron::get_mesh() const
{
  const unsigned TRI_EDGES = 3, QUAD_EDGES = 4;

  // first make sure that all faces are triangles
  #ifndef NDEBUG
  for (unsigned i=0; i< _faces.size(); i++)
    assert(_faces[i]->e.size() == TRI_EDGES || 
           _faces[i]->e.size() == QUAD_EDGES);
  #endif

  // setup vertices and mapping to vertices
  std::map<shared_ptr<Vertex>, unsigned> vertex_mapping;
  std::vector<Origin3d> vertices;
  for (unsigned i=0; i< _vertices.size(); i++)
  {
    vertices.push_back(_vertices[i]->o);
    vertex_mapping[_vertices[i]] = i;
  }

  // iterate over all faces
  std::vector<IndexedTri> faces;
  for (unsigned i=0; i< _faces.size(); i++)
  {
    // iterate over all vertices in the face
    VertexFaceIterator vfi(_faces[i], true);
    IndexedTri it;
    it.a = vertex_mapping[*vfi];
    assert(vfi.has_next());
    vfi.advance();
    it.b = vertex_mapping[*vfi];
    assert(vfi.has_next());
    vfi.advance();
    it.c = vertex_mapping[*vfi];
    faces.push_back(it);
    if (vfi.has_next())
    {
      IndexedTri it2;
      it2.a = it.c; 
      vfi.advance();
      it2.b = vertex_mapping[*vfi];
      it2.c = it.a;
      faces.push_back(it2);
    }
  }

  return IndexedTriArray(vertices.begin(), vertices.end(), faces.begin(), faces.end()); 
}

/// Determines whether the polyhedron is convex
void Polyhedron::determine_convexity()
{
  // set convexity to -inf to begin
  _convexity = -std::numeric_limits<double>::max();

  // iterate over each face
  for (unsigned i=0; i< _faces.size(); i++)
  {
    // get the plane for the face
    Plane p = _faces[i]->get_plane();

    // get the pose for the plane
    shared_ptr<const Pose3d> pose = p.get_pose();

    // get the maximum signed distance from any vertex to the face 
    for (unsigned j=0; j< _vertices.size(); j++)
      _convexity = std::max(_convexity, p.calc_signed_distance(Point3d(_vertices[j]->o, pose)));
  }

  // indicate that convexity has been calculated
  _convexity_computed = true;
}

/// Gets the plane containing a face
Plane Polyhedron::Face::get_plane() const
{
  if (e.size() < 3)
    throw std::runtime_error("Tried to get plane containing face when face has fewer than 3 edges!"); 

  // get this as non-constant
  boost::shared_ptr<const Polyhedron::Face> fthis_const(shared_from_this());
  boost::shared_ptr<Polyhedron::Face> fthis = boost::const_pointer_cast<Polyhedron::Face>(fthis_const);

  // iterate over vertices
  Polyhedron::VertexFaceIterator vfi(fthis, true);

  // need three vertices - get first two arbitrarily
  shared_ptr<Polyhedron::Vertex> v1 = *vfi;
  assert(vfi.has_next());
  vfi.advance();
  shared_ptr<Polyhedron::Vertex> v2 = *vfi;

  // now iterate through the remainder
  Origin3d normal(0.0, 0.0, 0.0);
  do 
  {
    // make sure that there is another vertex
    assert(vfi.has_next());

    // get the next vertex
    vfi.advance();
    shared_ptr<Polyhedron::Vertex> v3 = *vfi;

    // see whether the three make a non-zero cross product
    if ((normal = Origin3d::cross(v2->o - v1->o, v3->o - v2->o)).norm() > NEAR_ZERO)
      break;
  }
  while (vfi.has_next());

  // verify that the normal norm is ok
  if (normal.norm() > NEAR_ZERO)
    normal.normalize();
  {
    // redetermine the normal, now using the biggest norm
    double biggest_nrm = 0.0;
    normal = Origin3d(0.0, 0.0, 0.0);
    vfi = Polyhedron::VertexFaceIterator(fthis, true);
    do 
    {
      // make sure that there is another vertex
      assert(vfi.has_next());

      // get the next vertex
      vfi.advance();
      shared_ptr<Polyhedron::Vertex> v3 = *vfi;
  
      // see whether the three make a non-zero cross product
      Origin3d cand_normal = Origin3d::cross(v2->o - v1->o, v3->o - v2->o);
      if (cand_normal.norm() > biggest_nrm)
      {
        normal = cand_normal;
        biggest_nrm = normal.norm();
      }
    }
    while (vfi.has_next());

    // normalize the normal
    normal /= biggest_nrm;
  }

  // compute d
  double d = normal.dot(v1->o);

  return Plane(Vector3d(normal, GLOBAL), d);
}

/// Creates a minimum polyhedron
Polyhedron::Polyhedron()
{
  // indicate convexity has not been computed
  _convexity_computed = false;
}

/// Creates a polyhedron from the given vectors.
Polyhedron::Polyhedron(const vector<shared_ptr<Vertex>>& v,
                       const vector<shared_ptr<Edge>>& e,
                       const vector<shared_ptr<Face>>& f) {
  // Make copies.
  _vertices = v;
  _edges = e;
  _faces = f;

  // Indicate convexity has not been computed.
  _convexity_computed = false;

  // Compute the bounding box.
  calc_bounding_box();
}

/// Assignment operator
Polyhedron& Polyhedron::operator=(const Polyhedron& p)
{
  // copy easy things
  _bb_min = p._bb_min;
  _bb_max = p._bb_max;
  _convexity = p._convexity;
  _convexity_computed = p._convexity_computed;

  // clear vertices, edges, and faces
  _vertices.clear();
  _edges.clear();
  _faces.clear();

  // make a copy of vertices 
  std::map<shared_ptr<Polyhedron::Vertex>, shared_ptr<Polyhedron::Vertex> > v_copy;
  for (unsigned i=0; i< p._vertices.size(); i++)
  {
    shared_ptr<Polyhedron::Vertex> v(new Polyhedron::Vertex);
    v->o = p._vertices[i]->o;
    v_copy[p._vertices[i]] = v;
    v->data = p._vertices[i]->data;
    _vertices.push_back(v);
  }

  // make a copy of edges
  std::map<shared_ptr<Polyhedron::Edge>, shared_ptr<Polyhedron::Edge> > e_copy;
  for (unsigned i=0; i< p._edges.size(); i++)
  {
    shared_ptr<Polyhedron::Edge> e(new Polyhedron::Edge);
    e->v1 = v_copy[p._edges[i]->v1];
    e->v2 = v_copy[p._edges[i]->v2];
    e_copy[p._edges[i]] = e;
    e->data = p._edges[i]->data;
    _edges.push_back(e);
  }

  // make a copy of faces
  std::map<shared_ptr<Polyhedron::Face>, shared_ptr<Polyhedron::Face> > f_copy;
  for (unsigned i=0; i< p._faces.size(); i++)
  {
    shared_ptr<Polyhedron::Face> f(new Polyhedron::Face);
    BOOST_FOREACH(weak_ptr<Edge> e, p._faces[i]->e)
      f->e.push_back(e_copy[shared_ptr<Edge>(e)]);
    f->data = p._faces[i]->data;
    f_copy[p._faces[i]] = f;
    _faces.push_back(f);
  }

  // update edge pointers for vertices
  for (unsigned i=0; i< _vertices.size(); i++)
    BOOST_FOREACH(weak_ptr<Edge> e, p._vertices[i]->e)
      _vertices[i]->e.push_back(e_copy[shared_ptr<Edge>(e)]);

  // update face pointers for edges
  for (unsigned i=0; i< _edges.size(); i++)
  {
    _edges[i]->face1 = f_copy[p._edges[i]->face1];
    _edges[i]->face2 = f_copy[p._edges[i]->face2];
  }

  return *this;
}

/// Does a shallow copy of this polyhedron
Polyhedron Polyhedron::shallow_copy() const
{
  Polyhedron p;
  p._bb_min = _bb_min;
  p._bb_max = _bb_max;
  p._convexity = _convexity;
  p._convexity_computed = _convexity_computed;
  p._edges = _edges;
  p._vertices = _vertices;
  p._faces = _faces;

  return p;
}

/// Transforms a polyhedron
Polyhedron Polyhedron::transform(const Transform3d& T) const
{
  // copy this
  Polyhedron p = *this;

  // get the vertices
  std::vector<shared_ptr<Vertex> >& v = p._vertices;

  // transform the vertices
  for (unsigned i=0; i< v.size(); i++)
    v[i]->o = Origin3d(T.transform_point(Point3d(v[i]->o, T.source)));

  return p;    
}

/// Writes the polyhedron to Wavefront OBJ format 
void Polyhedron::write_to_obj(const std::string& filename) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // open the OBJ file for writing
  std::ofstream out(filename.c_str());

  // create a map from all vertices to indices, writing vertices at same time
  std::map<shared_ptr<Polyhedron::Vertex>, unsigned> vmap;
  for (unsigned i=0; i< _vertices.size(); i++)
  {
    out << "v " << _vertices[i]->o[X] << " " << _vertices[i]->o[Y] << " " << _vertices[i]->o[Z] << std::endl; 
    vmap[_vertices[i]] = i+1;
  }

  // iterate through each face 
  for (unsigned i=0; i< _faces.size(); i++)
  {
    out << "f";
    VertexFaceIterator vfi(_faces[i], true);
    assert(vmap.find(*vfi) != vmap.end());
    out << " " << vmap[*vfi];
    assert(vfi.has_next());
    vfi.advance(); 
    assert(vmap.find(*vfi) != vmap.end());
    out << " " << vmap[*vfi];
    assert(vfi.has_next());
    vfi.advance(); 
    assert(vmap.find(*vfi) != vmap.end());
    out << " " << vmap[*vfi];
    while (vfi.has_next())
    {
      vfi.advance();
      assert(vmap.find(*vfi) != vmap.end());
      out << " " << vmap[*vfi];
    }
    out << std::endl;
  }

  out.close();
}

/// Computes the Minkowski difference of two polyhedral primitives
/**
 * \return a polyhedron- the 'data' field of each vertex is of type 
 *         boost::shared_ptr<std::pair<int, int> >, where the first integer
 *         is the index of a vertex from the first polyhedron and the second
 *         integer is the index of a vertex from the second polyhedron.
 */
Polyhedron Polyhedron::calc_minkowski_diff(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const int DIM = 3;

  // get the vertices of A
  vector<Point3d> vA;
  pA->get_vertices(poseA, vA);

  // get the vertices of B
  vector<Point3d> vB;
  pB->get_vertices(poseB, vB);

  // convert A's vector of vertices to the global frame
  Transform3d wTposeA = Pose3d::calc_relative_pose(poseA, GLOBAL);
  for (unsigned i=0; i< vA.size(); i++)
    vA[i] = wTposeA.transform_point(vA[i]);

  // transform poses to global frame
  Transform3d wTposeB = Pose3d::calc_relative_pose(poseB, GLOBAL);
  for (unsigned i=0; i< vB.size(); i++)
    vB[i] = wTposeB.transform_point(vB[i]);

  // subtract vertices of B from vertices of A in the global frame, saving
  // original vertices that they came from
  const unsigned NVERTS = vA.size() * vB.size();
  vector<shared_ptr<Polyhedron::Vertex> > verts(NVERTS);
  for (unsigned i=0, k=0; i< vA.size(); i++)
    for (unsigned j=0; j< vB.size(); j++, k++)
    {
      verts[k] = shared_ptr<Polyhedron::Vertex>(new Polyhedron::Vertex);
      verts[k]->o = Origin3d(vA[i]) - Origin3d(vB[j]); 
      shared_ptr<std::pair<int, int> > intpair(new std::pair<int, int>(i, j));
      verts[k]->data = intpair;
    }

  // setup mapping for qhull
  map<coordT*, shared_ptr<Polyhedron::Vertex> > vertex_map;
  vector<coordT> qhull_points;
  qhull_points.resize(NVERTS*DIM);
  coordT* points_begin = &qhull_points.front();
  unsigned j=0;
  for (unsigned i=0, j=0; i< NVERTS; i++)
  {
    vertex_map[points_begin+j] = verts[i];
    qhull_points[j++] = verts[i]->o[X];
    qhull_points[j++] = verts[i]->o[Y];
    qhull_points[j++] = verts[i]->o[Z];
  }

  // ******************* prepare to calculate convex hull *******************
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";  // TODO: remove Fx option?
  FILE* outfile, * errfile;

  // setup constants for qhull
  const int N_POINTS = NVERTS; 
  const boolT IS_MALLOC = false;
  if (N_POINTS < 4)
    return Polyhedron();

  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

  // construct the convex hull
  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code)
  {
    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw NumericalException(); 
    return Polyhedron();
  }

  // create a new Polyhedron
  Polyhedron poly;

  // iterate through all vertices
  vertexT* vertex;
  FORALLvertices
  {
    boost::shared_ptr<Polyhedron::Vertex> v = vertex_map[vertex->point];
    poly._vertices.push_back(v);
  }

  if (LOGGING(LOG_COMPGEOM))
  {
    // triangulation may or may not be disabled...
    qh_triangulate();
    for (unsigned i=0; i< poly._vertices.size(); i++)
      FILE_LOG(LOG_COMPGEOM) << "vertex " << i << ": " << poly._vertices[i]->o << std::endl;
  }

  // triangulate
  qh_triangulate();

  // need maps for new edges created for simplicial and non-simplicial facets
  std::map<std::pair<vertexT*, vertexT*>, boost::shared_ptr<Polyhedron::Edge> > v_edges;
  std::map<ridgeT*, boost::shared_ptr<Polyhedron::Edge> > r_edges;

  // setup an iterator for simplicial facets
  std::map<std::pair<vertexT*, vertexT*>, boost::shared_ptr<Polyhedron::Edge> >::const_iterator vei;
  
  // iterate through facets
  for (facetT* facet=qh facet_list;facet && facet->next;facet=facet->next)
  {
    if (!facet->vertices)
      continue;

    // create a new facet
    boost::shared_ptr<Polyhedron::Face> f(new Polyhedron::Face);

    // see whether the facet is simplicial- it changes how we must process
    // edges
    if (facet->simplicial)
    {
      // see how the facet is oriented
      bool cw = (facet->toporient ^ qh_ORIENTclock);

      // edges will be between each vertex; get all vertices in the facet
      vertexT** vertex_pointer = (vertexT**)& ((facet->vertices)->e[0].p); 
      vertexT* v1 = *vertex_pointer++;
      vertexT* v2 = *vertex_pointer++;
      vertexT* v3 = *vertex_pointer;

      // setup three edges
      boost::shared_ptr<Polyhedron::Edge> e12, e23, e31;

      // create / lookup the three edges
      CREATE_LOOKUP(v1, v2, e12);
      CREATE_LOOKUP(v2, v3, e23);
      CREATE_LOOKUP(v3, v1, e31);

      // add all three edges to the face
      if (!cw)
      {
        f->e.push_back(e12);
        f->e.push_back(e23);
        f->e.push_back(e31);
      }
      else
      {
        f->e.push_back(e12);
        f->e.push_back(e31);
        f->e.push_back(e23);
      }
      assert(e12 != e23 && e12 != e31 && e23 != e31);
    }
    else
    {
      // setup lists of ccw vertices
      std::list<vertexT*> ccw_vertices;

      // facet is non-simplicial; iterate over the "ridges" (edges)
      ridgeT* ridge;    // for iterating...
      ridgeT** ridgep;  // ...over ridges
      FOREACHridge_( facet->ridges )
      {
        // setup the edge
        boost::shared_ptr<Polyhedron::Edge> e;

        // get whether the ridge is cw or ccw
        bool cw = ((ridge->top == facet) ^ qh_ORIENTclock);

        // see whether the ridge/edge is already in the map
        std::map<ridgeT*, boost::shared_ptr<Polyhedron::Edge> >::const_iterator new_edge_iter = r_edges.find(ridge);
        if (new_edge_iter == r_edges.end())
        {
          // create an edge and add it to the map
          e = boost::shared_ptr<Polyhedron::Edge>(new Polyhedron::Edge);
          poly._edges.push_back(e);
          r_edges[ridge] = e;

          // get the pointer to the vertices of the ridge
          vertexT** vertex_pointer = (vertexT**)& ((ridge->vertices)->e[0].p); 
          vertexT* vertex = *vertex_pointer;
          vertexT* v1 = *vertex_pointer;
          vertexT* v2 = *++vertex_pointer;

          // setup the vertices
          e->v1 = vertex_map[v1->point];
          e->v2 = vertex_map[v2->point];
          assert(e->v1 != e->v2);

          // add edge to the vertices
          e->v1->e.push_back(e);
          e->v2->e.push_back(e);

          // add the edge
          v_edges[std::make_pair(v1, v2)] = e;
          if (cw) 
            e->face1 = f;
          else
            e->face2 = f;
        }
        else
          e = new_edge_iter->second;

        // get the pointer to the vertices of the ridge
        vertexT** vertex_pointer = (vertexT**)& ((ridge->vertices)->e[0].p); 
        vertexT* vertex = *vertex_pointer;
        vertexT* v1 = *vertex_pointer;
        vertexT* v2 = *++vertex_pointer;

        // setup the vertices
        vertexT* mini_list[2];
        mini_list[0] = v1;
        mini_list[1] = v2;

        // add vertices to the list
        if (cw)
          std::swap(mini_list[0], mini_list[1]);
        ccw_vertices.push_back(mini_list[0]);
        ccw_vertices.push_back(mini_list[1]);

        // setup face- we'll use the convention that face1 is qhull's top
        // and face2 is qhull's bottom
        if (ridge->top == facet)
          e->face1 = f;
        else
          e->face2 = f;
      }

      // now setup the edge traversal 
      std::list<vertexT*>::const_iterator ccw_iter = ccw_vertices.begin();
      while (ccw_iter != ccw_vertices.end())
      {
        // get two vertices
        vertexT* v1 = *ccw_iter++;
        vertexT* v2 = *ccw_iter++;
 
        // setup the edge
        boost::shared_ptr<Polyhedron::Edge> e;

        // lookup the edge
        if ((vei = v_edges.find(std::make_pair(v1, v2))) != v_edges.end())
          e = vei->second;
        else if ((vei = v_edges.find(std::make_pair(v2, v1))) != v_edges.end())
          e = vei->second;
        else
          assert(false);

        // add it to the list
        f->e.push_back(e);
      }
    }

    // add the face to the polyhedron
    poly._faces.push_back(f);
  }

  // setup edge and vertex map for each face
  #ifndef NDEBUG
  std::map<boost::shared_ptr<Polyhedron::Edge>, unsigned> emap;
  for (unsigned i=0; i< poly._edges.size(); i++)
    emap[poly._edges[i]] = i;
  #endif

  // setup edge walk in each face
  for (unsigned i=0; i< poly._faces.size(); i++)
  {
    // get the edges in face
    boost::shared_ptr<Polyhedron::Face> f = poly._faces[i]; 
    FILE_LOG(LOG_COMPGEOM) << "examining face " << i << std::endl;

    // loop over all edges in the face 
    BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we, f->e)
    {
      // get the edge
      boost::shared_ptr<Polyhedron::Edge> e(we);
      #ifndef NDEBUG
      FILE_LOG(LOG_COMPGEOM) << "examining edge " << emap[e] << std::endl;
      #endif

      // loop over all other edges in the face
      BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we2, f->e)
      {
        // get the edge
        boost::shared_ptr<Polyhedron::Edge> e2(we2);

        // don't process same edge twice
        if (e == e2)
          continue;

        #ifndef NDEBUG
        FILE_LOG(LOG_COMPGEOM) << "  against edge " << emap[e2] << std::endl;
        #endif
      }
    }      
  }

  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  // mark the polyhedron as convex (here and in calc_convex_hull())
  poly._convexity_computed = true;
  poly._convexity = -1.0;

  // calculate the axis-aligned bounding box  
  poly.calc_bounding_box();

  return poly;
}

/// Calculates the bounding box
void Polyhedron::calc_bounding_box()
{
  const unsigned THREE_D = 3;
  
  if (_vertices.empty())
  {
    _bb_min.set_zero();
    _bb_max.set_zero();
  }
  else
  {
    // setup the bounding box
    _bb_min = _vertices.front()->o;
    _bb_max = _vertices.front()->o;
    for (unsigned i=1; i< _vertices.size(); i++)
      for (unsigned j=0; j< THREE_D; j++)
      {
        if (_vertices[i]->o[j] < _bb_min[j])
          _bb_min[j] = _vertices[i]->o[j];
        else if (_vertices[i]->o[j] > _bb_max[j])
          _bb_max[j] = _vertices[i]->o[j];
      }
  }
}

/// Finds the feature(s) of this polyhedron closest to the point
/**
 * \param p the query point
 * \param closest_features the closest features on return
 * \param inside whether the point is inside or outside the polyhedron on return
 * \return the distance
 */
/*
double Polyhedron::find_closest_features(const Origin3d& p, std::list<shared_ptr<Polyhedron::Feature> >& closest_features, bool& inside) const 
{
  boost::shared_ptr<const Pose2d> GLOBAL2D;
  std::set<shared_ptr<Polyhedron::Edge> > pedges;

  // setup closest features
  closest_features.clear();

  // setup closest distance
  double closest_dist = std::numeric_limits<double>::max();

  // set the point to inside initially
  inside = true;

  // iterate over all faces
  for (unsigned i=0; i< _faces.size(); i++)
  {
    // get the plane containing the face
    Plane plane = _faces[i]->get_plane();

    // get the distance of the point to the face
    Point3d pp(p, GLOBAL);
    double plane_dist = std::fabs(plane.calc_signed_distance(pp));

    // see whether the point is above the plane
    if (plane_dist > 0.0)
      inside = false;

    // if the distance is more than the closest distance, do not
    // check this facet further
    if (std::fabs(plane_dist) > closest_dist)
      continue; 

    // project the point onto the plane
    Ravelin::Vector2d p2d(plane.to_2D(pp), GLOBAL2D);

    // iterate over every edge in the ccw face-edge sequence; see whether 
    // point is always to the left side 
    bool strictly_inside = true;
    VertexFaceIterator vfi(_faces[i], true);
    shared_ptr<Polyhedron::Vertex> vi = *vfi;
    shared_ptr<Polyhedron::Vertex> vfirst = *vfi;
    do 
    {
      // get the next vertex
      assert(vfi.has_next());
      vfi.advance();
      shared_ptr<Polyhedron::Vertex> vj = *vfi;

      // project the two vertices of the edge to 2D
      Ravelin::Vector2d ea2d(plane.to_2D(Point3d(vi->o, GLOBAL)), GLOBAL2D);
      Ravelin::Vector2d eb2d(plane.to_2D(Point3d(vj->o, GLOBAL)), GLOBAL2D);

      // check whether point is strictly inside 
      CompGeom::OrientationType otype = CompGeom::area_sign(p2d, ea2d, eb2d);
      if (otype != CompGeom::eLeft)
      {
        strictly_inside = false;
        break;
      }

      // copy vj to vi
      vi = vj;
    }
    while (vfi.has_next());

    // do one more time to account for the last vertex
    if (strictly_inside)
    {
      shared_ptr<Polyhedron::Vertex> vj = vfirst;
    
      // project the two vertices of the edge to 2D
      Ravelin::Vector2d ea2d(plane.to_2D(Point3d(vi->o, GLOBAL)), GLOBAL2D);
      Ravelin::Vector2d eb2d(plane.to_2D(Point3d(vj->o, GLOBAL)), GLOBAL2D);

      // check whether point is strictly inside 
      CompGeom::OrientationType otype = CompGeom::area_sign(p2d, ea2d, eb2d);
      if (otype != CompGeom::eLeft)
        strictly_inside = false;
    }

    // if the projected point is strictly inside, then the closest feature
    // adjacent to this face is the interior of the face and the distance is
    // projected distance to the face
    if (strictly_inside)
    {
      // record closest feature
      if (std::fabs(plane_dist) < closest_dist)
      {
        closest_features.clear();
        closest_dist = std::fabs(plane_dist);
      }
      closest_features.push_back(_faces[i]);
    } 
    // proj. point on border or outside polygon; check distance from each edge 
    else
    {
      // iterate over all edges
      shared_ptr<Polyhedron::Edge> closest_edge;
      EdgeFaceIterator efi(_faces[i]);
      while (true) 
      {
        // get the next edge
        shared_ptr<Polyhedron::Edge> e = *efi;

        // don't process edge twice
        if (pedges.find(e) != pedges.end())
        {
          if (!efi.has_ccw())
            break;
          efi.advance_ccw(); 
          continue;
        }
        else
          pedges.insert(e);

        // compute the distance between the edge and the point 
        double t;
        LineSeg3 seg;
        Point3d closest;
        seg.first = Point3d(e->v1->o, GLOBAL);
        seg.second = Point3d(e->v2->o, GLOBAL);
        double dist = CompGeom::calc_dist(seg, pp, t, closest);

        // check to see whether the features should be added
        bool add = (dist <= closest_dist);

        if (dist < closest_dist)
        {
          closest_dist = dist;
          closest_features.clear();
        }
        if (add)
        {
          if (t < NEAR_ZERO)
            closest_features.push_back(e->v2);
          else if (std::fabs(t-1.0) < NEAR_ZERO)
            closest_features.push_back(e->v1);
          else
            closest_features.push_back(e);
        }

        // see whether to keep iterating
        if (!efi.has_ccw())
          break;
        efi.advance_ccw(); 
      }
    }
  }

  // finally, make the list of closest features unique
  closest_features.sort();
  closest_features.erase(std::unique(closest_features.begin(), closest_features.end()), closest_features.end());
  return closest_dist;
}
*/
/*
/// Finds the feature adjacent to this edge (including the edge itself) that is closest to the query point
void Polyhedron::Edge::find_closest_feature(const Origin3d& p, std::map<shared_ptr<Polyhedron::Feature>, double>& distances, double& closest_dist, std::list<shared_ptr<Polyhedron::Feature> >& closest_features)
{
  // compute the distance from the edge itself (if it hasn't already been computed)
  if (distances.find(shared_from_this()) == distances.end())
  {
  }

  // check left face
  if (distances.find(face1) == distances.end())
    face1->find_closest_feature(p, distances, closest_dist, closest_features);  

  // check right face 
  if (distances.find(face2) == distances.end())
    face2->find_closest_feature(p, distances, closest_dist, closest_features);

  
}

/// Finds the closest feature to this face that is closer to this point
void Polyhedron::Face::find_closest_feature(const Origin3d& p, std::map<shared_ptr<Polyhedron::Feature>, double>& distances, double& closest_dist, std::list<shared_ptr<Polyhedron::Feature> >& closest_features)
{
  boost::shared_ptr<const Pose2d> GLOBAL2D;

  // get the plane containing this face
  Plane plane = get_plane();

  // get the distance of the point to the face
  Point3d pp(p, GLOBAL);
  double plane_dist = std::fabs(plane.calc_signed_distance(pp));

  // project the point onto the plane
  Ravelin::Vector2d p2d(plane.to_2D(pp), GLOBAL2D);

  // iterate over every edge in the ccw face-edge sequence; see whether 
  // point is always to the left side 
  bool strictly_inside = true;
  shared_ptr<Polyhedron::Edge> e(shared_from_this()->e.front()); 
  do 
  {
    // project the two vertices of the edge to 2D
    Ravelin::Vector2d ea2d(plane.to_2D(Point3d(e->v1->o, GLOBAL)), GLOBAL2D);
    Ravelin::Vector2d eb2d(plane.to_2D(Point3d(e->v2->o, GLOBAL)), GLOBAL2D);

    // check whether point is strictly inside 
    CompGeom::OrientationType otype = CompGeom::area_sign(p2d, ea2d, eb2d);
    if (otype != CompGeom::eLeft)
    {
      strictly_inside = false;
      break;
    }

    // get the next edge
    e = shared_ptr<Polyhedron::Edge>((e->face1 == shared_from_this()) ? e->prevR : e->prevL);
  }
  while (e != shared_ptr<Polyhedron::Edge>(shared_from_this()->e.front()));

  // case 1: point is strictly inside the convex polygon
  // this means we need to check all adjacent edges 
  if (strictly_inside)
  {
    // store the distance for this face 
    assert(distances.find(shared_from_this()) == distances.end());
    distances[shared_from_this()] = plane_dist;

    // see whether this distance is closer than the closest
    if (closest_dist > plane_dist)
    {
      closest_dist = plane_dist;
      closest_features.clear();
      closest_features.push_back(shared_from_this());
    }
    else if (closest_dist == plane_dist)
      closest_features.push_back(shared_from_this());

    // iterate over all edges
    e = shared_ptr<Polyhedron::Edge>(shared_from_this()->e.front());
    do
    {
      // don't check the edge twice
      if (distances.find(e) != distances.end())
      {
        e = shared_ptr<Polyhedron::Edge>((e->face1 == shared_from_this()) ? e->prevR : e->prevL);  
        continue;
      }

      // check closest features from this edge 
      e->find_closest_feature(p, distances, closest_dist, closest_features);   

      // get the next edge
      e = shared_ptr<Polyhedron::Edge>((e->face1 == shared_from_this()) ? e->prevR : e->prevL);  
    }
    while (e != shared_ptr<Polyhedron::Edge>(shared_from_this()->e.front()));

    return; 
  } 

  // case 2: point is on the border or outside the polygon; check the distance 
  // from each edge
  double min_dist = std::numeric_limits<double>::max();
  shared_ptr<Polyhedron::Edge> closest_edge;
  e = shared_ptr<Polyhedron::Edge>(shared_from_this()->e.front());
  do 
  {
    // don't process edge twice
    if (distances.find(e) != distances.end())
    {
      e = shared_ptr<Polyhedron::Edge>((e->face1 == shared_from_this()) ? e->prevR : e->prevL);  
      continue;
    }

    // compute the distance between the edge and the point 
    double t;
    LineSeg3 seg;
    Point3d closest;
    seg.first = Point3d(e->v1->o, GLOBAL);
    seg.first = Point3d(e->v2->o, GLOBAL);
    double dist = CompGeom::calc_dist(seg, pp, t, closest);

    if (dist < min_dist)
    {
      min_dist = dist;
      closest_edge = e;
    }

    // get the next edge
    e = shared_ptr<Polyhedron::Edge>((e->face1 == shared_from_this()) ? e->prevR : e->prevL);  
  }
  while (e != shared_ptr<Polyhedron::Edge>(shared_from_this()->e.front()));

  // record the distance for this face 
  assert(distances.find(shared_from_this()) == distances.end());
  distances[shared_from_this()] = min_dist;

  // see whether this distance is closer than the closest
  if (closest_dist > min_dist)
  {
    closest_dist = min_dist;
    closest_features.clear();
    closest_features.push_back(closest_edge);
  }
  else if (closest_dist == min_dist)
    closest_features.push_back(closest_edge);

  // check features adjacent to closest edge
  closest_edge->find_closest_feature(p, distances, closest_dist, closest_features);   
}
*/

/// Sends this polyhedron to the specified stream using VRML
void Polyhedron::to_vrml(std::ostream& out, const Polyhedron& p, Origin3d diffuse_color, bool wireframe)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup a mapping from vertices to indices
  map<shared_ptr<Polyhedron::Vertex>, unsigned> mapping;
  for (unsigned i=0; i< p._vertices.size(); i++)
    mapping[p._vertices[i]] = i;
  
  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material { diffuseColor " << diffuse_color[0] << " " << diffuse_color[1] << " " << diffuse_color[2] << " } }" << std::endl;
  out << "  geometry ";
  if (!wireframe)
    out << "IndexedFaceSet {" << std::endl;
  else
    out << "IndexedLineSet {" << std::endl;
  out << "    coord Coordinate { point [ ";
  for (unsigned i=0; i< p._vertices.size(); i++)
    out << p._vertices[i]->o[X] << " " << p._vertices[i]->o[Y] << " " << p._vertices[i]->o[Z] << ", ";
  out << " ] }" << std::endl;
  out << "    coordIndex [";
  for (unsigned i=0; i< p._faces.size(); i++)
  {
    // iterate over the vertices
    Polyhedron::VertexFaceIterator vfi(p._faces[i], true);

    // get the first vertex if it's a solid
    shared_ptr<Polyhedron::Vertex> v1 = *vfi;
    shared_ptr<Polyhedron::Vertex> v = v1;

    // loop through the vertices
    while (true)
    {
      out << " " << mapping[v];
      if (vfi.has_next())
      {
        vfi.advance();
        v = *vfi;
      }
      else
        break;
    }

    // if it's a solid, print the first vertex again
    if (!wireframe)
      out << " " << mapping[v1];

    // close it out
    out << " -1";
  }

  out << " ] } }" << std::endl;
}

// prints out the polyhedron
std::ostream& Moby::operator<<(std::ostream& out, const Polyhedron& m)
{
  const std::vector<shared_ptr<Polyhedron::Face> > faces = m.get_faces();
  const std::vector<shared_ptr<Polyhedron::Vertex> > vertices = m.get_vertices();
  const std::vector<shared_ptr<Polyhedron::Edge> > edges = m.get_edges();

  // setup mapping from vertices to indices
  std::map<shared_ptr<Polyhedron::Vertex>, unsigned> vmap;
  for (unsigned i=0; i< vertices.size(); i++)
  {
    out << "vertex " << i << ": " << vertices[i]->o << std::endl;
    vmap[vertices[i]] = i;
  }

  // setup mapping from edges to indices
  std::map<shared_ptr<Polyhedron::Edge>, unsigned> emap;
  for (unsigned i=0; i< edges.size(); i++)
    emap[edges[i]] = i;

  // process all faces
  for (unsigned i=0; i< faces.size(); i++)
  {
    shared_ptr<Polyhedron::Face> f = faces[i];

    // outputting face f
    out << "face " << i << " edges: " << f->e.size() << std::endl;
 
    // output edges from the face
    BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we, f->e)
    {
      // get the real edge
      shared_ptr<Polyhedron::Edge> e(we);
 
      // see whether the edge is to the left or right
      bool leftT = (e->face1 == f) ? false : true;

      if (leftT)
        out << "edge " << emap[e] << " (L): " << vmap[e->v1] << ", " << vmap[e->v2] << std::endl;
      else
        out << "edge  " << emap[e] << " (R): " << vmap[e->v1] << ", " << vmap[e->v2] << std::endl;
    }

    // output all vertices in counter clockwise order
    out << "vertices (ccw):";
    Polyhedron::VertexFaceIterator ccw_iter(f, true);
    out << " " << (*ccw_iter)->o << "; ";
    while (ccw_iter.has_next())
    {
      ccw_iter.advance();
      out << " " << (*ccw_iter)->o << "; ";
    }
    out << std::endl;

    // output all vertices in clockwise order
    out << "vertices (cw):";
    Polyhedron::VertexFaceIterator cw_iter(f, false);
    out << " " << (*cw_iter)->o << "; ";
    while (cw_iter.has_next())
    {
      cw_iter.advance();
      out << " " << (*cw_iter)->o << "; ";
    }
    out << std::endl;
  }

  return out;
}

// prints out a vertex 
std::ostream& Moby::operator<<(std::ostream& out, const Polyhedron::Vertex& m)
{
  out << m.o;
  return out;
}

// prints out an edge
std::ostream& Moby::operator<<(std::ostream& out, const Polyhedron::Edge& m)
{
  out << m.v1->o << " " << m.v2->o;
  return out;
}

// prints out a face 
std::ostream& Moby::operator<<(std::ostream& out, const Polyhedron::Face& m)
{
  // get an iterator to the edges in the face
  std::list<boost::weak_ptr<Polyhedron::Edge> >::const_iterator ccw_iter = m.e.begin();
  shared_ptr<Polyhedron::Edge> e1(*ccw_iter++);
  shared_ptr<Polyhedron::Edge> e2(*ccw_iter++);

  // get the common face
  shared_ptr<Polyhedron::Face> e1_faces[2];
  shared_ptr<Polyhedron::Face> e2_faces[2];
  shared_ptr<Polyhedron::Face> common_face;
  e1_faces[0] = e1->face1;  e1_faces[1] = e1->face2;
  e2_faces[0] = e2->face1;  e2_faces[1] = e2->face2;
  std::sort(e1_faces, e1_faces+2);
  std::sort(e2_faces, e2_faces+2);
  std::vector<shared_ptr<Polyhedron::Face> > result;
  std::set_intersection(e1_faces, e1_faces+2, e2_faces, e2_faces+2, std::back_inserter(result));
  assert(result.size() == 1);
  common_face = result.front();

  // setup the face iterator
  Polyhedron::VertexFaceIterator vfi(common_face, true); 
  out << (*vfi)->o;
  while (vfi.has_next())
  {
    vfi.advance();
    out << " " << (*vfi)->o;
  }

  return out;
}

/// Constructs a vertex-face iterator
Polyhedron::VertexFaceIterator::VertexFaceIterator(shared_ptr<Polyhedron::Face> f, bool ccw)
{
  // save this face
  this->f = f;

  // save whether it is counter-clockwise iteration
  this->ccw = ccw;

  // setup the iterator
  if (ccw)
    ccw_iter = f->e.begin();
  else
    cw_iter = f->e.rbegin(); 

  // pick a vertex from the edge
  if (ccw)
  {
    // get the next edge
    list<weak_ptr<Edge> >::const_iterator ccw_iter2 = ccw_iter;
    ccw_iter2++;

    // get the vertices
    shared_ptr<Edge> e(*ccw_iter);
    shared_ptr<Edge> e2(*ccw_iter2);

    // the starting vertex should be the vertex of e that is shared with e2 
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v1;
    else
      v = e->v2;
  }
  else
  {
    // get the next edge
    list<weak_ptr<Edge> >::const_reverse_iterator cw_iter2 = cw_iter;
    cw_iter2++;

    // get the vertices
    shared_ptr<Edge> e(*cw_iter);
    shared_ptr<Edge> e2(*cw_iter2);

    // the starting vertex should be the vertex of e that is shared with e2 
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v1;
    else
      v = e->v2;
  }
}

/// Dereferences the iterator
shared_ptr<Polyhedron::Vertex> Polyhedron::VertexFaceIterator::operator*()
{
  return v;
}

/// Advances the iterator clockwise
void Polyhedron::VertexFaceIterator::advance()
{
  // pick a vertex from the edge
  if (ccw)
  {
    // advance the edge iterator
    ccw_iter++;

    // get the next edge
    list<weak_ptr<Edge> >::const_iterator ccw_iter2 = ccw_iter;
    ccw_iter2++;
    if (ccw_iter2 == f->e.end())
    {
      // we've come to the end, return the last vertex
      shared_ptr<Edge> e(*ccw_iter);
      if (v == e->v1)
        v = e->v2;
      else
        v = e->v1;
      return;
    }

    // get the vertices
    shared_ptr<Edge> e(*ccw_iter);
    shared_ptr<Edge> e2(*ccw_iter2);

    // we want to set v equal to the vertex that is shared with e2 
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v1;
    else
      v = e->v2;
  }
  else
  {
    // advance the edge iterator
    cw_iter++;

    // get the next edge
    list<weak_ptr<Edge> >::const_reverse_iterator cw_iter2 = cw_iter;
    cw_iter2++;
    if (cw_iter2 == f->e.rend())
    {
      // we've come to the end, return the last vertex
      shared_ptr<Edge> e(*cw_iter);
      if (v == e->v1)
        v = e->v2;
      else
        v = e->v1;
      return;
    }

    // get the vertices
    shared_ptr<Edge> e(*cw_iter);
    shared_ptr<Edge> e2(*cw_iter2);

    // we want to set v equal to the vertex that is shared with e2 
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v1;
    else
      v = e->v2;
  }
}

/// Checks to see whether the iterator can be advanced clockwise
bool Polyhedron::VertexFaceIterator::has_next()
{
  // see whether we are at the end of the list
  if (ccw)
  {
    list<weak_ptr<Edge> >::const_iterator ccw_iter2 = ccw_iter;    
    ccw_iter2++;
    return (ccw_iter2 != f->e.end()); 
  }
  else
  {
    list<weak_ptr<Edge> >::const_reverse_iterator cw_iter2 = cw_iter;    
    cw_iter2++;
    return (cw_iter2 != f->e.rend());
  }
}

/// Finds closest points between two closest features on two polyhedra.
/// @param fA the closest feature on polyhedron A
/// @param fB the closest feature on polyhedron B
/// @param wTa the transformation between polyhedron A and the world
/// @param wTb the transformation between polyhedron B and the world
/// @param[out] pA a closest point on A
/// @param[out] pB a closest point on B
/// @param[out] normalsA the summed normal(s) at pA
/// @param[out] normalsB the summed normal(s) at pB
void Polyhedron::find_closest_points(boost::shared_ptr<const Polyhedron::Feature> fA, boost::shared_ptr<const Polyhedron::Feature> fB, const Ravelin::Transform3d& wTa, const Ravelin::Transform3d& wTb, Point3d& pA, Point3d& pB, Vector3d& normalsA, Vector3d& normalsB) {
  if (dynamic_pointer_cast<const Polyhedron::Vertex>(fA)) {

    // fA is an vertex
    shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(fA);

    // Set normalsA.
    normalsA = vA->sum_coincident_normals();

    // Setup the point from A.
    Ravelin::Vector3d vAa(vA->o, wTa.source);
    pA = wTa.transform_point(vAa);

    if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {

      // fB is a vertex
      shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(fB);

      // Set normalsB.
      normalsB = vB->sum_coincident_normals();

      // Setup the point from B
      Ravelin::Vector3d vBb(vB->o, wTb.source);

      //Transforming point into world frame
      pB = wTb.transform_point(vBb);
    }
    else if (dynamic_pointer_cast<const Polyhedron::Edge>(fB)) {

      // fB is an edge
      shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(fB);

      Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), wTb.source);
      Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), wTb.source);

      // transform the edge into the world frame
      Ravelin::Vector3d vB1w = wTb.transform_point(vB1b);
      Ravelin::Vector3d vB2w = wTb.transform_point(vB2b);

      // Create line segment
      LineSeg3 line(vB1w,vB2w);

      // Compute distance and closest point
      Point3d p;
      double t;
      double dist = CompGeom::calc_dist(line,pA,t,p);
      pB = p;

      // Set normal.
      normalsB = eB->face1->get_plane().get_normal() +
          eB->face2->get_plane().get_normal();
    }
    else {
      // fB is a face
      boost::shared_ptr<const Pose3d> GLOBAL3D;

      // cast features to non-constant
      boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(fA);
      boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(fB);
      boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);

      // Transform point from A into B's frame.
      Ravelin::Vector3d vAb = wTb.inverse_transform_point(wTa.transform_point(vAa));
      vAb.pose = GLOBAL3D;                 // hack around plane being in B's frame

      // Find the minimum
      Plane planeB = faceB->get_plane();   // plane will be in B's frame
      double dist = planeB.calc_signed_distance(vAb);

      // project the point onto the plane
      Ravelin::Vector3d vAb_on_planeB = vAb - planeB.get_normal()*dist;
      vAb_on_planeB.pose = wTb.source;
      // this is correct because when v-clip ends and a face vertex case
      // the vertex will always be in the voronoi region, and therefore,
      // the vertex projection is always on face B.
      vAb_on_planeB.pose = wTb.source;
      pB = wTb.transform_point(vAb_on_planeB);

      // Set normal.
      normalsB = faceB->get_plane().get_normal();
    }
  }
  else if (dynamic_pointer_cast<const Polyhedron::Edge>(fA)) {

    // fA is an edge
    if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {

      // already implemented, just need to flip it.
      find_closest_points(fB, fA, wTb, wTa, pB, pA, normalsB, normalsA);
    }
    else if (dynamic_pointer_cast<const Polyhedron::Edge>(fB)) {
      // Features are two edges.

      // Cast pointers
      boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(fA);
      boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(fB);

      // Set normals.
      normalsA = eA->face1->get_plane().get_normal() +
          eA->face2->get_plane().get_normal();
      normalsB = eB->face1->get_plane().get_normal() +
          eB->face2->get_plane().get_normal();

      //create vectors
      Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), wTa.source);
      Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), wTa.source);

      Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), wTb.source);
      Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), wTb.source);

      //transform features to global frame
      Ravelin::Vector3d vA1w = wTa.transform_point(vA1a);
      Ravelin::Vector3d vA2w = wTa.transform_point(vA2a);
      Ravelin::Vector3d vB1w = wTb.transform_point(vB1b);
      Ravelin::Vector3d vB2w = wTb.transform_point(vB2b);

      //create line segment
      LineSeg3 lineA(vA1w,vA2w);
      LineSeg3 lineB(vB1w,vB2w);

      // compute distance and closest point
      Point3d p1;
      Point3d p2;
      double dist = CompGeom::calc_closest_points(lineA,lineB,p1,p2);
      pA = p1;
      pB = p2;
    }
    else {

      // should not be reached b/c V-Clip does not return face/edge.
      assert(false);
    }
  }
  else {
    if (dynamic_pointer_cast<const Polyhedron::Vertex>(fB)) {

      // already implemented, just need to flip it.
      find_closest_points(fB, fA, wTb, wTa, pB, pA, normalsB, normalsA);
    }
    else{
      // should not be reached
      assert(false);
    }
  }

  // Set poses.
  normalsA.pose = wTa.source;
  normalsB.pose = wTb.source;
}

/// Summary: Executes the V-Clip algorithm on two polyhedra, determining closest features and signed distance

/* 
 Decription The function takes in two Polyhedra the Pose of the two Polyhedra and two initial feature.
 The function will then executes the V-Clip algorithm. If the two polyhedra are separated, the function will
 return the closest feature and returns the distance between the two object.
 If two objects are interpenetrating, the function will try to find the penetration distance using a heuristic
 that works only when one face is penetrated. If the heuristic is not available, the algorithm will return -1.0
 which will signal the simulator to perform another operation to find the sign distance.
*/
 double Polyhedron::vclip(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
 // return -1.0;
  FeatureType fA, fB;
  const Polyhedron& polyA = pA->get_polyhedron();
  const Polyhedron& polyB = pB->get_polyhedron();

  // defining the maximum iteration based on the number of total features 
  // in the two polyhedra
  int feature_sum = polyA.get_faces().size()+polyA.get_vertices().size()+polyA.get_edges().size()
                  + polyB.get_faces().size()+polyB.get_vertices().size()+polyB.get_edges().size();
  int MAX_ITERATION = feature_sum;

  // get the transformation between A and B, and vice versa
  Transform3d aTb = Pose3d::calc_relative_pose(poseB, poseA);
  Transform3d bTa = aTb.inverse(); 

  FILE_LOG(LOG_COLDET) << "===========V-clip started========"<< std::endl;
  FILE_LOG(LOG_COLDET) << "poseA: "<< *poseA << std::endl;
  FILE_LOG(LOG_COLDET) << "poseB: "<< *poseB << std::endl;
  FILE_LOG(LOG_COLDET) << "aTb: " << aTb << std::endl;

  // if closest feature of A is null, pick features for A and B arbitrarily
  if(!closestA){
    std::vector<boost::shared_ptr<Face> >::const_iterator vi = pA->get_polyhedron().get_faces().begin();
    closestA = boost::shared_ptr<Feature>(*vi);
    vi = pB->get_polyhedron().get_faces().begin();
    closestB = boost::shared_ptr<Feature>(*vi);
  }

  // determine feature type for A
  if (dynamic_pointer_cast<const Polyhedron::Vertex>(closestA))
    fA = eVertex;
  else if (dynamic_pointer_cast<const Polyhedron::Edge>(closestA))
    fA = eEdge;
  else
    fA = eFace;

  // determine feature type for B 
  if (dynamic_pointer_cast<const Polyhedron::Vertex>(closestB))
    fB = eVertex;
  else if (dynamic_pointer_cast<const Polyhedron::Edge>(closestB))
    fB = eEdge;
  else
    fB = eFace;

  int iteration_count = 0;
  // iterate through the algorithm
  while (true)
  {
    iteration_count++;
    // handle vertex/vertex case
    if (fA == eVertex && fB == eVertex)
    {
      boost::shared_ptr<const Polyhedron::Vertex> a = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
      boost::shared_ptr<const Polyhedron::Vertex> b = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Vertex Vertex Case=====" << std::endl <<*a << std::endl << *b <<std::endl;

      Polyhedron::UpdateRule r = update_vertex_vertex(fA, fB, aTb, closestA, closestB);

      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        return dist; 
      }
    }
    // handle vertex/edge cases
    else if (fA == eVertex && fB == eEdge)
    {
      boost::shared_ptr<const Polyhedron::Vertex> a = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
      boost::shared_ptr<const Polyhedron::Edge> b = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Vertex Edge Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
      Polyhedron::UpdateRule r = update_vertex_edge(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        return dist; 
      }
    }
    else if (fB == eVertex && fA == eEdge)
    {
      boost::shared_ptr<const Polyhedron::Edge> a = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
      boost::shared_ptr<const Polyhedron::Vertex> b = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Edge Vertex Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
      Polyhedron::UpdateRule r = update_vertex_edge(fB, fA, bTa, closestB, closestA);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fB, fA, closestB, closestA, bTa);

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        return dist; 
      }
    }
    // handle edge/edge case
    else if (fA == eEdge && fB == eEdge)
    {
      boost::shared_ptr<const Polyhedron::Edge> a = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
      boost::shared_ptr<const Polyhedron::Edge> b = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
      FILE_LOG(LOG_COLDET)<< "=====Entering Edge Edge Case=====" << std::endl <<*a << std::endl << *b <<std::endl;

      Polyhedron::UpdateRule r = update_edge_edge(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        return dist; 
      }
    }
    // handle edge/face cases
    else if (fA == eEdge && fB == eFace)
    {
      boost::shared_ptr<const Polyhedron::Edge> a = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
      boost::shared_ptr<const Polyhedron::Face> b = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Edge Face Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
      Polyhedron::UpdateRule r = update_edge_face(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm

      if(r == eContinue && iteration_count > MAX_ITERATION)
      {
        FILE_LOG(LOG_COLDET) << "v-clip maximum iteration number has been reached; calculating distance" << std::endl;
        FILE_LOG(LOG_COLDET) << " using the current feature" << std::endl;
        r=eInterpenetrating;
      }

      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist;

      //start searching for penetration 

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        if(is_one_face_penetration(closestB, pA, closestA, bTa))
        { 
          FILE_LOG(LOG_COLDET) << "Quick Case" << std::endl;
          find_deepest_feature(closestB, closestA, fA, bTa);
          dist = calc_dist(fA, fB, closestA, closestB, aTb);
        }
        else
        {
          FILE_LOG(LOG_COLDET) << "Slow Case" << std::endl;          
          return -1.0;//minkowski_optimum_distance(pA,pB,aTb);
        }

        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        dist = calc_dist(fA, fB, closestA, closestB, aTb);
        return dist; 
      }
    }
    else if (fB == eEdge && fA == eFace)
    {
      boost::shared_ptr<const Polyhedron::Face> a = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
      boost::shared_ptr<const Polyhedron::Edge> b = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Face Edge Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
     
      Polyhedron::UpdateRule r = update_edge_face(fB, fA, bTa, closestB, closestA);

      if(r == eContinue && iteration_count > MAX_ITERATION)
      {
        FILE_LOG(LOG_COLDET) << "v-clip maximum iteration number has been reached; calculating distance" << std::endl;
        FILE_LOG(LOG_COLDET) << " using the current feature" << std::endl;
        r=eInterpenetrating;
      }

      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist;

      //start searching for penetration 

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        if(is_one_face_penetration(closestA, pB, closestB, aTb))
        { 
          FILE_LOG(LOG_COLDET) << "Quick Case" << std::endl;
          return -1.0;
          find_deepest_feature(closestA, closestB, fB, aTb);
          dist = calc_dist(fA, fB, closestA, closestB, aTb);
        }
        else
        {
          FILE_LOG(LOG_COLDET) << "Slow Case" << std::endl;
          return -1.0;//minkowski_optimum_distance(pA,pB,aTb);
        }

        FILE_LOG(LOG_COLDET)<< "penetrating" << std::endl;
        return -dist;
      }
      else
      {
        dist = calc_dist(fA, fB, closestA, closestB, aTb);
        return dist; 
      }
    }
    // handle face/face case
    else if (fA == eFace && fB == eFace)
    {
      // demote face b to an edge
      // cast pointer
      boost::shared_ptr<const Polyhedron::Face> faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
      std::list<boost::weak_ptr<Edge> >::const_iterator ei = (faceB->e).begin();
      closestB = boost::shared_ptr<const Polyhedron::Feature>(*ei);
      fB = eEdge;
      continue;
    }
    // handle vertex/face case
    else if (fA == eVertex && fB == eFace)
    {
      boost::shared_ptr<const Polyhedron::Vertex> a = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
      boost::shared_ptr<const Polyhedron::Face> b = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Vertex Face Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
     
      Polyhedron::UpdateRule r = update_vertex_face(fA, fB, aTb, closestA, closestB, pB->get_polyhedron());

      if(r == eContinue && iteration_count > MAX_ITERATION)
      {
        FILE_LOG(LOG_COLDET) << "v-clip maximum iteration number has been reached; calculating distance" << std::endl;
        FILE_LOG(LOG_COLDET) << " using the current feature" << std::endl;
        r=eInterpenetrating;
      }
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);
      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET) << "Slow Case" << std::endl;
        return -1.0;
      }
      else
      {
        return dist; 
      }

    }    //handle face/vertex case
    else if (fB == eVertex && fA == eFace)
    {
      boost::shared_ptr<const Polyhedron::Face> a = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
      boost::shared_ptr<const Polyhedron::Vertex> b = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
      FILE_LOG(LOG_COLDET) << "=====Entering Face Vertex Case=====" << std::endl <<*a << std::endl << *b <<std::endl;
     
      Polyhedron::UpdateRule r = update_vertex_face(fB, fA, bTa, closestB, closestA, pA->get_polyhedron());

      if(r == eContinue && iteration_count > MAX_ITERATION)
      {
        FILE_LOG(LOG_COLDET) << "v-clip maximum iteration number has been reached; calculating distance" << std::endl;
        FILE_LOG(LOG_COLDET) << " using the current feature" << std::endl;
        r=eInterpenetrating;
      }
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fB, fA, closestB, closestA, bTa);

      FILE_LOG(LOG_COLDET)<< "converged" << std::endl;
      if (r == eInterpenetrating)
      {
        FILE_LOG(LOG_COLDET) << "Slow Case" << std::endl;
        return -1.0;
      }
      else
      {
        return dist; 
      }
    }
  }
  }


/// Check if only one face of the two polyhedron is penetrated by checking if the edge of the penetrated face is penetrating
/// assuming no fully interpenetration
/// TODO: add fully interpenetration detection
bool Polyhedron::is_one_face_penetration(boost::shared_ptr<const Polyhedron::Feature>& fFace, boost::shared_ptr<const PolyhedralPrimitive> pE, boost::shared_ptr<const Polyhedron::Feature>& fEdge, Ravelin::Transform3d& fTe)
{
  // 1.Cast pointers 
  boost::shared_ptr<const Polyhedron::Face> f = boost::static_pointer_cast<const Polyhedron::Face>(fFace);
  boost::shared_ptr<const Polyhedron::Edge> e = boost::static_pointer_cast<const Polyhedron::Edge>(fEdge);

  // Creating triangles base on the vertices of the faces in pE
  // Assumption : All the faces are triangular
  // The code should have a better performance if the faces are sorted based on how close they are to the edge

  std::vector<boost::shared_ptr<Polyhedron::Face> > facesE = pE->get_polyhedron().get_faces();
  std::vector<Triangle> triangles;
  for(std::vector<boost::shared_ptr<Polyhedron::Face> >::const_iterator fEi = facesE.begin(); fEi != facesE.end(); ++fEi)
  {
    boost::shared_ptr<Polyhedron::Face> curFace(*fEi);
    Polyhedron::VertexFaceIterator vfi(curFace, true);
    std::vector<Ravelin::Vector3d> vface;
    FILE_LOG(LOG_COLDET)<<"Face: " << *curFace << std::endl;
    while (vfi.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> temp_vertex = *vfi;
      Ravelin::Vector3d vertex_vector_e(temp_vertex->o, fTe.source);
      Ravelin::Vector3d vertex_vector_f = fTe.transform_point(vertex_vector_e);
      vface.push_back(vertex_vector_f);
      vfi.advance();
    }

    //Adding the last vertex into the vector since the last vertex is not added
    boost::shared_ptr<Polyhedron::Vertex> temp_vertex = *vfi;
    Ravelin::Vector3d vertex_vector_e(temp_vertex->o, fTe.source);
    Ravelin::Vector3d vertex_vector_f = fTe.transform_point(vertex_vector_e);
    vface.push_back(vertex_vector_f);

    if(vface.size() == 3)
    {
      Triangle t(vface[0], vface[1], vface[2]);
      triangles.push_back(t);
    }
    else if(vface.size() > 3)
    {
      FILE_LOG(LOG_COLDET) << "faces are not triangles"<<std::endl;// TODO: Add triangulating so that the algorithm works for non-triangular faces
      CompGeom::triangulate_convex_polygon(vface.begin(),vface.end(),std::back_inserter(triangles));
    }
    
  }

  FILE_LOG(LOG_COLDET) << "finish adding triangles"<<std::endl;

  // 2.For all edge in fFace
  std::list<boost::weak_ptr<Polyhedron::Edge> > edgesF = f->e;

  for(std::list<boost::weak_ptr<Polyhedron::Edge> >::const_iterator eiF = edgesF.begin(); eiF != edgesF.end(); ++eiF)
  {
    // create segment
    boost::shared_ptr<Polyhedron::Edge> cur_edge(*eiF);
    Ravelin::Vector3d v1 (cur_edge->v1->o, fTe.target);
    Ravelin::Vector3d v2 (cur_edge->v2->o, fTe.target);
    LineSeg3 edge_seg(v1,v2);

    
    Point3d p1,p2;

    // check if the edge intersect with any triangles
    for(std::vector<Triangle>::const_iterator ti = triangles.begin(); ti != triangles.end(); ++ti)
    {
      if(CompGeom::intersect_seg_tri(edge_seg,*ti,p1,p2) != CompGeom::eSegTriNoIntersect)
        return false;
    }
  }

  //if it does not return true.
  return true;

}


// This function detects the feature on a polyhedron that has the most negative distance to a face from another polyhedron
// The function is used as a fast heuristics to find the penetration distance when only one face is penetrated.
void Polyhedron::find_deepest_feature(boost::shared_ptr<const Polyhedron::Feature>& face, boost::shared_ptr<const Polyhedron::Feature>& edge, Moby::Polyhedron::FeatureType& fE, Ravelin::Transform3d& fTe)
{
  // 1.cur_vertex = the vertices of the edge that is penetrating 

  boost::shared_ptr<Polyhedron::Vertex> cur_vertex; 
  boost::shared_ptr<const Polyhedron::Edge> e = boost::static_pointer_cast<const Polyhedron::Edge>(edge);
  boost::shared_ptr<const Polyhedron::Face> f = boost::static_pointer_cast<const Polyhedron::Face>(face);

  Plane p = f->get_plane();
  Ravelin::Vector3d norm = p.get_normal();
  norm.pose = fTe.target;
  p.set_normal(norm);
  boost::shared_ptr<Polyhedron::Vertex> v1 = e->v1;

  Ravelin::Vector3d v1_e(v1->o,fTe.source);
  Ravelin::Vector3d v1_f = fTe.transform_point(v1_e);

  // because it is already known that the edge is intersecting with the face, at
  // least one vertex should have a negative distance tward the face.
  if(p.calc_signed_distance(v1_f)<NEAR_ZERO)
  {
    cur_vertex = v1;
  }
  else
  {
    cur_vertex = e->v2;
  }


  while(true)
  {
    boost::shared_ptr<Polyhedron::Vertex> next_vertex;
    Ravelin::Vector3d cur_vector_e (cur_vertex->o, fTe.source);
    Ravelin::Vector3d cur_vector_f = fTe.transform_point(cur_vector_e);
    double cur_dist = p.calc_signed_distance(cur_vector_f);
    std::list<boost::weak_ptr<Polyhedron::Edge> > es = cur_vertex->e;

    for (std::list<boost::weak_ptr<Polyhedron::Edge> >::const_iterator ei = es.begin(); ei != es.end(); ++ei)
    {

      //find the vertex that is not cur_vertex
      boost::shared_ptr<const Polyhedron::Edge> cur_e(*ei);
      boost::shared_ptr<Polyhedron::Vertex> v = cur_e->v1;

      if(v == cur_vertex)
      {
        v = cur_e->v2;
      }

      FILE_LOG(LOG_COLDET) << "v: " << *v << std::endl;

      Ravelin::Vector3d v_e(v->o,fTe.source);
      Ravelin::Vector3d v_f = fTe.transform_point(v_e);
      double v_dist = p.calc_signed_distance(v_f);

      FILE_LOG(LOG_COLDET) << "dist: " << cur_dist - v_dist<<std::endl;

      if(cur_dist - v_dist > NEAR_ZERO)
      {
        next_vertex = v;
        cur_dist = v_dist;
      }
        
    }

    //nothing to update to
    if(!next_vertex)
    {
      break;
    }
    else
    {
      cur_vertex = next_vertex;
    }


  }

  edge = cur_vertex;
  fE=eVertex;
  return;
}


// this function is used to calculate the minimum translational distance to separate two polyhedron
// using the minkowski difference approach
double Polyhedron::minkowski_optimum_distance(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, Ravelin::Transform3d& aTb)
{
  shared_ptr<const Pose3d> GLOBAL3D;
  Ravelin::Origin3d o(0,0,0);
  Ravelin::Vector3d origin_vector(o,GLOBAL3D);
  Polyhedron mdiff = Polyhedron::calc_minkowski_diff(pA, pB, aTb.target, aTb.source);
  std::vector<boost::shared_ptr<Polyhedron::Face> > faces = mdiff.get_faces();
  double min_dist = std::numeric_limits<double>::max();
  
  for(std::vector<boost::shared_ptr<Polyhedron::Face> >::const_iterator fi = faces.begin(); fi != faces.end(); ++fi)
  {

    double dist = fabs((*fi)->get_plane().calc_signed_distance(origin_vector));
    FILE_LOG(LOG_COLDET) << dist << std::endl;
    if(min_dist - dist > NEAR_ZERO)
    {
      min_dist = dist;
    }
  }

  return min_dist;

}
/*// promotes the result features to the highest possible dimension
void promote_featrues(FeatureType& fA, FeatureType& fB, boost::shared_ptr<const Polyhedron::Feature>& closestA, boost::shared_ptr<const Polyhedron::Feature>& closestB, Ravelin::Transform3d& aTb)
{
  dist = calc_dist(fA, fB, closestA, closestB, aTb);
  if(fA == eVertex && fB == eVertex)
  {
    // Cast pointers
    boost::shared_ptr<const Polyhedron::Vertex> vertA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
    boost::shared_ptr<const Polyhedron::Vertex> vertB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);

    std::list<boost::weak_ptr<Edge> > edgesA = vertA->e;
    std::list<boost::weak_ptr<Edge> > edgesB = vertB->e;

    // look at the edges coincide to B, check if any of them are as much distant as B from A 
    for(std::list<boost::weak_ptr<Edge> >::const_iterator eiA = edgesA.begin(); eiA != edgesA.end(); ++eiA)
    {
      boost::shared_ptr<const Polyhedron::Edge> edgeA(*eiA);
      for(std::list<boost::weak_ptr<Edge> >::const_iterator eiB = edgesB.begin(); eiB != edgesB.end(); ++eiB)
      {
        boost::shared_ptr<const Polyhedron::Edge> edgeB(*eiB);
        boost::shared_ptr<const Polyhedron::Feature> vertB_2 = edgeB->v1;
        if(vertB_2 == vertB)
        {
          vertB_2 == edgeB->v2;
        }
        dist_edge = calc_dist(eEdge, eVertex, edgeA, vertB_2, aTb);

        if(dist_edge == dist){
          fA = eEdge;
          fB = eEdge;
          closestB = edgeB;
          closestA = edgeA;
          return promote_featrues(fA, fB, closestA, closestB, aTb);
        }
      }
    }

    // swap the role of A and B and repeat
    for(std::list<boost::weak_ptr<Edge> >::const_iterator eiB = edgesB.begin(); eiB != edgesB.end(); ++eiB)
    {
      boost::shared_ptr<const Polyhedron::Edge> edgeB(*eiB);
      for(std::list<boost::weak_ptr<Edge> >::const_iterator eiA = edgesA.begin(); eiA != edgesA.end(); ++eiA)
      {
        boost::shared_ptr<const Polyhedron::Edge> edgeA(*eiA);
        boost::shared_ptr<const Polyhedron::Feature> vertA_2 = edgeA->v1;
        if(vertA_2 == vertA)
        {
          vertA_2 == edgeA->v2;
        }
        dist_edge = calc_dist(eEdge, eVertex, edgeA, vertA_2, aTb);

        if(dist_edge == dist){
          fA = eEdge;
          fB = eEdge;
          closestB = edgeB;
          closestA = edgeA;
          return promote_featrues(fA, fB, closestA, closestB, aTb);
        }
      }
    }

    
  }
  else if(fA == eVertex && fB == eEdge)
  {

  }
  else if (fA == eVertex && fB == eFace)
  {

  }
  else if (fA == eEdge && fB == eEdge)
  {

  }
  else if (fA == eEdge && fB == eFace)
  {

  }
  // if A is in a higher dimension than B, swap A and B and call function again
  else
  {
   Ravelin::Transform3d bTa = aTb.inverse();
    promote_featrues(fB, fA, closestB, closestA, bTa);
  }

}
*/

/// Calculates the signed distance of the polyhedron from the point
/**
 * \param closest_facet the closest facet to the point on return
 */
double Polyhedron::calc_signed_distance(const Origin3d& p, unsigned& closest_facet) const
{
  throw std::runtime_error("Implement me!");
}

/// Finds the closest feature of the polyhedron to the point, given the closest facet
/**
 * \param closest_facet the closest facet to the point on return
 */
boost::shared_ptr<Polyhedron::Feature> Polyhedron::find_closest_feature(const Origin3d& p, unsigned closest_facet) const
{
  throw std::runtime_error("Implement me!");
}

/// Computes the distance between two features
double Polyhedron::calc_dist(FeatureType fA, FeatureType fB, boost::shared_ptr<const Polyhedron::Feature> closestA, boost::shared_ptr<const Polyhedron::Feature> closestB, Ravelin::Transform3d& aTb)
{

//  case1: vertex vs. vertex
  if (fA == eVertex && fB == eVertex)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
    boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
    
    // create vectors
    Ravelin::Vector3d vAa(vA->o, aTb.target);
    Ravelin::Vector3d vBb(vB->o, aTb.source);
    
    //Transforming B into frame A
    Ravelin::Vector3d vBa = aTb.transform_point(vBb);
    Ravelin::Vector3d diff = vAa-vBa;
    return diff.norm();

  } 
  else if (fA == eVertex && fB == eEdge)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
    boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);

    // create vectors for the vertex
    Ravelin::Vector3d vAa(vA->o, aTb.target);
    Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
    Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
    
    // transform fB into frame A
    Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
    Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

    //create line segment
    LineSeg3 line(vB1a,vB2a);
     
    //compute distance
    Point3d p;
    double t;
    double dist = CompGeom::calc_dist(line,vAa,t,p);
    return dist;
  } 
  else if (fA == eEdge && fB == eVertex)
  {
    // casting pinters
    boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
    boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);

    // create vectors for the vertex
    Ravelin::Vector3d vBb(vB->o, aTb.source);
    Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
    Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);
    
    //transform fB into frame A
    Ravelin::Vector3d vBa = aTb.transform_point(vBb);

    //create line segment
    LineSeg3 line(vA1a,vA2a);
     
    //computing distance
    Point3d p;
    double t;
    double dist = CompGeom::calc_dist(line,vBa,t,p);
    return dist;
  }
  else if(fA == eEdge && fB == eEdge)
  {
    //casting pointers
    boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
    boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
    
    //create vectors
    Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
    Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);

    Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
    Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
    
    //transform fB into frame A
    Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
    Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

    //create line segment
    LineSeg3 lineA(vA1a,vA2a);
    LineSeg3 lineB(vB1a,vB2a);

    // compute distance
    Point3d p1;
    Point3d p2;
    double dist = CompGeom::calc_closest_points(lineA,lineB,p1,p2);
    return dist;
  }
  else if (fA == eVertex && fB == eFace)
  {
    //creating null pointer for later use as place holder
    boost::shared_ptr<const Pose2d> GLOBAL2D;
    boost::shared_ptr<const Pose3d> GLOBAL3D;
    //cast features to non-constant
    boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
    boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);     

    FILE_LOG(LOG_COLDET) << "Closest A" << *vA <<std::endl;
    FILE_LOG(LOG_COLDET) << "Closest B" << *faceB <<std::endl;

    // create vector for A
    Ravelin::Vector3d vAa(vA->o, aTb.target);
  
    // transform A using the inverse of aTb
    Ravelin::Transform3d bTa = aTb.inverse();
    Ravelin::Vector3d vAb = bTa.transform_point(vAa);
    vAb.pose = GLOBAL3D;
    // find the minimum
    Plane planeB = faceB->get_plane();
    double dist = planeB.calc_signed_distance(vAb);

    //FILE_LOG(LOG_COLDET) << "Distance between plane and " << vAb<< " is "<< dist << std::endl;

    // project the point onto the plane
    Ravelin::Vector3d vAb_on_planeB = vAb - planeB.get_normal()*dist;

    //done using vAb, resetting the frame
    vAb.pose = aTb.source;

    // create the projection matrix to project to 2D
    Ravelin::Matrix3d R2D = CompGeom::calc_3D_to_2D_matrix(planeB.get_normal());
    Ravelin::Vector2d vAb_on_planeB_2d(CompGeom::to_2D(vAb_on_planeB,R2D),GLOBAL2D);

    Polyhedron::VertexFaceIterator vfiB(faceB,true);    
    std::vector<Ravelin::Vector2d> vB2d;
    
    //create the vertex vector
    while (vfiB.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v = *vfiB;
      vfiB.advance();
      Ravelin::Vector3d p(v->o, aTb.source);
      Ravelin::Origin2d o=CompGeom::to_2D(p,R2D);
      Ravelin::Vector2d v2d(o,GLOBAL2D);
      vB2d.push_back(v2d);
    }
    boost::shared_ptr<Polyhedron::Vertex> v=*vfiB;
    Ravelin::Vector3d p(v->o, aTb.source);
    Ravelin::Origin2d o=CompGeom::to_2D(p,R2D);
    Ravelin::Vector2d v2d(o,GLOBAL2D);
    vB2d.push_back(v2d);

    // find the relationship between the point and the face
    CompGeom::PolygonLocationType relation = CompGeom::polygon_location(vB2d.begin(), vB2d.end(), vAb_on_planeB_2d);

    // if the point is outside we need to find distance between edges
    if (relation==CompGeom::ePolygonOutside)
    {
      //create the vector for B
      std::list<boost::weak_ptr<Edge> > eBs = faceB->e;
  
      // initialize the minimum
      double min_dist = std::numeric_limits<double>::max();

      // finding the distance between edges
      for(std::list<boost::weak_ptr<Edge> >::iterator eBsi = eBs.begin(); eBsi!=eBs.end();++eBsi)
      {
        boost::shared_ptr<Edge> eB(*eBsi);
        Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
        Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
    
        //transform fB into frame A
        Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
        Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

        //create a line segment
        LineSeg3 line(vB1a,vB2a);

        //compute the distance and compare to minimum
        Point3d p;
        double t;
        double dist = CompGeom::calc_dist(line,vAa,t,p);
        //FILE_LOG(LOG_COLDET) << "Distance between edges and " << vAa<< " is "<< dist << std::endl;
        if(min_dist > dist)
          min_dist = dist;
      }

      return min_dist;
    }
    else
      return std::fabs(dist);
  }
  else if (fA == eFace && fB == eVertex)
  {
    // creating null pointer for later use as place holder
    boost::shared_ptr<const Pose2d> GLOBAL2D;
    boost::shared_ptr<const Pose3d> GLOBAL3D;
    // cast pointers
    boost::shared_ptr<const Polyhedron::Face> _faceA = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<Polyhedron::Face> faceA = boost::const_pointer_cast<Polyhedron::Face>(_faceA); 
    boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
    

    FILE_LOG(LOG_COLDET) << "Closest A" << *faceA <<std::endl;
    FILE_LOG(LOG_COLDET) << "Closest B" << *vB <<std::endl;

    //transform B
    Ravelin::Vector3d vBb(vB->o, aTb.source);
    Ravelin::Vector3d vBa = aTb.transform_point(vBb);

    //find distance between point and the plane the face is on
    Plane planeA = faceA->get_plane();
    vBa.pose = GLOBAL3D;
    double dist = planeA.calc_signed_distance(vBa);
    
    // project the point to the plane
    Ravelin::Vector3d vBa_on_planeA = vBa - planeA.get_normal()*dist;
    vBa.pose = aTb.target;


    // create the projection matrix to project shapes on to 2D
    Ravelin::Matrix3d R2D = CompGeom::calc_3D_to_2D_matrix(planeA.get_normal());
    Ravelin::Vector2d vBa_on_planeA_2d(CompGeom::to_2D(vBa_on_planeA,R2D),GLOBAL2D);

    Polyhedron::VertexFaceIterator vfiA(faceA,true);    
    std::vector<Ravelin::Vector2d> vA2d;
    
    // create the vertex vector
    while (vfiA.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v = *vfiA;
      vfiA.advance();
      Ravelin::Vector3d p(v->o, aTb.target);
      Ravelin::Origin2d o=CompGeom::to_2D(p,R2D);
      Ravelin::Vector2d v2d(o,GLOBAL2D);
      vA2d.push_back(v2d);
    }
    boost::shared_ptr<Polyhedron::Vertex> v = *vfiA;
    Ravelin::Vector3d p(v->o, aTb.target);
    Ravelin::Origin2d o=CompGeom::to_2D(p,R2D);
    Ravelin::Vector2d v2d(o,GLOBAL2D);
    vA2d.push_back(v2d);

    //find location of point w.r.t. polygon
    CompGeom::PolygonLocationType relation = CompGeom::polygon_location(vA2d.begin(), vA2d.end(), vBa_on_planeA_2d);

    //if the point is outside, vertex / edge distance needs to be calculated
    if (relation == CompGeom::ePolygonOutside)
    {
      // if outside plane, minimum distance is computed by comparing edges
      // set the minimum
      double min_dist=std::numeric_limits<double>::max();

      //creating vertexes for edges of a
      std::list<boost::weak_ptr<Edge> > eAs=faceA->e;
      for (std::list<boost::weak_ptr<Edge> >::iterator eAsi = eAs.begin(); eAsi!=eAs.end();++eAsi)
      {
        boost::shared_ptr<Edge> eA(*eAsi);
        Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
        Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);

        // create line segment
        LineSeg3 line(vA1a,vA2a);

        // compute distance and comparing to minimum
        Point3d p;
        double t;
        double dist = CompGeom::calc_dist(line,vBa,t,p);
        //FILE_LOG(LOG_COLDET) << "Distance between edges and" << vBa<< " is "<< dist << std::endl;
        if(min_dist>dist)
          min_dist=dist;
      }

      return min_dist;
    }
    else
      return std::fabs(dist);
  }
  else if (fA == eEdge && fB == eFace)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
    boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);
  
    // create lineseg for A
    Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
    Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);
    LineSeg3 lineA(vA1a,vA2a);
    FILE_LOG(LOG_COLDET) << *faceB <<std::endl;

    // transform B
    Polyhedron::VertexFaceIterator vfiBb(faceB,true);
    std::vector<Ravelin::Vector3d> vBa;
    while(vfiBb.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v=*vfiBb;
      vfiBb.advance();
      Ravelin::Vector3d p(v->o, aTb.source);
      Ravelin::Vector3d pa = aTb.transform_point(p);
      vBa.push_back(pa);
    }

    // add the last point in the VertexFaceIterator into the vector
    boost::shared_ptr<Polyhedron::Vertex> v = *vfiBb;
    Ravelin::Vector3d p(v->o, aTb.source);
    Ravelin::Vector3d pa = aTb.transform_point(p);
    vBa.push_back(pa);

    // triangulate the polygon
    std::vector<Triangle> triB;
    CompGeom::triangulate_convex_polygon(vBa.begin(),vBa.end(),std::back_inserter(triB));
      
    // find the minimum distance between line and all triangles
    Point3d p1,p2;
    double min_dist = std::numeric_limits<double>::max();
    for(std::vector<Triangle>::iterator t = triB.begin(); t!=triB.end();++t)
    {
      //FILE_LOG(LOG_COLDET)<< *t <<std::endl;
      double tmp=Triangle::calc_sq_dist(*t,lineA,p1,p2);
      //  double tmp=0;
      if (tmp < min_dist)
        min_dist=tmp;
    }

    return std::sqrt(min_dist);
  }
  else if (fA == eFace && fB == eEdge)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Face> faceA_const = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<Polyhedron::Face> faceA = boost::const_pointer_cast<Polyhedron::Face>(faceA_const);
    boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);

    // creating iterator for A
    Polyhedron::VertexFaceIterator vfiAa(faceA,true);
    std::vector<Ravelin::Vector3d> vAa;
    while(vfiAa.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v=*vfiAa;
      vfiAa.advance();
      Ravelin::Vector3d pa(v->o, aTb.target);
      vAa.push_back(pa);
    }
    boost::shared_ptr<Polyhedron::Vertex> v=*vfiAa;
    Ravelin::Vector3d pa(v->o, aTb.target);
    vAa.push_back(pa);    

    // triangulate the polygon
    std::vector<Triangle> triA;
    CompGeom::triangulate_convex_polygon(vAa.begin(),vAa.end(),std::back_inserter(triA));
      
    // transform and create lineSegB
    Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
    Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
    
    // transform fB into frame A
    Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
    Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

    // create segment
    LineSeg3 lineB(vB1a,vB2a);

    // find minimum distance between line and all triangles
    Point3d p1,p2;
    double min_dist = std::numeric_limits<double>::max();
    for (std::vector<Triangle>::iterator t = triA.begin(); t!=triA.end();++t)
    {
      double tmp=Triangle::calc_sq_dist(*t,lineB,p1,p2);
      if(tmp < min_dist)
        min_dist=tmp;
    }

    return std::sqrt(min_dist);
  }
  else 
  {
    assert(fA == eFace && fB == eFace);

    // cast pointers
    boost::shared_ptr<const Polyhedron::Face> faceA_const = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<Polyhedron::Face> faceA = boost::const_pointer_cast<Polyhedron::Face>(faceA_const);
    boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);

    // create iterator for A
    Polyhedron::VertexFaceIterator vfiAa(faceA,true);
    std::vector<Ravelin::Vector3d> vAa;
    while (vfiAa.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v = *vfiAa;
      vfiAa.advance();
      Ravelin::Vector3d pa(v->o, aTb.target);
      vAa.push_back(pa);
    }
    boost::shared_ptr<Polyhedron::Vertex> v = *vfiAa;
    Ravelin::Vector3d pa(v->o, aTb.target);
    vAa.push_back(pa);

    // triangulate
    std::vector<Triangle> triA;
    CompGeom::triangulate_convex_polygon(vAa.begin(),vAa.end(),std::back_inserter(triA));
  
    // create iterator for B
    Polyhedron::VertexFaceIterator vfiBb(faceB,true);
    std::vector<Ravelin::Vector3d> vBa;
    while(vfiBb.has_next())
    {
      boost::shared_ptr<Polyhedron::Vertex> v = *vfiBb;
      vfiBb.advance();
      Ravelin::Vector3d p(v->o, aTb.source);
      Ravelin::Vector3d pa = aTb.transform_point(p);
      vBa.push_back(pa);
    }
    v = *vfiBb;
    Ravelin::Vector3d p(v->o, aTb.source);
    pa = aTb.transform_point(p);
    vBa.push_back(pa);    
    
    // triangulate
    std::vector<Triangle> triB;
    CompGeom::triangulate_convex_polygon(vBa.begin(),vBa.end(),std::back_inserter(triB));
      
    // find minimum distance between all pairs of triangles
    Point3d p1,p2;
    double min_dist = std::numeric_limits<double>::max();
    for(std::vector<Triangle>::iterator tA = triA.begin(); tA!=triA.end();++tA)
      for(std::vector<Triangle>::iterator tB = triB.begin(); tB!=triB.end();++tB)
      {
        double tmp = Triangle::calc_sq_dist(*tA,*tB,p1,p2);
        if (tmp < min_dist)
          min_dist = tmp;
      }

    return std::sqrt(min_dist);
  }
}

// Calculate the voronoi plane of feature A with one of its neighbor features B
// Any point that has a negative distance to the plane is closer to B than it is 
// to A
// if feature A and feature B are swapped, the normal of the plane should flip 
boost::shared_ptr<Plane> Polyhedron::voronoi_plane (FeatureType fA, FeatureType fB, boost::shared_ptr<const Ravelin::Pose3d> pose, boost::shared_ptr<const Polyhedron::Feature>& featureA, boost::shared_ptr<const Polyhedron::Feature>& featureB)
{
  boost::shared_ptr<Plane> result;

  // handle VP(V,E) case
  if (fA==eVertex && fB==eEdge)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Vertex> vertexA = boost::static_pointer_cast<const Polyhedron::Vertex>(featureA);
    boost::shared_ptr<const Polyhedron::Edge> edgeB = boost::static_pointer_cast<const Polyhedron::Edge>(featureB);



    //creating the vector for the vertex and a random end of the second edge
    Ravelin::Vector3d pA(vertexA->o,pose);
    Ravelin::Vector3d pB(edgeB->v1->o,pose);
    
    // if the points are the same, create a new point to address normalization
    if((pA-pB).norm() < NEAR_ZERO)
      pB=Ravelin::Vector3d(edgeB->v2->o,pose);

    // calculate the normal
    Ravelin::Vector3d normal = pA - pB;
    normal.normalize();

    // setup the plane
    result = boost::shared_ptr<Plane>( new Plane(normal,pA));
  }
  // handle VP(E,V) case
  else if (fA==eEdge && fB==eVertex)
  {
    result = voronoi_plane(fB,fA,pose,featureB,featureA);

    result->set_normal(-(result->get_normal()));
    result->offset = -(result->offset);
  }
  // handle VP(E,F) case
  else if (fA==eEdge && fB==eFace)
  {
    // cast pointers
    boost::shared_ptr<const Polyhedron::Edge> edgeA = boost::static_pointer_cast<const Polyhedron::Edge>(featureA);
    boost::shared_ptr<const Polyhedron::Face> faceB_const = boost::static_pointer_cast<const Polyhedron::Face>(featureB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(faceB_const);

    // find the normal vector in the global frame of the plane the face is on 
    Ravelin::Vector3d face_normal = faceB->get_plane().get_normal();

    // find the vector representing the edge
    Ravelin::Vector3d v1(edgeA->v1->o, pose);
    Ravelin::Vector3d v2(edgeA->v2->o, pose);
    Ravelin::Vector3d edge_vector = v2 - v1;
    face_normal.pose = pose;
    // find the cross product between them to find the Voronoi plane normal
    Ravelin::Vector3d voronoi_normal = Ravelin::Vector3d::cross(edge_vector, face_normal);
    voronoi_normal.normalize();
    FILE_LOG(LOG_COLDET) << "v1: "<< v1 << std::endl << *faceB << std::endl << voronoi_normal << std::endl;
    // create the Voronoi plane    
    result = boost::shared_ptr<Plane>( new Plane(voronoi_normal,v1));

    // Since we don't know whether the normal is in the correct orientation,
    // we have to test it by plugging in one vertex of the face that is not the two vertices of the edge to test it out
    
    // create the VertexFaceIterator
    VertexFaceIterator vfi(faceB, true);

    // find a point that is not v1 or v2 
    boost::shared_ptr<Polyhedron::Vertex> test_vert = *vfi;
    Ravelin::Vector3d test_vect(test_vert->o, pose);

    if ((test_vect - v1).norm() < NEAR_ZERO ||
        (test_vect - v2).norm() < NEAR_ZERO)
    {
      while(vfi.has_next())
      {
        vfi.advance();
        test_vect= Ravelin::Vector3d((*vfi)->o, pose);
        if ((test_vect - v1).norm() > NEAR_ZERO &&
            (test_vect - v2).norm() > NEAR_ZERO)
        {
          test_vert = *vfi;
          break;
        }
      }
    }
    
    // the signed distance should be negative, so if the signed distance is positive, we have to reverse the vector
    if (result->calc_signed_distance(Ravelin::Vector3d(test_vert->o, pose)) > NEAR_ZERO)
    {
      result->set_normal(-(result->get_normal()));
      result->offset = -(result->offset);
    }
  }
  else if (fA==eFace && fB==eEdge)
  {
    // cast pointers
    result = voronoi_plane(fB,fA,pose,featureB,featureA);

    result->set_normal(-(result->get_normal()));
    result->offset = -(result->offset);
  }
  FILE_LOG(LOG_COLDET) << *result << std::endl;

  return result;
}


// Implements Mirtich's clipEdge algorithm
// The algorithm takes in an edge E from a polyhedron and multiple voronoi plane from another polyhedron
// and returns if the the edge intersects with the plane
// if the edge is intersected, the function also returns the portion of the edge the intersection is at
// and which neighbor feature the voronoi plane represents.
bool Polyhedron::clip_edge(boost::shared_ptr<const Polyhedron::Edge> edge, Transform3d fTe, double& min_lambda, double& max_lambda, boost::shared_ptr<const Polyhedron::Feature >& min_N, boost::shared_ptr<const Polyhedron::Feature >& max_N, const std::list<std::pair< boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > >& planes_neighbors) 
{  
  Ravelin::Vector3d t_e(edge->v1->o, fTe.source);
  Ravelin::Vector3d h_e(edge->v2->o, fTe.source);
  Ravelin::Vector3d t = fTe.transform_point(t_e);
  Ravelin::Vector3d h = fTe.transform_point(h_e);
  
  FILE_LOG(LOG_COLDET)<< "fTe " << fTe <<std::endl;
  FILE_LOG(LOG_COLDET)<< "t_e: " << t_e <<std::endl;
  FILE_LOG(LOG_COLDET)<< "h_e: " << h_e <<std::endl;
  FILE_LOG(LOG_COLDET)<< "t: " << t <<std::endl;
  FILE_LOG(LOG_COLDET)<< "h: " << h <<std::endl;

  //iterating through the pair list
  std::list< std::pair< boost::shared_ptr<const Polyhedron::Feature >, boost::shared_ptr<Plane> > >::const_iterator pni;
  for (pni=planes_neighbors.begin();pni!=planes_neighbors.end();++pni)
  {
    //getting the neighbor feature and the v-plane
    boost::shared_ptr<const Polyhedron::Feature> N=(*pni).first;
    boost::shared_ptr<Plane> P=(*pni).second;
    double lambda;
    FILE_LOG(LOG_COLDET)<< "Edge clipping plane: " << *P <<std::endl;
    //calculate the distance from the two end of the edge to the plane
    double dt = P->calc_signed_distance(t);
    double dh = P->calc_signed_distance(h);
    FILE_LOG(LOG_COLDET)<< "dt: " << dt <<std::endl;
    FILE_LOG(LOG_COLDET)<< "dh: " << dh <<std::endl;

    //if the edge is completely clipped
    if (dt<-NEAR_ZERO && dh<-NEAR_ZERO)
    {
      min_N=max_N=N;
      return false;
    }
    //if only one side is clipped
    else if(dt<-NEAR_ZERO)
    {  
      //Find lambda
      lambda = dt/(dt-dh);

      //if the current lambda is larger than the minimum possible lambda
      //We have to update the minimum lambda
      if (lambda-min_lambda > NEAR_ZERO)
      {
        min_lambda = lambda;
        min_N = N;
        //if the edge is completely clipped
        if (min_lambda-max_lambda > NEAR_ZERO)
          return false;
      }
    }
    else if (dh<-NEAR_ZERO)
    {
      // find lambda
      lambda = dt/(dt-dh);

      // if the current lambda is smaller than the maximum possible lambda
      // we have to update the minimum lambda
      if(lambda-max_lambda < -NEAR_ZERO)
      {
        max_lambda = lambda;
        max_N = N;

        //if the edge is completely clipped
        if(max_lambda - min_lambda < -NEAR_ZERO)
          return false;
      }
    }
  }
  return true;
}

// Does a post-clipping derivative check
// This test is executed after an edge is clipped by voronoi regions.
// The test finds whether the neighbour feature of X is closer to the edge than X is
// In order to do this, this method examine the derivative of distance wrt to lambda,
// the portion of the edge, when lambda = max_lambda or min_lambda. If the derivative is 
// positive at min_lambda or negative at max_lambda. It means the neighbor feature is
// closer. True is returned if one neighbour feature is closer, false is returned otherwise.

bool Polyhedron::post_clip_deriv_check(FeatureType& fX, boost::shared_ptr<const Polyhedron::Feature >& X , boost::shared_ptr<const Polyhedron::Edge> edge, Transform3d& xTe, double& min_lambda, double& max_lambda, boost::shared_ptr<const Polyhedron::Feature >& min_N, boost::shared_ptr<const Polyhedron::Feature >& max_N)
{

  boost::shared_ptr<Ravelin::Pose3d> GLOBAL3D;

  // create a vector pointing of the edge from tail to head(v1 to v2)
  Ravelin::Vector3d t_e(edge->v1->o, xTe.source);
  Ravelin::Vector3d h_e(edge->v2->o, xTe.source);
  Ravelin::Vector3d t = xTe.transform_point(t_e);
  Ravelin::Vector3d h = xTe.transform_point(h_e);
  Ravelin::Vector3d u = h - t;

  // placeholders for the sign of the derivative at min_lambda and max_lambda 
  double Ddot_min,Ddot_max;

  // Vertex case
  if (fX == eVertex)
  {
    boost::shared_ptr<const Polyhedron::Vertex > vX = boost::static_pointer_cast<const Polyhedron::Vertex>(X);
    Ravelin::Vector3d v(vX->o, xTe.target);
    Ddot_min = Ravelin::Vector3d::dot(u,(t+min_lambda*u-v));
    Ddot_max = Ravelin::Vector3d::dot(u,(t+max_lambda*u-v));
  }
  // Edge case
  else if (fX == eEdge)
  {
    // in Edge case, we evaluate the dervative wrt the neighbor feature instead of X

    if(min_N)
    {// calculate dDot_min
      if(boost::dynamic_pointer_cast<const Polyhedron::Vertex>(min_N))
      {
        boost::shared_ptr<const Polyhedron::Vertex > vX = boost::static_pointer_cast<const Polyhedron::Vertex>(min_N);
        Ravelin::Vector3d v(vX->o, xTe.target);
        Ddot_min = Ravelin::Vector3d::dot(u,(t+min_lambda*u-v));
      }
      // if the neighbor is not a vertex, it must be a face
      else
      {
        boost::shared_ptr<const Polyhedron::Face> faceX =  boost::static_pointer_cast<const Polyhedron::Face>(min_N);
        Plane p = faceX->get_plane();
        Ravelin::Vector3d n = p.get_normal();
        n.pose = xTe.target;

        // calculate dDot_min
        Ddot_min = Ravelin::Vector3d::dot(u,n);
        Ravelin::Vector3d v = t+u*min_lambda;
        v.pose = GLOBAL3D;
        // if the signed distance is negative, we need to reverse the sign
        if(p.calc_signed_distance(v)<-NEAR_ZERO)
        {
          Ddot_min = -Ddot_min;
        }
      }
    }

    // calculate dDot_max
    if(max_N)
    {
      if(boost::dynamic_pointer_cast<const Polyhedron::Vertex>(max_N))
      {
        boost::shared_ptr<const Polyhedron::Vertex > vX = boost::static_pointer_cast<const Polyhedron::Vertex>(max_N);
        Ravelin::Vector3d v(vX->o, xTe.target);
        Ddot_max = Ravelin::Vector3d::dot(u,(t+max_lambda*u-v));
      }
      // if the neighbor is not a vertex, it must be a face
      else
      {
        boost::shared_ptr<const Polyhedron::Face> faceX =  boost::static_pointer_cast<const Polyhedron::Face>(max_N);
        Plane p = faceX->get_plane();
        Ravelin::Vector3d n = p.get_normal();
        n.pose = xTe.target;

        // calculate dDot_min
        Ddot_max = Ravelin::Vector3d::dot(u,n);

        Ravelin::Vector3d v = t+u*max_lambda;

        v.pose = GLOBAL3D;
        // if the signed distance is negatve, we need to reverse the sign
        if(p.calc_signed_distance(v)<-NEAR_ZERO)
          Ddot_max = -Ddot_max;
      }
    }
  }
  // face case
  else if (fX == eFace)
  {
    boost::shared_ptr<const Polyhedron::Face> fX =  boost::static_pointer_cast<const Polyhedron::Face>(X);
    Plane p = fX->get_plane();

    // calculate dDot_min
    Ravelin::Vector3d n = p.get_normal();

    n.pose = xTe.target;

    // calculate dDot_min
    Ddot_min = Ravelin::Vector3d::dot(u,n);

    Ravelin::Vector3d vx = t+u*min_lambda;
    vx.pose = GLOBAL3D;
    // if the signed distance is negatve, we need to reverse the sign
    if(p.calc_signed_distance(vx)<-NEAR_ZERO)
      Ddot_min = -Ddot_min;

    // calculate dDot_max
    Ddot_max = Ravelin::Vector3d::dot(u,n);

    vx = t+u*max_lambda;
    vx.pose = GLOBAL3D;
    // if the signed distance is negatve, we need to reverse the sign
    if(p.calc_signed_distance(vx)<-NEAR_ZERO)
      Ddot_max = -Ddot_max;
  }

  FILE_LOG(LOG_COLDET) << "Ddot_max: " << Ddot_max <<std::endl;
  FILE_LOG(LOG_COLDET) << "Ddot_min: " << Ddot_min <<std::endl;
  FILE_LOG(LOG_COLDET) << "max: " << max_N <<std::endl;
  FILE_LOG(LOG_COLDET) << "maxV: " << boost::dynamic_pointer_cast<const Polyhedron::Vertex>(max_N) <<std::endl;
  FILE_LOG(LOG_COLDET) << "maxF: " << boost::dynamic_pointer_cast<const Polyhedron::Face>(max_N) <<std::endl;

  FILE_LOG(LOG_COLDET) << "min: " << min_N <<std::endl;
  FILE_LOG(LOG_COLDET) << "minV: " << boost::dynamic_pointer_cast<const Polyhedron::Vertex>(min_N) <<std::endl;
  FILE_LOG(LOG_COLDET) << "minF: " << boost::dynamic_pointer_cast<const Polyhedron::Face>(min_N) <<std::endl;


  // check whether it is posible to update X
  // the neighbour feature is closer iff the derivative is positive
  // at minimum lambda or negative at maximum lambda
  if(min_N && Ddot_min>NEAR_ZERO)
  {
    X = min_N;
    if (boost::dynamic_pointer_cast<const Polyhedron::Vertex>(min_N))
      fX = eVertex;  
    else if (boost::dynamic_pointer_cast<const Polyhedron::Face>(min_N))
      fX = eFace;
    else
      fX = eEdge;
    return true; 
  }
  else if(max_N && Ddot_max<-NEAR_ZERO)
  {
    X = max_N;
    if (boost::dynamic_pointer_cast<const Polyhedron::Vertex>(max_N))
      fX = eVertex;  
    else if (boost::dynamic_pointer_cast<const Polyhedron::Face>(max_N))
      fX = eFace;
    else 
      fX = eEdge;
    return true; 
  }

  return false;
}

// test if the vertex is inside the polygon
// if not, update the face to the face that is the most distant from vertex
Polyhedron::UpdateRule Polyhedron::handle_local_minimum(boost::shared_ptr<const Polyhedron::Vertex>& V, FeatureType& fF, boost::shared_ptr<const Polyhedron::Feature>& face,  const Polyhedron& face_poly, const Ravelin::Transform3d& fTv)
{
  // global pose
  boost::shared_ptr<Ravelin::Pose3d> GLOBAL3D; 

  // check whether the vertex has a negative distance with all faces in face_poly
  double d_max = -std::numeric_limits<double>::max();
  const std::vector<boost::shared_ptr<Face> >& faces = face_poly.get_faces();
  std::vector<boost::shared_ptr<Face> >::const_iterator fi;
  Ravelin::Vector3d v(V->o, fTv.source);
  Ravelin::Vector3d vf = fTv.transform_point(v);
  boost::shared_ptr<const Polyhedron::Feature> f0;

  // set the frame of vf to global for calculation
  vf.pose = GLOBAL3D;

  for(fi = faces.begin(); fi != faces.end() ; ++fi)
  {
    boost::shared_ptr<const Face> f = *fi;
    Plane p = f->get_plane();
    double d = p.calc_signed_distance(vf);
    FILE_LOG(LOG_COLDET) << "The distance between the vertex "<< vf << " and face " << *f << " is " << d << std::endl;
    if (d - d_max > NEAR_ZERO)
    {
      d_max = d;
      f0 = f;
    }
  }
  // check ends

  if (d_max < NEAR_ZERO)
    return eInterpenetrating;

  face = f0;
  fF = eFace;
  return eContinue;
}

/// update rule for vertex/vertex case
Polyhedron::UpdateRule Polyhedron::update_vertex_vertex(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{

  // search for Voronoi plane from three coincident edges to vertex A,
  // which vertex B violates
  boost::shared_ptr<const Polyhedron::Vertex> vertA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
  boost::shared_ptr<const Polyhedron::Vertex> vertB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);

  std::list<boost::weak_ptr<Edge> > es = vertA->e;
  std::list<boost::weak_ptr<Edge> >::const_iterator ei;
  Ravelin::Vector3d vectorB_B(vertB->o,aTb.source);
  Ravelin::Vector3d vectorB = aTb.transform_point(vectorB_B);

  FeatureType fEdge=eEdge;
  FeatureType fVertex=eVertex;

  FILE_LOG(LOG_COLDET)<< vectorB <<std::endl;

  for(ei=es.begin(); ei!=es.end(); ++ei)
  {
    boost::shared_ptr<const Feature> e(*ei);
    boost::shared_ptr<Plane> vp = voronoi_plane(fVertex,fEdge,aTb.target,closestA,e);

    FILE_LOG(LOG_COLDET)<< *vp <<std::endl<<vp->calc_signed_distance(vectorB)<<std::endl;

    if (vp->calc_signed_distance(vectorB) < -NEAR_ZERO)
    {
      // vertex B violates plane from this edge; update vA to eA
      closestA = e;
      fA=eEdge;
      return eContinue;
    }
  }

  // search for Voronoi plane from three coincident edges to vertex B,
  // which vertex A violates
  es = vertB->e;
  Ravelin::Vector3d vectorA_A(vertA->o,aTb.target);
  Transform3d bTa = aTb.inverse();
  Ravelin::Vector3d vectorA = bTa.transform_point(vectorA_A);

  FILE_LOG(LOG_COLDET)<< vectorA <<std::endl;

  for(ei=es.begin(); ei!=es.end(); ++ei)
  {
    boost::shared_ptr<const Feature> e(*ei);
    boost::shared_ptr<Plane> vp = voronoi_plane(fVertex,fEdge,aTb.source,closestB,e);

    FILE_LOG(LOG_COLDET)<< *vp <<std::endl<<vp->calc_signed_distance(vectorA)<<std::endl;

    if(vp->calc_signed_distance(vectorA) < -NEAR_ZERO)
    {
      // vertex A violates plane from this edge; update vB to eA
      closestB = e;
      fB=eEdge;
      return eContinue;
    }
  }

  // still here = no violations
  return eDone;
}

/// update rule for vertex/edge case
Polyhedron::UpdateRule Polyhedron::update_vertex_edge(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  boost::shared_ptr<const Polyhedron::Vertex> vertA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
  boost::shared_ptr<const Polyhedron::Edge> edgeB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
  const FeatureType F_FACE=eFace;
  const FeatureType F_EDGE=eEdge;
  const FeatureType F_VERTEX=eVertex;  


  // search for Voronoi plane from those coincident to eB that vA violates
  Ravelin::Vector3d vectA_A(vertA->o,aTb.target);
  Ravelin::Vector3d vectA = aTb.inverse().transform_point(vectA_A);
  FILE_LOG(LOG_COLDET) << vectA << std::endl;
  //V(E,V1)
  boost::shared_ptr<const Feature>  N = edgeB->v1;
  boost::shared_ptr<Plane> vp = voronoi_plane(F_EDGE, F_VERTEX,aTb.source,closestB,N);

  FILE_LOG(LOG_COLDET)<< *(edgeB->v1) << std::endl << *vp<< std::endl << vp->calc_signed_distance(vectA) << std::endl;
  double dist = vp->calc_signed_distance(vectA);
  if( dist< -NEAR_ZERO || fabs(dist) < NEAR_ZERO )
  {
    // vA violates plane; update edge to coincident plane
    closestB = N;
    fB = eVertex;
    return eContinue;
  }


  //V(E,V2)
  N = edgeB->v2;
  vp = voronoi_plane(F_EDGE,F_VERTEX,aTb.source,closestB,N);

  FILE_LOG(LOG_COLDET)<< *vp<< std::endl << vp->calc_signed_distance(vectA) << std::endl;
  dist = vp->calc_signed_distance(vectA);
  if( dist< -NEAR_ZERO || fabs(dist) < NEAR_ZERO )
  {
    // vA violates plane; update edge to coincident plane
    closestB = N;
    fB = eVertex;
    return eContinue;
  }

  // V(E,F1)
  N = edgeB->face1;
  vp = voronoi_plane(F_EDGE,F_FACE,aTb.source,closestB,N);

  FILE_LOG(LOG_COLDET)<< *vp<< std::endl << vp->calc_signed_distance(vectA) << std::endl;

  if(vp->calc_signed_distance(vectA) < -NEAR_ZERO)
  {
    // vA violates plane; update edge to coincident plane
    closestB = N;
    fB = eFace;
    return eContinue;
  }  

  // V(E,F1)
  N = edgeB->face2;
  vp = voronoi_plane(F_EDGE,F_FACE,aTb.source,closestB,N);

  FILE_LOG(LOG_COLDET)<< *vp<< std::endl << vp->calc_signed_distance(vectA) << std::endl;

  if (vp->calc_signed_distance(vectA) < -NEAR_ZERO)
  {
    // vA violates plane; update edge to coincident plane
    closestB = N;
    fB = eFace;
    return eContinue;
  }  

  FILE_LOG(LOG_COLDET)<< "Start clipping the edge"<< std::endl;
  // clip eB against the Voronoi region of vA
  double min_lambda=0;
  double max_lambda=1;
  boost::shared_ptr<const Polyhedron::Feature > min_N, max_N;
  
  // create vp-feature pair list
  std::list<std::pair<boost::shared_ptr<const Polyhedron::Feature> , boost::shared_ptr<Plane> > > planes_neighbors;
  std::list<boost::weak_ptr<Edge> > es = vertA->e;
  std::list<boost::weak_ptr<Edge> >::iterator ei;

  FeatureType fEdge=eEdge;
  FeatureType fVertex=eVertex;

  for(ei=es.begin(); ei!=es.end(); ++ei)
  {
    boost::shared_ptr<const Feature> e(*ei);
    boost::shared_ptr<Plane> const vp = voronoi_plane(fVertex,fEdge,aTb.target,closestA,e);
    std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > pn(e,vp);
    planes_neighbors.push_back(pn);
  }

  bool clip_result=clip_edge(edgeB, aTb, min_lambda, max_lambda, min_N, max_N, planes_neighbors);

  FILE_LOG(LOG_COLDET) << min_lambda << std::endl << max_lambda << std::endl;
  // check whether the edge is completely clipped by one feature
  if(min_N==max_N && min_N)
  {
    FILE_LOG(LOG_COLDET)<< "Completely Clipped"<< std::endl;
    closestA=min_N;
    fA=eEdge;
    return eContinue;
  }
  else
  {
    // check derivative and update the feature
    if(post_clip_deriv_check(fA, closestA, edgeB, aTb, min_lambda, max_lambda, min_N, max_N))
      return eContinue;
    // if V is not updated after all of porcess
    else
      return eDone;
  }
}

/// update rule for vertex/face case
Polyhedron::UpdateRule Polyhedron::update_vertex_face(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB, const Polyhedron& face_poly)
{
  // search for Voronoi plane from VP(F,E) that vertex A violate
  boost::shared_ptr<const Polyhedron::Vertex> vertA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
  boost::shared_ptr<const Polyhedron::Face> faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
  const FeatureType F_EDGE=eEdge;
  const FeatureType F_VERTEX=eVertex;
  const FeatureType F_FACE=eFace;

  std::list<boost::weak_ptr<Edge> > es = faceB->e;
  std::list<boost::weak_ptr<Edge> >::iterator ei;
  Ravelin::Vector3d vectorA(vertA->o,aTb.target);
  Ravelin::Vector3d vectorA_b = aTb.inverse().transform_point(vectorA);
  // since we are looking for violate distance, it will always be negative.
  double max_violate_distance = 0;
  boost::shared_ptr<const Polyhedron::Feature> max_violate_feature;

  // iterate through edge list
  for (ei=es.begin(); ei!=es.end(); ++ei)
  {
    boost::shared_ptr<const Feature> e(*ei);
    boost::shared_ptr<Plane> vp = voronoi_plane(F_FACE,F_EDGE,aTb.source,closestB,e);
    double dist = vp->calc_signed_distance(vectorA_b);
    boost::shared_ptr<const Edge> e_temp = boost :: static_pointer_cast<const Edge>(e);
    FILE_LOG(LOG_COLDET) << *e_temp << std::endl << *vp << std::endl<< "Vector: " << vectorA_b << std::endl << "dist: " << dist << std::endl;
    if ( dist < max_violate_distance || fabs(max_violate_distance-dist)<NEAR_ZERO)
    {
      // vertex B violates plane from this edge; update vA to eA
      max_violate_distance = dist;
      max_violate_feature = e;
    }
  }

  // if violated plane exists
  if (max_violate_feature)
  {
    closestB = max_violate_feature;
    fB = eEdge;
    return eContinue;
  }

  // check whether there are edges from V pointing toward F
  // Transforming plane B to A
  Plane p_b = faceB->get_plane();
  Ravelin::Vector3d normal_b = p_b.get_normal();
  normal_b.pose = aTb.source;
  Ravelin::Vector3d normal_a = aTb.transform_vector(normal_b);
  FILE_LOG(LOG_COLDET) << "The Face is " << *faceB <<std::endl << "Normal before transformation is " <<normal_b <<endl << "normal after transformation is " << normal_a << std::endl;

  // find a random vertex from face b
  boost::shared_ptr<Polyhedron::Face> faceB_non_const = boost::const_pointer_cast<Polyhedron::Face>(faceB);
  Polyhedron::VertexFaceIterator vfi(faceB_non_const, true);
  shared_ptr<Polyhedron::Vertex> v1 = *vfi;
  Ravelin::Vector3d v1_b (v1->o, aTb.source);
  Ravelin::Vector3d v1_a = aTb.transform_point(v1_b);
  normal_a.normalize();
  Plane p_a(normal_a, v1_a);

  double D_va = p_a.calc_signed_distance(vectorA);

  FILE_LOG(LOG_COLDET) << v1_b << " , " << v1_a <<std::endl << vectorA << std::endl << p_a << std::endl << D_va << std::endl;


  es = vertA->e;

  FILE_LOG(LOG_COLDET) << "edge list size: " <<es.size() <<std::endl;

  for(ei=es.begin(); ei!=es.end(); ++ei)
  {
    //find the vertex V' of e that is not vertA
    boost::shared_ptr<const Edge> e(*ei);
    boost::shared_ptr<Polyhedron::Vertex> v_prime = e->v1;

    if(vertA == v_prime)
      v_prime = e->v2;

    // calculate D(V')
    Ravelin::Vector3d vectorA_prime(v_prime->o, aTb.target);
    double D_v_prime = p_a.calc_signed_distance(vectorA_prime);
    FILE_LOG(LOG_COLDET) << "v_prime: " << *v_prime << std::endl << "D_v_prime: " << D_v_prime << std::endl << "D_va: " << D_va <<std::endl;
    if( D_va*(D_va - D_v_prime) > NEAR_ZERO)
    {
      // vertex B violates plane from this edge; update vA to eA
      closestA = e;
      fA=eEdge;
      return eContinue;
    }    
  }
  FILE_LOG(LOG_COLDET) << "edge test complete" << std::endl;
  // if there are no edge pointing to the face and the
  // vertex is absolutely out of the polyhedron. We are done 
    if (D_va > NEAR_ZERO)
      return eDone;
  // if the point has negative distance, it might be inside the polyhedron
  // handle_local_minimum will be used to test this and update.
  Transform3d bTa = aTb.inverse();
  return handle_local_minimum(vertA, fB, closestB, face_poly, bTa);
}

/// update rule for edge/edge case
Polyhedron::UpdateRule Polyhedron::update_edge_edge(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  
  //Casting
  boost::shared_ptr<const Polyhedron::Edge> edgeA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
  boost::shared_ptr<const Polyhedron::Edge> edgeB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
  const FeatureType F_FACE=eFace;
  const FeatureType F_VERTEX=eVertex;
  const FeatureType F_EDGE=eEdge;  

  {
    // clip eB against the Edge-Vertex Voronoi plane of EA
    double min_lambda=0;
    double max_lambda=1;
    boost::shared_ptr<const Polyhedron::Feature > min_N , max_N;
  
    // creating vp-feature pair list
    std::list<std::pair<boost::shared_ptr<const Polyhedron::Feature> , boost::shared_ptr<Plane> > > planes_neighbors;

    // creating Edge vs. Vertex1 vp and adding to the list
    boost::shared_ptr<const Polyhedron::Feature> V = edgeA->v1;
    boost::shared_ptr<Plane> vp = voronoi_plane(F_EDGE, F_VERTEX, aTb.target, closestA, V);
    std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > pn(V , vp);
    planes_neighbors.push_back(pn);

    // create edge vs. vertex2 vp and add to the list
    V = edgeA->v2;
    vp = voronoi_plane(F_EDGE, F_VERTEX, aTb.target, closestA, V);
    pn = std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> >(V , vp);
    planes_neighbors.push_back(pn);

    // clip edge
    bool clip_result = clip_edge(edgeB, aTb, min_lambda, max_lambda, min_N, max_N, planes_neighbors);

    // logging clip result for debugging purposes
    FILE_LOG(LOG_COLDET) << "min_lambda: " << min_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_lambda: " << max_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "min_N: " << min_N <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_N: " << max_N <<std::endl;


    // check if the edge is completely clipped by one features
    if(min_N==max_N && min_N)
    {
      closestA = min_N;
      fA = eVertex;
      return eContinue;
    }
    else
    {
      // check derivative and update the feature
      if (post_clip_deriv_check(fA, closestA, edgeB, aTb, min_lambda, max_lambda, min_N, max_N))
        return eContinue;
    }
  
    // emptying the list
    planes_neighbors.clear();

    // create edge vs. face1 vp and add it to the list
    boost::shared_ptr<const Polyhedron::Feature> F = edgeA->face1;
    vp = voronoi_plane(F_EDGE,F_FACE,aTb.target,closestA,F);
    pn = std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> >(F,vp);
    planes_neighbors.push_back(pn);

    // create edge vs. face2 vp and add it to the list
    F = edgeA->face2;
    vp = voronoi_plane(F_EDGE,F_FACE,aTb.target,closestA,F);
    pn = std::pair<boost::shared_ptr<const  Polyhedron::Feature>, boost::shared_ptr<Plane> >(F,vp);
    planes_neighbors.push_back(pn);

    clip_result=clip_edge(edgeB, aTb, min_lambda, max_lambda, min_N, max_N, planes_neighbors);
    FILE_LOG(LOG_COLDET) << "min_lambda: " << min_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_lambda: " << max_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "min_N: " << min_N <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_N: " << max_N <<std::endl;

    // check whether the edge is completely clipped by one feature
    if(min_N==max_N && min_N)
    {
      closestA = min_N;
      fA=eFace;
      return eContinue;
    }
    else
    {
      // check derivative and update the feature
      if (post_clip_deriv_check(fA, closestA, edgeB, aTb, min_lambda, max_lambda, min_N, max_N))
        return eContinue;
    }
  }  

  // swap role of A and B
  Transform3d bTa = aTb.inverse();
  {
    // clip eA against the Edge-Vertex Voronoi plane of EB
    double min_lambda = 0;
    double max_lambda = 1;
    boost::shared_ptr<const Polyhedron::Feature > min_N, max_N;
  
    // creating vp-feature pair list
    std::list<std::pair<boost::shared_ptr<const Polyhedron::Feature> , boost::shared_ptr<Plane> > > planes_neighbors;

    // create edge vs. vertex1 and add it to the list
    boost::shared_ptr<const Polyhedron::Feature> V = edgeB->v1;
    boost::shared_ptr<Plane> vp = voronoi_plane(F_EDGE, F_VERTEX, aTb.source, closestB, V);
    std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > pn(V , vp);
    planes_neighbors.push_back(pn);

    // create edge vs. vertex2 and add it to the list
    V = edgeB->v2;
    vp = voronoi_plane(F_EDGE, F_VERTEX, aTb.source, closestB, V);
    pn = std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > (V , vp);
    planes_neighbors.push_back(pn);

    bool clip_result = clip_edge(edgeA, bTa, min_lambda, max_lambda, min_N, max_N, planes_neighbors);
    FILE_LOG(LOG_COLDET) << "min_lambda: " << min_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_lambda: " << max_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "min_N: " << min_N <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_N: " << max_N <<std::endl;

    // check whether the edge is completely clipped by one feature
    if (min_N==max_N && min_N)
    {
      closestB = min_N;
      fB = eVertex;
      return eContinue;
    }
    else
    {
      //Check derivative and update the feature
      if(post_clip_deriv_check(fB, closestB, edgeA, bTa, min_lambda, max_lambda, min_N, max_N))
        return eContinue;
    }

    // clip eB against the Edge-Face Voronoi plane of EA
    // create vp-feature pair list
    planes_neighbors.clear();

    // create edge vs. vertex1 and add it to the list
    boost::shared_ptr<const Polyhedron::Feature> F = edgeB->face1;
    vp = voronoi_plane(F_EDGE,F_FACE,aTb.source,closestB,F);
    pn = std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> >(F,vp);
    planes_neighbors.push_back(pn);

    // create edge vs. vertex2 and add it to the list
    F = edgeB->face2;
    vp = voronoi_plane(F_EDGE,F_FACE,aTb.source,closestB,F);
    pn = std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> >(F,vp);
    planes_neighbors.push_back(pn);

    clip_result=clip_edge(edgeA, bTa, min_lambda, max_lambda, min_N, max_N, planes_neighbors);
    FILE_LOG(LOG_COLDET) << "min_lambda: " << min_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_lambda: " << max_lambda <<std::endl;
    FILE_LOG(LOG_COLDET) << "min_N: " << min_N <<std::endl;
    FILE_LOG(LOG_COLDET) << "max_N: " << max_N <<std::endl;

    // check whether the edge is completely clipped by one feature
    if(min_N==max_N && min_N)
    {
      closestB=min_N;
      fB=eFace;
      return eContinue;
    }
    else
    {
      // check derivative and update the feature
      if(post_clip_deriv_check(fB, closestB, edgeA, bTa, min_lambda, max_lambda, min_N, max_N))
        return eContinue;
    }
  } 
 
  return eDone;
}

/// Update rule for edge/face case
Polyhedron::UpdateRule Polyhedron::update_edge_face(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  // cast pointers
  boost::shared_ptr<const Polyhedron::Edge> edgeA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
  boost::shared_ptr<const Polyhedron::Face> faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
  const FeatureType F_FACE=eFace;
  const FeatureType F_VERTEX=eVertex;
  const FeatureType F_EDGE=eEdge; 
  boost::shared_ptr<Ravelin::Pose3d> GLOBAL3D; 

  Ravelin::Transform3d bTa = aTb.inverse(); 

  std::list<std::pair<boost::shared_ptr<const Polyhedron::Feature> , boost::shared_ptr<Plane> > > planes_neighbors;
  std::list<boost::weak_ptr<Edge> > es = faceB->e;
  std::list<boost::weak_ptr<Edge> >::const_iterator ei;

  for (ei=es.begin(); ei!=es.end(); ++ei)
  {
    boost::shared_ptr<const Feature> e(*ei);
    boost::shared_ptr<Plane> const vp = voronoi_plane(F_FACE,F_EDGE,aTb.source,closestB,e);
    std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > pn(e,vp);
    planes_neighbors.push_back(pn);
  }

  // clip edge
  double min_lambda = 0;
  double max_lambda = 1;
  boost::shared_ptr<const Polyhedron::Feature > min_N, max_N;
  bool clip_result = clip_edge(edgeA, bTa, min_lambda, max_lambda, min_N, max_N, planes_neighbors);

  FILE_LOG(LOG_COLDET)<< "min_lambda: " <<min_lambda <<std::endl;
  FILE_LOG(LOG_COLDET)<< "max_lambda: " <<max_lambda <<std::endl;  
  FILE_LOG(LOG_COLDET)<< "Face-Edge clipped" <<std::endl;  

  // check whether the edge is completely clipped
  if(!clip_result)
  {
    //if it is, we need find the closest feature of the face to the edge
    FILE_LOG(LOG_COLDET)<< "Edge completely excluded, looking for the closest feature of the face" <<std::endl;  
  
    //Initializing the start edge and the iterator
    boost::shared_ptr<const Polyhedron::Feature > cur_feature;
    boost::shared_ptr<const Polyhedron::Feature > prev_feature;

   /* 
    * hueristic: choose min_N or max_N, based on which
    * corresponding region contains more of edge being clipped.
    * NOTE: this heuristic is presented in the paper; however, after
    * testing, it seams like the heuristic is causing the algorithm to 
    * stuck in a infinite loop, and therefore not used.
    */
  
    // if(min_lambda+max_lambda>1.0)
    //   cur_feature = min_N;
    // else
    //   cur_feature = max_N;

    // Randomly choose from min_N and max_N (if both exist)
    if(min_N && max_N)
    {
      int n_choice = rand() % 2;
      if(n_choice)
        cur_feature = min_N;
      else
        cur_feature = max_N;
    }
    else if(!min_N)
    {
      cur_feature = max_N;
    }else
    {
      cur_feature = min_N;
    }

    // find the 
    for (ei = es.begin(); ei != es.end(); ++ei)
    {
      boost::shared_ptr<const Polyhedron::Edge > e(*ei);
      if(e == cur_feature)
        break;
    }
    
    // update feature B to the closest edge or vertex on F to E
    while (true)
    {

      // clip feature A with the two Edge-Vertex voronoi-plane of cur_edge
      boost::shared_ptr<const Polyhedron::Feature > next_feature;
      boost::shared_ptr<const Polyhedron::Edge> cur_edge = boost::static_pointer_cast<const Polyhedron::Edge>(cur_feature);
      planes_neighbors.clear();
      boost::shared_ptr<const Polyhedron::Feature> tail = cur_edge->v1;
      boost::shared_ptr<const Polyhedron::Feature> head = cur_edge->v2;
      boost::shared_ptr<Plane> const vp_t = voronoi_plane(F_EDGE, F_VERTEX, aTb.source, cur_feature, tail);
      boost::shared_ptr<Plane> const vp_h = voronoi_plane(F_EDGE, F_VERTEX, aTb.source, cur_feature, head);
      std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > np_t(tail,vp_t);
      std::pair<boost::shared_ptr<const Polyhedron::Feature>, boost::shared_ptr<Plane> > np_h(head,vp_h);
      planes_neighbors.push_back(np_t);
      planes_neighbors.push_back(np_h);
      double min_lambda = 0;
      double max_lambda = 1;
      boost::shared_ptr<const Polyhedron::Feature > min_N, max_N;
      bool clip_result = clip_edge(edgeA, bTa, min_lambda, max_lambda, min_N, max_N, planes_neighbors);

      // check derivative 
      Ravelin::Vector3d t_A(edgeA->v1->o, bTa.source);
      Ravelin::Vector3d t = bTa.transform_point(t_A);
      Ravelin::Vector3d h_A(edgeA->v2->o, bTa.source);
      Ravelin::Vector3d h = bTa.transform_point(h_A);
      Ravelin::Vector3d u = h - t;

      // placeholders for the derivative minimum/maximum
      double Ddot_min,Ddot_max=0;

      // calculate dDot_min using vertex formula
      if(min_N)
      {
        boost::shared_ptr<const Polyhedron::Vertex > vB = boost::static_pointer_cast<const Polyhedron::Vertex>(min_N);
        Ravelin::Vector3d v(vB->o, bTa.target);
        Ddot_min = Ravelin::Vector3d::dot(u,(t+min_lambda*u-v));
      }
      if(max_N)
      {
        boost::shared_ptr<const Polyhedron::Vertex > vB = boost::static_pointer_cast<const Polyhedron::Vertex>(max_N);
        Ravelin::Vector3d v(vB->o, bTa.target);
        Ddot_max = Ravelin::Vector3d::dot(u,(t+max_lambda*u-v));
      }

      FILE_LOG(LOG_COLDET) << "cur_edge: " << *cur_edge <<std::endl;
      FILE_LOG(LOG_COLDET) << "min_lambda: " << min_lambda <<std::endl;
      FILE_LOG(LOG_COLDET) << "max_lambda: " << max_lambda <<std::endl;
      FILE_LOG(LOG_COLDET) << "min_N: " << min_N <<std::endl;
      FILE_LOG(LOG_COLDET) << "max_N: " << max_N <<std::endl;
      FILE_LOG(LOG_COLDET) << "Ddot_min: " << Ddot_min <<std::endl;
      FILE_LOG(LOG_COLDET) << "Ddot_max: " << Ddot_max <<std::endl;

      // check which Edge to update to
      if(min_N && Ddot_min > NEAR_ZERO)
      {
        // find the neighbor edge that has tail as its vertex
        // get the one of the two neighbor edges
        std::list<boost::weak_ptr<Edge> >::const_iterator ei2 = ei;
        ei2++;
        if(ei2 == es.end())
          ei2 = es.begin();
        boost::shared_ptr<const Edge> neighbor_e1(*ei2);

        boost::shared_ptr<const Polyhedron::Vertex> minV = boost::static_pointer_cast<const Polyhedron::Vertex>(min_N);
        FILE_LOG(LOG_COLDET) << "minV: " << *minV <<std::endl;
        
        FILE_LOG(LOG_COLDET) << "neighbor_e1: " << *neighbor_e1 <<std::endl;


        // if tail is in neighbor_e1
        if (neighbor_e1->v1 == min_N|| neighbor_e1->v2 == min_N)
        {
          next_feature = neighbor_e1;
          ei = ei2;
        }
        else
        {
          // pointing ei2 to the other neighbor edge;
          if (ei != es.begin())
          {
            ei2 = ei;
            ei2--;
          }
          else
          {
            ei2 = es.end();
            ei2--;
          }
          next_feature = boost::shared_ptr<const Edge>(*ei2);

          ei = ei2;
        }
        // The closest feature is a point
        if(next_feature == prev_feature)
        {
          closestB = min_N;
          fB = eVertex;
          return eContinue;
        }

        prev_feature = cur_feature;
        cur_feature = next_feature;
        
        // the current edge is not the shortest distance one, do everything again with the updated edge
        continue; 
      }
      else if (max_N && Ddot_max < -NEAR_ZERO)
      {
        // find the neighboring edge that has head as its vertex
        // get the two neighboring edges
        std::list<boost::weak_ptr<Edge> >::const_iterator ei2 = ei;
        ei2++;

        boost::shared_ptr<const Polyhedron::Vertex> maxV = boost::static_pointer_cast<const Polyhedron::Vertex>(max_N);
        FILE_LOG(LOG_COLDET) << "maxV: " << *maxV <<std::endl;
        

        if (ei2 == es.end())
          ei2 = es.begin();
        boost::shared_ptr<Edge> neighbor_e1(*ei2);

        FILE_LOG(LOG_COLDET) << "neighbor_e1: " << *neighbor_e1 <<std::endl;

        // if head is in neighbor_e1
        if ((neighbor_e1->v1 == max_N)||(neighbor_e1->v2 == max_N))
        {
          next_feature = neighbor_e1;
          ei = ei2;
        }
        else
        {
          // pointing ei2 to the other neighbor edge;
          if(ei != es.begin())
          {
            ei2 = ei;
            ei2--;
          }
          else
          {
            ei2 = es.end();
            ei2--;
          }
          next_feature = boost::shared_ptr<Edge>(*ei2);

          ei = ei2;
        }

        // The closest feature is a point
        if(next_feature == prev_feature)
        {
          closestB = max_N;
          fB = eVertex;
          return eContinue;
        }

        prev_feature = cur_feature;
        cur_feature = next_feature;


        // the current edge is not the shortest distance one, do everything again with the updated edge
        continue; 
      }

      // if we are here, we are at the shortest distance edge
      closestB = cur_feature;
      fB = eEdge;
      break;
    }

    return eContinue;
  }

  // not fully clipped
  Ravelin::Vector3d t_a(edgeA->v1->o, aTb.target);
  Ravelin::Vector3d h_a(edgeA->v2->o, aTb.target);
  Ravelin::Vector3d t = bTa.transform_point(t_a);
  Ravelin::Vector3d h = bTa.transform_point(h_a);
  Ravelin::Vector3d u = h - t;
  Plane p = faceB->get_plane();

  Ravelin::Vector3d min_vx = t+u*min_lambda;
  Ravelin::Vector3d max_vx = t+u*max_lambda;
  min_vx.pose = GLOBAL3D;
  max_vx.pose = GLOBAL3D;


  double min_d = p.calc_signed_distance(min_vx);
  double max_d = p.calc_signed_distance(max_vx);

  FILE_LOG(LOG_COLDET) << "min_vx: " << min_vx << std::endl << "max_vx: " << max_vx <<std::endl; 
  FILE_LOG(LOG_COLDET) << "min_d: " << min_d << std::endl << "max_d: " << max_d <<std::endl; 
 
  //check for interpenetration
  if (min_d*max_d < NEAR_ZERO)
    return eInterpenetrating;
  
  Ravelin::Vector3d n = p.get_normal();

  n.pose = aTb.source;
  // calculate dDot_min
  double dDot_min = Ravelin::Vector3d::dot(u,n);

  // if the signed distance is negatve, we need to reverse the sign
  if (min_d < -NEAR_ZERO)
    dDot_min=-dDot_min;

  // return the vertex that is closer to the face
  if (dDot_min >= -NEAR_ZERO)
  {
    if (min_N)
    {
      closestB = min_N;
      fB = eEdge;
    }
    else
    {
      closestA = edgeA->v1;
      fA = eVertex;
    }
  }
  else
  {
    if (max_N)
    {
      closestB = max_N;
      fB = eEdge;
    }
    else
    {
      closestA = edgeA->v2;
      fA = eVertex;
    }
  }
  return eContinue;
}
