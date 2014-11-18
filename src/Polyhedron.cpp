/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <map>
#include <Moby/Types.h>
#include <Moby/Log.h>
#include <Moby/CompGeom.h>
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

// macro for creating an edge or finding a previously created edge
#define CREATE_LOOKUP(vA, vB, eAB) { \
if ((vei = v_edges.find(std::make_pair(vA, vB))) != v_edges.end()) \
{ \
  eAB = vei->second; \
  if (cw) \
  { \
    assert(!eAB->faceR); \
    eAB->faceR = f; \
  } \
  else \
  { \
    assert(!eAB->faceL); \
    eAB->faceL = f; \
  } \
} \
else if ((vei = v_edges.find(std::make_pair(vB, vA))) != v_edges.end()) \
{ \
  eAB = vei->second; \
  if (cw) \
  { \
    assert(!eAB->faceL); \
    eAB->faceL = f; \
  } \
  else \
  { \
    assert(!eAB->faceR); \
    eAB->faceR = f; \
  } \
} \
else \
{ \
  eAB = boost::shared_ptr<Polyhedron::Edge>(new Polyhedron::Edge); \
  v_edges[std::make_pair(vA, vB)] = eAB; \
  if (cw) \
    eAB->faceR = f; \
  else \
    eAB->faceL = f; \
  eAB->v1 = vertex_map[vA->point]; \
  eAB->v2 = vertex_map[vB->point]; \
  poly._edges.push_back(eAB); \
  eAB->v1->e.push_back(eAB); \
  eAB->v2->e.push_back(eAB); \
} \
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
  if (normal.norm() < NEAR_ZERO)
    throw std::runtime_error("Tried to get plane containing degenerate face!");
  else
    normal.normalize();

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

/// Assignment operator
Polyhedron& Polyhedron::operator=(const Polyhedron& p)
{
  _bb_min = p._bb_min;
  _bb_max = p._bb_max;
  _convexity = p._convexity;
  _convexity_computed = p._convexity_computed;
  _edges = p._edges;
  _vertices = p._vertices;
  _faces = p._faces;

  return *this;
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

  // subtract vertices of B from vertices of A in the global frame
  const unsigned NVERTS = vA.size() * vB.size();
  vector<shared_ptr<Polyhedron::Vertex> > vnew(NVERTS);
  for (unsigned i=0, k=0; i< vA.size(); i++)
    for (unsigned j=0; j< vB.size(); j++)
    {
      vnew[k] = shared_ptr<Polyhedron::Vertex>(new Polyhedron::Vertex);
      vnew[k]->o = Origin3d(vA[i]) - Origin3d(vB[j]); 
      shared_ptr<std::pair<int, int> > intpair(new std::pair<int, int>(i, j));
      vnew[k]->data = intpair;
    }

  // ******************* compute convex hull begins *******************
  const unsigned X = 0, Y = 1, Z = 2;
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";  // TODO: remove Fx option?
  FILE* outfile, * errfile;

  // setup constants for qhull
  const int DIM = 3;
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

  // setup the points
  std::map<coordT*, boost::shared_ptr<Polyhedron::Vertex> > vertex_map;
  std::vector<coordT> qhull_points;
  qhull_points.resize(N_POINTS*DIM);
  coordT* points_begin = &qhull_points.front();
  unsigned j=0;
  for (unsigned i = 0, j=0; i< NVERTS; i++) 
  {
    vertex_map[points_begin+j] = vnew[i];
    qhull_points[j++] = vnew[i]->o[X];
    qhull_points[j++] = vnew[i]->o[Y];
    qhull_points[j++] = vnew[i]->o[Z];
  }

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif  

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code != 0)
  {
    // points are not collinear.. unsure of the error...
    FILE_LOG(LOG_COMPGEOM) << "Polyhedron::calc_convex_hull() - unable to execute qhull on points:" << std::endl;
    for (unsigned i=0; i< NVERTS; i++)
      FILE_LOG(LOG_COMPGEOM) << "  " << vnew[i] << std::endl;

    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // release the mutex, since we're not using qhull anymore
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

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

    // see how the facet is oriented
    bool cw = (facet->toporient ^ qh_ORIENTclock);

    // see whether the facet is simplicial- it changes how we must process
    // edges
    if (facet->simplicial)
    {
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
      f->e.push_back(e12);
      f->e.push_back(e23);
      f->e.push_back(e31);
      assert(e12 != e23 && e12 != e31 && e23 != e31);
    }
    else
    {
      // facet is non-simplicial; iterate over the "ridges" (edges)
      ridgeT* ridge;    // for iterating...
      ridgeT** ridgep;  // ...over ridges
      FOREACHridge_( facet->ridges )
      {
        // setup the edge
        boost::shared_ptr<Polyhedron::Edge> e;

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

          // setup the vertices
          e->v1 = vertex_map[vertex->point];
          vertex = *++vertex_pointer;
          e->v2 = vertex_map[vertex->point];
          assert(e->v1 != e->v2);

          // add edge to the vertices
          e->v1->e.push_back(e);
          e->v2->e.push_back(e);
        }
        else
          e = new_edge_iter->second;

        // setup face- we'll use the convention that faceR is qhull's top
        // and faceL is qhull's bottom
        if (ridge->top == facet)
          e->faceR = f;
        else
          e->faceL = f;

        // add the edge to the face
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
      FILE_LOG(LOG_COMPGEOM) << "examining edge " << emap[e] << std::endl;

      // loop over all other edges in the face
      BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we2, f->e)
      {
        // get the edge
        boost::shared_ptr<Polyhedron::Edge> e2(we2);

        // don't process same edge twice
        if (e == e2)
          continue;

        FILE_LOG(LOG_COMPGEOM) << "  against edge " << emap[e2] << std::endl;

        // look for edges matching up
        if (e->faceL == f)
        {
          if (e2->faceL == f)
          {
            if (e->v2 == e2->v1)
            {
              e->nextL = e2;
              e2->prevL = e;
            }
            else if (e->v1 == e2->v2)
            {
              e->prevL = e2;
              e2->nextL = e;
            }
          }
          else
          {
            if (e->v2 == e2->v2)
            {
              e->nextL = e2;
              e2->nextR = e;
            }
            else if (e->v1 == e2->v1)
            {
              e->prevL = e2;
              e2->prevR = e;
            }
          }
        }
        else if (e->faceR == f)
        {
          if (e2->faceR == f)
          {
            if (e->v2 == e2->v1)
            {
              e->nextR = e2;
              e2->prevR = e;
            }
            else if (e->v1 == e2->v1)
            {
              e->prevR = e2;
              e2->nextR = e;
            }
          }
          else if (e2->faceL == f)
          {
            if (e->v2 == e2->v2)
            {
              e->nextR = e2;
              e2->nextL = e;
            }
            else if (e->v1 == e2->v1)
            {
              e->prevR = e2;
              e2->prevL = e;
            }
          }
        } 
      }
    }      
  }

  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif 

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  // mark the polyhedron as convex (here and in calc_convex_hull())
  poly._convexity_computed = true;
  poly._convexity = -1.0;

  // calculate the axis-aligned bounding box  
  poly.calc_bounding_box();

  // ******************* compute convex hull ends *******************

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

/*
/// Finds the feature adjacent to this edge (including the edge itself) that is closest to the query point
void Polyhedron::Edge::find_closest_feature(const Origin3d& p, std::map<shared_ptr<Polyhedron::Feature>, double>& distances, double& closest_dist, std::list<shared_ptr<Polyhedron::Feature> >& closest_features)
{
  // compute the distance from the edge itself (if it hasn't already been computed)
  if (distances.find(shared_from_this()) == distances.end())
  {
  }

  // check left face
  if (distances.find(faceR) == distances.end())
    faceR->find_closest_feature(p, distances, closest_dist, closest_features);  

  // check right face 
  if (distances.find(faceL) == distances.end())
    faceL->find_closest_feature(p, distances, closest_dist, closest_features);

  
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
    e = shared_ptr<Polyhedron::Edge>((e->faceR == shared_from_this()) ? e->prevR : e->prevL);
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
        e = shared_ptr<Polyhedron::Edge>((e->faceR == shared_from_this()) ? e->prevR : e->prevL);  
        continue;
      }

      // check closest features from this edge 
      e->find_closest_feature(p, distances, closest_dist, closest_features);   

      // get the next edge
      e = shared_ptr<Polyhedron::Edge>((e->faceR == shared_from_this()) ? e->prevR : e->prevL);  
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
      e = shared_ptr<Polyhedron::Edge>((e->faceR == shared_from_this()) ? e->prevR : e->prevL);  
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
    e = shared_ptr<Polyhedron::Edge>((e->faceR == shared_from_this()) ? e->prevR : e->prevL);  
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
  for (unsigned i=0; i< edges.size(); i++)
    out << "edge " << i << ":  prevL: " << emap[shared_ptr<Polyhedron::Edge>(edges[i]->prevL)] << " nextL: " << emap[shared_ptr<Polyhedron::Edge>(edges[i]->nextL)] << "  prevR: " << emap[shared_ptr<Polyhedron::Edge>(edges[i]->prevR)] << "  nextR: " << emap[shared_ptr<Polyhedron::Edge>(edges[i]->nextR)] << std::endl;

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
      bool leftT = (e->faceR == f) ? false : true;

      if (leftT)
        out << "edge " << emap[e] << " (L): " << vmap[e->v1] << ", " << vmap[e->v2] << std::endl;
      else
        out << "edge  " << emap[e] << " (R): " << vmap[e->v1] << ", " << vmap[e->v2] << std::endl;
    }

    // now do a counter-clockwise visit around the face
    out << "ccw visit:";
    shared_ptr<Polyhedron::Edge> e(f->e.front());
    do
    {
      out << " " << emap[e];
      if (e->faceR == f)
        e = shared_ptr<Polyhedron::Edge>(e->prevR);
      else
        e = shared_ptr<Polyhedron::Edge>(e->nextL);
    }
    while (e != shared_ptr<Polyhedron::Edge>(f->e.front()));
    out << std::endl;

    // now do a clockwise visit around the face
    out << "cw visit:";
    e = shared_ptr<Polyhedron::Edge>(f->e.front());
    do
    {
      out << " " << emap[e];
      if (e->faceR == f)
        e = shared_ptr<Polyhedron::Edge>(e->nextR);
      else
        e = shared_ptr<Polyhedron::Edge>(e->prevL);
    }
    while (e != shared_ptr<Polyhedron::Edge>(f->e.front()));
    out << std::endl;
  }

  return out;
}

Polyhedron::EdgeFaceIterator::EdgeFaceIterator(shared_ptr<Polyhedron::Face> f)
{
  // save this face
  this->f = f;

  // pick the first edge of the face to start from
  term = e = shared_ptr<Polyhedron::Edge>(f->e.front());
}

/// Dereferences the iterator
shared_ptr<Polyhedron::Edge> Polyhedron::EdgeFaceIterator::operator*()
{
  return e;
}

/// Advances the iterator clockwise
void Polyhedron::EdgeFaceIterator::advance_cw()
{
  if (e->faceR == f)
    e = shared_ptr<Polyhedron::Edge>(e->nextR);
  else
    e = shared_ptr<Polyhedron::Edge>(e->prevL);
}

/// Advances the iterator counter-clockwise
void Polyhedron::EdgeFaceIterator::advance_ccw()
{
  if (e->faceL == f)
    e = shared_ptr<Polyhedron::Edge>(e->nextL);
  else
    e = shared_ptr<Polyhedron::Edge>(e->prevR);
}

/// Checks to see whether the iterator can be advanced clockwise
bool Polyhedron::EdgeFaceIterator::has_cw()
{
  return ((e->faceL == f && shared_ptr<Polyhedron::Edge>(e->prevL) != term) || 
          (e->faceR == f && shared_ptr<Polyhedron::Edge>(e->nextR) != term));
}

/// Checks to see whether the iterator can be advanced counter-clockwise
bool Polyhedron::EdgeFaceIterator::has_ccw()
{
  return ((e->faceL == f && shared_ptr<Polyhedron::Edge>(e->nextL) != term) || 
          (e->faceR == f && shared_ptr<Polyhedron::Edge>(e->prevR) != term));
}

/// Constructs a vertex-face iterator
Polyhedron::VertexFaceIterator::VertexFaceIterator(shared_ptr<Polyhedron::Face> f, bool ccw)
{
  // save this face
  this->f = f;

  // save whether it is counter-clockwise iteration
  this->ccw = ccw;

  // pick the edge of the face to start from
  e = term = shared_ptr<Polyhedron::Edge>(f->e.front());

  // determine the vertex from cw / ccw
  if (ccw)
  {
    if (e->faceL == f)
      v = e->v1;
    else
      v = e->v2;
  }
  else
  {
    if (e->faceR == f)
      v = e->v2;
    else
      v = e->v1;
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
  if (ccw)
  {
    if (e->faceL == f)
    {
      e = shared_ptr<Polyhedron::Edge>(e->nextL);
      v = (e->faceL == f) ? e->v1 : e->v2;
    }
    else
    {
      assert(e->faceR == f);
      e = shared_ptr<Polyhedron::Edge>(e->prevR);
      v = (e->faceL == f) ? e->v1 : e->v2;
    }
  }
  else
  {
    if (e->faceR == f)
    {
      e = shared_ptr<Polyhedron::Edge>(e->nextR);
      v = (e->faceL == f) ? e->v1 : e->v2;
    }
    else
    {
      e = shared_ptr<Polyhedron::Edge>(e->prevL);
      v = (e->faceL == f) ? e->v1 : e->v2;
    }
  }
}

/// Checks to see whether the iterator can be advanced clockwise
bool Polyhedron::VertexFaceIterator::has_next()
{
  if (ccw)
  {
    if (e->faceL == f)
    {
      if (shared_ptr<Polyhedron::Edge>(e->nextL) == term)
        return false;
    }
    else
    {
      if (shared_ptr<Polyhedron::Edge>(e->prevR) == term)
        return false;
    }
  }
  else
  {
    if (e->faceR == f)
    {
      if (shared_ptr<Polyhedron::Edge>(e->nextR) == term)
        return false;
    }
    else
    {
      if (shared_ptr<Polyhedron::Edge>(e->prevL) == term)
        return false;
    }
  }

  return true;
}

