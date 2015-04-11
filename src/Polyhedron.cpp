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
/*
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
    for (unsigned j=0; j< vB.size(); j++, k++)
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

        // setup face- we'll use the convention that face1 is qhull's top
        // and face2 is qhull's bottom
        if (ridge->top == facet)
          e->face1 = f;
        else
          e->face2 = f;

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

        // look for edges matching up
        if (e->face2 == f)
        {
          if (e2->face2 == f)
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
        else if (e->face1 == f)
        {
          if (e2->face1 == f)
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
          else if (e2->face2 == f)
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
*/

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

    // see whether v1 or v2 is in the next edge
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v2;
    else
      v = e->v1;
  }
  else
  {
    // get the next edge
    list<weak_ptr<Edge> >::const_reverse_iterator cw_iter2 = cw_iter;
    cw_iter2++;

    // get the vertices
    shared_ptr<Edge> e(*cw_iter);
    shared_ptr<Edge> e2(*cw_iter2);

    // see whether v1 or v2 is in the next edge
    if (e->v1 == e2->v1 || e->v1 == e2->v2)
      v = e->v2;
    else
      v = e->v1;
  }

  // set even
  even = true;
}

/// Dereferences the iterator
shared_ptr<Polyhedron::Vertex> Polyhedron::VertexFaceIterator::operator*()
{
  return v;
}

/// Advances the iterator clockwise
void Polyhedron::VertexFaceIterator::advance()
{
  // for "even" cases, we just need to advance the vertex
  if (even)
  {
    if (ccw)
    {
      shared_ptr<Edge> e(*ccw_iter);
      if (v == e->v1)
        v = e->v2;
      else
      {
        assert(v == e->v2);
        v = e->v1; 
      }
    }
    else
    {
      shared_ptr<Edge> e(*cw_iter);
      if (v == e->v1)
        v = e->v2;
      else
      {
        assert(v == e->v2);
        v = e->v1; 
      }
    }
  }
  else
  {
    // for "odd" cases, we need to update the edge
    // pick a vertex from the edge
    if (ccw)
    {
      // get the next edge
      list<weak_ptr<Edge> >::const_iterator ccw_iter2 = ccw_iter;
      ccw_iter2++;

      // get the vertices
      shared_ptr<Edge> e(*ccw_iter);
      shared_ptr<Edge> e2(*ccw_iter2);

      // see whether v1 or v2 is in the next edge
      if (e->v1 == e2->v1 || e->v1 == e2->v2)
        v = e->v2;
      else
        v = e->v1;

      // advance the edge iterator
      ccw_iter++;
    }
    else
    {
      // get the next edge
      list<weak_ptr<Edge> >::const_reverse_iterator cw_iter2 = cw_iter;
      cw_iter2++;

      // get the vertices
      shared_ptr<Edge> e(*cw_iter);
      shared_ptr<Edge> e2(*cw_iter2);

      // see whether v1 or v2 is in the next edge
      if (e->v1 == e2->v1 || e->v1 == e2->v2)
        v = e->v2;
      else
        v = e->v1;

      // advance the edge iterator
      cw_iter++;
    }
  }

  // alter the even flag
  even = !even;
}

/// Checks to see whether the iterator can be advanced clockwise
bool Polyhedron::VertexFaceIterator::has_next()
{
  // if the case is even, we always have another
  if (even)
    return true;

  // if not, see whether we are at the end of the list
  if (ccw)
    return (ccw_iter != f->e.end()); 
  else
    return (cw_iter != f->e.rend());
}

/// Executes the V-Clip algorithm on two polyhedra, determining closest features and signed distance
double Polyhedron::vclip(shared_ptr<const PolyhedralPrimitive> pA, shared_ptr<const PolyhedralPrimitive> pB, shared_ptr<const Pose3d> poseA, shared_ptr<const Pose3d> poseB, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  FeatureType fA, fB;

  // get the transformation between A and B, and vice versa
  Transform3d aTb = Pose3d::calc_relative_pose(poseB, poseA);
  Transform3d bTa = aTb.inverse(); 

  // TODO: if closest feature of A is null, pick features for A and B arbitrarily

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

  // iterate through the algorithm
  while (true)
  {
    // handle vertex/vertex case
    if (fA == eVertex && fB == eVertex)
    {
      Polyhedron::UpdateRule r = update_vertex_vertex(fA, fB, aTb, closestA, closestB);

      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    // handle vertex/edge cases
    else if (fA == eVertex && fB == eEdge)
    {
      Polyhedron::UpdateRule r = update_vertex_edge(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    else if (fB == eVertex && fA == eEdge)
    {
      Polyhedron::UpdateRule r = update_vertex_edge(fB, fA, bTa, closestB, closestA);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fB, fA, closestB, closestA, bTa);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    // handle edge/edge case
    else if (fA == eEdge && fB == eEdge)
    {
      Polyhedron::UpdateRule r = update_edge_edge(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    // handle edge/face cases
    else if (fA == eEdge && fB == eFace)
    {
      Polyhedron::UpdateRule r = update_edge_face(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    else if (fB == eEdge && fA == eFace)
    {
      Polyhedron::UpdateRule r = update_edge_face(fB, fA, bTa, closestB, closestA);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fB, fA, closestB, closestA, bTa);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
    // handle face/face case
    else if (fA == eFace && fB == eFace)
    {
      Polyhedron::UpdateRule r = update_face_face(fA, fB, aTb, closestA, closestB);
      
      // look for continuing to run algorithm
      if (r == eContinue)
        continue;

      // otherwise, we have converged
      double dist = calc_dist(fA, fB, closestA, closestB, aTb);

      if (r == eInterpenetrating)
        return -dist;
      else
        return dist; 
    }
  }
}

//double Polyhedron::calc_dist();

/// Computes the distance between two features
double Polyhedron::calc_dist(FeatureType fA, FeatureType fB, boost::shared_ptr<const Polyhedron::Feature> closestA, boost::shared_ptr<const Polyhedron::Feature> closestB, Ravelin::Transform3d& aTb)
{
	//case1: Vertex vs. Vertex
	if(fA == eVertex && fB == eVertex){
		
		//casting pointers
		boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
		boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
		
		//Creating vectors
		Ravelin::Vector3d vAa(vA->o , aTb.target);
		Ravelin::Vector3d vBb(vB->o , aTb.source);
		
		//Transforming B into frame A
		Ravelin::Vector3d vBa = aTb.transform_point(vBb);
		Ravelin::Vector3d diff = vAa-vBa;
		return diff.norm();

	}else if(fA == eVertex && fB == eEdge){

		//casting pinters
		boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
		boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);

		//Better if the following process is put into a function (maybe?)
		//Creating vectors for the vertex
		Ravelin::Vector3d vAa(vA->o , aTb.target);
		Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
		Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
		
		//transforming fB into frame A
		Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
		Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

  		//creating segment
 		LineSeg3 line(vB1a,vB2a);
 		
 		//place holders
 		Point3d p;
 		double t;

 		//computing distance
 		double dist = CompGeom::calc_dist(line,vAa,t,p);
 		return dist;
	}else if(fA == eEdge && fB == eVertex){
		//casting pinters
		boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
		boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);

		//Better if the following process is put into a function (maybe?)
		//Creating vectors for the vertex
		Ravelin::Vector3d vBb(vB->o , aTb.source);
		Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
		Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);
		
		//transforming fB into frame A
		Ravelin::Vector3d vBa = aTb.transform_point(vBb);

  		//creating segment
 		LineSeg3 line(vA1a,vA2a);
 		
 		//place holders
 		Point3d p;
 		double t;

 		//computing distance
 		double dist = CompGeom::calc_dist(line,vBa,t,p);
 		return dist;
	}else if(fA == eEdge && fB == eEdge){

		//casting pointers
		boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
		boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
		
		//creating vectors
		Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
		Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);

		Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
		Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
		
		//transforming fB into frame A
		Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
		Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

  		//creating segment
 		LineSeg3 lineA(vA1a,vA2a);
 		LineSeg3 lineB(vB1a,vB2a);

 		//place holder
 		Point3d p1;
 		Point3d p2;
 		double dist = CompGeom::calc_closest_points(lineA,lineB,p1,p2);
 		return dist;
	}else if(fA == eVertex && fB == eFace){

		//casting
		boost::shared_ptr<const Polyhedron::Vertex> vA = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
		boost::shared_ptr<const Polyhedron::Face> faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);

			//creating vector for A
		Ravelin::Vector3d vAa(vA->o , aTb.target);

    //creating vector for B
    std::list<boost::weak_ptr<Edge> > eBs=faceB->e;

    //minimum finding
    Plane planeB=faceB->get_plane();
    double min=planeB.calc_signed_distance(vAa);
		for(std::list<boost::weak_ptr<Edge> >::iterator eBsi = eBs.begin(); eBsi!=eBs.end();++eBsi){

      boost::shared_ptr<Edge> eB(*eBsi);
      Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
      Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
    
      //transforming fB into frame A
      Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
      Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

      //creating segment
      LineSeg3 line(vB1a,vB2a);
    
      //place holders
      Point3d p;
      double t;

      //computing distance and comparing to minimum
      double dist = CompGeom::calc_dist(line,vAa,t,p);
      if(min>dist){
        min=dist;
      }
    }
  return min;


	}else if(fA == eFace && fB == eVertex){

    //Casting pointers
    boost::shared_ptr<const Polyhedron::Face> faceA = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<const Polyhedron::Vertex> vB = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
  
    //transforming B
    Ravelin::Vector3d vBb(vB->o , aTb.source);
    Ravelin::Vector3d vBa = aTb.transform_point(vBb);

    std::list<boost::weak_ptr<Edge> > eAs=faceA->e;
    //minimum finding
    Plane planeA=faceA->get_plane();
    double min=planeA.calc_signed_distance(vBa);
    for(std::list<boost::weak_ptr<Edge> >::iterator eAsi = eAs.begin(); eAsi!=eAs.end();++eAsi){

      boost::shared_ptr<Edge> eA(*eAsi);
      Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
      Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);
      //creating segment
      LineSeg3 line(vA1a,vA2a);
    
    //place holders
      Point3d p;
      double t;

    //computing distance and comparing to minimum
      double dist = CompGeom::calc_dist(line,vBa,t,p);
      if(min>dist){
        min=dist;
      }
    }
  return min;
}else if(fA == eEdge && fB == eFace){

		//Casting pointers
		boost::shared_ptr<const Polyhedron::Edge> eA = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
		boost::shared_ptr<const Polyhedron::Face> _faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(_faceB);
  
		//creating lineseg for A
		Ravelin::Vector3d vA1a(Ravelin::Origin3d(eA->v1->o), aTb.target);
		Ravelin::Vector3d vA2a(Ravelin::Origin3d(eA->v2->o), aTb.target);
		LineSeg3 lineA(vA1a,vA2a);

		//Transforming B
		Polyhedron::VertexFaceIterator vfiBb(faceB,true);
 		std::vector<Ravelin::Vector3d> vBa;
		while(vfiBb.has_next()){
    		boost::shared_ptr<Polyhedron::Vertex> v=*vfiBb;
    		vfiBb.advance();
    		Ravelin::Vector3d p(v->o, aTb.source);
    		Ravelin::Vector3d pa=aTb.transform_point(p);
    		vBa.push_back(pa);
  		}

		//Triangulate
  		std::vector<Triangle> triB;
  		CompGeom::triangulate_convex_polygon(vBa.begin(),vBa.end(),triB.begin());
  		
  		//finding minimum distance between line and all triangles
  		Point3d p1,p2;
  		double min = DBL_MAX;
  		for(std::vector<Triangle>::iterator t = triB.begin(); t!=triB.end();++t){
  			double tmp=Triangle::calc_sq_dist(*t,lineA,p1,p2);
  			if(tmp<min){
  				min=tmp;
  			}
  		}
    	return min;
	}else if(fA == eFace && fB == eEdge){

		//Casting pointers
		boost::shared_ptr<const Polyhedron::Face> _faceA = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<Polyhedron::Face> faceA = boost::const_pointer_cast<Polyhedron::Face>(_faceA);
		boost::shared_ptr<const Polyhedron::Edge> eB = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);

		//creating Iterator for A
		Polyhedron::VertexFaceIterator vfiAa(faceA,true);
 		std::vector<Ravelin::Vector3d> vAa;
		while(vfiAa.has_next()){
    		boost::shared_ptr<Polyhedron::Vertex> v=*vfiAa;
    		vfiAa.advance();
    		Ravelin::Vector3d pa(v->o, aTb.target);
    		vAa.push_back(pa);
  		}

  		//Triangulating
		std::vector<Triangle> triA;
  		CompGeom::triangulate_convex_polygon(vAa.begin(),vAa.end(),triA.begin());
  		
		//Transforming and creating lineSegB
		Ravelin::Vector3d vB1b(Ravelin::Origin3d(eB->v1->o), aTb.source);
		Ravelin::Vector3d vB2b(Ravelin::Origin3d(eB->v2->o), aTb.source);
		
		//transforming fB into frame A
		Ravelin::Vector3d vB1a = aTb.transform_point(vB1b);
		Ravelin::Vector3d vB2a = aTb.transform_point(vB2b);

  		//creating segment
 		LineSeg3 lineB(vB1a,vB2a);

  		//finding minimum distance between line and all triangles
  		Point3d p1,p2;
  		double min = DBL_MAX;
		for(std::vector<Triangle>::iterator t = triA.begin(); t!=triA.end();++t){
  			double tmp=Triangle::calc_sq_dist(*t,lineB,p1,p2);
  			if(tmp<min){
  				min=tmp;
  			}
  		}
    	return min;
	}else if(fA == eFace && fB == eFace){
		//Casting pointers
		boost::shared_ptr<const Polyhedron::Face> _faceA = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
    boost::shared_ptr<Polyhedron::Face> faceA = boost::const_pointer_cast<Polyhedron::Face>(_faceA);
		boost::shared_ptr<const Polyhedron::Face> _faceB = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
    boost::shared_ptr<Polyhedron::Face> faceB = boost::const_pointer_cast<Polyhedron::Face>(_faceB);

		//creating Iterator for A
		Polyhedron::VertexFaceIterator vfiAa(faceA,true);
 		std::vector<Ravelin::Vector3d> vAa;
		while(vfiAa.has_next()){
    		boost::shared_ptr<Polyhedron::Vertex> v=*vfiAa;
    		vfiAa.advance();
    		Ravelin::Vector3d pa(v->o, aTb.target);
    		vAa.push_back(pa);
  		}

  		//Triangulating
		std::vector<Triangle> triA;
  		CompGeom::triangulate_convex_polygon(vAa.begin(),vAa.end(),triA.begin());
	
		//creating iterator for B
  		Polyhedron::VertexFaceIterator vfiBb(faceB,true);
 		std::vector<Ravelin::Vector3d> vBa;
		while(vfiBb.has_next()){
    		boost::shared_ptr<Polyhedron::Vertex> v=*vfiBb;
    		vfiBb.advance();
    		Ravelin::Vector3d p(v->o, aTb.source);
    		Ravelin::Vector3d pa=aTb.transform_point(p);
    		vBa.push_back(pa);
  		}

		//Triangulate
  		std::vector<Triangle> triB;
  		CompGeom::triangulate_convex_polygon(vBa.begin(),vBa.end(),triB.begin());
  		
  		//finding minimum distance between all pairs of triangles
  		Point3d p1,p2;
  		double min = DBL_MAX;
		for(std::vector<Triangle>::iterator tA = triA.begin(); tA!=triA.end();++tA){
  			for(std::vector<Triangle>::iterator tB = triB.begin(); tB!=triB.end();++tB){
  				double tmp=Triangle::calc_sq_dist(*tA,*tB,p1,p2);
  				if(tmp<min){
  					min=tmp;
  				}
  			}
  		}
  		return min;
	}
}

/// Does the case of update vertex/vertex
Polyhedron::UpdateRule Polyhedron::update_vertex_vertex(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  // search for Voronoi plane from three coincident edges to vertex A,
  // which vertex B violates

    // vertex B violates plane from this edge; update vA to eA

    return eContinue;

  // search for Voronoi plane from three coincident edges to vertex B,
  // which vertex A violates

    // vertex A violates plane from this edge; update vB to eB

    return eContinue;

  
  // still here = no violations
  return eDone;
}

/// Does the case of update vertex/edge
Polyhedron::UpdateRule Polyhedron::update_vertex_edge(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
  // search for Voronoi plane from those coincident to eB that vA violates

    // vA violates plane; update edge to coincident plane

    return eContinue;

  // clip eB against the Voronoi region of vA

  // ... to be continued... 
}

/// Does the case of update vertex/face
Polyhedron::UpdateRule Polyhedron::update_vertex_face(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
}

/// Does the case of update edge/edge
Polyhedron::UpdateRule Polyhedron::update_edge_edge(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
}

/// Does the case of update edge/face
Polyhedron::UpdateRule Polyhedron::update_edge_face(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
}

/// Does the case of update face/face
Polyhedron::UpdateRule Polyhedron::update_face_face(FeatureType& fA, FeatureType& fB, Transform3d& aTb, shared_ptr<const Polyhedron::Feature>& closestA, shared_ptr<const Polyhedron::Feature>& closestB)
{
}

