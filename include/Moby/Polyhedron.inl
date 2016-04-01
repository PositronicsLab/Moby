/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

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

/// Computes the convex hull for a polyhedron
template <class ForwardIterator>
Polyhedron Polyhedron::calc_convex_hull(ForwardIterator begin, ForwardIterator end)
{
  // ******************* compute convex hull begins *******************
  const unsigned X = 0, Y = 1, Z = 2;
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";  // TODO: remove Fx option?
  FILE* outfile, * errfile;

  // setup vertices for the convex hull
  std::vector<boost::shared_ptr<Polyhedron::Vertex> > verts;
  while (begin != end)
  {
    verts.push_back(boost::shared_ptr<Polyhedron::Vertex>(new Polyhedron::Vertex));
    verts.back()->o = Ravelin::Origin3d(*begin);
    begin++;
  }

  // get number of vertices
  const unsigned NVERTS = verts.size();
  
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
    vertex_map[points_begin+j] = verts[i];
    qhull_points[j++] = verts[i]->o[X];
    qhull_points[j++] = verts[i]->o[Y];
    qhull_points[j++] = verts[i]->o[Z];
  }

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code != 0)
  {
    // points are not collinear.. unsure of the error...
    FILE_LOG(LOG_COMPGEOM) << "Polyhedron::calc_convex_hull() - unable to execute qhull on points:" << std::endl;
    for (unsigned i=0; i< NVERTS; i++)
      FILE_LOG(LOG_COMPGEOM) << "  " << verts[i] << std::endl;

    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw Ravelin::NumericalException(); 
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

  // TODO: remove this when we can handle non-simplicial facets
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

  // make convexity on the polyhedron
  poly._convexity_computed = true;
  poly._convexity = -1.0;

  // compute the bounding box
  poly.calc_bounding_box();

  return poly; 
}

// undefine CREATE_LOOKUP
#undef CREATE_LOOKUP


