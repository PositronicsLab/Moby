/// Determines contact data between two geometries that are touching or interpenetrating
template <class OutputIterator>
OutputIterator CCD::find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{
  // look for special cases
  PrimitivePtr pA = cgA->get_geometry();
  PrimitivePtr pB = cgB->get_geometry();
  if (boost::dynamic_pointer_cast<SpherePrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_sphere_sphere(cgA, cgB, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return find_contacts_sphere_plane(cgA, cgB, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<BoxPrimitive>(pB))
      return find_contacts_box_sphere(cgB, cgA, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pB))
      return find_contacts_sphere_heightmap(cgA, cgB, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<BoxPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return find_contacts_plane_generic(cgB, cgA, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pB))
      return find_contacts_polyhedron_polyhedron(cgA, cgB, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_box_sphere(cgA, cgB, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_sphere_heightmap(cgB, cgA, output_begin, TOL);
    else if (pB->is_convex())
      return find_contacts_convex_heightmap(cgB, cgA, output_begin, TOL);
    else
      return find_contacts_heightmap_generic(cgA, cgB, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<PlanePrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_sphere_plane(cgB, cgA, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<CylinderPrimitive>(pB))
      return find_contacts_cylinder_plane(cgB, cgA, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<TorusPrimitive>(pB))
      return find_contacts_torus_plane(cgB, cgA, output_begin, TOL);
    else
      return find_contacts_plane_generic(cgA, cgB, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<CylinderPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return find_contacts_cylinder_plane(cgA, cgB, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<PolyhedralPrimitive>(pB))
      return find_contacts_polyhedron_polyhedron(cgA, cgB, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return find_contacts_plane_generic(cgB, cgA, output_begin, TOL);
  }
  else if (boost::dynamic_pointer_cast<TorusPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
    {
      return find_contacts_torus_plane(cgA, cgB, output_begin, TOL);
    }
  }
  else // no special case for A
  {
    if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pB))
    {
      if (pA->is_convex())
        return find_contacts_convex_heightmap(cgA, cgB, output_begin, TOL);
      else
        return find_contacts_heightmap_generic(cgB, cgA, output_begin, TOL);
    }
    else if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
    {
      return find_contacts_plane_generic(cgB, cgA, output_begin, TOL);
    }
  }

  // still here? just use the generic contact finder
  return find_contacts_generic(cgA, cgB, output_begin, TOL);
}

inline std::vector<Ravelin::Vector3d> CCD::create_edge_vector(const std::vector<boost::shared_ptr<Polyhedron::Edge> > &edges, Ravelin::Transform3d &wTe)
{
  std::vector<Ravelin::Vector3d > edge_vectors;
  BOOST_FOREACH(boost::shared_ptr<Polyhedron::Edge> edge, edges)
  {
    Ravelin::Vector3d v1_e(edge->v1->o,wTe.source);
    Ravelin::Vector3d v2_e(edge->v2->o,wTe.source);

    Ravelin::Vector3d v1_w = wTe.transform_point(v1_e);
    Ravelin::Vector3d v2_w = wTe.transform_point(v2_e);
    
    Ravelin::Vector3d edge_v = (v2_w - v1_w);
    edge_v.normalize();
    edge_vectors.push_back(edge_v);
  }

  return edge_vectors;
}

inline void CCD::project(std::vector<Ravelin::Vector3d> vectors, Ravelin::Vector3d axis, double &min_dot, double &max_dot, int &min_index, int &max_index)
{
  for(int i = 0 ; i < vectors.size(); ++i)
  {
    double value = axis.dot(vectors[i]);
    if (!min_dot || value < min_dot)
    {
      min_dot = value;
      min_index = i;
    }

    if(!max_dot || value > max_dot)
    {
      max_dot = value;
      max_index = i;
    }

  }
}

inline void CCD::create_convex_hull_list(boost::shared_ptr<Polyhedron::Vertex> start_vert, const Ravelin::Vector3d& axis, const Ravelin::Transform3d& wTv, std::vector<Ravelin::Vector3d>& ch_vectors)
{
  std::vector <boost::shared_ptr<Polyhedron::Vertex> > visited_vertices;
  std::queue <boost::shared_ptr<Polyhedron::Vertex> > search_list;

  //adding the first projection into the map and list
  search_list.push(start_vert);
  Ravelin::Vector3d v(start_vert->o,wTv.source);
  Ravelin::Vector3d vw = wTv.transform_point(v);
  visited_vertices.push_back(start_vert);
  ch_vectors.push_back(vw);

  //create testing plane
  Plane test_plane(axis,vw);

  while(!search_list.empty())
  {
    boost::shared_ptr<Polyhedron::Vertex> cur_vertex = search_list.front();
    search_list.pop();
    std::list<boost::weak_ptr<Polyhedron::Edge> > edges = cur_vertex->e;
    for (std::list<boost::weak_ptr<Polyhedron::Edge> >::iterator ei = edges.begin(); ei!=edges.end(); ei++)
    {
      boost::shared_ptr<Polyhedron::Edge> neighbor_e(*ei);
      boost::shared_ptr<Polyhedron::Vertex> next_vertex = neighbor_e->v1;
      if( next_vertex == cur_vertex)
      {
        next_vertex = neighbor_e->v2;
      }

      std::vector <boost::shared_ptr<Polyhedron::Vertex> >::iterator vvi = std::find(visited_vertices.begin(), visited_vertices.end(), next_vertex);
      if(vvi != visited_vertices.end())
      {
          // we have visit this 
        continue;
      }
      else
      {
        visited_vertices.push_back(cur_vertex);
        Ravelin::Vector3d v(cur_vertex->o,wTv.source);
        Ravelin::Vector3d vw = wTv.transform_point(v);
        if (fabs(test_plane.calc_signed_distance(vw)) < NEAR_ZERO)
        {
          ch_vectors.push_back(vw);
          search_list.push(next_vertex);
        }
      }
    }
  }
}

/// Finds contacts between two polyhedra
template <class OutputIterator>
OutputIterator CCD::find_contacts_polyhedron_polyhedron(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{
  const double INF = std::numeric_limits<double>::max();
  const unsigned TRI_VERTS = 3;
  std::vector<std::pair<Ravelin::Vector3d, double> > hs;
  enum FeatureType { eNone, eVertex, eEdge, eFace };
  FeatureType featA = eNone, featB = eNone;
  std::vector<Point3d> verts;
  boost::shared_ptr<Polyhedron::Vertex> vA, vB;
  boost::shared_ptr<Polyhedron::Edge> eA, eB;
  boost::shared_ptr<Polyhedron::Face> fA, fB;
  TessellatedPolyhedronPtr tpoly;
  boost::shared_ptr<Ravelin::Pose2d> GLOBAL_2D;

  // get the two primitives
  boost::shared_ptr<const PolyhedralPrimitive> pA = boost::dynamic_pointer_cast<const PolyhedralPrimitive>(cgA->get_geometry());
  boost::shared_ptr<const PolyhedralPrimitive> pB = boost::dynamic_pointer_cast<const PolyhedralPrimitive>(cgB->get_geometry());

  // get the two polyhedra
  const Polyhedron& polyA = pA->get_polyhedron();
  const Polyhedron& polyB = pB->get_polyhedron();

  // get the two poses
  boost::shared_ptr<const Ravelin::Pose3d> poseA = pA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> poseB = pB->get_pose(cgB);

  // get transforms to global frame
  Ravelin::Transform3d wTa = Ravelin::Pose3d::calc_relative_pose(poseA, GLOBAL);
  Ravelin::Transform3d wTb = Ravelin::Pose3d::calc_relative_pose(poseB, GLOBAL);

  // call v-clip
  boost::shared_ptr<const Polyhedron::Feature> closestA;
  boost::shared_ptr<const Polyhedron::Feature> closestB;
  double dist = Polyhedron::vclip(pA, pB, poseA, poseB, closestA, closestB);
  FILE_LOG(LOG_COLDET) << "v-clip reports distance of " << dist << std::endl;

  // see whether to generate contacts
  if (dist > TOL)
    return output_begin;

  // The normal, which will be the vector of translation, and the contact plane.
  Ravelin::Vector3d normal;
  Plane contact_plane;

  // Determine the vector to translate the geometries. This vector defines
  // the contact plane. Once the geometries are translated, all vertices on the
  // contact plane are treated as being on the contact manifold.
  if (dist > 0.0)
  {
    // Setup the normal in the global frame
    Point3d closestAw, closestBw;
    find_closest_points(closestA, closestB, wTa, wTb, closestAw, closestBw);
    normal =  closestAw - closestBw;
    normal.normalize();
    contact_plane.set_normal(normal);
    contact_plane.offset = 0.5*normal.dot(closestAw - closestBw);
  }
  else
  {
    // Compute first set of testing vectors from face normals
    std::vector<Ravelin::Vector3d> test_vectors;
    const std::vector<boost::shared_ptr<Polyhedron::Face> >& fA = polyA.get_faces();
    const std::vector<boost::shared_ptr<Polyhedron::Face> >& fB = polyB.get_faces();
    for (unsigned i=0; i< fA.size(); i++)
      test_vectors.push_back(wTa.transform_vector(Ravelin::Vector3d(fA[i]->get_plane().get_normal().data(), poseA)));
    for (unsigned i=0; i< fB.size(); i++)
      test_vectors.push_back(wTb.transform_vector(Ravelin::Vector3d(fB[i]->get_plane().get_normal().data(), poseB)));

    // create testing axes (cross-products of edges from A and B)
    const std::vector<boost::shared_ptr<Polyhedron::Edge> >& edgesA = polyA.get_edges();
    const std::vector<boost::shared_ptr<Polyhedron::Edge> >& edgesB = polyB.get_edges();
    std::vector<Ravelin::Vector3d> evA = create_edge_vector(edgesA, wTa);
    std::vector<Ravelin::Vector3d> evB = create_edge_vector(edgesB, wTb);
    for (std::vector<Ravelin::Vector3d>::iterator evAi = evA.begin(); evAi != evA.end(); ++evAi) {
      for (std::vector<Ravelin::Vector3d>::iterator evBi = evB.begin(); evBi != evB.end(); ++evBi) {
        Ravelin::Vector3d xv = Ravelin::Vector3d::cross(*evAi, *evBi);
        double nrm = xv.norm();
        if (nrm > NEAR_ZERO) {
          xv /= nrm;
          test_vectors.push_back(xv);
        }
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
      int min_index_a, min_index_b, max_index_a, max_index_b;

      // project vertices onto the candidate axis
      // NOTE: this could be done in lg n time of the number of vertices rather
      // than 
      project(vector_a, *test_i, min_a, max_a, min_index_a, max_index_a);
      project(vector_b, *test_i, min_b, max_b, min_index_b, max_index_b);

      // Compute the amount of overlap.
      double o1 = max_a - min_b;
      double o2 = max_b - min_a;

      if (o1 > NEAR_ZERO && o2 > NEAR_ZERO) {
        // there is an overlap
        double overlap = std::min(o1, o2);
        boost::shared_ptr <Polyhedron::Vertex> v1, v2;

        if (min_overlap - overlap > NEAR_ZERO) {
          min_overlap = overlap;
          min_axis = *test_i;
          if (fabs(overlap - o1) > NEAR_ZERO) {
            direction = -1;
            a_vertex = vAa[max_index_a];
            b_vertex = vBb[min_index_b];
          } else {
            direction = 1;
            b_vertex = vBb[max_index_b];
            a_vertex = vAa[min_index_a];
          }
        }
      }
    }

    // ensure that the distance is negated
    assert(min_overlap >= 0.0);
    dist = -min_overlap;
    Point3d closestAw = wTa.transform_point(Point3d(a_vertex->o, poseA));
    Point3d closestBw = wTb.transform_point(Point3d(b_vertex->o, poseB));
    normal =  closestAw - closestBw;
    normal.normalize();
    contact_plane.set_normal(normal);
    contact_plane.offset = 0.5*normal.dot(closestAw - closestBw);
  }

  // TODO: need to look for the case where there are two minimally overlapping axes
  // TODO: do we need to call v-clip to find closest points (for when bodies are intersecting), but after they have been pushed apart to a kissing configuration?
  // TODO: setup the contact plane for intersecting bodies
  // TODO: figure out what is the case for distance zero
  // TODO: make the normal offset method work properly for intersecting bodies
  const double HALF_DIST = dist * 0.5;

  // Determine the vertices from A that are on the contact plane.
  std::vector<Point3d> vertsA_on_contact;
  const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vertsA = polyA.get_vertices();
  for (unsigned i=0; i< vertsA.size(); i++)
  {
    // project the point toward the plane
    Point3d point = wTa.transform_point(Point3d(vertsA[i]->o, poseA)) - normal*HALF_DIST;
    if (contact_plane.calc_signed_distance(point) < NEAR_ZERO)
      vertsA_on_contact.push_back(point);
  }

  // Determine the vertices from B that are on the contact plane
  std::vector<Point3d> vertsB_on_contact;
  const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vertsB = polyB.get_vertices();
  for (unsigned i=0; i< vertsB.size(); i++)
  {
    // project the point toward the plane
    Point3d point = wTb.transform_point(Point3d(vertsB[i]->o, poseB)) + normal*HALF_DIST;
    if (contact_plane.calc_signed_distance(point) < NEAR_ZERO)
      vertsB_on_contact.push_back(point);
  }

  // Setup the points of intersection.
  std::vector<Ravelin::Vector3d> convhull_a, convhull_b;

  // Try to calculate convex hull of A
  if (vertsA_on_contact.size() >= 3)  // minimum necessary for 2D hull
  {
    // Hull might still fail
    try {
      // computing the convex hull of the vertices from A
      CompGeom::calc_convex_hull(vertsA_on_contact.begin(),
                                 vertsA_on_contact.end(),
                                 normal,
                                 std::back_inserter(convhull_a));
    }
    catch (Ravelin::NumericalException e) {}

    // check for ccw
    if (!CompGeom::ccw(convhull_a.begin(), convhull_a.end(),normal))
      std::reverse(convhull_a.begin(), convhull_a.end());
  }

  // Try to calculate convex hull of B
  if (vertsB_on_contact.size() >= 3)  // minimum necessary for 2D hull
  {
    // Hull might still fail
    try {
      // computing the convex hull of the vertices from b
      CompGeom::calc_convex_hull(vertsB_on_contact.begin(),
                                 vertsB_on_contact.end(),
                                 normal,
                                 std::back_inserter(convhull_b));
    }
    catch (Ravelin::NumericalException e) {}

    // Make the ordering of the polygon ccw with respect to the normal.
    if (!CompGeom::ccw(convhull_b.begin(), convhull_b.end(),normal))
      std::reverse(convhull_b.begin(), convhull_b.end());
  }

  // Now determine the type of contact shape.
  if (!convhull_a.empty()) {
    // Look for case that points from B correspond to a polygon too
    if (!convhull_b.empty()) {

      // Intersect the polygons
      std::vector <Point3d> isect;
      CompGeom::intersect_convex_polygons(convhull_a.begin(),
                                          convhull_a.end(),
                                          convhull_b.begin(),
                                          convhull_b.end(),
                                          normal,
                                          std::back_inserter(isect));

      // Output the points.
      for (unsigned i = 0; i < isect.size(); i++)
        *output_begin++ = create_contact(cgA, cgB, isect[i], normal, dist);
      return output_begin;
    } else {
      // We won't have a polygon from B. vertsB either represents a point or
      // a line segment. Look for both cases by first finding line endpoints.
      LineSeg3 endpoints;
      CompGeom::determine_seg_endpoints(vertsB_on_contact.begin(),
                                        vertsB_on_contact.end(), endpoints);
      if ((endpoints.first - endpoints.second).norm() < NEAR_ZERO) {
        // It's effectively just a single point. We'll use it. Note that we
        // could try something like seeing whether the point is inside the
        // polygon. We know the bodies are in contact and we know the normal.
        // We just need a point.
        *output_begin++ = create_contact(cgA, cgB, endpoints.first, normal,
                                         dist);
        return output_begin;
      }
      else
      {
        // Two separate points: intersect the line segment with the polygon.
        double te, tl;
        if (!CompGeom::intersect_seg_convex_polygon(convhull_a.begin(),
                                                   convhull_a.end(),
                                                   normal,
                                                   endpoints, te, tl)) {
          // In the case of no intersection, we _should_ log an error message so
          // that this behavior can be examined in more detail. Instead, using
          // the reasoning above, we'll choose to believe that having a somewhat
          // wrong point (as would be the case if there were no intersection) is
          // better than having no point at all.
          *output_begin++ = create_contact(cgA, cgB, endpoints.first, normal,
                                           dist);
          return output_begin;
        }

        Point3d p1 = (1-tl)*endpoints.first + tl*endpoints.second;
        Point3d p2 = (1-te)*endpoints.first + te*endpoints.second;
        *output_begin++ = create_contact(cgA, cgB, p1, normal, dist);
        *output_begin++ = create_contact(cgA, cgB, p2, normal, dist);
        return output_begin;
      }
    }
  }
  else  // vertsA does not represent a polygon
  {
    // vertsA represents either a point or a line segment. Look for both
    // cases by first finding line endpoints.
    LineSeg3 endpointsA;
    CompGeom::determine_seg_endpoints(vertsA_on_contact.begin(),
                                      vertsA_on_contact.end(), endpointsA);

    // Look for the case that points from B correspond to a polygon.
    if (!convhull_b.empty()) {

      // Determine the topology of points from A.
      if ((endpointsA.first - endpointsA.second).norm() < NEAR_ZERO) {
        // We've found a point. See discussion on identical case above.
        *output_begin++ = create_contact(cgA, cgB, endpointsA.first,
                                         normal, dist);
        return output_begin;
      } else {
        // A represents a line seg; intersect the line segment with the polygon.
        double te, tl;
        if (!CompGeom::intersect_seg_convex_polygon(convhull_a.begin(),
                                               convhull_a.end(),
                                               normal,
                                               endpointsA, te, tl)) {
          // In the case of no intersection, we _should_ log an error message so
          // that this behavior can be examined in more detail. Instead, using
          // the reasoning above, we'll choose to believe that having a somewhat
          // wrong point (as would be the case if there were no intersection) is
          // better than having no point at all.
            *output_begin++ = create_contact(cgA, cgB, endpointsA.first, normal,
                                             dist);
            return output_begin;
        }

        Point3d p1 = (1-tl)*endpointsA.first + tl*endpointsA.second;
        Point3d p2 = (1-te)*endpointsA.first + te*endpointsA.second;
        *output_begin++ = create_contact(cgA, cgB, p1, normal, dist);
        *output_begin++ = create_contact(cgA, cgB, p2, normal, dist);
        return output_begin;
      }

      // Should not still be here- we just handled case that B is a polygon.
      assert(false);
    } else // vertics from B do not correspond to a polygon either.
    {
      // If we're here, then we could have a point from each shape, a
      // point and a line segment from the shapes, or a line segment
      // from each shape. The last is the one that requires the most
      // conversation because we then need to intersect the line segments.
      // For the same reason as was discussed in the point-on-polygon case
      // above, we'll simply return the single point from a shape otherwise.

      // Determine topology of points from B.
      LineSeg3 endpointsB;
      CompGeom::determine_seg_endpoints(vertsB_on_contact.begin(),
                                        vertsB_on_contact.end(), endpointsB);

      // Look for the easy cases
      if ((endpointsA.first - endpointsA.second).norm() < NEAR_ZERO) {
        *output_begin++ = create_contact(cgA, cgB, endpointsA.first, normal,
                                         dist);
        return output_begin;
      }
      if ((endpointsB.first - endpointsB.second).norm() < NEAR_ZERO) {
        *output_begin++ = create_contact(cgA, cgB, endpointsB.first, normal,
                                         dist);
        return output_begin;
      }

      // If we're still here- it's the line intersection case. Urgh.
      // Rather than intersect the line segments properly, which is a PITA,
      // let's just find a closest point between them.
      Point3d cpA, cpB;
      CompGeom::calc_closest_points(endpointsA, endpointsB, cpA, cpB);
      *output_begin++ =
          create_contact(cgA, cgB, (cpA + cpB) * 0.5, normal, dist);
      return output_begin;
    } // end else for points from B not belonging to a polygon
  } // end else for dist > 0.0
    /*

else if (vertsA_on_contact.size() == 1) {

  // When there is only one vertex of polyhedron A
  // The intersection will be that point
  isect.push_back(vertsA_on_contact.front());

} else if (vertsB_on_contact.size() == 1) {

  // When there is only one vertex of polyhedron A
  // The intersection will be that point
  isect.push_back(vertsB_on_contact.front());

} else if (vertsA_on_contact.size() == 2 && vertsB_on_contact.size() == 2 ) {

  // When both of them are segments
  // We then intersects the segments to find intersection of them.

  LineSeg3 seg_a(Point3d(vertsA_on_contact[0].data(), GLOBAL), Point3d(vertsA_on_contact[1].data(), GLOBAL));
  LineSeg3 seg_b(Point3d(vertsA_on_contact[0].data(), GLOBAL), Point3d(vertsA_on_contact[1].data(), GLOBAL));

  Point3d isect1,isect2;

  CompGeom::SegSegIntersectType result = CompGeom::intersect_segs(seg_a, seg_b, isect1, isect2);

  // Maybe a switch on the result should be used(?)

  if(result != CompGeom::eSegSegNoIntersect){
    isect.push_back(isect1);
    if(result == CompGeom::eSegSegEdge){
      isect.push_back(isect2);
    }
  }

} else if (vertsA_on_contact.size() == 2){

  // segment A vs polygon B
  // we then tries to intersect the segment a with polygon B


  // projecting all points to 2d
  Ravelin::Matrix3d _2DT3D = CompGeom::calc_3D_to_2D_matrix(normal);
  double offset_3d = CompGeom::determine_3D_to_2D_offset(vertsA_on_contact[0], _2DT3D);

  std::vector <Ravelin::Origin2d> verts_2D_a, verts_2D_b;

  for (std::vector<Ravelin::Vector3d>::iterator vAoci = vertsA_on_contact.begin();
        vAoci != vertsA_on_contact.end(); ++vAoci) {
    verts_2D_a.push_back(CompGeom::to_2D(*vAoci, _2DT3D));
  }

  for (std::vector<Ravelin::Vector3d>::iterator vBoci = vertsB_on_contact.begin();
      vBoci != vertsB_on_contact.end(); ++vBoci) {

    verts_2D_b.push_back(CompGeom::to_2D(*vBoci, _2DT3D));

  }

  // create segment A
  LineSeg2 seg_a(Point2d(verts_2D_a[0], GLOBAL_2D),
                  Point2d(verts_2D_a[1], GLOBAL_2D));

  // create the polygon by computing convex hull
  std::vector <Point2d> convhull_b_2D;
  std::vector <LineSeg2> isect2;
  CompGeom::calc_convex_hull(verts_2D_b.begin(), verts_2D_b.end(),
                              convhull_b_2D.begin());

  // make sure the list is counter-clock-wise
  if (!CompGeom::ccw(convhull_b_2D.begin(), convhull_b_2D.end())) {
    std::reverse(convhull_b_2D.begin(), convhull_b_2D.end());
  }

  // intersecting the line and the polygon
  CompGeom::intersect_seg_polygon(convhull_b_2D.begin(), convhull_b_2D.end(),
                                  seg_a, isect2.begin());

  // projecting the points back to 3d and add them to the list
  for (std::vector<LineSeg2>::iterator isecti = isect2.begin();
        isecti != isect2.end(); ++isecti) {

        LineSeg2 l = *isecti;

        Ravelin::Vector2d v = l.first;
        Ravelin::Origin2d o(v.x(), v.y());

        isect.push_back(Ravelin::Vector3d(
                        CompGeom::to_3D(o,_2DT3D.inverse(), offset_3d),
                        GLOBAL));

        v = l.second;
        o = Ravelin::Origin2d(v.x(), v.y());

        isect.push_back(Ravelin::Vector3d(
                        CompGeom::to_3D(o, _2DT3D.inverse(), offset_3d),
                        GLOBAL));
  }

} else if (vertsB_on_contact.size() == 2) {

  // segment B vs polygon A
  // we then tries to intersect the segment B with polygon A


  // projecting all points to 2d
  Ravelin::Matrix3d _2DT3D = CompGeom::calc_3D_to_2D_matrix(normal);
  double offset_3d = CompGeom::determine_3D_to_2D_offset(vertsA_on_contact[0], _2DT3D);

  std::vector <Ravelin::Origin2d> verts_2D_a, verts_2D_b;

  for (std::vector<Ravelin::Vector3d>::iterator vAoci = vertsA_on_contact.begin();
        vAoci != vertsA_on_contact.end(); ++vAoci) {
    verts_2D_a.push_back(CompGeom::to_2D(*vAoci, _2DT3D));
  }

  for (std::vector<Ravelin::Vector3d>::iterator vBoci = vertsB_on_contact.begin();
      vBoci != vertsB_on_contact.end(); ++vBoci) {

    verts_2D_b.push_back(CompGeom::to_2D(*vBoci, _2DT3D));

  }

  // create segment b
  LineSeg2 seg_b(Point2d(verts_2D_b[0], GLOBAL_2D),
                  Point2d(verts_2D_b[1], GLOBAL_2D));

  // create the polygon by computing convex hull
  std::vector <Point2d> convhull_a_2D;
  std::vector <LineSeg2> isect2;
  CompGeom::calc_convex_hull(verts_2D_a.begin(), verts_2D_a.end(),
                              convhull_a_2D.begin());

  // make sure the list is counter-clock-wise
  if (!CompGeom::ccw(convhull_a_2D.begin(), convhull_a_2D.end())) {
    std::reverse(convhull_a_2D.begin(), convhull_a_2D.end());
  }

  // intersecting the line and the polygon
  CompGeom::intersect_seg_polygon(convhull_a_2D.begin(), convhull_a_2D.end(),
                                  seg_b, isect2.begin());

  // projecting the points back to 3d and add them to the list
  for (std::vector<LineSeg2>::iterator isecti = isect2.begin();
        isecti != isect2.end(); ++isecti) {

        LineSeg2 l = *isecti;

        Ravelin::Vector2d v = l.first;
        Ravelin::Origin2d o(v.x(), v.y());

        isect.push_back(Ravelin::Vector3d(
                        CompGeom::to_3D(o,_2DT3D.inverse(), offset_3d),
                        GLOBAL));

        v = l.second;
        o = Ravelin::Origin2d(v.x(), v.y());

        isect.push_back(Ravelin::Vector3d(
                        CompGeom::to_3D(o, _2DT3D.inverse(), offset_3d),
                        GLOBAL));
  }
}





CompGeom::calc_convex_hull(verts.begin(), verts.end(), normal, std::back_inserter(isect));

  // Output the points.
for (unsigned i = 0; i < isect.size(); i++)
  *output_begin++ = create_contact(cgA, cgB, isect[i], normal, dist);
*/
  /*
  // case #1: attempt to use volume of intersection
  if (dist <= 0.0) {
    // initialize
    std::vector <boost::shared_ptr<Polyhedron::Edge> >
        edgesA = polyA.get_edges();
    std::vector <boost::shared_ptr<Polyhedron::Edge> >
        edgesB = polyB.get_edges();
    std::vector <boost::shared_ptr<Polyhedron::Vertex> >
        vAa = polyA.get_vertices();
    std::vector <boost::shared_ptr<Polyhedron::Vertex> >
        vBb = polyB.get_vertices();

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

    // create testing axes (edges from A, edges from B, and their cross product)
    std::vector <Ravelin::Vector3d> evA = create_edge_vector(edgesA, wTa);
    std::vector <Ravelin::Vector3d> evB = create_edge_vector(edgesB, wTb);
    std::vector <Ravelin::Vector3d> test_vectors;
    for (std::vector<Ravelin::Vector3d>::iterator evAi = evA.begin();
         evAi != evA.end(); ++evAi) {
      for (std::vector<Ravelin::Vector3d>::iterator evBi = evB.begin();
           evBi != evB.end(); ++evBi) {
        Ravelin::Vector3d xv = Ravelin::Vector3d::cross(*evAi, *evBi);
        if (xv.norm() > NEAR_ZERO) {
          xv.normalize();
          test_vectors.push_back(xv);
        }
      }
    }
    test_vectors.insert(test_vectors.end(), evA.begin(), evA.end());
    test_vectors.insert(test_vectors.end(), evB.begin(), evB.end());

    // ***********************************************************************
    // find the minimum overlap
    // ***********************************************************************
    double min_overlap = std::numeric_limits<double>::max();
    Ravelin::Vector3d min_axis;
    boost::shared_ptr <Polyhedron::Vertex> a_vertex, b_vertex;
    int direction = 1;
    for (std::vector<Ravelin::Vector3d>::iterator test_i = test_vectors.begin();
         test_i != test_vectors.end(); ++test_i) {
      double min_a, max_a, min_b, max_b;
      boost::shared_ptr <Polyhedron::Vertex> minV_a, minV_b, maxV_a, maxV_b;
      int min_index_a, min_index_b, max_index_a, max_index_b;

      // projecting shapes to axis
      project(vector_a, *test_i, min_a, max_a, min_index_a, max_index_a);
      project(vector_b, *test_i, min_b, max_b, min_index_b, max_index_b);

      double o1 = max_a - min_b;
      double o2 = max_b - min_a;

      if (o1 > NEAR_ZERO && o2 > NEAR_ZERO) {
        // there is an overlap
        double overlap = std::min(o1, o2);
        boost::shared_ptr <Polyhedron::Vertex> v1, v2;

        if (min_overlap - overlap > NEAR_ZERO) {
          min_overlap = overlap;
          min_axis = *test_i;
          if (fabs(overlap - o1) > NEAR_ZERO) {
            a_vertex = vAa[max_index_a];
            b_vertex = vBb[min_index_b];
            direction = 1;
          } else {
            b_vertex = vBb[max_index_b];
            a_vertex = vAa[min_index_a];
            direction = -1;
          }
        }
      }
    }

    // feature search for two polyhedron 
    min_axis = min_axis * direction;
    std::vector <Ravelin::Vector3d> ch_vectors_a, ch_vectors_b;
    create_convex_hull_list(a_vertex, min_axis, wTa, ch_vectors_a);

    // project b (to what?)
    Ravelin::Vector3d trans_v = (min_axis * min_overlap);
    Ravelin::Origin3d trans_o(trans_v.x(), trans_v.y(), trans_v.z());
    Ravelin::Transform3d sep_trans(trans_o);
    Ravelin::Transform3d w_1Tb = sep_trans * wTb;
    create_convex_hull_list(b_vertex, min_axis, w_1Tb, ch_vectors_b);

    std::vector <Ravelin::Vector3d> ch_a, ch_b;
    std::vector <Ravelin::Vector3d> isect;
    if (ch_vectors_a.size() == 1) {
      isect.push_back(ch_vectors_a.front());
    } else if (ch_vectors_a.size() == 2) {

      Ravelin::Matrix3d _2DT3D = CompGeom::calc_3D_to_2D_matrix(min_axis);
      double offset_3d =
          CompGeom::determine_3D_to_2D_offset(ch_vectors_a[0], _2DT3D);
      std::vector <Ravelin::Origin2d> ch_2D_a, ch_2D_b;
      for (std::vector<Ravelin::Vector3d>::iterator chvi = ch_vectors_a.begin();
           chvi != ch_vectors_a.end(); ++chvi) {
        ch_2D_a.push_back(CompGeom::to_2D(*chvi, _2DT3D));
      }
      for (std::vector<Ravelin::Vector3d>::iterator chvi = ch_vectors_b.begin();
           chvi != ch_vectors_b.end(); ++chvi) {
        ch_2D_b.push_back(CompGeom::to_2D(*chvi, _2DT3D));
      }

      LineSeg2
          seg_a(Point2d(ch_2D_a[0], GLOBAL_2D), Point2d(ch_2D_a[1], GLOBAL_2D));

      std::vector <Point2d> ch_b;
      std::vector <LineSeg2> isect2;
      CompGeom::calc_convex_hull(ch_2D_b.begin(), ch_2D_b.end(), ch_b.begin());

      if (!CompGeom::ccw(ch_b.begin(), ch_b.end())) {
        std::reverse(ch_b.begin(), ch_b.end());
      }
      CompGeom::intersect_seg_polygon(ch_b.begin(),
                                      ch_b.end(),
                                      seg_a,
                                      isect2.begin());

      for (std::vector<LineSeg2>::iterator chvi = isect2.begin();
           chvi != isect2.end(); ++chvi) {
        LineSeg2 l = *chvi;
        Ravelin::Vector2d v = l.first;
        Ravelin::Origin2d o(v.x(), v.y());
        isect.push_back(Ravelin::Vector3d(CompGeom::to_3D(o,
                                                          _2DT3D.inverse(),
                                                          offset_3d), GLOBAL));
        v = l.second;
        o = Ravelin::Origin2d(v.x(), v.y());
        isect.push_back(Ravelin::Vector3d(CompGeom::to_3D(o,
                                                          _2DT3D.inverse(),
                                                          offset_3d), GLOBAL));
      }

    } else {
      if (ch_vectors_b.size() == 1) {
        isect.push_back(ch_vectors_b[0]);
      } else if (ch_vectors_b.size() == 2) {

        Ravelin::Matrix3d _2DT3D = CompGeom::calc_3D_to_2D_matrix(min_axis);
        double offset_3d =
            CompGeom::determine_3D_to_2D_offset(ch_vectors_b[0], _2DT3D);
        std::vector <Ravelin::Origin2d> ch_2D_a, ch_2D_b;
        for (std::vector<Ravelin::Vector3d>::iterator
                 chvi = ch_vectors_a.begin(); chvi != ch_vectors_a.end();
             ++chvi) {
          ch_2D_a.push_back(CompGeom::to_2D(*chvi, _2DT3D));
        }
        for (std::vector<Ravelin::Vector3d>::iterator
                 chvi = ch_vectors_b.begin(); chvi != ch_vectors_b.end();
             ++chvi) {
          ch_2D_b.push_back(CompGeom::to_2D(*chvi, _2DT3D));
        }

        LineSeg2 seg_b
            (Point2d(ch_2D_b[0], GLOBAL_2D), Point2d(ch_2D_b[1], GLOBAL_2D));

        std::vector <Point2d> ch_a;
        std::vector <LineSeg2> isect2;
        CompGeom::calc_convex_hull(ch_2D_a.begin(),
                                   ch_2D_a.end(),
                                   ch_a.begin());

        if (!CompGeom::ccw(ch_a.begin(), ch_a.end())) {
          std::reverse(ch_a.begin(), ch_a.end());
        }
        CompGeom::intersect_seg_polygon(ch_a.begin(),
                                        ch_a.end(),
                                        seg_b,
                                        isect2.begin());

        for (std::vector<LineSeg2>::iterator chvi = isect2.begin();
             chvi != isect2.end(); ++chvi) {
          LineSeg2 l = *chvi;
          Ravelin::Vector2d v = l.first;
          Ravelin::Origin2d o(v.x(), v.y());
          isect.push_back(Ravelin::Vector3d(CompGeom::to_3D(o,
                                                            _2DT3D.inverse(),
                                                            offset_3d),
                                            GLOBAL));
          v = l.second;
          o = Ravelin::Origin2d(v.x(), v.y());
          isect.push_back(Ravelin::Vector3d(CompGeom::to_3D(o,
                                                            _2DT3D.inverse(),
                                                            offset_3d),
                                            GLOBAL));
        }

      } else {
        CompGeom::calc_convex_hull(ch_vectors_a.begin(),
                                   ch_vectors_a.end(),
                                   min_axis,
                                   std::back_inserter(ch_a));
        CompGeom::calc_convex_hull(ch_vectors_b.begin(),
                                   ch_vectors_b.end(),
                                   min_axis,
                                   std::back_inserter(ch_b));
        CompGeom::intersect_convex_polygons(ch_vectors_a.begin(),
                                            ch_vectors_a.end(),
                                            ch_vectors_b.begin(),
                                            ch_vectors_b.end(),
                                            min_axis,
                                            isect.begin());
      }
    }

    for (unsigned i = 0; i < isect.size(); i++)
      *output_begin++ =
          create_contact(cgA, cgB, isect[i], min_axis, -min_overlap);
  }
*/
  return output_begin;
}

/// Finds contacts between two polyhedra
template <class OutputIterator>
OutputIterator CCD::find_contacts_polyhedron_polyhedron_soft(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{
  const double INF = std::numeric_limits<double>::max();
  std::vector<std::pair<Ravelin::Vector3d, double> > hs;
  enum FeatureType { eNone, eVertex, eEdge, eFace };
  FeatureType featA = eNone, featB = eNone;
  boost::shared_ptr<Polyhedron::Vertex> vA, vB;
  boost::shared_ptr<Polyhedron::Edge> eA, eB;
  boost::shared_ptr<Polyhedron::Face> fA, fB;
  TessellatedPolyhedronPtr tpoly;
  boost::shared_ptr<Ravelin::Pose2d> GLOBAL_2D;

  // get the two primitives
  boost::shared_ptr<const PolyhedralPrimitive> pA = boost::dynamic_pointer_cast<const PolyhedralPrimitive>(cgA->get_geometry());
  boost::shared_ptr<const PolyhedralPrimitive> pB = boost::dynamic_pointer_cast<const PolyhedralPrimitive>(cgB->get_geometry());

  // get the two polyhedra
  const Polyhedron& polyA = pA->get_polyhedron();
  const Polyhedron& polyB = pB->get_polyhedron();

  // get the two poses
  boost::shared_ptr<const Ravelin::Pose3d> poseA = pA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> poseB = pB->get_pose(cgB);

  // get transforms to global frame
  Ravelin::Transform3d wTa = Ravelin::Pose3d::calc_relative_pose(poseA, GLOBAL);
  Ravelin::Transform3d wTb = Ravelin::Pose3d::calc_relative_pose(poseB, GLOBAL);

  // call v-clip
  boost::shared_ptr<const Polyhedron::Feature> closestA;
  boost::shared_ptr<const Polyhedron::Feature> closestB;
  double dist = Polyhedron::vclip(pA, pB, poseA, poseB, closestA, closestB);
  FILE_LOG(LOG_COLDET) << "v-clip reports distance of " << dist << std::endl;

  // see whether to generate contacts
  if (dist > TOL)
    return output_begin;

  if (dist < 0.0)
  {
    // NOTE: this is Evan's LP based implementation
    // get the halfspaces
    PolyhedralPrimitive::get_halfspaces(polyA, poseA, wTa, std::back_inserter(hs));
    PolyhedralPrimitive::get_halfspaces(polyB, poseB, wTb, std::back_inserter(hs));

    // find the interior point
    Ravelin::Origin3d ip;
    dist = -CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);
    FILE_LOG(LOG_COLDET) << "LP distance: " << dist << std::endl;

    // calculate the half-space intersection
    try
    {
      tpoly = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
    }
    catch (Ravelin::NumericalException e)
    {
      // case #2: polyhedra are kissing
      FILE_LOG(LOG_COLDET) << "qhull unable to compute intersection *volume*" << std::endl;

      // setup the normal
      std::vector<std::pair<Plane, bool> > planes;

      // init sets of vertices
      std::set<boost::shared_ptr<Polyhedron::Vertex> > svertsA, svertsB;

      // get the interior point in poseA and poseB
      Point3d ipA = Ravelin::Pose3d::transform_point(poseA, Point3d(ip, GLOBAL));
      Point3d ipB = Ravelin::Pose3d::transform_point(poseB, Point3d(ip, GLOBAL));

      // get all faces on the interior point plane from polyhedron A
      for (unsigned i=0; i< polyA.get_faces().size(); i++)
      {
        // get the plane containing the face; face thinks it is in global frame
        // but it really is in poseA's frame
        Plane plane = polyA.get_faces()[i]->get_plane();
        Ravelin::Vector3d nnew = plane.get_normal();
        nnew.pose = poseA;
        plane.set_normal(nnew);

        // compute the signed distance from the interior point
        // if signed distance is greater than zero, don't store the vertices 
        double sdist = plane.calc_signed_distance(ipA);
        if (std::fabs(sdist) > NEAR_ZERO)
          continue;

        // if the normal hasn't been set, set it
        if (planes.empty())
          planes.push_back(std::make_pair(-plane.transform(wTa), true));
        else
        {
          // normal has already been set; make sure normal is aligned with
          // set normal
          Ravelin::Vector3d normal_cand = -wTa.transform_vector(plane.get_normal());
          bool aligned_with_one = false;
          for (unsigned j=0; j< planes.size(); j++)
            if (planes[j].first.get_normal().dot(normal_cand) > 1.0 - NEAR_ZERO)
              aligned_with_one = true;
            if (!aligned_with_one)
              planes.push_back(std::make_pair(-plane.transform(wTa), true)); 
          } 
        }

      // loop over B's faces
        for (unsigned i=0; i< polyB.get_faces().size(); i++)
        {
        // get the plane containing the face; face thinks it is in global frame
        // but it really is in poseB's frame
          Plane plane = polyB.get_faces()[i]->get_plane();
          Ravelin::Vector3d nnew = plane.get_normal();
          nnew.pose = poseB;
          plane.set_normal(nnew);

        // compute the signed distance from the interior point
        // if signed distance is greater than zero, don't store the vertices 
          double sdist = plane.calc_signed_distance(ipB);
          if (std::fabs(sdist) > NEAR_ZERO)
            continue;

        // if the normal hasn't been set, set it
          if (planes.empty())
            planes.push_back(std::make_pair(plane.transform(wTb), false));
          else
          {
          // normal has already been set; make sure normal is aligned with
          // set normal
            Ravelin::Vector3d normal_cand = wTb.transform_vector(plane.get_normal());
            bool aligned_with_one = false;
            for (unsigned j=0; j< planes.size(); j++)
              if (planes[j].first.get_normal().dot(normal_cand) > 1.0 - NEAR_ZERO)
                aligned_with_one = true;
              if (!aligned_with_one)
                planes.push_back(std::make_pair(plane.transform(wTb), false)); 
            } 
          }

      // get all vertices from A and B
          const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vertsA = polyA.get_vertices();
          const std::vector<boost::shared_ptr<Polyhedron::Vertex> >& vertsB = polyB.get_vertices();

      // now we need to find *exactly* one normal
          if (planes.empty())
          {
            FILE_LOG(LOG_COLDET) << "unable to find normal!" << std::endl;
            return output_begin;
          }
          else if (planes.size() > 1)
          {
        // look for planes with non-zero distances on both sides 
            for (unsigned i=0; i< planes.size(); i++)
            {
          // get the plane
              const Plane& plane = planes[i].first;

          // see whether all vertices from the other polyhedron are all on
          // one side of the plane 
              double min_neg = 0.0;
              double max_pos = 0.0;

          // see which polyhedron plane was taken from and process the other
              if (!planes[i].second)
              { 
            // loop through all vertices from A
                for (unsigned j=0; j< vertsA.size(); j++)
                {
                  double sdist = plane.calc_signed_distance(wTa.transform_point(Point3d(vertsA[j]->o, poseA)));
                  if (sdist < min_neg)
                    min_neg = sdist;
                  if (sdist > max_pos)
                    max_pos = sdist;
                  if (-min_neg > NEAR_ZERO && max_pos > NEAR_ZERO)
                  {
                    planes[i] = planes.back();
                    planes.pop_back();
                    i--;
                    break;
                  }
                }
              }
              else
              {            
            // loop through all vertices from B
                for (unsigned j=0; j< vertsB.size(); j++)
                {
                  double sdist = plane.calc_signed_distance(wTb.transform_point(Point3d(vertsB[j]->o, poseB)));
                  if (sdist < min_neg)
                    min_neg = sdist;
                  if (sdist > max_pos)
                    max_pos = sdist;
                  if (-min_neg > NEAR_ZERO && max_pos > NEAR_ZERO)
                  {
                    planes[i] = planes.back();
                    planes.pop_back();
                    i--;
                    break;
                  }
                }
              }
            } 
          } 

      // if we still have more than one plane, normal is indeterminate 
          if (planes.size() > 1)
          {
            FILE_LOG(LOG_COLDET) << "normal is indeterminate!" << std::endl;
            return output_begin;
          }

      // get the normal
          const Ravelin::Vector3d& n = planes.front().first.get_normal(); 

      // compute the convex hulls of the two sets of vertices
          std::vector<Point3d> voA, voB, hullA, hullB;
          BOOST_FOREACH(boost::shared_ptr<Polyhedron::Vertex> v, vertsA)
          voA.push_back(wTa.transform_point(Point3d(v->o, poseA)));
          BOOST_FOREACH(boost::shared_ptr<Polyhedron::Vertex> v, vertsB)
          voB.push_back(wTb.transform_point(Point3d(v->o, poseB)));
          CompGeom::calc_convex_hull(voA.begin(), voA.end(), n, std::back_inserter(hullA));
          CompGeom::calc_convex_hull(voB.begin(), voB.end(), n, std::back_inserter(hullB));

      // compute the intersection of the convex hulls 
          std::vector<Point3d> isect;
          CompGeom::intersect_convex_polygons(hullA.begin(), hullA.end(), 
            hullB.begin(), hullB.end(), n, 
            std::back_inserter(isect));

      // create the contacts
          for (unsigned i=0; i< isect.size(); i++)
            *output_begin++ = create_contact(cgA, cgB, isect[i], n, 0.0);

          return output_begin;
        }

    // setup minimum distance, contact normal, and the contact plane offset
        double min_dist = INF;
        Ravelin::Vector3d normal(GLOBAL);
        double offset = INF;

    // get the vertices from the polyhedron
        const std::vector<Ravelin::Origin3d>& vertices = tpoly->get_vertices();

    // compute the skip distance, used to compute the normal 
        double skip_dist = std::numeric_limits<double>::max();
        const std::vector<IndexedTri>& facets = tpoly->get_facets();
        for (unsigned i=0; i< facets.size(); i++)
        {
      // setup a triangle
          Triangle tri(Point3d(vertices[facets[i].a], GLOBAL), 
           Point3d(vertices[facets[i].b], GLOBAL),
           Point3d(vertices[facets[i].c], GLOBAL));

      // get the maximum distance of all vertices
          double max_vert_dist = -std::numeric_limits<double>::max();

      // ensure that the plane containing this triangle is from 
      // polyhedron B by checking that signed distances from all vertices of B 
      // are negative
      // NOTE: this is currently an O(n) operation, but it could be turned into
      //       an O(lg N) one
          for (unsigned j=0; j< polyB.get_vertices().size(); j++)
          {
            Point3d v = wTb.transform_point(Point3d(polyB.get_vertices()[j]->o, poseB));
            double tri_dist = tri.calc_signed_dist(v);
            max_vert_dist = std::max(tri_dist, max_vert_dist);
          }

      // update the skip distance
          skip_dist = std::min(skip_dist, max_vert_dist);
        }

    // make the skip distance tolerant
        skip_dist = std::max(skip_dist, NEAR_ZERO);

    // determine the normal; the normal will be from the face that
    // yields the minimum maximum distance
        for (unsigned i=0; i< facets.size(); i++)
        {
      // setup a triangle
          Triangle tri(Point3d(vertices[facets[i].a], GLOBAL), 
           Point3d(vertices[facets[i].b], GLOBAL),
           Point3d(vertices[facets[i].c], GLOBAL));

      // ensure that the plane containing this triangle is from 
      // polyhedron B by checking that signed distances from all vertices of B 
      // are negative
      // NOTE: this is currently an O(n) operation, but it could be turned into
      //       an O(lg N) one
          bool skip = false;
          for (unsigned j=0; j< polyB.get_vertices().size(); j++)
          {
            Point3d v = wTb.transform_point(Point3d(polyB.get_vertices()[j]->o, poseB));
            double tri_dist = tri.calc_signed_dist(v);
            if (tri_dist > skip_dist + NEAR_ZERO)
            {
              skip = true;
              break;
            }
          }
          if (skip)
            continue;

      // get the reverse of the normal 
          Ravelin::Vector3d ncand = -tri.calc_normal();

      // get the offset
          double offset_cand = tri.calc_offset(ncand);

      // get the extremal point in the direction of the inverse normal
      // NOTE: this is currently an O(n) operation, but it could be turned into
      //       an O(lg N) one
          Point3d p(tpoly->find_extreme_vertex(Ravelin::Origin3d(ncand)), GLOBAL);

      // compute the distance of the extremal point from the face
          double dist = ncand.dot(p) - offset_cand;
          FILE_LOG(LOG_COLDET) << "candidate normal: " << tri.calc_normal() << " signed dist: " << dist << std::endl;
          assert(dist > -NEAR_ZERO);
          if (dist < min_dist)
          {
        // if the absolute distance is less than the minimum, we've found our
        // normal
            min_dist = dist;
            normal = -ncand;
            offset = -offset_cand; 
          }
        }
        FILE_LOG(LOG_COLDET) << "contact normal: " << normal << std::endl;

    // get each vertex, creating a contact point
        for (unsigned i=0; i< vertices.size(); i++)
        {
      // compute the interpenetration depth
          double depth = normal.dot(Point3d(vertices[i], GLOBAL)) - offset;

      // create the contact
          *output_begin++ = create_contact(cgA, cgB, Point3d(vertices[i], GLOBAL), normal, depth);
        }

        return output_begin;
      }

  // case #3: use closest features
  // get the type of the first feature
      if (boost::dynamic_pointer_cast<const Polyhedron::Vertex>(closestA))
      {
        featA = eVertex;
        boost::shared_ptr<const Polyhedron::Vertex> vA_const = boost::static_pointer_cast<const Polyhedron::Vertex>(closestA);
        vA = boost::const_pointer_cast<Polyhedron::Vertex>(vA_const);
      }
      else if (boost::dynamic_pointer_cast<const Polyhedron::Edge>(closestA))
      {
        featA = eEdge;
        boost::shared_ptr<const Polyhedron::Edge> eA_const = boost::static_pointer_cast<const Polyhedron::Edge>(closestA);
        eA = boost::const_pointer_cast<Polyhedron::Edge>(eA_const);
      }
      else 
      {
        featA = eFace;
        boost::shared_ptr<const Polyhedron::Face> fA_const = boost::static_pointer_cast<const Polyhedron::Face>(closestA);
        fA = boost::const_pointer_cast<Polyhedron::Face>(fA_const);
      }

  // get the type of the second feature
      if (boost::dynamic_pointer_cast<const Polyhedron::Vertex>(closestB))
      {
        featB = eVertex;
        boost::shared_ptr<const Polyhedron::Vertex> vB_const = boost::static_pointer_cast<const Polyhedron::Vertex>(closestB);
        vB = boost::const_pointer_cast<Polyhedron::Vertex>(vB_const);
      }
      else if (boost::dynamic_pointer_cast<const Polyhedron::Edge>(closestB))
      {
        featB = eEdge;
        boost::shared_ptr<const Polyhedron::Edge> eB_const = boost::static_pointer_cast<const Polyhedron::Edge>(closestB);
        eB = boost::const_pointer_cast<Polyhedron::Edge>(eB_const);
      }
      else 
      {
        featB = eFace;
        boost::shared_ptr<const Polyhedron::Face> fB_const = boost::static_pointer_cast<const Polyhedron::Face>(closestB);
        fB = boost::const_pointer_cast<Polyhedron::Face>(fB_const);
      }

  // find *a* pair of closest points on the two closest features
      Point3d pA0, pB0;
      if (featA == eVertex)
      {
    // vertex serves fine
        Point3d pA = Point3d(vA->o, poseA);
        pA0 = wTa.transform_point(pA);
        if (featB == eVertex)
        {
          Point3d pB = Point3d(vB->o, poseB);
          pB0 = wTb.transform_point(pB);
        }
        else if (featB == eEdge)
        {
      // find the closest point on the line segment to pA
          LineSeg3 seg;
          seg.first = wTb.transform_point(Point3d(eB->v1->o, poseB));
          seg.second = wTb.transform_point(Point3d(eB->v2->o, poseB));
          double t;
          CompGeom::calc_dist(seg, pA0, t, pB0);
        }
        else
        {
      // find the closest point on the face to pA
          Triangle t = get_triangle(fB, wTb);
          Triangle::calc_sq_dist(t, pA0, pB0);
        }
      }
      else if (featA == eEdge)
      {
        if (featB == eVertex)
        {
      // get the point in the global frame
          pB0 = wTb.transform_point(Point3d(vB->o, poseB));

      // find the closest point on the line segment to pB
          LineSeg3 seg;
          seg.first = wTa.transform_point(Point3d(eA->v1->o, poseA));
          seg.second = wTa.transform_point(Point3d(eA->v2->o, poseA));
          double t;
          CompGeom::calc_dist(seg, pB0, t, pA0);
        }
        else if (featB == eEdge)
        {
      // find two closest points on the two line segments
          LineSeg3 segA, segB;
          segA.first = wTa.transform_point(Point3d(eA->v1->o, poseA));
          segA.second = wTa.transform_point(Point3d(eA->v2->o, poseA));
          segB.first = wTb.transform_point(Point3d(eB->v1->o, poseB));
          segB.second = wTb.transform_point(Point3d(eB->v2->o, poseB));
          CompGeom::calc_closest_points(segA, segB, pA0, pB0);
        }
        else
        {
      // find a closest point on the edge to a closest point on the face
          LineSeg3 seg;
          seg.first = wTa.transform_point(Point3d(eA->v1->o, poseA));
          seg.second = wTa.transform_point(Point3d(eA->v2->o, poseA));

      // find the closest point on the face to the line segment 
          Triangle t = get_triangle(fB, wTb);
          Triangle::calc_sq_dist(t, seg, pB0, pA0);
        }
      }
      else
      {
        if (featB == eVertex)
        {
          pB0 = wTb.transform_point(Point3d(vB->o, poseB));

      // find the closest point on the face to pB
          Triangle t = get_triangle(fA, wTa);
          Triangle::calc_sq_dist(t, pB0, pA0);
        }
        else if (featB == eEdge)
        {
      // find a closest point on the edge to a closest point on the face
          LineSeg3 seg;
          seg.first = wTb.transform_point(Point3d(eB->v1->o, poseB));
          seg.second = wTb.transform_point(Point3d(eB->v2->o, poseB));

      // find the closest point on the face to the line segment 
          Triangle t = get_triangle(fA, wTa);
          Triangle::calc_sq_dist(t, seg, pA0, pB0);
        }
        else
        {
      // find closest points on each face 
          Triangle tA = get_triangle(fA, wTa);
          Triangle tB = get_triangle(fB, wTb);
          Triangle::calc_sq_dist(tA, tB, pA0, pB0);
        }
      }

  // calculate the vector between *a* pair of closest points on the two 
  // features
      Ravelin::Vector3d d = Ravelin::Vector3d::normalize(pA0 - pB0);

  // determine the two plane offsets
      double o1 = -d.dot(pA0);
      double o2 = d.dot(pB0);

  // setup points from each polyhedron
      std::vector<Point3d> pointsA, pointsB;

  // get all vertices from the first polyhedron whose projection along the
  // normal is as close to the closest point to the second polyhedron
      std::vector<boost::shared_ptr<Polyhedron::Vertex> > vAa = polyA.get_vertices();
      for (unsigned i=0; i< vAa.size(); i++)
      {
    // get the vertex in the global frame
        Point3d p = wTa.transform_point(Point3d(vAa[i]->o, poseA));

    // calculate the signed distance
        double sdist = -d.dot(p) - o1;
        if (sdist > -NEAR_ZERO)
          pointsA.push_back(p);
      }

  // get all vertices from the second polyhedron whose projection along the
  // normal is as close to the closest point to the first polyhedron
      std::vector<boost::shared_ptr<Polyhedron::Vertex> > vBb = polyB.get_vertices();
      for (unsigned i=0; i< vBb.size(); i++)
      {
    // get the vertex in the global frame
        Point3d p = wTb.transform_point(Point3d(vBb[i]->o, poseB));

    // calculate the signed distance
        double sdist = d.dot(p) - o2;
        if (sdist > -NEAR_ZERO)
          pointsB.push_back(p);
      }

  // if only one closest point on a polyhedron, no intersections are necessary 
      if (pointsA.size() == 1 || pointsB.size() == 1)
      {
    // create a contact point midway between pA0 and pB0 
        pA0 += pB0; 
        pA0 *= 0.5;
        *output_begin++ = create_contact(cgA, cgB, pA0, d, dist);
        return output_begin;
      }

  // ** from here on out, we'll have to project points to 2D
  // determine the 3D to 2D projection matrix
      Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(d);

  // determine two offsets and average them
      o1 = CompGeom::determine_3D_to_2D_offset(pointsA.front(), R);  
      o2 = CompGeom::determine_3D_to_2D_offset(pointsB.front(), R);  
      double o = (o1 + o2)*0.5;

  // project all points to the same plane
      std::vector<Ravelin::Origin2d> pointsA2, pointsB2;
      CompGeom::to_2D(pointsA.begin(), pointsA.end(), R, std::back_inserter(pointsA2));
      CompGeom::to_2D(pointsB.begin(), pointsB.end(), R, std::back_inserter(pointsB2));

  // if two closest points on each polyhedron, possible edge/edge contact;
  // compute the intersection
      if (pointsA.size() == 2 && pointsB.size() == 2)
      {
    // setup the two line segments
        LineSeg2 seg1(Point2d(pointsA2.front(), GLOBAL_2D), Point2d(pointsA2.back(), GLOBAL_2D));
        LineSeg2 seg2(Point2d(pointsB2.front(), GLOBAL_2D), Point2d(pointsB2.back(), GLOBAL_2D));

    // compute the line segment intersection
        Point2d isect, isect2;
        if (CompGeom::intersect_segs(seg1, seg2, isect, isect2) == CompGeom::eSegSegEdge)
        {
      // prepare to project back to 3D
          Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

      // project the intersections back to 3D and setup contact points there 
          Point3d p(CompGeom::to_3D(Ravelin::Origin2d(isect), RT, o), GLOBAL);
          Point3d p2(CompGeom::to_3D(Ravelin::Origin2d(isect), RT, o), GLOBAL);
          *output_begin++ = create_contact(cgA, cgB, p, d, dist);
          *output_begin++ = create_contact(cgA, cgB, p2, d, dist);
          return output_begin;
        }
        else
        {
      // treat this as a point contact
      // create a contact point midway between pA0 and pB0 
          pA0 += pB0; 
          pA0 *= 0.5;
          *output_begin++ = create_contact(cgA, cgB, pA0, d, dist);
          return output_begin;
        }
      }

  // if two closest points on one polyhedron, and three or more on another,
  // possible edge/face contact; compute the intersection 
      if (pointsA.size() == 2)
      {
    // setup the line segment
        LineSeg2 seg(Point2d(pointsA2.front(), GLOBAL_2D), Point2d(pointsA2.back(), GLOBAL_2D));

    // determine the polygon corresponding to the other
        std::vector<Ravelin::Origin2d> pointsB2_cvx;
        CompGeom::calc_convex_hull(pointsB2.begin(), pointsB2.end(), std::back_inserter(pointsB2_cvx));

    // copy as Vector2d type
        std::vector<Ravelin::Vector2d> cvx(pointsB2_cvx.size());
        for (unsigned i=0; i< pointsB2_cvx.size(); i++)
          cvx[i] = pointsB2_cvx[i];

    // intersect the segment with a convex polygon
        double te, tl;
        if (CompGeom::intersect_seg_convex_polygon(cvx.begin(), cvx.end(), seg, te, tl))
        {
      // prepare to project back to 3D
          Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

      // compute the two endpoints
          Point2d pe2 = seg.first*(1.0-te) + seg.second*te;
          Point2d pl2 = seg.first*(1.0-tl) + seg.second*tl;

      // project the points to 3D
          Point3d pe(CompGeom::to_3D(Ravelin::Origin2d(pe2), RT, o), GLOBAL);
          Point3d pl(CompGeom::to_3D(Ravelin::Origin2d(pl2), RT, o), GLOBAL);

      // create contact points
          *output_begin++ = create_contact(cgA, cgB, pe, d, dist);
          *output_begin++ = create_contact(cgA, cgB, pl, d, dist);
          return output_begin;
        }
        else
        {
      // treat this as a point contact
      // create a contact point midway between pA0 and pB0 
          pA0 += pB0; 
          pA0 *= 0.5;
          *output_begin++ = create_contact(cgA, cgB, pA0, d, dist);
          return output_begin;
        }
      }
      else if (pointsB.size() == 2)
      {
    // setup the line segment
        LineSeg2 seg(Point2d(pointsB2.front(), GLOBAL_2D), Point2d(pointsB2.back(), GLOBAL_2D));

    // determine the polygon corresponding to the other
        std::vector<Ravelin::Origin2d> pointsA2_cvx;
        CompGeom::calc_convex_hull(pointsA2.begin(), pointsA2.end(), std::back_inserter(pointsA2_cvx));

    // copy as Vector2d type
        std::vector<Ravelin::Vector2d> cvx(pointsA2_cvx.size());
        for (unsigned i=0; i< pointsA2_cvx.size(); i++)
          cvx[i] = pointsA2_cvx[i];

    // intersect the segment with a convex polygon
        double te, tl;
        if (CompGeom::intersect_seg_convex_polygon(cvx.begin(), cvx.end(), seg, te, tl))
        {
      // prepare to project back to 3D
          Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

      // compute the two endpoints
          Point2d pe2 = seg.first*(1.0-te) + seg.second*te;
          Point2d pl2 = seg.first*(1.0-tl) + seg.second*tl;

      // project the points to 3D
          Point3d pe(CompGeom::to_3D(Ravelin::Origin2d(pe2), RT, o), GLOBAL);
          Point3d pl(CompGeom::to_3D(Ravelin::Origin2d(pl2), RT, o), GLOBAL);

      // create contact points
          *output_begin++ = create_contact(cgA, cgB, pe, d, dist);
          *output_begin++ = create_contact(cgA, cgB, pl, d, dist);
          return output_begin;
        }
        else
        {
      // treat this as a point contact
      // create a contact point midway between pA0 and pB0 
          pA0 += pB0; 
          pA0 *= 0.5;
          *output_begin++ = create_contact(cgA, cgB, pA0, d, dist);
          return output_begin;
        }
      }

  // *** if we're here, we have at least three contacts on each polyhedron ***

  // compute their convex hulls
      std::vector<Ravelin::Origin2d> pointsA2_cvx, pointsB2_cvx;
      CompGeom::calc_convex_hull(pointsA2.begin(), pointsA2.end(), std::back_inserter(pointsA2_cvx));
      CompGeom::calc_convex_hull(pointsB2.begin(), pointsB2.end(), std::back_inserter(pointsB2_cvx));

  // copy as Vector2d types
      std::vector<Ravelin::Vector2d> cvxA(pointsA2_cvx.size()), cvxB(pointsB2_cvx.size());
      for (unsigned i=0; i< pointsA2_cvx.size(); i++)
        cvxA[i] = pointsA2_cvx[i];
      for (unsigned i=0; i< pointsB2_cvx.size(); i++)
        cvxB[i] = pointsB2_cvx[i];

  // compute their intersection
      std::vector<Ravelin::Vector2d> isect;
      CompGeom::intersect_convex_polygons(cvxA.begin(), cvxA.end(), cvxB.begin(), cvxB.end(), std::back_inserter(isect));

  // project back to 3D and compute contact points
      Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);
      for (unsigned i=0; i< isect.size(); i++)
      {
        Point3d p(CompGeom::to_3D(Ravelin::Origin2d(isect[i]), RT, o), GLOBAL);
        *output_begin++ = create_contact(cgA, cgB, p, d, dist);
      }
      return output_begin;

/*
  // check possibilities here
  if (pointsA.size() == 1)
  {
    // vertex from A; average all points from B
    Point3d p = pointsB.front();
    for (unsigned i=1; i< pointsB.size(); i++)
      p += pointsB[i];
    p /= pointsB.size();

    // create a contact point midway between the vertex from A and p 
    p += pointsA.front();
    p *= 0.5;
    *output_begin++ = create_contact(cgA, cgB, p, d, dist);
    return output_begin;
  }
  else if (pointsA.size() == 2)
  {
  }

      // create a contact midway between the midpoint of the edge and the vertex
      Point3d pe = (pointsB.front() + pointsB.back())*0.5;
      Point3d p = (pointsA.front() + 

  // TODO: 'for' loop here

    // get the vertex
    boost::shared_ptr<Polyhedron::Vertex> v=*vfiB;
    vfiB.advance();

    // get the projected distance from the contact plane

    Ravelin::Vector3d p(v->o, cgB->get_pose());
    Ravelin::Vector3d p0=wTB.transform_point(p);
    v2d.push_back(Point2d(CompGeom::to_2D(p0,R2D), GLOBAL_2D));
  }


  // TODO: determine the true intersection type (i.e., most general) based on 
  // this

  // now handle on a case by case basis
  if (featA == eVertex)
  {
    if (featB == eVertex)
      return find_contacts_vertex_vertex(cgA, cgB, vA, vB, dist, output_begin);
    else if (featB == eEdge)
      return find_contacts_vertex_edge(cgA, cgB, vA, eB, dist, output_begin);
    else
      return find_contacts_vertex_face(cgA, cgB, vA, fB, dist, output_begin);
  }
  else if (featA == eEdge)
  {
    if (featB == eVertex)
      return find_contacts_vertex_edge(cgB, cgA, vB, eA, dist, output_begin);
    else if (featB == eEdge)
      return find_contacts_edge_edge(cgA, cgB, eA, eB, dist, output_begin);
    else
      return find_contacts_edge_face(cgA, cgB, eA, fB, dist, output_begin);
  }
  else
  {
    if (featB == eVertex)
      return find_contacts_vertex_face(cgB, cgA, vB, fA, dist, output_begin);
    else if (featB == eEdge)
      return find_contacts_edge_face(cgB, cgA, eB, fA, dist, output_begin);
    else
      return find_contacts_face_face(cgA, cgB, fA, fB, dist, output_begin);
  }
*/
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_vertex_vertex(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> v1, boost::shared_ptr<Polyhedron::Vertex> v2, double signed_dist, OutputIterator output_begin){

  //TODO: Implement this(maybe)
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_vertex_edge(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> v, boost::shared_ptr<Polyhedron::Edge> e, double signed_dist, OutputIterator output_begin){

  //TODO: Implement this(maybe)
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_vertex_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> vA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin)
{
  // NOTE: the contact point will be *at* the vertex
  //convert the vertex point to a Vector3d
  Ravelin::Vector3d p(vA->o, cgA->get_pose());

  //converted 
  //normal will point away from B
  Ravelin::Vector3d normal(Ravelin::Origin3d(fB->get_plane().get_normal()), cgB->get_pose());

  *output_begin++ = create_contact(cgA, cgB, p, normal, signed_dist);
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_edge_edge(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> e1, boost::shared_ptr<Polyhedron::Edge> e2, double signed_dist, OutputIterator output_begin)
{
  //TODO: Implement this(maybe)
  return output_begin;
}

/// Intersects the edge against the face, assuming that both are located in the same plane
/**
 * \note contact points are created along the edge
 */
template <class OutputIterator>
 OutputIterator CCD::find_contacts_edge_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> eA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin)
 {
  // set global frame (in 2D)
  const boost::shared_ptr<const Ravelin::Pose2d> GLOBAL_2D;

  //1. calc normal of face
  Ravelin::Vector3d normal(Ravelin::Origin3d(fB->get_plane().get_normal()), cgB->get_pose());

  //We calculate the transform from cgB to GLOBAL for later use
  Ravelin::Transform3d wTB = Ravelin::Pose3d::calc_relative_pose(cgB->get_pose(),GLOBAL);

  //transform the normal to the global frame
  Ravelin::Vector3d normal0=wTB.transform_vector(normal);
  Ravelin::Matrix3d R2D = CompGeom::calc_3D_to_2D_matrix(normal0);

  //2. get all vertex from face
  Polyhedron::VertexFaceIterator vfi(fB,true);
  std::vector<Point2d> v2d;
  while(vfi.has_next()){
    boost::shared_ptr<Polyhedron::Vertex> v=*vfi;
    vfi.advance();
    Ravelin::Vector3d p(v->o, cgB->get_pose());
    Ravelin::Vector3d p0=wTB.transform_point(p);
    v2d.push_back(Point2d(CompGeom::to_2D(p0,R2D), GLOBAL_2D));
  }

  //checking if the points are CCW. If not reverse it
  if(!CompGeom::ccw(v2d.begin(),v2d.end()))
    std::reverse(v2d.begin(),v2d.end());

  //Find transform matrix for the edge
  Ravelin::Transform3d wTA = Ravelin::Pose3d::calc_relative_pose(cgA->get_pose(),GLOBAL);

  //Transforming v1
  Ravelin::Vector3d v1(Ravelin::Origin3d(eA->v1->o), cgA->get_pose());
  Ravelin::Vector3d v10=wTA.transform_point(v1);
  Ravelin::Origin2d o1_2d=CompGeom::to_2D(v10,R2D);

  //Transforming v2
  Ravelin::Vector3d v2(Ravelin::Origin3d(eA->v2->o), cgA->get_pose());
  Ravelin::Vector3d v20=wTA.transform_point(v1);
  Ravelin::Origin2d o2_2d=CompGeom::to_2D(v20,R2D);

  //creating segment
  Ravelin::Vector2d v1_2d(o1_2d, GLOBAL_2D);
  Ravelin::Vector2d v2_2d(o2_2d, GLOBAL_2D);  
  LineSeg2 line(v1_2d,v2_2d);

  //creacting result holders
  double te,tl;

  //intersecting line to face
  CompGeom::intersect_seg_convex_polygon(v2d.begin(),v2d.end(),line,te,tl);

  //5. calc 2 contact points
  Ravelin::Vector3d c1=v10*(1-te)+v20*te;
  Ravelin::Vector3d c2=v10*(1-tl)+v20*tl;

  //adding output
  *output_begin++ = create_contact(cgA, cgB, c1, normal0, signed_dist);
  *output_begin++ = create_contact(cgA, cgB, c2, normal0, signed_dist);
  return output_begin;
}

/**
 * \note it is not clear which plane (that containing fA or that containing
 *       fB) the contact points will lie upon
 */
template <class OutputIterator>
OutputIterator CCD::find_contacts_face_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Face> fA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin)
{
  //Face A:
  //calc normal of face
  Ravelin::Vector3d normal(Ravelin::Origin3d(fA->get_plane().get_normal()), cgA->get_pose());

  //We calculate the transform from cgB to GLOBAL for later use
  Ravelin::Transform3d wTA = Ravelin::Pose3d::calc_relative_pose(cgA->get_pose(),GLOBAL);

  //transform the normal to the global frame
  Ravelin::Vector3d normal_0=wTA.transform_vector(normal);

  //get all vertex from face
  Polyhedron::VertexFaceIterator vfiA(fA,true);
  std::vector<Ravelin::Vector3d> v3dA;
  while(vfiA.has_next()){
    boost::shared_ptr<Polyhedron::Vertex> v=*vfiA;
    vfiA.advance();
    Ravelin::Vector3d p(v->o, cgA->get_pose());
    Ravelin::Vector3d p0=wTA.transform_point(p);
    v3dA.push_back(p0);
  }

  //Face B:
  //We calculate the transform from cgB to GLOBAL for later use
  Ravelin::Transform3d wTB = Ravelin::Pose3d::calc_relative_pose(cgB->get_pose(),GLOBAL);

  //2. get all vertex from face
  Polyhedron::VertexFaceIterator vfiB(fB,true);
  std::vector<Ravelin::Vector3d> v3dB;
  while(vfiB.has_next()){
    boost::shared_ptr<Polyhedron::Vertex> v=*vfiB;
    vfiB.advance();
    Ravelin::Vector3d p(v->o, cgB->get_pose());
    Ravelin::Vector3d p0=wTB.transform_point(p);
    v3dB.push_back(p0);
  }

  //intersectiong 2 faces
  std::vector<Point3d> isect;
  CompGeom::intersect_convex_polygons(v3dA.begin(),v3dA.end(),v3dB.begin(),v3dB.end(),normal_0,std::back_inserter(isect));

  //adding output
  for(unsigned i=0;i<isect.size();i++)
    *output_begin++ = create_contact(cgA, cgB, isect[i], -normal_0, signed_dist);  
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{
  std::vector<Point3d> vA, vB;
  double dist;
  std::vector<Ravelin::Vector3d> n;

  // get the vertices from A and B
  cgA->get_vertices(vA);
  cgB->get_vertices(vB);

  // examine all points from A against B
  for (unsigned i=0; i< vA.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgB->calc_dist_and_normal(vA[i], n)) <= TOL)
    {
      // add the contact points
      for (unsigned j=0; j< n.size(); j++)
        *output_begin++ = create_contact(cgA, cgB, vA[i], n[j], dist);
    }
  }

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgA->calc_dist_and_normal(vB[i], n)) <= TOL)
    {
      // add the contact points
      for (unsigned j=0; j< n.size(); j++)
        *output_begin++ = create_contact(cgA, cgB, vB[i], -n[j], dist);
    }
  }

  return output_begin;
}


// find the contacts between a plane and a generic shape
template <class OutputIterator>
OutputIterator CCD::find_contacts_cylinder_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
{
  Ravelin::Vector3d normal;
  Point3d p; // this is plane

  // Set intitial value of distance to contact
  double d = std::numeric_limits<double>::infinity();

  // get the two primitives
  boost::shared_ptr<CylinderPrimitive> pA = boost::dynamic_pointer_cast<CylinderPrimitive>(cgA->get_geometry());
  boost::shared_ptr<PlanePrimitive> pB = boost::dynamic_pointer_cast<PlanePrimitive>(cgB->get_geometry());

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_cylinder_plane() entered with tolerance " << TOL << std::endl;
  FILE_LOG(LOG_COLDET) << " body A: " << cgA->get_single_body()->body_id << std::endl;
  FILE_LOG(LOG_COLDET) << " body B: " << cgB->get_single_body()->body_id << std::endl;

  // get the pose for the plane primitive
  boost::shared_ptr<const Ravelin::Pose3d> Pplane = pB->get_pose(cgB);
  // get the pose for the cylinder
  boost::shared_ptr<const Ravelin::Pose3d> Pcyl = pA->get_pose(cgA);

  ///////////////
  const double R = pA->get_radius();
  const double H = pA->get_height();
  const unsigned Y = 1;

  Ravelin::Transform3d pPc = Ravelin::Pose3d::calc_relative_pose(Pcyl,Pplane);

  // cN is the cylinder axis with respect to the plane
  Ravelin::Vector3d cN = Ravelin::Vector3d(
   Ravelin::Matrix3d(pPc.q).get_column(Y),
   Pplane);
  cN.normalize();

  // Plane Normal
  Ravelin::Vector3d n(0,1,0,Pplane);

  // Contact normal is always plane normal
  normal = Ravelin::Pose3d::transform_vector(Moby::GLOBAL,n);

  // cylinder o w.r.t. plane
  Point3d c0(pPc.x.data(),Pplane);

  double n_dot_cN = n.dot(cN);

  Ravelin::Vector3d axial_dir = cN;
  if(n_dot_cN > 0)
    axial_dir = -axial_dir;
  axial_dir.normalize();

  if(fabs(n_dot_cN) > 1.0-1e-8){
    FILE_LOG(LOG_COLDET) << " -- Cylinder axis is perpendicular to Plane" << std::endl;

    Point3d x = (H/2.0)*axial_dir + c0;

    d = x.dot(n);

    if (d > TOL)
      return o;

    int res = 4;
    for(int i=0;i<res;i++){
      Ravelin::Vector3d tan1,tan2;
      Ravelin::Vector3d::determine_orthonormal_basis(n,tan1,tan2);
      double t = (M_PI*2.0)/res * (double)i;
      Point3d p_cylinder = x + R*cos(t)*tan1 + R*sin(t)*tan2;
      p = Ravelin::Pose3d::transform_point(Moby::GLOBAL,p_cylinder);

      *o++ = create_contact(cgA, cgB, Ravelin::Pose3d::transform_point(GLOBAL, p), normal, d);
    }

  } else if(fabs(n_dot_cN) < 1e-8){
    FILE_LOG(LOG_COLDET) << " -- Cylinder axis is parallel to Plane"<< std::endl;

    Point3d x = c0 - R*n;

    d = x.dot(n);
    if (d > TOL)
      return o;

    double res[2] = {-1.0,1.0};
    for(int i=0;i<2;i++){
      double t = res[i]*H/2.0;
      Point3d p_cylinder = x + axial_dir*t;
      p = Ravelin::Pose3d::transform_point(Moby::GLOBAL,p_cylinder);

      *o++ = create_contact(cgA, cgB, Ravelin::Pose3d::transform_point(GLOBAL, p), normal, d);
    }

  } else {
    FILE_LOG(LOG_COLDET) << " -- Cylinder edge is closest to plane"<< std::endl;
    //(axis_cylinder x (n_plane x axis_cylinder))
    Ravelin::Vector3d radial_dir =
    Ravelin::Vector3d::cross(
      axial_dir,
      Ravelin::Vector3d::cross(axial_dir,n)
      );
    radial_dir.normalize();

    Point3d x = (H/2.0)*axial_dir + R*radial_dir + c0;

    d = x.dot(n);

    if (d > TOL)
      return o;

    p =  Ravelin::Pose3d::transform_point(Moby::GLOBAL,x);
//    Point3d pP = x + d*n;
//    pthis =  Ravelin::Pose3d::transform_point(Moby::GLOBAL,pP);
    *o++ = create_contact(cgA, cgB, Ravelin::Pose3d::transform_point(GLOBAL, p), normal, d);

  }
  return o;
}
// find the contacts between a plane and a generic shape
template <class OutputIterator>
OutputIterator CCD::find_contacts_sphere_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two primitives
  boost::shared_ptr<SpherePrimitive> pA = boost::dynamic_pointer_cast<SpherePrimitive>(cgA->get_geometry());
  boost::shared_ptr<PlanePrimitive> pB = boost::dynamic_pointer_cast<PlanePrimitive>(cgB->get_geometry());

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_sphere_plane() entered with tolerance " << TOL << std::endl;
  FILE_LOG(LOG_COLDET) << " body A: " << cgA->get_single_body()->body_id << std::endl;
  FILE_LOG(LOG_COLDET) << " body B: " << cgB->get_single_body()->body_id << std::endl;

  // get the poses for the two primitives
  boost::shared_ptr<const Ravelin::Pose3d> sph_pose = pA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> plane_pose = pB->get_pose(cgB);

  // get the sphere in the plane pose
  Point3d sph_c(0.0, 0.0, 0.0, sph_pose);
  Point3d sph_c_plane = Ravelin::Pose3d::transform_point(plane_pose, sph_c);

  // get the lowest point on the sphere
  double dist = sph_c_plane[Y] - pA->get_radius();

  // check the tolerance
  if (dist > TOL)
    return o;

  // setup the contact point
  Point3d p(sph_c_plane[X], 0.5*(sph_c_plane[Y] - pA->get_radius()), sph_c_plane[Z], plane_pose);

  // setup the normal
  Ravelin::Vector3d n(0.0, 1.0, 0.0, plane_pose);
  n = Ravelin::Pose3d::transform_vector(GLOBAL, n);

  // create the contact 
  *o++ = create_contact(cgA, cgB, Ravelin::Pose3d::transform_point(GLOBAL, p), n, dist);

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_sphere_plane() exited" << std::endl;

  // copy points to o
  return o;
}

// find the contacts between a plane and a generic shape
template <class OutputIterator>
OutputIterator CCD::find_contacts_plane_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
{
  std::vector<Point3d> vB;
  double dist;
  std::vector<Ravelin::Vector3d> n;

  // get the plane primitive
  boost::shared_ptr<PlanePrimitive> pA = boost::dynamic_pointer_cast<PlanePrimitive>(cgA->get_geometry());

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_plane_generic() entered with tolerance " << TOL << std::endl;
  FILE_LOG(LOG_COLDET) << " body A: " << cgA->get_single_body()->body_id << std::endl;
  FILE_LOG(LOG_COLDET) << " body B: " << cgB->get_single_body()->body_id << std::endl;

  // get the vertices from B
  cgB->get_vertices(vB);

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // compute the distance
    double dist = cgA->calc_dist_and_normal(vB[i], n);

    // see whether the point is inside the primitive
    FILE_LOG(LOG_COLDET) << "point " << vB[i] << " distance: " << dist << std::endl;
    if (dist <= TOL)
    {
      // add the contact point
      for (unsigned j=0; j< n.size(); j++)
        *o++ = create_contact(cgA, cgB, vB[i], -n[j], dist);
    }
  }

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_plane_generic() exited" << std::endl;

  // copy points to o
  return o;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_heightmap_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
{
  std::vector<Point3d> vA, vB;
  double dist;
  std::vector<Ravelin::Vector3d> n;

  // get the heightmap primitive
  boost::shared_ptr<HeightmapPrimitive> hmA = boost::dynamic_pointer_cast<HeightmapPrimitive>(cgA->get_geometry());

  // get the bounding volume for cgB
  PrimitivePtr pB = cgB->get_geometry();
  BVPtr bvB = pB->get_BVH_root(cgB);

  // get the vertices from A and B
  hmA->get_vertices(bvB, hmA->get_pose(cgA), vA);
  cgB->get_vertices(vB);

  // examine all points from A against B
  for (unsigned i=0; i< vA.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgB->calc_dist_and_normal(vA[i], n)) <= TOL)
    {
      // add the contact points
      for (unsigned j=0; j< n.size(); j++)
        *o++ = create_contact(cgA, cgB, vA[i], -n[j], dist);
    }
  }

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgA->calc_dist_and_normal(vB[i], n)) <= TOL)
    {
      // add the contact point
      for (unsigned j=0; j< n.size(); j++)
        *o++ = create_contact(cgA, cgB, vB[i], n[j], dist);
    }
  }

  // copy points to o
  return o;
}

/// Finds contacts between a sphere and a heightmap
template <class OutputIterator>
OutputIterator CCD::find_contacts_sphere_heightmap(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{

  const unsigned X = 0, Z = 2;

  // get the output iterator
  OutputIterator o = output_begin;

  // setup a vector of contacts
  std::vector<Constraint> contacts;

  // get the sphere and heightmap
  boost::shared_ptr<SpherePrimitive> sA = boost::dynamic_pointer_cast<SpherePrimitive>(cgA->get_geometry());
  boost::shared_ptr<HeightmapPrimitive> hmB = boost::dynamic_pointer_cast<HeightmapPrimitive>(cgB->get_geometry());

  // get the two poses for the primitives
  boost::shared_ptr<const Ravelin::Pose3d> pA = sA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> pB = hmB->get_pose(cgB);

  // get the transform from the sphere pose to the heightmap
  Ravelin::Transform3d T = Ravelin::Pose3d::calc_relative_pose(pA, pB);

  // transform the sphere center to the height map space
  Point3d ps_c(0.0, 0.0, 0.0, pA);
  Point3d ps_c_B = T.transform_point(ps_c);

  // get the lowest point on the sphere (toward the heightmap)
  Ravelin::Vector3d vdir(0.0, -1.0*sA->get_radius(), 0.0, pB);

  // get the lowest point on the sphere
  Point3d sphere_lowest = ps_c_B + vdir;

  // get the height of the lowest point on the sphere above the heightmap
  double min_sphere_dist = hmB->calc_height(sphere_lowest);
  if (min_sphere_dist < TOL)
  {
    // setup the contact point
    Point3d point = Ravelin::Pose3d::transform_point(GLOBAL, ps_c_B);

    // setup the normal
    Ravelin::Vector3d normal = Ravelin::Vector3d(0.0, 1.0, 0.0, pB);
    if (min_sphere_dist >= 0.0)
    {
      double gx, gz;
      hmB->calc_gradient(Ravelin::Pose3d::transform_point(pB, ps_c_B), gx, gz);
      normal = Ravelin::Vector3d(-gx, 1.0, -gz, pB);
      normal.normalize();
    }
    normal = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
    contacts.push_back(create_contact(cgA, cgB, point, normal, min_sphere_dist));
  }

  // get the corners of the bounding box in pB pose
  Point3d bv_lo = ps_c_B;
  Point3d bv_hi = ps_c_B;
  bv_lo[X] -= sA->get_radius();
  bv_hi[X] += sA->get_radius();
  bv_lo[Z] -= sA->get_radius();
  bv_hi[Z] += sA->get_radius();

  // get the heightmap width, depth, and heights
  double width = hmB->get_width();
  double depth = hmB->get_depth();
  const Ravelin::MatrixNd& heights = hmB->get_heights();

  // get the lower i and j indices
  unsigned lowi = constrain_unsigned((bv_lo[X]+width*0.5)*(heights.rows()-1)/width,heights.rows()-1);
  unsigned lowj = constrain_unsigned((bv_lo[Z]+depth*0.5)*(heights.columns()-1)/depth,heights.columns()-1);

  // get the upper i and j indices
  unsigned upi = constrain_unsigned(((bv_hi[X]+width*0.5)*(heights.rows()-1)/width)+1,heights.rows()-1);
  unsigned upj = constrain_unsigned(((bv_hi[Z]+depth*0.5)*(heights.columns()-1)/depth)+1,heights.columns()-1);

  // iterate over all points in the bounding region
  for (unsigned i=lowi; i< upi; i++)
    for (unsigned j=lowj; j< upj; j++)
    {
      // compute the point on the heightmap
      double x = -width*0.5+width*i/(heights.rows()-1);
      double z = -depth*0.5+depth*j/(heights.columns()-1);
      Point3d p(x, heights(i,j), z, pB);

      // get the distance from the primitive
      Point3d p_A = Ravelin::Pose3d::transform_point(pA, p);
      double dist = sA->calc_signed_dist(p_A);

      // ignore distance if it isn't sufficiently close
      if (dist > TOL)
        continue;

      // setup the contact point
      Point3d point = Ravelin::Pose3d::transform_point(GLOBAL, p_A);

      // setup the normal
      Ravelin::Vector3d normal = Ravelin::Vector3d(0.0, 1.0, 0.0, pB);
      if (dist >= 0.0)
      {
        double gx, gz;
        hmB->calc_gradient(Ravelin::Pose3d::transform_point(pB, p_A), gx, gz);
        normal = Ravelin::Vector3d(-gx, 1.0, -gz, pB);
        normal.normalize();
      }
      normal = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
      contacts.push_back(create_contact(cgA, cgB, point, normal, dist));
    }

  // create the normal pointing from B to A
    return std::copy(contacts.begin(), contacts.end(), o);
  }

/// Finds contacts for a convex shape and a heightmap
template <class OutputIterator>
  OutputIterator CCD::find_contacts_convex_heightmap(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
  {
    const unsigned X = 0, Y = 1, Z = 2;

  // get the output iterator
    OutputIterator o = output_begin;

  // setup a vector of contacts
    std::vector<Constraint> contacts;

  // get the convex primitive and heightmap
    PrimitivePtr sA = cgA->get_geometry();
    boost::shared_ptr<HeightmapPrimitive> hmB = boost::dynamic_pointer_cast<HeightmapPrimitive>(cgB->get_geometry());

  // get the two poses for the primitives
    boost::shared_ptr<const Ravelin::Pose3d> pA = sA->get_pose(cgA);
    boost::shared_ptr<const Ravelin::Pose3d> pB = hmB->get_pose(cgB);

  // get the transform from the primitive pose to the heightmap
    Ravelin::Transform3d T = Ravelin::Pose3d::calc_relative_pose(pA, pB);

  // intersect vertices from the convex primitive against the heightmap
    std::vector<Point3d> cverts;
    sA->get_vertices(pA, cverts);
    for (unsigned i=0; i< cverts.size(); i++)
    {
      Point3d pt = T.transform_point(cverts[i]);
      const double HEIGHT = hmB->calc_height(pt);
      if (HEIGHT < TOL)
      {
      // setup the contact point
        Point3d point = Ravelin::Pose3d::transform_point(GLOBAL, pt);

      // setup the normal
        Ravelin::Vector3d normal = Ravelin::Vector3d(0.0, 1.0, 0.0, pB);
        if (HEIGHT >= 0.0)
        {
          double gx, gz;
          hmB->calc_gradient(pt, gx, gz);
          normal = Ravelin::Vector3d(-gx, 1.0, -gz, pB);
          normal.normalize();
        }
        normal = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
        contacts.push_back(create_contact(cgA, cgB, point, normal, HEIGHT));
      }
    }

  // get the bounding volume for the primitive
    BVPtr bv = sA->get_BVH_root(cgA);

  // get the lower and upper bounds of the BV
    Point3d bv_lo = bv->get_lower_bounds();
    Point3d bv_hi = bv->get_upper_bounds();

  // create an OBB
    OBB obb;
    obb.l[X] = (bv_hi[X] - bv_lo[X])*0.5;
    obb.l[Y] = (bv_hi[Y] - bv_lo[Y])*0.5;
    obb.l[Z] = (bv_hi[Z] - bv_lo[Z])*0.5;
    obb.R = T.q;
    obb.center = Point3d(T.x, T.target);

  // get the AABB points in heightmap space
    bv_lo = obb.get_lower_bounds();
    bv_hi = obb.get_upper_bounds();

  // get the heightmap width, depth, and heights
    double width = hmB->get_width();
    double depth = hmB->get_depth();
    const Ravelin::MatrixNd& heights = hmB->get_heights();

  // get the lower i and j indices
    unsigned lowi = constrain_unsigned((bv_lo[X]+width*0.5)*(heights.rows()-1)/width,heights.rows()-1);
    unsigned lowj = constrain_unsigned((bv_lo[Z]+depth*0.5)*(heights.columns()-1)/depth,heights.columns()-1);

  // get the upper i and j indices
    unsigned upi = constrain_unsigned(((bv_hi[X]+width*0.5)*(heights.rows()-1)/width)+1,heights.rows()-1);
    unsigned upj = constrain_unsigned(((bv_hi[Z]+depth*0.5)*(heights.columns()-1)/depth)+1,heights.columns()-1);

  // iterate over all points in the bounding region
    for (unsigned i=lowi; i<= upi; i++)
      for (unsigned j=lowj; j< upj; j++)
      {
      // compute the point on the heightmap
        double x = -width*0.5+width*i/(heights.rows()-1);
        double z = -depth*0.5+depth*j/(heights.columns()-1);
        Point3d p(x, heights(i,j), z, pB);

      // get the distance from the primitive
        Point3d p_A = Ravelin::Pose3d::transform_point(pA, p);
        double dist = sA->calc_signed_dist(p_A);

      // ignore distance if it isn't sufficiently close
        if (dist > TOL)
          continue;

      // setup the contact point
        Point3d point = Ravelin::Pose3d::transform_point(GLOBAL, p_A);

      // setup the normal
        Ravelin::Vector3d normal = Ravelin::Vector3d(0.0, 1.0, 0.0, pB);
        if (dist >= 0.0)
        {
          double gx, gz;
          hmB->calc_gradient(Ravelin::Pose3d::transform_point(pB, p_A), gx, gz);
          normal = Ravelin::Vector3d(-gx, 1.0, -gz, pB);
          normal.normalize();
        }
        normal = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
        contacts.push_back(create_contact(cgA, cgB, point, normal, dist));
      }

  // create the normal pointing from B to A
      return std::copy(contacts.begin(), contacts.end(), o);
    }

/// Finds contacts for two spheres (one piece of code works for both separated and non-separated spheres)
template <class OutputIterator>
    OutputIterator CCD::find_contacts_sphere_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
    {
  // get the output iterator
      OutputIterator o = output_begin;

  // get the two spheres
      boost::shared_ptr<SpherePrimitive> sA = boost::dynamic_pointer_cast<SpherePrimitive>(cgA->get_geometry());
      boost::shared_ptr<SpherePrimitive> sB = boost::dynamic_pointer_cast<SpherePrimitive>(cgB->get_geometry());

  // setup new pose for primitive A that refers to the underlying geometry
      boost::shared_ptr<Ravelin::Pose3d> PoseA(new Ravelin::Pose3d(*sA->get_pose()));
      PoseA->rpose = cgA->get_pose();

  // setup new pose for primitive B that refers to the underlying geometry
      boost::shared_ptr<Ravelin::Pose3d> PoseB(new Ravelin::Pose3d(*sB->get_pose()));
      PoseB->rpose = cgB->get_pose();

  // get the two sphere centers in the global frame
      PoseA->update_relative_pose(GLOBAL);
      PoseB->update_relative_pose(GLOBAL);
      Point3d cA0(PoseA->x, GLOBAL);
      Point3d cB0(PoseB->x, GLOBAL);

  // determine the distance between the two spheres
      Ravelin::Vector3d d = cA0 - cB0;
      double dist = d.norm() - sA->get_radius() - sB->get_radius();
      if (dist > TOL)
        return o;

  // get the closest points on the two spheres
      Ravelin::Vector3d n = Ravelin::Vector3d::normalize(d);
      Point3d closest_A = cA0 - n*sA->get_radius();
      Point3d closest_B = cB0 + n*sB->get_radius();

  // create the contact point halfway between the closest points
      Point3d p = (closest_A + closest_B)*0.5;

  // create the normal pointing from B to A
      *o++ = create_contact(cgA, cgB, p, n, dist);

      return o;
    }

/// Gets contact points between a box and a sphere
template <class OutputIterator>
    OutputIterator CCD::find_contacts_box_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
    {
  // get the box and the sphere
      boost::shared_ptr<BoxPrimitive> bA = boost::dynamic_pointer_cast<BoxPrimitive>(cgA->get_geometry());
      boost::shared_ptr<SpherePrimitive> sB = boost::dynamic_pointer_cast<SpherePrimitive>(cgB->get_geometry());

  // get the relevant poses for both
      boost::shared_ptr<const Ravelin::Pose3d> box_pose = bA->get_pose(cgA);
      boost::shared_ptr<const Ravelin::Pose3d> sphere_pose = sB->get_pose(cgB);

  // find closest points
      Point3d psph(sphere_pose), pbox(box_pose);
      double dist = bA->calc_closest_points(sB, pbox, psph);
      if (dist > TOL)
        return o;

  // NOTE: we aren't actually finding the deepest point of interpenetration
  // from the sphere into the box...

  // if the distance between them is greater than zero, return the midpoint
  // of the two points as the contact point
      Point3d p;
      Ravelin::Vector3d normal;
      if (dist > 0.0)
      {
        Point3d psph_global = Ravelin::Pose3d::transform_point(GLOBAL, psph);
        Point3d pbox_global = Ravelin::Pose3d::transform_point(GLOBAL, pbox);
        p = (psph_global + pbox_global)*0.5;
        normal = pbox_global - psph_global;
        double nrm = normal.norm();
        if (nrm > NEAR_ZERO)
          normal /= nrm;
        else
        {
          normal = Ravelin::Pose3d::transform_vector(GLOBAL, psph);
          normal.normalize();
        }
      }
      else
      {
        p = Ravelin::Pose3d::transform_point(GLOBAL, psph);
        normal = Ravelin::Pose3d::transform_vector(GLOBAL, psph);
        normal.normalize();
      }

  // create the contact
      *o++ = create_contact(cgA, cgB, p, normal, dist);

      return o;
    }

/// Gets contact points between a torus and a plane 
template <class OutputIterator>
    OutputIterator CCD::find_contacts_torus_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
    {
      const unsigned Z = 2;
      const double EPS = NEAR_ZERO * 100.0;
      const unsigned MAX_CONTACTS = 4;

      FILE_LOG(LOG_COLDET) << "CCD::find_contacts_torus_plane(.) entered" << std::endl;

  // get the plane and torus primitives
      PrimitivePtr plane_geom = boost::dynamic_pointer_cast<Primitive>(cgB->get_geometry());
      boost::shared_ptr<TorusPrimitive> torus_geom = boost::dynamic_pointer_cast<TorusPrimitive>(cgA->get_geometry());

  // get the major and minor radii
      const double R = torus_geom->get_major_radius(); 
      const double r = torus_geom->get_minor_radius(); 

  // get the poses for the primitives
      boost::shared_ptr<const Ravelin::Pose3d> Pplane = plane_geom->get_pose(cgB);
      boost::shared_ptr<const Ravelin::Pose3d> Ptorus = torus_geom->get_pose(cgA); 

  // get the transformation from the torus's space to the plane's space
      Ravelin::Transform3d tPp = Ravelin::Pose3d::calc_relative_pose(Pplane, Ptorus);

  // Z column of rotation matrix (plane to torus)
  // is plane normal in torus frame
      Ravelin::Vector3d n_plane = tPp.transform_vector(Ravelin::Vector3d(0.0, 1.0, 0.0, Pplane));
      n_plane.normalize();

  // Convert to global coords for output
      Ravelin::Vector3d normal = Ravelin::Pose3d::transform_vector(GLOBAL,n_plane);

  // Torus axis is z-axis in torus frame
      Ravelin::Vector3d k(0,0,1,Ptorus);

  // Set intitial value of distance to contact
      double d = std::numeric_limits<double>::infinity();

  // if Torus is aligned with plane:
  // Return distance torus origin to
  // closest point on plane less pipe r
      double n_dot_k = n_plane.dot(k);
      if (std::fabs(n_dot_k) > 1.0-EPS){
    // d = depth
    // p0 = plane origin, p = plane normal
    // l0 = line origin, l = line direction

    // plane origin: plane origin in torus frame
    // line origin: torus origin in torus frame
        Point3d p0(tPp.x,Ptorus), l0(0,0,0,Ptorus);

    // plane normal: plane normal in torus frame
    // line direction: torus k axis
        Ravelin::Vector3d n = n_plane,l = k;

    // distance torus to closest point on plane is:
    // distance torus origin to closest point on plane
    // - distance torus edge to torus origin
        d = (p0 - l0).dot(n)/(l.dot(n)) - r;
        if (d > TOL)
          return o;

    // iterate over maximum number of contacts 
        for (unsigned i=0; i< MAX_CONTACTS; i++)
        { 
          double t = (double) i/MAX_CONTACTS * (2.0 * M_PI) - M_PI;
          Point3d p_torus(R*std::cos(t),R*std::sin(t),-r,Ptorus);
          Point3d point = Ravelin::Pose3d::transform_point(GLOBAL,p_torus);
          FILE_LOG(LOG_COLDET) << " -- Torus is parallel to plane"<< std::endl;
          FILE_LOG(LOG_COLDET) << "Point: "<<  point << std::endl;
          FILE_LOG(LOG_COLDET) << "Normal: "<<  normal << std::endl;
          FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
          FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;

      // create the contact
          *o++ = create_contact(cgA, cgB, point, normal, d);
        }

        return o;
      }

  //((n_plane x axis_torus) x axis_torus)
      Ravelin::Vector3d d_ring = Ravelin::Vector3d::cross(Ravelin::Vector3d::cross(n_plane,k), k);
      d_ring.normalize();

  // if Torus is _|_ with plane:
  // Return distance torus to plane less pipe r and ring R
      if(std::fabs(n_dot_k) < EPS){
    // d = depth
    // p0 = plane origin, p = plane normal
    // l0 = line origin, l = line direction

    // plane origin: plane origin in torus frame
    // line origin: torus origin in torus frame
        Point3d p0 = tPp.transform_point(Point3d(0.0, 0.0, 0.0, Pplane));
        Point3d l0(0,0,0,Ptorus);

    // plane normal: plane normal in torus frame
    // line direction: on xy-plane of torus
    //   parallel to plane normal in torus frame
        Ravelin::Vector3d n = n_plane,l = d_ring;
        d = (p0 - l0).dot(n)/(l.dot(n)) - (r+R);
        if (d > TOL)
          return o;

        Point3d p_torus = l*(d + r + R) + l0;
        Point3d point = Ravelin::Pose3d::transform_point(Moby::GLOBAL,p_torus);
        FILE_LOG(LOG_COLDET) << " -- Torus is perpendicular to plane"<< std::endl;
        FILE_LOG(LOG_COLDET) << "Point: "<<  point << std::endl;
        FILE_LOG(LOG_COLDET) << "Normal: "<<  normal << std::endl;
        FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
        FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;

    // create the contact
        *o++ = create_contact(cgA, cgB, point, normal, d);

        return o;
      }

  //   ((d_ring x axis_torus) x n_plane ) x (d_ring x axis_torus)
  //    ^ tangent to pipe   ^
  //   ^ _|_ to plane normal            ^
  //   ^ toward plane on torus pipe                             ^
      Ravelin::Vector3d d_pipe = Ravelin::Vector3d::cross(
        Ravelin::Vector3d::cross(Ravelin::Vector3d::cross(d_ring,k),n_plane),
        Ravelin::Vector3d::cross(-d_ring,k)
        );
      d_pipe.normalize();

  // d = depth
  // p0 = plane origin, p = plane normal
  // l0 = line origin, l = line direction

  // plane origin: plane origin in torus frame
  // line origin: torus origin in torus frame
      Point3d p0(tPp.x,Ptorus), l0 = R * d_ring;

  // plane normal: plane normal in torus frame
  // line direction: on xy-plane of torus
  //   parallel to plane normal in torus frame
      Ravelin::Vector3d n = n_plane,l = d_pipe;
      d = (p0 - l0).dot(n)/(l.dot(n)) - r;
      if (d > TOL)
        return o;

  //point on torus closest to plane;
      Point3d p_torus = R * d_ring + r * d_pipe;
      p_torus.pose = Ptorus;
  // TODO: find the point in the torus's space such that
  //       tPp.transform_point(.) results in the value of y closest to
  //       negative infinity
      Point3d point = Ravelin::Pose3d::transform_point(Moby::GLOBAL,p_torus);

  // create the contact
      *o++ = create_contact(cgA, cgB, point, normal, d);

      FILE_LOG(LOG_COLDET) << "Point: "<<  point << std::endl;
      FILE_LOG(LOG_COLDET) << "Normal: "<<  normal << std::endl;
      FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
      FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;

      return o;
    }

/// Does insertion sort -- custom comparison function not supported (uses operator<)
template <class BidirectionalIterator>
    void CCD::insertion_sort(BidirectionalIterator first, BidirectionalIterator last)
    {
  // safety check; exit if nothing to do
      if (first == last)
        return;

      BidirectionalIterator min = first;

  // loop
      BidirectionalIterator i = first;
      i++;
      for (; i != last; i++)
        if (*i < *min)
          min = i;

  // swap the iterators
        std::iter_swap(first, min);
        while (++first != last)
          for (BidirectionalIterator j = first; *j < *(j-1); --j)
            std::iter_swap((j-1), j);
        }

