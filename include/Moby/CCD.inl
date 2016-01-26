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

/// Finds contacts between two polyhedra
template <class OutputIterator>
OutputIterator CCD::find_contacts_polyhedron_polyhedron(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin, double TOL)
{
  const double INF = std::numeric_limits<double>::max();
  std::vector<std::pair<Ravelin::Vector3d, double> > hs;
  enum FeatureType { eNone, eVertex, eEdge, eFace };
  FeatureType featA = eNone, featB = eNone;
  boost::shared_ptr<Polyhedron::Vertex> vA, vB;
  boost::shared_ptr<Polyhedron::Edge> eA, eB;
  boost::shared_ptr<Polyhedron::Face> fA, fB;
  TessellatedPolyhedronPtr tpoly;
 
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

  // see whether to generate contacts
  if (dist > TOL)
    return output_begin; 

  // case #1: use volume of intersection
  if (dist <= 0.0)
  {
    // get the halfspaces
    PolyhedralPrimitive::get_halfspaces(polyA, poseA, wTa, std::back_inserter(hs));
    PolyhedralPrimitive::get_halfspaces(polyB, poseB, wTb, std::back_inserter(hs));

    // find the interior point
    Ravelin::Origin3d ip;
    dist = -CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);

    // calculate the half-space intersection
    try
    {
      tpoly = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
    }
    catch (NumericalException e)
    {
      // case #2: polyhedra are kissing

      // setup the normal
      Ravelin::Vector3d normal;

      // init sets of vertices
      std::set<boost::shared_ptr<Polyhedron::Vertex> > vertsA, vertsB;

      // get the interior point in poseA and poseB
      Point3d ipA = Ravelin::Pose3d::transform_point(poseA, Point3d(ip, GLOBAL));
      Point3d ipB = Ravelin::Pose3d::transform_point(poseB, Point3d(ip, GLOBAL));

      // get all faces on the interior point plane from polyhedron B; pick
      // normal from this too
      bool normal_unset = true, bad_normal_A = false;

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

        // add all vertices of the face to the set
        BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we, polyA.get_faces()[i]->e)
        {
          boost::shared_ptr<Polyhedron::Edge> e(we);
          vertsA.insert(e->v1);
          vertsA.insert(e->v2);
        }

        // if the normal hasn't been set, set it
        if (normal_unset)
        {
          normal_unset = false;
          normal = -wTa.transform_vector(plane.get_normal());
        }
        else
        {
          // normal has already been set; make sure normal is aligned with
          // set normal
          if (normal.dot(-wTa.transform_vector(plane.get_normal())) < 1.0 - NEAR_ZERO)
            bad_normal_A = true; 
        } 
      }

      // see whether we need to try to reset the normal
      if (bad_normal_A)
        normal_unset = true;

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
        if (normal_unset)
        {
          normal_unset = false;
          normal = wTb.transform_vector(plane.get_normal());
        }
        else
        {
          // normal has already been set; make sure normal is aligned with
          // set normal
          if (normal.dot(wTb.transform_vector(plane.get_normal())) < 1.0 - NEAR_ZERO)
            return output_begin; 
        } 

        // add all vertices of the face to the set
        BOOST_FOREACH(boost::weak_ptr<Polyhedron::Edge> we, polyB.get_faces()[i]->e)
        {
          boost::shared_ptr<Polyhedron::Edge> e(we);
          vertsB.insert(e->v1);
          vertsB.insert(e->v2);
        }
      }

      // if normal isn't set, quit
      if (normal_unset)
       return output_begin;

      // compute the convex hulls of the two sets of vertices
      std::vector<Point3d> voA, voB, hullA, hullB;
      BOOST_FOREACH(boost::shared_ptr<Polyhedron::Vertex> v, vertsA)
        voA.push_back(wTa.transform_point(Point3d(v->o, poseA)));
      BOOST_FOREACH(boost::shared_ptr<Polyhedron::Vertex> v, vertsB)
        voB.push_back(wTb.transform_point(Point3d(v->o, poseB)));
      CompGeom::calc_convex_hull(voA.begin(), voA.end(), normal, std::back_inserter(hullA));
      CompGeom::calc_convex_hull(voB.begin(), voB.end(), normal, std::back_inserter(hullB));

      // compute the intersection of the convex hulls 
      std::vector<Point3d> isect;
      CompGeom::intersect_convex_polygons(hullA.begin(), hullA.end(), 
                                          hullB.begin(), hullB.end(), normal, 
                                          std::back_inserter(isect));

      // create the contacts
      for (unsigned i=0; i< isect.size(); i++)
        *output_begin++ = create_contact(cgA, cgB, isect[i], normal, 0.0);

      return output_begin;
    }

    // setup minimum distance, contact normal, and the contact plane offset
    double min_dist = INF;
    Ravelin::Vector3d normal(GLOBAL);
    double offset = INF;

    // get the vertices from the polyhedron
    const std::vector<Ravelin::Origin3d>& vertices = tpoly->get_vertices();
 
    // determine the normal; the normal will be from the face that
    // yields the minimum maximum distance
    const std::vector<IndexedTri>& facets = tpoly->get_facets();
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
        if (tri.calc_signed_dist(v) > NEAR_ZERO)
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

    // get each vertex, creating a contact point
    for (unsigned i=0; i< vertices.size(); i++)
    {
      // compute the interpenetration depth
      double depth = normal.dot(Point3d(vertices[i], GLOBAL)) - offset;
      assert(depth < NEAR_ZERO);

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
OutputIterator CCD::find_contacts_vertex_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Vertex> vA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin){
  //convert the vertex point to a Vector3d
  Ravelin::Vector3d p(vA->o, cgA->get_pose());
  //converted 
  //normal will point away from B
  Ravelin::Vector3d normal(Ravelin::Origin3d(fB->get_plane().get_normal()), cgB->get_pose());

  *output_begin++ = create_contact(cgA, cgB, p, normal, signed_dist);
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_edge_edge(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> e1, boost::shared_ptr<Polyhedron::Edge> e2, double signed_dist, OutputIterator output_begin){

  //TODO: Implement this(maybe)
  return output_begin;
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_edge_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Edge> eA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin){

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
  if(!CompGeom::ccw(v2d.begin(),v2d.end())){
    std::reverse(v2d.begin(),v2d.end());
  }

  //Finding transform matrix for the edge
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


template <class OutputIterator>
OutputIterator CCD::find_contacts_face_face(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, boost::shared_ptr<Polyhedron::Face> fA, boost::shared_ptr<Polyhedron::Face> fB, double signed_dist, OutputIterator output_begin){
  
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
  for(unsigned i=0;i<isect.size();i++){
    *output_begin++ = create_contact(cgA, cgB, isect[i], -normal_0, signed_dist);  
  }
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
  std::vector<UnilateralConstraint> contacts;

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
  std::vector<UnilateralConstraint> contacts;

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
    normal = Ravelin::Vector3d::normalize(pbox_global - psph_global);
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

