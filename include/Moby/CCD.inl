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
    if (boost::dynamic_pointer_cast<BoxPrimitive>(pB))
      return find_contacts_box_box(cgA, cgB, output_begin, TOL);
    else if (boost::dynamic_pointer_cast<PlanePrimitive>(pB))
      return find_contacts_plane_generic(cgB, cgA, output_begin, TOL); 
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
    else 
      return find_contacts_plane_generic(cgA, cgB, output_begin, TOL); 
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
OutputIterator CCD::find_contacts_sphere_plane(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two primitives
  boost::shared_ptr<SpherePrimitive> pA = boost::dynamic_pointer_cast<SpherePrimitive>(cgA->get_geometry());
  boost::shared_ptr<PlanePrimitive> pB = boost::dynamic_pointer_cast<PlanePrimitive>(cgB->get_geometry());

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_sphere_plane() entered with tolerance " << TOL << std::endl;
  FILE_LOG(LOG_COLDET) << " body A: " << cgA->get_single_body()->id << std::endl;
  FILE_LOG(LOG_COLDET) << " body B: " << cgB->get_single_body()->id << std::endl;

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
 
  // check tolerance
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

  // get the bounding volume for cgB
  PrimitivePtr pB = cgB->get_geometry();
  BVPtr bvB = pB->get_BVH_root(cgB);

  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_plane_generic() entered with tolerance " << TOL << std::endl;
  FILE_LOG(LOG_COLDET) << " body A: " << cgA->get_single_body()->id << std::endl;
  FILE_LOG(LOG_COLDET) << " body B: " << cgB->get_single_body()->id << std::endl;

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
  unsigned lowi = (unsigned) ((bv_lo[X]+width*0.5)*(heights.rows()-1)/width);
  unsigned lowj = (unsigned) ((bv_lo[Z]+depth*0.5)*(heights.columns()-1)/depth);

  // get the upper i and j indices
  unsigned upi = (unsigned) ((bv_hi[X]+width*0.5)*(heights.rows()-1)/width)+1;
  unsigned upj = (unsigned) ((bv_hi[Z]+depth*0.5)*(heights.columns()-1)/depth)+1;

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
  unsigned lowi = (unsigned) ((bv_lo[X]+width*0.5)*(heights.rows()-1)/width);
  unsigned lowj = (unsigned) ((bv_lo[Z]+depth*0.5)*(heights.columns()-1)/depth);

  // get the upper i and j indices
  unsigned upi = (unsigned) ((bv_hi[X]+width*0.5)*(heights.rows()-1)/width)+1;
  unsigned upj = (unsigned) ((bv_hi[Z]+depth*0.5)*(heights.columns()-1)/depth)+1;

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
  if (dist < TOL)
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

/// Gets the distance of this box from a sphere
template <class OutputIterator>
OutputIterator CCD::find_contacts_box_box(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o, double TOL) 
{
  // get the two boxes 
  boost::shared_ptr<BoxPrimitive> bA = boost::dynamic_pointer_cast<BoxPrimitive>(cgA->get_geometry());
  boost::shared_ptr<BoxPrimitive> bB = boost::dynamic_pointer_cast<BoxPrimitive>(cgB->get_geometry());

  // get the relevant poses for both 
  boost::shared_ptr<const Ravelin::Pose3d> bA_pose = bA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> bB_pose = bB->get_pose(cgB);

  // find closest points
  Point3d pboxA(bA_pose), pboxB(bB_pose);
  double dist = CP::find_cpoint(bA, bB, bA_pose, bB_pose, pboxA, pboxB);
  if (dist > TOL)
    return o;

  // if the distance between them is greater than zero, return the midpoint
  // of the two points as the contact point
  Point3d p;
  if (dist > NEAR_ZERO)
  {
    Ravelin::Vector3d normal;
    Point3d pboxA_global = Ravelin::Pose3d::transform_point(GLOBAL, pboxA);
    Point3d pboxB_global = Ravelin::Pose3d::transform_point(GLOBAL, pboxB); 
    p = (pboxA_global + pboxB_global)*0.5; 
    normal = Ravelin::Vector3d::normalize(pboxB_global - pboxA_global);

    // create the contact
    *o++ = create_contact(cgA, cgB, p, normal, dist);
  }
  else
  {
    std::vector<Ravelin::Vector3d> normals;
    p = Ravelin::Pose3d::transform_point(GLOBAL, pboxA);
    bA->calc_dist_and_normal(pboxA, normals);

    // create the contacts
    for (unsigned i=0; i< normals.size(); i++)
      *o++ = create_contact(cgA, cgB, p, normals[i], dist);
  }
 
  // call generic find contacts find
  return find_contacts_generic(cgA, cgB, o, TOL);
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

