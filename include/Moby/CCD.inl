/// Determines contact data between two geometries that are touching or interpenetrating 
template <class OutputIterator>
OutputIterator CCD::find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
{
  std::vector<Point3d> vA, vB;
  double dist;
  Ravelin::Vector3d n;

  // look for special cases
  PrimitivePtr pA = cgA->get_geometry();
  PrimitivePtr pB = cgB->get_geometry();
  if (boost::dynamic_pointer_cast<SpherePrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_sphere_sphere(cgA, cgB, output_begin);
    if (boost::dynamic_pointer_cast<BoxPrimitive>(pB))
      return find_contacts_box_sphere(cgB, cgA, output_begin);
    if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pB))
      return find_contacts_sphere_heightmap(cgA, cgB, output_begin);
  }
  else if (boost::dynamic_pointer_cast<BoxPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_box_sphere(cgA, cgB, output_begin);
  }
  else if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pA))
  {
    if (boost::dynamic_pointer_cast<SpherePrimitive>(pB))
      return find_contacts_sphere_heightmap(cgB, cgA, output_begin);
    else
      return find_contacts_heightmap_generic(cgA, cgB, output_begin); 
  }
  else // no special case for A
  {
    if (boost::dynamic_pointer_cast<HeightmapPrimitive>(pB))
      return find_contacts_heightmap_generic(cgB, cgA, output_begin); 
  }

  // setup list of added events
  std::list<Event> e;

  // setup the closest distance found
  double min_dist = std::numeric_limits<double>::max();

  // get the vertices from A and B
  cgA->get_vertices(vA);
  cgB->get_vertices(vB);

  // examine all points from A against B  
  for (unsigned i=0; i< vA.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgB->calc_dist_and_normal(vA[i], n))-NEAR_ZERO <= min_dist)
    {
      // see whether to throw out the old points
      if (dist-NEAR_ZERO < min_dist && min_dist > 0.0)
        e.clear();

      // setup the new minimum distance
      min_dist = std::min(min_dist, std::max(0.0, dist));

      // add the contact point
      e.push_back(create_contact(cgA, cgB, vA[i], n)); 
    }
  }

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgA->calc_dist_and_normal(vB[i], n))-NEAR_ZERO <= min_dist)
    {
      // see whether to throw out the old points
      if (dist-NEAR_ZERO < min_dist && min_dist > 0.0)
        e.clear();

      // setup the new minimum distance
      min_dist = std::min(min_dist, std::max(0.0, dist));
     
      // add the contact point
      e.push_back(create_contact(cgA, cgB, vB[i], -n)); 
    }
  }

  // copy points to o
  return std::copy(e.begin(), e.end(), output_begin);
}
      
template <class OutputIterator>
OutputIterator CCD::find_contacts_heightmap_generic(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
{
  std::vector<Point3d> vA, vB;
  double dist;
  Ravelin::Vector3d n;

  // setup list of added events
  std::list<Event> e;

  // setup the closest distance found
  double min_dist = std::numeric_limits<double>::max();

  // get the heightmap primitive
  boost::shared_ptr<HeightmapPrimitive> hmA = boost::dynamic_pointer_cast<HeightmapPrimitive>(cgA);

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
    if ((dist = cgB->calc_dist_and_normal(vA[i], n))-NEAR_ZERO <= min_dist)
    {
      // see whether to throw out the old points
      if (dist-NEAR_ZERO < min_dist && min_dist > 0.0)
        e.clear();

      // setup the new minimum distance
      min_dist = std::min(min_dist, std::max(0.0, dist));

      // add the contact point
      e.push_back(create_contact(cgA, cgB, vA[i], n)); 
    }
  }

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // see whether the point is inside the primitive
    if ((dist = cgA->calc_dist_and_normal(vB[i], n))-NEAR_ZERO <= min_dist)
    {
      // see whether to throw out the old points
      if (dist-NEAR_ZERO < min_dist && min_dist > 0.0)
        e.clear();

      // setup the new minimum distance
      min_dist = std::min(min_dist, std::max(0.0, dist));
     
      // add the contact point
      e.push_back(create_contact(cgA, cgB, vB[i], -n)); 
    }
  }

  // copy points to o
  return std::copy(e.begin(), e.end(), output_begin);
}

/// Finds contacts for a sphere and a heightmap 
template <class OutputIterator>
OutputIterator CCD::find_contacts_sphere_heightmap(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
{
  const unsigned Y = 1;

  // get the output iterator
  OutputIterator o = output_begin; 

  // get the sphere and heightmap
  boost::shared_ptr<SpherePrimitive> sA = boost::dynamic_pointer_cast<SpherePrimitive>(cgA->get_geometry());
  boost::shared_ptr<HeightmapPrimitive> hmB = boost::dynamic_pointer_cast<HeightmapPrimitive>(cgB->get_geometry());

  // get the two poses for the primitives
  boost::shared_ptr<const Ravelin::Pose3d> pA = sA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> pB = hmB->get_pose(cgB);

  // get the transform from sphere to the heightmap
  Ravelin::Transform3d T = Ravelin::Pose3d::calc_relative_pose(pB, pA);

  // transform the sphere center to the height map space
  Point3d ps_c(0.0, 0.0, 0.0, pA);
  Point3d ps_c_B = T.transform_point(ps_c);

  // compute the distance
  double y = hmB->calc_height(ps_c_B);
  double d = y - sA->get_radius();

  // see whether we can exit early
  if (d > NEAR_ZERO)
    return o;

  // setup closest point on heightmap
  Point3d pheightmap = ps_c_B;
  pheightmap[Y] = y; 

  // setup closest point on sphere
  Ravelin::Vector3d xlat(0.0, -1.0, 0.0, pB);
  if (d > 0.0)
    xlat *= sA->get_radius();
  else
    xlat *= sA->get_radius()+d;
  Point3d psphere = ps_c + Ravelin::Pose3d::transform_vector(pA, xlat);

  // setup the normal at the heightmap
  double gx, gz;
  hmB->calc_gradient(pheightmap, gx, gz);
  Ravelin::Vector3d normal(gx, 1.0, gz, pB);
  normal = Ravelin::Pose3d::transform_vector(GLOBAL, normal); 

  // create the contact point at the heightmap 
  Point3d p = Ravelin::Pose3d::transform_point(GLOBAL, pheightmap);

  // create the normal pointing from B to A
  *o++ = create_contact(cgA, cgB, p, normal); 

  return o;    
}

/// Finds contacts for two spheres (one piece of code works for both separated and non-separated spheres)
template <class OutputIterator>
OutputIterator CCD::find_contacts_sphere_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
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
  if (d.norm() - sA->get_radius() - sB->get_radius())
    return o;  

  // get the closest points on the two spheres
  Ravelin::Vector3d n = Ravelin::Vector3d::normalize(d);
  Point3d closest_A = cA0 - n*sA->get_radius();
  Point3d closest_B = cB0 + n*sB->get_radius();

  // create the contact point halfway between the closest points
  Point3d p = (closest_A + closest_B)*0.5;

  // create the normal pointing from B to A
  *o++ = create_contact(cgA, cgB, p, n); 

  return o;    
}

/// Gets the distance of this box from a sphere
template <class OutputIterator>
OutputIterator CCD::find_contacts_box_sphere(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator o) 
{
  // get the box and the sphere 
  boost::shared_ptr<BoxPrimitive> bA = boost::dynamic_pointer_cast<BoxPrimitive>(cgA->get_geometry());
  boost::shared_ptr<SpherePrimitive> sB = boost::dynamic_pointer_cast<SpherePrimitive>(cgB->get_geometry());

  // get the relevant poses for both 
  boost::shared_ptr<const Ravelin::Pose3d> box_pose = bA->get_pose(cgA);
  boost::shared_ptr<const Ravelin::Pose3d> sphere_pose = sB->get_pose(cgB);

  // get the sphere center in a's frame 
  Point3d sph_c(0.0, 0.0, 0.0, sphere_pose);
  Point3d sph_c_A = Ravelin::Pose3d::transform_point(box_pose, sph_c);

  // get the closest point
  Point3d pbox(box_pose);
  double dist = bA->calc_closest_point(sph_c_A, pbox) - sB->get_radius();
  FILE_LOG(LOG_COLDET) << "CCD::find_contacts_box_sphere(): distance is " << dist << std::endl;
  if (dist > NEAR_ZERO)
    return o;

  // compute farthest interpenetration of box inside sphere
  Point3d box_c(0.0, 0.0, 0.0, box_pose);
  Ravelin::Vector3d normal = Ravelin::Pose3d::transform_point(GLOBAL, box_c);
  normal.normalize();

  // determine closest point on the sphere 
  Point3d psph = Ravelin::Pose3d::transform_point(GLOBAL, sph_c);
  psph += normal * (sB->get_radius() + std::min(dist, 0.0));

  // if the distance is greater than zero, return the midpoint
  Point3d p;
  if (dist > 0.0)
    p = (psph+Ravelin::Pose3d::transform_point(GLOBAL, pbox))*0.5;
  else
    p = psph;

  // create the contact
  *o++ = create_contact(cgA, cgB, p, normal);

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

