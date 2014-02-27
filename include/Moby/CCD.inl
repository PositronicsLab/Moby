/// Finds two contact events between two dynamic bodies
template <class OutputIterator>
OutputIterator CCD::find_contacts(DynamicBodyPtr dbA, DynamicBodyPtr dbB, OutputIterator output_begin) 
{
  // drill down to pairs of rigid bodies
  ArticulatedBodyPtr abA = boost::dynamic_pointer_cast<ArticulatedBody>(dbA);
  ArticulatedBodyPtr abB = boost::dynamic_pointer_cast<ArticulatedBody>(dbB);
  if (abA && abB)
  {
    BOOST_FOREACH(RigidBodyPtr rbA, abA->get_links())
      BOOST_FOREACH(RigidBodyPtr rbB, abB->get_links())
        output_begin = find_contacts(rbA, rbB, output_begin);

    return output_begin;
  }
  else if (abA)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rbB = boost::dynamic_pointer_cast<RigidBody>(dbB);
    BOOST_FOREACH(RigidBodyPtr rbA, abA->get_links())
      output_begin = find_contacts(rbA, rbB, output_begin);

    return output_begin;
  }
  else if (abB)
  {
    double dt = std::numeric_limits<double>::max();
    RigidBodyPtr rbA = boost::dynamic_pointer_cast<RigidBody>(dbA);
    BOOST_FOREACH(RigidBodyPtr rbB, abB->get_links())
      output_begin = find_contacts(rbA, rbB, output_begin);

    return output_begin;
  }
  else
  {
    RigidBodyPtr rbA = boost::dynamic_pointer_cast<RigidBody>(dbA);
    RigidBodyPtr rbB = boost::dynamic_pointer_cast<RigidBody>(dbB);
    return find_contacts(rbA, rbB, output_begin); 
  }
}

/// Finds two contact events between two rigid bodies 
template <class OutputIterator>
OutputIterator CCD::find_contacts(RigidBodyPtr rbA, RigidBodyPtr rbB, OutputIterator output_begin) 
{
  if (rbA == rbB || (!rbA->is_enabled() && !rbB->is_enabled()))
    return output_begin;

  BOOST_FOREACH(CollisionGeometryPtr cgA, rbA->geometries)
    BOOST_FOREACH(CollisionGeometryPtr cgB, rbB->geometries)
    {
      Point3d pA, pB;
      double dist = CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB);
      if (dist <= 0.0)
        output_begin = find_contacts_not_separated(cgA, cgB, output_begin);
      else if (dist < NEAR_ZERO)
        output_begin = find_contacts_separated(cgA, cgB, dist, output_begin);
    }

  return output_begin; 
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_separated(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, double min_dist, OutputIterator output_begin)
{
  std::vector<Point3d> vA, vB;
  Point3d pA, pB;
  Ravelin::Vector3d n;

  // look for special cases
  PrimitivePtr primA = cgA->get_geometry();
  PrimitivePtr primB = cgB->get_geometry();
  if (typeid(primA) == typeid(boost::shared_ptr<SpherePrimitive>))
  {
    if (typeid(primB) == typeid(boost::shared_ptr<SpherePrimitive>))
      return find_contacts_sphere_sphere(cgA, cgB, output_begin);
  }

  // get the output iterator
  OutputIterator o = output_begin; 

  // get the vertices from A and B
  cgA->get_vertices(vA);
  cgB->get_vertices(vB);

  // examine all points from A against B  
  for (unsigned i=0; i< vA.size(); i++)
  {
    // get the distance from the point to the primitive 
    double dist = cgB->calc_dist_and_normal(vA[i], n);

    // see whether the distance is comparable to the minimum distance
    if (dist - NEAR_ZERO <= min_dist)
      *o++ = create_contact(cgA, cgB, vA[i], n); 
  }    

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // get the distance from the point to the primitive 
    double dist = cgA->calc_dist_and_normal(vB[i], n);

    // see whether the distance is comparable to the minimum distance
    if (dist - NEAR_ZERO <= min_dist)
      *o++ = create_contact(cgA, cgB, vB[i], -n); 
  }

  return o;    
}

/// Determines contact data between two geometries that are touching or interpenetrating 
template <class OutputIterator>
OutputIterator CCD::find_contacts_not_separated(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
{
  std::vector<Point3d> vA, vB;
  Ravelin::Vector3d n;

  // look for special cases
  PrimitivePtr pA = cgA->get_geometry();
  PrimitivePtr pB = cgB->get_geometry();
  if (typeid(pA) == typeid(boost::shared_ptr<SpherePrimitive>))
  {
    if (typeid(pB) == typeid(boost::shared_ptr<SpherePrimitive>))
      return find_contacts_sphere_sphere(cgA, cgB, output_begin);
  }

  // copy the output iterator
  OutputIterator o = output_begin;

  // indicate, as of yet, no contacts have been added
  bool added = false;

  // get the vertices from A and B
  cgA->get_vertices(vA);
  cgB->get_vertices(vB);

  // examine all points from A against B  
  for (unsigned i=0; i< vA.size(); i++)
  {
    // see whether the point is inside the primitive
    if (cgB->calc_dist_and_normal(vA[i], n) <= 0.0)
    {
      *o++ = create_contact(cgA, cgB, vA[i], n); 
      added = true;
    }
  }

  // examine all points from B against A
  for (unsigned i=0; i< vB.size(); i++)
  {
    // see whether the point is inside the primitive
    if (cgA->calc_dist_and_normal(vB[i], n) <= 0.0)
    {
      *o++ = create_contact(cgA, cgB, vB[i], -n); 
      added = true;
    }
  }

  // if no contacts have been found, use the separated distance function
  if (!added)
  {
    // examine all points from A against B  
    for (unsigned i=0; i< vA.size(); i++)
    {
      // get the distance from the point to the primitive 
      double dist = cgB->calc_dist_and_normal(vA[i], n);

      // see whether the distance is comparable to the minimum distance
      if (dist <= NEAR_ZERO)
        *o++ = create_contact(cgA, cgB, vA[i], n); 
    }    

    // examine all points from B against A
    for (unsigned i=0; i< vB.size(); i++)
    {
      // get the distance from the point to the primitive 
      double dist = cgA->calc_dist_and_normal(vB[i], n);

      // see whether the distance is comparable to the minimum distance
      if (dist <= NEAR_ZERO)
        *o++ = create_contact(cgA, cgB, vB[i], -n); 
    }    
  }

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
  
  // get the closest points on the two spheres
  Ravelin::Vector3d d = cA0 - cB0;
  Ravelin::Vector3d n = Ravelin::Vector3d::normalize(d);
  Point3d closest_A = cA0 - n*sA->get_radius();
  Point3d closest_B = cB0 + n*sB->get_radius();

  // create the contact point halfway between the closest points
  Point3d p = (closest_A + closest_B)*0.5;

  // create the normal pointing from B to A
  *o++ = create_contact(cgA, cgB, p, n); 

  return o;    
}

