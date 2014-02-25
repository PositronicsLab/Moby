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
  BOOST_FOREACH(CollisionGeometryPtr cgA, rbA->geometries)
    BOOST_FOREACH(CollisionGeometryPtr cgB, rbB->geometries)
    {
      Point3d pA, pB;
      if (CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB) <= 0.0)
        output_begin = find_contacts_not_separated(cgA, cgB, output_begin);
    }

  return output_begin; 
}

template <class OutputIterator>
OutputIterator CCD::find_contacts_separated(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, OutputIterator output_begin)
{
/*
  // distance functions
  // 1) computing separating distance between bodies (for CA)
  // 2) computing interpenetration between bodies (for constraint vio)
  // 3) computing separating distance between bodies and normal (for contacts)
  //    -- for interpenetrating bodies, this includes finding points
  //       of interpenetration
 */
  std::vector<Point3d> vA, vB;
  Point3d pA, pB;
  Ravelin::Vector3d n;

  // get the output iterator
  OutputIterator o = output_begin; 

  // get the distance between the closest points 
  double min_dist = CollisionGeometry::calc_signed_dist(cgA, cgB, pA, pB);

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


