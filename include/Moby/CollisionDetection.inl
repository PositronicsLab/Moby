/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Gets the set of single bodies in the collision detector
template <class OutputIterator>
OutputIterator CollisionDetection::get_single_bodies(OutputIterator output_begin) const
{
  std::list<SingleBodyPtr> sbs;
  BOOST_FOREACH(CollisionGeometryPtr cg, _geoms)
    sbs.push_back(cg->get_single_body());
  sbs.sort();
  return std::copy(sbs.begin(), std::unique(sbs.begin(), sbs.end()), output_begin);
}

/// Gets the set of rigid bodies in the collision detector
template <class OutputIterator>
OutputIterator CollisionDetection::get_rigid_bodies(OutputIterator output_begin) const
{
  std::list<RigidBodyPtr> rbs;
  BOOST_FOREACH(CollisionGeometryPtr cg, _geoms)
  {
    RigidBodyPtr rb = boost::dynamic_pointer_cast<RigidBody>(cg->get_single_body());
    if (rb)
      rbs.push_back(rb);
  }
  rbs.sort();
  return std::copy(rbs.begin(), std::unique(rbs.begin(), rbs.end()), output_begin);
}

/// Gets the set of dynamic bodies in the collision detector
template <class OutputIterator>
OutputIterator CollisionDetection::get_dynamic_bodies(OutputIterator output_begin) const
{
  std::list<DynamicBodyPtr> dbs;
  BOOST_FOREACH(CollisionGeometryPtr cg, _geoms)
  {
    SingleBodyPtr sb = cg->get_single_body();
    ArticulatedBodyPtr ab = sb->get_articulated_body();
    if (ab)
      dbs.push_back(boost::dynamic_pointer_cast<DynamicBody>(ab));
    else
      dbs.push_back(sb);
  }
  dbs.sort();
  return std::copy(dbs.begin(), std::unique(dbs.begin(), dbs.end()), output_begin);
}


