/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

template <class OutputIterator>
OutputIterator RigidBody::get_parent_links(OutputIterator begin) const
{
  BOOST_FOREACH(const InnerJointData& i, _inner_joints)
    *begin++ = RigidBodyPtr(i.parent);

  return begin;
}

template <class OutputIterator>
OutputIterator RigidBody::get_child_links(OutputIterator begin) const
{
  BOOST_FOREACH(const OuterJointData& o, _outer_joints)
    *begin++ = RigidBodyPtr(o.child);

  return begin;
}

/// Gets all collision geometries of this rigid body (including all descendants of collision geometries)
template <class OutputIterator>
OutputIterator RigidBody::get_all_collision_geometries(OutputIterator begin) const
{
  std::list<CollisionGeometryPtr> descendants;

  // process each collision geometry
  BOOST_FOREACH(CollisionGeometryPtr g, geometries)
  {
    // add the top-level i'th child geometry
    *begin++ = g;

    // add all descendants of this geometry
    g->get_sub_geometries(std::back_inserter(descendants));
    while (!descendants.empty())
    {
      *begin++ = descendants.front();
      descendants.pop_front();
    }
  }
  
  return begin;
}


