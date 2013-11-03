/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

template <class OutputIterator>
OutputIterator RigidBody::get_parent_links(OutputIterator begin) const
{
  BOOST_FOREACH(JointPtr j, _inner_joints)
    *begin++ = get_parent_link(j);

  return begin;
}

template <class OutputIterator>
OutputIterator RigidBody::get_child_links(OutputIterator begin) const
{
  BOOST_FOREACH(JointPtr j, _outer_joints)
    *begin++ = get_child_link(j);

  return begin;
}

