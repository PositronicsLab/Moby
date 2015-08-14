/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

template <class OutputIterator>
OutputIterator RigidBody::get_parent_links(OutputIterator begin) const
{
  BOOST_FOREACH(boost::shared_ptr<Ravelin::Jointd> j, _inner_joints)
    *begin++ = get_parent_link(j);

  return begin;
}

template <class OutputIterator>
OutputIterator RigidBody::get_child_links(OutputIterator begin) const
{
  BOOST_FOREACH(boost::shared_ptr<Ravelin::Jointd> j, _outer_joints)
    *begin++ = get_child_link(j);

  return begin;
}


