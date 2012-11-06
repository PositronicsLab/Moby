template <class BidirectionalIterator>
void Event::insertion_sort(BidirectionalIterator first, BidirectionalIterator last)
{
  // exit if nothing to do
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

template <class OutputIterator>
OutputIterator Event::get_super_bodies(OutputIterator begin) const
{
  // look for empty event
  if (event_type == Event::eNone)
    return begin;

  // look for limit event
  if (event_type == Event::eLimit)
  {
    RigidBodyPtr outboard = limit_joint->get_outboard_link();
    *begin++ = outboard->get_articulated_body();
  }
  else if (event_type == Event::eConstraint)
  {
    RigidBodyPtr outboard = constraint_joint->get_outboard_link();
    *begin++ = outboard->get_articulated_body();
  }  
  else if (event_type == Event::eContact)
  {
    SingleBodyPtr sb1 = contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_geom2->get_single_body();
    ArticulatedBodyPtr ab1 = sb1->get_articulated_body();
    ArticulatedBodyPtr ab2 = sb2->get_articulated_body();
    if (ab1)
      *begin++ = ab1;
    else
    {
      if (sb1->is_enabled())
        *begin++ = sb1;
    }
    if (ab2)
      *begin++ = ab2;
    else
    {
      if (sb2->is_enabled())
        *begin++ = sb2;
    }
  }
  else
    assert(false);

  return begin;
}

