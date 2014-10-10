template <class BidirectionalIterator>
void UnilateralConstraint::insertion_sort(BidirectionalIterator first, BidirectionalIterator last)
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
OutputIterator UnilateralConstraint::get_super_bodies(OutputIterator begin) const
{
  DynamicBodyPtr db1, db2;
  unsigned nb = get_super_bodies(db1, db2);
  switch (nb)
  {
    case 0: break;
    case 1: if (db1) *begin++ = db1; break;
    case 2: if (db1) *begin++ = db1; if (db2) *begin++ = db2; break;
    default:
      assert(false);
  }
  return begin;
} 


