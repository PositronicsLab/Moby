/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Gets all BV nodes
/**
 * The output is ordered by levels in the hierarchy.
 */
template <class OutputIterator>
OutputIterator BV::get_all_BVs(OutputIterator begin) const
{
  std::queue<BVPtr> q;
  q.push(boost::const_pointer_cast<BV>(shared_from_this()));

  while (!q.empty())
  {
    BVPtr o = q.front();
    q.pop();
    *begin++ = o;
    BOOST_FOREACH(BVPtr child, o->children)
      q.push(child);
  }

  return begin;
}

/// Gets all leaf nodes
template <class OutputIterator>
OutputIterator BV::get_all_leafs(OutputIterator begin) const
{
  std::stack<BVPtr> s;
  BOOST_FOREACH(BVPtr child, children)
    s.push(child);

  while (!s.empty())
  {
    BVPtr o = s.top();
    s.pop();
    if (o->is_leaf())
      *begin++ = o;
    else
      BOOST_FOREACH(BVPtr child, o->children)
        s.push(o);
  }

  return begin;
}

/// Intersects two BV trees; returns list of all leaf-level intersecting BVs
/**
 * \param a the root of the first BV tree
 * \param b the root of the second BV tree
 * \param aTb the transform from b's frame to a's frame (i.e., 
 *        inverse(transform(a)) * transform(b))
 * \param bTa the transform from a's frame to b's frame (i.e., 
 *        inverse(transform(b)) * transform(a))
 * \param output_begin iterator to the beginning of a list of pairs of type 
 *        pair<shared_ptr<BV>, shared_ptr<BV> >; first element in each pair
 *        comes from the first BV tree, second element comes from the second
 *        tree
 * \return iterator to the end of a list of pairs of type 
 *        pair<shared_ptr<BV>, shared_ptr<BV> >; first element in each pair
 *        comes from the first BV tree, second element comes from the second
 *        tree
 */
template <class OutputIterator>
OutputIterator BV::intersect_BV_trees(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, const Ravelin::Transform3d& bTa, OutputIterator output_begin)
{
  typedef boost::tuple<BVPtr, BVPtr, bool> BVTuple;
  std::queue<BVTuple> q;

  #ifdef _DEBUG_BV_
  std::cout << "BV::intersect_BV_trees() entered" << std::endl;
  #endif

  // add a and b to the queue, in that order
  q.push(boost::make_tuple(a, b, false));

  // drill down alternatingly until both trees are exhausted
  while (!q.empty())
  {
    // get the two nodes off of the front of the queue
    BVTuple& obbt = q.front();
    a = obbt.get<0>();
    b = obbt.get<1>();
    bool reversed = obbt.get<2>();
    const Ravelin::Transform3d& iT = (!reversed) ? aTb : bTa;
    q.pop();

    // check to see whether they intersect
    if (!BV::intersects(a.get(), b.get(), iT))
      continue;

    // they do intersect; if they are both leaves, add them to the output
    if (a->is_leaf() && b->is_leaf())
      *output_begin++ = (!reversed) ? std::make_pair(a, b) : std::make_pair(b, a);
    else
    {
      // they are not both leafs; attempt to check children of a
      if (!a->is_leaf())
        BOOST_FOREACH(BVPtr child, a->children)
          q.push(boost::make_tuple(b, child, !reversed));
      else // b _must_ not be a leaf
        BOOST_FOREACH(BVPtr child, b->children)
          q.push(boost::make_tuple(child, a, !reversed));
    }
  }

  #ifdef _DEBUG_BV_
  std::cout << "BV::intersect_BV_trees() exited" << std::endl;
  #endif

  return output_begin;
}

