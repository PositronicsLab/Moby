/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Constructs the object from a collection of dynamic bodies
template <class InputIterator>
GeneralizedCCD::GeneralizedCCD(InputIterator begin, InputIterator end) 
{
  eps_tolerance = std::sqrt(std::numeric_limits<double>::epsilon()); 
  while (begin != end) 
    add_dynamic_body(*begin++); 
}

/// Intersects two BV leafs together and returns collision data (if any)
template <class OutputIterator>
OutputIterator GeneralizedCCD::intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const
{
  // NOTE: if we want to speed up this static collision check (slightly),
  // we could institute an object-level map from BV leafs to triangles

  // get the primitives 
  PrimitivePtr a_primitive = geom_a->get_geometry(); 
  PrimitivePtr b_primitive = geom_b->get_geometry();

  // get the mesh data from the plugins
  const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& mdata_a = a_primitive->get_sub_mesh(a);
  const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& mdata_b = b_primitive->get_sub_mesh(b); 

  // get the meshes
  const IndexedTriArray& a_mesh = *mdata_a.first;
  const IndexedTriArray& b_mesh = *mdata_b.first;

  // get the triangles
  const std::list<unsigned>& a_tris = mdata_a.second;
  const std::list<unsigned>& b_tris = mdata_b.second;

  // do intersections
  BOOST_FOREACH(unsigned a_idx, a_tris)
  {
    // get the triangle
    Triangle ta = a_mesh.get_triangle(a_idx, aTb.target);

    // loop over all triangles in b
    BOOST_FOREACH(unsigned b_idx, b_tris)
    {
      // get the untransformed second triangle
      Triangle utb = b_mesh.get_triangle(b_idx, aTb.source);

      // transform second triangle
      Triangle tb = Triangle::transform(utb, aTb);

      // see whether triangles intersect
      if (!CompGeom::query_intersect_tri_tri(ta, tb))
        continue;

      // they do intersect, add colliding pair
      CollidingTriPair cp;
      cp.geom1 = geom_a;
      cp.geom2 = geom_b;
      cp.mesh1 = &a_mesh;
      cp.mesh2 = &b_mesh;
      cp.tri1 = a_idx;
      cp.tri2 = b_idx;
      *output_begin++ = cp;

      // see whether to stop
      if (this->mode == eFirstContact)
        return output_begin; 
    }
  }

  return output_begin;
} 

/// Does insertion sort -- custom comparison function not supported (uses operator<)
template <class BidirectionalIterator>
void GeneralizedCCD::insertion_sort(BidirectionalIterator first, BidirectionalIterator last)
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

