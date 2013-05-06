/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Intersects two BV leafs together and returns collision data (if any)
template <class OutputIterator>
OutputIterator MeshDCD::intersect_BV_leafs(BVPtr a, BVPtr b, const std::pair<Ravelin::Quatd, Ravelin::Origin3d>& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const
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
    Triangle ta = a_mesh.get_triangle(a_idx);

    // loop over all triangles in b
    BOOST_FOREACH(unsigned b_idx, b_tris)
    {
      // get the untransformed second triangle
      Triangle utb = b_mesh.get_triangle(b_idx);

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


