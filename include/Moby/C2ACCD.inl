/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Constructs the object from a collection of dynamic bodies
template <class InputIterator>
C2ACCD::C2ACCD(InputIterator begin, InputIterator end) 
{
  while (begin != end) 
    add_dynamic_body(*begin++); 
}

/// Intersects two BV leafs together and returns collision data (if any)
template <class OutputIterator>
OutputIterator C2ACCD::intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const
{
  // NOTE: if we want to speed up this static collision check (slightly),
  // we could institute an object-level map from BV leafs to triangles

  // get the primitives 
  boost::shared_ptr<const Primitive> a_primitive = geom_a->get_geometry(); 
  boost::shared_ptr<const Primitive> b_primitive = geom_b->get_geometry();

  // get the mesh data from the plugins
  const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& mdata_a = _meshes.find(a)->second;
  const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& mdata_b = _meshes.find(b)->second; 

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

/// Gets the vertex indices from selected facets of an IndexedTriArray
/**
 * \param tris the triangle mesh
 * \param fselect_begin iterator to a container of type unsigned (facet indices)
 * \param fselect_end iterator to a container of type unsigned (facet indices)
 * \param output beginning iterator to a container of type Ravelin::Vector3; unique
 *        vertices are copied here on return
 * \return ending iterator to a container of type Ravelin::Vector3; unique vertices
 *         are copied here on return
 */
template <class InputIterator, class OutputIterator>
OutputIterator C2ACCD::get_vertices(const IndexedTriArray& tris, InputIterator fselect_begin, InputIterator fselect_end, OutputIterator output)
{
  // init a list of vertices  
  std::list<unsigned> verts;

  // get facets from the triangle array
  const std::vector<IndexedTri>& facets = tris.get_facets();
  while (fselect_begin != fselect_end)
  {
    const IndexedTri& f = facets[*fselect_begin++];
    verts.push_back(f.a);
    verts.push_back(f.b);
    verts.push_back(f.c);
  }

  // sort the list so that we don't process repeated vertices
  verts.sort();
  std::list<unsigned>::const_iterator new_end = std::unique(verts.begin(), verts.end());

  // copy vertices to output
  const std::vector<Ravelin::Origin3d>& vertices = tris.get_vertices();
  for (std::list<unsigned>::const_iterator i = verts.begin(); i != new_end; i++)
    *output++ = vertices[*i];

  return output;
}

