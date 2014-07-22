/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Gets the vertex indices from selected facets of an IndexedTriArray
/**
 * \param tris the triangle mesh
 * \param fselect_begin iterator to a container of type unsigned (facet indices)
 * \param fselect_end iterator to a container of type unsigned (facet indices)
 * \param output beginning iterator to a container of type Point3d; unique
 *        vertices are copied here on return
 * \return ending iterator to a container of type Point3d; unique vertices
 *         are copied here on return
 */
template <class InputIterator, class OutputIterator>
OutputIterator TriangleMeshPrimitive::get_vertices(const IndexedTriArray& tris, InputIterator fselect_begin, InputIterator fselect_end, OutputIterator output, boost::shared_ptr<const Ravelin::Pose3d> P)
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
    *output++ = Point3d(vertices[*i], P);

  return output;
}

