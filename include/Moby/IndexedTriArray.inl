/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

inline std::ostream& operator<<(std::ostream& out, const IndexedTriArray& mesh)
{
   const std::vector<Ravelin::Origin3d>& verts = mesh.get_vertices();
   const std::vector<IndexedTri>& facets = mesh.get_facets();
   BOOST_FOREACH(const IndexedTri& tri, facets)
   {
      out << "  : " << tri.a << " " << tri.b << " " << tri.c << std::endl;
      out << "  : " << verts[tri.a] << " " << verts[tri.b] << " " << verts[tri.c] << std::endl;
   }

  return out;
}

/// Creates an indexed triangle mesh from containers of vertices and facets
/**
 * \param vertices an iterator to the beginning of a container of Point3d
 *        objects
 * \param verts_end an iterator to the end of a container of Point3d objects
 * \param facets_begin an iterator to the beginning of a container of 
 *        IndexedTri objects
 * \param facets_end an iterator to the end of a container of IndexedTri objects
 */
template <class ForwardIterator1, class ForwardIterator2>
IndexedTriArray::IndexedTriArray(ForwardIterator1 verts_begin, ForwardIterator1 verts_end, ForwardIterator2 facets_begin, ForwardIterator2 facets_end)
{
  _vertices = boost::shared_ptr<const std::vector<Ravelin::Origin3d> >(new std::vector<Ravelin::Origin3d>(verts_begin, verts_end));
  _facets = boost::shared_ptr<const std::vector<IndexedTri> >(new std::vector<IndexedTri>(facets_begin, facets_end));

  // validate that indices are within range
  validate();

  // calculate incident facets and determine coplanar features
  calc_incident_facets();
  determine_coplanar_features();
}

/// Converts an indexed triangle mesh to a container of triangles
/*
 * \param output_begin an iterator to the beginning of a container of Triangle 
 *        objects
 * \return an iterator to the end of a container of Triangle objects
 */
template <class OutputIterator>
OutputIterator IndexedTriArray::get_tris(OutputIterator output_begin, boost::shared_ptr<const Ravelin::Pose3d> P) const 
{
  // get facets and vertices
  const std::vector<Ravelin::Origin3d>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();

  for (unsigned i=0; i< facets.size(); i++)
    *output_begin++ = Triangle(Point3d(vertices[facets[i].a], P), Point3d(vertices[facets[i].b], P), Point3d(vertices[facets[i].c], P));

  return output_begin;
}

/// Intersects two meshes together and returns indices of intersecting triangles
/**
 * \param mesh_a the first mesh
 * \param mesh_b the second mesh
 * \param output_begin an iterator pointing to the beginning of a container of 
 *        std::pair<unsigned, unsigned> objects (on return, container will
 *        hold indices of intersecting triangles of mesh_a and mesh_b)
 * \param exit_early if <b>true</b>, exits after first intersection detected 
 * \return an iterator pointing to the end of a container of 
 *        std::pair<unsigned, unsigned> objects (on return, container will
 *        hold indices of intersecting triangles of mesh_a and mesh_b)
 */
template <class OutputIterator>
OutputIterator IndexedTriArray::intersect(const IndexedTriArray& mesh_a, const IndexedTriArray& mesh_b, OutputIterator output_begin, boost::shared_ptr<const Ravelin::Pose3d> Pa, boost::shared_ptr<const Ravelin::Pose3d> Pb, bool exit_early)
{
  #ifdef _DEBUG_COLDET_
  std::cout << "IndexedTriArray::intersect() entered" << std::endl;
  #endif

  // form vectors of triangles
  std::vector<Triangle> tris_a, tris_b;
  mesh_a.get_tris(std::back_inserter(tris_a), Pa);
  mesh_b.get_tris(std::back_inserter(tris_b), Pb);

  // get the transform between Pa and Pb
  Ravelin::Transform3d aTb = Ravelin::Pose3d::calc_relative_pose(Pb, Pa);

  // case 1: identity transform between Pa and Pb
  if (aTb.is_identity())
  {
    // process all pairs
    for (unsigned i=0; i< tris_a.size(); i++)
    {
      // get the i'th triangle of a
      const Triangle& ta = tris_a[i];

      // loop over all triangles of b
      for (unsigned j=0; j< tris_b.size(); j++)
      {
        // get the j'th triangle of b
        const Triangle& tb = tris_b[j];

        #ifdef _DEBUG_COLDET_
        std::cout << "checking tri/tri intersection: ";
        std::cout << "  triangle 1: " << ta << std::endl;
        std::cout << "  triangle 2: " << tb << std::endl;
        std::cout << "  result? ";
        std::cout << query_intersect_tri_tri(ta, tb) << std::endl;
        #endif

        // process them against each other
        if (query_intersect_tri_tri(ta, tb))
        {
          // add the indices to the pairs
          *output_begin++ = std::make_pair(i,j);  

          // exit early if warranted
          if (exit_early)
          {
            #ifdef _DEBUG_COLDET_
            std::cout << "(exiting early since intersection detected)" << std::endl;
            std::cout << "IndexedTriArray::intersect() exited" << std::endl;
            #endif

            return output_begin;
          }
        }
      }
    }
  }
  else
  {
    // process all pairs
    for (unsigned i=0; i< tris_a.size(); i++)
    {
      // get the i'th triangle of a
      const Triangle& ta = tris_a[i];

      // loop over all triangles of b
      for (unsigned j=0; j< tris_b.size(); j++)
      {
        // get the transformed j'th triangle of b
        Triangle tb = Triangle::transform(tris_b[j], aTb);

        #ifdef _DEBUG_COLDET_
        std::cout << "checking tri/tri intersection: ";
        std::cout << "  triangle 1: " << ta << std::endl;
        std::cout << "  triangle 2: " << tb << std::endl;
        std::cout << "  result? ";
        std::cout << query_intersect_tri_tri(ta, tb) << std::endl;
        #endif

        // process them against each other
        if (query_intersect_tri_tri(ta, tb))
        {
          // add the indices to the pairs
          *output_begin++ = std::make_pair(i,j);  

          // exit early if warranted
          if (exit_early)
          {
            #ifdef _DEBUG_COLDET_
            std::cout << "(exiting early since intersection detected)" << std::endl;
            std::cout << "IndexedTriArray::intersect() exited" << std::endl;
            #endif

            return output_begin;
          }
        }
      }
    }
  }

  #ifdef _DEBUG_COLDET_
  std::cout << "IndexedTriArray::intersect() exited" << std::endl;
  #endif

  return output_begin;
}

