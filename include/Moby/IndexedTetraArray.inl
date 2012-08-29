/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

inline std::ostream& operator<<(std::ostream& out, const IndexedTetraArray& mesh)
{
   const std::vector<Vector3>& verts = mesh.get_vertices();
   const std::vector<IndexedTetra>& tetra = mesh.get_tetra();
   BOOST_FOREACH(const IndexedTetra& tet, tetra)
   {
      out << "  : " << tet.a << " " << tet.b << " " << tet.c << " " << tet.d << std::endl;
      out << "  : " << verts[tet.a] << " " << verts[tet.b] << " " << verts[tet.c] << " " << verts[tet.d] << std::endl;
   }

  return out;
}

/// Creates an indexed tetrahedral mesh from containers of vertices and tetra
/**
 * \param vertices an iterator to the beginning of a container of Vector3
 *        objects
 * \param verts_end an iterator to the end of a container of Vector3 objects
 * \param tetra_begin an iterator to the beginning of a container of 
 *        IndexedTetra objects
 * \param tetra_end an iterator to the end of a container of IndexedTetra objects
 */
template <class InputIterator1, class InputIterator2>
IndexedTetraArray::IndexedTetraArray(InputIterator1 verts_begin, InputIterator1 verts_end, InputIterator2 tetra_begin, InputIterator2 tetra_end)
{
  _vertices = boost::shared_ptr<const std::vector<Vector3> >(new std::vector<Vector3>(verts_begin, verts_end));
  _tetra = boost::shared_ptr<const std::vector<IndexedTetra> >(new std::vector<IndexedTetra>(tetra_begin, tetra_end));

  // validate that indices are within range and tetrahedra oriented correctly
  validate();
}


