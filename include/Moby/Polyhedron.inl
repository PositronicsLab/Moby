/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Constructs a Polyhedron from iterators to a container holding Point3d objects and iterators to a container holding IndexedTri types
template <class InputIterator1, class InputIterator2>
Polyhedron::Polyhedron(InputIterator1 verts_begin, InputIterator1 verts_end, InputIterator2 facets_begin, InputIterator2 facets_end)
{
  // indicate that convexity has not been computed
  _convexity_computed = false;

  // create the mesh
  _mesh = IndexedTriArray(verts_begin, verts_end, facets_begin, facets_end);

  // calculate the bounding box
  calc_bounding_box();
}

