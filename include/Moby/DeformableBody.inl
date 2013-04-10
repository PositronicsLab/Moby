/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Gets the vertex indices from selected facets of an IndexedTriArray
/**
 * \param tselect_begin iterator to a container of type unsigned (tetra indices)
 * \param tselect_end iterator to a container of type unsigned (tetra indices)
 * \param output beginning iterator to a container of type Ravelin::Vector3; unique
 *        vertices are copied here on return
 * \return ending iterator to a container of type Ravelin::Vector3; unique vertices
 *         are copied here on return
 */
template <class InputIterator, class OutputIterator>
OutputIterator DeformableBody::get_vertices(InputIterator tselect_begin, InputIterator tselect_end, OutputIterator output) const
{
  // init a list of vertices
  std::list<unsigned> verts;

  // get tetrahedron from the tetrahedron array
  while (tselect_begin != tselect_end)
  {
    const IndexedTetra& t = _tetrahedra[*tselect_begin++];
    verts.push_back(t.a);
    verts.push_back(t.b);
    verts.push_back(t.c);
    verts.push_back(t.d);
  }

  // sort the list so that we don't process repeated vertices
  verts.sort();
  std::list<unsigned>::const_iterator new_end = std::unique(verts.begin(), verts.end());

  // copy vertices to output
  for (std::list<unsigned>::const_iterator i = verts.begin(); i != new_end; i++)
    *output++ = _nodes[*i]->x;

  return output;
}

/// Gets the bounding boxes of all tetrahedra that a point may be within
/**
 * \param p the point
 * \param output_begin the beginning of a container of type unsigned that will
 *        hold the IDs of the tetrahedra that _may_ contain p
 */
template <class OutputIterator>
OutputIterator DeformableBody::get_tetrahedra(const Ravelin::Point3d& p, OutputIterator output_begin) const
{
  const double TOL = 1e-5;
  std::stack<BVPtr> S;
  Ravelin::Point3d closest, farthest;
  bool output = false;

  FILE_LOG(LOG_BV) << "DeformableBody::get_tetrahedra() entered" << std::endl;
  FILE_LOG(LOG_BV) << "  query point: " << p << std::endl;

  S.push(_hroot);
  while (!S.empty())
  {
    BVPtr bv = S.top();
    AABB* aabb = (AABB*) bv.get();
    S.pop();

    FILE_LOG(LOG_BV) << "  examining AABB: " << *aabb << std::endl;
    
    // if point is within, process further
    if (!AABB::outside(*aabb, p, TOL))
    {
      // if AABB is a leaf, store all tetrahedra contained within
      if (bv->is_leaf())
      {
        FILE_LOG(LOG_BV) << " -- AABB is leaf node; adding all contained tetrahedra" << std::endl;

        std::map<BVPtr, std::list<unsigned> >::const_iterator map_iter = _aabb_tetra_map.find(bv);
        assert(map_iter != _aabb_tetra_map.end());
        for (std::list<unsigned>::const_iterator i = map_iter->second.begin(); i != map_iter->second.end(); i++)
        {
          *output_begin++ = *i;
          output = true;
        }
      }
      else
      {
        FILE_LOG(LOG_BV) << " -- AABB is internal node; checking children" << std::endl;

        // process children
        BOOST_FOREACH(BVPtr bv_child, bv->children)
          S.push(bv_child);
      }
    }
    else
      FILE_LOG(LOG_BV) << " -- point is outside of AABB" << std::endl;
  }

  // if we have some output, we are done
  if (output)
  {
    FILE_LOG(LOG_BV) << "DeformableBody::get_tetrahedra() exited" << std::endl;
    return output_begin;
  }

  FILE_LOG(LOG_BV) << " -- point is outside all bounding boxes; finding closest leafs" << std::endl;

  // if the point isn't inside any of the bounding boxes, then we have to get
  // the closest leaf bounding boxes
  S.push(_hroot);
  while (!S.empty())
  {
    BVPtr bv = S.top();
    AABB* aabb = (AABB*) bv.get();
    S.pop();

    FILE_LOG(LOG_BV) << "  examining AABB: " << *aabb << std::endl;

    // if AABB is a leaf, store all tetrahedra contained within
    if (aabb->is_leaf())
    {
      FILE_LOG(LOG_BV) << " -- AABB is leaf node; adding all contained tetrahedra" << std::endl;

      std::map<BVPtr, std::list<unsigned> >::const_iterator map_iter = _aabb_tetra_map.find(bv);
      assert(map_iter != _aabb_tetra_map.end());
      for (std::list<unsigned>::const_iterator i = map_iter->second.begin(); i != map_iter->second.end(); i++)
        *output_begin++ = *i;
    }
    else
    {
      FILE_LOG(LOG_BV) << " -- AABB is internal node; checking children" << std::endl;

      // make an AABB vector from the children
      std::vector<BVPtr> bvs;
      BOOST_FOREACH(BVPtr bv_child, bv->children)
        bvs.push_back(bv_child);

      // mark all children as fit for processing initially
      std::vector<bool> process(bv->children.size(), true);

      // now, process all pairs of children
      for (unsigned i=0; i< bvs.size(); i++)
      {
        // no need reprocessing...
        if (!process[i])
          continue;

        // get the aabb
        AABB* aabb_i = (AABB*) bvs[i].get();

        // get the farthest point on AABB i to p and its squared distance
        double dist_farthest = AABB::get_farthest_point(*aabb_i, p, farthest);

        // check against all other boxes
        for (unsigned j=0; j< bvs.size(); j++)
        {
          if (!process[j] || i == j)
            continue;

          // get the aabb
          AABB* aabb_j = (AABB*) bvs[j].get();

          // get the closest point on AABB j to p
          AABB::get_closest_point(*aabb_j, p, closest);

          // get the distance from the closest point to p
          double dist_closest = (closest - p).norm_sq();

          // if the closest point on AABB j is farther than the farthest point
          // on AABB i, we don't have to process this box
          if (dist_farthest < dist_closest)
            process[j] = false;
        }
      }

      // process children
      for (unsigned i=0; i< bvs.size(); i++)
        if (process[i])
          S.push(bvs[i]);
    }
  }

  FILE_LOG(LOG_BV) << "DeformableBody::get_tetrahedra() exited" << std::endl;
  return output_begin;
}

