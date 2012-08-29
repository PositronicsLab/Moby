/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Gets all descendant collision geometries (unordered)
template <class OutputIterator>
OutputIterator CollisionGeometry::get_sub_geometries(OutputIterator begin) const
{
  std::stack<CollisionGeometryPtr> cgstack;
  
  // add all children of this to the stack
  for (unsigned i=0; i< _children.size(); i++)
    cgstack.push(_children[i]);
  
  // process them..
  while (!cgstack.empty())
  {
    // get the geometry off of the top of the stack
    CollisionGeometryPtr cg = cgstack.top();
    cgstack.pop();
    
    // save it
    *begin++ = cg;
    
    // add all children to the stack
    BOOST_FOREACH(CollisionGeometryPtr child, cg->children)
      cgstack.push(child);
  }
  
  return begin;
} // end method


