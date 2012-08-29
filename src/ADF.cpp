/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <fstream>
#include <set>
#include <cmath>
#include <queue>
#include <limits>
#include <list>
#include <stack>

#ifdef USE_INVENTOR
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoCoordinate3.h>
#endif

#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/CompGeom.h>
#include <Moby/Polyhedron.h>
#include <Moby/ADF.h>

using namespace boost;
using namespace Moby;

/// Constructs an adaptively-sampled distance field (ADF) with a recursion level of 0
/**
 * The bounds of the ADF are set to -/+ infinity and the distances are set to 
 * positive infinity.
 */
ADF::ADF()
{
  // set the bounds to +/- infinity
  set_bounds(ONES_3 * -std::numeric_limits<Real>::max(), ONES_3 * std::numeric_limits<Real>::max());
  _distances = std::vector<Real>(BOX_VERTICES, std::numeric_limits<Real>::max());
}

/// Constructs an octree with the specified parent, recursion level, and bounds
/**
 * The bounds of the ADF are set to -/+ infinity and the distances are set to 
 * positive infinity.
 */
ADF::ADF(shared_ptr<ADF> parent, const Vector3& lo_bounds, const Vector3& hi_bounds)
{
  set_bounds(lo_bounds, hi_bounds);
  _distances = std::vector<Real>(BOX_VERTICES, std::numeric_limits<Real>::max());
  _parent = parent;
}

/// Constructs an ADF with the specified parent and vertices
/**
 * The vertices are expected to be in the following order:
 *        6---7
 *       /|  /|
 *      3---5 |
 *      | 2-|-4
 *      |/  |/
 *      0---1
 */
ADF::ADF(shared_ptr<ADF> parent, const std::vector<Vector3ConstPtr>& vertices)
{
  _parent = parent;
  _lo_bounds = *vertices.front();
  _hi_bounds = *vertices.back();
  _vertices = vertices;
}

/// Returns the index of the subvolume being occupied
/**
 * \return an index in the range [0,7]; assertion fails if the point is outside
 *         of the volume
 */
unsigned ADF::get_sub_volume_idx(const Vector3& point) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the bounds of the ADF
  Vector3 half = 0.5*(_lo_bounds + _hi_bounds);

  // Subdivision order
  //      6---7
  //     /|  /|
  //    3---5 |
  //    | 2-|-4
  //    |/  |/
  //    0---1
  
  if (point[X] < half[X])
  {
    if (point[Y] < half[Y])
    {
      if (point[Z] < half[Z])
        return 0;
      else
        return 2;
    }
    else
    {
      if (point[Z] < half[Z])
        return 3;
      else
        return 6;
    }
  }
  else
  {
    if (point[Y] < half[Y])
    {
      if (point[Z] < half[Z])
        return 1;
      else
        return 4;
    }
    else
    {
      if (point[Z] < half[Z])
        return 5;
      else
        return 7;
    }
  }

  assert(false);
  return std::numeric_limits<unsigned>::max();
}

/// Returns NULL if empty, ADF of leaf if occupied
shared_ptr<ADF> ADF::is_cell_occupied(const Vector3& point) const
{
  const unsigned THREE_D = 3;

  // first test whether point is in box
  if (_parent.expired())
    for (unsigned i=0; i< THREE_D; i++)
      if (point[i] < _lo_bounds[i] || point[i] > _hi_bounds[i])
      {
        FILE_LOG(LOG_ADF) << "ADF::is_cell_occupied() - " << point << " is outside bounding box!" << std::endl;
        FILE_LOG(LOG_ADF) << "  lower bound: " << _lo_bounds << std::endl;
        FILE_LOG(LOG_ADF) << "  upper bound: " << _hi_bounds << std::endl;

        return shared_ptr<ADF>();
      }

  // if there are no children, return this
  if (_children.empty())
    return const_pointer_cast<ADF>(shared_from_this());
  else
  {
    unsigned subvolume = get_sub_volume_idx(point);

    FILE_LOG(LOG_ADF) << "ADF::is_cell_occupied() - point is inside sub volume: " << subvolume << std::endl;

    return _children[subvolume]->is_cell_occupied(point);
  }
}

/// Gets the bounds for this ADF cell
void ADF::get_bounds(Vector3& lo, Vector3& hi) const
{
  lo = _lo_bounds;
  hi = _hi_bounds;
}

/// (Re)sets the bounds for this ADF cell
void ADF::set_bounds(const Vector3& lo_bound, const Vector3& hi_bound)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // set the bounds
  _lo_bounds = lo_bound;
  _hi_bounds = hi_bound;
  
  // set the corner vertices of this ADF cell
  // The corner vertices are set as follows, where '0' and '7' are the vertices
  // corresponding to the lower and upper bounds of the cell, respectively.
  //      6---7
  //     /|  /|
  //    3---5 |
  //    | 2-|-4
  //    |/  |/
  //    0---1

  // reset the vertices 
  _vertices.clear();
  _vertices.push_back(Vector3Ptr(new Vector3(_lo_bounds[X], _lo_bounds[Y], _lo_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_hi_bounds[X], _lo_bounds[Y], _lo_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_lo_bounds[X], _lo_bounds[Y], _hi_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_lo_bounds[X], _hi_bounds[Y], _lo_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_hi_bounds[X], _lo_bounds[Y], _hi_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_hi_bounds[X], _hi_bounds[Y], _lo_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_lo_bounds[X], _hi_bounds[Y], _hi_bounds[Z])));
  _vertices.push_back(Vector3Ptr(new Vector3(_hi_bounds[X], _hi_bounds[Y], _hi_bounds[Z])));

  // reset the ADF
  reset();
}

//// Resets the octree
void ADF::reset()
{
  // recursively destroy all child nodes
  for (unsigned i=0; i< _children.size(); i++)
    _children[i]->reset();
  _children.clear();
}

/// Sets the vector of distances for this cell
/**
 * \param distances a 8-dimension vector of distances corresponding to each
 *        vertex of the cell.
 * The distances must correspond to the following:
 *
 *       6---7
 *     /|  /|
 *    3---5 |
 *    | 2-|-4
 *    |/  |/
 *    0---1
 *
 * where '0' and '7' correspond to the bounds of this cell.
 */
void ADF::set_distances(const std::vector<Real>& distances)
{
  assert(distances.size() == BOX_VERTICES);
  _distances = distances;
}

/// Distance function for a triangle mesh
Real ADF::trimesh_distance_function(const Vector3& pt, void* data)
{
  // get the Polyhedron 
  Polyhedron& poly = *(Polyhedron*) data;

  // compute the distance
  return poly.calc_signed_distance(pt);
}

/// Builds an ADF from a polyhedron
shared_ptr<ADF> ADF::build_ADF(Polyhedron& poly, unsigned max_recursion, Real epsilon, Real max_pos_dist, Real max_neg_dist)
{
  // determine the bounding box for the triangle mesh
  std::pair<Vector3, Vector3> bb = poly.get_bounding_box_corners();

  // build the ADF
  return build_ADF(bb.first, bb.second, &trimesh_distance_function, max_recursion, epsilon, max_pos_dist, max_neg_dist, (void*) &poly);
}

/// Builds an ADF using a bounding box and distance function
/**
 * \param lo the lower bounds
 * \param hi the upper bounds
 * \param dfn the distance function 
 * \param max_recursion the maximum recursion level within the ADF octree
 * \param epsilon the tolerance below which subdivision stops
 * \param max_pos_dist the maximum positive (external) distance to build the 
 *         ADF; if max_pos_dist is negative, then the maximum positive distance
 *         is computed to be 1% of the diagonal of the bounding box
 * \param max_neg_dist the maximum negative (internal) distance to build the 
 *         ADF; default value is infinity
 * \return a shared pointer to the constructed ADF octree
 */
shared_ptr<ADF> ADF::build_ADF(const Vector3& lo, const Vector3& hi, Real (*dfn)(const Vector3&, void*), unsigned max_recursion, Real epsilon, Real max_pos_dist, Real max_neg_dist, void* data)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real COMPUTED_EXTRA = 0.01;

  FILE_LOG(LOG_ADF) << "building ADF with focus on iso-surface" << std::endl;
  FILE_LOG(LOG_ADF) << "  triangle mesh bounds: " << lo << " / " << hi << std::endl;

  // if the maximum positive distance is negative, compute the maximum positive
  // distance
  if (max_pos_dist < 0.0)
  {
    Vector3 diag(hi[X] - lo[X], hi[Y] - lo[Y], hi[Z] - lo[Z]);
    max_pos_dist = diag.norm() * COMPUTED_EXTRA;
  }

  // we'll make each dimension of the box slightly bigger to better represent 
  // the iso-surface
  const Real INV_SQRT_3 = 1.0/std::sqrt(3.0);
  Vector3 enlarge = Vector3(1.0, 1.0, 1.0)*INV_SQRT_3*max_pos_dist;
  Vector3 new_lo = lo - enlarge;
  Vector3 new_hi = hi + enlarge;

  FILE_LOG(LOG_ADF) << "  root ADF bounds: " << lo << " / " << hi << std::endl;

  // create the root of the ADF
  shared_ptr<ADF> root(new ADF);
  root->set_bounds(new_lo, new_hi);
  root->set_distances(dfn, data);

  // setup a queue for processing
  std::queue<shared_ptr<ADF> > q;
  q.push(root);

  while (!q.empty())
  {
    // get the element off of the front of the queue
    shared_ptr<ADF> cell = q.front();
    q.pop();

    // if it's not possible to subdivide the current cell any more, continue
    if (cell->get_recursion_level() == max_recursion)
      continue;

    // if cell is completely inside, continue
    bool completely_inside = true;
    for (unsigned i=0; i< BOX_VERTICES; i++)
      if (cell->get_distances()[i] > -max_neg_dist)
      {
        completely_inside = false;
        break;
      }
    if (completely_inside)
      continue;

    // get samples over this cell
    std::vector<Vector3> samples;
    cell->get_samples(samples);

    // otherwise, check whether the cell's distances within desired tolerance
    Real mean_err = 0;
    bool within_tol = true;
    for (unsigned i=0; i< samples.size(); i++)
    {              
      Real calc_dist = cell->calc_signed_distance(samples[i]);
      Real true_dist = dfn(samples[i], data);
      FILE_LOG(LOG_ADF) << " sample: " << samples[i] << std::endl;
      FILE_LOG(LOG_ADF) << "   true distance: " << true_dist << std::endl;
      FILE_LOG(LOG_ADF) << "   calculated distance: " << calc_dist << std::endl;
      mean_err += std::fabs(calc_dist - true_dist);
      if (std::fabs(calc_dist - true_dist) > epsilon)
      {
        within_tol = false;
        break;
      }
    }

    mean_err /= samples.size();
    FILE_LOG(LOG_ADF) << " ADF cell at level " << cell->get_recursion_level() << " mean error: " << mean_err << std::endl;

    // if the cell is within tolerance, don't need to subdivide
    if (within_tol)
      continue;

    // subdivide the cell
    cell->subdivide(dfn, data);

    // add children to the queue for processing
    const std::vector<shared_ptr<ADF> >& children = cell->get_children();
    for (unsigned i=0; i< OCT_CHILDREN; i++)
      q.push(children[i]);
  
    FILE_LOG(LOG_ADF) << "subdivided cell: " << std::endl << *cell;
    FILE_LOG(LOG_ADF) << "children: " << std::endl;
    for (unsigned i=0; i< OCT_CHILDREN; i++)
    {
      FILE_LOG(LOG_ADF) << "  child " << i << ": " << std::endl;
      FILE_LOG(LOG_ADF) << *children[i] << std::endl;
    }
  }

  return root;
}

/// Determines whether the given point is within this ADF's bounding box
bool ADF::contains(const Vector3& point) const
{
  const unsigned THREE_D = 3;

  for (unsigned i=0; i< THREE_D; i++)
    if (point[i] < _lo_bounds[i] || point[i] > _hi_bounds[i])
      return false;

  return true;
}

/// Intersects two ADFs and returns the new ADF
/**
 * \param adf1 the first ADF
 * \param adf2 the second ADF
 * \param epsilon the tolerance to the true distance below which further recursion
 *        is terminated
 * \param recursion_limit the maximum number of levels of the ADF octree
 * \return a smart pointer to the new ADF
 * \note closest triangles to samples in ADF cells are not preserved
 */
shared_ptr<ADF> ADF::intersect(shared_ptr<ADF> adf1, shared_ptr<ADF> adf2, Real epsilon, unsigned recursion_limit)
{
  const unsigned THREE_D = 3;

  // setup a map of processed vertices
  std::map<Vector3ConstPtr, Real> processed;

  // determine the intersection of the bounding box for the two ADFs
  Vector3 lo, hi, lo1, hi1, lo2, hi2;
  adf1->get_bounds(lo1, hi1);
  adf2->get_bounds(lo2, hi2);
  for (unsigned i=0; i< THREE_D; i++)
  {
    lo[i] = std::max(lo1[i], lo2[i]);
    hi[i] = std::min(hi1[i], hi2[i]);
  }

  FILE_LOG(LOG_ADF) << "ADF::intersect() called" << std::endl;
  FILE_LOG(LOG_ADF) << "  bounds for ADF 1: " << lo1 << " / " << hi1 << std::endl;
  FILE_LOG(LOG_ADF) << "  bounds for ADF 2: " << lo2 << " / " << hi2 << std::endl;
  FILE_LOG(LOG_ADF) << "  bounds for new ADF: " << lo << " / " << hi << std::endl;

  // create the root of the new ADF
  shared_ptr<ADF> root(new ADF);
  root->set_bounds(lo, hi);

  // set the distances on the root node
  std::vector<Real> distances(BOX_VERTICES);
  const std::vector<Vector3ConstPtr>& root_vertices = root->get_vertices();
  for (unsigned i=0; i< BOX_VERTICES; i++)
  {
    // compute the distance
    distances[i] = calc_max_distance(adf1, adf2, *root_vertices[i]);

    // save the distance
    processed[root_vertices[i]] = distances[i];
  }
  root->set_distances(distances);

  // setup a queue for processing
  std::queue<shared_ptr<ADF> > q;
  q.push(root);

  // process until no more processing necessary
  while (!q.empty())
  {
    // get the ADF cell off of the front of the queue
    shared_ptr<ADF> cell = q.front();
    q.pop();

    // if it's not possible to subdivide the current cell any more, continue
    if (cell->get_recursion_level() == recursion_limit)
      continue;

    // get the bounds on this cell
    Vector3 lo, hi;
    cell->get_bounds(lo, hi);

    // get samples over this cell
    std::vector<Vector3> samples;
    cell->get_samples(samples);

    // otherwise, check whether the cell's distances are within tolerance
    bool within_tol = true;
    for (unsigned i=0; within_tol && i< samples.size(); i++)
      if (std::fabs(cell->calc_signed_distance(samples[i]) - calc_max_distance(adf1, adf2, samples[i])) > epsilon)
        within_tol = false;

    // if the cell is within tolerance, don't need to subdivide
    if (within_tol)
      continue;

    // otherwise, subdivide the cell
    cell->subdivide();

    // compute distances for all of the cell's children and add the cell's
    // children to the queue for processing
    const std::vector<shared_ptr<ADF> >& children = cell->get_children();
    for (unsigned i=0; i< OCT_CHILDREN; i++)
    {
      // get the vertices for the i'th child
      const std::vector<Vector3ConstPtr>& vertices = children[i]->get_vertices();

      // process all vertices
      for (unsigned j=0; j< BOX_VERTICES; j++)
      {
        // if the j'th vertex has been processed already, save the distance
        std::map<Vector3ConstPtr, Real>::iterator k = processed.find(vertices[j]);
        if (k != processed.end())
          distances[j] = k->second;
        else
          distances[j] = calc_max_distance(adf1, adf2, *vertices[j]);

        // save the computed distances
        processed[vertices[j]] = distances[j];
      }

      // set the distances for the i'th child
      children[i]->set_distances(distances);

      // add the i'th child to the queue for processing
      q.push(children[i]);
    }
  }

  return root;
}

/// Determines the normal to the surface at a point
/**
 * \note this uses the method defined by Frisken et al. [2000]
 */
Vector3 ADF::determine_normal(const Vector3& point) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the cell that contains this point
  shared_ptr<ADF> cell = is_cell_occupied(point);
  assert(cell);

  // get the bounds of the cell
  Vector3 lo, hi;
  cell->get_bounds(lo, hi);

  // calculate the signed distance at the point projected to all six facets
  // of the cell
  Real left = cell->calc_signed_distance(Vector3(lo[X], point[Y], point[Z]));
  Real right = cell->calc_signed_distance(Vector3(hi[X], point[Y], point[Z]));
  Real up = cell->calc_signed_distance(Vector3(point[X], hi[Y], point[Z]));
  Real down = cell->calc_signed_distance(Vector3(point[X], lo[Y], point[Z]));
  Real front = cell->calc_signed_distance(Vector3(point[X], point[Y], hi[Z]));
  Real back = cell->calc_signed_distance(Vector3(point[X], point[Y], lo[Z]));

  // compute the gradient at this point
  return Vector3::normalize(Vector3(right - left, up - down, front - back));
}

/// Gets all leaf nodes in the sub-tree rooted at this cell
/**
 * \note the list is not cleared before this operation
 * \note cells are not inserted into the list in any order
 */
void ADF::get_all_leaf_nodes(std::list<shared_ptr<ADF> >& leafs) const
{
  // if this is a leaf node, add itself to the list
  if (is_leaf())
    leafs.push_back(const_pointer_cast<ADF>(shared_from_this()));
  else
    for (unsigned i=0; i< OCT_CHILDREN; i++)
      _children[i]->get_all_leaf_nodes(leafs);
}

/// Gets all cells in the sub-tree rooted at this cell
/**
 * \note the list is not cleared before this operation
 * \note the list is ordered as a pre-order
 */
void ADF::get_all_cells(std::list<shared_ptr<ADF> >& cells) const
{
  // add this node to the list
  cells.push_back(const_pointer_cast<ADF>(shared_from_this()));

  // process children recursively, if this is not a leaf
  if (!is_leaf())
    for (unsigned i=0; i< OCT_CHILDREN; i++)
      _children[i]->get_all_cells(cells);
}

/// Determines a sample on the iso-surface of the ADF
bool ADF::generate_iso_sample(Vector3& sample, Real epsilon) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // if this is not a leaf node, find a leaf node
  if (!is_leaf())
  {
    // get the vector of children, and permute it randomly
    std::vector<shared_ptr<ADF> > children_copy = _children;
    std::random_shuffle(children_copy.begin(), children_copy.end());

    // attempt to sample from each child
    for (unsigned i=0; i< OCT_CHILDREN; i++)
      if (children_copy[i]->generate_iso_sample(sample, epsilon))
        return true;

    // if we made it here, we were unable to sample
    return false;
  }

  // get the distances and vertices, and shuffle them
  std::vector<Vector3ConstPtr> vertices = _vertices;
  std::vector<Real> distances = _distances;
  std::random_shuffle(vertices.begin(), vertices.end());
  std::random_shuffle(distances.begin(), distances.end());

  // this is a leaf node, see whether it crosses the iso-surface
  bool positive_found = false, negative_found = false;
  for (unsigned i=0; i< BOX_VERTICES && !(positive_found && negative_found); i++)
  {
    if (distances[i] < 0)
      negative_found = true;
    else if (distances[i] > 0)
      positive_found = true;
    else
    {
      // if we're here, the distance is zero, and we've found the iso-surface!
      sample = *vertices[i];
      return true;
    }
  }

  // if there is not a positive / negative crossing, this leaf doesn't cross the
  // iso-surface
  if (!(positive_found && negative_found))
    return false;

  // setup some things for mathematica generated code
  Real q1 = _distances[0];
  Real q2 = _distances[1];
  Real q3 = _distances[2];
  Real q4 = _distances[3];
  Real q5 = _distances[4];
  Real q6 = _distances[5];
  Real q7 = _distances[6];
  Real q8 = _distances[7];
  Real x1x = _lo_bounds[X];
  Real x1y = _lo_bounds[Y];
  Real x1z = _lo_bounds[Z];
  Real x8x = _hi_bounds[X];
  Real x8y = _hi_bounds[Y];
  Real x8z = _hi_bounds[Z];

  // it does cross the iso-surface, solve for the point randomly
  Real px, py, pz;

  // set all values random initially, we'll analytically determine the
  // dependent variable 
  px = (Real) rand() * (_hi_bounds[X] - _lo_bounds[X])/RAND_MAX + _lo_bounds[X];
  py = (Real) rand() * (_hi_bounds[Y] - _lo_bounds[Y])/RAND_MAX + _lo_bounds[Y];
  pz = (Real) rand() * (_hi_bounds[Z] - _lo_bounds[Z])/RAND_MAX + _lo_bounds[Z];

  // determine the denominator of px (solving analytically)
  Real denom = py*pz*q1 - py*pz*q2 - py*pz*q3 - py*pz*q4 + 
       py*pz*q5 + py*pz*q6 + py*pz*q7 - py*pz*q8 - 
       pz*q1*x1y + pz*q2*x1y + pz*q3*x1y - pz*q5*x1y - 
       py*q1*x1z + py*q2*x1z + py*q4*x1z - py*q6*x1z + 
       q1*x1y*x1z - q2*x1y*x1z + pz*q4*x8y - pz*q6*x8y - 
       pz*q7*x8y + pz*q8*x8y - q4*x1z*x8y + q6*x1z*x8y + 
       py*q3*x8z - py*q5*x8z - py*q7*x8z + py*q8*x8z - 
       q3*x1y*x8z + q5*x1y*x8z + q7*x8y*x8z - q8*x8y*x8z;

  // if the denominator is not zero, finish solving for px
  if (std::fabs(denom) > 1e-12)
  {
    // we've found valid values for py and pz; determine the numerator
    Real num = py*pz*q1*x1x - py*pz*q3*x1x - 
        py*pz*q4*x1x + py*pz*q7*x1x - pz*q1*x1x*x1y + 
        pz*q3*x1x*x1y - py*q1*x1x*x1z + py*q4*x1x*x1z + 
        q1*x1x*x1y*x1z - py*pz*q2*x8x + py*pz*q5*x8x + 
        py*pz*q6*x8x - py*pz*q8*x8x + pz*q2*x1y*x8x - 
        pz*q5*x1y*x8x + py*q2*x1z*x8x - py*q6*x1z*x8x - 
        q2*x1y*x1z*x8x + pz*q4*x1x*x8y - pz*q7*x1x*x8y - 
        q4*x1x*x1z*x8y - pz*q6*x8x*x8y + pz*q8*x8x*x8y + 
        q6*x1z*x8x*x8y + py*q3*x1x*x8z - py*q7*x1x*x8z - 
        q3*x1x*x1y*x8z - py*q5*x8x*x8z + py*q8*x8x*x8z + 
        q5*x1y*x8x*x8z + q7*x1x*x8y*x8z - q8*x8x*x8y*x8z;

    // determine px
    px = num/denom;
  }
  else
  {
    // failed to solve for px, try solving for py
    Real denom = px*pz*q1 - px*pz*q2 - px*pz*q3 - px*pz*q4 + 
        px*pz*q5 + px*pz*q6 + px*pz*q7 - px*pz*q8 - 
        pz*q1*x1x + pz*q3*x1x + pz*q4*x1x - pz*q7*x1x - 
        px*q1*x1z + px*q2*x1z + px*q4*x1z - px*q6*x1z + 
        q1*x1x*x1z - q4*x1x*x1z + pz*q2*x8x - pz*q5*x8x - 
        pz*q6*x8x + pz*q8*x8x - q2*x1z*x8x + q6*x1z*x8x + 
        px*q3*x8z - px*q5*x8z - px*q7*x8z + px*q8*x8z - 
        q3*x1x*x8z + q7*x1x*x8z + q5*x8x*x8z - q8*x8x*x8z;

    if (std::fabs(denom) > 1e-12)
    {
      // we've found valid values for px and pz; determine the numerator
      Real num = px*pz*q1*x1y - px*pz*q2*x1y - 
        px*pz*q3*x1y + px*pz*q5*x1y - pz*q1*x1x*x1y + 
        pz*q3*x1x*x1y - px*q1*x1y*x1z + px*q2*x1y*x1z + 
        q1*x1x*x1y*x1z + pz*q2*x1y*x8x - pz*q5*x1y*x8x - 
        q2*x1y*x1z*x8x - px*pz*q4*x8y + px*pz*q6*x8y + 
        px*pz*q7*x8y - px*pz*q8*x8y + pz*q4*x1x*x8y - 
        pz*q7*x1x*x8y + px*q4*x1z*x8y - px*q6*x1z*x8y - 
        q4*x1x*x1z*x8y - pz*q6*x8x*x8y + pz*q8*x8x*x8y + 
        q6*x1z*x8x*x8y + px*q3*x1y*x8z - px*q5*x1y*x8z - 
        q3*x1x*x1y*x8z + q5*x1y*x8x*x8z - px*q7*x8y*x8z + 
        px*q8*x8y*x8z + q7*x1x*x8y*x8z - q8*x8x*x8y*x8z;

      // determine py
      py = num/denom;
    }
    else
    {
      // last possibility: solve for pz
      Real denom = px*py*q1 - px*py*q2 - px*py*q3 - px*py*q4 + 
        px*py*q5 + px*py*q6 + px*py*q7 - px*py*q8 - 
        py*q1*x1x + py*q3*x1x + py*q4*x1x - py*q7*x1x - 
        px*q1*x1y + px*q2*x1y + px*q3*x1y - px*q5*x1y + 
        q1*x1x*x1y - q3*x1x*x1y + py*q2*x8x - py*q5*x8x - 
        py*q6*x8x + py*q8*x8x - q2*x1y*x8x + q5*x1y*x8x + 
        px*q4*x8y - px*q6*x8y - px*q7*x8y + px*q8*x8y - 
        q4*x1x*x8y + q7*x1x*x8y + q6*x8x*x8y - q8*x8x*x8y;

      // the denominator must not be equal to zero
      assert(std::fabs(denom) > 1e-12);

      // determine the numerator
      Real num = px*py*q1*x1z - px*py*q2*x1z - 
        px*py*q4*x1z + px*py*q6*x1z - py*q1*x1x*x1z + 
        py*q4*x1x*x1z - px*q1*x1y*x1z + px*q2*x1y*x1z + 
        q1*x1x*x1y*x1z + py*q2*x1z*x8x - py*q6*x1z*x8x - 
        q2*x1y*x1z*x8x + px*q4*x1z*x8y - px*q6*x1z*x8y - 
        q4*x1x*x1z*x8y + q6*x1z*x8x*x8y - px*py*q3*x8z + 
        px*py*q5*x8z + px*py*q7*x8z - px*py*q8*x8z + 
        py*q3*x1x*x8z - py*q7*x1x*x8z + px*q3*x1y*x8z - 
        px*q5*x1y*x8z - q3*x1x*x1y*x8z - py*q5*x8x*x8z + 
        py*q8*x8x*x8z + q5*x1y*x8x*x8z - px*q7*x8y*x8z + 
        px*q8*x8y*x8z + q7*x1x*x8y*x8z - q8*x8x*x8y*x8z;

      // determine pz
      pz = num/denom;
    }
  }

  // set the sample
  sample = Vector3(px, py, pz);

  // check the distance
  if (std::fabs(calc_signed_distance(sample)) > epsilon)
  {
    FILE_LOG(LOG_ADF) << "sampling failed: signed distance = " << calc_signed_distance(sample) << std::endl;
    return false;
  }

  return true;
}

/// Checks distance to a point using two ADFs and returns the maximum signed distance
/**
 * If the point is inside both ADFs, the maximum distance is returned.  If it is
 * only inside one ADF, that distance is returned.  If it is inside neither
 * ADF, +infinity is returned.
 */
Real ADF::calc_max_distance(shared_ptr<ADF> adf1, shared_ptr<ADF> adf2, const Vector3& point)
{
  // check to see whether the point is inside the ADFs
  bool inside1 = adf1->contains(point);
  bool inside2 = adf2->contains(point);

  // if it's not inside either, return +infinity
  if (!inside1 && !inside2)
    return std::numeric_limits<Real>::max();

  // if it's only inside one, return that one
  if (inside1 && !inside2)
    return adf1->calc_signed_distance(point);
  else if (inside2 && !inside1)
    return adf2->calc_signed_distance(point);
  else
    return std::max(adf1->calc_signed_distance(point), adf2->calc_signed_distance(point));
}

/// Simplifies an ADF by coalescing cells
void ADF::simplify(Real epsilon)
{
  // if there are no children, cannot simplify
  if (_children.empty())
    return;  

  // call simplify on the child nodes first
  if (!_children.empty())
    for (unsigned i=0; i< _children.size(); i++)
      _children[i]->simplify(epsilon);

  // if one of the children is not a leaf, cannot simplify further
  for (unsigned i=0; i< _children.size(); i++)
    if (!_children[i]->is_leaf())
      return;

  // we can possibly coalesce the children; check distances at 19 points
  std::map<Vector3ConstPtr, Real> check;
  for (unsigned i=0; i< OCT_CHILDREN; i++)
  {
    const std::vector<Vector3ConstPtr>& vertices = _children[i]->get_vertices();
    const std::vector<Real>& distances = _children[i]->get_distances();
    for (unsigned j=0; j< BOX_VERTICES; j++)
      check[vertices[j]] = distances[j];
  }

  // get the set of vertices of this cell (sorted)
  std::set<Vector3ConstPtr> vertex_set(_vertices.begin(), _vertices.end());
  
  // check distances at all points that don't already exist in this cell
  for (std::map<Vector3ConstPtr, Real>::iterator i = check.begin(); i != check.end(); i++)
  {
    // if the vertex exists in this cell, don't check it
    if (vertex_set.find(i->first) != vertex_set.end())
      continue;

    // if the distance is above the tolerance, we aren't able to coalesce 
    if (std::fabs(tri_linear_interp(_vertices, _distances, *i->first) - i->second) > epsilon)
      return; 
  }

  // if we made it here, we can coalesce!
  _children.clear();
}

/// Computes the signed distance using this ADF cell using trilinear interpolation
Real ADF::calc_signed_distance(const Vector3& point) const
{
  // get the leaf cell that contains this point
  shared_ptr<ADF> leaf = is_cell_occupied(point);

  // verify that the leaf is non-NULL
  assert(leaf);

  // get the vector of vertices
  const std::vector<Vector3ConstPtr>& vertices = leaf->get_vertices();  

  // get the vector of distances
  const std::vector<Real>& distances = leaf->get_distances();

  // perform trilinear interpolation
  return tri_linear_interp(vertices, distances, point);
}

/// Gets 19 points sampled over this cell
/**
 * The 19 points are: the center of the cell, the center of each face (6), and
 * the center of each edge (12).
 * \note the passed vector is not cleared
 */
void ADF::get_samples(std::vector<Vector3>& samples) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // add the center to the samples
  samples.push_back(Vector3(_lo_bounds + _hi_bounds)/2);

  // compute midpoints of the bounds
  Real mid_x = (_lo_bounds[X] + _hi_bounds[X])/2;
  Real mid_y = (_lo_bounds[Y] + _hi_bounds[Y])/2;
  Real mid_z = (_lo_bounds[Z] + _hi_bounds[Z])/2;

  // add the center of each face
  samples.push_back(Vector3(mid_x, mid_y, _lo_bounds[Z]));
  samples.push_back(Vector3(mid_x, mid_y, _hi_bounds[Z]));
  samples.push_back(Vector3(_lo_bounds[X], mid_y, mid_z));
  samples.push_back(Vector3(_hi_bounds[X], mid_y, mid_z));
  samples.push_back(Vector3(mid_x, _lo_bounds[Y], mid_z));
  samples.push_back(Vector3(mid_x, _hi_bounds[Y], mid_z));

  // add the center of each edge
  samples.push_back(Vector3(mid_x,_lo_bounds[Y],_lo_bounds[Z]));  
  samples.push_back(Vector3(mid_x,_lo_bounds[Y],_hi_bounds[Z]));  
  samples.push_back(Vector3(mid_x,_hi_bounds[Y],_lo_bounds[Z]));  
  samples.push_back(Vector3(mid_x,_hi_bounds[Y],_hi_bounds[Z]));  
  samples.push_back(Vector3(_lo_bounds[X],mid_y,_lo_bounds[Z]));  
  samples.push_back(Vector3(_lo_bounds[X],mid_y,_hi_bounds[Z]));  
  samples.push_back(Vector3(_hi_bounds[X],mid_y,_lo_bounds[Z]));  
  samples.push_back(Vector3(_hi_bounds[X],mid_y,_hi_bounds[Z]));  
  samples.push_back(Vector3(_lo_bounds[X],_lo_bounds[Y],mid_z));  
  samples.push_back(Vector3(_lo_bounds[X],_hi_bounds[Y],mid_z));  
  samples.push_back(Vector3(_hi_bounds[X],_lo_bounds[Y],mid_z));  
  samples.push_back(Vector3(_hi_bounds[X],_hi_bounds[Y],mid_z));  
}

/// Gets the OpenInventor scenegraph representation of the ADF as a point set
SoSeparator* ADF::render() const
{
  #ifdef USE_INVENTOR
  SoSeparator* separator = new SoSeparator;
  render(separator);
  return separator;
  #else
  return NULL;
  #endif
}

/// Gets the OpenInventor scenegraph representation of the ADF as a point set
void ADF::render(SoSeparator* separator) const
{
  #ifdef USE_INVENTOR
  const unsigned POINTS_PER_DIM = 10;
  const unsigned X = 0, Y = 1, Z = 2;

  // if this node is not a leaf, render all children...
  if (!_children.empty())
    for (unsigned i=0; i< _children.size(); i++)
      _children[i]->render(separator);
  else
  {
    // if the distances do not have both positive and negative numbers, quit
    bool positive_found = false, negative_found = false;
    for (unsigned i=0; i< _distances.size(); i++)
      if (_distances[i] > 0.0)
        positive_found = true;
      else if (_distances[i] < 0.0)
        negative_found = true;
      else if (_distances[i] == 0.0)
      {
        positive_found = true;
        negative_found = true;
        break;
      }
    if (!(positive_found && negative_found))
      return;

    unsigned idx = 0;
    SoCoordinate3* coords = new SoCoordinate3;
    Real px = _lo_bounds[X], py = _lo_bounds[Y], pz = _lo_bounds[Z];
    const Real xlen = (*_vertices[1])[X] - (*_vertices[0])[X];
    const Real ylen = (*_vertices[3])[Y] - (*_vertices[0])[Y];
    const Real zlen = (*_vertices[2])[Z] - (*_vertices[0])[Z];
    const Real xinc = xlen / (POINTS_PER_DIM-1);
    const Real yinc = ylen / (POINTS_PER_DIM-1);
    const Real zinc = zlen / (POINTS_PER_DIM-1);
    for (unsigned i=0; i< POINTS_PER_DIM; i++, px += xinc)
    {
      py = _lo_bounds[Y];
      for (unsigned j=0; j< POINTS_PER_DIM; j++, py += yinc)
      {
        pz = _lo_bounds[Z];
        for (unsigned k=0; k< POINTS_PER_DIM; k++, pz+= zinc)
        {
          px = (Real) rand() / RAND_MAX * (_hi_bounds[X] - _lo_bounds[X]) + _lo_bounds[X];
          py = (Real) rand() / RAND_MAX * (_hi_bounds[Y] - _lo_bounds[Y]) + _lo_bounds[Y];
          pz = (Real) rand() / RAND_MAX * (_hi_bounds[Z] - _lo_bounds[Z]) + _lo_bounds[Z];
          if (std::fabs(calc_signed_distance(Vector3(px, py, pz))) < 1e-2)
            coords->point.set1Value(idx++, px, py, pz);
        }
      }
    }

    // setup the point set (if any)
    if (idx == 0)
      return;

    // create the point set
    SoPointSet* ps = new SoPointSet;
    ps->numPoints = idx;

    // add the coordinates and the point set to the separator
    SoSeparator* newsep = new SoSeparator;
    newsep->addChild(coords);
    newsep->addChild(ps);
    separator->addChild(newsep);

    SoCoordinate3* lcoords = new SoCoordinate3;
    lcoords->point.set1Value(0, (*_vertices[0])[X], (*_vertices[0])[Y], (*_vertices[0])[Z]);
    lcoords->point.set1Value(1, (*_vertices[1])[X], (*_vertices[1])[Y], (*_vertices[1])[Z]);
    lcoords->point.set1Value(2, (*_vertices[4])[X], (*_vertices[4])[Y], (*_vertices[4])[Z]);
    lcoords->point.set1Value(3, (*_vertices[2])[X], (*_vertices[2])[Y], (*_vertices[2])[Z]);
    lcoords->point.set1Value(4, (*_vertices[0])[X], (*_vertices[0])[Y], (*_vertices[0])[Z]);
    lcoords->point.set1Value(5, (*_vertices[3])[X], (*_vertices[3])[Y], (*_vertices[3])[Z]);
    lcoords->point.set1Value(6, (*_vertices[5])[X], (*_vertices[5])[Y], (*_vertices[5])[Z]);
    lcoords->point.set1Value(7, (*_vertices[1])[X], (*_vertices[1])[Y], (*_vertices[1])[Z]);
    lcoords->point.set1Value(8, (*_vertices[4])[X], (*_vertices[4])[Y], (*_vertices[4])[Z]);
    lcoords->point.set1Value(9, (*_vertices[7])[X], (*_vertices[7])[Y], (*_vertices[7])[Z]);
    lcoords->point.set1Value(10, (*_vertices[5])[X], (*_vertices[5])[Y], (*_vertices[5])[Z]);
    lcoords->point.set1Value(11, (*_vertices[7])[X], (*_vertices[7])[Y], (*_vertices[7])[Z]);
    lcoords->point.set1Value(12, (*_vertices[6])[X], (*_vertices[6])[Y], (*_vertices[6])[Z]);
    lcoords->point.set1Value(13, (*_vertices[3])[X], (*_vertices[3])[Y], (*_vertices[3])[Z]);
    lcoords->point.set1Value(14, (*_vertices[6])[X], (*_vertices[6])[Y], (*_vertices[6])[Z]);
    lcoords->point.set1Value(15, (*_vertices[2])[X], (*_vertices[2])[Y], (*_vertices[2])[Z]);
    lcoords->point.set1Value(16, (*_vertices[4])[X], (*_vertices[4])[Y], (*_vertices[4])[Z]);
    SoLineSet* ls = new SoLineSet;
    ls->numVertices = 17;

    SoSeparator* sep = new SoSeparator;
    SoBaseColor* rndcolor = new SoBaseColor;
    Real r = (Real) rand() / RAND_MAX;
    Real g = (Real) rand() / RAND_MAX;
    Real b = (Real) rand() / RAND_MAX;
    rndcolor->rgb.setValue(SbColor(r,g,b));
    sep->addChild(rndcolor);
    sep->addChild(lcoords);
    sep->addChild(ls);
    separator->addChild(sep);
  }
  #endif
}

/// Performs trilinear interpolation at point p, given a vector of 8 three-dimensional points and their associated functional values
/**
 * The structure of the points must be setup as a box as follows:
 *       6---7
 *     /|  /|
 *    3---5 |
 *    | 2-|-4
 *    |/  |/
 *    0---1
 */
Real ADF::tri_linear_interp(const std::vector<Vector3ConstPtr>& x, const std::vector<Real>& q, const Vector3& p)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // make sure that there are exactly eight values
  assert(x.size() == BOX_VERTICES);
  assert(q.size() == BOX_VERTICES);

  // get the lengths in the x, y, and z directions
  Real xlen = (*x[1])[X] - (*x[0])[X];
  Real ylen = (*x[3])[Y] - (*x[0])[Y];
  Real zlen = (*x[2])[Z] - (*x[0])[Z];

  // get distance from the bounds of the cube in each direction 
  Real dHIx = (*x[1])[X] - p[X];
  Real dLOx = p[X] - (*x[0])[X];
  Real dHIy = (*x[3])[Y] - p[Y];
  Real dLOy = p[Y] - (*x[0])[Y];
  Real dHIz = (*x[2])[Z] - p[Z];
  Real dLOz = p[Z] - (*x[0])[Z];

  // compute the interpolated value
  Real value = 0;
  value += q[0] * dHIx * dHIy * dHIz;
  value += q[1] * dLOx * dHIy * dHIz;
  value += q[2] * dHIx * dHIy * dLOz;
  value += q[3] * dHIx * dLOy * dHIz;
  value += q[4] * dLOx * dHIy * dLOz;
  value += q[5] * dLOx * dLOy * dHIz;
  value += q[6] * dHIx * dLOy * dLOz;
  value += q[7] * dLOx * dLOy * dLOz;

  // scale the value by the lengths of the sides
  value /= xlen * ylen * zlen;

  FILE_LOG(LOG_ADF) << "ADF::tri_linear_interp() entered" << std::endl;
  FILE_LOG(LOG_ADF) << "  evaluating point: " << p << std::endl;
  FILE_LOG(LOG_ADF) << "  vertices and q values: " << std::endl;
  for (unsigned i=0; i< BOX_VERTICES; i++)
    FILE_LOG(LOG_ADF) << "    " << *x[i] << " --> " << q[i] << std::endl;
  FILE_LOG(LOG_ADF) << "ADF::tri_linear_interp() exited" << std::endl;

  return value;
} 

/// Counts the number of cells within this ADF
unsigned ADF::count_cells() const
{
  unsigned count = 0;

  // count all children recursively
  for (unsigned i=0; i< _children.size(); i++)
    count += _children[i]->count_cells();

  // count this cell too
  return count+1;
}

/// Subdivides this ADF cell 
void ADF::subdivide()
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (!_children.empty())
    return;

  // compute the side length
  Vector3 slen = (_hi_bounds - _lo_bounds)/2;

  // allocate the vector for the children 
  _children = std::vector<shared_ptr<ADF> >(OCT_CHILDREN);

  // Subdivision diagram
  //
  // (low y)       (mid y)     (high y)
  // 2---e---4    f---g---i    6---q---7
  // |   |   |    |   |   |    |   |   |
  // |   |   |    |   |   |    |   |   |
  // b---c---d    j---k---m    r---s---t
  // |   |   |    |   |   |    |   |   |
  // |   |   |    |   |   |    |   |   |
  // 0---a---1    n---o---p    3---u---5

  // the boxes
  // lower front left: 0
  // lower front right: 1
  // lower back left: 2
  // lower back right: 4
  // upper front left: 3
  // upper front right: 5
  // upper back left: 6
  // upper back right: 7 

  // each box is specified as follows, where 0 and 7 are the lower and upper bounds
  //      6---7
  //     /|  /|
  //    3---5 |
  //    | 2-|-4
  //    |/  |/
  //    0---1

  // setup midpoints for x, y, and z
  Real midx = _lo_bounds[X] + slen[X];
  Real midy = _lo_bounds[Y] + slen[Y];
  Real midz = _lo_bounds[Z] + slen[Z];

  // setup vectors for sub bounds
  Vector3 sub_lo_bounds, sub_hi_bounds;

  // setup the new vertices (a..u)
  Vector3ConstPtr a(new Vector3(midx, _lo_bounds[Y], _lo_bounds[Z]));
  Vector3ConstPtr b(new Vector3(_lo_bounds[X], _lo_bounds[Y], midz));
  Vector3ConstPtr c(new Vector3(midx, _lo_bounds[Y], midz));
  Vector3ConstPtr d(new Vector3(_hi_bounds[X], _lo_bounds[Y], midz));
  Vector3ConstPtr e(new Vector3(midx, _lo_bounds[Y], _hi_bounds[Z]));
  Vector3ConstPtr f(new Vector3(_lo_bounds[X], midy, _hi_bounds[Z]));
  Vector3ConstPtr g(new Vector3(midx, midy, _hi_bounds[Z]));
  Vector3ConstPtr i(new Vector3(_hi_bounds[X], midy, _hi_bounds[Z]));
  Vector3ConstPtr j(new Vector3(_lo_bounds[X], midy, midz));
  Vector3ConstPtr k(new Vector3(midx, midy, midz));
  Vector3ConstPtr m(new Vector3(_hi_bounds[X], midy, midz));
  Vector3ConstPtr n(new Vector3(_lo_bounds[X], midy, _lo_bounds[Z]));
  Vector3ConstPtr o(new Vector3(midx, midy, _lo_bounds[Z]));
  Vector3ConstPtr p(new Vector3(_hi_bounds[X], midy, _lo_bounds[Z]));
  Vector3ConstPtr q(new Vector3(midx, _hi_bounds[Y], _hi_bounds[Z]));
  Vector3ConstPtr r(new Vector3(_lo_bounds[X], _hi_bounds[Y], midz));
  Vector3ConstPtr s(new Vector3(midx, _hi_bounds[Y], midz));
  Vector3ConstPtr t(new Vector3(_hi_bounds[X], _hi_bounds[Y], midz));
  Vector3ConstPtr u(new Vector3(midx, _hi_bounds[Y], _lo_bounds[Z]));
  
  // box 0
  std::vector<Vector3ConstPtr> box_0_verts;
  box_0_verts.push_back(_vertices[0]);
  box_0_verts.push_back(a);
  box_0_verts.push_back(b);
  box_0_verts.push_back(n);
  box_0_verts.push_back(c);
  box_0_verts.push_back(o);
  box_0_verts.push_back(j);
  box_0_verts.push_back(k);
  _children[0] = shared_ptr<ADF>(new ADF(shared_from_this(), box_0_verts));
  
  // box 1
  std::vector<Vector3ConstPtr> box_1_verts;
  box_1_verts.push_back(a);
  box_1_verts.push_back(_vertices[1]);
  box_1_verts.push_back(c);
  box_1_verts.push_back(o);
  box_1_verts.push_back(d);
  box_1_verts.push_back(p);
  box_1_verts.push_back(k);
  box_1_verts.push_back(m);
  _children[1] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_1_verts));
  
  // box 2
  std::vector<Vector3ConstPtr> box_2_verts;
  box_2_verts.push_back(b);
  box_2_verts.push_back(c);
  box_2_verts.push_back(_vertices[2]);
  box_2_verts.push_back(j);
  box_2_verts.push_back(e);
  box_2_verts.push_back(k);
  box_2_verts.push_back(f);
  box_2_verts.push_back(g);
_children[2] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_2_verts));
  
  // box 4
  std::vector<Vector3ConstPtr> box_4_verts;
  box_4_verts.push_back(c);
  box_4_verts.push_back(d);
  box_4_verts.push_back(e);
  box_4_verts.push_back(k);
  box_4_verts.push_back(_vertices[4]);
  box_4_verts.push_back(m);
  box_4_verts.push_back(g);
  box_4_verts.push_back(i);
  _children[4] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_4_verts));

  // box 3
  std::vector<Vector3ConstPtr> box_3_verts;
  box_3_verts.push_back(n);
  box_3_verts.push_back(o);
  box_3_verts.push_back(j);
  box_3_verts.push_back(_vertices[3]);
  box_3_verts.push_back(k);
  box_3_verts.push_back(u);
  box_3_verts.push_back(r);
  box_3_verts.push_back(s);
  _children[3] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_3_verts));
  
  // box 5
  std::vector<Vector3ConstPtr> box_5_verts;
  box_5_verts.push_back(o);
  box_5_verts.push_back(p);
  box_5_verts.push_back(k);
  box_5_verts.push_back(u);
  box_5_verts.push_back(m);
  box_5_verts.push_back(_vertices[5]);
  box_5_verts.push_back(s);
  box_5_verts.push_back(t);
  _children[5] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_5_verts));
  
  // box 6
  std::vector<Vector3ConstPtr> box_6_verts;
  box_6_verts.push_back(j);
  box_6_verts.push_back(k);
  box_6_verts.push_back(f);
  box_6_verts.push_back(r);
  box_6_verts.push_back(g);
  box_6_verts.push_back(s);
  box_6_verts.push_back(_vertices[6]);
  box_6_verts.push_back(q);
  _children[6] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_6_verts));
  
  // box 7
  std::vector<Vector3ConstPtr> box_7_verts;
  box_7_verts.push_back(k);
  box_7_verts.push_back(m);
  box_7_verts.push_back(g);
  box_7_verts.push_back(s);
  box_7_verts.push_back(i);
  box_7_verts.push_back(t);
  box_7_verts.push_back(q);
  box_7_verts.push_back(_vertices[7]);
  _children[7] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_7_verts));

  // clear distances of this cell to save memory
  _distances.clear();
}


/// Subdivides this ADF cell and computes distance quickly for children using mesh
void ADF::subdivide(Real (*dfn)(const Vector3&, void*), void* data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (!_children.empty())
    return;

  // compute the side length
  Vector3 slen = (_hi_bounds - _lo_bounds)/2;

  // allocate the vector for the children 
  _children = std::vector<shared_ptr<ADF> >(OCT_CHILDREN);

  // Subdivision diagram
  //
  // (low y)       (mid y)     (high y)
  // 2---e---4    f---g---i    6---q---7
  // |   |   |    |   |   |    |   |   |
  // |   |   |    |   |   |    |   |   |
  // b---c---d    j---k---m    r---s---t
  // |   |   |    |   |   |    |   |   |
  // |   |   |    |   |   |    |   |   |
  // 0---a---1    n---o---p    3---u---5

  // the boxes
  // lower front left: 0
  // lower front right: 1
  // lower back left: 2
  // lower back right: 4
  // upper front left: 3
  // upper front right: 5
  // upper back left: 6
  // upper back right: 7 

  // each box is specified as follows, where 0 and 7 are the lower and upper bounds
  //      6---7
  //     /|  /|
  //    3---5 |
  //    | 2-|-4
  //    |/  |/
  //    0---1

  // setup midpoints for x, y, and z
  Real midx = _lo_bounds[X] + slen[X];
  Real midy = _lo_bounds[Y] + slen[Y];
  Real midz = _lo_bounds[Z] + slen[Z];

  // setup vectors for sub bounds
  Vector3 sub_lo_bounds, sub_hi_bounds;

  // setup the new vertices (a..u)
  Vector3ConstPtr a(new Vector3(midx, _lo_bounds[Y], _lo_bounds[Z]));
  Vector3ConstPtr b(new Vector3(_lo_bounds[X], _lo_bounds[Y], midz));
  Vector3ConstPtr c(new Vector3(midx, _lo_bounds[Y], midz));
  Vector3ConstPtr d(new Vector3(_hi_bounds[X], _lo_bounds[Y], midz));
  Vector3ConstPtr e(new Vector3(midx, _lo_bounds[Y], _hi_bounds[Z]));
  Vector3ConstPtr f(new Vector3(_lo_bounds[X], midy, _hi_bounds[Z]));
  Vector3ConstPtr g(new Vector3(midx, midy, _hi_bounds[Z]));
  Vector3ConstPtr i(new Vector3(_hi_bounds[X], midy, _hi_bounds[Z]));
  Vector3ConstPtr j(new Vector3(_lo_bounds[X], midy, midz));
  Vector3ConstPtr k(new Vector3(midx, midy, midz));
  Vector3ConstPtr m(new Vector3(_hi_bounds[X], midy, midz));
  Vector3ConstPtr n(new Vector3(_lo_bounds[X], midy, _lo_bounds[Z]));
  Vector3ConstPtr o(new Vector3(midx, midy, _lo_bounds[Z]));
  Vector3ConstPtr p(new Vector3(_hi_bounds[X], midy, _lo_bounds[Z]));
  Vector3ConstPtr q(new Vector3(midx, _hi_bounds[Y], _hi_bounds[Z]));
  Vector3ConstPtr r(new Vector3(_lo_bounds[X], _hi_bounds[Y], midz));
  Vector3ConstPtr s(new Vector3(midx, _hi_bounds[Y], midz));
  Vector3ConstPtr t(new Vector3(_hi_bounds[X], _hi_bounds[Y], midz));
  Vector3ConstPtr u(new Vector3(midx, _hi_bounds[Y], _lo_bounds[Z]));
  
  // compute necessary distances
  Real a_d = dfn(*a, data);
  Real b_d = dfn(*b, data);
  Real c_d = dfn(*c, data);
  Real d_d = dfn(*d, data);
  Real e_d = dfn(*e, data);
  Real f_d = dfn(*f, data);
  Real g_d = dfn(*g, data);
  Real i_d = dfn(*i, data);
  Real j_d = dfn(*j, data);
  Real k_d = dfn(*k, data);
  Real m_d = dfn(*m, data);
  Real n_d = dfn(*n, data);
  Real o_d = dfn(*o, data);
  Real p_d = dfn(*p, data);
  Real q_d = dfn(*q, data);
  Real r_d = dfn(*r, data);
  Real s_d = dfn(*s, data);
  Real t_d = dfn(*t, data);
  Real u_d = dfn(*u, data);

  // box 0
  std::vector<Vector3ConstPtr> box_0_verts;
  box_0_verts.push_back(_vertices[0]);
  box_0_verts.push_back(a);
  box_0_verts.push_back(b);
  box_0_verts.push_back(n);
  box_0_verts.push_back(c);
  box_0_verts.push_back(o);
  box_0_verts.push_back(j);
  box_0_verts.push_back(k);
  std::vector<Real> box_0_distances;
  box_0_distances.push_back(_distances[0]);
  box_0_distances.push_back(a_d);
  box_0_distances.push_back(b_d);
  box_0_distances.push_back(n_d);
  box_0_distances.push_back(c_d);
  box_0_distances.push_back(o_d);
  box_0_distances.push_back(j_d);
  box_0_distances.push_back(k_d);
  _children[0] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_0_verts));
  _children[0]->_distances = box_0_distances;
  
  // box 1
  std::vector<Vector3ConstPtr> box_1_verts;
  box_1_verts.push_back(a);
  box_1_verts.push_back(_vertices[1]);
  box_1_verts.push_back(c);
  box_1_verts.push_back(o);
  box_1_verts.push_back(d);
  box_1_verts.push_back(p);
  box_1_verts.push_back(k);
  box_1_verts.push_back(m);
  std::vector<Real> box_1_distances;
  box_1_distances.push_back(a_d);
  box_1_distances.push_back(_distances[1]);
  box_1_distances.push_back(c_d);
  box_1_distances.push_back(o_d);
  box_1_distances.push_back(d_d);
  box_1_distances.push_back(p_d);
  box_1_distances.push_back(k_d);
  box_1_distances.push_back(m_d);
  _children[1] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_1_verts));
  _children[1]->_distances = box_1_distances;
  
  // box 2
  std::vector<Vector3ConstPtr> box_2_verts;
  box_2_verts.push_back(b);
  box_2_verts.push_back(c);
  box_2_verts.push_back(_vertices[2]);
  box_2_verts.push_back(j);
  box_2_verts.push_back(e);
  box_2_verts.push_back(k);
  box_2_verts.push_back(f);
  box_2_verts.push_back(g);
  std::vector<Real> box_2_distances;
  box_2_distances.push_back(b_d);
  box_2_distances.push_back(c_d);
  box_2_distances.push_back(_distances[2]);
  box_2_distances.push_back(j_d);
  box_2_distances.push_back(e_d);
  box_2_distances.push_back(k_d);
  box_2_distances.push_back(f_d);
  box_2_distances.push_back(g_d);
_children[2] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_2_verts));
  _children[2]->_distances = box_2_distances;
  
  // box 4
  std::vector<Vector3ConstPtr> box_4_verts;
  box_4_verts.push_back(c);
  box_4_verts.push_back(d);
  box_4_verts.push_back(e);
  box_4_verts.push_back(k);
  box_4_verts.push_back(_vertices[4]);
  box_4_verts.push_back(m);
  box_4_verts.push_back(g);
  box_4_verts.push_back(i);
  std::vector<Real> box_4_distances;
  box_4_distances.push_back(c_d);
  box_4_distances.push_back(d_d);
  box_4_distances.push_back(e_d);
  box_4_distances.push_back(k_d);
  box_4_distances.push_back(_distances[4]);
  box_4_distances.push_back(m_d);
  box_4_distances.push_back(g_d);
  box_4_distances.push_back(i_d);
  _children[4] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_4_verts));
  _children[4]->_distances = box_4_distances;

  // box 3
  std::vector<Vector3ConstPtr> box_3_verts;
  box_3_verts.push_back(n);
  box_3_verts.push_back(o);
  box_3_verts.push_back(j);
  box_3_verts.push_back(_vertices[3]);
  box_3_verts.push_back(k);
  box_3_verts.push_back(u);
  box_3_verts.push_back(r);
  box_3_verts.push_back(s);
  std::vector<Real> box_3_distances;
  box_3_distances.push_back(n_d);
  box_3_distances.push_back(o_d);
  box_3_distances.push_back(j_d);
  box_3_distances.push_back(_distances[3]);
  box_3_distances.push_back(k_d);
  box_3_distances.push_back(u_d);
  box_3_distances.push_back(r_d);
  box_3_distances.push_back(s_d);
  _children[3] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_3_verts));
  _children[3]->_distances = box_3_distances;
  
  // box 5
  std::vector<Vector3ConstPtr> box_5_verts;
  box_5_verts.push_back(o);
  box_5_verts.push_back(p);
  box_5_verts.push_back(k);
  box_5_verts.push_back(u);
  box_5_verts.push_back(m);
  box_5_verts.push_back(_vertices[5]);
  box_5_verts.push_back(s);
  box_5_verts.push_back(t);
  std::vector<Real> box_5_distances;
  box_5_distances.push_back(o_d);
  box_5_distances.push_back(p_d);
  box_5_distances.push_back(k_d);
  box_5_distances.push_back(u_d);
  box_5_distances.push_back(m_d);
  box_5_distances.push_back(_distances[5]);
  box_5_distances.push_back(s_d);
  box_5_distances.push_back(t_d);
  _children[5] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_5_verts));
  _children[5]->_distances = box_5_distances;
  
  // box 6
  std::vector<Vector3ConstPtr> box_6_verts;
  box_6_verts.push_back(j);
  box_6_verts.push_back(k);
  box_6_verts.push_back(f);
  box_6_verts.push_back(r);
  box_6_verts.push_back(g);
  box_6_verts.push_back(s);
  box_6_verts.push_back(_vertices[6]);
  box_6_verts.push_back(q);
  std::vector<Real> box_6_distances;
  box_6_distances.push_back(j_d);
  box_6_distances.push_back(k_d);
  box_6_distances.push_back(f_d);
  box_6_distances.push_back(r_d);
  box_6_distances.push_back(g_d);
  box_6_distances.push_back(s_d);
  box_6_distances.push_back(_distances[6]);
  box_6_distances.push_back(q_d);
  _children[6] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_6_verts));
  _children[6]->_distances = box_6_distances;
  
  // box 7
  std::vector<Vector3ConstPtr> box_7_verts;
  box_7_verts.push_back(k);
  box_7_verts.push_back(m);
  box_7_verts.push_back(g);
  box_7_verts.push_back(s);
  box_7_verts.push_back(i);
  box_7_verts.push_back(t);
  box_7_verts.push_back(q);
  box_7_verts.push_back(_vertices[7]);
  std::vector<Real> box_7_distances;
  box_7_distances.push_back(k_d);
  box_7_distances.push_back(m_d);
  box_7_distances.push_back(g_d);
  box_7_distances.push_back(s_d);
  box_7_distances.push_back(i_d);
  box_7_distances.push_back(t_d);
  box_7_distances.push_back(q_d);
  box_7_distances.push_back(_distances[7]);
  _children[7] = shared_ptr<ADF>(new ADF(shared_from_this(),  box_7_verts));
  _children[7]->_distances = box_7_distances;
}

/// Saves this ADF to a file
void ADF::save_to_file(const std::string& filename) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // compose the map of vertices
  std::map<Vector3ConstPtr, unsigned> vertices;
  std::list<Vector3ConstPtr> verts;
  std::stack<shared_ptr<const ADF> > s;
  s.push(shared_from_this());
  while (!s.empty())
  {
    // get the top cell off of the stack
    shared_ptr<const ADF> cell = s.top();
    s.pop();

    // save the vertices
    for (unsigned i=0; i< BOX_VERTICES; i++)
      if (vertices.find(cell->_vertices[i]) == vertices.end())
      {
        vertices[cell->_vertices[i]] = verts.size();
        verts.push_back(cell->_vertices[i]);
      }

    // add all children to the stack
    for (unsigned i=0; i< cell->_children.size(); i++)
      s.push(cell->_children[i]);
  }

  // open the file
  std::ofstream out(filename.c_str());

  // write the number of vertices
  out << verts.size() << std::endl;

  // write the vertices
  for (std::list<Vector3ConstPtr>::const_iterator i = verts.begin(); i != verts.end(); i++)
    out << (**i)[X] << " " << (**i)[Y] << " " << (**i)[Z] << std::endl;

  // write each cell
  unsigned written = 0;
  s.push(shared_from_this());
  while (!s.empty())
  {
    // get the top cell off of the stack
    shared_ptr<const ADF> cell = s.top();
    s.pop();
    
    // write the vertex indices
    for (unsigned i=0; i< BOX_VERTICES; i++)
    {
      assert(vertices.find(cell->_vertices[i]) != vertices.end());
      out << vertices[cell->_vertices[i]] << std::endl;
    }

    // write the distances
    for (unsigned i=0; i< BOX_VERTICES; i++)
      out << cell->_distances[i] << " " << std::endl;

    // write the bounds
    out << cell->_lo_bounds[X] << " " << cell->_lo_bounds[Y] << " " << cell->_lo_bounds[Z] << std::endl;
    out << cell->_hi_bounds[X] << " " << cell->_hi_bounds[Y] << " " << cell->_hi_bounds[Z] << std::endl;

    // write whether this is an internal node
    out << !cell->_children.empty() << std::endl;

    // add all children of the cell to the stack for processing
    for (unsigned i=0; i< cell->_children.size(); i++)
      s.push(cell->_children[i]);

    // add to the count of written cells
    written++;
  }

  // close the file
  out.close();

  FILE_LOG(LOG_ADF) << written << " cells written" << std::endl;
}

/// Reads a ADF tree from the given file
shared_ptr<ADF> ADF::load_from_file(const std::string& filename)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // open the file
  std::ifstream in(filename.c_str());

  // read in the number of vertices
  unsigned verts;
  in >> verts;

  // read the vertices
  std::map<unsigned, Vector3ConstPtr> vertices;
  for (unsigned i=0; i< verts; i++)
  {
    Vector3 vec;
    in >> vec[X];
    in >> vec[Y];
    in >> vec[Z];
    Vector3ConstPtr vc(new Vector3(vec));
    vertices[i] = vc;
  }

  // create the root cell and add it to the stack
  shared_ptr<ADF> root(new ADF);
  std::stack<shared_ptr<ADF> > s;
  s.push(root);
  unsigned read = 0;

  // keep reading until the stack is empty
  while (!s.empty())
  {
    // get the cell off of the top of the stack
    shared_ptr<ADF> cell = s.top();
    s.pop();

    // read the cell vertices
    cell->_vertices.clear();
    for (unsigned i=0; i< BOX_VERTICES; i++)
    {
      unsigned idx;
      in >> idx;
      cell->_vertices.push_back(vertices[idx]);
    }

    // read the distances
    cell->_distances = std::vector<Real>(BOX_VERTICES);
    for (unsigned i=0; i< BOX_VERTICES; i++)
      in >> cell->_distances[i];

    // read the bounds
    in >> cell->_lo_bounds[X];
    in >> cell->_lo_bounds[Y];
    in >> cell->_lo_bounds[Z];
    in >> cell->_hi_bounds[X];
    in >> cell->_hi_bounds[Y];
    in >> cell->_hi_bounds[Z];

    // read whether this cell is a leaf node or not
    bool internal;
    in >> internal;

    // if this is not a leaf, create and add all children to the stack
    if (internal)
    {
      cell->_children.clear();
      for (unsigned i=0; i< BOX_VERTICES; i++)
      {
        shared_ptr<ADF> child_cell(new ADF);
        child_cell->_parent = cell;
        cell->_children.push_back(child_cell);
        s.push(child_cell);
      }
    }

    // update the count of read cells
    read++;
  }

  // close the file
  in.close();

  FILE_LOG(LOG_ADF) << read << " cells read" << std::endl;

  return root;
}

/// Sets the distances of all vertices of this cell
void ADF::set_distances(Real (*dfn)(const Vector3&, void*), void* data)
{
  _distances = std::vector<Real>(BOX_VERTICES);
  for (unsigned i=0; i< BOX_VERTICES; i++)
    _distances[i] = dfn(*_vertices[i], data);
}

/// Gets the recursion level of this ADF
unsigned ADF::get_recursion_level() const
{
  if (_parent.expired())
    return 0;
  else
  {
    shared_ptr<ADF> parent(_parent);
    return parent->get_recursion_level()+1;
  }
}

/// Prints out the stats on this ADF cell
std::ostream& Moby::operator<<(std::ostream& out, const ADF& adf)
{
  const std::vector<Vector3ConstPtr>& vertices = adf.get_vertices();
  const std::vector<Real>& distances = adf.get_distances();
  Vector3 lo, hi;
  adf.get_bounds(lo, hi);

  out << "cell level: " << adf.get_recursion_level() << std::endl;
  out << "leaf? " << adf.is_leaf() << std::endl;
  out << "bounds: " << lo << " / " << hi << std::endl;
  out << "cell vertices and distances: " << std::endl;
  for (unsigned i=0; i< vertices.size(); i++)
    out << "  " << *vertices[i] << " --> " << distances[i] << std::endl;

  return out;
}

