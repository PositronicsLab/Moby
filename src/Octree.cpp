/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <queue>
#include <limits>
#include <Moby/Octree.h>

using namespace Moby;

/// Constructs an octree with zero minimum resolution
/**
 * The bounds of the octree are set to -/+ infinity.
 */
Octree::Octree()
{
	_bounds_lo = Vector3(-std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max());
	_bounds_hi = Vector3(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max());
	_num_points = 0;
	_minres = 0;
}

/// Constructs an octree with the specified minimum resolution
/**
 * The bounds of the octree are set to -/+ infinity.
 */
Octree::Octree(Real minres)
{
	_bounds_lo = Vector3(-std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max());
	_bounds_hi = Vector3(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max());
	_num_points = 0;
	_minres = minres;
}

/// Constructs an octree with the specified minimum resolution and bounds
Octree::Octree(Real minres, const Vector3& lo_bounds, const Vector3& hi_bounds)
{
	_bounds_lo = lo_bounds; 
	_bounds_hi = hi_bounds; 
	_num_points = 0;
	_minres = minres;
}

/// Constructs an octree with the specified parent, minimum resolution, and bounds
Octree::Octree(OctreePtr parent, Real minres, const Vector3& lo_bounds, const Vector3& hi_bounds)
{
	_bounds_lo = lo_bounds; 
	_bounds_hi = hi_bounds; 
	_num_points = 0;
	_minres = minres;
	_parent = parent;
}

/// Inserts a point into the octree
void Octree::insert(const Vector3& point)
{
	const unsigned X = 0;

	#ifndef NDEBUG
	// verify that the point is not outside of the bounds of this octree node
	const unsigned THREE_D = 3;
	for (unsigned i=0; i< THREE_D; i++)
		assert(point[i] >= _bounds_lo[i] && point[i] <= _bounds_hi[i]);
	#endif

	// increment the number of points in this node
	_num_points++;

	// don't subdivide if we are below minimum resolution
	if ((_bounds_hi[X] - _bounds_lo[X])/2 < _minres)
		return;

	// subdivide, if possible
	if (subdivide())
	{
		unsigned subvolume_idx = get_sub_volume_idx(point);
		_children[subvolume_idx]->insert(point);
	}
}

/// Returns the index of the subvolume being occupied
/**
 * \return an index in the range [0,7]; assertion fails if the point is outside
 *         of the volume
 */
unsigned Octree::get_sub_volume_idx(const Vector3& point) const
{
	const unsigned X = 0, Y = 1, Z = 2;

	Vector3 half = 0.5*(_bounds_lo + _bounds_hi);

	// Subdivision order
	// 	   1---5
  //     /|  /|
  //    0---4 |
  //    | 3-|-7
  //    |/  |/
  //    2---6
	
	if (point[X] < half[X])
	{
		if (point[Y] < half[Y])
		{
			if (point[Z] < half[Z])
				return 2;
			else
				return 3;
		}
		else
		{
			if (point[Z] < half[Z])
				return 0;
			else
				return 1;
		}
	}
	else
	{
		if (point[Y] < half[Y])
		{
			if (point[Z] < half[Z])
				return 6;
			else
				return 7;
		}
		else
		{
			if (point[Z] < half[Z])
				return 4;
			else
				return 5;
		}
	}

	assert(false);
	return std::numeric_limits<unsigned>::max();
}

/// Determines (quickly) whether any part of a rectangular region of space is occupied
bool Octree::is_occupied(const Vector3& lo_bounds, const Vector3& hi_bounds) const
{
	// only should be called on top-level bounding box
	assert(_parent.expired());

	// verify that regions are correct
	#ifndef NDEBUG
	for (unsigned i=0; i< 3; i++)
		assert(lo_bounds[i] <= hi_bounds[i]);
	#endif

	// if no points, return false
	if (_num_points == 0)
		return false;

	// setup a queue, with this as the first element
	std::queue<OctreePtr> q;
	q.push(boost::const_pointer_cast<Octree>(shared_from_this()));

	while (!q.empty())
	{
		// get the node at the front of the queue
		OctreePtr node = q.front();
		q.pop();

		// process the children, if any, for this octree
		const std::vector<OctreePtr>& children = node->get_children();
		for (std::vector<OctreePtr>::const_iterator i = children.begin(); i != children.end(); i++)
		{
			// if there are no hits for this child, continue
			if ((*i)->get_num_points() == 0)
				continue;

			// get the bounding region for this child
			Vector3 child_lo_bounds, child_hi_bounds;
			(*i)->get_bounds(child_lo_bounds, child_hi_bounds);
			
			// determine how much of the bounding box is within the bounding region
			BBQueryType query = within_bounding_box(std::make_pair(child_lo_bounds, child_hi_bounds), std::make_pair(lo_bounds, hi_bounds));

			// if the cell is entirely within the bounding box, return now
			if (query == eFullyInside)
				return true;
 
			// if the cell lies partially within the bounding box, add it to the queue
			if (query == ePartiallyInside)
				q.push(*i);
		}
	}

	// if still here, region is clear
	return false;
}

/// Gets the bounds for this Octree cell
void Octree::get_bounds(Vector3& lo, Vector3& hi) const
{
	lo = _bounds_lo;
	hi = _bounds_hi;
}

/// Clears the point from the Octree
/**
 * \return <b>true</b> if the point was found in the octree (and deleted),
 *         and <b>false</b> otherwise
 */
bool Octree::clear_cell(const Vector3& point)
{
	const unsigned THREE_D = 3;

	// for parent, make sure cell is occupied
	if (_parent.expired() && !is_cell_occupied(point))
		return false;

	// check whether point is out of bounds of this node
	for (unsigned i=0; i< THREE_D; i++)
		if (_bounds_lo[i] > point[i] || _bounds_hi[i] < point[i])
			return false;

	// if there are no points in this cell, return false
	if (_num_points == 0)
		return false;

	// recursively remove point from occupying cell
	_num_points--;

	// if there are no more points left in this branch of the octree, remove
	// all child nodes
	if (_num_points == 0)
		_children.clear();
	else if (!_children.empty())
	{
		unsigned subvolume_idx = get_sub_volume_idx(point);
		_children[subvolume_idx]->clear_cell(point);
	}

	return true;
}

/// Determines whether a bounding box is within another bounding box
/**
 * \return <b>eFullyInside</b> if query is fully inside box, 
 * <b>ePartiallyInside</b> if query intersects box, and <b>eOutside</b> if
 * the two bounding boxes are separate.
 */
Octree::BBQueryType Octree::within_bounding_box(const std::pair<Vector3, Vector3>& query, const std::pair<Vector3, Vector3>& reference)
{
	const unsigned X = 0, Y = 1, Z = 2;

	// setup all points of the query
	std::list<Vector3> points;
	points.push_back(Vector3(query.first[X], query.first[Y], query.first[Z]));
	points.push_back(Vector3(query.first[X], query.first[Y], query.second[Z]));
	points.push_back(Vector3(query.first[X], query.second[Y], query.first[Z]));
	points.push_back(Vector3(query.first[X], query.first[Y], query.second[Z]));
	points.push_back(Vector3(query.second[X], query.first[Y], query.first[Z]));
	points.push_back(Vector3(query.second[X], query.first[Y], query.second[Z]));
	points.push_back(Vector3(query.second[X], query.second[Y], query.first[Z]));
	points.push_back(Vector3(query.second[X], query.first[Y], query.second[Z]));
	
	// set indicators to false initially
	bool outside_one = false, inside_one = false;
	for (std::list<Vector3>::iterator i = points.begin(); i != points.end(); i++)
	{
		// check the i'th point
		if (within_bounding_box(*i, reference))
			inside_one = true;
		else
			outside_one = true;

		// check whether we can exit prematurely
		if (inside_one && outside_one)
			return ePartiallyInside;
	}

	// return the appropriate value
	return (inside_one) ? eFullyInside : eOutside;
}

/// Determines whether a query point is within the bounding box
bool Octree::within_bounding_box(const Vector3& query, const std::pair<Vector3, Vector3>& reference)
{
	const unsigned THREE_D = 3;

	for (unsigned i=0; i< THREE_D; i++)
		if (query[i] < reference.first[i] || query[i] > reference.second[i])
			return false;

	return true;
}

/// Returns NULL if empty, Octree of leaf if occupied
OctreePtr Octree::is_cell_occupied(const Vector3& point) const
{
	const unsigned THREE_D = 3;

	// first test whether point is in box
	if (_parent.expired())
		for (unsigned i=0; i< THREE_D; i++)
			if (point[i] < _bounds_lo[i] || point[i] > _bounds_hi[i])
				return OctreePtr();

	// if there are no points in this leaf, return NULL
	if (_num_points == 0)
		return OctreePtr();

	// if there are no children, return this
	if (_children.empty())
		return boost::const_pointer_cast<Octree>(shared_from_this());
	else
	{
		unsigned subvolume = get_sub_volume_idx(point);
		return _children[subvolume]->is_cell_occupied(point);
	}

	return OctreePtr();
}

/// (Re)sets the bounds for this octree
void Octree::set_bounds(const Vector3& lo_bound, const Vector3& hi_bound)
{
	_bounds_lo = lo_bound;
	_bounds_hi = hi_bound;
	reset();
}

//// Resets the octree
void Octree::reset()
{
	// clear the number of points
	_num_points = 0;

	// recursively destroy all child nodes
	for (unsigned i=0; i< _children.size(); i++)
		_children[i]->reset();
	_children.clear();
}

/// Subdivides this octree (if possible)
/**
 * \return <b>true</b> if successfully subdivides, <b>false</b> otherwise
 */
bool Octree::subdivide()
{
	const unsigned X = 0, Y = 1, Z = 2;

	if (!_children.empty())
		return true;

	// compute the side length
	Vector3 sidelen = (_bounds_hi - _bounds_lo)/2;

	// check that the minimum resolution is not exceeded (we're assuming cell
	// is a cube)
	if (sidelen[X] < _minres)
		return false;

	// allocate the vector for octree children
	_children = std::vector<OctreePtr>(OCT_CHILDREN); 

	// Subdivision order
  // 	   1---5
  //     /|  /|
  //    0---4 |
  //    | 3-|-7
  //    |/  |/
  //    2---6

	// setup vectors for sub bounds
	Vector3 sub_bounds_lo, sub_bounds_hi;

	// box 0
	sub_bounds_lo[X] = _bounds_lo[X];
	sub_bounds_lo[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z];
	sub_bounds_hi[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_hi[Y] = _bounds_hi[Y];
	sub_bounds_hi[Z] = _bounds_lo[Z] + sidelen[Z];
	_children[0] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 1
	sub_bounds_lo[X] = _bounds_lo[X];
	sub_bounds_lo[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z] + sidelen[Z];
	sub_bounds_hi[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_hi[Y] = _bounds_hi[Y];
	sub_bounds_hi[Z] = _bounds_hi[Z];
	_children[1] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 2
	sub_bounds_lo[X] = _bounds_lo[X];
	sub_bounds_lo[Y] = _bounds_lo[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z];
	sub_bounds_hi[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_hi[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_hi[Z] = _bounds_lo[Z] + sidelen[Z];
	_children[2] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 3
	sub_bounds_lo[X] = _bounds_lo[X];
	sub_bounds_lo[Y] = _bounds_lo[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z] + sidelen[Z];
	sub_bounds_hi[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_hi[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_hi[Z] = _bounds_hi[Z];
	_children[3] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 4
	sub_bounds_lo[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_lo[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z];
	sub_bounds_hi[X] = _bounds_hi[X];
	sub_bounds_hi[Y] = _bounds_hi[Y];
	sub_bounds_hi[Z] = _bounds_lo[Z] + sidelen[Z];
	_children[4] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 5
	sub_bounds_lo[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_lo[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z] + sidelen[Z];
	sub_bounds_hi[X] = _bounds_hi[X];
	sub_bounds_hi[Y] = _bounds_hi[Y];
	sub_bounds_hi[Z] = _bounds_hi[Z];
	_children[5] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 6
	sub_bounds_lo[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_lo[Y] = _bounds_lo[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z];
	sub_bounds_hi[X] = _bounds_hi[X];
	sub_bounds_hi[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_hi[Z] = _bounds_lo[Z] + sidelen[Z];
	_children[6] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));
	
	// box 7
	sub_bounds_lo[X] = _bounds_lo[X] + sidelen[X];
	sub_bounds_lo[Y] = _bounds_lo[Y];
	sub_bounds_lo[Z] = _bounds_lo[Z] + sidelen[Z];
	sub_bounds_hi[X] = _bounds_hi[X];
	sub_bounds_hi[Y] = _bounds_lo[Y] + sidelen[Y];
	sub_bounds_hi[Z] = _bounds_hi[Z];
	_children[7] = OctreePtr(new Octree(shared_from_this(), _minres, sub_bounds_lo, sub_bounds_hi));

	return true;
}

