/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_OCTREE_H_
#define _MOBY_OCTREE_H_

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <list>
#include <Moby/Types.h>
#include <Moby/Vector3.h>

namespace Moby {

/// A generic octree
class Octree : public boost::enable_shared_from_this<Octree>
{
	public:
		Octree();
		Octree(Real minres);
		Octree(Real minres, const Vector3& lo_bounds, const Vector3& hi_bounds);
		Octree(boost::shared_ptr<Octree> parent, Real minres, const Vector3& lo_bounds, const Vector3& hi_bounds);
		void insert(const Vector3& point);
		bool clear_cell(const Vector3& point);
		OctreePtr is_cell_occupied(const Vector3& point) const;
		bool is_occupied(const Vector3& lo_bounds, const Vector3& hi_bounds) const;
		void set_bounds(const Vector3& lo_bounds, const Vector3& hi_bounds);
		void get_bounds(Vector3& lo_bounds, Vector3& hi_bounds) const;
		void reset();

		/// Gets the vector of children of this octree node
		const std::vector<OctreePtr>& get_children() const { return _children; }

		/// Sets the parent of this octree
		void set_parent(OctreePtr parent) { _parent = parent; }

		/// Gets the parent of this octree
		OctreePtr get_parent() const { return (_parent.expired()) ? OctreePtr() : OctreePtr(_parent); }

		/// Gets the number of points in this octree node
		unsigned get_num_points() const { return _num_points; }

		/// Gets the minimum resolution for this octree
		Real get_min_res() const { return _minres; }

		/// Sets the minimum resolution for this octree
		void set_min_res(Real minres) { _minres = minres; }

		/// Determines whether this octree node is a leaf
		bool is_leaf() const { return _children.empty(); }

	private:
		enum BBQueryType { eOutside, ePartiallyInside, eFullyInside };
		static const unsigned OCT_CHILDREN = 8;
		std::vector<OctreePtr> _children;
		boost::weak_ptr<Octree> _parent;
		Vector3 _bounds_lo, _bounds_hi;
		unsigned _num_points;
		Real _minres;

		static BBQueryType within_bounding_box(const std::pair<Vector3, Vector3>& query, const std::pair<Vector3, Vector3>& reference);
		static bool within_bounding_box(const Vector3& query, const std::pair<Vector3, Vector3>& reference);
		unsigned get_sub_volume_idx(const Vector3& point) const;
		bool subdivide();
}; // end class

} // end namespace

#endif

