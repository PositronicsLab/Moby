/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU General Public 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GENERALIZED_CCD_PLUGIN_H
#define _GENERALIZED_CCD_PLUGIN_H

#include <vector>
#include <Moby/Types.h>
#include <Moby/OBB.h>

namespace Moby {

/// A plugin for the generalized continuous collision detector
class GeneralizedCCDPlugin
{
  public:
    /// Returns the root bounding volume for intersection testing 
    virtual BVPtr get_BVH_root() = 0;

    /// Gets vertices corresponding to the bounding volume
    virtual void get_vertices(BVPtr bv, std::vector<Point3d>& vertices) = 0; 

    /// Determines whether a point is inside/on the geometry
    /**
     * \param gs_BV the bounding volume in which the point resides
     * \param p the query point
     * \param normal the normal to the query point, if the point is inside/on
     *        the geometry
     * \return <b>true</b> if the point is inside or on the geometry, 
     *         <b>false</b> otherwise
     */
    virtual bool point_inside(BVPtr bv, const Point3d& p, Ravelin::Vector3d& normal) = 0;

    /// Determines whether a line segment and the shape intersect
    /**
     * \param bv a bounding volume (to speed intersection testing)
     * \param seg the line segment
     * \param isect the point of intersection, on return (if any)
     * \param normal the normal to the shape at the point of intersection, 
     *          on return (if any)
     * \return the parameter,t, of the intersection seg.first + seg.second*t
     *         or -1 if no intersection
     */
    virtual double intersect_seg(BVPtr bv, const LineSeg3& seg, Point3d& isect, Ravelin::Vector3d& normal) = 0;

    /// Gets mesh data for the geometry with the specified bounding volume
    /**
     * \param bv the bounding data from which the corresponding mesh data will
     *           be taken
     * \return a pair consisting of an IndexedTriArray and a list of triangle
     *         indices encapsulated by bv
     */
    virtual std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_mesh(BVPtr bv) = 0;
}; // end class

} // end namespace

#endif

