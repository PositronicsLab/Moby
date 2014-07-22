/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_DUMMY_BV_H_
#define _MOBY_DUMMY_BV_H_

#include <Moby/BV.h>

namespace Moby {

/// An abstract bounding volume 
/**
 * \note the BV is generally constructed such that its frame is aligned with
 *       that of the underlying rigid body (or collision geometry).  
 *       Therefore, the center of the BV is computed relative to the 
 *       center-of-mass of the body (or the center of the geometry).  The
 *       orientation of the BV will always be identical to the orientation of
 *       the rigid body (or collision geometry).
 */
class DummyBV : public BV 
{
  public:
    virtual ~DummyBV() {}

    /// Nothing will be output
    std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const { return out; }

    /// Gets the associated pose for this bounding volume
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const { return GLOBAL; }

    /// Determines whether a point is outside the bounding volume
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const { return false; }

    /// Determines whether a line segment intersects the bounding volume
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const { q = seg.first*tmin + seg.second*((double) 1.0 - tmin); return true; }

    /// Virtual function that calculates a velocity-expanded BV
    /**
     * \param g the geometry that this bounding volume represents
     * \param dt the time step
     * \param v the velocity
     * \return the velocity-expanded bounding volume
     */ 
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const { return boost::const_pointer_cast<BV>(get_this()); }

    /// Volume will be zero
    virtual double calc_volume() const { return 0.0; }

    /// Gets the lower bounds
    virtual Point3d get_lower_bounds() const { double INF = std::numeric_limits<double>::max(); return Point3d(-INF, -INF, -INF); }

    /// Gets the upper bounds
    virtual Point3d get_upper_bounds() const { double INF = std::numeric_limits<double>::max(); return Point3d(INF, INF, INF); }

    /// Transforms the dummy BV (does nothing)
    virtual void transform(const Ravelin::Transform3d& T, BV* result) const {}
}; // end class

} // end namespace

#endif

