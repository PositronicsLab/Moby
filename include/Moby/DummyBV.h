/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
    std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const { return out; }

    /// Determines whether a point is outside the bounding volume
    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) { return false; }

    /// Determines whether a line segment intersects the bounding volume
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { q = seg.first*tmin + seg.second*((Real) 1.0 - tmin); return true; }

    /// Virtual function that calculates a velocity-expanded BV
    /**
     * \param g the geometry that this bounding volume represents
     * \param dt the time step
     * \param lv the linear velocity
     * \param av the angular velocity
     * \return the velocity-expanded bounding volume
     */ 
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const { return boost::const_pointer_cast<BV>(get_this()); }

    /// Volume will be zero
    virtual Real calc_volume() const { return 0.0; }

    /// Gets the lower bounds
    virtual Vector3 get_lower_bounds(const Matrix4& T) { Real INF = std::numeric_limits<Real>::max(); return Vector3(-INF, -INF, -INF); }

    /// Gets the upper bounds
    virtual Vector3 get_upper_bounds(const Matrix4& T) { Real INF = std::numeric_limits<Real>::max(); return Vector3(INF, INF, INF); }

}; // end class

} // end namespace

#endif

