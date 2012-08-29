/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_AABB_H_
#define _MOBY_AABB_H_

#include <Moby/OBB.h>
#include <Moby/BV.h>

namespace Moby {

/// An axis-aligned bounding box 
class AABB : public BV
{
  public:
    AABB() {}

    template <class InputIterator>
    AABB(InputIterator begin, InputIterator end);

    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) const { return AABB::outside(*this, point, tol); }
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { return AABB::intersects(*this, seg, tmin, tmax, q); }
    virtual std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const;
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const;
    virtual Real calc_volume() const;
    virtual Vector3 get_lower_bounds(const Matrix4& T);
    virtual Vector3 get_upper_bounds(const Matrix4& T);
    OBB get_OBB() const;
    static bool intersects(const AABB& a, const AABB& b);
    static bool intersects(const AABB& a, const AABB& b, const Matrix4& aTb);
    static bool outside(const AABB& a, const Vector3& point, Real tol = NEAR_ZERO);
    static bool intersects(const AABB& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q);
    static void get_closest_point(const AABB& a, const Vector3& p, Vector3& closest);
    static Real get_farthest_point(const AABB& a, const Vector3& p, Vector3& farthest);

    /// The lower corner of the AABB
    Vector3 minp;

    /// The upper corner of the AABB;
    Vector3 maxp;
}; // end class

#include "AABB.inl"

} // end namespace

#endif

