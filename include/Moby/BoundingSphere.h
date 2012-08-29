/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_BSPHERE_H_
#define _MOBY_BSPHERE_H_

#include <Moby/BV.h>
#include <Moby/DummyBV.h>

namespace Moby {

/// A sphere used for bounding geometry
class BoundingSphere : public BV 
{
  public:
    BoundingSphere();
    BoundingSphere(const BoundingSphere& bsphere) { operator=(bsphere); }
    BoundingSphere(const Vector3& center, Real radius);
    BoundingSphere& operator=(const BoundingSphere& bsphere);
    virtual std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const;
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const;
    static Real calc_dist(const BoundingSphere& s1, const BoundingSphere& s2);
    static bool intersects(const BoundingSphere& a, const BoundingSphere& b);
    static bool intersects(const BoundingSphere& a, const BoundingSphere& b, const Matrix4& T);
    static bool intersects(const BoundingSphere& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q);
    static bool outside(const BoundingSphere& a, const Vector3& point, Real tol = NEAR_ZERO);
    boost::shared_ptr<const BoundingSphere> get_this() const { return boost::dynamic_pointer_cast<const BoundingSphere>(shared_from_this()); }
    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) const { return BoundingSphere::outside(*this, point, tol); }
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { return BoundingSphere::intersects(*this, seg, tmin, tmax, q); }

    template <class ForwardIterator>
    BoundingSphere(ForwardIterator begin, ForwardIterator end);

    /// Gets the lower bounds on the bounding sphere
    virtual Vector3 get_lower_bounds(const Matrix4& T) { return T.get_translation()+Vector3(center[0]-radius, center[1]-radius, center[2] - radius); }

    /// Gets the upper bounds on the bounding sphere
    virtual Vector3 get_upper_bounds(const Matrix4& T) { return T.get_translation()+Vector3(center[0]+radius, center[1]+radius, center[2] + radius); }

    /// Center of the bounding box
    Vector3 center;

    /// The radius of the bounding sphere (we use a float b/c accuracy here not so important)
    Real radius;

    /// Calculates the volume of this bounding volume
    virtual Real calc_volume() const { return M_PI * radius * radius * radius; }
}; // end class

// inline functions
#include "BoundingSphere.inl"

} // end namespace

#endif

