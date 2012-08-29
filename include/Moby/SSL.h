/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_SSL_H_
#define _MOBY_SSL_H_

#include <stack>
#include <iostream>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Moby/AAngle.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix4.h>

namespace Moby {

/// A sphere-swept line (SSL) that optionally allows building an SSL tree
class SSL : public BV 
{
  public:
    SSL();
    SSL(const SSL& obb) { operator=(obb); }
    SSL(const Vector3& p1, const Vector3& p2, Real radius);
    SSL(const SSL& s, const Vector3& v);
    void operator=(const SSL& s);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const;
    static Real calc_dist(const SSL& o, const Vector3& p);
    static Real calc_dist(const SSL& a, const SSL& b, Vector3& cpa, Vector3& cpb);
    static Real calc_dist(const SSL& a, const SSL& b, const Matrix4& aTb, Vector3& cpa, Vector3& cpb);
    static bool intersects(const SSL& a, const SSL& b);
    static bool intersects(const SSL& a, const SSL& b, const Matrix4& T);
    static bool intersects(const SSL& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q);
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { return SSL::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSL& a, const Vector3& point, Real tol = NEAR_ZERO);
    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) const { return SSL::outside(*this, point, tol); }
    boost::shared_ptr<SSL> get_this() { return boost::dynamic_pointer_cast<SSL>(shared_from_this()); }
    boost::shared_ptr<const SSL> get_this() const { return boost::dynamic_pointer_cast<const SSL>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const;
    unsigned calc_size() const;
    virtual Vector3 get_lower_bounds(const Matrix4& T);
    virtual Vector3 get_upper_bounds(const Matrix4& T);

    /// Calculates (approximate?) volume of the SSL 
    virtual Real calc_volume() const { return (p1 - p2).norm() * M_PI * radius * radius * radius; }

    /// The first point of the line segment 
    Vector3 p1;

    /// The second point of the line segment
    Vector3 p2;

    /// Radius of the spherical addition 
    Real radius;

  private:
    static Real calc_sq_dist(const LineSeg3& seg, const Vector3& q, Vector3& cp);
    static Real calc_sq_dist(const LineSeg3& seg0, const LineSeg3& seg1, Vector3& cp0, Vector3& cp1);
    static Real calc_sq_dist(const Vector3& origin, const Vector3& dir, const LineSeg3& seg, Vector3& cp_line, Vector3& cp_seg, Real& line_param);

}; // end class

// include inline functions
#include "SSL.inl"

} // end namespace

#endif

