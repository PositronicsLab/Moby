/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_SSR_H_
#define _MOBY_SSR_H_

#include <stack>
#include <iostream>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Moby/AAngle.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/FastThreadable.h>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix4.h>

namespace Moby {

/// A sphere-swept rectangle (SSR) that optionally allows building an SSR tree
/**
 * \note the SSR is generally constructed such that its frame is aligned with
 *       that of the underlying rigid body (or collision geometry).  
 *       Therefore, the center of the SSR is computed relative to the 
 *       center-of-mass of the body (or the center of the geometry).  The
 *       orientation of the SSR will always be identical to the orientation of
 *       the rigid body (or collision geometry).
 */
class SSR : public BV 
{
  public:
    SSR();
    SSR(const SSR& obb) { operator=(obb); }
    SSR(const Vector3& center, const Matrix3& R, const Vector2& l, Real radius);
    SSR(const SSR& s, const Vector3& v);
    void operator=(const SSR& s);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const;
    static Real calc_dist(const SSR& a, const Vector3& p);
    static Real calc_dist(const SSR& a, const LineSeg3& s);
    static Real calc_dist(const SSR& a, const SSR& b, Vector3& cpa, Vector3& cpb);
    static Real calc_dist(const SSR& a, const SSR& b, const Matrix4& aTb, Vector3& cpa, Vector3& cpb);
    static bool intersects(const SSR& a, const SSR& b);
    static bool intersects(const SSR& a, const SSR& b, const Matrix4& T);
    static bool intersects(const SSR& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q);
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { return SSR::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSR& a, const Vector3& point, Real tol = NEAR_ZERO);
    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) const { return SSR::outside(*this, point, tol); }
    boost::shared_ptr<SSR> get_this() { return boost::dynamic_pointer_cast<SSR>(shared_from_this()); }
    boost::shared_ptr<const SSR> get_this() const { return boost::dynamic_pointer_cast<const SSR>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const;
    unsigned calc_size() const;
    virtual Vector3 get_lower_bounds(const Matrix4& T);
    virtual Vector3 get_upper_bounds(const Matrix4& T);
    void get_rect_verts(Vector3 rect_verts[4]) const;

    template <class ForwardIterator>
    void expand_to_fit(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    SSR(ForwardIterator begin, ForwardIterator end);

    /// Calculates (approximate?) volume of the SSR 
    virtual Real calc_volume() const { return std::sqrt(l[0]*l[0] + l[1]*l[1]) + 2*radius; }

    /// Center of the volume
    Vector3 center;

    /// Lengths of the rectangle sides 
    Vector2 l;

    /// Orientation of this SSR
    Matrix3 R;

    /// Radius of the spherical addition 
    Real radius;

  private:
    static Real calc_sq_dist(const Vector3& p, const Vector3& rect_center, const Vector3& axis1, const Vector3& axis2, const Vector2& lengths, Vector3& cp_rect);
    static Real calc_sq_dist(const Vector3& origin, const Vector3& dir, const LineSeg3& seg, Vector3& cp_line, Vector3& cp_seg, Real& line_param);
    static Real calc_sq_dist(const Vector3& origin, const Vector3& dir, const Vector3& rect_center, const Vector3& axis1, const Vector3& axis2, const Vector2& lengths, Vector3& cp_line, Vector3& cp_rect, Real& line_param);
    static Real calc_sq_dist(const LineSeg3& seg, const Vector3& rect_center, const Vector3& axis1, const Vector3& axis2, const Vector2& lengths, Vector3& cp_seg, Vector3& cp_rect);
    static Real calc_sq_dist(const Vector3& a_center, const Vector3& aaxis1, const Vector3& aaxis2, const Vector2& alengths, const Vector3& b_center, const Vector3& baxis1, const Vector3& baxis2, const Vector2& blengths, Vector3& cpa, Vector3& cpb);

    template <class ForwardIterator>
    void calc_lengths_and_radius(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void align(ForwardIterator begin, ForwardIterator end, const Vector3& d1, Vector3& d2);
}; // end class

// include inline functions
#include "SSR.inl"

} // end namespace

#endif

