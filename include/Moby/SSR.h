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
#include <Ravelin/LinAlgd.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Point3d.h>
#include <Ravelin/AAngled.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/FastThreadable.h>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>

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
    SSR(const Ravelin::Point3d& center, const Ravelin::Matrix3d& R, const Ravelin::Vector2d& l, double radius);
    SSR(const SSR& s, const Ravelin::Vector3d& v);
    void operator=(const SSR& s);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, double dt, const Ravelin::Vector3d& lv, const Ravelin::Vector3d& av) const;
    static double calc_dist(const SSR& a, const Ravelin::Point3d& p);
    static double calc_dist(const SSR& a, const LineSeg3& s);
    static double calc_dist(const SSR& a, const SSR& b, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static double calc_dist(const SSR& a, const SSR& b, const std::pair<Ravelin::Quatd, Ravelin::Origin3d>& aTb, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static bool intersects(const SSR& a, const SSR& b);
    static bool intersects(const SSR& a, const SSR& b, const std::pair<Ravelin::Quatd, Ravelin::Origin3d>& T);
    static bool intersects(const SSR& a, const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q);
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q) const { return SSR::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSR& a, const Ravelin::Point3d& point, double tol = NEAR_ZERO);
    virtual bool outside(const Ravelin::Point3d& point, double tol = NEAR_ZERO) const { return SSR::outside(*this, point, tol); }
    boost::shared_ptr<SSR> get_this() { return boost::dynamic_pointer_cast<SSR>(shared_from_this()); }
    boost::shared_ptr<const SSR> get_this() const { return boost::dynamic_pointer_cast<const SSR>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    unsigned calc_size() const;
    virtual Ravelin::Point3d get_lower_bounds(const Ravelin::Pose3d& T);
    virtual Ravelin::Point3d get_upper_bounds(const Ravelin::Pose3d& T);
    void get_rect_verts(Ravelin::Point3d rect_verts[4]) const;

    template <class ForwardIterator>
    void expand_to_fit(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    SSR(ForwardIterator begin, ForwardIterator end);

    /// Calculates (approximate?) volume of the SSR 
    virtual double calc_volume() const { return std::sqrt(l[0]*l[0] + l[1]*l[1]) + 2*radius; }

    /// Center of the volume
    Ravelin::Point3d center;

    /// Lengths of the rectangle sides 
    Ravelin::Vector2d l;

    /// Orientation of this SSR
    Ravelin::Matrix3d R;

    /// Radius of the spherical addition 
    double radius;

  private:
    static double calc_sq_dist(const Ravelin::Point3d& p, const Ravelin::Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Ravelin::Point3d& cp_rect);
    static double calc_sq_dist(const Ravelin::Point3d& origin, const Ravelin::Vector3d& dir, const LineSeg3& seg, Ravelin::Point3d& cp_line, Ravelin::Point3d& cp_seg, double& line_param);
    static double calc_sq_dist(const Ravelin::Point3d& origin, const Ravelin::Vector3d& dir, const Ravelin::Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Ravelin::Point3d& cp_line, Ravelin::Point3d& cp_rect, double& line_param);
    static double calc_sq_dist(const LineSeg3& seg, const Ravelin::Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Ravelin::Point3d& cp_seg, Ravelin::Point3d& cp_rect);
    static double calc_sq_dist(const Ravelin::Point3d& a_center, const Ravelin::Vector3d& aaxis1, const Ravelin::Vector3d& aaxis2, const Ravelin::Vector2d& alengths, const Ravelin::Point3d& b_center, const Ravelin::Vector3d& baxis1, const Ravelin::Vector3d& baxis2, const Ravelin::Vector2d& blengths, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);

    template <class ForwardIterator>
    void calc_lengths_and_radius(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void align(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& d1, Ravelin::Vector3d& d2);
}; // end class

// include inline functions
#include "SSR.inl"

} // end namespace

#endif

