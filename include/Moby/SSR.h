/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
    SSR(const Point3d& center, const Ravelin::Matrix3d& R, const Ravelin::Vector2d& l, double radius);
    SSR(const SSR& s, const Ravelin::Vector3d& v);
    SSR& operator=(const SSR& s);
    virtual BVPtr calc_swept_BV(CollisionGeometryPtr g, const Ravelin::SVelocityd& v) const;
    static double calc_dist(const SSR& a, const Point3d& p);
    static double calc_dist(const SSR& a, const LineSeg3& s);
    static double calc_dist(const SSR& a, const SSR& b, Point3d& cpa, Point3d& cpb);
    static double calc_dist(const SSR& a, const SSR& b, const Ravelin::Transform3d& aTb, Point3d& cpa, Point3d& cpb);
    static bool intersects(const SSR& a, const SSR& b);
    static bool intersects(const SSR& a, const SSR& b, const Ravelin::Transform3d& T);
    static bool intersects(const SSR& a, const LineSeg3& seg, double& tmin, double tmax, Point3d& q);
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Point3d& q) const { return SSR::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSR& a, const Point3d& point, double tol = NEAR_ZERO);
    virtual bool outside(const Point3d& point, double tol = NEAR_ZERO) const { return SSR::outside(*this, point, tol); }
    boost::shared_ptr<SSR> get_this() { return boost::dynamic_pointer_cast<SSR>(shared_from_this()); }
    boost::shared_ptr<const SSR> get_this() const { return boost::dynamic_pointer_cast<const SSR>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    unsigned calc_size() const;
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_relative_pose() const { return center.pose; }
    virtual void transform(const Ravelin::Transform3d& T, BV* result) const;
    virtual Point3d get_lower_bounds() const;
    virtual Point3d get_upper_bounds() const;
    void get_rect_verts(Point3d rect_verts[4]) const;

    template <class ForwardIterator>
    void expand_to_fit(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    SSR(ForwardIterator begin, ForwardIterator end);

    /// Calculates (approximate?) volume of the SSR 
    virtual double calc_volume() const { return std::sqrt(l[0]*l[0] + l[1]*l[1]) + 2*radius; }

    /// Center of the volume
    Point3d center;

    /// Lengths of the rectangle sides 
    Ravelin::Vector2d l;

    /// Orientation of this SSR
    Ravelin::Matrix3d R;

    /// Radius of the spherical addition 
    double radius;

  private:
    static double calc_sq_dist(const Point3d& p, const Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Point3d& cp_rect);
    static double calc_sq_dist(const Point3d& origin, const Ravelin::Vector3d& dir, const LineSeg3& seg, Point3d& cp_line, Point3d& cp_seg, double& line_param);
    static double calc_sq_dist(const Point3d& origin, const Ravelin::Vector3d& dir, const Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Point3d& cp_line, Point3d& cp_rect, double& line_param);
    static double calc_sq_dist(const LineSeg3& seg, const Point3d& rect_center, const Ravelin::Vector3d& axis1, const Ravelin::Vector3d& axis2, const Ravelin::Vector2d& lengths, Point3d& cp_seg, Point3d& cp_rect);
    static double calc_sq_dist(const Point3d& a_center, const Ravelin::Vector3d& aaxis1, const Ravelin::Vector3d& aaxis2, const Ravelin::Vector2d& alengths, const Point3d& b_center, const Ravelin::Vector3d& baxis1, const Ravelin::Vector3d& baxis2, const Ravelin::Vector2d& blengths, Point3d& cpa, Point3d& cpb);

    template <class ForwardIterator>
    void calc_lengths_and_radius(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void align(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& d1, Ravelin::Vector3d& d2);
}; // end class

// include inline functions
#include "SSR.inl"

} // end namespace

#endif

