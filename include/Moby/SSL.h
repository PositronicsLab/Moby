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
#include <Ravelin/Vector3d.h>
#include <Ravelin/Point3d.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/AAngled.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>

namespace Moby {

/// A sphere-swept line (SSL) that optionally allows building an SSL tree
class SSL : public BV 
{
  public:
    SSL();
    SSL(const SSL& obb) { operator=(obb); }
    SSL(const Ravelin::Point3d& p1, const Ravelin::Point3d& p2, double radius);
    SSL(const SSL& s, const Ravelin::Vector3d& v);
    void operator=(const SSL& s);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, double dt, const Ravelin::Twistd& v) const;
    static double calc_dist(const SSL& o, const Ravelin::Point3d& p);
    static double calc_dist(const SSL& a, const SSL& b, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static double calc_dist(const SSL& a, const SSL& b, const std::pair<Ravelin::Quatd, Ravelin::Origin3d>& aTb, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static bool intersects(const SSL& a, const SSL& b);
    static bool intersects(const SSL& a, const SSL& b, const std::pair<Ravelin::Quatd, Ravelin::Origin3d>& aTb);
    static bool intersects(const SSL& a, const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q);
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q) const { return SSL::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const SSL& a, const Ravelin::Point3d& point, double tol = NEAR_ZERO);
    virtual bool outside(const Ravelin::Point3d& point, double tol = NEAR_ZERO) const { return SSL::outside(*this, point, tol); }
    boost::shared_ptr<SSL> get_this() { return boost::dynamic_pointer_cast<SSL>(shared_from_this()); }
    boost::shared_ptr<const SSL> get_this() const { return boost::dynamic_pointer_cast<const SSL>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    unsigned calc_size() const;
    virtual Ravelin::Point3d get_lower_bounds(const Ravelin::Pose3d& T);
    virtual Ravelin::Point3d get_upper_bounds(const Ravelin::Pose3d& T);

    /// Calculates (approximate?) volume of the SSL 
    virtual double calc_volume() const { return (p1 - p2).norm() * M_PI * radius * radius * radius; }

    /// The first point of the line segment 
    Ravelin::Point3d p1;

    /// The second point of the line segment
    Ravelin::Point3d p2;

    /// Radius of the spherical addition 
    double radius;

  private:
    static double calc_sq_dist(const LineSeg3& seg, const Ravelin::Point3d& q, Ravelin::Point3d& cp);
    static double calc_sq_dist(const LineSeg3& seg0, const LineSeg3& seg1, Ravelin::Point3d& cp0, Ravelin::Point3d& cp1);
    static double calc_sq_dist(const Ravelin::Point3d& origin, const Ravelin::Vector3d& dir, const LineSeg3& seg, Ravelin::Point3d& cp_line, Ravelin::Point3d& cp_seg, double& line_param);

}; // end class

// include inline functions
#include "SSL.inl"

} // end namespace

#endif

