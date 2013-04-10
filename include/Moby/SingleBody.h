/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SINGLE_BODY_H
#define _SINGLE_BODY_H

#include <Moby/DynamicBody.h>

namespace Moby {

/// Superclass for both rigid and deformable bodies 
class SingleBody : public DynamicBody 
{
  public:
    virtual ~SingleBody() {}
    virtual DynamicBodyPtr get_super_body() const;

    /// Integrates this body
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator)
    {
      if (!is_enabled())
        return;
      else
        DynamicBody::integrate(t, h, integrator);
    }

    /// Gets the position of the center-of-mass of the body
    virtual Ravelin::Point3d get_position() const = 0;

    /// Gets the velocity of the body
    virtual const Ravelin::Twistd& get_velocity() const = 0;

    /// Adds a wrench to the body
    virtual void add_wrench(const Ravelin::Wrenchd& w) = 0;

    /// Applies an impulse at a point on the body
    virtual void apply_impulse(const Ravelin::Wrenchd& w) = 0;

    /// Calculates the mass of the body
    virtual double calc_mass() const = 0;

    /// Gets the articulated body that this body is a part of (if any)
    virtual ArticulatedBodyPtr get_articulated_body() const = 0;

    /// Determines whether the body is enabled
    virtual bool is_enabled() const = 0;

    /// Calculates the velocity at a point on the body
    virtual Ravelin::Vector3d calc_point_vel(const Ravelin::Point3d& point) const = 0;

    /// Calculates the velocity at a point on the body in a given direction
    double calc_point_vel(const Ravelin::Point3d& point, const Ravelin::Vector3d& dir) { return calc_point_vel(point).dot(dir); }

  private:

    virtual double get_aspeed() const { return get_velocity().get_angular().norm(); }
}; // end class

} // end namespace

#endif

