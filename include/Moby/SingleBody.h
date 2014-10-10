/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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

    /// Gets the computation frame for the body
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_computation_frame() const = 0;

    /// Gets the mass of the body (for gravity calculation)
    virtual double get_mass() const = 0;

    /// Gets the pose of the body
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_pose() const = 0;

    /// Gets the acceleration of the body
    virtual const Ravelin::SAcceld& get_accel() = 0;

    /// Gets the velocity of the body
    virtual const Ravelin::SVelocityd& get_velocity() = 0;

    /// Applies an impulse at a point on the body
    virtual void apply_impulse(const Ravelin::SMomentumd& w) = 0;

    /// Calculates the mass of the body
    virtual double calc_mass() const = 0;

    /// Gets the articulated body that this body is a part of (if any)
    virtual ArticulatedBodyPtr get_articulated_body() const = 0;

    /// Determines whether the body is enabled
    virtual bool is_enabled() const = 0;

    /// Calculates the velocity at a point on the body *in the body frame*
    virtual Ravelin::Vector3d calc_point_vel(const Point3d& point) const = 0;

    /// Calculates the velocity at a point on the body in a given direction
    double calc_point_vel(const Point3d& point, const Ravelin::Vector3d& dir);

    /// Gets the maximum angular speed of this body (useful for collision detection)
    virtual double get_aspeed() 
    { 
      const Ravelin::SVelocityd& v = get_velocity();
      boost::shared_ptr<const Ravelin::Pose3d> F = get_pose();
      Ravelin::SVelocityd vi = Ravelin::Pose3d::transform(F, v);
      return vi.get_angular().norm(); 
    }

}; // end class

} // end namespace

#endif

