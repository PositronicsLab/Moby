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

    /// Integrates this body
    virtual void integrate(Real t, Real h, boost::shared_ptr<Integrator<VectorN> > integrator)
    {
      if (!is_enabled())
        return;
      else
        DynamicBody::integrate(t, h, integrator);
    }

    /// Gets the position of the center-of-mass of the body
    virtual const Vector3& get_position() const = 0;

    /// Gets the linear velocity of the center-of-mass of the body
    virtual const Vector3& get_lvel() const = 0;

    /// Gets the angular velocity of the body
    virtual const Vector3& get_avel() const = 0;

    /// Adds a force to the body
    virtual void add_force(const Vector3& f) = 0;

    /// Adds a force at a point on the body
    virtual void add_force(const Vector3& f, const Vector3& point) = 0;

    /// Applies an impulse at a point on the body
    virtual void apply_impulse(const Vector3& j, const Vector3& point) = 0;

    /// Calculates the mass of the body
    virtual Real calc_mass() const = 0;

    /// Gets the articulated body that this body is a part of (if any)
    virtual ArticulatedBodyPtr get_articulated_body() const = 0;

    /// Calculates the velocity at a point on the body
    virtual Vector3 calc_point_vel(const Vector3& point) const = 0;

    /// Calculates the velocity at a point on the body in a given direction
    Real calc_point_vel(const Vector3& point, const Vector3& dir) { return calc_point_vel(point).dot(dir); }

  private:

    virtual Real get_aspeed() const { return get_avel().norm(); }
}; // end class

} // end namespace

#endif

