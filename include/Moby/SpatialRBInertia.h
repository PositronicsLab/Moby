/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_RB_INERTIA_H
#define _SPATIAL_RB_INERTIA_H

#include <Moby/SVector6.h>
#include <Moby/SMatrix6N.h>

namespace Moby {

/// A 6x6 spatial algebra matrix used for dynamics calculations
/** 
 * The matrix is represented by:
 * | J - hx*hx*m  hx*m |
 * | -hx*m        I*m  |
 * where hx is the skew symmetric matrix determined by h and I is the identity
 * matrix.
 */
class SpatialRBInertia
{
  public:
    SpatialRBInertia();
    SpatialRBInertia(Real m, const Vector3& h, const Matrix3& J);
    SpatialRBInertia(const SpatialRBInertia& source) { operator=(source); }
    void set_zero();
    static SpatialRBInertia zero() { SpatialRBInertia m; m.set_zero(); return m; }
    void to_matrix(Real data[]) const;
    SMatrix6N& mult(const SMatrix6N& m, SMatrix6N& result) const;
    SVector6 inverse_mult(const SVector6& v) const;
    SpatialRBInertia& operator=(const SpatialRBInertia& source);
    SpatialRBInertia& operator+=(const SpatialRBInertia& m);
    SpatialRBInertia& operator-=(const SpatialRBInertia& m);
    SpatialRBInertia& operator*=(const SpatialRBInertia& m) { return *this = operator*(m); }
    SpatialRBInertia& operator*=(Real scalar);
    SpatialRBInertia& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    SpatialRBInertia operator+(const SpatialRBInertia& m) const;
    SpatialRBInertia operator-(const SpatialRBInertia& m) const;
    SpatialRBInertia operator*(const SpatialRBInertia& m) const;
    SpatialRBInertia operator*(Real scalar) const { SpatialRBInertia m = *this; m *= scalar; return m; }
    SpatialRBInertia operator/(Real scalar) const { return operator*(1.0/scalar); }
    SVector6 operator*(const SVector6& v) const;
    SpatialRBInertia operator-() const;

    /// The rigid body mass
    Real m;

    /// The rigid body offset vector
    Vector3 h;

    /// The rigid body moment of inertia matrix
    Matrix3 J;
}; // end class

std::ostream& operator<<(std::ostream& out, const SpatialRBInertia& m);

} // end namespace

#endif

