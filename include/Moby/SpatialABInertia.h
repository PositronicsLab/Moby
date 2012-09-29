/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_AB_INERTIA_H
#define _SPATIAL_AB_INERTIA_H

#include <Moby/SVector6.h>
#include <Moby/SpatialRBInertia.h>
#include <Moby/Matrix3.h>

namespace Moby {

/// A 6x6 spatial algebra matrix typically used for dynamics calculations
class SpatialABInertia
{
  public:
    SpatialABInertia();
    SpatialABInertia(const Matrix3& M, const Matrix3& H, const Matrix3& J);
    SpatialABInertia(const SpatialABInertia& source) { operator=(source); }
    SpatialABInertia(const SpatialRBInertia& source) { operator=(source); }
    SpatialABInertia(const MatrixN& M);
    void set_zero();
    static SpatialABInertia zero() { SpatialABInertia m; m.set_zero(); return m; }
    SMatrix6N& mult(const SMatrix6N& m, SMatrix6N& result) const;
    SVector6 inverse_mult(const SVector6& v) const;
    SpatialABInertia& operator=(const SpatialRBInertia& source);
    SpatialABInertia& operator=(const SpatialABInertia& source);
    SpatialABInertia& operator+=(const SpatialABInertia& m);
    SpatialABInertia& operator-=(const SpatialABInertia& m);
    SpatialABInertia& operator*=(const SpatialABInertia& m) { return *this = operator*(m); }
    SpatialABInertia& operator*=(Real scalar);
    SpatialABInertia& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    SpatialABInertia operator+(const SpatialRBInertia& m) const;
    SpatialABInertia& operator+=(const SpatialRBInertia& m) { return *this = *this + m; }
    SpatialABInertia operator+(const SpatialABInertia& m) const;
    SpatialABInertia operator-(const SpatialABInertia& m) const;
    SpatialABInertia operator*(const SpatialABInertia& m) const;
    SpatialABInertia operator*(Real scalar) const { SpatialABInertia m = *this; m *= scalar; return m; }
    SpatialABInertia operator/(Real scalar) const { return operator*(1.0/scalar); }
    SVector6 operator*(const SVector6& v) const;
    SpatialABInertia operator-() const;
    static SpatialABInertia inverse_inertia(const SpatialABInertia& I);    
    static SpatialABInertia mult(const MatrixN& M, const SpatialABInertia& I);

    /// The upper left / lower right hand matrix 'H'
    Matrix3 H;

    /// The lower left matrix
    Matrix3 J;

    /// The upper right matrix
    Matrix3 M;
}; // end class

std::ostream& operator<<(std::ostream& out, const SpatialABInertia& m);

} // end namespace

#endif

