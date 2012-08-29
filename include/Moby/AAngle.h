/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _AXIS_ANGLE_H
#define _AXIS_ANGLE_H

#include <string>
#include <cmath>
#include <boost/shared_array.hpp>
#include <Moby/Constants.h>
#include <Moby/Types.h>

namespace Moby {

class Quat;
class Matrix3;
class Matrix4;
class Vector3;
class Vector4;

/// Class for representation of orientation by an angle around an axis
class AAngle
{
  public:
    AAngle();
    AAngle(const AAngle& source);
    AAngle(boost::shared_array<const Real> array);
    AAngle(const Quat* q);
    AAngle(Real x, Real y, Real z, Real angle);
    AAngle(const Vector3*  v, Real angle);
    AAngle(const VectorN* v);
    AAngle(const Matrix3* m);
    AAngle(const Matrix3* m, const Vector3* axis);
    AAngle(const Matrix4* m);
    bool epsilon_equals(const AAngle& a, Real epsilon) const;
    static bool epsilon_equals(const AAngle& a1, const AAngle& a2, Real epsilon);
    void set(boost::shared_array<const Real>  array);
    void set(Real x, Real y, Real z, Real angle);
    void set(const VectorN* v);
    void set(const Vector3* v, Real angle);
    void set(const Matrix3* m);
    void set(const Matrix3* m, const Vector3* axis);
    void set(const Matrix4* m);
    void set(const Quat* q);
    void operator=(const AAngle& source);
    AAngle operator*(const AAngle& a) const;
    void operator*=(const AAngle& a);

    /// Does a "safe" square root - input values sufficiently close to -0 are returned as 0
    Real safe_sqrt(Real x) { return (x < (Real) 0.0 && x > -NEAR_ZERO) ? (Real) 0.0 : std::sqrt(x); }

    Real angle;
    Real x;
    Real y;
    Real z;
};

std::ostream& operator<<(std::ostream& out, const AAngle& a);

}

#endif

