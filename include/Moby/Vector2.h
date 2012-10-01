/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR2_H
#define _VECTOR2_H

#include <boost/shared_array.hpp>
#include <limits>
#include <iostream>
#include <cmath>
#include <Moby/Types.h>

namespace Moby {

class Matrix2;

/// A two-dimensional floating point vector used for computational geometry calculations
class Vector2
{
  public:
    Vector2() {}
    Vector2(Real x, Real y);
    Vector2(const Real* array);
    Real dot(const Vector2& v) const { return v[0]*_data[0] + v[1]*_data[1]; }
    static Real dot(const Vector2& v1, const Vector2& v2) { return v1[0]*v2[0] + v1[1]*v2[1]; }
    Real norm_inf() const { return std::max(std::fabs(_data[0]), std::fabs(_data[1])); }
    Real norm() const { return std::sqrt(norm_sq()); }
    Real norm_sq() const { return dot(*this, *this); }
    void normalize() { assert(norm() > std::numeric_limits<Real>::epsilon()); operator/=(norm()); }
    static Vector2 normalize(const Vector2& v) { Vector2 w = v; w.normalize(); return w; }
    unsigned size() const { return 2; }
    static Real norm(const Vector2& v) { return std::sqrt(norm_sq(v)); }
    static Real norm_sq(const Vector2& v) { return v.dot(v); }
    Matrix2* outer_prod(const Vector2& v, Matrix2* m) const { return outer_prod(*this, v, m); }
    static Matrix2* outer_prod(const Vector2& v1, const Vector2& v2, Matrix2* m);
    void set_zero() { _data[0] = _data[1] = 0.0; }
    void set_one() { _data[0] = _data[1] = 1.0; }
    static Vector2 zero() { return Vector2(0.0, 0.0); }
    Vector2& operator=(const Vector2& v) { _data[0] = v[0]; _data[1] = v[1]; return *this; }
    Vector2 operator+(const Vector2& v) const { return Vector2(_data[0] + v[0], _data[1] + v[1]); }
    Vector2 operator-(const Vector2& v) const { return Vector2(_data[0] - v[0], _data[1] - v[1]); }
    Vector2& operator+=(const Vector2& v) { _data[0] += v[0]; _data[1] += v[1];  return *this; }
    Vector2& operator-=(const Vector2& v) { _data[0] -= v[0]; _data[1] -= v[1];  return *this; }
    Vector2& operator*=(Real scalar) { _data[0] *= scalar; _data[1] *= scalar; return *this; }
    Vector2& operator/=(Real scalar) { _data[0] /= scalar; _data[1] /= scalar; return *this; }
    Vector2 operator*(Real scalar) const { Vector2 v = *this; v *= scalar; return v; }
    Vector2 operator/(Real scalar) const { Vector2 v = *this; v /= scalar; return v; }
    Vector2 operator-() const { return Vector2(-_data[0], -_data[1]); }
    Real& operator[](const unsigned i) { assert(i < 2); return _data[i]; }
    Real operator[](const unsigned i) const { assert(i < 2); return _data[i]; }
    Real* data() { return _data; }
    const Real* data() const { return _data; }
    unsigned rows() const { return 2; }
    unsigned columns() const { return 1; }
    Vector2& resize(unsigned m, unsigned n) { assert(m == 2 && n == 1); return *this; }

    /// Computes the "perp" operator
    Vector2 perp() const { return Vector2(_data[1], -_data[0]); }

    /// Computes the dot "perp" product
    Real dot_perp(const Vector2& v) const { return dot(v.perp()); }

    /// Get iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Get iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Get iterator to the end of the data
    Real* end() { return &_data[0] + 2; }

    /// Get iterator to the end of the data
    const Real* end() const { return &_data[0] + 2; }

    /// Copies from the source vector
    Vector2& copy_from(const Vector2& source) { return operator=(source); }

  private:
    Real _data[2];
}; // end class

inline Vector2 operator*(Real scalar, const Vector2& v) { return v * scalar; }

/// Writes a Vector2 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const Vector2& v)
{
  out << '[' << v[0] << ',' << ' ' << v[1] << ']' << ' ';
  return out;
};

} // end namespace

#endif
