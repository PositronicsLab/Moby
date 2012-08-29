/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR3_H
#define _VECTOR3_H

#include <stdexcept>
#include <limits>
#include <cmath>
#include <iostream>
#include <Moby/Types.h>

namespace Moby {

/// A three-dimensional floating point vector
class Vector3
{
  public:
    Vector3() {}
    Vector3(Real x, Real y, Real z);
    Vector3(const Real* array);
    Vector3(const Vector3& source) { operator=(source); }
    void swap(Vector3& source);
    Real dot(const Vector3& v) const { return v[0]*_data[0] + v[1]*_data[1] + v[2]*_data[2]; }
    static Real dot(const Vector3& v1, const Vector3& v2) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
    void normalize() { assert(norm() > std::numeric_limits<Real>::epsilon()); operator/=(norm()); }
    void normalize_or_zero() { Real nrm = norm(); if (nrm > std::numeric_limits<Real>::epsilon()) operator/=(nrm); else _data[0] = _data[1] = _data[2] = (Real) 0.0; }
    static Vector3 normalize(const Vector3& v) { Vector3 w = v; w.normalize(); return w; }
    unsigned size() const { return 3; }
    Vector3& resize(unsigned N, bool keep = true) { if (N != 3) throw std::runtime_error("Can't resize a Vector3 to size other than 3!"); }
    bool is_finite() const;
    Real norm_inf() const { return std::max(std::fabs(_data[0]), std::max(std::fabs(_data[1]), std::fabs(_data[2]))); }
    Real norm() const { return std::sqrt(norm_sq()); }
    Real norm_sq() const { return dot(*this, *this); }
    static Real norm(const Vector3& v) { return std::sqrt(norm_sq(v)); }
    static Real norm_sq(const Vector3& v) { return v.dot(v); }
    Matrix3* outer_prod(const Vector3& v, Matrix3* m) const { return outer_prod(*this, v, m); }
    static Matrix3* outer_prod(const Vector3& v1, const Vector3& v2, Matrix3* m);
    void set_zero() { _data[0] = _data[1] = _data[2] = 0.0; }
    void set_one() { _data[0] = _data[1] = _data[2] = 1.0; }
    static Vector3 zero() { return Vector3(0.0, 0.0, 0.0); }
    static Vector3 one() { return Vector3(1.0, 1.0, 1.0); }
    bool operator<(const Vector3& v) const;
    Vector3& operator=(const Vector3& v) { _data[0] = v[0]; _data[1] = v[1]; _data[2] = v[2]; return *this; }
    Vector3 operator+(const Vector3& v) const { return Vector3(_data[0] + v[0], _data[1] + v[1], _data[2] + v[2]); }
    Vector3 operator-(const Vector3& v) const { return Vector3(_data[0] - v[0], _data[1] - v[1], _data[2] - v[2]); }
    Vector3& operator+=(const Vector3& v) { _data[0] += v[0]; _data[1] += v[1]; _data[2] += v[2]; return *this; }
    Vector3& operator-=(const Vector3& v) { _data[0] -= v[0]; _data[1] -= v[1]; _data[2] -= v[2]; return *this; }
    Vector3& operator*=(Real scalar) { _data[0] *= scalar; _data[1] *= scalar; _data[2] *= scalar; return *this; }
    Vector3 operator*(Real scalar) const { Vector3 v = *this; v *= scalar; return v; }
    Vector3 operator/(Real scalar) const { Vector3 v = *this; v /= scalar; return v; }
    Vector3& operator/=(Real scalar) { _data[0] /= scalar; _data[1] /= scalar; _data[2] /= scalar; return *this; }
    Vector3 operator-() const { return Vector3(-_data[0], -_data[1], -_data[2]); }
    static Vector3 cross(const Vector3& v1, const Vector3& v2);
    static Vector3 determine_orthogonal_vec(const Vector3& v);
    static void determine_orthonormal_basis(const Vector3& v1, Vector3& v2, Vector3& v3);
    Real& operator[](const unsigned i) { assert(i < 3); return _data[i]; }
    Real operator[](const unsigned i) const { assert(i < 3); return _data[i]; }
    Real* data() { return _data; }
    const Real* data() const { return _data; }

    /// Get iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Get iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Get iterator to the end of the data
    Real* end() { return &_data[0] + 3; }

    /// Get iterator to the end of the data
    const Real* end() const { return &_data[0] + 3; }

    /// Copies from the source vector
    Vector3& copy_from(const Vector3& source) { return operator=(source); }

  private:
    static bool rel_equal(Real x, Real y);
    Real _data[3];
}; // end class

inline Vector3 operator*(Real scalar, const Vector3& v) { return v * scalar; }

/// Writes a Vector3 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const Vector3& v)
{
  out << "[" << v[0] << ", " << v[1] << ", " << v[2] << "] ";
  return out;
};

} // end namespace

#endif

