/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR6_H
#define _VECTOR6_H

#include <Moby/Vector3.h>
#include <Moby/VectorN.h>

namespace Moby {

/// A 6-dimensional floating-point vector for use with spatial algebra
/**
 * Note that spatial algebra defines the dot product in an unusual manner: if vector x = [a; b] and 
 * vector y = [c; d] then x'y = [b'; a'][c d] = dot(b,c) + dot(a,d).
 */
class SVector6
{
  public:
    SVector6() {}
    SVector6(Real x, Real y, Real z, Real a, Real b, Real c);
    SVector6(const Real* array);
    SVector6(const Vector3& upper, const Vector3& lower);
    SVector6(const VectorN& v);
    unsigned size() const { return 6; }
    static SVector6 spatial_cross(const SVector6& v1, const SVector6& v2);
    static SVector6 zero() { return SVector6(0,0,0,0,0,0); }
    Real dot(const SVector6& v) const { return dot(*this, v); }
    static Real dot(const SVector6& v1, const SVector6& v2);
    void set_lower(const Vector3& lower);
    void set_upper(const Vector3& upper);
    Vector3 get_lower() const;
    Vector3 get_upper() const;
    SVector6& operator=(const SVector6& source);
    void transpose();
    static SVector6 transpose(const SVector6& v);
    Real& operator[](const unsigned i) { assert(i < 6); return _data[i]; }
    Real operator[](const unsigned i) const { assert(i < 6); return _data[i]; }
    Real* data() { return _data; }
    const Real* data() const { return _data; }
    SVector6 operator-() const;
    SVector6 operator*(Real scalar) const { SVector6 v = *this; return v*= scalar; }
    SVector6 operator/(Real scalar) const { SVector6 v = *this; return v/= scalar; }
    SVector6& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    SVector6& operator*=(Real scalar);
    SVector6& operator-=(const SVector6& v);
    SVector6& operator+=(const SVector6& v);
    SVector6 operator+(const SVector6& v) const;
    SVector6 operator-(const SVector6& v) const;

    /// Get iterator to the start of the data
    Real* begin() { return _data; }

    /// Get iterator to the start of the data
    const Real* begin() const { return _data; }

    /// Get iterator to the end of the data
    Real* end() { return _data + 6; }

    /// Get iterator to the end of the data
    const Real* end() const { return _data + 6; }

  private:
    Real _data[6];
}; // end class
} // end namespace

#endif

