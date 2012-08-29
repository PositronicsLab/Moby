/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX6_H
#define _MATRIX6_H

#include <Moby/SVector6.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/MatrixN.h>

namespace Moby {

class Quat;
class AAngle;
class SMatrix6N;

/// A 6x6 spatial algebra matrix typically used for dynamics calculations
class SMatrix6
{
  public:
    SMatrix6();
    SMatrix6(const SMatrix6& source) { operator=(source); }
    SMatrix6(const Real* source);
    SMatrix6(const MatrixN& source);
    unsigned size() const { return 6; }
    unsigned rows() const { return 6; }
    unsigned columns() const { return 6; }
    Matrix3 get_upper_left() const;
    Matrix3 get_upper_right() const;
    Matrix3 get_lower_left() const;
    Matrix3 get_lower_right() const;
    void set_upper_left(const Matrix3& m);
    void set_upper_right(const Matrix3& m);
    void set_lower_left(const Matrix3& m);
    void set_lower_right(const Matrix3& m);
    void set_identity();
    void set_zero();
    static SMatrix6 zero() { SMatrix6 m; m.set_zero(); return m; }
    static SMatrix6 identity() { SMatrix6 m; m.set_identity(); return m; }
    void transpose();
    SMatrix6 transpose(const SMatrix6& m);
    SMatrix6N* mult(const SMatrix6N* m, SMatrix6N* result) const;
    MatrixN& mult(const SMatrix6N* m, MatrixN& result) const;
    SMatrix6& operator=(const SMatrix6& source);
    SMatrix6& operator+=(const SMatrix6& m);
    SMatrix6& operator-=(const SMatrix6& m);
    SMatrix6& operator*=(const SMatrix6& m) { return *this = operator*(m); }
    SMatrix6& operator*=(Real scalar);
    SMatrix6& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    SMatrix6 operator+(const SMatrix6& m) const;
    SMatrix6 operator-(const SMatrix6& m) const;
    SMatrix6 operator*(const SMatrix6& m) const;
    SMatrix6 operator*(Real scalar) const { SMatrix6 m = *this; m *= scalar; return m; }
    SMatrix6 operator/(Real scalar) const { return operator*(1.0/scalar); }
    SVector6 operator*(const SVector6& v) const;
    SMatrix6 operator-() const;
    static SMatrix6 spatial_cross(const SVector6& v);
    static SMatrix6 calc_spatial_transform(const Matrix4& source, const Matrix4& target);
    static SMatrix6 calc_spatial_transform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix3& targetR, const Vector3& targetX);
    static SMatrix6 inverse_inertia(const SMatrix6& I);    
    Real& operator()(const unsigned i, const unsigned j) { assert(i < 6 && j < 6); return _data[j*6+i]; }
    Real operator()(const unsigned i, const unsigned j) const { assert(i < 6 && j < 6); return _data[j*6+i]; }

    /// Dummy method
    SMatrix6& resize(unsigned rows, unsigned columns, bool preserve = false) { assert(rows == 6 && columns == 6); return *this; }

    /// Gets iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Gets iterator to the end of the data
    Real* end() { return &_data[0] + 36; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return &_data[0] + 36; }

    /// Gets a constant pointer to the beginning of the matrix array 
    const Real* data() const { return _data; }

    /// Gets a pointer to the beginning of the matrix array
    Real* data() { return _data; }

  private:
    Real _data[36];
}; // end class

std::ostream& operator<<(std::ostream& out, const SMatrix6& m);

} // end namespace

#endif

