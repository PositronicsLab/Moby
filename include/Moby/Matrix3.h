/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX3_H
#define _MATRIX3_H

#include <iostream>
#include <Moby/Vector3.h>

namespace Moby {

class Quat;
class AAngle;
class MatrixN;

/// A 3x3 matrix that may be used for orientation, inertia tensors, etc.
class Matrix3
{
  public:
    Matrix3() { }
    Matrix3(const Real* array);
    Matrix3(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22);
    Matrix3(const Quat* q);
    Matrix3(const Matrix3& source) { operator=(source); }
    Matrix3(const AAngle* a);
    Real norm_inf() const;
    unsigned size() const { return 3; }
    unsigned rows() const { return 3; }
    unsigned columns() const { return 3; }
    bool epsilon_equals(const Matrix3& m, Real epsilon) const;
    static bool epsilon_equals(const Matrix3& m1, const Matrix3& m2, Real epsilon);
    bool is_symmetric(Real tolerance) const;
    bool orthonormalize();
    bool is_orthonormal() const;    
    Real det() const;
    void inverse();
    static Vector3 calc_differential(const Matrix3& R1, const Matrix3& R2);
    static Matrix3 inverse(const Matrix3& m);
    void set_rot_X(Real angle);
    static Matrix3 rot_X(Real angle);
    void set_rot_Y(Real angle);
    static Matrix3 rot_Y(Real angle);
    void set_rot_Z(Real angle);
    static Matrix3 rot_Z(Real angle);
    void set(const AAngle* a);
    void set(const Quat* q);
    static Matrix3 skew_symmetric(Real a, Real b, Real c);
    static Matrix3 skew_symmetric(const Vector3& v);
    static Vector3 inverse_skew_symmetric(const Matrix3& R);
    static Matrix3 transpose(const Matrix3& m);
    void transpose();
    static bool valid_rotation(const Matrix3& R);
    static bool valid_rotation_scale(const Matrix3& R);
    static Matrix3 identity() { Matrix3 m; m.set_identity(); return m; }
    static Matrix3 zero() { Matrix3 m; m.set_zero(); return m; }
    void set_row(unsigned i, const Vector3& v);
    void set_row(unsigned i, const Real* source);
    void set_column(unsigned i, const Vector3& v);
    void set_column(unsigned i, const Real* source);
    Vector3 get_row(unsigned i) const;
    void get_row(unsigned i, Real* source) const;
    void get_row(unsigned i, Vector3& result) const;
    void get_column(unsigned i, Real* source) const;
    Vector3 get_column(unsigned i) const;
    void get_column(unsigned i, Vector3& result) const;
    void set_identity();
    void set_zero();
    Matrix3 mult(const Matrix3& m) const;
    Vector3 mult(const Vector3& v) const;
    Vector3 transpose_mult(const Vector3& v) const;
    Matrix3 mult_transpose(const Matrix3& m) const;
    Matrix3 transpose_mult_transpose(const Matrix3& m) const;
    Matrix3 transpose_mult(const Matrix3& m) const;
    static MatrixN* mult(const Matrix3& m1, const MatrixN* m2, MatrixN* result) { return m1.mult(m2, result); }
    MatrixN* mult(const MatrixN* m, MatrixN* result) const;
    MatrixN* mult_transpose(const MatrixN* m, MatrixN* result) const;
    MatrixN* transpose_mult_transpose(const MatrixN* m, MatrixN* result) const;
    MatrixN* transpose_mult(const MatrixN* m, MatrixN* result) const;
    Matrix3& operator=(const Matrix3& source);
    Matrix3& operator+=(const Matrix3& m);
    Matrix3& operator-=(const Matrix3& m);
    Matrix3& operator*=(const Matrix3& m) { return *this = *this * m; }
    Matrix3& operator*=(Real scalar);
    Matrix3& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    Vector3 operator*(const Vector3& v) const { return mult(v); } 
    Matrix3 operator+(const Matrix3& m) const { Matrix3 n = *this; n += m; return n; }
    Matrix3 operator-(const Matrix3& m) const { Matrix3 n = *this; n -= m; return n; }
    Matrix3 operator*(const Matrix3& m) const { return mult(m); }
    Matrix3 operator*(Real scalar) const { Matrix3 m = *this; m *= scalar; return m; }
    Matrix3 operator/(Real scalar) const { return operator*(1.0/scalar); }
    Matrix3 operator-() const; 
    Real& operator()(unsigned i, unsigned j) { assert(i < 3 && j < 3); return _data[j*3+i]; }
    Real operator()(unsigned i, unsigned j) const { assert(i < 3 && j < 3); return _data[j*3+i]; }

    /// Dummy method
    Matrix3& resize(unsigned rows, unsigned columns, bool preserve = false) { assert(rows == 3 && columns == 3); return *this; }

    /// Gets iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Gets iterator to the end of the data
    Real* end() { return &_data[0] + 9; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return &_data[0] + 9; }

    /// Gets constant pointer to the beginning of the matrix array
    const Real* data() const { return _data; }

    /// Gets pointer to the beginning of the matrix array
    Real* data() { return _data; }

  private:
    static bool rel_equal(Real x, Real y, Real tolerance);
    bool orthonormalize(Vector3& a, Vector3& b, Vector3& c);
    Real _data[9];
}; // end class

std::ostream& operator<<(std::ostream& out, const Matrix3& m);

} // end namespace


#endif
