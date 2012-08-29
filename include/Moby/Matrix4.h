/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX4_H
#define _MATRIX4_H

#include <Moby/Vector3.h>

class SbMatrix;

namespace Moby {

class Quat;
class AAngle;
class Matrix3;

/// A 4x4 transformation matrix.
/**
 * This transformation matrix supports rotation and translation, but not
 * scaling; if the rotation component is set to scale also, undefined
 * behavior will result.
 * The underlying data is stored in a column-major 
 * (e.g., the element at row 1, column 0 is element "1") though the transform
 * itself is non-OpenGL style (i.e., translation is in the fourth column).
 */
class Matrix4
{
  public:
    Matrix4();
    Matrix4(const Matrix4& source) { operator=(source); }
    Matrix4(const Real* array);
    Matrix4(const AAngle* a);
    Matrix4(const Matrix3* m);
    Matrix4(const Quat* q);
    Matrix4(const AAngle* a, const Vector3* v);
    Matrix4(const Matrix3* m, const Vector3* v);
    Matrix4(const Quat* q, const Vector3* v);
    Matrix4(const Vector3* v);
    Matrix4(Real m00, Real m01, Real m02, Real m03, Real m10, Real m11, Real m12, Real m13, Real m20, Real m21, Real m22, Real m23);
    static Matrix4 identity() { Matrix4 T; T.set_identity(); return T; }
    unsigned size() const { return 4; }
    unsigned rows() const { return 4; }
    unsigned columns() const { return 4; }
    bool epsilon_equals(const Matrix4& m, Real epsilon) const;
    static bool epsilon_equals(const Matrix4& m1, const Matrix4& m2, Real epsilon);
    static Matrix4 interpolate(const Matrix4& m1, const Matrix4& m2, Real t);
    Matrix3 get_rotation() const;
    void get_rotation(Matrix3* m) const; 
    Vector3 get_translation() const;
    Vector3 mult_point(const Vector3& v) const;
    Vector3 mult_vector(const Vector3& v) const;
    Vector3 inverse_mult_point(const Vector3& v) const;
    Vector3 transpose_mult_vector(const Vector3& v) const;
    void set_identity();
    void inverse_transform();
    static Matrix4 inverse_transform(const Matrix4& m);
    void get_translation(Real* array) const;
    void get_translation(Vector3& v) const;
    void get_translation(Real& x, Real& y, Real& z) const;
    void set(const Real* array);
    void set(const AAngle* a, const Vector3*  v);
    void set(const Matrix3* m, const Vector3*  v);
    void set(const Quat* q, const Vector3*  v);
    void set_rotation(const AAngle* a);
    void set_rotation(const Quat* q);
    void set_rotation(const Matrix3* m);
    void set_translation(Real x, Real y, Real z);
    void set_translation(const Real* array);
    void set_translation(const Vector3& v);
    Matrix4& operator=(const Matrix4& source);
    static bool valid_transform(const Matrix4& m);
    Real& operator()(unsigned i, unsigned j) { assert(i < 3 && j < 4); return _data[j*3+i]; }
    Real operator()(unsigned i, unsigned j) const { assert(i < 3 && j < 4); return _data[j*3+i]; }
    Matrix4 operator*(const Matrix4& m) const;
    const Real* data() const { return _data; }
    Real* data() { return _data; }

    /// Gets iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Gets iterator to the end of the data
    Real* end() { return &_data[0] + 12; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return &_data[0] + 12; }

  private:
    Real _data[4*3]; 
}; // end class

std::ostream& operator<<(std::ostream& out, const Matrix4& m);

} // end namespace

#endif

