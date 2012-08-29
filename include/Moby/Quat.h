/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_QUAT_H
#define _MOBY_QUAT_H

#include <cmath>
#include <Moby/VectorN.h>
#include <Moby/Vector3.h>
#include <Moby/MatrixN.h>

namespace Moby {

class Matrix3;
class Matrix4;
class AAngle;

/// Quaternion class used for orientation
/**
 * This class is used to represent orientation via unit quaternions.  Note, however, that the user is responsible
 * for normalizing the quaternions; normalization is not automatic.
 */
class Quat
{
  public:
    Quat();
    Quat(Real x, Real y, Real z, Real w);
    Quat(boost::shared_array<const Real>  q);
    Quat(const Quat& q);
    Quat(const VectorN& v);
    Quat(const Matrix3* m);
    Quat(const Matrix4* m);
    Quat(const AAngle* a);
    bool unit() const;
    static Quat zero();
    void conjugate();
    static Quat conjugate(const Quat& q);
    bool epsilon_equals(const Quat& q, Real epsilon);
    static bool epsilon_equals(const Quat& q1, const Quat& q2, Real epsilon);
    void slerp(const Quat& q, Real alpha);
    void lerp(const Quat& q, Real alpha);
    static Quat slerp(const Quat& q1, const Quat& q2, Real alpha);
    static Quat lerp(const Quat& q1, const Quat& q2, Real alpha);
    void inverse();
    static Quat inverse(const Quat& q);
    Real norm_sq() const;
    Real norm() const { return magnitude(); }
    Real norm_inf() const;
    void normalize();
    static Quat normalize(const Quat& q); 
    void set(const VectorN& v);
    void set(const AAngle* a);
    void set(const Matrix3* m);
    void set(const Matrix4* m);
    void calc_generalized_force(const Vector3& t, Real gt[4]) const;
    unsigned size() const { return 4; }
    Real& operator[](unsigned i);
    Real operator[](unsigned i) const;
    Quat operator-() const;
    Quat operator-(const Quat& q) const;
    Quat& operator-=(const Quat& q);
    Quat operator+(const Quat& q) const;
    Quat& operator+=(const Quat& q);
    Quat& operator=(const Quat& q);
    Quat operator*(const Quat& q) const;
    Quat operator/(const Quat& q) const;
    Vector3 operator*(const Vector3& v) const;
    Quat operator*(Real scalar) const;
    Quat operator/(Real scalar) const { return operator*(1.0/scalar); }
    Quat& operator*=(const Quat& q);
    Quat& operator*=(Real scalar);
    Quat& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    Real magnitude() const;
    static Quat deriv(const Quat& q, const Vector3& w);
    static Quat dderiv(const Quat& q, const Vector3& omega, const Vector3& alpha);
    static Vector3 to_omega(const Quat& q, const Quat& qd);
    static VectorN to_vect(const Quat& q);
    Vector3 G_mult(Real qw, Real qx, Real qy, Real qz) const;
    Quat G_transpose_mult(const Vector3& v) const;
    Quat L_transpose_mult(const Vector3& v) const;
    MatrixN& determine_G(MatrixN& G) const;
    MatrixN& determine_L(MatrixN& L) const;
    Vector3 L_mult(const Quat& q) const;
    Real* data() { return &x; }
    const Real* data() const { return &x; }
    
    /// First quaternion component
    Real x;

    /// Second quaternion component
    Real y;

    /// Third quaterion component
    Real z; 

    /// Fourth quaternion component
    Real w;

    /// Copies from the source quaternion
    Quat& copy_from(const Quat& source) { return operator=(source); }

  private:
    static Real sgn(Real x) { return (x >= (Real) 0.0) ? (Real) 1.0 : (Real) -1.0; }

    /// Computes a "safe" square root
    /**
     * If the input is negative, safe_sqrt() will return 0.
     * \note the result will be incorrect if the magnitude of the input is
     *       effectively greater than zero.
     */
    static Real safe_sqrt(Real x) { return (x < (Real) 0.0) ? (Real) 0.0 : std::sqrt(x); }
};

Quat operator*(Real scalar, const Quat& q);
std::ostream& operator<<(std::ostream& out, const Quat& q);

}
#endif
