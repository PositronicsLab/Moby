/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX2_H
#define _MATRIX2_H

#include <iostream>
#include <Moby/cblas.h>
#include <Moby/Vector2.h>
#include <Moby/MissizeException.h>

namespace Moby {

/// A general 2x2 matrix
class Matrix2
{
  public:
    Matrix2() {}
    Matrix2(const Real* array);
    Matrix2(Real m00, Real m01, Real m10, Real m11);
    Matrix2(const Matrix2& source) { operator=(source); }
    unsigned size() const { return 2; }
    unsigned rows() const { return 2; }
    unsigned columns() const { return 2; }
    bool is_symmetric(Real tolerance) const;
    bool orthonormalize();
    bool is_orthonormal() const;    
    Real det() const;
    void inverse();
    static Matrix2 inverse(const Matrix2& m);
    void set_rot_Z(Real angle);
    static Matrix2 rot_Z(Real angle);
    static Matrix2 transpose(const Matrix2& m);
    void transpose();
    static bool valid_rotation(const Matrix2& R);
    static Matrix2 identity() { Matrix2 m; m.set_identity(); return m; }
    static Matrix2 zero() { Matrix2 m; m.set_zero(); return m; }
    void set_row(unsigned i, const Vector2& v);
    void set_row(unsigned i, const Real* source);
    void set_column(unsigned i, const Vector2& v);
    void set_column(unsigned i, const Real* source);
    Vector2 get_row(unsigned i) const;
    void get_row(unsigned i, Real* source) const;
    void get_row(unsigned i, Vector2& result) const;
    void get_column(unsigned i, Real* source) const;
    Vector2 get_column(unsigned i) const;
    void get_column(unsigned i, Vector2& result) const;
    void set_identity();
    void set_zero();
    Matrix2 mult(const Matrix2& m) const;
    Vector2 mult(const Vector2& v) const;
    Vector2 transpose_mult(const Vector2& v) const;
    Matrix2 mult_transpose(const Matrix2& m) const;
    Matrix2 transpose_mult_transpose(const Matrix2& m) const;
    Matrix2 transpose_mult(const Matrix2& m) const;
    Matrix2& operator=(const Matrix2& source);
    Matrix2& operator+=(const Matrix2& m);
    Matrix2& operator-=(const Matrix2& m);
    Matrix2& operator*=(const Matrix2& m) { return *this = *this * m; }
    Matrix2& operator*=(Real scalar);
    Matrix2& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    Vector2 operator*(const Vector2& v) const { return mult(v); } 
    Matrix2 operator+(const Matrix2& m) const { Matrix2 n = *this; n += m; return n; }
    Matrix2 operator-(const Matrix2& m) const { Matrix2 n = *this; n -= m; return n; }
    Matrix2 operator*(const Matrix2& m) const { return mult(m); }
    Matrix2 operator*(Real scalar) const { Matrix2 m = *this; m *= scalar; return m; }
    Matrix2 operator/(Real scalar) const { return operator*(1.0/scalar); }
    Matrix2 operator-() const; 
    Real& operator()(unsigned i, unsigned j) { assert(i < 2 && j < 2); return _data[j*2+i]; }
    Real operator()(unsigned i, unsigned j) const { assert(i < 2 && j < 2); return _data[j*2+i]; }

    /// Dummy method
    Matrix2& resize(unsigned rows, unsigned columns, bool preserve = false) { assert(rows == 2 && columns == 2); return *this; }

    /// Gets iterator to the start of the data
    Real* begin() { return _data; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return _data; }

    /// Gets iterator to the end of the data
    Real* end() { return _data + 4; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return _data + 4; }

    /// Gets a constant pointer to the beginning of the matrix array 
    const Real* data() const { return _data; }

    /// Gets a pointer to the beginning of the matrix array
    Real* data() { return _data; }

    template <class T>
    T& transpose_mult_transpose(const T& x, T& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (cols != this->rows())
        throw MissizeException();
      y.resize(this->columns(), rows);
      CBLAS::gemm(CblasTrans, CblasTrans, 2, rows, 2, *this, 2, x, rows, (Real) 1.0, (Real) 0.0, y, 2); 
      return y;
    }

    template <class T>
    T& mult_transpose(const T& x, T& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (cols != this->columns())
        throw MissizeException();
      y.resize(this->rows(), rows);
      CBLAS::gemm(CblasNoTrans, CblasTrans, 2, rows, 2, *this, x, rows, (Real) 1.0, (Real) 0.0, y, 2); 
      return y;
    }

    template <class T>
    T& transpose_mult(const T& x, T& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (rows != this->columns())
        throw MissizeException();
      y.resize(this->columns(), cols);
      if (cols > 1)
        CBLAS::gemm(CblasTrans, CblasNoTrans, 2, cols, 2, *this, 2, x, rows, (Real) 1.0, (Real) 0.0, y, 2); 
      else
        CBLAS::gemv(CblasTrans, 2, 2, *this, 2, x, 1, (Real) 1.0, (Real) 0.0, y, 1);
      return y;
    }

    template <class T>
    T& mult(const T& x, T& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (rows != this->columns())
        throw MissizeException();
      y.resize(this->rows(), cols);
      if (cols > 1)
        CBLAS::gemm(CblasNoTrans, CblasNoTrans, 2, cols, 2, *this, x, (Real) 1.0, (Real) 0.0, y); 
      else
        CBLAS::gemv(CblasNoTrans, 2, 2, *this, 2, x, 1, (Real) 1.0, (Real) 0.0, y, 1);
      return y;
    }

  private:
    static bool rel_equal(Real x, Real y, Real tolerance);
    bool orthonormalize(Vector2& a, Vector2& b);
    Real _data[4];
}; // end class

std::ostream& operator<<(std::ostream& out, const Matrix2& m);

} // end namespace


#endif
