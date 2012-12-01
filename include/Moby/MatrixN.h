/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIXN_H
#define _MATRIXN_H

#include <Moby/cblas.h>
#include <Moby/VectorN.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/BlockIterator.h>
#include <Moby/InvalidIndexException.h>

namespace Moby {
        
/// A generic, possibly non-square matrix.  
/**
 * The underlying data is stored in column-major format (e.g., the element at row 1, column 0 is the second element in the flat array).
 */
class MatrixN
{
  public:
    MatrixN();
    MatrixN(unsigned rows, unsigned columns);
    MatrixN(const MatrixN& source);
    MatrixN(MatrixN& source);
    MatrixN(unsigned rows, unsigned columns, const Real* array);
    MatrixN(unsigned rows, unsigned columns, boost::shared_array<Real> array);
    MatrixN(const VectorN& v, bool transpose=false);
    MatrixN(const Matrix3& m);
    MatrixN(const Matrix4& m);
    MatrixN& set_identity();
    MatrixN& set_identity(unsigned sz);
    bool is_symmetric(Real tolerance = -1.0) const;
    static MatrixN identity(unsigned dim);
    Real norm_inf() const;
    static MatrixN construct_variable(unsigned rows, unsigned cols, ...);
    virtual ~MatrixN() { }
    MatrixN& zero_upper_triangle();
    MatrixN& zero_lower_triangle();
    void set_sub_row_block(unsigned row, ...);
    static VectorN diag_mult(const VectorN& d, const VectorN& v);
    static MatrixN diag_mult(const VectorN& d, const MatrixN& m);
    static MatrixN diag_mult_transpose(const VectorN& d, const MatrixN& m);
    static VectorN& diag_mult(const VectorN& d, const VectorN& v, VectorN& result);
    static MatrixN& diag_mult(const VectorN& d, const MatrixN& m, MatrixN& result);
    static MatrixN& diag_mult_transpose(const VectorN& d, const MatrixN& m, MatrixN& result);
    MatrixN select_square(const std::vector<bool>& indices) const;
    MatrixN& select_square(const std::vector<bool>& indices, MatrixN& result) const;
    MatrixN& set(const VectorN& v, bool transpose=false);
    MatrixN& set(unsigned rows, unsigned columns, const Real* array);
    MatrixN& set(unsigned rows, unsigned columns, boost::shared_array<Real> array);
    MatrixN& set(const Matrix3& m);
    MatrixN& set(const Matrix4& m);
    BlockIterator block_start(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end);
    BlockIterator block_end(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end);
    VectorN get_row(unsigned i) const;
    VectorN get_column(unsigned i) const;
    VectorN& get_column(unsigned i, VectorN& result) const;
    void get_column(unsigned i, Real* output) const;
    VectorN& get_row(unsigned i, VectorN& result) const;
    void get_row(unsigned i, Real* output) const;
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    virtual MatrixN& resize(unsigned rows, unsigned columns, bool preserve = false);
    MatrixN& remove_row(unsigned i);
    MatrixN& remove_column(unsigned i);
    MatrixN& set_row(unsigned i, const VectorN& v);
    MatrixN& set_row(unsigned i, const Vector3& v);
    MatrixN& set_row(unsigned i, const Real* source);
    MatrixN& set_column(unsigned i, const VectorN& v);
    MatrixN& set_column(unsigned i, const Vector3& v);
    MatrixN& set_column(unsigned i, const Real* source);
    MatrixN& set_sub_mat(unsigned row_start, unsigned col_start, const VectorN& v, bool transpose = false);
    MatrixN get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end) const;
    VectorN& get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end, VectorN& v) const;
    bool epsilon_equals(const MatrixN& m, Real epsilon);
    static bool epsilon_equals(const MatrixN& m1, const MatrixN& m2, Real epsilon);
    static MatrixN zero(unsigned rows, unsigned columns);
    MatrixN& negate();
    MatrixN& set_zero();
    virtual MatrixN& transpose();
    static MatrixN transpose(const MatrixN& m);
    static MatrixN& transpose(const MatrixN& m, MatrixN& result);
    virtual MatrixN& operator=(const MatrixN& source);
    virtual MatrixN& operator=(MatrixN& source);
    VectorN operator*(const VectorN& v) const { return mult(v); }
    static MatrixN& mult(const MatrixN& m1, const MatrixN& m2, MatrixN& result) { return m1.mult(m2, result); }
    MatrixN mult(const MatrixN& m) const;
    MatrixN mult_transpose(const MatrixN& m) const;
    MatrixN transpose_mult(const MatrixN& m) const;
    MatrixN transpose_mult_transpose(const MatrixN& m) const;
    VectorN mult(const VectorN& v) const;
    VectorN transpose_mult(const VectorN& v) const;
    MatrixN operator*(const MatrixN& m) const { return mult(m); }
    MatrixN operator*(Real scalar) const;
    MatrixN operator/(Real scalar) const { return operator*(1.0/scalar); }
    MatrixN operator+(const MatrixN& m) const;
    MatrixN& operator+=(const MatrixN& m);
    MatrixN operator-(const MatrixN& m) const;
    MatrixN operator-() const;
    MatrixN& operator-=(const MatrixN& m);
    MatrixN& operator/=(Real scalar);
    MatrixN& operator*=(Real scalar);
    MatrixN& operator*=(const MatrixN& m);
    const Real& operator()(unsigned i, unsigned j) const { assert(i < _rows && j < _columns); return _data[j*_rows+i]; }
    Real& operator()(const unsigned i, const unsigned j) { assert(i < _rows && j < _columns); return _data[j*_rows+i]; }
//    Real operator()(const unsigned i, const unsigned j) const { assert(i < _rows && j < _columns); return _data[j*_rows+i]; }
    Real* data() { return _data.get(); }
    const Real* data() const { return _data.get(); }    

    /// Sets this to a m x n sized zero matrix
    MatrixN& set_zero(unsigned m, unsigned n) { return resize(m,n).set_zero(); }

    /// Gets iterator to the start of the data
    Real* begin() { return &_data[0]; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return &_data[0]; }

    /// Gets iterator to the end of the data
    Real* end() { return &_data[0] + _rows*_columns; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return &_data[0] + _rows*_columns; }

    template <class ForwardIterator>
    MatrixN select_rows(ForwardIterator row_start, ForwardIterator row_end) const;

    template <class ForwardIterator>
    MatrixN& select_rows(ForwardIterator row_start, ForwardIterator row_end, MatrixN& m) const;

    template <class ForwardIterator>
    MatrixN select_columns(ForwardIterator col_start, ForwardIterator col_end) const;

    template <class ForwardIterator>
    MatrixN& select_columns(ForwardIterator col_start, ForwardIterator col_end, MatrixN& m) const;

    template <class ForwardIterator1, class ForwardIterator2>
    MatrixN select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end) const;

    template <class ForwardIterator1, class ForwardIterator2>
    MatrixN& select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end, MatrixN& m) const;

    template <class ForwardIterator1, class ForwardIterator2>
    VectorN select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end) const;

    template <class ForwardIterator1, class ForwardIterator2>
    VectorN& select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end, VectorN& v) const;

    template <class ForwardIterator>
    MatrixN select_square(ForwardIterator start, ForwardIterator end) const;

    template <class ForwardIterator>
    MatrixN& select_square(ForwardIterator start, ForwardIterator end, MatrixN& m) const;

    template <class M>
    MatrixN& set_sub_mat(unsigned row_start, unsigned col_start, const M& m, bool transpose = false);

    template <class M>
    M& get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end, M& m, bool transpose = false) const;

    /// Swaps this matrix with another
    void swap(MatrixN& source) { std::swap(_data, source._data); std::swap(_rows, source._rows); std::swap(_columns, source._columns); }

    /// Copies FROM the source matrix
    MatrixN& copy_from(const MatrixN& source) { operator=(source); return *this;}

    template <class T, class U>
    U& transpose_mult_transpose(const T& x, U& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (cols != this->rows())
        throw MissizeException();
      y.resize(this->columns(), rows);
      if (_rows == 0 || _columns == 0 || rows == 0)
      {
        y.set_zero();
        return y;
      }
      CBLAS::gemm(CblasTrans, CblasTrans, _columns, rows, _rows, *this, _rows, x, rows, (Real) 1.0, (Real) 0.0, y, _columns); 
      return y;
    }

    template <class T, class U>
    U& mult_transpose(const T& x, U& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (cols != this->columns())
        throw MissizeException();
      y.resize(this->rows(), rows);
      if (_rows == 0 || _columns == 0 || rows == 0)
      {
        y.set_zero();
        return y;
      }
      CBLAS::gemm(CblasNoTrans, CblasTrans, _rows, rows, _columns, *this, _rows, x, rows, (Real) 1.0, (Real) 0.0, y, _rows); 
      return y;
    }

    template <class T, class U>
    U& transpose_mult(const T& x, U& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (rows != this->rows())
        throw MissizeException();
      y.resize(this->columns(), cols);
      if (_rows == 0 || _columns == 0 || cols == 0)
      {
        y.set_zero();
        return y;
      }
      if (cols > 1)
        CBLAS::gemm(CblasTrans, CblasNoTrans, _columns, cols, _rows, *this, _rows, x, rows, (Real) 1.0, (Real) 0.0, y, _columns); 
      else
        CBLAS::gemv(CblasTrans, _rows, _columns, *this, _rows, x, 1, (Real) 1.0, (Real) 0.0, y, 1);
      return y;
    }

    template <class T, class U>
    U& mult(const T& x, U& y) const
    {
      unsigned rows = x.rows();
      unsigned cols = x.columns();
      if (rows != this->columns())
        throw MissizeException();
      y.resize(this->rows(), cols);
      if (_rows == 0 || _columns == 0 || cols == 0)
      {
        y.set_zero();
        return y;
      }
      if (cols > 1)
        CBLAS::gemm(CblasNoTrans, CblasNoTrans, _rows, cols, _columns, *this, _rows, x, rows, (Real) 1.0, (Real) 0.0, y, _rows); 
      else
        CBLAS::gemv(CblasNoTrans, _rows, _columns, *this, _rows, x, 1, (Real) 1.0, (Real) 0.0, y, 1);
      return y;
    }

  protected:
    boost::shared_array<Real> _data;
    unsigned _rows;
    unsigned _columns;
    unsigned _capacity;
}; // end class

std::ostream& operator<<(std::ostream& out, const MatrixN& m);
std::istream& operator>>(std::istream& in, MatrixN& m);
// include inline functions
#include "MatrixN.inl"

} // end namespace

#endif

