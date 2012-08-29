/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_MATRIX_N_H_
#define _SPARSE_MATRIX_N_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/SparseVectorN.h>
#include <Moby/VectorN.h>

namespace Moby {

/// A sparse matrix represented in 'CSR' format
class SparseMatrixN
{
  public:
    SparseMatrixN();
    SparseMatrixN(unsigned m, unsigned n, const std::map<std::pair<unsigned, unsigned>, Real>& values);
    SparseMatrixN(unsigned m, unsigned n, boost::shared_array<unsigned> ptr, boost::shared_array<unsigned> indices, boost::shared_array<Real> data);
    SparseMatrixN(const MatrixN& m);
    static SparseMatrixN identity(unsigned n);
    VectorN& mult(const VectorN& x, VectorN& result) const;
    VectorN& transpose_mult(const VectorN& x, VectorN& result) const;
    MatrixN& mult(const MatrixN& m, MatrixN& result) const;
    MatrixN mult(const MatrixN& m) const;
    MatrixN& mult_transpose(const MatrixN& m, MatrixN& result) const;
    MatrixN mult_transpose(const MatrixN& m) const;
    MatrixN& transpose_mult(const MatrixN& m, MatrixN& result) const;
    MatrixN transpose_mult(const MatrixN& m) const;
    MatrixN& transpose_mult_transpose(const MatrixN& m, MatrixN& result) const;
    MatrixN transpose_mult_transpose(const MatrixN& m) const;
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    SparseMatrixN get_sub_mat(unsigned rstart, unsigned rend, unsigned cstart, unsigned cend) const;
    SparseVectorN& get_row(unsigned i, SparseVectorN& row) const;
    SparseVectorN get_row(unsigned i) const;
    SparseVectorN& get_column(unsigned i, SparseVectorN& column) const;
    SparseVectorN get_column(unsigned i) const;
    boost::shared_array<unsigned> get_indices() const { return _indices; }
    boost::shared_array<unsigned> get_ptr() const { return _ptr; }
    boost::shared_array<Real> get_data() const { return _data; }
    SparseMatrixN& operator-=(const SparseMatrixN& m);
    SparseMatrixN& operator+=(const SparseMatrixN& m);
    SparseMatrixN& operator*=(Real scalar);
    SparseMatrixN& negate();
    static SparseMatrixN& outer_square(const VectorN& g, SparseMatrixN& result);
    static SparseMatrixN& outer_square(const SparseVectorN& v, SparseMatrixN& result);
    MatrixN to_dense() const;
    MatrixN& to_dense(MatrixN& m) const;

  protected:
    SparseMatrixN(unsigned m, unsigned n) { _rows = m; _columns = n; }
    boost::shared_array<unsigned> _indices;  // column indices
    boost::shared_array<unsigned> _ptr;      // starting indices for row i 
    boost::shared_array<Real> _data;         // actual data
    unsigned _rows;
    unsigned _columns;

  private:
    void set(unsigned rows, unsigned columns, const std::map<std::pair<unsigned, unsigned>, Real>& values);
}; // end class

std::ostream& operator<<(std::ostream& out, const SparseMatrixN& s);

} // end namespace

#endif

