/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SMATRIX6N_H
#define _SMATRIX6N_H

#include <Moby/MatrixN.h>

namespace Moby {

class SpatialABInertia;
class SpatialRBInertia;
        
/// A generic, possibly non-square 6xn spatial algebra matrix 
/**
 * The underlying data is stored in a column-major (e.g., the element at row 1, column 0 is element "1").
 */
class SMatrix6N : public MatrixN 
{
  public:
    SMatrix6N();
    SMatrix6N(unsigned columns);
    SMatrix6N(const MatrixN& source);
    SMatrix6N(MatrixN& source);
    SMatrix6N(const SMatrix6N& source);
    SMatrix6N(SMatrix6N& source);
    SMatrix6N(unsigned columns, const Real* array);
    MatrixN transpose() const;
    MatrixN& mult(const MatrixN& m, SMatrix6N& result) const { return MatrixN::mult(*this, m, result); }
    static MatrixN& transpose(const SMatrix6N& m, MatrixN& result);
    SMatrix6N& copy_from(const SMatrix6N& source) { MatrixN::operator=(source); return *this; }
    void swap(SMatrix6N& source) { std::swap(_data, source._data); }
    SMatrix6N& operator=(const SMatrix6N& source) { MatrixN::operator=(source); return *this; }
    SMatrix6N& operator=(SMatrix6N& source) { MatrixN::operator=(source); return *this; }
    SVector6 mult(const VectorN& v) const; 
    VectorN& mult(const SVector6& v, VectorN& result) const;
    virtual MatrixN& resize(unsigned rows, unsigned columns, bool preserve = false);
    virtual MatrixN& transpose_mult(const MatrixN& m, MatrixN& result) const;
    virtual MatrixN& transpose_mult(const SpatialRBInertia* m, MatrixN& result) const;
    virtual MatrixN& transpose_mult(const SpatialABInertia* m, MatrixN& result) const;
    virtual MatrixN& transpose_mult_transpose(const MatrixN& m, MatrixN& result) const;
    virtual VectorN& transpose_mult(const VectorN& v, VectorN& result) const;
//    static SMatrix6N& mult(const SMatrix6& m1, const SMatrix6N& m2, SMatrix6N& result);
    static MatrixN& mult_transpose(const MatrixN& m1, const SMatrix6N& m2, MatrixN& result);
    SVector6 get_column(unsigned i) const;
    void set_column(unsigned i, const SVector6& v);
    VectorN& transpose_mult(const SVector6& v, VectorN& result) const;
    VectorN transpose_mult(const SVector6& v) const;
    SMatrix6N operator*(const MatrixN& m) const { return SMatrix6N(MatrixN::operator*(m)); }
    SVector6 operator*(const VectorN& v) const { return mult(v); }
}; // end class

/*
inline SMatrix6N operator*(const SMatrix6& m1, const SMatrix6N& m2)
{
  SMatrix6N result;
  SMatrix6N::mult(m1, m2, result);
  return result;
}
*/

} // end namespace

#endif

