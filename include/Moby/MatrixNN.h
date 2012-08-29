/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIXNN_H
#define _MATRIXNN_H

#include <Moby/VectorN.h>
#include <Moby/MatrixN.h>

namespace Moby {
        
/// A generic square matrix
class MatrixNN : public MatrixN
{
  public:
    MatrixNN();
    MatrixNN(unsigned dim);
    MatrixNN(const MatrixN& source);
    MatrixNN(const MatrixNN& source);
    MatrixNN(unsigned int dim, const Real* array);
    MatrixNN(unsigned int dim, boost::shared_array<Real> array);
    MatrixNN& inverse();
    virtual ~MatrixNN() {}
    bool is_symmetric(Real tolerance) const;
    static MatrixNN inverse(const MatrixNN& m);
    static MatrixNN zero(unsigned dim);
    static MatrixNN identity(unsigned dim);
    MatrixNN& set_identity();
    virtual MatrixN& resize(unsigned m, unsigned n, bool preserve_data = false);
    MatrixNN& resize(unsigned n, bool preserve_data = false) { resize(n,n,preserve_data); return *this; }
    virtual MatrixNN& operator=(const MatrixN& source);
    virtual MatrixNN& operator=(const MatrixNN& source) { MatrixN::operator=(source); return *this; }
    void swap(MatrixNN& source) { std::swap(_data, source._data); std::swap(_rows, source._rows); std::swap(_columns, source._columns); }
    static VectorN diag_mult(const VectorN& d, const VectorN& v);
    static MatrixN diag_mult(const VectorN& d, const MatrixN& m);
    static MatrixN diag_mult_transpose(const VectorN& d, const MatrixN& m);
    static VectorN& diag_mult(const VectorN& d, const VectorN& v, VectorN& result);
    static MatrixN& diag_mult(const VectorN& d, const MatrixN& m, MatrixN& result);
    static MatrixN& diag_mult_transpose(const VectorN& d, const MatrixN& m, MatrixN& result);
    MatrixNN operator*(Real scalar) const { return MatrixNN(MatrixN::operator*(scalar)); }
    MatrixN operator*(const MatrixN& m) const { return MatrixN::operator*(m); }
    VectorN operator*(const VectorN& v) const { return MatrixN::operator*(v); }
    Vector2 operator*(const Vector2& v) const;
    Vector3 operator*(const Vector3& v) const;
    MatrixNN& set_zero() { MatrixN::set_zero(); return *this; }
    MatrixNN select(const std::vector<bool>& indices) const;
    MatrixNN& select(const std::vector<bool>& indices, MatrixNN& result) const;

    /// Sets this to a n-sized zero matrix
    MatrixNN& set_zero(unsigned n) { return resize(n).set_zero(); }

    /// Sets this to a n-sized identity matrix
    MatrixNN& set_identity(unsigned n) { return resize(n).set_identity(); }

    /// Copies FROM the source matrix
    MatrixNN& copy_from(const MatrixNN& source) { operator=(source); return *this; }

    /// Returns the number of rows / columns of this matrix
    unsigned size() const { return _rows; }

    template <class ForwardIterator>
    MatrixNN select(ForwardIterator start, ForwardIterator end) const;

    template <class ForwardIterator>
    MatrixNN& select(ForwardIterator start, ForwardIterator end, MatrixNN& m) const;

    template <class ForwardIterator1, class ForwardIterator2>
    MatrixN& select(ForwardIterator1 rstart, ForwardIterator1 rend, ForwardIterator2 cstart, ForwardIterator2 cend, MatrixN& m) const { return MatrixN::select(rstart, rend, cstart, cend, m); }

    template <class ForwardIterator1, class ForwardIterator2>
    MatrixN select(ForwardIterator1 rstart, ForwardIterator1 rend, ForwardIterator2 cstart, ForwardIterator2 cend) const { return MatrixN::select(rstart, rend, cstart, cend); }
}; // end class

// include inline functions
#include "MatrixNN.inl"

} // end namespace

#endif

