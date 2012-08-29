/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_VECTOR_N_H_
#define _SPARSE_VECTOR_N_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/VectorN.h>

namespace Moby {

/// A sparse vector represented in 'CSR' format
class SparseVectorN
{
  public:
    SparseVectorN() { _size = _nelm = 0; }
    SparseVectorN(unsigned n, const std::map<unsigned, Real>& values);
    SparseVectorN(unsigned n, unsigned nnz, boost::shared_array<unsigned> indices, boost::shared_array<Real> data);
    SparseVectorN(const VectorN& v);
    Real dot(const VectorN& x) const;
    Real square() const;
    unsigned size() const { return _size; }
    unsigned num_elements() const { return _nelm; }
    boost::shared_array<unsigned> get_indices() const { return _indices; }
    boost::shared_array<Real> get_data() const { return _data; }
    VectorN& to_dense(VectorN& result) const;
    VectorN to_dense() const { VectorN result; to_dense(result); return result; }
    SparseVectorN operator-() const;
    SparseVectorN& operator*=(Real scalar);
    SparseVectorN& mult(Real scalar, SparseVectorN& result) const;

  protected:
    boost::shared_array<unsigned> _indices;   // indices of the data
    boost::shared_array<Real> _data;          // 
    unsigned _size;
    unsigned _nelm;
}; // end class

std::ostream& operator<<(std::ostream& out, const SparseVectorN& s);

} // end namespace

#endif

