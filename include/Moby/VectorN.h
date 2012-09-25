/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTORN_H
#define _VECTORN_H

#include <cstdarg>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <Moby/cblas.h>
#include <Moby/Types.h>
#include <Moby/Vector2.h>
#include <Moby/Vector3.h>
#include <Moby/SVector6.h>
#include <Moby/MissizeException.h>
#include <Moby/InvalidIndexException.h>

namespace Moby {

class MatrixN;
  
/// A generic N-dimensional floating point vector
class VectorN
{
  public:
    template <class ForwardIterator>
    VectorN(ForwardIterator begin, ForwardIterator end);

    VectorN();
    VectorN(unsigned N);
    VectorN(const VectorN& source);
    VectorN(VectorN& source);
    VectorN(const Vector2& v);
    VectorN(const Vector3& v);
    VectorN(const SVector6& v);
    VectorN(unsigned N, const Real* array);
    VectorN(unsigned N, boost::shared_array<Real> array);
    static VectorN construct_variable(unsigned N, ...);
    virtual ~VectorN() {}
    virtual Real dot(const VectorN& v) const { return dot(*this, v); }
    static Real dot(const VectorN& v1, const VectorN& v2);
    VectorN& normalize() { assert(norm() > 0.0); operator*=(1.0/norm()); return *this; }
    static Real calc_abs_err(const VectorN& vapprox, const VectorN& vtrue);
    static Real calc_rel_err(const VectorN& vapprox, const VectorN& vtrue);
    static VectorN normalize(const VectorN& v);
    unsigned size() const { return _len; }
    static Real norm1(const VectorN& v);
    static Real norm(const VectorN& v);
    static Real norm_sq(const VectorN& v) { return VectorN::dot(v, v); }
    Real norm_inf() const;
    Real norm1() const { return norm1(*this); }
    Real norm() const { return norm(*this); }
    Real norm_sq() const { return norm_sq(*this); }
    static VectorN one(unsigned N);
    MatrixN* outer_prod(const VectorN& v, MatrixN* m) const { return outer_prod(*this, v, m); }
    static MatrixN* outer_prod(const VectorN& v1, const VectorN& v2, MatrixN* m);
    static VectorN concat(const VectorN& v1, const VectorN& v2);
    static VectorN& concat(const VectorN& v1, const VectorN& v2, VectorN& result);
    VectorN& augment(const VectorN& v);
    VectorN& resize(unsigned N, bool preserve = false);
    VectorN& set_zero();
    VectorN& set_one();
    VectorN get_sub_vec(unsigned start_idx, unsigned end_idx) const;
    VectorN& get_sub_vec(unsigned start_idx, unsigned end_idx, VectorN& v) const;
    VectorN& set_sub_vec(unsigned start_idx, const VectorN& v);
    static VectorN zero(unsigned n);
    bool is_finite() const;
    VectorN& negate();
    VectorN& operator=(Real r);
    VectorN& operator=(const Vector3& source);
    VectorN& operator=(VectorN& source);
    VectorN& operator=(const VectorN& source);
    VectorN operator+(const VectorN& v) const;
    VectorN& operator+=(const VectorN& v);
    VectorN operator-(const VectorN& v) const;
    VectorN& operator-=(const VectorN& v);
    VectorN operator*(Real scalar) const;
    VectorN& operator*=(Real scalar);
    VectorN operator/(Real scalar) const;
    VectorN& operator/=(Real scalar) { return operator*=(1.0/scalar); }
    VectorN operator-() const;
    Real& operator[](const unsigned i) { assert(i < _len); return _data[i]; }
    Real operator[](const unsigned i) const { assert(i < _len); return _data[i]; }
    Real* data() { return _data.get(); }
    const Real* data() const { return _data.get(); }
    bool operator<(const VectorN& v) const;
    bool operator==(const VectorN& v) const;
    static bool epsilon_equals(const VectorN& v1, const VectorN& v2, Real epsilon);
    bool epsilon_equals(const VectorN& v, Real epsilon) const;
    static VectorN parse(const std::string& s);
    VectorN select(const std::vector<bool>& indices) const;
    VectorN& select(const std::vector<bool>& indices, VectorN& result) const;
    VectorN& resize(unsigned m, unsigned n) { assert(n == 1); resize(m, false); return *this; }

    template <class ForwardIterator>
    VectorN select(ForwardIterator idx_begin, ForwardIterator idx_end) const;

    template <class ForwardIterator>
    VectorN& select(ForwardIterator idx_begin, ForwardIterator idx_end, VectorN& v) const;

    template <class ForwardIterator>
    VectorN& set(ForwardIterator idx_begin, ForwardIterator idx_end, const VectorN& v);

    /// Sets this to a n-length zero vector
    VectorN& set_zero(unsigned n) { return resize(n).set_zero(); } 

    /// Sets this to a n-length ones vector
    VectorN& set_one(unsigned n) { return resize(n).set_one(); }

    /// Gets iterator to the start of the data
    Real* begin() { return (_len > 0) ? &_data[0] : NULL; }

    /// Gets constant iterator to the start of the data
    const Real* begin() const { return (_len > 0) ? &_data[0] : NULL; }

    /// Gets iterator to the end of the data
    Real* end() { return (_len > 0) ? &_data[0] + _len : NULL; }

    /// Gets constant iterator to the end of the data
    const Real* end() const { return (_len > 0) ? &_data[0] + _len : NULL; }

    /// Copies the source vector    
    VectorN& copy_from(const VectorN& source) { return operator=(source); }

    /// Swaps this vector with another
    void swap(VectorN& source) { std::swap(_data, source._data); std::swap(_len, source._len); std::swap(_capacity, source._capacity); }

    template <class V>
    V& get_sub_vec(unsigned start, unsigned end, V& v) const;

    template <class V>
    VectorN& set_sub_vec(unsigned start, const V& v);

    unsigned rows() const { return _len; }
    unsigned columns() const { return 1; }

  protected:
    static bool rel_equal(Real x, Real y);
    boost::shared_array<Real> _data;
    unsigned _len;
    unsigned _capacity;
}; // end class


std::ostream& operator<<(std::ostream& out, const VectorN& v);
std::istream& operator>>(std::istream& in, VectorN& v);

// include inline functions
#include "VectorN.inl"

} // end namespace

#endif
