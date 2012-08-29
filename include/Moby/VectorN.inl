/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

inline VectorN operator*(Real scalar, const VectorN& v) { return v * scalar; }

/// Gets a sub-vector from this vector
template <class V>
V& VectorN::get_sub_vec(unsigned start, unsigned end, V& v) const
{
  if (start > end || end > _len)
    throw InvalidIndexException();

  // resize the vector
  const unsigned SZ = end - start;
  v.resize(SZ);

  // see whether to exit now
  if (_len == 0 || SZ == 0)
    return v;

  // copy using BLAS
  CBLAS::copy(SZ, _data.get()+start, 1, v.data(), 1);

  return v;
}

/// Sets a sub-vector of this vector
template <class V>
VectorN& VectorN::set_sub_vec(unsigned start, const V& v)
{
  // get size of v
  const unsigned SZ = v.size();

  // see whether index is valid
  if (start + SZ > _len)
    throw InvalidIndexException();

  // check whether to exit now
  if (SZ == 0)
    return *this;

  // copy using BLAS
  CBLAS::copy(SZ, v.data(), 1, _data.get()+start, 1);

  return *this;
}

/// Sets a vector from iterators
template <class ForwardIterator>
VectorN::VectorN(ForwardIterator begin, ForwardIterator end)
{
  // setup the vector
  _len = std::distance(begin, end);
  _capacity = _len;
  _data = boost::shared_array<Real>(new Real[_len]);

  // copy the elements
  std::copy(begin, end, this->begin());
}

/// Gets a subvector (not necessarily contiguous)
template <class ForwardIterator>
VectorN& VectorN::select(ForwardIterator idx_begin, ForwardIterator idx_end, VectorN& v) const
{
  // setup the vector
  unsigned len = std::distance(idx_begin, idx_end);
  if (len > _len)
    throw MissizeException();
  v.resize(len);

  // copy
  unsigned vidx = 0;
  for (ForwardIterator i = idx_begin; i != idx_end; i++)
    v[vidx++] = _data[*i];

  return v;
}

/// Gets a subvector (not necessarily contiguous)
template <class ForwardIterator>
VectorN VectorN::select(ForwardIterator idx_begin, ForwardIterator idx_end) const
{
  VectorN v;
  return select(idx_begin, idx_end, v);
}  

/// Sets a subvector; other components are unchanged
template <class ForwardIterator>
VectorN& VectorN::set(ForwardIterator idx_begin, ForwardIterator idx_end, const VectorN& v)
{
  unsigned len = std::distance(idx_begin, idx_end);
  if (len != v.size() || len > _len)
    throw MissizeException();

  // copy
  unsigned vidx = 0;
  for (ForwardIterator i = idx_begin; i != idx_end; i++)
    _data[*i] = v[vidx++];

  return *this;
}

