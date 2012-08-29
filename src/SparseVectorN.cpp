/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/MissizeException.h>
#include <Moby/Constants.h>
#include <Moby/SparseVectorN.h>

using std::map;
using boost::shared_array;
using namespace Moby;

/// Creates a sparse vector from a dense vector
SparseVectorN::SparseVectorN(const VectorN& x)
{
  // setup the vector size
  _size = x.size();

  // determine the nonzero elements
  _nelm = 0;
  shared_array<bool> nz_elms(new bool[x.size()]);
  for (unsigned i=0; i< x.size(); i++)
    if (std::fabs(x[i]) < NEAR_ZERO)
    {
      _nelm++;
      nz_elms[i] = true;
    }
    else
      nz_elms[i] = false;

  // declare memory
  _indices = shared_array<unsigned>(new unsigned[_nelm]);
  _data = shared_array<Real>(new Real[_nelm]);

  // setup data
  for (unsigned i=0, j=0; i< x.size(); i++)
    if (nz_elms[i])
    {
      _indices[j] = i;
      _data[j] = x[i];
      j++;
    }
}

SparseVectorN::SparseVectorN(unsigned n, const map<unsigned, Real>& values)
{
  // declare memory
  _nelm = values.size();
  _size = n;
  _indices = shared_array<unsigned>(new unsigned[_nelm]);
  _data = shared_array<Real>(new Real[_nelm]);
  
  unsigned j=0;
  for (map<unsigned, Real>::const_iterator i = values.begin(); i != values.end(); i++)
  {
    _indices[j] = i->first;
    _data[j] = i->second;
    j++;
  }
}

SparseVectorN::SparseVectorN(unsigned n, unsigned nelms, shared_array<unsigned> indices, shared_array<Real> data)
{
  _size = n;
  _nelm = nelms;
  _indices = indices;
  _data = data;
}

/// Computes the dot product between a sparse vector and itself
Real SparseVectorN::square() const
{
  Real result = 0;

  for (unsigned i=0; i< _nelm; i++)
    result += _data[i] * _data[i];

  return result;
}

/// Computes the dot product between a sparse vector and a dense vector
Real SparseVectorN::dot(const VectorN& x) const
{
  Real result = 0;

  if (x.size() != _size)
    throw MissizeException();

  for (unsigned i=0; i< _nelm; i++)
    result += x[_indices[i]] * _data[i];

  return result;
}

/// Gets the dense version of the vector 
VectorN& SparseVectorN::to_dense(VectorN& result) const
{
  result.set_zero(_size);
  for (unsigned i=0; i< _nelm; i++)
    result[_indices[i]] = _data[i];

  return result;
}

/// Multiplies this sparse vector by a scalar (in place)
SparseVectorN& SparseVectorN::operator*=(Real scalar)
{
  for (unsigned i=0; i< _nelm; i++)
    _data[i] *= scalar;
  return *this;
}

/// Multiplies this sparse vector by a scalar and returns the result in a new vector
SparseVectorN& SparseVectorN::mult(Real scalar, SparseVectorN& result) const
{
  result._size = this->_size;
  result._nelm = this->_nelm;
  result._indices = shared_array<unsigned>(new unsigned[this->_nelm]);
  result._data = shared_array<Real>(new Real[this->_nelm]);
  for (unsigned i=0; i< this->_nelm; i++)
  {
    result._indices[i] = this->_indices[i];
    result._data[i] = this->_data[i] * scalar;
  }

  return result;
}

/// Negates this sparse vector and returns the result in a new vector
SparseVectorN SparseVectorN::operator-() const
{
  SparseVectorN s;
  s._size = this->_size;
  s._nelm = this->_nelm;
  s._indices = shared_array<unsigned>(new unsigned[this->_nelm]);
  s._data = shared_array<Real>(new Real[this->_nelm]);
  for (unsigned i=0; i< this->_nelm; i++)
  {
    s._indices[i] = this->_indices[i];
    s._data[i] = -this->_data[i];
  }
  return s;
}

