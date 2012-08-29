/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <cmath>
#include <Moby/LinAlg.h>
#include <Moby/MissizeException.h>
#include <Moby/MatrixNN.h>

using namespace Moby;
using boost::shared_array;
using std::vector;

/// Default constructor - constructs an empty matrix
MatrixNN::MatrixNN() : MatrixN()
{
}

/// Constructs a MatrixNN object from a MatrixN object
MatrixNN::MatrixNN(const MatrixN& source) : MatrixN(source)
{
  assert (source.rows() == source.columns());
}

/// Constructs a dim*dim matrix 
MatrixNN::MatrixNN(unsigned dim) : MatrixN(dim, dim)
{
}

/// Constructs a matrix from an array
/**
 * \param the dimension of the matrix rows and columns
 * \param array an array of dim*dim Real values in row-major format
 */
MatrixNN::MatrixNN(unsigned dim, const Real* array) : MatrixN(dim, dim, array)
{
}

/// Constructs a matrix from an array
/**
 * \param size the number of size (rows/columns) of the matrix
 * \param array an array of rows*columns Real values in row-major format;
 * \note MatrixNN will modify the contents of array until the MatrixNN object is
 *       destructed
 */
MatrixNN::MatrixNN(unsigned size, shared_array<Real> array) : MatrixN(size, size, array)
{
}

/// Copy constructor
MatrixNN::MatrixNN(const MatrixNN& source) : MatrixN(source)
{
}

/// Sets this equal to the MatrixN source
/**
 * \throws MissizeException if source is not a square matrix
 */
MatrixNN& MatrixNN::operator=(const MatrixN& source) 
{ 
  if (source.rows() != source.columns())
    throw MissizeException();
  MatrixN::operator=(source);
  return *this;
}

/// Multiplies a 2x2 MatrixNN by a vector
/**
 * \throws MissizeException if this is not a 2x2 matrix
 */
Vector2 MatrixNN::operator*(const Vector2& v) const 
{ 
  if (size() != v.size())
    throw MissizeException();
  return Vector2(_data[0]*v[0] + _data[2]*v[1], _data[1]*v[0] + _data[3]*v[1]); 
}

/// Multiplies a 3x3 MatrixNN by a vector
/**
 * \throws MissizeException if this is not a 2x2 matrix
 */
Vector3 MatrixNN::operator*(const Vector3& v) const 
{ 
  if (size() != v.size())
    throw MissizeException();
  return Vector3(_data[0]*v[0] + _data[3]*v[1] + _data[6]*v[2], _data[1]*v[0] + _data[4]*v[1] + _data[7]*v[2], _data[2]*v[0] + _data[5]*v[1] + _data[8]*v[2]); 
}

/// Checks whether the given matrix is symmetric to the specified tolerance
bool MatrixNN::is_symmetric(Real tolerance) const
{
  // make sure that tolerance is non-negative
  tolerance = std::fabs(tolerance);

  // check symmetry
  for (unsigned i=0; i< _rows; i++)
    for (unsigned j=i+1; j< _columns; j++)
      if (std::fabs(operator()(i,j) - operator()(j,i)) > tolerance)
        return false;

  return true;
}

/// Sets this matrix to its inverse
MatrixNN& MatrixNN::inverse()
{
  LinAlg::inverse(*this);
  return *this;
}

/// Determines the inverse of the given matrix
MatrixNN MatrixNN::inverse(const MatrixNN& m)
{
  MatrixNN n(m);
  n.inverse();
  return n;  
}

/// Returns the identity matrix
/**
 * \param n the number of rows / columns of the matrix
 * \return a NxN identity matrix
 */
MatrixNN MatrixNN::identity(unsigned dim)
{
  MatrixNN m(dim);
  m.set_identity();
  
  return m;
}

/// Sets this to the identity matrix 
MatrixNN& MatrixNN::set_identity()
{
  set_zero();
  for (unsigned i=0; i< _rows; i++)
    _data[i*_rows+i] = 1.0;
  return *this;
}

/// Returns the zero matrix
/**
 * \param n the number of rows / columns of the matrix
 * \return a NxN zero matrix
 */
MatrixNN MatrixNN::zero(unsigned n)
{
  MatrixNN m(n);
  m.set_zero();

  return m;
}

/// Resizes this matrix
/**
 * \note memory is not reallocated, unless the matrix is growing in size
 */
MatrixN& MatrixNN::resize(unsigned m, unsigned n, bool preserve_data)
{
  if (m != n)
    throw std::runtime_error("Cannot resize a square matrix to non-square!");

  return MatrixN::resize(m, n, preserve_data);
}

/// Multiplies the diagonal matrix formed from d by the matrix m
MatrixN MatrixNN::diag_mult(const VectorN& d, const MatrixN& m)
{
  MatrixN result;
  diag_mult(d, m, result);
  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix transpose(m)
MatrixN MatrixNN::diag_mult_transpose(const VectorN& d, const MatrixN& m)
{
  MatrixN result;
  diag_mult_transpose(d, m, result);
  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix m
MatrixN& MatrixNN::diag_mult(const VectorN& d, const MatrixN& m, MatrixN& result)
{
  if (d.size() != m.rows())
    throw MissizeException();

  result.resize(d.size(), m.columns());
  for (unsigned i=0; i< m.columns(); i++)
    std::transform(d.begin(), d.end(), m.begin()+m.rows()*i, result.begin()+result.rows()*i, std::multiplies<Real>());

  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix transpose(m)
MatrixN& MatrixNN::diag_mult_transpose(const VectorN& d, const MatrixN& m, MatrixN& result)
{
  if (d.size() != m.columns())
    throw MissizeException();

  // copy the transpose of m to the result
  MatrixN::transpose(m, result);

  // do the specified number of times
  for (unsigned i=0; i< m.rows(); i++)
    CBLAS::scal(d.size(), d[i], result.begin()+i, result.rows());

  return result;
}

/// Multiplies the diagonal matrix formed from d by the vector v
VectorN& MatrixNN::diag_mult(const VectorN& d, const VectorN& v, VectorN& result)
{
  if (d.size() != v.size())
    throw MissizeException();

  result.resize(d.size());
  std::transform(d.begin(), d.end(), v.begin(), result.begin(), std::multiplies<Real>());
  return result;
}

/// Multiplies the diagonal matrix formed from d by the vector v
VectorN MatrixNN::diag_mult(const VectorN& d, const VectorN& v)
{
  VectorN result;
  return diag_mult(d, v, result);
}

/// Selects a submatrix from this
MatrixNN MatrixNN::select(const vector<bool>& indices) const
{
  MatrixNN result;
  return select(indices, result);
}

/// Selects a submatrix from this
MatrixNN& MatrixNN::select(const vector<bool>& indices, MatrixNN& result) const
{
  const unsigned n = size();
  assert(n == indices.size());

  // determine how many indices are selected and resize result vector
  unsigned nselect = 0;
  for (unsigned i=0; i< n; i++)
    if (indices[i])
      nselect++;
  result.resize(nselect);

  // get the source and target data
  const Real* source = data();
  Real* target = result.data();
 
  // copy 
  for (unsigned i=0, s=0, t=0; i< n; i++)
    for (unsigned j=0; j< n; j++, s++)
      if (indices[i] && indices[j])
        target[t++] = source[s];

  return result;
}
