/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <Moby/FastThreadable.h>
#include <Moby/cblas.h>
#include <Moby/MissizeException.h>
#include <Moby/Constants.h>
#include <Moby/SpatialABInertia.h>
#include <Moby/SMatrix6N.h>

using namespace Moby;

/// Default constructor - constructs an empty matrix
SMatrix6N::SMatrix6N() : MatrixN()
{
}

/// Constructs a 6 x columns dimensional matrix 
/**
 * Sets matrix to the 6 x columns dimensional zero matrix
 */
SMatrix6N::SMatrix6N(unsigned columns) : MatrixN(6, columns)
{
}

/// Constructs a matrix from an array
/**
 * \param rows the number of rows of the matrix
 * \param columns the number of columns of the matrix
 * \param array an array of rows*columns Real values in row-major format
 */
SMatrix6N::SMatrix6N(unsigned columns, const Real* array) : MatrixN(6, columns, array)
{
}

/// Copy constructor
SMatrix6N::SMatrix6N(const MatrixN& source) : MatrixN(source) 
{
}

/// Assignment constructor
/**
 * \note source is destroyed on return
 */
SMatrix6N::SMatrix6N(MatrixN& source) : MatrixN(source) 
{
  if (source.rows() != 6)
    throw std::runtime_error("SMatrix6N::SMatrix6N() - source matrix does not have six rows!");
}

/// Copy constructor
SMatrix6N::SMatrix6N(const SMatrix6N& source) : MatrixN(source) 
{
}

/// Assignment constructor
/**
 * \note source is destroyed on return
 */
SMatrix6N::SMatrix6N(SMatrix6N& source) : MatrixN(source) 
{
}

/// Gets the desired column of this matrix
SVector6 SMatrix6N::get_column(unsigned i) const
{
  const unsigned SPATIAL_DIM = 6;

  // make sure i is valid
  if (i >= columns())
    throw InvalidIndexException();

  // get the vector
  SVector6 v;
  for (unsigned j=0, k=SPATIAL_DIM*i; j< SPATIAL_DIM; j++, k++)
    v[j] = _data[k];

  return v;
}

/// Sets the desired column of this matrix
void SMatrix6N::set_column(unsigned i, const SVector6& v)
{
  const unsigned SPATIAL_DIM = 6;

  // make sure i is valid
  if (i >= columns())
    throw InvalidIndexException();

  // set the values
  for (unsigned j=0, k=SPATIAL_DIM*i; j< SPATIAL_DIM; j++, k++)
    _data[k] = v[j];
}

/// Returns the transpose of the given matrix
MatrixN& SMatrix6N::transpose(const SMatrix6N& m, MatrixN& result)
{
  unsigned X = 0, Z = 2, A = 3, C = 5;

  // check for empty matrix
  if (m.columns() == 0)
    return result.set_zero(0, 6);

  // setup the transposed matrix
  result.resize(m.columns(), m.rows());

  // handle the top part first
  SAFESTATIC FastThreadable<MatrixN> work, workT;
  m.get_sub_mat(X, Z+1, 0, m.columns(), work());
  MatrixN::transpose(work(), workT());
  result.set_sub_mat(0, A, workT());

  // handle the bottom part now 
  m.get_sub_mat(A, C+1, 0, m.columns(), work());
  MatrixN::transpose(work(), workT());
  result.set_sub_mat(0, X, workT());

  return result;
}

/// Returns the transpose of this matrix
MatrixN SMatrix6N::transpose() const
{  
  MatrixN result;
  SMatrix6N::transpose(*this, result);
  return result;
}

/// Resizes this matrix
MatrixN& SMatrix6N::resize(unsigned rows, unsigned columns, bool preserve)
{
  if (rows != 6)
    throw std::runtime_error("SMatrix6N::resize() - attempt to resize spatial matrix rows when != 6");

  return MatrixN::resize(rows, columns, preserve);
}

/// Multiplies a standard matrix by the transpose of a spatial matrix
MatrixN& SMatrix6N::mult_transpose(const MatrixN& m1, const SMatrix6N& m2, MatrixN& result)
{
  if (m1.columns() != m2.columns())
    throw MissizeException();

  // resize the result matrix
  result.resize(m1.rows(), m2.rows());

  // look for empty result
  if (m1.rows() == 0 || m2.rows() == 0)
  {
    result.set_zero();
    return result;
  }

  // setup parameters for GEMM
  const unsigned M = m1.rows();
  const unsigned N = 3;
  const unsigned K = m1.columns();
  const Real ALPHA = 1.0;
  const Real BETA = 0.0;
  const unsigned LDA = M;
  const unsigned LDB = 6;
  const unsigned LDC = M;

  // we first want to multiply m1 by the bottom half of m2
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, M, N, K, ALPHA, m1.begin(), LDA, m2.begin()+N, LDB, BETA, result.begin(), LDC);

  // now multiply m1 by the top half of m2
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, M, N, K, ALPHA, m1.begin(), LDA, m2.begin(), LDB, BETA, result.begin()+result.rows()*N, LDC);

  return result;
}

/// Multiplies this spatial matrix by a spatial vector
VectorN& SMatrix6N::mult(const SVector6& v, VectorN& result) const
{
  if (_columns != v.size())
    throw MissizeException();

  // resize the result vector
  result.resize(_rows);

  // look for empty result
  if (_columns == 0)
  {
    result.set_zero();
    return result;
  }

  // do the multiplication
  CBLAS::gemv(CblasNoTrans, *this, v, (Real) 1.0, (Real) 0.0, result); 

  return result;
}

/// Multiplies the transpose of this matrix by a non-spatial vector
VectorN& SMatrix6N::transpose_mult(const VectorN& v, VectorN& result) const
{
  if (_rows != v.size())
    throw MissizeException();

  SVector6 v6(v.begin());
  return transpose_mult(v6, result);
}

/// Multiplies the transpose of this matrix by a non-spatial vector
VectorN SMatrix6N::transpose_mult(const SVector6& v) const
{
  VectorN result;
  transpose_mult(v, result);
  return result;
}

/// Multiplies the transpose of this matrix by a spatial vector
VectorN& SMatrix6N::transpose_mult(const SVector6& v, VectorN& result) const
{
  // resize the result
  result.resize(_columns);

  // look for empty result
  if (_columns == 0)
    return result;

  // setup parameters for GEMV
  const unsigned M = 3;
  const unsigned N = _columns;
  const Real ALPHA = 1.0;
  const Real BETA1 = (Real) 0.0;
  const Real BETA2 = (Real) 1.0;
  const unsigned LDA = 6;
  const unsigned INCX = 1, INCY = 1;

  // we first want to multiply the bottom half of this by v
  CBLAS::gemv(CblasColMajor, CblasTrans, M, N, ALPHA, begin()+M, LDA, v.begin(), INCX, BETA1, result.begin(), INCY);

  // now multiply the top half of this by v
  CBLAS::gemv(CblasColMajor, CblasTrans, M, N, ALPHA, begin(), LDA, v.begin()+M, INCX, BETA2, result.begin(), INCY);

  return result;
}

/// Multiplies the transpose of this matrix by a spatial rigid body inertia matrix
MatrixN& SMatrix6N::transpose_mult(const SpatialRBInertia* I, MatrixN& result) const
{
  const unsigned SPATIAL_DIM = 6;

  if (_rows != SPATIAL_DIM)
    throw MissizeException();

  // resize the result matrix
  result.resize(_columns, SPATIAL_DIM);

  // look for empty result
  if (_columns == 0)
  {
    result.set_zero();
    return result;
  }

  // get the data from result
  Real* rdata = result.data();
  const unsigned LDR = result.rows();

  // loop over all columns
  for (unsigned i=0, j=0; i< _columns; i++, j+= SPATIAL_DIM)
  {
    // get the vectors corresponding to this column
    Vector3 bottom(_data[j], _data[j+1], _data[j+2]);
    Vector3 top(_data[j+3], _data[j+4], _data[j+5]);

    // do the computation
    Vector3 h_x_bottom = Vector3::cross(I->h, bottom);
    Vector3 top_x_h = Vector3::cross(top, I->h);
    Vector3 left = h_x_bottom + I->J*top;
    Vector3 right = bottom*I->m + top_x_h;
    
    unsigned kk=i;
    rdata[kk] = left[0]; kk += LDR;
    rdata[kk] = left[1]; kk += LDR;
    rdata[kk] = left[2]; kk += LDR;
    rdata[kk] = right[0]; kk += LDR;
    rdata[kk] = right[1]; kk += LDR;
    rdata[kk] = right[2];
  }

  return result;
}
 
/// Multiplies the transpose of this matrix by a spatial articulated body inertia matrix
MatrixN& SMatrix6N::transpose_mult(const SpatialABInertia* I, MatrixN& result) const
{
  const unsigned SPATIAL_DIM = 6;

  if (_rows != SPATIAL_DIM)
    throw MissizeException();

  // resize the result matrix
  result.resize(_columns, SPATIAL_DIM);

  // look for empty result
  if (_columns == 0)
  {
    result.set_zero();
    return result;
  }

  // get the transpose of all matrices in I
  const Matrix3 HT = Matrix3::transpose(I->H);

  // get the data from result
  Real* rdata = result.data();
  const unsigned LDR = result.rows();

  // loop over all columns
  for (unsigned i=0, j=0; i< _columns; i++, j+= SPATIAL_DIM)
  {
    // get the vectors corresponding to this column
    Vector3 top(_data[j], _data[j+1], _data[j+2]);
    Vector3 bottom(_data[j+3], _data[j+4], _data[j+5]);

    // do the computation
    Vector3 left = I->H*bottom + I->J*top;
    Vector3 right = I->M*bottom + HT*top;
    
    unsigned kk=i;
    rdata[kk] = left[0]; kk += LDR;
    rdata[kk] = left[1]; kk += LDR;
    rdata[kk] = left[2]; kk += LDR;
    rdata[kk] = right[0]; kk += LDR;
    rdata[kk] = right[1]; kk += LDR;
    rdata[kk] = right[2];
  }

  return result;
}
 
/// Multiplies the transpose of this matrix by a non-spatial matrix
MatrixN& SMatrix6N::transpose_mult(const MatrixN& m, MatrixN& result) const
{
  if (_rows != m.rows())
    throw MissizeException();

  // resize the result matrix
  result.resize(_columns, m.columns());

  // look for empty result
  if (_columns == 0 || m.columns() == 0)
    return result.set_zero();

  // setup parameters for GEMM
  const unsigned M = _columns;
  const unsigned N = m.columns();
  const unsigned K = 3;
  const Real ALPHA = 1.0;
  const Real BETA1 = 0.0;
  const Real BETA2 = 1.0;
  const unsigned LDA = 6;
  const unsigned LDB = 6;
  const unsigned LDC = M;

  // we first want to multiply the bottom half of this by the top half of m
  CBLAS::gemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, ALPHA, begin()+K, LDA, m.begin(), LDB, BETA1, result.begin(), LDC);

  // now multiply the top half of this by the bottom half of m, making sure
  // to add the result (this is why BETA = 1.0)
  CBLAS::gemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, ALPHA, begin(), LDA, m.begin()+K, LDB, BETA2, result.begin(), LDC);

  return result;
}

/// Multiplies the transpose of this matrix by the transpose of a non-spatial matrix
MatrixN& SMatrix6N::transpose_mult_transpose(const MatrixN& m, MatrixN& result) const
{
  if (_rows != m.columns())
    throw MissizeException();

  // resize the result matrix
  result.resize(_columns, m.rows());

  // look for empty result
  if (_columns == 0 || m.rows() == 0)
  {
    result.set_zero();
    return result;
  }

  // setup parameters for GEMM
  const unsigned M = _columns;
  const unsigned N = m.rows();
  const unsigned K = 3;
  const Real ALPHA = 1.0;
  const Real BETA1 = 0.0;
  const Real BETA2 = 1.0;
  const unsigned LDA = 6;
  const unsigned LDB = N;
  const unsigned LDC = M;

  // we first want to multiply the bottom half of this by the left half of m
  CBLAS::gemm(CblasColMajor, CblasTrans, CblasTrans, M, N, K, ALPHA, begin()+K, LDA, m.begin(), LDB, BETA1, result.begin(), LDC);

  // now multiply the top half of this by the right half of m, making sure
  // to add the result (this is why BETA2 = 1.0)
  CBLAS::gemm(CblasColMajor, CblasTrans, CblasTrans, M, N, K, ALPHA, begin(), LDA, m.begin()+K*m.rows(), LDB, BETA2, result.begin(), LDC);

  return result;
}

/*
/// Multiplies a 6x6 spatial matrix by a 6xn spatial matrix
SMatrix6N& SMatrix6N::mult(const SMatrix6& m1, const SMatrix6N& m2, SMatrix6N& result)
{
  if (m2.rows() != m1.columns())
    throw MissizeException();

  // resize the result
  result.resize(m1.rows(), m2.columns());

  // look for zero result
  if (m1.rows() == 0 || m1.columns() == 0 || m2.columns() == 0)
  {
    result.set_zero();
    return result;
  } 

  // carry out multiplication
  CBLAS::gemm(CblasNoTrans, CblasNoTrans, m1, m2, (Real) 1.0, (Real) 0.0, result); 
  return result;
}
*/

/// Multiplies a 6xN matrix by a Nx1 vector and returns a 6-dimensional vector
SVector6 SMatrix6N::mult(const VectorN& v) const
{
  if (columns() != v.size())
    throw MissizeException();

  // look for empty result
  if (_columns == 0)
    return ZEROS_6;

  // setup the result vector
  SVector6 result;

  // carry out multiplication
  CBLAS::gemv(CblasNoTrans, *this, v, (Real) 1.0, (Real) 0.0, result);

  return result;
}

