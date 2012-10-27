/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <cmath>
#include <iomanip>
#include <Moby/FastThreadable.h>
#include <Moby/MissizeException.h>
#include <Moby/MatrixN.h>

using std::vector;
using namespace Moby;
using boost::shared_array;

/// Default constructor - constructs an empty matrix
MatrixN::MatrixN()
{
  _rows = _columns = 0;
  _capacity = 0;
}

/// Constructs a rows x columns dimensional (unitialized) matrix 
MatrixN::MatrixN(unsigned rows, unsigned columns)
{
  _rows = rows;
  _columns = columns;
  _capacity = rows*columns;
  if (rows > 0 && columns > 0)
    _data = shared_array<Real>(new Real[rows*columns]);
}

/// Constructs a matrix from a vector
/**
 * \param v the vector
 * \param transpose if true the matrix will be transposed (1xn instead of nx1)
 */
MatrixN::MatrixN(const VectorN& v, bool transpose)
{
  _rows = _columns = _capacity = 0;
  set(v, transpose);
}

/// Constructs a matrix from an array
/**
 * \param rows the number of rows of the matrix
 * \param columns the number of columns of the matrix
 * \param array an array of rows*columns Real values in row-major format
 */
MatrixN::MatrixN(unsigned rows, unsigned columns, const Real* array)
{
  _rows = _columns = _capacity = 0;
  set(rows, columns, array);
}

/// Constructs a matrix from an array
/**
 * \param rows the number of rows of the matrix
 * \param columns the number of columns of the matrix
 * \param array an array of rows*columns Real values in row-major format;
 * \note MatrixN will modify the contents of array until the MatrixN object is
 *       destructed
 */
MatrixN::MatrixN(unsigned rows, unsigned columns, shared_array<Real> array)
{
  _capacity = rows*columns;
  _rows = rows;
  _columns = columns;
  _data = array;
}

/// Constructs a matrix from a Matrix3
MatrixN::MatrixN(const Matrix3& m)
{
  _rows = _columns = _capacity = 0;
  set(m);
}

/// Constructs a matrix from a Matrix3
MatrixN::MatrixN(const Matrix4& m)
{
  _rows = _columns = _capacity = 0;
  set(m);
}

/// Copy constructor
MatrixN::MatrixN(const MatrixN& source)
{
  _rows = _columns = _capacity = 0;
  MatrixN::operator=(source);
}

/// Copy (assignment) constructor
/**
 * \note destroys source on return
 */
MatrixN::MatrixN(MatrixN& source)
{
  _rows = source._rows;
  _columns = source._columns;
  _capacity = source._capacity;
  std::swap(_data, source._data);
  source._rows = source._columns = 0;
  source._capacity = 0;
}

/// Creates an iterator to the start of a block
BlockIterator MatrixN::block_start(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  unsigned i = row_start;
  unsigned j = col_start;
  unsigned rows = row_end - row_start;
  unsigned cols = col_end - col_start;

  if (i+rows > _rows || j+cols > _columns)
    throw MissizeException();

  BlockIterator b;
  b._count = 0;
  b._sz = rows*cols;
  b._data_start = (_data) ? _data.get() : NULL;
  b._current_data = (_data) ? &_data[j*_rows+i] : NULL;
  b._matrix_rows = _rows;
  b._matrix_columns = _columns;
  b._block_rows = rows;
  b._block_columns = cols;

  return b;
}

/// Creates an iterator to the end of a block
BlockIterator MatrixN::block_end(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  unsigned rows = row_end - row_start;
  unsigned cols = col_end - col_start;
  return block_start(row_start, row_end, col_start, col_end) + (rows*cols);
}

/// Computes the l-infinity norm of this matrix
Real MatrixN::norm_inf() const
{
  Real nrm = (Real) 0.0;
  const unsigned NELM = _rows * _columns;
  for (unsigned i=0; i< NELM; i++)
    nrm = std::max(nrm, std::fabs(_data[i]));

  return nrm;
}

/// Sets a matrix from a Matrix3
MatrixN& MatrixN::set(const Matrix3& m)
{
  const unsigned SZ = 3;
  resize(SZ,SZ);
  CBLAS::copy(SZ*SZ, m.data(), 1, _data.get(), 1);
  return *this;
}

/// Sets a matrix from a Matrix4
MatrixN& MatrixN::set(const Matrix4& m)
{
  // resize the matrix
  const unsigned SZ = 4;
  resize(SZ,SZ);

  // copy
  const unsigned MCOL1 = 0, MCOL2 = 3, MCOL3 = 6, MCOL4 = 9;
  const unsigned COL1 = 0, COL2 = 4, COL3 = 8, COL4 = 12;
  const unsigned WX = 3, WY = 7, WZ = 11, WW = 15;
  const Real* mdata = m.begin();
  for (unsigned i=0; i< 3; i++)
  {
    _data[COL1 + i] = mdata[MCOL1 + i];
    _data[COL2 + i] = mdata[MCOL2 + i];
    _data[COL3 + i] = mdata[MCOL3 + i];
    _data[COL4 + i] = mdata[MCOL4 + i];
  }

  // set last row
  _data[WX] = _data[WY] = _data[WZ] = (Real) 0.0;
  _data[WW] = (Real) 1.0;
  return *this;
}

/// Sets a matrix from a vector
MatrixN& MatrixN::set(const VectorN& v, bool transpose)
{
  // resize the matrix
  if (!transpose)
    resize(v.size(), 1);
  else
    resize(1, v.size());

  if (_rows > 0 && _columns > 0)
    CBLAS::copy(v.size(), v.data(), 1, _data.get(), 1);

  return *this;
}

/// Sets a matrix from an array
MatrixN& MatrixN::set(unsigned rows, unsigned columns, const Real* array)
{
  // resize the matrix
  resize(rows, columns);

  // copy the data
  if (_rows > 0 && _columns > 0)
    CBLAS::copy(_rows*_columns, array, 1, _data.get(), 1);

  return *this;
}

/// Sets a matrix from an array
MatrixN& MatrixN::set(unsigned rows, unsigned columns, shared_array<Real> array)
{
  // just use this array
  _rows = rows;
  _columns = columns;
  _capacity = rows*columns;
  _data = array;
  return *this;
}

/// Sets up a block matrix using a variable number of matrices
/**
 * \note undefined behavior will result if the passed matrices do not span
 *       all of the columns of the matrix!
 */
void MatrixN::set_sub_row_block(unsigned row, ...)
{
  // initialize the variable argument list
  std::va_list lst;
  va_start(lst, row);

  // determine how many rows from the first matrix
  unsigned col = 0;
  const MatrixN& m = *va_arg(lst, const MatrixN*);
  const unsigned NROWS = m.rows();
  set_sub_mat(row, col, m);
  col += m.columns();

  // continue reading blocks
  while (col < columns())
  {
    const MatrixN& m = *va_arg(lst, const MatrixN*);
    assert(NROWS == m.rows());
    set_sub_mat(row, col, m);
    col += m.columns();
  }

  // do the variable argument end
  va_end(lst);
}

/// Constructs a MatrixN using a variable number of double values
/**
 * \note the values are given row by row
 * \note There is no means in C++ to check the types of a list of variable
 * arguments.  If the variable arguments are not of type double, then
 * unexpected values will result. Constructing a matrix using the
 * statement MatrixN::construct_variable(2, 2, 1.0, 1.0, 1.0, 0) is incorrect 
 * because the programmer has assumed that the integer 0 will be converted to 
 * a double type; this is not the case.
 */
MatrixN MatrixN::construct_variable(unsigned rows, unsigned cols, ...)
{
  MatrixN m(rows, cols);
  std::va_list lst;
  va_start(lst, cols);
  for (unsigned i=0; i< rows; i++)
    for (unsigned j=0; j< cols; j++)
      m(i,j) = (Real) va_arg(lst, double);
  va_end(lst);
  return m;
}

/// Removes a row from the matrix
/**
 * \note downsizes the matrix -- does not reallocate memory
 */
MatrixN& MatrixN::remove_row(unsigned i)
{
  VectorN tmp(_columns);
  for (unsigned j=i+1; j< _rows; j++)
  {
    get_row(j, tmp);
    set_row(j-1, tmp);
  }

  // downsize the matrix
  _rows--;
  return *this;
}

/// Removes a column from the matrix
/**
 * \note downsizes the matrix -- does not reallocate memory
 */
MatrixN& MatrixN::remove_column(unsigned i)
{
  VectorN tmp(_rows);
  for (unsigned j=i+1; j< _columns; j++)
  {
    get_column(j, tmp);
    set_column(j-1, tmp);
  }

  // downsize the matrix
  _columns--;
  return *this;
}

/// Sets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param v a vector of the same length as the number of columns of <b>this</b>
 */
MatrixN& MatrixN::set_row(unsigned i, const VectorN& v)
{
  if (i >= _rows)
    throw InvalidIndexException();
  if (v.size() != _columns)
    throw MissizeException();

  // use BLAS for copying
  if (_columns > 0)
    CBLAS::copy(_columns, v.data(), 1, _data.get()+i, _rows);

  return *this;
}

/// Sets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param v a 3D vector
 * \note the number of columns of this must be three!
 */
MatrixN& MatrixN::set_row(unsigned i, const Vector3& v)
{
  if (i >= _rows)
    throw InvalidIndexException();
  if (_columns != v.size())
    throw MissizeException();
 
  // use BLAS for copying
  if (_columns > 0)
    CBLAS::copy(_columns, v.begin(), 1, _data.get()+i, _rows);

  return *this;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
VectorN MatrixN::get_row(unsigned i) const
{
  VectorN v;
  get_row(i, v);
  return v;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param target an array of size columns() or greater
 */
void MatrixN::get_row(unsigned i, Real* target) const
{
  if (i >= _rows)
    throw InvalidIndexException();
 
 if (_columns > 0)
    CBLAS::copy(_columns, _data.get()+i, _rows, target, 1);
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
VectorN& MatrixN::get_row(unsigned i, VectorN& result) const
{
  if (i >= _rows)
    throw InvalidIndexException();

  // resize the vector
  result.resize(_columns);

  // use BLAS for copying
  if (_columns > 0)
    CBLAS::copy(_columns, _data.get()+i, _rows, result.data(), 1);

  return result;
}

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param v a vector of the same length as the number of rows of <b>this</b> 
 */
MatrixN& MatrixN::set_column(unsigned i, const VectorN& v)
{
  if (i >= _columns)
    throw InvalidIndexException();
  if (v.size() != _rows)
    throw MissizeException();

  // use BLAS for copying
  if (_rows > 0)
    CBLAS::copy(_rows, v.data(), 1, _data.get()+i*_rows, 1);

  return *this;
}

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column 
 * \param v a 3D vector
 * \note the number of rows of this must be three!
 */
MatrixN& MatrixN::set_column(unsigned i, const Vector3& v)
{
  if (i >= _columns)
    throw InvalidIndexException();
  if (_rows != 3)
    throw MissizeException();

  // use BLAS for copying
  if (_columns > 0)
    CBLAS::copy(_columns, v.data(), 1, _data.get()+i, _rows);

  return *this;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 */
VectorN MatrixN::get_column(unsigned i) const
{
  VectorN v;
  get_column(i, v);
  return v;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param target an array of size rows() or greater
 */
void MatrixN::get_column(unsigned i, Real* target) const
{
  if (i >= _columns)
    throw InvalidIndexException();

  if (_rows > 0)
    CBLAS::copy(_rows, _data.get()+i*_rows, 1, target, 1);
}

/// Gets the specified column of the matrix
/*
 * \param i the 0-index of the column
 */
VectorN& MatrixN::get_column(unsigned i, VectorN& result) const
{
  if (i >= _columns)
    throw InvalidIndexException();

  // resize the vector
  result.resize(_rows);

  // use BLAS for copying
  if (_rows > 0)
    CBLAS::copy(_rows, _data.get()+i*_rows, 1, result.data(), 1);

  return result;
}

/// Determines whether this matrix is equal to another within the specified tolerance
bool MatrixN::epsilon_equals(const MatrixN& m, Real epsilon)
{
  return epsilon_equals(*this, m, epsilon);
}

/// Determines whether two matrices are equal within the specified tolerance
bool MatrixN::epsilon_equals(const MatrixN& m1, const MatrixN& m2, Real epsilon)
{
  // get the dimensions M and N
  unsigned M = m1.rows();
  unsigned N = m1.columns();

  // check that the two matrices are of the same dimension
  if (M != m2.rows() || N != m2.columns())
    throw MissizeException();

  // get the two arrays
  const Real* x1 = m1.data();
  const Real* x2 = m2.data();
      
  for (unsigned i=0; i< M*N; i++)
    if (std::fabs(x1[i] - x2[i]) > epsilon)  
      return false;
  return true;
}

/// Returns the zero matrix
/**
 * \param rows the number of rows
 * \param columns the number of columns of the matrix
 * \return a rows x columns zero matrix
 */
MatrixN MatrixN::zero(unsigned int rows, unsigned int columns)
{
  // create the new matrix
  MatrixN m(rows, columns);

  // get the array for the matrix, and set all values to zero
  Real* x = m.data();
  for (unsigned i=0; i< rows*columns; i++)
    x[i] = (Real) 0.0;

  return m; 
}

/// Resizes this matrix, optionally preserving its existing elements
/**
 * \note this method keeps from reallocating memory unless absolutely
 *       necessary (i.e., if the matrix grows or preserve=true, then memory
 *       will need to be reallocated.
 */
MatrixN& MatrixN::resize(unsigned rows, unsigned columns, bool preserve)
{
  shared_array<Real> newdata;

  // if the matrix is already the proper size, exit
  if (_rows == rows && _columns == columns)
    return *this;

  // if we can downsize, do that..
  if (rows*columns <= _capacity && (_rows == rows || !preserve))
  {
    _rows = rows;
    _columns = columns;
    return *this;
  }

  // create a new array
  if (rows > 0 && columns > 0)
    newdata = shared_array<Real>(new Real[rows*columns]);

  // preserve existing elements, if desired
  if (preserve && _rows > 0 && _columns > 0)
  {
    const unsigned n = std::min(_rows, rows);
    for (unsigned i=0; i< _columns; i++)
      CBLAS::copy(n, _data.get()+_rows*i,1,newdata.get()+rows*i,1);
  }

  // set the new data
  _data = newdata;
  _rows = rows;
  _columns = columns;
  _capacity = rows*columns;
  return *this;
}

/// Sets the matrix to the zero matrix
MatrixN& MatrixN::set_zero()
{
  for (unsigned i=0; i< _rows*_columns; i++)
    _data[i] = (Real) 0.0;
  return *this;
}

/// Sets this matrix to its transpose
MatrixN& MatrixN::transpose()
{
  // do fastest transpose first (if possible)
  if (_rows == 1 || _columns == 1)
  {
    std::swap(_rows, _columns);
    return *this;
  }  

  // do second fastest transpose, if possible
  if (_rows == _columns)
  {
    for (unsigned i=0; i< _rows; i++)
      for (unsigned j=i+1; j< _rows; j++)
        std::swap(_data[i*_columns+j], _data[j*_rows+i]);

    return *this;
  }

  // do slowest transpose operation
  SAFESTATIC FastThreadable<MatrixN> n;
  n().resize(_columns, _rows);
  Real* ndata = n().data();
  for (unsigned i=0; i< _rows; i++)
    for (unsigned j=0; j< _columns; j++)
      ndata[i*_columns + j] = _data[j*_rows + i];
  this->copy_from(n()); 

  return *this;
}

/// Determines the transpose of a matrix and stores the result in a given matrix
MatrixN& MatrixN::transpose(const MatrixN& m, MatrixN& result)
{
  // resize the result
  result.resize(m.columns(), m.rows());

  const Real* mdata = m.data();
  Real* rdata = result.data();
  for  (unsigned i=0; i< m.rows(); i++)
    for (unsigned j=0; j< m.columns(); j++)
      rdata[i*m.columns()+j] = mdata[j*m.rows()+i];
  
  return result;
}

/// Determines the transpose of the matrix
MatrixN MatrixN::transpose(const MatrixN& m)
{
  MatrixN n = m;
  return n.transpose();
}

/// Sets the specified sub vector
/**
 * \param row_start the row to start (inclusive)
 * \param col_start the column to start (inclusive)
 * \param v the source vector
 * \param transpose set to true if the vector is to be transposed
 * \note fails assertion if v is too large to insert into this
 */
MatrixN& MatrixN::set_sub_mat(unsigned row_start, unsigned col_start, const VectorN& v, bool transpose)
{
  if ((!transpose && (row_start + v.size() > _rows || col_start >= _columns)) ||
      (transpose && (row_start >= _rows || col_start + v.size() > _columns)))
    throw MissizeException();

  // copy only column of v using BLAS
  if (!transpose)
    CBLAS::copy(v.size(), v.data(), 1, _data.get()+row_start+col_start*_rows, 1);
  else
    CBLAS::copy(v.size(), v.data(), 1, _data.get()+row_start+col_start*_rows, _rows);

  return *this;
}

/// Gets the specified sub matrix as a vector
/**
 * \param row_start the row to start (inclusive)
 * \param row_end the row to end (exclusive)
 * \param col_start the column to start (inclusive)
 * \param col_end the column to end (exclusive)
 * \return a (row_end - row_start) x (col_end - col_start) sized matrix
 */
VectorN& MatrixN::get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end, VectorN& result) const
{
  // check whether we can exit early
  if (row_end == row_start)
  {
    result.set_zero(0);
    return result;
  }
  else if (col_end == col_start)
  {
    result.set_zero(0); 
    return result;
  }

  // error checking: indices too big
  if (row_end - row_start > 1 && col_end - col_start > 1)
    throw std::runtime_error("MatrixN::get_sub_mat() called with vector argument but matrix is required for these indices!");
 
  // error checking: invalid indices
  if (row_end > rows() || col_end > columns() || row_end < row_start || col_end < col_start)
    throw InvalidIndexException();

  // determine whether we are doing the vector or the transpose
  if (row_end - row_start == 1)
  {
    // doing the transpose
    result.resize(col_end - col_start);

    // do the copy
    CBLAS::copy(result.size(), begin()+col_start*rows()+row_start, rows(), result.begin(), 1);
  }
  else
  {
    // no transpose
    result.resize(row_end - row_start);

    // do the copy
    CBLAS::copy(result.size(), begin()+col_start*rows()+row_start, 1, result.begin(), 1); 
  }

  return result;
}

/// Gets the specified sub matrix
/**
 * \param row_start the row to start (inclusive)
 * \param row_end the row to end (exclusive)
 * \param col_start the column to start (inclusive)
 * \param col_end the column to end (exclusive)
 * \return a (row_end - row_start) x (col_end - col_start) sized matrix
 */
MatrixN MatrixN::get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end) const
{
  MatrixN m;
  return get_sub_mat(row_start, row_end, col_start, col_end, m);
}

/// Sets this matrix to that of another
MatrixN& MatrixN::operator=(const MatrixN& m)
{
  resize(m.rows(), m.columns());

  // copy m over this using BLAS
  if (_rows > 0 && _columns > 0)
    CBLAS::copy(m.rows()*m.columns(), m.data(), 1, _data.get(), 1);

  return *this;
}

/// Assigns this method to another
/**
 * \note destroys source on return
 */
MatrixN& MatrixN::operator=(MatrixN& source)
{
  _rows = source._rows;
  _columns = source._columns;
  source._rows = source._columns = 0;
  std::swap(_capacity, source._capacity);
  std::swap(_data, source._data);
  return *this;
}

/// Multiplies this matrix by a scalar
MatrixN MatrixN::operator*(Real scalar) const
{
  MatrixN this_copy(*this);
  this_copy *= scalar;
  return this_copy;
}

/// Multiplies this matrix by a vector
VectorN MatrixN::mult(const VectorN& v) const
{
  VectorN result;
  mult(v, result);
  return result;
}

/// Multiplies the transpose of this matrix by a vector
VectorN MatrixN::transpose_mult(const VectorN& v) const
{
  VectorN result;
  transpose_mult(v, result);
  return result;
}

/// Multiplies this matrix by another matrix
MatrixN MatrixN::mult(const MatrixN& m) const
{
  MatrixN result;
  mult(m, result);
  return result;
}

/// Multiplies the transpose of this matrix by another matrix
MatrixN MatrixN::transpose_mult(const MatrixN& m) const
{
  MatrixN result;
  transpose_mult(m, result);
  return result;
}

/// Multiplies this matrix by the transpose of another matrix
MatrixN MatrixN::mult_transpose(const MatrixN& m) const
{
  MatrixN result;
  mult_transpose(m, result);
  return result;
}

/// Multiplies the transpose of this matrix by the transpose of another matrix
MatrixN MatrixN::transpose_mult_transpose(const MatrixN& m) const
{
  MatrixN result;
  transpose_mult_transpose(m, result);
  return result;
}

/// Multiplies this matrix by another in place
MatrixN& MatrixN::operator*=(Real scalar)
{
  // call BLAS scaling function
  if (_rows > 0 && _columns > 0)
    CBLAS::scal(_rows*_columns, scalar, _data.get(), 1);
  return *this;
}

/// Divides this matrix by a scalar in place
MatrixN& MatrixN::operator/=(Real scalar)
{
  // call BLAS scaling function
  if (_rows > 0 && _columns > 0)
    CBLAS::scal(_rows*_columns, (Real) 1.0/scalar, _data.get(), 1);
  return *this;
}

/// Multiplies this matrix by another in place
MatrixN& MatrixN::operator*=(const MatrixN& m)
{
  operator=(mult(m));
  return *this;
}

/// Adds two matrices together
MatrixN MatrixN::operator+(const MatrixN& m) const
{
  MatrixN n(*this);
  n += m;
  return n;
}

/// Subtracts m2 from m1
MatrixN MatrixN::operator-(const MatrixN& m) const
{
  MatrixN n(*this);
  n -= m;
  return n;
}

/// Returns the negation of this matrix
MatrixN MatrixN::operator-() const
{
  MatrixN m(_rows, _columns);
  if (_rows > 0 && _columns > 0) 
    std::transform(begin(), end(), m.begin(), std::negate<Real>());
  return m;
}

/// Negates this matrix in place
MatrixN& MatrixN::negate()
{
  unsigned n = _rows * _columns;
  for (unsigned i=0; i< n; i++)
    _data[i] = -_data[i];
  return *this;
}

/// Adds m to *this in place
MatrixN& MatrixN::operator+=(const MatrixN& m)
{
  if (_rows != m.rows() || _columns != m.columns())
    throw MissizeException(); 

  if (_rows > 0 && _columns > 0) 
    std::transform(begin(), end(), m.begin(), begin(), std::plus<Real>());
  return *this;
}

/// Subtracts m from *this in place
MatrixN& MatrixN::operator-=(const MatrixN& m)
{
  if (_rows != m.rows() || _columns != m.columns())
    throw MissizeException(); 
  
  if (_rows > 0 && _columns > 0) 
    std::transform(begin(), end(), m.begin(), begin(), std::minus<Real>());
  return *this;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const MatrixN& m)
{
  const unsigned OUTPUT_PRECISION = 8;

  if (m.rows() == 0 || m.columns() == 0)
  {
    out << "(empty)" << std::endl;
    return out;
  }
  
  for (unsigned i=0; i< m.rows(); i++)
  {
    for (unsigned j=0; j< m.columns()-1; j++)
      out << std::setprecision(OUTPUT_PRECISION) << m(i,j) << " ";
    out << std::setprecision(OUTPUT_PRECISION) << m(i,m.columns()-1) << std::endl;
  }
   
  return out;
}

/// Reads a matrix from the stream
std::istream& Moby::operator>>(std::istream& in, MatrixN& m)
{
  unsigned rows, columns;

  // read in the size
  in >> rows;
  in >> columns;

  // create the matrix
  m.resize(rows, columns);
  for (unsigned i=0; i< rows; i++)
    for (unsigned j=0; j< columns; j++) 
      in >> m(i,j);

  return in;
} 

/// Checks whether the given matrix is symmetric to the specified tolerance
bool MatrixN::is_symmetric(Real tolerance) const
{
  if (_rows != _columns)
    return false;

  // make sure that tolerance is positive 
  if (tolerance <= (Real) 0.0)
    tolerance = std::numeric_limits<Real>::epsilon() * norm_inf() * _rows;

  // check symmetry
  const unsigned LD = _rows;

  // loop over columns
  for (unsigned i=0, ii=0; i< _rows; i++, ii+= LD)
    for (unsigned j=0, jj=0; j< i; j++, jj+= LD)
      if (std::fabs(_data[jj+i] - _data[ii+j]) > tolerance)
        return false;

  return true;
}

/// Zeros the upper triangle of the matrix
MatrixN& MatrixN::zero_upper_triangle()
{
  // check for easy exit
  if (_rows == 0 || _columns == 0)
    return *this;

  // zero the upper triangle
  for (unsigned i=0, s=_rows; i< _rows; i++, s+= _rows+1)
    for (unsigned j=i+1, r=0; j< _columns; j++, r+= _rows)
      _data[r+s] = (Real) 0.0;

  return *this;
}

/// Zeros the lower triangle of the matrix
MatrixN& MatrixN::zero_lower_triangle()
{
  // check for easy exit
  if (_rows == 0 || _columns == 0)
    return *this;

  // zero the lower triangle
  for (unsigned i=1, s=1; i< _rows; i++, s++)
    for (unsigned j=0, r=0; j< std::min(i, _columns); j++, r+= _rows)
      _data[r+s] = (Real) 0.0;

  return *this;
}


/// Returns an identity matrix
MatrixN MatrixN::identity(unsigned n)
{
  MatrixN m;
  m.set_identity(n);
  return m;
}

/// Sets this matrix to the identity matrix
MatrixN& MatrixN::set_identity(unsigned i)
{
  resize(i,i);
  set_zero();
  for (unsigned i=0, j=0; i< _rows; i++, j+= _rows+1)
    _data[j] = (Real) 1.0;

  return *this;
}

/// Sets this matrix to the identity matrix
MatrixN& MatrixN::set_identity()
{
  if (_rows != _columns)
    throw MissizeException();
  set_identity(_rows);
}

/// Multiplies the diagonal matrix formed from d by the matrix m (slow, convenient method)
MatrixN MatrixN::diag_mult(const VectorN& d, const MatrixN& m) 
{
  MatrixN result;
  diag_mult(d, m, result);
  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix transpose(m) (slow, convenient method)
MatrixN MatrixN::diag_mult_transpose(const VectorN& d, const MatrixN& m)
{
  MatrixN result;
  diag_mult_transpose(d, m, result);
  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix m
MatrixN& MatrixN::diag_mult(const VectorN& d, const MatrixN& m, MatrixN& result)
{
  if (d.size() != m.rows())
    throw MissizeException();

  result.resize(d.size(), m.columns());
  for (unsigned i=0; i< m.columns(); i++)
    std::transform(d.begin(), d.end(), m.begin()+m.rows()*i, result.begin()+result.rows()*i, std::multiplies<Real>());

  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix transpose(m)
MatrixN& MatrixN::diag_mult_transpose(const VectorN& d, const MatrixN& m, MatrixN& result)
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
VectorN& MatrixN::diag_mult(const VectorN& d, const VectorN& v, VectorN& result)
{
  if (d.size() != v.size())
    throw MissizeException();

  result.resize(d.size());
  std::transform(d.begin(), d.end(), v.begin(), result.begin(), std::multiplies<Real>());
  return result;
}

/// Multiplies the diagonal matrix formed from d by the vector v (slow, convenient method)
VectorN MatrixN::diag_mult(const VectorN& d, const VectorN& v)
{
  VectorN result;
  return diag_mult(d, v, result);
}

/// Selects a submatrix from this (slow, convenient method)
MatrixN MatrixN::select_square(const vector<bool>& indices) const
{
  MatrixN result;
  return select_square(indices, result);
}

/// Selects a submatrix from this (square matrix)
MatrixN& MatrixN::select_square(const vector<bool>& indices, MatrixN& result) const
{
  const unsigned n = rows();
  assert(n == indices.size());

  // determine how many indices are selected and resize result vector
  unsigned nselect = 0;
  for (unsigned i=0; i< n; i++)
    if (indices[i])
      nselect++;
  result.resize(nselect, nselect);

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

