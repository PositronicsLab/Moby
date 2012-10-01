/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Gets the specified sub matrix
/**
 * \param row_start the row to start (inclusive)
 * \param row_end the row to end (exclusive)
 * \param col_start the column to start (inclusive)
 * \param col_end the column to end (exclusive)
 * \param transpose if true, stores the transpose of the submatrix into m
 * \return a (row_end - row_start) x (col_end - col_start) sized matrix
 */
template <class M>
M& MatrixN::get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end, M& m, bool transpose) const
{
  if (row_start > row_end || row_end > _rows || col_start > col_end || col_end > _columns)
    throw InvalidIndexException();

  // resize the matrix
  if (!transpose)
    m.resize(row_end - row_start, col_end - col_start);
  else
    m.resize(col_end - col_start, row_end - row_start);

  // see whether to exit now
  if (_rows == 0 || _columns == 0)
    return m;

  // copy each column using BLAS
  if (!transpose)
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), _data.get()+row_start+(col_start+i) * _rows, 1, m.data()+i*m.rows(), 1);
  else
    for (unsigned i=0; i< m.rows(); i++)
      CBLAS::copy(m.columns(), _data.get()+row_start+(col_start+i) * _rows, 1, m.data()+i, m.rows());
  
  return m;
}

/// Sets the specified sub matrix
/**
 * \param row_start the row to start (inclusive)
 * \param col_start the column to start (inclusive)
 * \param m the source matrix
 * \param transpose set to true if the matrix is to be transposed
 * \note fails assertion if m is too large to insert into this
 */
template <class M>
MatrixN& MatrixN::set_sub_mat(unsigned row_start, unsigned col_start, const M& m, bool transpose)
{
  if (!transpose)
  {
    if (row_start + m.rows() > _rows || col_start + m.columns() > _columns)
      throw MissizeException();
  }
  else if (row_start + m.columns() > _rows || col_start + m.rows() > _columns)
    throw MissizeException();

  // make sure there is data to copy
  if (m.rows() == 0 || m.columns() == 0)
    return *this;

  // copy each column of m using BLAS
  if (!transpose)
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), m.data()+i*m.rows(), 1, _data.get()+row_start+(col_start+i) * _rows, 1);
  else
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), m.data()+i*m.rows(), 1, _data.get()+row_start+col_start*_rows + i, _rows);

  return *this;
}

/// Gets a submatrix of columns (not necessarily a block)
template <class ForwardIterator>
MatrixN& MatrixN::select_columns(ForwardIterator col_start, ForwardIterator col_end, MatrixN& m) const
{
  // setup vectors of selections
  const unsigned ncols = std::distance(col_start, col_end);

  // resize matrix 
  m.resize(_rows, ncols);

  // make sure there is data to copy
  if (_rows == 0 || ncols == 0)
    return m;

  // populate m
  unsigned mi;
  ForwardIterator i;
  for (i=col_start, mi=0; i != col_end; i++, mi++)
    CBLAS::copy(_rows, begin()+_rows*(*i), 1, m.begin()+_rows*mi, 1);

  return m;
}

/// Gets a submatrix of rows (not necessarily a block)
template <class ForwardIterator>
MatrixN& MatrixN::select_rows(ForwardIterator row_start, ForwardIterator row_end, MatrixN& m) const
{
  // setup vectors of selections
  const unsigned nrows = std::distance(row_start, row_end);

  // resize matrix 
  m.resize(nrows, _columns);

  // make sure there is data to copy
  if (nrows == 0 || _columns == 0)
    return m;

  // populate m
  unsigned mi;
  ForwardIterator i;
  for (i=row_start, mi=0; i != row_end; i++, mi++)
    CBLAS::copy(_columns, begin()+*i, _rows, m.begin()+mi, nrows);

  return m;
}

/// Gets a submatrix of rows (not necessarily a block)
template <class ForwardIterator>
MatrixN MatrixN::select_columns(ForwardIterator col_start, ForwardIterator col_end) const
{
  MatrixN m;
  select_columns(col_start, col_end, m);
  return m;
}

/// Gets a submatrix of rows (not necessarily a block)
template <class ForwardIterator>
MatrixN MatrixN::select_rows(ForwardIterator row_start, ForwardIterator row_end) const
{
  MatrixN m;
  select_rows(row_start, row_end, m);
  return m;
}

/// Gets a submatrix (not necessarily a block)
template <class ForwardIterator1, class ForwardIterator2>
MatrixN& MatrixN::select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end, MatrixN& m) const
{
  // setup vectors of selections
  const unsigned nrows = std::distance(row_start, row_end);
  const unsigned ncols = std::distance(col_start, col_end);

  // resize matrix 
  m.resize(nrows, ncols);

  // make sure there is data to copy
  if (nrows == 0 || ncols == 0)
    return m;

  // populate m
  unsigned mi;
  ForwardIterator1 i;
  for (i=row_start, mi=0; i != row_end; i++, mi++)
  {
    unsigned mj;
    ForwardIterator2 j;

    // copy all flagged elements in column of A
    for (j=col_start, mj=0; j != col_end; j++, mj++)
      m(mi, mj) = operator()(*i,*j);
  }

  return m;
}

/// Gets a submatrix (not necessarily a block)
template <class ForwardIterator1, class ForwardIterator2>
MatrixN MatrixN::select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end) const
{
  MatrixN m;
  return select(row_start, row_end, col_start, col_end, m);
}

template <class ForwardIterator>
MatrixN& MatrixN::select_square(ForwardIterator start, ForwardIterator end, MatrixN& m) const
{
  const unsigned n = std::distance(start, end);

  // resize matrix
  m.resize(n,n);

  // populate m
  unsigned mi;
  ForwardIterator i;
  for (i=start, mi=0; i != end; i++, mi++)
  {
    unsigned mj;
    ForwardIterator j;

    // copy all flagged elements in column of A
    for (j=start, mj=0; j != end; j++, mj++)
      m(mi, mj) = operator()(*i,*j);
  }

  return m;
}

template <class ForwardIterator>
MatrixN MatrixN::select_square(ForwardIterator start, ForwardIterator end) const
{
  // create new matrix
  MatrixN m;
  return select_square(start, end, m);
}

