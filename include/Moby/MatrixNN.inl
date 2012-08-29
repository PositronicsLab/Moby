/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

template <class ForwardIterator>
MatrixNN& MatrixNN::select(ForwardIterator start, ForwardIterator end, MatrixNN& m) const
{
  const unsigned n = std::distance(start, end);

  // resize matrix
  m.resize(n);

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
MatrixNN MatrixNN::select(ForwardIterator start, ForwardIterator end) const
{
  // create new matrix
  MatrixNN m;
  return select(start, end, m);
}

