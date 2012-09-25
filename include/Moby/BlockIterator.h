/****************************************************************************
 * Copyright 2012 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _BLOCK_ITERATOR_H
#define _BLOCK_ITERATOR_H

#include <stdexcept>

namespace Moby {

/// A construct for iterating over a rectangular block of a matrix
class BlockIterator : public std::iterator<std::random_access_iterator_tag, Real>
{
  friend class MatrixN;

  public:
    BlockIterator()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _block_rows = 0;
      _block_columns = 0;
      _matrix_rows = 0;
      _matrix_columns = 0;
    }

    BlockIterator& operator+=(int n) 
    {
      assert(n >= 0);
      const unsigned NSKIP = _matrix_rows - _block_rows;
      for (int i=0; i< n; i++)
      {
        _count++;
        _current_data++;
        if (_count % _block_rows == 0)
          _current_data += NSKIP;
      }

      return *this; 
    }

    BlockIterator& operator-=(int n) 
    { 
      assert(n >= 0);
      const unsigned NSKIP = _matrix_rows - _block_rows;
      for (int i=0; i< n; i++)  
      {
        _count--;
        _current_data--;
        if (_count % _block_rows == 0)
          _current_data -= NSKIP;
      }

      return *this; 
    }

    Real& operator[](unsigned i) const
    {
      if (i > _sz)
        throw std::runtime_error("Data outside of scope!");
      const unsigned NSKIP = _matrix_rows - _block_rows;
      Real* data = _data_start;
      for (unsigned j=0; j< i; )
      {
        data++;
        if (++j % _block_rows == 0)
          data += NSKIP;
      }

      return *data;
    }

    BlockIterator operator+(int n) const
    {
      BlockIterator b = *this;
      b += n;
      return b;
    }

    BlockIterator operator-(int n) const
    {
      BlockIterator b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const BlockIterator& j) const
    {
      return _count < j._count;
    }

    bool operator>(const BlockIterator& j) const
    {
      return _count > j._count;
    }

    Real& operator*() const 
    {
      if (_count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const BlockIterator& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const BlockIterator& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _matrix_rows == b._matrix_rows &&
             _matrix_columns == b._matrix_columns &&
             _block_rows == b._block_rows &&
             _block_columns == b._block_columns);

      return (_count == b._count);
    }

    bool operator!=(const BlockIterator& j) const { return !operator==(j); }

    // prefix--
    BlockIterator& operator--() 
    { 
      _count--; 
      _current_data--;
      if (_count % _block_rows == 0)
        _current_data -= (_matrix_rows - _block_rows);

       return *this; 
    }

    // prefix++
    BlockIterator& operator++() 
    { 
      _count++;
      _current_data++;
      if (_count % _block_rows == 0)
        _current_data += (_matrix_rows - _block_rows);

      return *this;
    }

    // postfix--
    BlockIterator operator--(int n) 
    { 
      BlockIterator b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    BlockIterator operator++(int n) 
    {
      BlockIterator b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    BlockIterator& operator=(const BlockIterator& i)
    {
      _current_data = i._current_data;
      _count = i._count;
      _sz = i._sz;
      _data_start = i._data_start;
      _matrix_rows = i._matrix_rows;
      _matrix_columns = i._matrix_columns;
      _block_rows = i._block_rows;
      _block_columns = i._block_columns;
      return *this;
    }

  protected:
    int _count, _sz;
    Real* _data_start;
    Real* _current_data;
    unsigned _matrix_rows;
    unsigned _matrix_columns;
    unsigned _block_rows;
    unsigned _block_columns;
}; // end class

} // end namespace

#endif

