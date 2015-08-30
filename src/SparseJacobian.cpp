/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/MissizeException.h>
#include <Moby/SparseJacobian.h>

using std::vector;
using namespace Ravelin;
using namespace Moby;

/// Multiplies this sparse Jacobian by a vector
VectorNd& SparseJacobian::mult(const VectorNd& x, VectorNd& result) const
{
  VectorNd tmp;

  // check for proper size
  if (cols != x.size())
    throw MissizeException();

  // set the result size
  result.set_zero(rows);

  // look for number of blocks 
  if (blocks.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // get the relevant part of the segments
    SharedConstVectorNd x_segment =  x.segment(blocks[i].st_col_idx, blocks[i].st_col_idx+blocks[i].columns());
    SharedVectorNd result_segment = result.segment(blocks[i].st_row_idx, blocks[i].st_row_idx+blocks[i].rows());

    // do the computation
    blocks[i].block.mult(x_segment, tmp);
    result_segment += tmp;
  }

  return result;
}

/// Multiplies this sparse Jacobian by a matrix 
MatrixNd& SparseJacobian::mult(const MatrixNd& x, MatrixNd& result) const
{
  MatrixNd tmp;

  // check for proper size
  if (cols != x.rows())
    throw MissizeException();

  // set the result size
  result.set_zero(rows, x.columns());

  // look for number of blocks 
  if (blocks.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // assume block i is of size r x c
    // then block of x must be of size c x x.columns() 
    // and result must be of size r x x.columns()
    const unsigned R = blocks[i].rows();
    const unsigned C = blocks[i].columns();

    // get the relevant blocks
    SharedConstMatrixNd x_block =  x.block(blocks[i].st_col_idx, blocks[i].st_col_idx+C, 0, x.columns());
    SharedMatrixNd result_block = result.block(blocks[i].st_row_idx, blocks[i].st_row_idx+R, 0, x.columns());

    // do the computation
    blocks[i].block.mult(x_block, tmp);
    result_block += tmp;
  }

  return result;
}

/// Multiplies the transpose of this sparse Jacobian by a matrix 
MatrixNd& SparseJacobian::transpose_mult(const MatrixNd& x, MatrixNd& result) const
{
  MatrixNd tmp;

  // check for proper size
  if (rows != x.rows())
    throw MissizeException();

  // set the result size
  result.set_zero(cols, x.columns());

  // look for number of blocks 
  if (blocks.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // assume block i is of size r x c
    // then block of x must be of size c x x.columns() 
    // and result must be of size r x x.columns()
    const unsigned R = blocks[i].columns();
    const unsigned C = blocks[i].rows();

    // get the relevant blocks
    SharedConstMatrixNd x_block =  x.block(blocks[i].st_row_idx, blocks[i].st_row_idx+C, 0, x.columns());
    SharedMatrixNd result_block = result.block(blocks[i].st_col_idx, blocks[i].st_col_idx+R, 0, x.columns());

    // do the computation
    blocks[i].block.transpose_mult(x_block, tmp);
    result_block += tmp;
  }

  return result;
}

/// Multiples this sparse Jacobian by the transpose of another sparse Jacobian
MatrixNd& SparseJacobian::mult_transpose(const SparseJacobian& M, MatrixNd& result) const
{
  MatrixNd tmp;

  // set the result size
  result.set_zero(rows, M.rows);

  // get the blocks from M
  const std::vector<MatrixBlock>& x = M.blocks;

  // look for number of blocks 
  if (blocks.size() == 0 || x.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // assume block i is of size r x c
    // then input block of x must be of size c x d 
    // and result must be of size r x d
    const unsigned RSTART = blocks[i].st_col_idx;
    const unsigned CSTART = blocks[i].st_row_idx;
    const unsigned REND = RSTART + blocks[i].rows();
    const unsigned CEND = CSTART + blocks[i].columns();

    // loop over each input block
    for (unsigned j=0; j< x.size(); j++)
    {
      // get x column start and end
      const unsigned X_CSTART = x[j].st_col_idx;
      const unsigned X_CEND = X_CSTART + x[j].columns();

      // see whether the two blocks overlap 
      if (X_CEND < CSTART || X_CSTART >= CEND)
        continue;

      // get x row start and end
      const unsigned X_RSTART = x[j].st_row_idx;
      const unsigned X_REND = X_RSTART + x[j].rows();

      // get the common columns start and end
      const unsigned C_CSTART = std::max(CSTART, X_CSTART);
      const unsigned C_CEND = std::min(CEND, X_CEND);

      // get the result block - it's easiest
      SharedMatrixNd result_block = result.block(RSTART, REND, X_RSTART, X_REND);

      // get the appropriate blocks of the input matrices
      SharedConstMatrixNd i_block =  blocks[i].block.block(0, blocks[i].block.rows(), C_CSTART - CSTART, C_CEND - CSTART);
      SharedConstMatrixNd x_block =  x[j].block.block(0, x[j].rows(), C_CSTART - X_CSTART, C_CEND - X_CSTART);

      // do the computation
      i_block.mult_transpose(x_block, tmp);
      result_block += tmp;
    }
  }


  // resize the result matrix
  return result;
}

/// Multiplies this sparse Jacobian by a block diagonal matrix 
MatrixNd& SparseJacobian::mult(const vector<MatrixBlock>& x, unsigned result_cols, MatrixNd& result) const
{
  MatrixNd tmp;

  // set the result size
  result.set_zero(rows, result_cols);

  // look for number of blocks 
  if (blocks.size() == 0 || x.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // assume block i is of size r x c
    // then input block of x must be of size c x d 
    // and result must be of size r x d
    const unsigned RSTART = blocks[i].st_col_idx;
    const unsigned CSTART = blocks[i].st_row_idx;
    const unsigned REND = RSTART + blocks[i].rows();
    const unsigned CEND = CSTART + blocks[i].columns();

    // loop over each input block
    for (unsigned j=0; j< x.size(); j++)
    {
      // get x row start and end
      const unsigned X_RSTART = x[j].st_row_idx;
      const unsigned X_REND = X_RSTART + x[j].rows();

      // see whether the two blocks overlap 
      if (X_REND < CSTART || X_RSTART >= CEND)
        continue;

      // get x column start and end
      const unsigned X_CSTART = x[j].st_col_idx;
      const unsigned X_CEND = X_CSTART + x[j].columns();

      // get the common columns start and end
      const unsigned C_CSTART = std::max(CSTART, X_RSTART);
      const unsigned C_CEND = std::min(CEND, X_REND);

      // get the result block - it's easiest
      SharedMatrixNd result_block = result.block(RSTART, REND, X_CSTART, X_CEND);

      // get the appropriate blocks of the input matrices
      SharedConstMatrixNd i_block =  blocks[i].block.block(0, blocks[i].block.rows(), C_CSTART - CSTART, C_CEND - CSTART);
      SharedConstMatrixNd x_block =  x[j].block.block(C_CSTART - X_RSTART, C_CEND - X_RSTART, 0, x[j].columns());

      // do the computation
      i_block.mult(x_block, tmp);
      result_block += tmp;
    }
  }

  return result;
}

/// Multiplies this sparse Jacobian by a block diagonal matrix 
MatrixNd& SparseJacobian::mult(const vector<MatrixNd>& x, MatrixNd& result) const
{
  MatrixNd tmp;
  vector<unsigned> st_row_idx(x.size());

  // get the total number of columns
  st_row_idx.push_back(0);
  for (unsigned i=0; i< x.size(); i++)
  {
    assert(x[i].rows() == x[i].columns());
    st_row_idx.push_back(st_row_idx.back() + x[i].rows());
  }

  // set the result size
  result.set_zero(rows, st_row_idx.back());

  // look for number of blocks 
  if (blocks.size() == 0 || x.size() == 0)
    return result;

  // loop over each block
  for (unsigned i=0; i< blocks.size(); i++)
  {
    // assume block i is of size r x c
    // then input block of x must be of size c x d 
    // and result must be of size r x d
    const unsigned R = blocks[i].rows();
    const unsigned C = blocks[i].columns();

    // loop over each input block
    for (unsigned j=0, k=0; j< x.size(); j++)
    {
      // see whether the two correspond
      if (blocks[i].st_col_idx != k)
      {
        k += x[j].rows();
        continue;
      }

      // get D
      const unsigned D = x[j].columns(); 

      // verify that the blocks are the appropriate size
      if (C != D)
        throw MissizeException();

      // get the relevant block of the result
      SharedMatrixNd result_block = result.block(blocks[i].st_row_idx, blocks[i].st_row_idx+R, st_row_idx[j], st_row_idx[j]+D);

      // do the computation
      blocks[i].block.mult(x[j], tmp);
      result_block += tmp;
      
      // update k
      k += x[j].rows();
    }
  }

  return result;
}

/// Converts a sparse Jacobian to a dense matrix (for debugging purposes)
MatrixNd& SparseJacobian::to_dense(MatrixNd& M) const
{
  M.set_zero(rows, cols);

  for (unsigned i=0; i< blocks.size(); i++)
  {
    const unsigned R = blocks[i].rows();
    const unsigned C = blocks[i].columns();
    M.block(blocks[i].st_row_idx, blocks[i].st_row_idx+R, blocks[i].st_col_idx, blocks[i].st_col_idx+C) = blocks[i].block;
  }

  return M;
}

