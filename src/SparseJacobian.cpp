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
VectorNd& SparseJacobian::mult(const VectorNd& x, VectorNd& result)
{
  VectorNd tmp;

  // check for proper size
  if (_cols != x.size())
    throw MissizeException();

  // set the result size
  result.set_zero(_rows);

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
MatrixNd& SparseJacobian::mult(const MatrixNd& x, MatrixNd& result)
{
  MatrixNd tmp;

  // check for proper size
  if (_cols != x.rows())
    throw MissizeException();

  // set the result size
  result.set_zero(_rows, x.columns());

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

/// Multiplies this sparse Jacobian by a block diagonal matrix 
MatrixNd& SparseJacobian::mult(const vector<MatrixBlock>& x, unsigned result_cols, MatrixNd& result)
{
  MatrixNd tmp;

  // set the result size
  result.set_zero(_rows, result_cols);

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
    for (unsigned j=0; i< x.size(); j++)
    {
      // see whether the two correspond
      if (blocks[i].st_col_idx != x[j].st_row_idx)
        continue;

      // get D
      const unsigned D = x[j].columns(); 

      // verify that the blocks are the appropriate size
      if (C != x[j].rows())
        throw MissizeException();

      // get the relevant block of the result
      SharedMatrixNd result_block = result.block(blocks[i].st_row_idx, blocks[i].st_row_idx+R, x[j].st_col_idx, x[j].st_col_idx+D);

      // do the computation
      blocks[i].block.mult(x[j].block, tmp);
      result_block += tmp;
    }
  }

  return result;
}




