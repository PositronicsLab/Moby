#include <Moby/SparseJacobian.h>
#include "gtest/gtest.h"

using std::vector;
using namespace Ravelin;
using namespace Moby;

TEST(MultInertia, Regular)
{
  
}

TEST(MultInertia, Block)
{
  MatrixNd J;
  SparseJacobian Js;

  // setup the sparse Jacobian
  J.set_zero(10, 20);
  Js.rows = 10;
  Js.cols = 20;
  for (unsigned i=0; i< 10; i++)
  {
    Js.blocks.push_back(MatrixBlock());
    Js.blocks.back().block.resize(1,10);
    Js.blocks.back().st_row_idx = i;
    Js.blocks.back().st_col_idx = i;
    for (unsigned j=0; j< 10; j++)
    {
      Js.blocks.back().block(0,j) = (double) rand()/RAND_MAX;
      J(i,j+i) = Js.blocks.back().block(0,j);
    }
  }

  // now setup 10 2x2 matrices
  MatrixNd dense = MatrixNd::zero(20,20);
  vector<MatrixBlock> x;
  for (unsigned i=0; i< 10; i++)
  {
    x.push_back(MatrixBlock());
    x[i].st_row_idx = i*2;
    x[i].st_col_idx = i*2;
    x[i].block.resize(2,2);
    for (unsigned j=0; j< 2; j++)
      for (unsigned k=0; k< 2; k++)
        x[i].block(j,k) = (double) rand() / RAND_MAX;
    dense.block(i*2, i*2+2, i*2, i*2+2) = x[i].block;
  }

  // do the dense version
  MatrixNd result_dense;
  J.mult(dense, result_dense);

  // do the sparse multiplication 
  MatrixNd result_sparse;
  Js.mult(x, 20, result_sparse);

  // check the results
  for (unsigned i=0; i< result_dense.rows(); i++)
    for (unsigned j=0; j< result_dense.columns(); j++)
      EXPECT_NEAR(result_dense(i,j), result_sparse(i,j), 1e-6);
}

TEST(Mult, Vector)
{
  MatrixNd J;
  SparseJacobian Js;

  // setup the sparse Jacobian
  J.set_zero(10, 20);
  Js.rows = 10;
  Js.cols = 20;
  for (unsigned i=0; i< 10; i++)
  {
    Js.blocks.push_back(MatrixBlock());
    Js.blocks.back().block.resize(1,10);
    Js.blocks.back().st_row_idx = i;
    Js.blocks.back().st_col_idx = i;
    for (unsigned j=0; j< 10; j++)
    {
      Js.blocks.back().block(0,j) = (double) rand()/RAND_MAX;
      J(i,j+i) = Js.blocks.back().block(0,j);
    }
  }

  // randomly setup the vector
  VectorNd v(J.columns());
  for (unsigned i=0; i< v.size(); i++)
    v[i] = (double) rand() / RAND_MAX;

  // do the sparse Jacobian
  VectorNd result_sparse;
  Js.mult(v, result_sparse);

  // do the dense version
  VectorNd result_dense;
  J.mult(v, result_dense);

  // check the results
  for (unsigned i=0; i< result_dense.size(); i++)
    EXPECT_NEAR(result_dense[i], result_sparse[i], 1e-6);
}

TEST(Mult, Matrix)
{
  MatrixNd J;
  SparseJacobian Js;

  // setup the sparse Jacobian
  J.set_zero(10, 20);
  Js.rows = 10;
  Js.cols = 20;
  for (unsigned i=0; i< 10; i++)
  {
    Js.blocks.push_back(MatrixBlock());
    Js.blocks.back().block.resize(1,10);
    Js.blocks.back().st_row_idx = i;
    Js.blocks.back().st_col_idx = i;
    for (unsigned j=0; j< 10; j++)
    {
      Js.blocks.back().block(0,j) = (double) rand()/RAND_MAX;
      J(i,j+i) = Js.blocks.back().block(0,j);
    }
  }

  // randomly setup the vector
  MatrixNd V(J.columns(), 5);
  for (unsigned j=0; j< 5; j++)
    for (unsigned i=0; i< V.rows(); i++)
      V(i,j) = (double) rand() / RAND_MAX;

  // do the sparse Jacobian
  MatrixNd result_sparse;
  Js.mult(V, result_sparse);

  // do the dense version
  MatrixNd result_dense;
  J.mult(V, result_dense);

  // check the results
  for (unsigned i=0; i< result_dense.rows(); i++)
    for (unsigned j=0; j< result_dense.columns(); j++)
      EXPECT_NEAR(result_dense(i,j), result_sparse(i,j), 1e-6);
}

TEST(TransposeMult, Matrix)
{
  MatrixNd J;
  SparseJacobian Js;

  // setup the sparse Jacobian
  J.set_zero(10, 20);
  Js.rows = 10;
  Js.cols = 20;
  for (unsigned i=0; i< 10; i++)
  {
    Js.blocks.push_back(MatrixBlock());
    Js.blocks.back().block.resize(1,10);
    Js.blocks.back().st_row_idx = i;
    Js.blocks.back().st_col_idx = i;
    for (unsigned j=0; j< 10; j++)
    {
      Js.blocks.back().block(0,j) = (double) rand()/RAND_MAX;
      J(i,j+i) = Js.blocks.back().block(0,j);
    }
  }

  // randomly setup the vector
  MatrixNd V(J.rows(), 5);
  for (unsigned j=0; j< 5; j++)
    for (unsigned i=0; i< V.rows(); i++)
      V(i,j) = (double) rand() / RAND_MAX;

  // do the sparse Jacobian
  MatrixNd result_sparse;
  Js.transpose_mult(V, result_sparse);

  // do the dense version
  MatrixNd result_dense;
  J.transpose_mult(V, result_dense);

  // check the results
  for (unsigned i=0; i< result_dense.rows(); i++)
    for (unsigned j=0; j< result_dense.columns(); j++)
      EXPECT_NEAR(result_dense(i,j), result_sparse(i,j), 1e-6);
}

TEST(MultTranspose, Matrix)
{
  MatrixNd J;
  SparseJacobian Js;

  // setup the sparse Jacobian
  J.set_zero(10, 20);
  Js.rows = 10;
  Js.cols = 20;
  for (unsigned i=0; i< 10; i++)
  {
    Js.blocks.push_back(MatrixBlock());
    Js.blocks.back().block.resize(1,10);
    Js.blocks.back().st_row_idx = i;
    Js.blocks.back().st_col_idx = i;
    for (unsigned j=0; j< 10; j++)
    {
      Js.blocks.back().block(0,j) = (double) rand()/RAND_MAX;
      J(i,j+i) = Js.blocks.back().block(0,j);
    }
  }

  // do the sparse Jacobian
  MatrixNd result_sparse;
  Js.mult_transpose(Js, result_sparse);

  // do the dense version
  MatrixNd result_dense;
  J.mult_transpose(J, result_dense);

  // check the results
  for (unsigned i=0; i< result_dense.rows(); i++)
    for (unsigned j=0; j< result_dense.columns(); j++)
      EXPECT_NEAR(result_dense(i,j), result_sparse(i,j), 1e-6);
}

