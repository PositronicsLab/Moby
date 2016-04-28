#include <iostream>
#include <Moby/qpOASES.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>

int main( void ) {
  Moby::QPOASES _qp;

  const double INF = 1e+29;
  
  // Setup H matrix
  Ravelin::MatrixNd H;
  H.resize(2, 2);
  H(0, 0) = 4;
  H(0, 1) = 1;
  H(1, 0) = 1;
  H(1, 1) = 2;

  // linear coefs
  Ravelin::VectorNd p;
  p.resize(2);
  p[0] = 1;
  p[1] = 1;

  // Inequality constraints
  Ravelin::MatrixNd M;
  M.resize(2, 2);
  M(0, 0) = 1;
  M(0, 1) = 0;
  M(1, 1) = 1;
  M(1, 0) = 0;

  Ravelin::VectorNd q;
  q.resize(2);
  q[0] = 0;
  q[1] = 0;

  // Setup equality constraint
  Ravelin::MatrixNd A;
  A.resize(1, 2);
  A(0, 0) = 1;
  A(0, 1) = 1;

  Ravelin::VectorNd b;
  b.resize(1);
  b[0] = 1;

  // Setup the x bounds
  Ravelin::VectorNd lb, ub;
  lb.resize(2);
  ub.resize(2);
  lb[0] = 0;
  lb[1] = 0;
  ub[0] = INF;
  ub[1] = INF;

  // Output
  Ravelin::VectorNd z;
  z.resize(2);
  z[0] = 0;
  z[1] = 0;

  if (!_qp.qp_activeset(H, p, lb, ub, M, q, A, b, z)){
    std::cout << "QP failed to find feasible point" << std::endl;
  } else {
    std::cout << "QP solved successfully: " << z << std::endl;
  }
  return 0;
}

