#include <Ravelin/LinAlgd.h>
#include <Moby/Newton.h>

using namespace Ravelin;
using namespace Moby;

/// Does Newton's method for unconstrained minimization
void Newton::newton(VectorNd& x, unsigned max_iter, double f_tol, void* data, double (*f)(const Ravelin::VectorNd&, void*), void (*grad)(const Ravelin::VectorNd&, void*, Ravelin::VectorNd&), void (*hess)(const Ravelin::VectorNd&, void*, Ravelin::MatrixNd&))
{
  const double ALPHA = 0.05;  // parameter for backtracking line search
  const double BETA = 0.5;  // reduction parameter for backtracking line search
  MatrixNd H, cH;
  VectorNd g, dx, xprime;

  // compute f(x)
  double f_naught = f(x, data);

  // loop up to the specified number of iterations
  for (unsigned i=0; i< max_iter; i++)
  {
    // compute the gradient
    grad(x, data, g);

    // compute the Hessian
    hess(x, data, H);

    // factor the Hessian
    const double REG_UPDATE = 10.0;
    double reg = 1e-16;
    cH = H;
    while (true)
    {
      // attempt to factor
      if (LinAlgd::factor_chol(cH))
        break;
      
      // factorization failed; add regularization
      cH = H;
      RowIteratord ri = cH.row_iterator_begin();
      for (unsigned j=0; j< H.rows(); j++, ri += H.rows()+1)
        *ri += reg;

      // update regularization factor
      reg *= REG_UPDATE;
    }

    // solve for dx
    dx = g;
    dx.negate();
    LinAlgd::solve_chol_fast(cH, dx);

    // prepare to do backtracking line search
    const double grad_dot_dx = g.dot(dx);
    double s = 1.0;
    ((xprime = dx) *= s) += x;

    // do backtracking line search
    double fprime = f(xprime, data);
    while (fprime > f_naught + ALPHA * s * grad_dot_dx)
    {
      // update s and compute f'(x')  
      s *= BETA;
      ((xprime = dx) *= s) += x;
      fprime = f(xprime, data);
    }

    // see whether the change in f is sufficiently small to quit
    if (f_naught - f_tol < fprime)
      return;

    // update x
    x = xprime;

    // update f_naught
    f_naught = fprime;
  }   
}

