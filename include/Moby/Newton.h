#ifndef _NEWTON_H_
#define _NEWTON_H_ 

#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>

namespace Moby {

/// Class for doing unconstrained Newton minimization
class Newton
{
  public:
    static void newton(Ravelin::VectorNd& x, unsigned max_iter, double f_tol, void* data, double (*f)(const Ravelin::VectorNd&, void*), void (*grad)(const Ravelin::VectorNd&, void*, Ravelin::VectorNd&), void (*hess)(const Ravelin::VectorNd&, void*, Ravelin::MatrixNd&));
}; // end class def

} // end namespace

#endif

