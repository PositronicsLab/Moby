#include <Moby/VectorN.h>
#include <Moby/MatrixNN.h>
#include <Moby/Constants.h>

#ifndef _PATH_LCP_SOLVER_H_
#define _PATH_LCP_SOLVER_H_

namespace Moby {

/// Class for solving LCP problems using PATH's solver
class PathLCPSolver
{
  public:
    static bool solve_lcp(const MatrixNN& MM, const VectorN& q, VectorN& z, Real tol = NEAR_ZERO);
}; // end class PathLCPSolver

} // end namespace

#endif

