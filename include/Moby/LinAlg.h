/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_LINALG_H
#define _MOBY_LINALG_H

#include <Moby/Constants.h>
#include <Moby/MatrixNN.h>
#include <Moby/SparseMatrixN.h>
#include <Moby/VectorN.h>

namespace Moby {

/// Linear algebra routines
/**
LinAlg is a set of static routines that interface to LAPACK.  I have included only very few routines here, however they should be some of the most utilized: SVD, (SVD-based) pseudo-inverse, linear equation solving, and matrix inverse.
*/
class LinAlg
{
  public:
    static MatrixN nullspace(const MatrixN& A, Real tol = -1.0) { MatrixN ns; nullspace(A, ns, tol); return ns; }
    static MatrixN& nullspace(const MatrixN& A, MatrixN& ns, Real tol = -1.0);
    static Real cond(const MatrixN& x);
    static MatrixN to_matrix(const std::string& s);
    static VectorN to_vector(const std::string& s);
    static bool factor_chol(MatrixNN& M);
    static void factor_QR(MatrixN& AR, MatrixN& Q);
    static void factor_QR(MatrixN& AR, MatrixN& Q, std::vector<int>& PI);
    static bool factor_LU(MatrixN& M, std::vector<int>& IPIV);
    static void factor_LDL(MatrixNN& M, std::vector<int>& IPIV);
    static MatrixNN& inverse_LU(MatrixNN& M, const std::vector<int>& IPIV);
    static VectorN solve_chol(const MatrixNN& M, const VectorN& b);
    static MatrixN solve_chol(const MatrixNN& M, const MatrixN& B);
    static VectorN& solve_chol_fast(const MatrixNN& M, VectorN& xb);
    static MatrixN& solve_chol_fast(const MatrixNN& M, MatrixN& XB);
    static VectorN& solve_LDL_fast(const MatrixNN& M, const std::vector<int>& IPIV, VectorN& xb);
    static MatrixN& solve_LDL_fast(const MatrixNN& M, const std::vector<int>& IPIV, MatrixN& XB);
    static VectorN solve_LU(const MatrixNN& M, bool transpose, const std::vector<int>& IPIV, const VectorN& b);
    static VectorN& solve_LU_fast(const MatrixNN& M, bool transpose, const std::vector<int>& IPIV, VectorN& xb);
    static MatrixN solve_LU(const MatrixNN& M, bool transpose, const std::vector<int>& IPIV, const MatrixN& B);
    static MatrixN& solve_LU_fast(const MatrixNN& M, bool transpose, const std::vector<int>& IPIV, MatrixN& XB);
    static MatrixN solve(const MatrixNN& A, const MatrixN& B);
    static VectorN solve(const MatrixNN& A, const VectorN& b);
    static VectorN& solve_tri_fast(const MatrixNN& A, bool utri, bool transpose_A, VectorN& xb);
    static MatrixN& solve_tri_fast(const MatrixNN& A, bool utri, bool transpose_A, MatrixN& XB);
    static VectorN solve_tri(const MatrixNN& A, bool utri, bool transpose_A, const VectorN& b);
    static MatrixN solve_tri(const MatrixNN& A, bool utri, bool transpose_A, const MatrixN& B);
    static VectorN& solve_LS(const SparseMatrixN& A, const VectorN& b, VectorN& x, Real damp = (Real) 0.0, unsigned max_iter = std::numeric_limits<unsigned>::max(), Real tol = (Real) 0.0);
    static MatrixN solve_LS(const MatrixN& A, const MatrixN& B, Real tol = -1.0);
    static VectorN solve_LS(const MatrixN& A, const VectorN& b, Real tol = -1.0);
    static MatrixN solve_symmetric(const MatrixNN& A, const MatrixN& B);
    static VectorN solve_symmetric(const MatrixNN& A, const VectorN& b);
    static MatrixN solve_SPD(const MatrixNN& A, const MatrixN& B);
    static VectorN solve_SPD(const MatrixNN& A, const VectorN& b);
    static MatrixN& solve_fast(MatrixNN& A, MatrixN& XB);
    static VectorN& solve_fast(MatrixNN& A, VectorN& xb);
    static MatrixN& solve_LS_fast(MatrixN& A, MatrixN& XB, Real tol = -1.0);
    static VectorN& solve_LS_fast(MatrixN& A, VectorN& xb, Real tol = -1.0);
    static VectorN& solve_LS_fast2(MatrixN& A, VectorN& xb);
    static MatrixN& solve_symmetric_fast(MatrixNN& A, MatrixN& XB);
    static VectorN& solve_symmetric_fast(MatrixNN& A, VectorN& xb);
    static MatrixN& solve_SPD_fast(MatrixNN& A, MatrixN& XB);
    static VectorN& solve_SPD_fast(MatrixNN& A, VectorN& xb);
    static void solve_iterative(const MatrixNN& A, const VectorN& b, VectorN& x, unsigned iter);
    static void svd(const MatrixN& A, MatrixNN& U, VectorN& S, MatrixNN& V);
    static MatrixNN& inverse(MatrixNN& mat);
    static MatrixNN& inverse_symmetric(MatrixNN& mat);
    static MatrixNN& inverse_chol(MatrixNN& mat);
    static MatrixNN& inverse_SPD(MatrixNN& mat);
    static MatrixNN pseudo_inverse(const MatrixNN& A, Real tol=-1.0);
    static MatrixN pseudo_inverse(const MatrixN& A, Real tol=-1.0);
    static MatrixNN& pseudo_inverse(const MatrixNN& A, MatrixNN& result, Real tol=-1.0);
    static MatrixN& pseudo_inverse(const MatrixN& A, MatrixN& result, Real tol=-1.0);
    static unsigned calc_rank(const MatrixN& x, Real tol=-1.0);
    static void eig_symm(const MatrixNN& A, VectorN& evals);
    static void eig_symm(const MatrixNN& A, VectorN& evals, MatrixNN& evecs);
    static bool is_SPSD(const MatrixNN& A, Real tolerance);
    static bool is_SPD(const MatrixNN& A, Real tolerance);
    static MatrixN& gauss_elim(MatrixN& A, Real tol = (Real) -1.0);
}; // end class

} // end namespace

#endif
