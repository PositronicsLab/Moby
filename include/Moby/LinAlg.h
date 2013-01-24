/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_LINALG_H
#define _MOBY_LINALG_H

#include <Moby/Constants.h>
#include <Moby/MatrixN.h>
#include <Moby/Matrix2.h>
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
    static MatrixN& nullspace(MatrixN& A, MatrixN& ns, Real tol = -1.0);
    static Real cond(MatrixN& A);
    static MatrixN to_matrix(const std::string& s);
    static VectorN to_vector(const std::string& s);
    static bool factor_chol(MatrixN& M);
    static void factor_QR(MatrixN& AR, MatrixN& Q);
    static void factor_QR(MatrixN& AR, MatrixN& Q, std::vector<int>& PI);
    static bool factor_LU(MatrixN& M, std::vector<int>& IPIV);
    static void factor_LDL(MatrixN& M, std::vector<int>& IPIV);
    static MatrixN& inverse_LU(MatrixN& M, const std::vector<int>& IPIV);
    static VectorN solve_chol(const MatrixN& M, const VectorN& b);
    static MatrixN solve_chol(const MatrixN& M, const MatrixN& B);
    static VectorN& solve_chol_fast(const MatrixN& M, VectorN& xb);
    static MatrixN& solve_chol_fast(const MatrixN& M, MatrixN& XB);
    static VectorN& solve_LDL_fast(const MatrixN& M, const std::vector<int>& IPIV, VectorN& xb);
    static MatrixN& solve_LDL_fast(const MatrixN& M, const std::vector<int>& IPIV, MatrixN& XB);
    static VectorN solve_LU(const MatrixN& M, bool transpose, const std::vector<int>& IPIV, const VectorN& b);
    static VectorN& solve_LU_fast(const MatrixN& M, bool transpose, const std::vector<int>& IPIV, VectorN& xb);
    static MatrixN solve_LU(const MatrixN& M, bool transpose, const std::vector<int>& IPIV, const MatrixN& B);
    static MatrixN& solve_LU_fast(const MatrixN& M, bool transpose, const std::vector<int>& IPIV, MatrixN& XB);
    static MatrixN solve(const MatrixN& A, const MatrixN& B);
    static VectorN solve(const MatrixN& A, const VectorN& b);
    static VectorN& solve_tridiagonal_fast(VectorN& dlo, VectorN& diag, VectorN& dup, VectorN& xb); 
    static MatrixN& solve_tridiagonal_fast(VectorN& dlo, VectorN& diag, VectorN& dup, MatrixN& XB); 
    static VectorN& solve_tri_fast(const MatrixN& A, bool utri, bool transpose_A, VectorN& xb);
    static MatrixN& solve_tri_fast(const MatrixN& A, bool utri, bool transpose_A, MatrixN& XB);
    static VectorN& solve_LS(const SparseMatrixN& A, const VectorN& b, VectorN& x, Real damp = (Real) 0.0, unsigned max_iter = std::numeric_limits<unsigned>::max(), Real tol = (Real) 0.0);
    static MatrixN solve_LS(const MatrixN& A, const MatrixN& B, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol = -1.0);
    static VectorN solve_LS(const MatrixN& A, const VectorN& b, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol = -1.0);
    static MatrixN solve_symmetric(const MatrixN& A, const MatrixN& B);
    static VectorN solve_symmetric(const MatrixN& A, const VectorN& b);
    static MatrixN solve_SPD(const MatrixN& A, const MatrixN& B);
    static VectorN solve_SPD(const MatrixN& A, const VectorN& b);
    static MatrixN& solve_fast(MatrixN& A, MatrixN& XB);
    static VectorN& solve_fast(MatrixN& A, VectorN& xb);
    static MatrixN& solve_LS_fast1(MatrixN& A, MatrixN& XB, Real tol = -1.0) { solve_LS_fast(A, XB, svd1, tol); }
    static VectorN& solve_LS_fast1(MatrixN& A, VectorN& xb, Real tol = -1.0) { solve_LS_fast(A, xb, svd1, tol); }
    static MatrixN& solve_LS_fast2(MatrixN& A, MatrixN& XB, Real tol = -1.0) { solve_LS_fast(A, XB, svd2, tol); }
    static VectorN& solve_LS_fast2(MatrixN& A, VectorN& xb, Real tol = -1.0) { solve_LS_fast(A, xb, svd2, tol); }
    static MatrixN& solve_LS_fast(MatrixN& A, MatrixN& XB, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol = -1.0);
    static VectorN& solve_LS_fast(MatrixN& A, VectorN& xb, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol = -1.0);
    static MatrixN& solve_symmetric_fast(MatrixN& A, MatrixN& XB);
    static VectorN& solve_symmetric_fast(MatrixN& A, VectorN& xb);
    static MatrixN& solve_SPD_fast(MatrixN& A, MatrixN& XB);
    static VectorN& solve_SPD_fast(MatrixN& A, VectorN& xb);
    static void solve_iterative(const MatrixN& A, const VectorN& b, VectorN& x, unsigned iter);
    static void svd(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V);
    static void svd1(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V);
    static void svd2(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V);
    static MatrixN& inverse(MatrixN& mat);
    static MatrixN& inverse_symmetric(MatrixN& mat);
    static MatrixN& inverse_chol(MatrixN& mat);
    static MatrixN& inverse_SPD(MatrixN& mat);
    static MatrixN& pseudo_inverse(MatrixN& A, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol=-1.0);
    static unsigned calc_rank(MatrixN& x, Real tol=-1.0);
    static void eig_symm(MatrixN& A, VectorN& evals);
    static void eig_symm_plus(MatrixN& A_evecs, VectorN& evals);
    static bool is_SPSD(const MatrixN& A, Real tolerance);
    static bool is_SPD(const MatrixN& A, Real tolerance);
    static MatrixN& gauss_elim(MatrixN& A, Real tol = (Real) -1.0);
    static void givens(Real a, Real b, Real& c, Real& s);
    static Matrix2 givens(Real c, Real s);
    static void householder(Real alpha, const VectorN& x, Real& tau, VectorN& v);
    static void update_QR_rank1(MatrixN& Q, MatrixN& R, const VectorN& u, const VectorN& v);
    static void update_QR_delete_cols(MatrixN& Q, MatrixN& R, unsigned k, unsigned p);
    static void update_QR_insert_cols(MatrixN& Q, MatrixN& R, MatrixN& U, unsigned k);
    static void update_QR_insert_rows(MatrixN& Q, MatrixN& R, MatrixN& U, unsigned k);
    static void update_QR_delete_rows(MatrixN& Q, MatrixN& R, unsigned k, unsigned p);
}; // end class

} // end namespace

#endif
