/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <stdexcept>
#include <Moby/cblas.h>

template <>
void CBLAS::rot(const int N, double *X, const int incX,
                double *Y, const int incY, const double c, const double s)
{
  cblas_drot(N, X, incX, Y, incY, c, s);
}

template <>
void CBLAS::rotg(double& a, double& b, double& c, double& s)
{
  cblas_drotg(&a, &b, &c, &s);
}

template <>
void CBLAS::scal(int N, double alpha, double* X, int incX)
{
  cblas_dscal(N, alpha, X, incX);
}

template <>
void CBLAS::scal(int N, float alpha, float* X, int incX)
{
  cblas_sscal(N, alpha, X, incX);
}

template <>
void CBLAS::copy(int N, const double* X, int incX, double* Y, int incY)
{
  cblas_dcopy(N, X, incX, Y, incY);
}

template <>
void CBLAS::copy(int N, const float* X, int incX, float* Y, int incY)
{
  cblas_scopy(N, X, incX, Y, incY);
}

template <>
void CBLAS::axpy(int N, double alpha, const double* X, int incX, double* Y, int incY)
{
  cblas_daxpy(N, alpha, X, incX, Y, incY);
}

template <>
void CBLAS::axpy(int N, float alpha, const float* X, int incX, float* Y, int incY)
{
  cblas_saxpy(N, alpha, X, incX, Y, incY);
}

template <>
double CBLAS::dot(int N, const double* X, int incX, const double* Y, int incY)
{
  return cblas_ddot(N, X, incX, Y, incY);
}

template <>
void CBLAS::rot(const int N, float *X, const int incX,
                float *Y, const int incY, const float c, const float s)
{
  cblas_srot(N, X, incX, Y, incY, c, s);
}

template <>
void CBLAS::rotg(float& a, float& b, float& c, float& s)
{
  cblas_srotg(&a, &b, &c, &s);
}

template <>
float CBLAS::dot(int N, const float* X, int incX, const float* Y, int incY)
{
  return cblas_sdot(N, X, incX, Y, incY);
}

template <>
double CBLAS::nrm2(int N, const double* X, int incX)
{
  return cblas_dnrm2(N, X, incX);
}

template <>
float CBLAS::nrm2(int N, const float* X, int incX)
{
  return cblas_snrm2(N, X, incX);
}

template <>
void CBLAS::ger(enum CBLAS_ORDER order, int M, int N, double alpha, const double* X, int incX, const double* Y, int incY, double* A, int lda)
{
  cblas_dger(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

template <>
void CBLAS::trsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int m, int n, double alpha, const double* A, int lda, double* B, int ldb)
{
  cblas_dtrsm(CblasColMajor, side, uplo, transA, CblasNonUnit, m, n, alpha, A, lda, B, ldb);  
}

template <>
void CBLAS::trsv(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int n, const double* A, int lda, double* x, int incx)
{
  cblas_dtrsv(CblasColMajor, uplo, transA, CblasNonUnit, n, A, lda, x, incx);  
}

template <>
void CBLAS::ger(enum CBLAS_ORDER order, int M, int N, float alpha, const float* X, int incX, const float* Y, int incY, float* A, int lda)
{
  cblas_sger(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

template <>
void CBLAS::gemv(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N, double alpha, const double* A, int lda, const double* X, int incX, double beta, double* Y, int incY)
{
  #ifndef ADDRESS_ATLAS_BUG
  cblas_dgemv(order, transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  #else
  #define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))
  if (transA == CblasNoTrans)
    cblas_dgemv(order, transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  else
  {
    // NOTE: this code adapted from GSL
    int i, j;
    int lenX, lenY;

    const int Trans = (transA != CblasConjTrans) ? transA : CblasTrans;

    if (M == 0 || N == 0)
      return;

    if (alpha == 0.0 && beta == 1.0)
      return;

    if (Trans == CblasNoTrans) {
      lenX = N;
      lenY = M;
    } else {
      lenX = M;
      lenY = N;
    }

    /* form  y := beta*y */
    if (beta == 0.0) {
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        Y[iy] = 0.0;
        iy += incY;
      }
    } else if (beta != 1.0) {
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        Y[iy] *= beta;
        iy += incY;
      }
    }

    if (alpha == 0.0)
      return;

    if ((order == CblasRowMajor && Trans == CblasNoTrans)
        || (order == CblasColMajor && Trans == CblasTrans)) {
      /* form  y := alpha*A*x + y */
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        double temp = 0.0;
        int ix = OFFSET(lenX, incX);
        for (j = 0; j < lenX; j++) {
          temp += X[ix] * A[lda * i + j];
          ix += incX;
        }
        Y[iy] += alpha * temp;
        iy += incY;
      }
    }   else if ((order == CblasRowMajor && Trans == CblasTrans)
                 || (order == CblasColMajor && Trans == CblasNoTrans)) {
      /* form  y := alpha*A'*x + y */
      int ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        const double temp = alpha * X[ix];
        if (temp != 0.0) {
          int iy = OFFSET(lenY, incY);
          for (i = 0; i < lenY; i++) {
            Y[iy] += temp * A[lda * j + i];
            iy += incY;
          }
        }
        ix += incX;
      }
    } else {
      throw std::runtime_error("unrecognized operation");
    }
  }
  #undef OFFSET
  #endif
}

template <>
void CBLAS::gemv(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N, float alpha, const float* A, int lda, const float* X, int incX, float beta, float* Y, int incY)
{
  cblas_sgemv(order, transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template <>
void CBLAS::gemm(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB, int M, int N, int K, double alpha, const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc)
{
  assert(ldc >= 1 && ldc >= M);
  cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <>
void CBLAS::gemm(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB, int M, int N, int K, float alpha, const float* A, int lda, const float* B, int ldb, float beta, float* C, int ldc)
{
  cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <>
void CBLAS::trsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int m, int n, float alpha, const float* A, int lda, float* B, int ldb)
{
  cblas_strsm(CblasColMajor, side, uplo, transA, CblasNonUnit, m, n, alpha, A, lda, B, ldb);  
}

template <>
void CBLAS::trsv(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int n, const float* A, int lda, float* x, int incx)
{
  cblas_strsv(CblasColMajor, uplo, transA, CblasNonUnit, n, A, lda, x, incx);  
}


// *****************************************************************
// arbitrary precision routines
// *****************************************************************

#ifdef BUILD_ARBITRARY_PRECISION
template <>
mpfr::mpreal CBLAS::nrm2(int N, const mpfr::mpreal* X, int incX)
{
  return cblas_anrm2(N, X, incX);
}

template <>
void CBLAS::ger(enum CBLAS_ORDER order, int M, int N, mpfr::mpreal alpha, const mpfr::mpreal* X, int incX, const mpfr::mpreal* Y, int incY, mpfr::mpreal* A, int lda)
{
  cblas_ager(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

template <>
void CBLAS::gemv(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N, mpfr::mpreal alpha, const mpfr::mpreal* A, int lda, const mpfr::mpreal* X, int incX, mpfr::mpreal beta, mpfr::mpreal* Y, int incY)
{
  cblas_agemv(order, transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template <>
void CBLAS::gemm(enum CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB, int M, int N, int K, mpfr::mpreal alpha, const mpfr::mpreal* A, int lda, const mpfr::mpreal* B, int ldb, mpfr::mpreal beta, mpfr::mpreal* C, int ldc)
{
  cblas_agemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <>
mpfr::mpreal CBLAS::dot(int N, const mpfr::mpreal* X, int incX, const mpfr::mpreal* Y, int incY)
{
  return cblas_adot(N, X, incX, Y, incY);
}

template <>
void CBLAS::axpy(int N, mpfr::mpreal alpha, const mpfr::mpreal* X, int incX, mpfr::mpreal* Y, int incY)
{
  cblas_aaxpy(N, alpha, X, incX, Y, incY);
}

template <>
void CBLAS::scal(int N, mpfr::mpreal alpha, mpfr::mpreal* X, int incX)
{
  cblas_ascal(N, alpha, X, incX);
}

template <>
void CBLAS::copy(int N, const mpfr::mpreal* X, int incX, mpfr::mpreal* Y, int incY)
{
  cblas_acopy(N, X, incX, Y, incY);
}

template <>
static void CBLAS::trsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int m, int n, mpfr::mpreal alpha, const mpfr::mpreal* A, int lda, mpfr::mpreal* B, int ldb)
{
  cblas_atrsm(CblasColMajor, side, uplo, transA, CblasNonUnit, m, n, alpha, A.begin(), lda, B.begin(), ldb);  
}

template <>
static void CBLAS::trsv(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transA, int n, const mpfr::mpreal* A, int lda, mpfr::mpreal* x, int incx)
{
  cblas_atrsv(CblasColMajor, uplo, transA, CblasNonUnit, n, A.begin(), lda, x.begin(), incx);  
}

#endif


