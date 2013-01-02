/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/algorithm/minmax.hpp>
#include <cstring>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <boost/algorithm/minmax.hpp>
#include <Moby/FastThreadable.h>
#include <Moby/cblas.h>
#include "clapack.h"
#include <Moby/MissizeException.h>
#include <Moby/NonsquareMatrixException.h>
#include <Moby/NumericalException.h>
#include <Moby/SingularException.h>
#include <Moby/Log.h>
#include <Moby/LinAlg.h>

using namespace Moby;
using std::vector;
using boost::shared_array;
using std::endl;

#ifdef BUILD_SINGLE
// ***********************************************************************
// single precision routines start here
// ***********************************************************************

float log2(float x)
{
  return std::log(x)/std::log(2.0f);
}

/// Calls LAPACK function for Givens rotation
static void lartg_(float* F, float* G, float* CS, float* SN, float* R)
{
  slartg_(F, G, CS, SN, R);
}

/// Calls LAPACK function for solving tridiagonal systems of linear equations
static void gtsv_(INTEGER* N, INTEGER* NRHS, SINGLE* DL, SINGLE* D, SINGLE* DU, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sgtsv_(N, NRHS, DL, D, DU, B, LDB, INFO);
}

/// Calls LAPACK function for solving triangular systems of linear equations
static void trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, SINGLE* AP, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  char DIAG = 'N';
  strtrs_(UPLO, TRANS, &DIAG, N, NRHS, AP, LDA, B, LDB, INFO);
}

/// Calls LAPACK function *ormqr_
static void ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, SINGLE* A, INTEGER* LDA, SINGLE* TAU, SINGLE* C, INTEGER* LDC, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  SINGLE WORK_QUERY;
  sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK_QUERY, &LWORK, INFO);
  assert(INFO == 0);

  // declare memory
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<SINGLE> > WORK;
  WORK().resize(LWORK);

  // do the real call now
  sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK()[0], &LWORK, INFO);
  assert(INFO == 0);
}

/// Calls LAPACK function for doing LDL' factorization of a matrix
static void sptrf_(char* UPLO, INTEGER* N, SINGLE* AP, INTEGER* IPIV, INTEGER* INFO)
{
  ssptrf_(UPLO, N, AP, IPIV, INFO);
}

/// Calls LAPACK function for solution to Ax=b using LDL' factorization
static void sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* AP, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  ssptrs_(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for least-squares solution to Ax=b using SVD
static void gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, SINGLE* RCOND, INTEGER* INFO)
{
  // create array to hold singular values
  INTEGER min_mn = std::min(*M, *N);
  SAFESTATIC FastThreadable<vector<float> > S;
  S().resize(min_mn);

   // necessary for computing WORK sizes
  INTEGER ISPEC = (INTEGER) 9;
  INTEGER TMP = (INTEGER) 0;
  const char* NAME = "SGELSD";
  const char* OPTS = " ";
  INTEGER smlsiz = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, NRHS, &TMP, strlen(NAME), strlen(OPTS));
  assert(smlsiz > (INTEGER) 0);
  INTEGER NLVL = std::max((INTEGER) 0, (INTEGER) (log2((float) min_mn/(float) (smlsiz+1))+1));
  INTEGER LIWORK = std::max((INTEGER) 1, 3*min_mn*NLVL + 11*min_mn);
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  IWORK().resize(LIWORK);
  INTEGER RANK;

  // first do WORK query
  float WORK_QUERY;
  INTEGER LWORK = -1;
  sgelsd_(M, N, NRHS, A, M, B, LDB, &S().front(), RCOND, &RANK, &WORK_QUERY, &LWORK, &IWORK().front(), INFO);
  assert(*INFO == 0);

  // setup LWORK 
  LWORK = (INTEGER) WORK_QUERY;

  // setup WORK array
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(std::max((INTEGER) 1, LWORK));

  // compute
  sgelsd_(M, N, NRHS, A, M, B, LDB, &S.front(), RCOND, &RANK, &WORK().front(), &LWORK, &IWORK().front(), INFO);
}

/// Calls appropriate LAPACK function for solving systems of linear equations Ax=b, where A is symmetric indefinite
static void sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{ 
  INTEGER LWORK = -1;
  float WORK_QUERY;

  // first, determine workspace size
  ssysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return;

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK);

  // call LAPACK
  ssysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for svd (divide and conquer)
static void gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* S, SINGLE* U, INTEGER* LDU, SINGLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  float WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  IWORK().resize(8*minmn);

  // call LAPACK to determine the optimal workspace size
  sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, IWORK(), INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK); 

  // call LAPACK once again
  sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK.front(), &LWORK, &IWORK().front(), INFO);
}

/// Calls LAPACK function for svd 
static void gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* S, SINGLE* U, INTEGER* LDU, SINGLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  float WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;

  // call LAPACK to determine the optimal workspace size
  sgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK); 

  // call LAPACK once again
  sgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK.front(), &LWORK, INFO);
}

/// Calls LAPACK function for computing eigenvalues and eigenvectors
static void syevd_(char* JOBZ, char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* EVALS, INTEGER* INFO)
{
  // do work query first
  float WORK_QUERY;
  INTEGER LWORK = -1;
  INTEGER IWORK_QUERY;
  INTEGER LIWORK = -1;
  ssyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK_QUERY, &LWORK, &IWORK_QUERY, &LIWORK, INFO);

  // set array sizes
  LWORK = (INTEGER) WORK_QUERY;
  LIWORK = IWORK_QUERY;
  SAFESTATIC FastThreadable<vector<float> > WORK;
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  WORK().resize(LWORK);
  IWORK().resize(LIWORK);
  ssyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK().front(), &LWORK, &IWORK().front(), &LIWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations using LU factorization
static void gesv_(INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* X, INTEGER* LDX, INTEGER* INFO)
{
  sgesv_(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

/// Calls LAPACK function for solving a system of linear equations from a Cholesky factorization
static void potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  spotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
}

/// Calls LAPACK function for computing matrix inverse using Cholesky factorization
static void potri_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* INFO)
{
  spotri_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for Cholesky factorization
static void potrf_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* INFO)
{
  spotrf_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for solving system of equations Ax=b, where A is PSD
static void posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sposv_(UPLO, N, NRHS, A, N, B, LDB, INFO);
}

/// Calls LAPACK function for inverting symmetric indefinite matrix using factorization
static void sytri_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(*N);
  ssytri_(UPLO, N, A, LDA, IPIV, &WORK().front(), INFO);
}

/// Calls LAPACK function for factorizing symmetric, indefinite matrix 
static void sytrf_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO) 
{
  // perform workspace query for factorization
  float WORK_QUERY;
  INTEGER LWORK = -1;
  ssytrf_(UPLO, N, A, LDA, IPIV, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup WORK array
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK);

  // perform the necessary factorization
  ssytrf_(UPLO, N, A, LDA, IPIV, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
static void geqp3_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* JPVT, SINGLE* TAU, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  SINGLE work_query;
  sgeqp3_(M, N, A, LDA, JPVT, TAU, &work_query, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) work_query;
  SAFESTATIC FastThreadable<vector<SINGLE> > WORK;
  WORK().resize(LWORK);

  // do QR factorization
  sgeqp3_(M, N, A, LDA, JPVT, TAU, &WORK()[0], &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
static void geqrf_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* TAU, INTEGER* INFO)
{
  // determine LWORK
  const char* NAME = "SGEQRF";
  const char* OPTS = " ";
  INTEGER ISPEC = 1L;
  INTEGER TMP = -1;
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, &TMP, &TMP, strlen(NAME), strlen(OPTS));
  assert(NB >= 0);

  // setup WORK vectors
  INTEGER LWORK = NB*(*N);
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK);
  sgeqrf_(M, N, A, LDA, TAU, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for LU factorization 
static void getrf_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  sgetrf_(M, N, A, LDA, IPIV, INFO);
}

/// Calls LAPACK function for matrix inversion using LU factorization
static void getri_(INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{ 
  // compute block size
  INTEGER TMP = -1;
  INTEGER ISPEC = 1;
  const char* NAME = "SGETRI";
  const char* OPTS = " ";
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, N, &TMP, &TMP, &TMP, strlen(NAME), strlen(OPTS));

  // setup LAPACK parameters
  INTEGER LWORK = (*N)*NB;

  // setup work vector
  SAFESTATIC FastThreadable<vector<float> > WORK;
  WORK().resize(LWORK);
  
  sgetri_(N, A, LDA, IPIV, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations with a given LU factorization
static void getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sgetrs_(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for forming Q from a QR factorization
static void orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, SINGLE* A, INTEGER* LDA, SINGLE* TAU, INTEGER* INFO)
{
  // do a workspace query
  INTEGER LWORK = -1;
  SINGLE WORK_QUERY;
  sorgqr_(M, N, K, A, LDA, TAU, &WORK_QUERY, &LWORK, INFO);

  // initialize the work array
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<SINGLE> > WORK;
  WORK().resize(LWORK);

  // call the function for real
  sorgqr_(M, N, K, A, LDA, TAU, &WORK().front(), &LWORK, INFO);
}

#else
#ifdef BUILD_DOUBLE
// ***********************************************************************
// double precision routines start here
// ***********************************************************************

double log2(double x)
{
  return std::log(x)/std::log(2.0);
}

/// Calls LAPACK function for Givens rotation
static void lartg_(double* F, double* G, double* CS, double* SN, double* R)
{
  dlartg_(F, G, CS, SN, R);
}

/// Calls LAPACK function for solving tridiagonal systems of linear equations
static void gtsv_(INTEGER* N, INTEGER* NRHS, DOUBLE* DL, DOUBLE* D, DOUBLE* DU, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dgtsv_(N, NRHS, DL, D, DU, B, LDB, INFO);
}

/// Calls LAPACK function for solving triangular systems of linear equations
static void trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, DOUBLE* AP, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  char DIAG = 'N';
  dtrtrs_(UPLO, TRANS, &DIAG, N, NRHS, AP, LDA, B, LDB, INFO);
}

/// Calls LAPACK function *ormqr_
static void ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, DOUBLE* C, INTEGER* LDC, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  DOUBLE WORK_QUERY;
  dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO == 0);

  // declare memory
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);

  // do the real call now
  dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK()[0], &LWORK, INFO);
  assert(*INFO == 0);
}

/// Calls LAPACK function for doing LDL' factorization of a matrix
static void sptrf_(char* UPLO, INTEGER* N, DOUBLE* AP, INTEGER* IPIV, INTEGER* INFO)
{
  dsptrf_(UPLO, N, AP, IPIV, INFO);
}

/// Calls LAPACK function for solution to Ax=b using LDL' factorization
static void sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* AP, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dsptrs_(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for least-squares solution to Ax=b using SVD
static void gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, DOUBLE* RCOND, INTEGER* INFO)
{
  // create array to hold singular values
  INTEGER min_mn = std::min(*M, *N);
  SAFESTATIC FastThreadable<vector<double> > S;
  S().resize(min_mn);

   // necessary for computing WORK sizes
  INTEGER ISPEC = (INTEGER) 9;
  INTEGER TMP = (INTEGER) 0;
  const char* NAME = "DGELSD";
  const char* OPTS = " ";
  INTEGER smlsiz = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, NRHS, &TMP, strlen(NAME), strlen(OPTS));
  assert(smlsiz > (INTEGER) 0);
  INTEGER NLVL = std::max((INTEGER) 0, (INTEGER) (log2((double) min_mn/(double) (smlsiz+1))+1));
  INTEGER LIWORK = std::max((INTEGER) 1, 3*min_mn*NLVL + 11*min_mn);
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  IWORK().resize(LIWORK);
  INTEGER RANK;

  // first do WORK query
  double WORK_QUERY;
  INTEGER LWORK = -1;
  dgelsd_(M, N, NRHS, A, M, B, LDB, &S().front(), RCOND, &RANK, &WORK_QUERY, &LWORK, &LIWORK, INFO);
  assert(*INFO == 0);

  // setup LWORK 
  LWORK = (INTEGER) WORK_QUERY;

  // setup WORK array
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(std::max((INTEGER) 1, LWORK));

  // compute
  dgelsd_(M, N, NRHS, A, M, B, LDB, &S().front(), RCOND, &RANK, &WORK().front(), &LWORK, &IWORK().front(), INFO);
}

/// Calls appropriate LAPACK function for solving systems of linear equations Ax=b, where A is symmetric indefinite
static void sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{ 
  INTEGER LWORK = -1;
  double WORK_QUERY;

  // first, determine workspace size
  dsysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDA, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);

  // call LAPACK
  dsysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for svd (divide and conquer)
static void gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* S, DOUBLE* U, INTEGER* LDU, DOUBLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  double WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  IWORK().resize(8*minmn);

  // call LAPACK to determine the optimal workspace size
  dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, &IWORK().front(), INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK); 

  // call LAPACK once again
  dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK().front(), &LWORK, &IWORK().front(), INFO);
}

/// Calls LAPACK function for svd
static void gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* S, DOUBLE* U, INTEGER* LDU, DOUBLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  double WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;

  // call LAPACK to determine the optimal workspace size
  dgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK); 

  // call LAPACK once again
  dgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for computing eigenvalues and eigenvectors
static void syevd_(char* JOBZ, char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* EVALS, INTEGER* INFO)
{
  // do work query first
  double WORK_QUERY;
  INTEGER LWORK = -1;
  INTEGER IWORK_QUERY;
  INTEGER LIWORK = -1;
  dsyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK_QUERY, &LWORK, &IWORK_QUERY, &LIWORK, INFO);

  // set array sizes
  LWORK = (INTEGER) WORK_QUERY;
  LIWORK = IWORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  SAFESTATIC FastThreadable<vector<INTEGER> > IWORK;
  WORK().resize(LWORK);
  IWORK().resize(LIWORK);

  dsyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK().front(), &LWORK, &IWORK().front(), &LIWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations using LU factorization
static void gesv_(INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* X, INTEGER* LDX, INTEGER* INFO)
{
  dgesv_(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

/// Calls LAPACK function for solving a system of linear equation from a Cholesky factorization
static void potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dpotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
}

/// Calls LAPACK function for computing matrix inverse using Cholesky factorization
static void potri_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* INFO)
{
  dpotri_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for Cholesky factorization
static void potrf_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* INFO)
{
  dpotrf_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for solving system of equations Ax=b, where A is PSD
static void posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dposv_(UPLO, N, NRHS, A, N, B, LDB, INFO);
}

/// Calls LAPACK function for inverting symmetric indefinite matrix using factorization
static void sytri_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(*N);
  dsytri_(UPLO, N, A, LDA, IPIV, &WORK().front(), INFO);
}

/// Calls LAPACK function for factorizing symmetric, indefinite matrix 
static void sytrf_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER *INFO) 
{
  // perform workspace query for factorization
  double WORK_QUERY;
  INTEGER LWORK = -1;
  dsytrf_(UPLO, N, A, LDA, IPIV, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup WORK array
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);

  // perform the necessary factorization
  dsytrf_(UPLO, N, A, LDA, IPIV, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
static void geqp3_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* JPVT, DOUBLE* TAU, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  DOUBLE work_query;
  dgeqp3_(M, N, A, LDA, JPVT, TAU, &work_query, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) work_query;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);

  // do QR factorization
  dgeqp3_(M, N, A, LDA, JPVT, TAU, &WORK()[0], &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
static void geqrf_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, INTEGER* INFO)
{
  // determine LWORK
  const char* NAME = "DGEQRF";
  const char* OPTS = " ";
  INTEGER ISPEC = 1L;
  INTEGER TMP = -1;
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, &TMP, &TMP, strlen(NAME), strlen(OPTS));
  assert(NB >= 0);

  // setup WORK vectors
  INTEGER LWORK = NB*(*N);
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);
  dgeqrf_(M, N, A, LDA, TAU, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for LU factorization 
static void getrf_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  dgetrf_(M, N, A, LDA, IPIV, INFO);
}

/// Calls LAPACK function for matrix inversion using LU factorization
static void getri_(INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{ 
  // compute block size
  INTEGER TMP = -1;
  INTEGER ISPEC = 1;
  const char* NAME = "DGETRI";
  const char* OPTS = " ";
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, N, &TMP, &TMP, &TMP, strlen(NAME), strlen(OPTS));

  // setup LAPACK parameters
  INTEGER LWORK = (*N)*NB;

  // setup work vector 
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);
  
  dgetri_(N, A, LDA, IPIV, &WORK().front(), &LWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations with a given LU factorization
static void getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dgetrs_(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for forming Q from a QR factorization
static void orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, INTEGER* INFO)
{
  // do a workspace query
  INTEGER LWORK = -1;
  DOUBLE WORK_QUERY;
  dorgqr_(M, N, K, A, LDA, TAU, &WORK_QUERY, &LWORK, INFO);

  // initialize the work array
  LWORK = (INTEGER) WORK_QUERY;
  SAFESTATIC FastThreadable<vector<DOUBLE> > WORK;
  WORK().resize(LWORK);

  // call the function for real
  dorgqr_(M, N, K, A, LDA, TAU, &WORK().front(), &LWORK, INFO);
}
#else
#ifdef BUILD_ARBITRARY_PRECISION
// ***********************************************************************
// artibrary precision routines start here
// ***********************************************************************

/// Calls LAPACK function for Givens rotation
static void lartg_(mpfr::mpreal* F, mpfr::mpreal* G, mpfr::mpreal* CS, mpfr::mpreal* SN, mpfr::mpreal* R)
{
  alartg_(F, G, CS, SN, R);
}

/// Calls LAPACK function for solving tridiagonal systems of linear equations
static void gtsv_(INTEGER* N, INTEGER* NRHS, mpfr::mpreal* DL, mpfr::mpreal* D, mpfr::mpreal* DU, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  agtsv_(N, NRHS, DL, D, DU, B, LDB, INFO);
}

/// Calls LAPACK function for solving triangular systems of linear equations
static void trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* AP, INTEGER* LDA, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  char DIAG = 'N';
  atrtrs_(UPLO, TRANS, &DIAG, N, NRHS, AP, LDA, B, LDB, INFO);
}

/// Calls LAPACK function *ormqr_
static void ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* TAU, mpfr::mpreal* C, INTEGER* LDC, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  mpfr::mpreal WORK_QUERY;
  aormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO == 0);

  // declare memory
  LWORK = (INTEGER) WORK_QUERY;
  vector<mpfr::mpreal> WORK(LWORK);

  // do the real call now
  aormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK[0], &LWORK, INFO);
  assert(*INFO == 0);
}

/// Calls LAPACK function for doing LDL' factorization of a matrix
static void sptrf_(char* UPLO, INTEGER* N, mpfr::mpreal* AP, INTEGER* IPIV, INTEGER* INFO)
{
  asptrf_(UPLO, N, AP, IPIV, INFO);
}

/// Calls LAPACK function for solution to Ax=b using LDL' factorization
static void sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* AP, INTEGER* IPIV, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  asptrs_(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for least-squares solution to Ax=b using SVD
static void gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* B, INTEGER* LDB, mpfr::mpreal* RCOND, INTEGER* INFO)
{
  // create array to hold singular values
  INTEGER min_mn = std::min(*M, *N);
  mpfr::mpreal* S = new mpfr::mpreal[min_mn];

   // necessary for computing WORK sizes
  INTEGER ISPEC = (INTEGER) 9;
  INTEGER TMP = (INTEGER) 0;
  const char* NAME = "DGELSD";
  const char* OPTS = " ";
  INTEGER smlsiz = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, NRHS, &TMP, strlen(NAME), strlen(OPTS));
  assert(smlsiz > (INTEGER) 0);
  const mpfr::mpreal l2c = (mpfr::mpreal) 1.0/std::log((mpfr::mpreal) 2.0);
  INTEGER NLVL = std::max((INTEGER) 0, (INTEGER) (std::log(min_mn/(smlsiz+1))*l2c+1));
  INTEGER LIWORK = 3*min_mn*NLVL + 11*min_mn;
  INTEGER* IWORK = new INTEGER[std::max((INTEGER) 1,LIWORK)];
  INTEGER RANK;

  // first do WORK query
  mpfr::mpreal WORK_QUERY;
  INTEGER LWORK = -1;
  agelsd_(M, N, NRHS, A, M, B, LDB, S, RCOND, &RANK, &WORK_QUERY, &LWORK, IWORK, INFO);
  assert(*INFO == 0);

  // setup LWORK 
  LWORK = (INTEGER) WORK_QUERY;

  // setup WORK array
  mpfr::mpreal* WORK = new mpfr::mpreal[std::max((INTEGER) 1, LWORK)];

  // compute
  agelsd_(M, N, NRHS, A, M, B, LDB, S, RCOND, &RANK, WORK, &LWORK, IWORK, INFO);

  delete [] IWORK;
  delete [] S;
  delete [] WORK;
}

/// Calls appropriate LAPACK function for solving systems of linear equations Ax=b, where A is symmetric indefinite
static void sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{ 
  INTEGER LWORK = -1;
  mpfr::mpreal WORK_QUERY;

  // first, determine workspace size
  asysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDA, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK];

  // call LAPACK
  asysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &LWORK, INFO);

  delete [] WORK;
}

/// Calls LAPACK function for svd (standard)
static void gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* S, mpfr::mpreal* U, INTEGER* LDU, mpfr::mpreal* V, INTEGER* LDVT, INTEGER* INFO)
{
  mpfr::mpreal WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;

  // call LAPACK to determine the optimal workspace size
  agesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK];

  // call LAPACK once again
  agesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, WORK, &LWORK, INFO);

  delete [] IWORK;
  delete [] WORK;
}

/// Calls LAPACK function for svd (divide and conquer)
static void gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* S, mpfr::mpreal* U, INTEGER* LDU, mpfr::mpreal* V, INTEGER* LDVT, INTEGER* INFO)
{
  mpfr::mpreal WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;
  INTEGER* IWORK = new INTEGER[8*minmn];

  // call LAPACK to determine the optimal workspace size
  agesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, IWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK];

  // call LAPACK once again
  agesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, WORK, &LWORK, IWORK, INFO);

  delete [] IWORK;
  delete [] WORK;
}

/// Calls LAPACK function for computing eigenvalues and eigenvectors
static void syevd_(char* JOBZ, char* UPLO, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* EVALS, INTEGER* INFO)
{
  // do work query first
  mpfr::mpreal WORK_QUERY;
  INTEGER LWORK = -1;
  INTEGER IWORK_QUERY;
  INTEGER LIWORK = -1;
  asyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK_QUERY, &LWORK, &IWORK_QUERY, &LIWORK, INFO);

  // set array sizes
  LWORK = (INTEGER) WORK_QUERY;
  LIWORK = IWORK_QUERY;
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK];
  INTEGER* IWORK = new INTEGER[LIWORK];
  asyevd_(JOBZ, UPLO, N, A, LDA, EVALS, WORK, &LWORK, IWORK, &LIWORK, INFO);

  delete [] WORK;
  delete [] IWORK;
}

/// Calls LAPACK function for solving system of linear equations using LU factorization
static void gesv_(INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, mpfr::mpreal* X, INTEGER* LDX, INTEGER* INFO)
{
  agesv_(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

/// Calls LAPACK function for solving a system of linear equation from a Cholesky factorization
static void potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  apotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
}

/// Calls LAPACK function for computing matrix inverse using Cholesky factorization
static void potri_(char* UPLO, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* INFO)
{
  apotri_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for Cholesky factorization
static void potrf_(char* UPLO, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* INFO)
{
  apotrf_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for solving system of equations Ax=b, where A is PSD
static void posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  aposv_(UPLO, N, NRHS, A, N, B, LDB, INFO);
}

/// Calls LAPACK function for inverting symmetric indefinite matrix using factorization
static void sytri_(char* UPLO, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  mpfr::mpreal* WORK = new mpfr::mpreal[*N];
  asytri_(UPLO, N, A, LDA, IPIV, WORK, INFO);
  delete [] WORK;
}

/// Calls LAPACK function for factorizing symmetric, indefinite matrix 
static void sytrf_(char* UPLO, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, INTEGER *INFO) 
{
  // perform workspace query for factorization
  mpfr::mpreal WORK_QUERY;
  INTEGER LWORK = -1;
  asytrf_(UPLO, N, A, LDA, IPIV, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup WORK array
  LWORK = (INTEGER) WORK_QUERY;
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK]; 

  // perform the necessary factorization
  asytrf_(UPLO, N, A, LDA, IPIV, WORK, &LWORK, INFO);
  delete [] WORK;
}

/// Calls LAPACK function for QR factorization
static void geqrf_(INTEGER* M, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* TAU, INTEGER* INFO)
{
  // determine LWORK
  const char* NAME = "DGEQRF";
  const char* OPTS = " ";
  INTEGER ISPEC = 1L;
  INTEGER TMP = -1;
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, &TMP, &TMP, strlen(NAME), strlen(OPTS));
  assert(NB >= 0);

  // setup WORK vectors
  INTEGER LWORK = NB*(*N);
  vector<mpfr::mpreal> WORK(LWORK);
  ageqrf_(M, N, A, LDA, TAU, &WORK.front(), &LWORK, INFO);
}

/// Calls LAPACK function for LU factorization 
static void getrf_(INTEGER* M, INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  agetrf_(M, N, A, LDA, IPIV, INFO);
}

/// Calls LAPACK function for matrix inversion using LU factorization
static void getri_(INTEGER* N, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{ 
  // compute block size
  INTEGER TMP = -1;
  INTEGER ISPEC = 1;
  const char* NAME = "DGETRI";
  const char* OPTS = " ";
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, N, &TMP, &TMP, &TMP, strlen(NAME), strlen(OPTS));

  // setup LAPACK parameters
  INTEGER LWORK = (*N)*NB;
  
  mpfr::mpreal* WORK = new mpfr::mpreal[LWORK];
  agetri_(N, A, LDA, IPIV, WORK, &LWORK, INFO);
  delete [] WORK;
}

/// Calls LAPACK function for solving system of linear equations with a given LU factorization
static void getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, mpfr::mpreal* A, INTEGER* LDA, INTEGER* IPIV, mpfr::mpreal* B, INTEGER* LDB, INTEGER* INFO)
{
  agetrs_(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for forming Q from a QR factorization
static void orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, mpfr::mpreal* A, INTEGER* LDA, mpfr::mpreal* TAU, INTEGER* INFO)
{
  // do a workspace query
  INTEGER LWORK = -1;
  mpfr::mpreal WORK_QUERY;
  aorgqr_(M, N, K, A, LDA, TAU, &WORK_QUERY, &LWORK, INFO);

  // initialize the work array
  LWORK = (INTEGER) WORK_QUERY;
  shared_array<mpfr::mpreal> WORK(new mpfr::mpreal[LWORK]);

  // call the function for real
  aorgqr_(M, N, K, A, LDA, TAU, WORK.get(), &LWORK, INFO);
}
#endif
#endif
#endif

/// Solves a tridiagonal system
/**
 * \param dl the (n-1) elements on the subdiagonal (destroyed on return)
 * \param d  the n elements on the diagonal (destroyed on return)
 * \param du the n elements on the superdiagonal (destroyed on return)
 * \param xb the right hand side on input, the solution on return
 * \return the solution
 * throws SingularException if the matrix is singular
 */
VectorN& LinAlg::solve_tridiagonal_fast(VectorN& dl, VectorN& d, VectorN& du, VectorN& xb)
{
  // determine N
  INTEGER N = d.size();
  if (xb.rows() != N)
    throw MissizeException();

  // call the tridiagonal solver
  INTEGER LDB = xb.rows();
  INTEGER NRHS = 1;
  INTEGER INFO;
  gtsv_(&N, &NRHS, dl.data(), d.data(), du.data(), xb.data(), &LDB, &INFO);

  // see whether solution successful
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return xb;
}

/// Solves a tridiagonal system
/**
 * \param dl the (n-1) elements on the subdiagonal (destroyed on return)
 * \param d  the n elements on the diagonal (destroyed on return)
 * \param du the n elements on the superdiagonal (destroyed on return)
 * \param XB the right hand side on input, the solution on return
 * \return the solution
 * throws SingularException if the matrix is singular
 */
MatrixN& LinAlg::solve_tridiagonal_fast(VectorN& dl, VectorN& d, VectorN& du, MatrixN& XB)
{
  // determine N
  INTEGER N = d.size();
  if (XB.rows() != N)
    throw MissizeException();

  // call the tridiagonal solver
  INTEGER LDB = XB.rows();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;
  gtsv_(&N, &NRHS, dl.data(), d.data(), du.data(), XB.data(), &LDB, &INFO);

  // see whether solution successful
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Performs a LDL' factorization of a symmetric, indefinite matrix
/**
 * \param A the matrix A on input; the factorized matrix on output
 */
void LinAlg::factor_LDL(MatrixN& A, vector<int>& IPIV)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
    return;

  // verify that A is symmetric
  #ifndef NDEBUG
  if (!A.is_symmetric())
    FILE_LOG(LOG_LINALG) << "LinAlg::factor_LDL() warning: matrix is assymetrical!" << endl;
  #endif

  // resize IPIV
  IPIV.resize(A.rows());

  // setup LAPACK args for factorization
  char UPLO = 'L';
  INTEGER N = A.rows();
  INTEGER INFO;

  // get A's data -- we're going to modify it directly
  Real* data = A.begin();

  // alter matrix to put into packed format
  for (unsigned j=0, k=0, r=0; j< A.rows(); j++, r++)
    for (unsigned i=0, s=0; i< A.columns(); i++, s+= A.rows())
      if (i >= j)
         data[k++] = data[r+s];

  // perform the factorization
  sptrf_(&UPLO, &N, A.data(), &IPIV.front(), &INFO);
  assert(INFO >= 0);
}

/// Performs the Cholesky factorization of a matrix
/**
 * \param A the matrix A on input; the factorized (upper triangular) matrix on output
 * \return <b>true</b> if matrix factored successfully, <b>false</b> otherwise
 */
bool LinAlg::factor_chol(MatrixN& A)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
    return true;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER INFO;

  // perform the Cholesky factorization
  potrf_(&UPLO, &N, A.data(), &N, &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    return false;

  // make the matrix upper triangular
  A.zero_lower_triangle();
  return true;
}

/// Performs the QR factorization of a matrix with column pivoting
/**
 * Factorizes A*P = Q*R
 * \param AQ the matrix A on input; the matrix R on output
 * \param Q the matrix Q on output
 * \param PI the column pivots on output
 */
void LinAlg::factor_QR(MatrixN& AR, MatrixN& Q, vector<int>& PI)
{
  // determine LAPACK parameters
  INTEGER M = AR.rows();
  INTEGER N = AR.columns();
  unsigned min_mn = (unsigned) std::min(M,N); 
 
  // setup tau vector
  SAFESTATIC FastThreadable<VectorN> TAU;
  TAU().resize(min_mn);

  // setup PI for entry
  PI.resize(N);
  for (int i=0; i< N; i++)
    PI[i] = 0;

  // call LAPACK
  INTEGER LDA = AR.rows();
  INTEGER INFO;
  geqp3_(&M, &N, AR.begin(), &LDA, &PI[0], TAU().begin(), &INFO);
  assert(INFO == 0);

  // correct indices for PI
  for (int i=0; i< N; i++)
    PI[i]--;

  // resize AR
  AR.resize(std::min(AR.rows(), AR.columns()), AR.columns(), true);

  // make R upper triangular
  AR.zero_lower_triangle();
}

/// Performs the QR factorization of a matrix
/**
 * \param AQ the m x n matrix A on input; the matrix min(m,n) x n R on output 
 * \param Q the m x min(m,n) matrix Q on output
 */
void LinAlg::factor_QR(MatrixN& AR, MatrixN& Q)
{
  SAFESTATIC FastThreadable<MatrixN> WORKM;
  SAFESTATIC FastThreadable<VectorN> TAU;

  // get matrix/vector
  VectorN& tau = TAU();
  MatrixN& workM = WORKM();

  // setup constants
  const unsigned m = AR.rows();
  const unsigned n = AR.columns();

  // check for zero sized matrix
  if (m == 0 || n == 0)
    return;

  // determine LAPACK parameters
  INTEGER M = AR.rows();
  INTEGER N = AR.columns();
  INTEGER MINMN = std::min(M, N);
  const unsigned min_mn = (unsigned) MINMN; 
 
  // setup tau vector
  tau.resize(min_mn);

  // call LAPACK
  INTEGER LDA = M;
  INTEGER INFO;
  geqrf_(&M, &N, AR.begin(), &M, tau.begin(), &INFO);
  assert(INFO == 0);

  // setup Q
  Q.resize(m,m);
  std::copy(AR.begin(), AR.begin()+m*min_mn, Q.begin());
  orgqr_(&M, &MINMN, &MINMN, Q.begin(), &M, tau.begin(), &INFO);
  
  // make R triangular
  AR.zero_lower_triangle();

  // note: R is m x n, so we don't have to resize 
}

/// Performs the LU factorization of a matrix
/**
 * \param A the m x n matrix to be factored; on exit, the factors L and U
 *        from the factorization M = P*L*U (unit diagonal elements of L are
 *        not stored)
 * \param IPIV on output, contains the pivot indices (for 1 <= i <= min(m,n),
 *        row i of A was interchanged with row IPIV[i])
 * \return <b>false</b> if A is singular, <b>true</b> otherwise
 */
bool LinAlg::factor_LU(MatrixN& A, vector<int>& IPIV)
{
  if (A.rows() == 0 || A.columns() == 0)
    return true;

  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  
  // resize IPIV if necessary
  if (IPIV.size() < std::min(A.rows(), A.columns()))
    IPIV.resize(std::min(A.rows(),A.columns()));

  // call LAPACK
  INTEGER INFO;
  getrf_(&M, &N, A.begin(), &M, &IPIV.front(), &INFO);
  return INFO == 0;
}

/// Computes a matrix inverse using the factorization determined via factor_LU()
/**
 * \param M the LU-factored matrix
 * \param IPIV the pivots determined by factor_LU()
 * \note throws SingularException if matrix is singular
 */
MatrixN& LinAlg::inverse_LU(MatrixN& M, const vector<int>& IPIV)
{
  if (M.rows() == 0 || M.columns() == 0)
    return M;

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // call lapack
  INTEGER INFO;
  INTEGER N = M.rows();
  getri_(&N, M.begin(), &N, (INTEGER*) &IPIV.front(), &INFO);
  if (INFO > 0)
    throw SingularException();

  return M;
}

/// Solves a triangular system of linear equations
/**
 * \param A the matrix
 * \param utri if <b>true</b> A is upper triangular (lower triangular otherwise)
 * \param transpose_A if <b>true</b>, solves A'*x = b
 * \param xb contains b on entry, x on return
 * \return reference to xb
 */
VectorN& LinAlg::solve_tri_fast(const MatrixN& A, bool utri, bool transpose_A, VectorN& xb)
{
  if (A.rows() != xb.size())
    throw MissizeException();

  if (A.rows() == A.columns())
    throw NonsquareMatrixException();

  if (A.rows() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char TRANS = (transpose_A) ? 'T' : 'N';
  char UPLO = (utri) ? 'U' : 'L';
  INTEGER N = A.rows();
  INTEGER LDA = A.rows();
  INTEGER LDB = xb.size();
  INTEGER NRHS = 1;
  INTEGER INFO;
  trtrs_(&UPLO, &TRANS, &N, &NRHS, (Real*) A.begin(), &LDA, xb.begin(), &LDB, &INFO);

  return xb;
}

/// Solves a triangular system of linear equations
/**
 * \param A the matrix
 * \param utri if <b>true</b> A is upper triangular (lower triangular otherwise)
 * \param transpose_A if <b>true</b>, solves A'*X = B
 * \param XB contains B on entry, X on return
 * \return reference to xb
 */
MatrixN& LinAlg::solve_tri_fast(const MatrixN& A, bool utri, bool transpose_A, MatrixN& XB)
{
  if (A.rows() != XB.rows())
    throw MissizeException();

  if (A.rows() == A.columns())
    throw NonsquareMatrixException();

  if (A.rows() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char TRANS = (transpose_A) ? 'T' : 'N';
  char UPLO = (utri) ? 'U' : 'L';
  INTEGER N = A.rows();
  INTEGER LDA = A.rows();
  INTEGER LDB = XB.rows();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;
  trtrs_(&UPLO, &TRANS, &N, &NRHS, (Real*) A.begin(), &LDA, XB.begin(), &LDB, &INFO);

  return XB;
}

/// Solves a system of linear equations using the factorization determined via factor_LDL()
/**
 * \param M the factorization performed by factor_LDL()
 * \param xb the right hand side (b) on input, the vector x on return
 */
VectorN& LinAlg::solve_LDL_fast(const MatrixN& M, const vector<int>& IPIV, VectorN& xb)
{
  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  if (M.rows() != xb.size())
    throw MissizeException();

  if (M.rows() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'L';
  INTEGER N = M.rows();
  INTEGER NRHS = 1;
  INTEGER LDB = xb.size();
  INTEGER INFO;

  // call the solver routine
  sptrs_(&UPLO, &N, &NRHS, (Real*) M.begin(), (int*) &IPIV.front(), xb.begin(), &LDB, &INFO);
  assert(INFO == 0);

  return xb;
}

/// Solves systems of linear equations using the factorization determined via factor_LDL()
/**
 * \param M the factorization performed by factor_LDL()
 * \param XB the right hand sides on input, the vectors X on return
 */
MatrixN& LinAlg::solve_LDL_fast(const MatrixN& M, const vector<int>& IPIV, MatrixN& XB)
{
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // check for empty matrix
  if (M.rows() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.rows();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = XB.rows();
  INTEGER INFO;

  // call the solver routine
  sptrs_(&UPLO, &N, &NRHS, (Real*) M.begin(), (int*) &IPIV.front(), XB.begin(), &LDB, &INFO);
  assert(INFO == 0);

  return XB;
}

/// Solves a system of linear equations using the factorization determined via factor_chol()
/**
 * \param M the Cholesky decomposition performed by factor_chol()
 * \param b the right hand side
 * \return the solution x to Mx = b
 */
VectorN LinAlg::solve_chol(const MatrixN& M, const VectorN& b)
{
  VectorN x;
  x = b;
  solve_chol_fast(M, x);
  return x;
}

/// Solves a system of linear equations using the factorization determined via factor_chol()
/**
 * \param M the Cholesky decomposition performed by factor_chol()
 * \param xb the right hand side on input, the vector x on return
 */
VectorN& LinAlg::solve_chol_fast(const MatrixN& M, VectorN& xb)
{
  if (M.rows() != xb.size())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  if (M.rows() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.rows();
  INTEGER NRHS = 1;
  INTEGER LDA = M.rows();
  INTEGER LDB = xb.size();
  INTEGER INFO;

  // call the solver routine
  potrs_(&UPLO, &N, &NRHS, (Real*) M.begin(), &LDA, xb.begin(), &LDB, &INFO);
  assert(INFO == 0);

  return xb;
}

/// Solves a system of linear equations using the factorization determined via factor_chol()
/**
 * \param M the Cholesky decomposition performed by factor_chol()
 * \param B the right hand side
 * \return the solution X to MX = B
 */
MatrixN LinAlg::solve_chol(const MatrixN& M, const MatrixN& B)
{
  MatrixN X;
  X = B;
  solve_chol_fast(M, X);
  return X;
}

/// Solves a system of linear equations using the factorization determined via factor_chol()
/**
 * \param M the Cholesky decomposition performed by factor_chol()
 * \param XB the right hand sides on input, the vectors X on return
 */
MatrixN& LinAlg::solve_chol_fast(const MatrixN& M, MatrixN& XB)
{
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // check for empty matrices
  if (M.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.rows();
  INTEGER NRHS = XB.columns();
  INTEGER LDA = M.rows();
  INTEGER LDB = XB.rows();
  INTEGER INFO;

  // call the solver routine
  potrs_(&UPLO, &N, &NRHS, (Real*) M.begin(), &LDA, XB.begin(), &LDB, &INFO);
  assert(INFO == 0);

  return XB;
}

/// Solves a system of linear equations using the factorization determined via factor_LU()
/**
 * \param M the LU factorization performed by factor_LU()
 * \param transpose if <b>true</b>, solves M'x = b
 * \param IPIV pivots computed by factor_LU()
 * \param b the right hand side
 * \return inv(M)*b or inv(M)'*b
 */
VectorN LinAlg::solve_LU(const MatrixN& M, bool transpose, const vector<int>& IPIV, const VectorN& b)
{
  VectorN x;
  x = b;
  return solve_LU_fast(M, transpose, IPIV, x);
}

/// Solves a system of linear equations using the factorization determined via factor_LU()
/**
 * \param M the LU factorization performed by factor_LU()
 * \param xb the right hand side on input, the vector x on return
 * \param transpose if <b>true</b>, solves M'x = b
 * \param IPIV pivots computed by factor_LU()
 */
VectorN& LinAlg::solve_LU_fast(const MatrixN& M, bool transpose, const vector<int>& IPIV, VectorN& xb)
{
  if (M.rows() == 0)
    return xb.set_zero();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  if (M.rows() != xb.size())
    throw MissizeException();

  // setup parameters to LAPACK
  INTEGER N = M.rows();
  INTEGER NRHS = 1;
  char TRANS = (transpose) ? 'T' : 'N';
  
  // call LAPACK
  INTEGER INFO;
  getrs_(&TRANS, &N, &NRHS, (Real*) M.begin(), &N, (INTEGER*) &IPIV.front(), xb.begin(), &N, &INFO);

  return xb;
}

/// Solves a system of linear equations using the factorization determined via factor_LU()
/**
 * \param M the LU factorization performed by factor_LU()
 * \param B the right hand side
 * \param transpose if <b>true</b>, solves M'X = B
 * \param IPIV pivots computed by factor_LU()
 * \return inv(M)*B or inv(M)'*B
 */
MatrixN LinAlg::solve_LU(const MatrixN& M, bool transpose, const vector<int>& IPIV, const MatrixN& B)
{
  MatrixN X;
  X = B;
  return solve_LU_fast(M, transpose, IPIV, X);
}

/// Solves a system of linear equations using the factorization determined via factor_LU()
/**
 * \param M the LU factorization performed by factor_LU()
 * \param XB the right hand side on input, the matrix X on return
 * \param transpose if <b>true</b>, solves M'X = B
 * \param IPIV pivots computed by factor_LU()
 */
MatrixN& LinAlg::solve_LU_fast(const MatrixN& M, bool transpose, const vector<int>& IPIV, MatrixN& XB)
{
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();

  // check for empty matrix
  if (M.rows() == 0)
    return XB.set_zero();

  // setup parameters to LAPACK
  INTEGER N = M.rows();
  INTEGER NRHS = XB.columns();
  char TRANS = (transpose) ? 'T' : 'N';
  
  // call LAPACK
  INTEGER INFO;
  getrs_(&TRANS, &N, &NRHS, (Real*) M.begin(), &N, (INTEGER*) &IPIV.front(), XB.begin(), &N, &INFO);

  return XB;
}

/// Calculates the rank of a matrix
unsigned LinAlg::calc_rank(MatrixN& A, Real tol)
{
  // look for easy out
  if (A.rows() == 0 || A.columns() == 0)
    return 0;

  // compute the SVD of A
  SAFESTATIC FastThreadable<MatrixN> U, V;
  SAFESTATIC FastThreadable<VectorN> S;
  MatrixN& Ux = U();
  MatrixN& Vx = V();
  VectorN& Sx = S();
  svd(A, Ux, Sx, Vx);

  // get the dimensions of A
  unsigned m = A.rows();
  unsigned n = A.columns();
  unsigned maxmn = std::max(m,n);
  unsigned minmn = std::min(m,n);

  // get the # of singular values
  Real tolerance = (tol < 0.0) ? S()[0] * maxmn * std::numeric_limits<Real>::epsilon() : tol;
  unsigned ns = 0;
  for (unsigned i=Sx.size()-1; i > 0; i--, ns++)
    if (Sx[i] > tolerance)
      break;

  assert(minmn > ns);
  return minmn - ns;
}

/// Computes the nullspace of a matrix
/**
 * \note A is destroyed on return
 */
MatrixN& LinAlg::nullspace(MatrixN& A, MatrixN& nullspace, Real tol)
{
  SAFESTATIC FastThreadable<MatrixN> U, V;
  SAFESTATIC FastThreadable<VectorN> S;

  // look for fast way out
  if (A.rows() == 0)
  {
    nullspace.set_zero(A.columns(), A.columns());
    for (unsigned i=0; i< A.columns(); i++)
      nullspace(i,i) = (Real) 1.0;
    return nullspace;
  }

  // compute the SVD of A
  MatrixN& Ux = U();
  VectorN& Sx = S();
  svd(A, Ux, Sx, nullspace);

  // get the dimensions of A
  unsigned m = A.rows();
  unsigned n = A.columns();
  boost::tuple<unsigned, unsigned> min_max = boost::minmax(m, n);
  unsigned minmn = min_max.get<0>();
  unsigned maxmn = min_max.get<1>();

  // get the # of singular values
  unsigned ns = 0;
  Real tolerance = (tol < 0.0) ? Sx[0] * maxmn * std::numeric_limits<Real>::epsilon() : tol;
  for (unsigned i=Sx.size()-1; i > 0; i--, ns++)
    if (Sx[i] > tolerance)
      break;

  // add in # of singular values from non-square matrices
  for (unsigned i=minmn; i < n; i++)
    ns++;

  // shift nullspace (if necesary)
  if (ns == nullspace.columns())
    return nullspace;
  else
  {
    BlockIterator bi = nullspace.block_start(0, nullspace.rows(), nullspace.columns()-ns, nullspace.columns());
    std::copy(bi, bi+ns*nullspace.rows(), nullspace.begin());
    return nullspace;
  }
}

/// Computes the condition number of a matrix
Real LinAlg::cond(MatrixN& A)
{
  SAFESTATIC FastThreadable<MatrixN> U, V;
  SAFESTATIC FastThreadable<VectorN> S;
  svd(A, U(), S(), V());
  return S()[0] / S()[S().size()-1];
}

/// Computes the psuedo-inverse of a matrix
MatrixN& LinAlg::pseudo_inverse(MatrixN& A, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol)
{
  // get the dimensionality of A
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy out
  if (m == 0)
  {
    A.resize(n,0);
    return A;
  }
  else if (n == 0)
  {
    A.resize(0, m);
    return A;
  }

  // init matrices for svd computation 
  SAFESTATIC FastThreadable<MatrixN> U, V, workM;
  SAFESTATIC FastThreadable<VectorN> S;

  // compute the svd
  MatrixN& Ux = U();
  MatrixN& Vx = V();
  VectorN& Sx = S();
  MatrixN& workMx = workM();
  svd(A, Ux, Sx, Vx);
  
  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
    tol = Sx[0] * std::max(m,n) * std::numeric_limits<Real>::epsilon();

  // compute 1/S
  unsigned S_len = Sx.size();

  // A is m x n, B is m x k
  // (L -> R, scaling V)    n^2 + n*min(n,m)*m + nmk [n < m < k, n < k < m]
  // (L -> R, scaling U')   m^2 + n*min(n,m)*m + nmk [m < n < k]
  // (R -> L, scaling U')   m^2 + m^2k + nmk + n^2k  [k < n < m]
  // (R -> L, scaling U'*B) m^2k + min(n,m)*k + n*min(m,n)*k [k < m < n, m < k < n]

  // compute inv(s) 
  for (unsigned i=0; i< S_len; i++)
    Sx[i] = (std::fabs(Sx[i]) > tol) ? (Real) 1.0/Sx[i] : (Real) 0.0;

  // scale U' or V, depending on whichever is smaller, m or n
  if (m < n)
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(m, Sx[i], &Ux(0,i), 1);
  else
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx[i], &Vx(0,i), 1);

  // size the result properly
  A.resize(n, m);

  // do the multiplication (V * U')
  CBLAS::gemm(CblasNoTrans, CblasTrans, n, m, minmn, Vx, n, Ux, m, (Real) 1.0, (Real) 0.0, A, n); 
  return A;
}

/// Constructs a vector from a string
/**
 * \param str a string representing the vector; elements are separated using whitespace, comma, or semicolon
 */
VectorN LinAlg::to_vector(const std::string& str)
{
  std::vector<Real> TMP;
  size_t idx = 0;

  // get the first 
  while ((idx = str.find_first_not_of(" \t\r\n;,", idx)) != std::string::npos)
  {
    size_t firstsep = str.find_first_of(" \t\r\n;,", idx);
    unsigned n = (firstsep != std::string::npos) ? firstsep - idx : str.size() - idx;
    TMP.push_back((Real) std::atof(str.substr(idx, n).c_str()));
    idx = firstsep;
  }

  // form the vector
  VectorN v(TMP.size());
  for (unsigned i=0; i< TMP.size(); i++)
    v[i] = TMP[i];

  return v;  
}

/// Constructs a matrix from a string
/** 
 * Elements in a row are separated using whitespace or commas.  Rows are
 * separated using semicolons.
 */
MatrixN LinAlg::to_matrix(const std::string& str)
{
  std::vector<std::vector<Real> > TMP;
  std::vector<std::string> strings;
  size_t idx = 0;

  // process each row independently
  int last = 0;
  for (unsigned int i=0; i< str.size(); i++)
  {
    if (str[i] == ';')
    {
      if (i - last <= 0)
        continue;
      strings.push_back(str.substr(last, i-last));
      last = i+1;
    }
    else if (i == str.size()-1 && i - last > 0)
      strings.push_back(str.substr(last, i-last+1));
  }
    
  // process..
  for (unsigned int i=0; i< strings.size(); i++)
  {
    // reset the index
    idx = 0;

    // add a new vector
    TMP.push_back(std::vector<Real>());
    
    // get the first 
    while ((idx = strings[i].find_first_not_of(" \t\r\n,", idx)) != std::string::npos)
    {
      size_t firstsep = strings[i].find_first_of(" \t\r\n,", idx);
      unsigned n = (firstsep != std::string::npos) ? firstsep - idx : strings[i].size() - idx;
      TMP[i].push_back((Real) std::atof(strings[i].substr(idx, n).c_str()));
      idx = firstsep;
    }

    // check for row size mismatch
    if (i > 0 && TMP[i].size() != TMP[i-1].size())
    {
      std::cerr << "LinAlg::to_matrix() - row size mismatch!" << std::endl;
      for (unsigned j=0; j<= i; j++)
        std::cerr << "Elements in row " << j << ": " << TMP[j].size() << std::endl;
      return MatrixN();
    }
  }
  
  // form the matrix
  MatrixN m(TMP.size(), TMP.front().size());
  for (unsigned int i=0; i< TMP.size(); i++)
    for (unsigned int j=0; j< TMP.front().size(); j++)
      m(i,j) = TMP[i][j];

  return m;
}

/// Determines whether a symmetric matrix is positive semi-definite
bool LinAlg::is_SPSD(const MatrixN& m, Real tol)
{
  SAFESTATIC FastThreadable<VectorN> evals_f;
  SAFESTATIC FastThreadable<MatrixN> m_copy;
  VectorN& evals = evals_f();
  eig_symm(m_copy(), evals);

  // make tolerance positive, if it is not already
  if (tol < (Real) 0.0)
    tol = std::fabs(evals[evals.size()-1]) * m.rows() * std::numeric_limits<Real>::max();

  // check whether all eigenvalues are non-negative to numerical tolerance
  for (unsigned i=0; i< evals.size(); i++)
    if (evals[i] < -tol)
      return false;

  return true;
}

/// Determines whether a matrix is positive-definite
bool LinAlg::is_SPD(const MatrixN& m, Real tol)
{
  // get the eigenvalues of the matrix
  SAFESTATIC FastThreadable<VectorN> evals_f;
  SAFESTATIC FastThreadable<MatrixN> m_copy;
  VectorN& evals = evals_f();
  eig_symm(m_copy(), evals);

  // make tolerance positive, if it is not already
  if (tol < (Real) 0.0)
    tol = std::fabs(evals[evals.size()-1]) * m.rows() * std::numeric_limits<Real>::max();

  // check whether all eigenvalues are positive to numerical tolerance
  for (unsigned i=0; i< evals.size(); i++)
    if (evals[i] < tol)
      return false;

  return true;
}

/// Computes the eigenvalues of the matrix A
/**
 * \param A a matrix
 * \param evals on return, the eigenvalues will be stored here in ascending order
 */
void LinAlg::eig_symm(MatrixN& A, VectorN& evals)
{
  // make sure A is not zero sized
  if (A.rows() == 0)
  {
    evals.resize(0);
    return;
  }

  // verify that A is square
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that the matrix is symmetric
  #ifndef NDEBUG
  if (!A.is_symmetric())
    std::cerr << "LinAlg::eig_symm() - matrix does not appear to be symmetric!" << std::endl;
  #endif

  // make sure that the eigenvalues array is the proper size
  evals.resize(A.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'N';
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.rows();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, A.data(), &LDA, evals.data(), &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

/// Computes the eigenvalues and eigenvectors of the matrix A
/**
 * \param A a square symmetric matrix on input, eigenvectors corresponding to eigenvalues on return
 * \param evals on return, the eigenvalues will be stored here in ascending order
 */
void LinAlg::eig_symm_plus(MatrixN& A_evecs, VectorN& evals)
{
  // make sure that A is not zero sized
  if (A_evecs.rows() == 0 || A_evecs.columns() == 0)
  {
    evals.resize(0);
    A_evecs.resize(0,0);
    return;
  }

  if (A_evecs.rows() != A_evecs.columns())
    throw NonsquareMatrixException();

  // verify that the matrix is symmetric
  #ifndef NDEBUG
  if (!A_evecs.is_symmetric())
    std::cerr << "LinAlg::eig_symm() - matrix does not appear to be symmetric!" << std::endl;
  #endif

  // make sure that the eigenvalues array is the proper size
  evals.resize(A_evecs.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'V';
  char UPLO = 'U';
  INTEGER N = A_evecs.rows();
  INTEGER LDA = A_evecs.rows();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, A_evecs.begin(), &LDA, evals.begin(), &INFO);
  assert(INFO == 0);
  
  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

void LinAlg::svd(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V)
{
  SAFESTATIC FastThreadable<MatrixN> A_backup;

  // copy A
  A_backup().copy_from(A);

  try
  {
    svd1(A, U, S, V);
  }
  catch (NumericalException e)
  {
    svd2(A_backup(), U, S, V); 
  }
}

/// Does an 'in place' SVD (destroying A)
/**
 * The singular value decomposition of A is U*S*V' (' is the transpose 
 * operator); to recompose A, it will be necessary to transpose V before
 * multiplication (i.e., V is returned by the algorithm, not V').
 * Note: passed matrices and vectors U, S, and V are resized as necessary. 
 * \param A the matrix on which the SVD will be performed (destroyed on return)
 * \param U on output, a A.rows() x A.rows() orthogonal matrix
 * \param S on output, a min(A.rows(), A.columns()) length vector of singular values
 * \param V on output, a A.columns() x A.columns() orthogonal matrix
 */
void LinAlg::svd1(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V)
{
  // make sure that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
  {
    U.set_zero(A.rows(), A.rows());
    S.resize(0);
    V.set_zero(A.columns(), A.columns());
    return;
  } 

  // setup U
  if (U.rows() != A.rows() || U.columns() != A.rows())
    U.resize(A.rows(), A.rows());

  // setup S
  unsigned minmn = std::min(A.rows(), A.columns());
  if (S.size() != minmn)
    S.resize(minmn);

  // setup V
  if (V.rows() != A.columns() || V.columns() != A.columns())
    V.resize(A.columns(), A.columns());

  // setup call to LAPACK
  char JOB = 'A';
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER LDA = A.rows();
  INTEGER LDU = U.rows();
  INTEGER LDVT = V.rows();
  INTEGER INFO;

  // call LAPACK 
  gesvd_(&JOB, &JOB, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Singular value decomposition failed to converge");

  // transpose V
  V.transpose();
}

/// Does an 'in place' SVD (destroying A), using divide and conquer algorithm
/**
 * The singular value decomposition of A is U*S*V' (' is the transpose 
 * operator); to recompose A, it will be necessary to transpose V before
 * multiplication (i.e., V is returned by the algorithm, not V').
 * Note: passed matrices and vectors U, S, and V are resized as necessary. 
 * \param A the matrix on which the SVD will be performed (destroyed on return)
 * \param U on output, a A.rows() x A.rows() orthogonal matrix
 * \param S on output, a min(A.rows(), A.columns()) length vector of singular values
 * \param V on output, a A.columns() x A.columns() orthogonal matrix
 */
void LinAlg::svd2(MatrixN& A, MatrixN& U, VectorN& S, MatrixN& V)
{
  // make sure that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
  {
    U.set_zero(A.rows(), A.rows());
    S.resize(0);
    V.set_zero(A.columns(), A.columns());
    return;
  } 

  // setup U
  if (U.rows() != A.rows() || U.columns() != A.rows())
    U.resize(A.rows(), A.rows());

  // setup S
  unsigned minmn = std::min(A.rows(), A.columns());
  if (S.size() != minmn)
    S.resize(minmn);

  // setup V
  if (V.rows() != A.columns() || V.columns() != A.columns())
    V.resize(A.columns(), A.columns());

  // setup call to LAPACK
  char JOBZ = 'A';
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER LDA = A.rows();
  INTEGER LDU = U.rows();
  INTEGER LDVT = V.rows();
  INTEGER INFO;

  // call LAPACK 
  gesdd_(&JOBZ, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Singular value decomposition failed to converge");

  // transpose V
  V.transpose();
}
/// Solves a symmetric, indefinite square matrix
/**
 * \param A the matrix to be solved; the matrix is destroyed on return
 * \param xb the RHS b (A*x = b) on input; the solution, x, on return
 */
VectorN& LinAlg::solve_symmetric_fast(MatrixN& A, VectorN& xb)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify A and b are compatible 
  if (A.columns() != xb.size())
    throw MissizeException();

  // make sure that A is not zero sized
  if (A.rows() == 0)
    return xb.resize(0);

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = 1;
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK
  sysv_(&UPLO, &N, &NRHS, A.data(), &N, &IPIV().front(), xb.data(), &N, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return xb;
}

/// Solves a symmetric, indefinite square matrix
/**
 * \param A the matrix to be solved; the matrix is destroyed on return
 * \param XB the RHS B (A*X = B) on input; the solution, X, on return
 */
MatrixN& LinAlg::solve_symmetric_fast(MatrixN& A, MatrixN& XB)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify A and b are compatible 
  if (A.columns() != XB.rows())
    throw MissizeException();

  // make sure that A is not zero sized
  if (A.rows() == 0)
    return XB.set_zero(0, 0);

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = XB.columns();
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK
  sysv_(&UPLO, &N, &NRHS, A.data(), &N, &IPIV().front(), XB.data(), &N, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Solves a symmetric, indefinite square matrix (slow but convenient method)
VectorN LinAlg::solve_symmetric(const MatrixN& A, const VectorN& b)
{
  // copy A and b
  MatrixN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call the main solver
  solve_symmetric_fast(A_copy, x);

  return x;
}

/// Solves multiple symmetric, indefinite systems of linear equations (slow but convenient method)
MatrixN LinAlg::solve_symmetric(const MatrixN& A, const MatrixN& B)
{
  // copy A and B 
  MatrixN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call the main solver
  solve_symmetric_fast(A_copy, X);

  return X;
}

/// Inverts the symmetric, indefinite matrix A
/**
 * \param A a square, symmetric indefinite matrix; inverse will be contained
 *        here on return
 */
MatrixN& LinAlg::inverse_symmetric(MatrixN& A)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero size
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER INFO;

  // setup pivot array
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);

  // perform the necessary factorization
  sytrf_(&UPLO, &N, A.data(), &N, &IPIV().front(), &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // perform the inversion
  sytri_(&UPLO, &N, A.data(), &N, &IPIV().front(), &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // now, make the matrix symmetric
  Real* data = A.data();
  for (unsigned i=1, ii=A.rows(); i< N; i++, ii += A.rows())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.rows())
      data[jj+i] = data[ii+j];

  return A;
}

/// Solves a system of equations A*x = b using a symmetric, positive-definite square matrix
/**
 * \param A the matrix coefficients; this matrix will be destroyed on return
 * \param xb on input b, on output x
 */
VectorN& LinAlg::solve_SPD_fast(MatrixN& A, VectorN& xb)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A and b are proper size
  if (A.columns() != xb.size())
    throw MissizeException();

  // verify that A is not zero size
  if (A.rows() == 0)
    return xb.resize(0);

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = 1;
  INTEGER INFO;

  // call LAPACK
  posv_(&UPLO, &N, &NRHS, A.data(), &N, xb.data(), &N, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return xb;
}

/// Solves a system of equations A*X = B using a symmetric, positive-definite square matrix
/**
 * \param A the matrix coefficients; this matrix will be destroyed on return
 * \param XB on input B, on output X
 */
MatrixN& LinAlg::solve_SPD_fast(MatrixN& A, MatrixN& XB)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A and b are proper size
  if (A.columns() != XB.rows())
    throw MissizeException();

  // verify that A is not zero size
  if (A.rows() == 0)
    return XB.set_zero(A.rows(), XB.columns());

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;

  // call LAPACK
  posv_(&UPLO, &N, &NRHS, A.data(), &N, XB.data(), &N, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Solves a system of equations A*X = B using a symmetric, positive-definite square matrix
MatrixN LinAlg::solve_SPD(const MatrixN& A, const MatrixN& B)
{
  // copy A and B
  MatrixN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call the main solver
  solve_SPD_fast(A_copy, X);

  return X;
}

/// Solves a system of equations A*x = b using a symmetric, positive-definite square matrix
VectorN LinAlg::solve_SPD(const MatrixN& A, const VectorN& b)
{
  // copy A and b
  MatrixN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call the main solver
  solve_SPD_fast(A_copy, x);

  return x;
}

/// Inverts the symmetric, positive-definite matrix A using Cholesky factorization
/**
 * \param A the Cholesky factorization of a matrix; contains the inverse
 *        on return
 */
MatrixN& LinAlg::inverse_chol(MatrixN& A)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero sized
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER INFO;

  // perform the inverse
  potri_(&UPLO, &N, A.data(), &N, &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // now, make the matrix symmetric
  Real* data = A.data();
  for (unsigned i=1, ii=A.rows(); i< N; i++, ii += A.rows())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.rows())
      data[jj+i] = data[ii+j];


  return A;
}

/// Inverts the symmetric, positive-definite matrix A using Cholesky factorization
/**
 * \param A a square, symmetric positive-definite matrix; contains the inverse
 *        on return
 */
MatrixN& LinAlg::inverse_SPD(MatrixN& A)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero sized
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER INFO;

  // perform the Cholesky factorization
  potrf_(&UPLO, &N, A.data(), &N, &INFO);
  assert(INFO >= 0);

  // perform the inverse
  potri_(&UPLO, &N, A.data(), &N, &INFO);
  if (INFO > 0)
    throw SingularException();

  // now, make the matrix symmetric
  Real* data = A.data();
  for (unsigned i=1, ii=A.rows(); i< N; i++, ii += A.rows())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.rows())
      data[jj+i] = data[ii+j];

  return A;
}

/// Inverts the matrix A using LU factorization
/**
 * \param A a square matrix; contains the inverse on return
 */
MatrixN& LinAlg::inverse(MatrixN& A)
{
  SAFESTATIC FastThreadable<vector<int> > pivots;
  factor_LU(A, pivots());
  inverse_LU(A, pivots());

  return A;
}

/// Less robust least squares solver (solves Ax = b)
/**
 * \note this method does not work!
 */
/*
VectorN& LinAlg::solve_LS_fast2(MatrixN& A, VectorN& xb)
{
  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
    return xb;

  // verify that A and b are appropriate sizes
  if (A.rows() != xb.size())
    throw MissizeException();

  // do QR factorization 
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER min_mn = std::min(A.rows(), A.columns());
 
  // setup tau vector
  SAFESTATIC FastThreadable<VectorN> TAU;
  TAU().resize(min_mn);

  // call LAPACK
  INTEGER INFO;
  geqrf_(&M, &N, A.begin(), &M, TAU().begin(), &INFO);
  assert(INFO == 0);

  // call ORMQR for work query
  char SIDE = 'L';
  char TRANS = 'T';
  INTEGER CM = xb.size();
  INTEGER CN = 1;
  INTEGER K = TAU().size();
  INTEGER LDA = A.rows();
  INTEGER LDC = CM;
  ormqr_(&SIDE, &TRANS, &CM, &CN, &K, A.begin(), &LDA, TAU().begin(), xb.begin(), &LDC, &INFO);

  // solve triangular system
  char UPLO = 'U';
  INTEGER LDB = xb.size();
  INTEGER NRHS = 1;
  trtrs_(&UPLO, &TRANS, &N, &NRHS, A.begin(), &LDA, xb.begin(), &LDB, &INFO);

  return xb; 
} 
*/

/// Sparse least squares solver (solves Ax = b)
/**
 * Solves rank-deficient and undertermined (minimum norm solution) systems.
 * Computes least-squares solutio nto overdetermined systems.
 * \param A the coefficient matrix
 * \param b the r.h.s.
 * \param x on output, the least-squares solution to the system
 * \return a reference to x
 */
VectorN& LinAlg::solve_LS(const SparseMatrixN& A, const VectorN& b, VectorN& x, Real damp, unsigned max_iter, Real tol)
{
  // init variables for LSQR
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  Real dampsq = damp*damp;
  Real ddnorm = (Real) 0.0;
  Real res2 = (Real) 0.0;
  Real anorm = (Real) 0.0;
  Real xnorm = (Real) 0.0;
  Real xxnorm = (Real) 0.0;
  Real z = (Real) 0.0;
  Real cs2 = (Real) -1.0;
  Real sn2 = (Real) 0.0;

  // setup necessary vectors
  VectorN u, v, w, tmp, dk;

  // setup first vectors u and v for bidiagonalization.  These satisfy 
  // beta*u = b, alpha*v = A'u
  b.get_sub_vec(0, m, u);
  x.set_zero(n);
  Real alpha = (Real) 0.0;
  Real beta = u.norm();
  if (beta > (Real) 0.0)
  {
    u /= beta;
    A.transpose_mult(u, v);
    alpha = v.norm();
  }
  if (alpha > (Real) 0.0)
  {
    v/= alpha;
    w.copy_from(v);
  }

  // look for solutio x = 0
  Real arnorm = alpha * beta;
  if (arnorm == (Real) 0.0)
    return x;

  Real rhobar = alpha;
  Real phibar = beta;
  Real bnorm = beta;
  Real rnorm = beta;
  Real r1norm = rnorm;
  Real r2norm = rnorm;

  // main iteration loop
  for (unsigned iter = 0; iter < max_iter; iter++)
  {
    // perform the next step of the bidiagonalization to obtain the next beta,
    // u, alpha, v.  These satisfy the relations
    // beta*u = a*v - alpha*u
    // alpha*v = A'*u - beta*v
    tmp.copy_from(u);
    u *= alpha;
    A.mult(v, u) -= tmp;
    beta = u.norm();
    if (beta > (Real) 0.0)
    {
      u /= beta;
      anorm = std::sqrt(anorm*anorm + alpha*alpha + beta*beta + damp*damp);
      tmp.copy_from(v);
      tmp *= beta;
      A.transpose_mult(u, v) -= tmp;
      alpha = v.norm();
      if (alpha > (Real) 0.0)
        v /= alpha;
    }

    // use a plane rotation to eliminate the damping parameter.  This alters
    // the diagonal (rhobar) of the lower bidiagonal matrix
    Real rhobar1 = std::sqrt(rhobar*rhobar + damp*damp);
    Real cs1 = rhobar / rhobar1;
    Real sn1 = damp / rhobar1;
    Real psi = sn1 * phibar;
    phibar = cs1 * phibar;

    // use a plane rotation to eliminate the subdiagonal element (beta) of the
    // lower bidiagonal matrix, giving an upper bidiagonal matrix
    Real rho = std::sqrt(rhobar1*rhobar1 + beta*beta);
    Real cs = rhobar1/rho;
    Real sn = beta/rho;
    Real theta = sn*alpha;
    rhobar = -cs*alpha;
    Real phi = cs*phibar;
    phibar = sn*phibar;
    Real tau = sn*phi;

    // update x and w
    Real t1 = phi / rho;
    Real t2 = -theta / rho;
    dk.copy_from(w);
    dk /= rho;
    tmp.copy_from(w) *= t1;
    x += tmp;
    tmp.copy_from(w) *= t2;
    w.copy_from(v) += tmp;
    ddnorm += dk.norm_sq();

    // use a plane rotation on the right to eliminate the super-diagonal
    // element (theta) of the upper bidiagonal matrix.  Then use the result
    // to estimate norm(x)
    Real delta = sn2 * rho;
    Real gambar = -cs2 * rho;
    Real rhs = phi - delta * z;
    Real zbar = rhs / gambar;
    xnorm = std::sqrt(xxnorm + zbar*zbar);
    Real gamma = std::sqrt(gambar*gambar + theta*theta);
    cs2 = gambar / gamma;
    sn2 = theta / gamma;
    Real z = rhs / gamma;
    xxnorm = xxnorm + z*z;

    // test for convergence.  first, estimate the norms of rbar and Abar'rbar
    Real res1 = phibar*phibar;
    res2 = res2 + psi*psi;
    rnorm = std::sqrt(res1+res2);
    arnorm = alpha*std::fabs(tau);
    Real r1sq = rnorm*rnorm - dampsq*xxnorm;
    r1norm = std::sqrt(std::fabs(r1sq));
    if (r1sq < 0)
      r1norm = -r1norm;
    r2norm = rnorm;

    // now use these norms to estimate certain other quantities, some of which
    // will be small near a solution
    Real test1 = rnorm / bnorm;
    Real test2 = arnorm/(anorm*rnorm);
    t1 = test1 / ((Real) 1.0 + anorm*xnorm/bnorm);
    Real rtol = tol + tol*anorm*xnorm / bnorm;

    // following tests guard against extremely small values of tol
    if ((Real) 1.0 + test2 <= (Real) 1.0)
      break;
    if ((Real) 1.0 + t1 <= (Real) 1.0)
      break;

    // allow for tolerances set by the user
    if (test2 <= tol || test1 <= rtol)
      break;
  }

  return x;
}
 
/// Most robust system of linear equations solver (solves Ax = b)
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param A the coefficient matrix (destroyed on return)
 * \param xb the vector b on input, the vector x on return
 * \param tol the tolerance for determining the rank of A; if tol < 0.0,
 *        tol is computed using machine epsilon
 */
VectorN& LinAlg::solve_LS_fast(MatrixN& A, VectorN& xb, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol)
{
  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
  {
    xb.resize(0);
    return xb;
  }

  // verify that A and b are appropriate sizes
  if (A.rows() != xb.size())
    throw MissizeException();

  // get the dimensionality of A
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy out
  if (m == 0 || n == 0)
  {
    xb.set_zero(n);
    return xb;
  }

  // init matrices for svd computation 
  SAFESTATIC FastThreadable<MatrixN> U, V;
  SAFESTATIC FastThreadable<VectorN> S, workv;

  // compute the svd
  MatrixN& Ux = U();
  MatrixN& Vx = V();
  VectorN& Sx = S();
  VectorN& workvx = workv();
  svd(A, Ux, Sx, Vx);
  
  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
    tol = Sx[0] * std::max(m,n) * std::numeric_limits<Real>::epsilon();

  // compute 1/S
  unsigned S_len = Sx.size();

  // compute inv(s) 
  for (unsigned i=0; i< S_len; i++)
    Sx[i] = (std::fabs(Sx[i]) > tol) ? (Real) 1.0/Sx[i] : (Real) 0.0;

  // compute U' * xb (yielding m x 1 vector)
  Ux.transpose_mult(xb, workvx);

  // scale this vector by inv(s) (yielding n x 1 vector, that is really only
  // minmn x 1)
  std::transform(Sx.begin(), Sx.begin()+minmn, workvx.begin(), Sx.begin(), std::multiplies<Real>());

  // multiply V * this vector (yielding n x 1 vector
  xb.resize(n);
  CBLAS::gemv(CblasNoTrans, n, minmn, Vx, n, Sx, 1, (Real) 1.0, (Real) 0.0, xb, 1);
  return xb;

// NOTE: this is disabled b/c it does not work as well...
/*
  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER NRHS = 1;
  INTEGER LDA = M;
  INTEGER LDB = std::max(M,N);
  INTEGER INFO;

  // b will need to contain sufficient elements to store x
  Real* rhs = new Real[std::max(M,N)];
  std::copy(A.begin(), A.end(), rhs);

  // compute
  gelsd_(&M, &N, &NRHS, A.data(), &LDA, rhs, &LDB, &stol, &INFO);

  // don't check success - SVD algorithm is very stable.. 

  // need to resize x
  xb.resize((unsigned) N, true);
  std::copy(rhs, rhs+N, xb.begin());

  // mark A as destroyed and free memory
  A.resize(0,0);
  delete [] rhs;

  return xb;
*/
}

/// Most robust system of linear equations solver (solves AX = B)
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param A the coefficient matrix (destroyed on return)
 * \param XB the matrix B on input, the matrix X on return
 * \param tol the tolerance for determining the rank of A; if tol < 0.0,
 *        tol is computed using machine epsilon
 */
MatrixN& LinAlg::solve_LS_fast(MatrixN& A, MatrixN& XB, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol)
{
  // verify that A and B are appropriate sizes
  if (A.rows() != XB.rows())
    throw MissizeException();

  // get the dimensionality of A
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  const unsigned k = XB.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy out
  if (m == 0 || n == 0)
  {
    XB.set_zero(n, XB.columns());
    return XB;
  }

  // init matrices for svd computation 
  SAFESTATIC FastThreadable<MatrixN> U, V, workM, workM2;
  SAFESTATIC FastThreadable<VectorN> S;

  // compute the svd
  MatrixN& Ux = U();
  MatrixN& Vx = V();
  VectorN& Sx = S();
  MatrixN& workMx = workM();
  svd(A, Ux, Sx, Vx);
  
  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
    tol = Sx[0] * std::max(m,n) * std::numeric_limits<Real>::epsilon();

  // compute 1/S
  unsigned S_len = Sx.size();

  // A is m x n, B is m x k
  // (L -> R, scaling V)    n^2 + n*min(n,m)*m + nmk [n < m < k, n < k < m]
  // (L -> R, scaling U')   m^2 + n*min(n,m)*m + nmk [m < n < k]
  // (R -> L, scaling U')   m^2 + m^2k + nmk + n^2k  [k < n < m]
  // (R -> L, scaling U'*B) m^2k + min(n,m)*k + n*min(m,n)*k [k < m < n, m < k < n]

  // compute inv(s) 
  for (unsigned i=0; i< S_len; i++)
    Sx[i] = (std::fabs(Sx[i]) > tol) ? (Real) 1.0/Sx[i] : (Real) 0.0;

  // check cases
  // case 1: n is smallest
  if (n < m && n < k)
  {
    // scale n columns of V
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx[i], &Vx(0,i), 1);

    // multiply scaled V by U' = workM
    workMx.resize(n, m);
    CBLAS::gemm(CblasNoTrans, CblasTrans, n, m, n, Vx, n, Ux, m, (Real) 1.0, (Real) 0.0, workMx, n);

    // multiply workM * XB
    Vx.resize(n,k);
    CBLAS::gemm(CblasNoTrans, CblasNoTrans, n, k, m, workMx, n, XB, m, (Real) 1.0, (Real) 0.0, Vx, n);
    XB.copy_from(Vx);
  }
  // case 2: m < n < k
  else if (m < n && n < k)
  {
    // scale columns of U
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(m, Sx[i], &Ux(0,i), 1);

    // multiply V by scaled U' = workM
    workMx.resize(n,m);
    CBLAS::gemm(CblasNoTrans, CblasTrans, n, m, m, Vx, n, Ux, m, (Real) 1.0, (Real) 0.0, workMx, n);

    // multiply workM * XB
    Vx.resize(n,k);
    CBLAS::gemm(CblasNoTrans, CblasNoTrans, n, k, m, workMx, n, XB, m, (Real) 1.0, (Real) 0.0, Vx, n);
    XB.copy_from(Vx);
  }
  // case 3: k < n < m
  else if (k < n && n < m)
  {
    // scale columns of U
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx[i], &Ux(0,i), 1);

    // multiply U' * XB (resulting in n x k matrix)
    workMx.resize(n,k);
    CBLAS::gemm(CblasTrans, CblasNoTrans, n, k, m, Ux, m, XB, m, (Real) 1.0, (Real) 0.0, workMx, n);

    // multiply V * workM
    Vx.mult(workMx, XB);
  }
  // case 4: n is largest
  else
  {
    assert(n >= m && n >= k);

    // scale m columns of V
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(n, Sx[i], &Vx(0,i), 1);

    // multiply U' * XB (resulting in m x k matrix)
    Ux.transpose_mult(XB, workMx);

    // multiply V * workM
    XB.resize(n,k);
    CBLAS::gemm(CblasNoTrans, CblasNoTrans, n, k, m, Vx, n, workMx, m, (Real) 1.0, (Real) 0.0, XB, n);
  }

  return XB;

// NOTE: this is disabled b/c it does not work as well...

/*
  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = std::max(M,N);
  INTEGER INFO;

  // allocate storage for solving the RHS
  const INTEGER XB_rows = std::max(M,NRHS);
  const INTEGER XB_cols = std::max(N,NRHS);
  Real* rhs = new Real[XB_rows*XB_cols];
  std::copy(XB.begin(), XB.end(), rhs);

  // compute
  gelsd_(&M, &N, &NRHS, A.data(), &M, rhs, &LDB, &stol, &INFO);

  // don't check success - SVD algorithm is very stable.. 

  // copy solution into XB
  XB.resize(N, NRHS);
  std::copy(rhs, rhs+(M*NRHS), XB.begin());

  // mark A as destroyed and free memory
  A.resize(0,0);
  delete [] rhs;

  return XB;
*/
}

/// Most robust system of linear equations solver
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param stol the tolerance for determining the rank of A; if stol < 0.0,
 *        stol is computed using machine epsilon
 */
VectorN LinAlg::solve_LS(const MatrixN& A, const VectorN& b, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real tol)
{
  // copy A and b
  MatrixN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call faster solver
  solve_LS_fast(A_copy, x, svd, tol);

  return x;
}

/// Most robust system of linear equations solver
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param stol the tolerance for determining the rank of A; if stol < 0.0,
 *        stol is computed using machine epsilon
 */
MatrixN LinAlg::solve_LS(const MatrixN& A, const MatrixN& B, void (*svd)(MatrixN&, MatrixN&, VectorN&, MatrixN&), Real stol)
{
  MatrixN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call faster solver
  solve_LS_fast(A_copy, X, svd, stol);

  return X;
}

/// Solves the general system Ax = b
/**
 * \param A a square matrix (destroyed on return)
 * \param xb the vector b on input, the vector x on return
 */
VectorN& LinAlg::solve_fast(MatrixN& A, VectorN& xb)
{
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
    return xb;

  // verify that A and b are compatible
  if (A.columns() != xb.size())
    throw MissizeException();

  // setup LAPACK parameters
  INTEGER N = A.rows();
  INTEGER NRHS = 1;
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK (use solving routine that uses LU factorization)
  gesv_(&N, &NRHS, A.data(), &N, &IPIV().front(), xb.data(), &N, &INFO);
 
 // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return xb;
}

/// Solves the general system AX = B
/**
 * \param A a square matrix (destroyed on return)
 * \param XB the matrix B on input, the matrix X on return
 */
MatrixN& LinAlg::solve_fast(MatrixN& A, MatrixN& XB)
{  
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A is not zero size
  if (A.rows() == 0)
    return XB;

  // verify that A and b are compatible
  if (A.columns() != XB.rows())
    throw MissizeException();

  // setup LAPACK parameters
  INTEGER N = A.rows();
  INTEGER NRHS = XB.columns();
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK (use solving routine that uses LU factorization)
  gesv_(&N, &NRHS, A.data(), &N, &IPIV().front(), XB.data(), &N, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Solves the general system Ax = b
/**
 * \param A a square matrix
 * \param b a vector of length A.rows()
 * \return the solution vector of length b.size()
 */
VectorN LinAlg::solve(const MatrixN& A, const VectorN& b)
{  
  // copy A and b
  MatrixN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call faster solver
  solve_fast(A_copy, x);

  return x;
}

/// Solves the general system AX = B
/**
 * \param A a square matrix
 * \param B a matrix with A.rows() rows
 * \return the solution matrix of dimension A.rows() x B.columns()
 */
MatrixN LinAlg::solve(const MatrixN& A, const MatrixN& B)
{  
  // copy A and B
  MatrixN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call faster solver
  solve_fast(A_copy, X);

  return X;
}

/// Computes a givens rotation
void LinAlg::givens(Real a, Real b, Real& c, Real& s)
{
  // setup LAPACK parameters
  Real UNUSED;

  // call LAPACK (use solving routine that uses LU factorization)
  lartg_(&a, &b, &c, &s, &UNUSED);

  // reverse s
  s = -s;
}

/// Computes the givens matrix given a c and s
Matrix2 LinAlg::givens(Real c, Real s)
{
  return Matrix2(c, s, -s, c);
}

/// Computes a householder vector
void LinAlg::householder(Real alpha, const VectorN& x, Real& tau, VectorN& v)
{
  Real s = x.norm_sq();
  v.copy_from(x);
  if (s < (Real) 0.0)
    tau = (Real) 0.0;
  else
  {
    Real t = std::sqrt(alpha*alpha + s);
    Real v_one = (alpha <= (Real) 0.0) ? (alpha-t) : -s/(alpha+t);
    tau = 2*v_one*v_one/(s + v_one*v_one);
    v /= v_one;
  }   
}

/// Updates a QR factorization by a rank-1 update
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 */
void LinAlg::update_QR_rank1(MatrixN& Q, MatrixN& R, const VectorN& u, const VectorN& v)
{
  // apply Givens rotations
}

/// Updates a QR factorization by deleting p columns starting at column idx k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param k the column index to start deleting at
 * \parma p the number of columns to delete 
 */
void LinAlg::update_QR_delete_cols(MatrixN& Q, MatrixN& R, unsigned k, unsigned p)
{
  SAFESTATIC VectorN workv, workv2;
  SAFESTATIC MatrixN workM, workM2, workM3;
  SAFESTATIC vector<Real> c, s, tau;
  SAFESTATIC vector<unsigned> select;
  SAFESTATIC vector<VectorN> V;
  
  const int m = Q.rows();
  const int n = R.columns();
  const int lim = std::min(m-1,n-(int) p);
  assert(k + p <= n);

  // make k one indexed to work with our 1-indexed algorithm
  k++;

  // simplest case
  if (k == n-p+1)
  {
    R.get_sub_mat(0,m,0,k-1,workM);
    R.copy_from(workM);
    Q.get_sub_mat(0,m,0,m,workM);
    Q.copy_from(workM);
    return;
  } 

  // next simplest case
  if (k > std::min(m-1,(int) (n-p)))
  {
    // get relevant indices of R
    select.clear();
    for (int i=1; i<= n; i++)
      if (i < k || i >= k+p)
        select.push_back(i-1);
    R.select_columns(select.begin(), select.end(), workM);
    R.copy_from(workM);
    Q.get_sub_mat(0,m,0,R.rows(),workM);
    Q.copy_from(workM);
    return;
  }

  // setup c, s
  c.resize(std::max(m,n)+1);
  s.resize(std::max(m,n)+1);

  // shift R
  R.get_sub_mat(0,m,k+p-1,n,workM);
  R.set_sub_mat(0,k-1,workM);

  // third simplest case: p = 1 and m >= n
  if (p == 1 && m >= n)
  {
    for (int j=k; j<= n-1; j++)
    {
      // compute Givens rotation
      givens(R(j-1,j-1), R(j,j-1), c[j], s[j]);
      
      // update R
      R(j-1,j-1) = c[j]*R(j-1,j-1) - s[j]*R(j,j-1);
      R.get_sub_mat(j-1,j+1,j-1,n-1,workM);
      givens(c[j], s[j]).transpose_mult(workM, workM2);
      R.set_sub_mat(j-1,j-1,workM2);
    }

    // compute upper triangular part of R
    R.get_sub_mat(0,m,0,n-1,workM);
    for (int i=1; i<= workM.rows(); i++)
      for (int j=1; j<= std::min(i-1,(int) workM.columns()); j++)
        workM(i-1,j-1) = (Real) 0.0;
    R.copy_from(workM);

    // update Q
    for (int j=k; j<= n-1; j++)
    {
      Q.get_sub_mat(0,m,j-1,j+1,workM);
      workM.mult(givens(c[j], s[j]), workM2);
      Q.set_sub_mat(0,j-1, workM2);
    }

    return;
  }

  // householder case (p > 1)
  tau.resize(lim+1);
  V.resize(lim+1);
  for (int j=k; j<= lim; j++)
  {
    // compute householder vector/scalar
    int last = std::min(j+(int) p, m);
    R.get_sub_mat(j, last, j-1, j, workv2);
    householder(R(j-1,j-1), workv2, tau[j], V[j]);

    // update R
    R(j-1,j-1) -= tau[j]*(R(j-1,j-1) + V[j].dot(workv2));
    if (j < n-p)
    {
      R.get_sub_mat(j-1,last,j,n-p, workM);
      workv.resize(V[j].size()+1);
      workv.set_sub_vec(1, V[j]);
      workv[0] = (Real) 1.0;
      workM.transpose_mult(workv, workv2);
      workv *= tau[j];
      VectorN::outer_prod(workv, workv2, &workM2);
      workM -= workM2;
      R.set_sub_mat(j-1,j,workM);
    }
  }

  // setup upper triangular R
  R.get_sub_mat(0,m,0,n-p,workM);
  R.copy_from(workM);
  for (int i=1; i<= workM.rows(); i++)
    for (int j=1; j<= std::min(i-1,(int) workM.columns()); j++)
      R(i-1,j-1) = (Real) 0.0;

  // setup Q
  for (int j=k; j<= lim; j++)
  {
    int last = std::min(j+(int) p,m);
    Q.get_sub_mat(0,m,j-1,last,workM);
    workv.resize(V[j].size()+1);
    workv.set_sub_vec(1, V[j]);
    workv[0] = (Real) 1.0;
    workv2.copy_from(workv) *= tau[j];
    VectorN::outer_prod(workv, workv2, &workM2);
    workM.mult(workM2, workM3);
    workM -= workM3;
    Q.set_sub_mat(0,j-1,workM);    
  }
//  Q.get_sub_mat(0,m,0,std::min(lim+(int) p, m), workM);
  Q.get_sub_mat(0,m,0,m, workM);
  Q.copy_from(workM);
}

/// Updates a QR factorization by inserting one or more columns at column idx k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param U a m x p matrix, destroyed on return
 */
void LinAlg::update_QR_insert_cols(MatrixN& Q, MatrixN& R, MatrixN& U, unsigned k)
{
  SAFESTATIC VectorN workv;
  SAFESTATIC MatrixN workM, workM2, Qu;
  SAFESTATIC MatrixN c, s;

  const int m = Q.rows();
  const int n = R.columns();
  const int p = U.columns();
  assert(U.rows() == m);

  // make k one indexed to work with our 1-indexed algorithm
  k++;

  // setup c, s
  c.resize(std::max(m,n)+1,std::max(m,n)+1);
  s.resize(std::max(m,n)+1,std::max(m,n)+1);

  // setup U
  Q.transpose_mult(U, workM);
  U.copy_from(workM);
  if (m > n+1)
  {
    // do a QR factorization
    U.get_sub_mat(n,m,0,p,workM);
    LinAlg::factor_QR(workM, Qu);
  }

  if (k <= n)
  {
    // zero out the rest with givens, stop at last column of U or last row that 
    // is reached first
    int jstop = std::min((int) p,m-(int) k-2);
    for (int j=1; j<= jstop; j++)
    {
      int istart = std::min(n+j,m);
      int upfirst = std::max(istart-j,1);
      for (int i=istart; i> j; i--)
      {
        givens(U(i-2,j-1),U(i-1,j-1), c(i,j), s(i,j));

        // update U
        U(i-2,j-1) = c(i,j)*U(i-2,j-1) - s(i,j)*U(i-1,j-1);
        if (j < p)
        {
          U.get_sub_mat(i-2,i,j,p,workM);
          givens(c(i,j), s(i,j)).transpose_mult(workM, workM2);
          U.set_sub_mat(i-2,j,workM2);
        }

        // update R
        R.get_sub_mat(i-2,i,upfirst-1,n,workM);
        givens(c(i,j), s(i,j)).transpose_mult(workM, workM2);
        R.set_sub_mat(i-2,upfirst-1,workM2);

        // update one more column next i step
        upfirst--;
      }
    }
  }

  // finish R
  workM.resize(U.rows(), U.columns()+R.columns());
  if (k == 1)
  {
    workM.set_sub_mat(0, 0, U);
    workM.set_sub_mat(0, U.columns(), R);
  }
  else if (k == n+1)
  {
    workM.set_sub_mat(0, 0, R);
    workM.set_sub_mat(0, R.columns(), U);
  }
  else
  {
    R.get_sub_mat(0,m,0,k-1,workM2);
    workM.set_sub_mat(0,0,workM2);
    workM.set_sub_mat(0,workM2.columns(),U);
    R.get_sub_mat(0,m,k-1,n,workM2);
    workM.set_sub_mat(0,U.columns()+k-1,workM2);
  }
  R.copy_from(workM);

  // finally, make R upper triangular
  for (int i=1; i<= (int) R.rows(); i++)
    for (int j=1; j<= std::min((int) R.columns(),i-1); j++)
      R(i-1,j-1) = (Real) 0.0; 

  // compute Q
  if (m > n+1)
  {
    Q.get_sub_mat(0,m,n,m,workM);
    workM.mult(Qu, workM2);
    Q.set_sub_mat(0,n,workM2);
  }
  if (k <= n)
  {
    int jstop = std::min((int) p,m-(int) k-2);
    for (int j=1; j<= jstop; j++)
    {
      int istart = std::min(n+j,m);
      for (int i=istart; i>= j+1; i--)
      {
        Q.get_sub_mat(0,m,i-2,i,workM);
        workM.mult(givens(c(i,j), s(i,j)), workM2);
        Q.set_sub_mat(0,i-2,workM2);
      }
    }
  }
}

/// Updates a QR factorization by inserting a block of rows, starting at index k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param U a p x n matrix (destroyed on return)
 * \param k the index to insert at
 */
void LinAlg::update_QR_insert_rows(MatrixN& Q, MatrixN& R, MatrixN& U, unsigned k)
{
  SAFESTATIC VectorN workv, workv2, workv3, workv4;
  SAFESTATIC MatrixN workM, workM2, workM3, Qu, Ru;
  SAFESTATIC vector<Real> c, s, tau;
  SAFESTATIC vector<VectorN> V;

  const int m = Q.rows();
  const int n = R.columns();
  const int lim = std::min(m,n);
  const int p = U.rows();

  // make k one indexed
  k++;

  // verify that U is the correct size
  assert(U.columns() == n);

  // simplest case: inserting a single row, m >= n
  if (m >= n && p == 1)
  {
    // resize c and s
   c.resize(n+1);
   s.resize(n+1);

    for (int j=1; j<= n; j++)
    {
      // compute Givens rotation
      givens(R(j-1,j-1), U(0, j-1), c[j], s[j]);
      
      // update R
      R(j-1,j-1) = c[j]*R(j-1,j-1) - s[j]*U(0,j-1);

      // update jth row of R and u
      R.get_sub_mat(j-1, j, j, n, workv);
      U.get_sub_mat(0, 1, j, n, workv2);
      workv3.copy_from(workv) *= c[j];
      workv4.copy_from(workv2) *= s[j];
      workv3 -= workv4;
      R.set_sub_mat(j-1,j, workv3, true);
      workv3.copy_from(workv) *= s[j];
      workv4.copy_from(workv2) *= c[j];
      workv3 += workv4;
      U.set_sub_mat(0,j, workv3, true); 
    }

    // setup new R
    workM.resize(R.rows() + 1, R.columns());
    workM.set_sub_mat(0, 0, R);
    BlockIterator b = workM.block_start(R.rows(), R.rows()+1, 0, R.columns());
    std::fill_n(b, R.columns(), (Real) 0.0);
    R.copy_from(workM);

    // compute new Q
    workM.resize(Q.rows()+1, Q.columns()+1);
    workM.set_sub_mat(0,0,Q);
    b = workM.block_start(Q.rows(), Q.rows()+1, 0, Q.columns());
    std::fill_n(b, Q.columns(), (Real) 0.0);
    b = workM.block_start(0, Q.rows(), Q.columns(), Q.columns()+1);
    std::fill_n(b, Q.rows(), (Real) 0.0);
    workM(Q.rows(), Q.columns()) = (Real) 1.0; 
    Q.copy_from(workM);
    if (k != m+1)
    {
      // permute Q
      workM2.resize(Q.rows(), Q.columns());
      Q.get_sub_mat(0,k-1,0,m+1,workM);
      workM2.set_sub_mat(0,0,workM);
      Q.get_row(m, workv);
      workM2.set_row(workM.rows(), workv);
      Q.get_sub_mat(k-1,m,0,m+1, workM);
      workM2.set_sub_mat(k, 0, workM);
      Q.copy_from(workM2);
    }
    for (int j=1; j<= n; j++)
    {
      Q.get_sub_mat(0,m+1,j-1,j,workv);
      Q.get_sub_mat(0,m+1,m,m+1,workv2);
      workv3.copy_from(workv) *= c[j];
      workv4.copy_from(workv2) *= s[j];
      workv3 -= workv4;
      Q.set_sub_mat(0,j-1,workv3,false);
      workv3.copy_from(workv) *= s[j];
      workv4.copy_from(workv2) *= c[j];
      workv3 += workv4;
      Q.set_sub_mat(0,m,workv3,false);
    }

    return;
  }

  // householder case (p > 1 or m < n)
  V.resize(lim+1);
  tau.resize(lim+1);

  for (int j=1; j<= lim; j++)
  {
    // compute householder
    U.get_sub_mat(0,p,j-1,j,workv);
    householder(R(j-1,j-1), workv, tau[j], V[j]);

    // remember old jth row of R
    R.get_sub_mat(j-1,j,j,n,workv);

    // update jth row of R
    R.get_sub_mat(j-1,j,j-1,n,workv2); // workv2 = R(j,j:n)
    workv2 *= ((Real) 1.0 - tau[j]);
    U.get_sub_mat(0,p,j-1,n,workM);   // workM = U(1:p,j:n)
    workM.transpose_mult(V[j], workv3) *= tau[j];
    workv2 -= workv3;
    R.set_sub_mat(j-1,j-1,workv2,true); 

    // update trailing part if U
    if (j < n)
    {
      U.get_sub_mat(0,p,j,n,workM);                 // get X = U(1:p,j+1:n)
      workv2.copy_from(V[j]) *= tau[j];             // tau*V
      VectorN::outer_prod(workv2, V[j], &workM2);   // tau*V*V'
      workM2.mult(workM, workM3);                   // tau*V*V'*X
      workM -= workM3;                              // X -= tau*V*V'*X
      VectorN::outer_prod(workv2, workv, &workM2);  // Y = tau*V*Rj'
      workM -= workM2;                              // X -= Y
      U.set_sub_mat(0,j,workM);
    }
  }

  // update R
  workM.resize(R.rows()+p,n);
  workM.set_sub_mat(0, 0, R);
  BlockIterator b = workM.block_start(R.rows(), workM.rows(), 0, n);
  std::fill_n(b, p*n, (Real) 0.0);
  R.copy_from(workM);
  if (m < n)
  {
    U.get_sub_mat(0, U.rows(), m,n, Ru);
    factor_QR(Ru, Qu);
    R.set_sub_mat(m,m,Ru);
  }

  // update Q
  workM.resize(Q.rows()+p, Q.columns()+p);
  workM.set_sub_mat(0,0,Q);
  b = workM.block_start(0, m, Q.columns(), Q.columns()+p);
  std::fill_n(b, m*p, (Real) 0.0);
  b = workM.block_start(Q.rows(), Q.rows()+p, 0, Q.columns());
  std::fill_n(b, p*Q.columns(), (Real) 0.0);
  b = workM.block_start(Q.rows(), Q.rows()+p, Q.columns(), Q.columns()+p);
  std::fill_n(b, p*p, (Real) 0.0);
  for (unsigned i=0, jj=m, kk=p; i< p; i++, jj++, kk++)
    workM(jj-1,kk-1) = (Real) 1.0;
  Q.copy_from(workM);
  if (k != m+1)
  {
    // permute Q
    Q.get_sub_mat(0,k-1,0,m+p,workM);
    Q.get_sub_mat(m,m+p,0,m+p,workM2);
    Q.get_sub_mat(k-1,m,0,m+p,workM3);
    Q.set_sub_mat(0,0,workM);
    Q.set_sub_mat(workM.rows(),0,workM2);
    Q.set_sub_mat(workM.rows()+workM2.rows(),0,workM3);
  }
  for (int j=1; j<= lim; j++)
  {
    // remember jth column of Q
    Q.get_sub_mat(0,m+p,j-1,j,workv);

    // update jth column
    workv2.copy_from(workv);
    workv2 *= ((Real) 1.0 - tau[j]);
    Q.get_sub_mat(0,m+p,m,m+p,workM);
    workM.mult(V[j], workv3) *= tau[j];
    workv2 -= workv3;
    Q.set_sub_mat(0,j-1,workv2);

    // update m+1:p columns of Qhat
    workv2.copy_from(V[j]) *= tau[j];                 // s = V*tau
    VectorN::outer_prod(workv, workv2, &workM);       // Y = Qhatk * V' * tau
    Q.get_sub_mat(0,m+p,m,m+p, workM2);
    workM -= workM2;
    workM2.mult(V[j], workv);                         // r = X*V
    VectorN::outer_prod(workv, workv2, &workM2);      // Y = X*V*V'*tau
    workM += workM2;
    workM.negate();
    Q.set_sub_mat(0,m,workM);
    if (m < n)
    {
      Q.get_sub_mat(0,m+p,m,m+p,workM);
      workM.mult(Qu, workM2);
      Q.set_sub_mat(0,m,workM2);
    }
  }
}

/// Updates a QR factorization by deleting a block of rows 
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param k the index to start deleting at
 * \param p the number of rows to delete
 */
void LinAlg::update_QR_delete_rows(MatrixN& Q, MatrixN& R, unsigned k, unsigned p)
{
  SAFESTATIC vector<Real> c, s;
  SAFESTATIC VectorN workv, workv2;
  SAFESTATIC MatrixN W, workM, workM2, cc, ss;

  const int m = Q.rows();
  const int n = R.columns();

  assert(k+p <= m);

  // make k 1-indexed
  k++;

  // simplest case: p = 1
  if (p == 1)
  {
    // resize c and s
    c.resize(m);
    s.resize(m);

    Q.get_sub_mat(k-1, k, 0, m, workv);
    for (int j=m-1; j>= 1; j--)
    {
      // compute givens rotation and update q
      givens(workv[j-1], workv[j], c[j], s[j]);
      workv[j-1] = c[j]*workv[j-1] - s[j]*workv[j];
    
      // update R if there is a nonzero row
      if (j <= n)
      {
        R.get_sub_mat(j-1,j+1,j-1,n,workM);
        givens(c[j], s[j]).transpose_mult(workM, workM2);
        R.set_sub_mat(j-1,j-1,workM2);
      }
    }
    R.get_sub_mat(1,m,0,n,workM);
    R.copy_from(workM);

    // compute Q
    if (k != 1)
    {
      Q.get_sub_mat(0,k-1,0,m,workM);
      Q.set_sub_mat(1,0,workM);
    }

    for (int j=m-1; j>= 2; j--)
    {
      Q.get_sub_mat(1,m,j-1,j+1,workM);
      workM.mult(givens(c[j], s[j]), workM2);
      Q.set_sub_mat(1,j-1,workM2);
    }

    // do not need to update 1st column of Q
    Q.get_sub_mat(1,m,0,1,workv); 
    Q.get_sub_mat(1,m,1,2,workv2);
    workv *= s[1];
    workv2 *= c[1];
    workv += workv2;
    Q.set_sub_mat(1,1,workv);
    Q.get_sub_mat(1,m,1,m,workM);
    Q.copy_from(workM);
    return;
  }

  // "standard" case
  cc.resize(p+1,m);
  ss.resize(p+1,m);
  Q.get_sub_mat(k-1,k+p-1,0,m,W);
  for (int i=1; i<= p; i++)
  {
    for (int j=m-1; j>= i; j--)
    {
      givens(W(i-1,j-1), W(i-1,j), cc(i,j), ss(i,j));
      W(i-1,j-1) = W(i-1,j-1) * cc(i,j) - W(i-1,j) * ss(i,j);
      W.get_sub_mat(i,p,j-1,j+1,workM);
      workM.mult(givens(cc(i,j), ss(i,j)), workM2);
      W.set_sub_mat(i,j-1, workM2);

      // update R if there is a nonzero row
      if (j <= n+i-1)
      {
        R.get_sub_mat(j-1,j+1,j-i,n,workM);
        givens(cc(i,j), ss(i,j)).transpose_mult(workM, workM2);
        R.set_sub_mat(j-1,j-i,workM2);
      }
    }
  }

  // compute the new R
  R.get_sub_mat(p,m,0,n,workM);
  R.copy_from(workM);

  // compute the new Q
  if (k != 1)
  {
    Q.get_sub_mat(0,k-1,0,m,workM);
    Q.set_sub_mat(p,0,workM);
  }
  for (int i=1; i<= p; i++)
  {
    for (int j=m-1; j>= i+1; j--)
    {
      Q.get_sub_mat(p,m,j-1,j+1,workM);
      workM.mult(givens(cc(i,j), ss(i,j)), workM2);
      Q.set_sub_mat(p,j-1,workM2);
    }

    Q.get_sub_mat(p,m,i-1,i, workv) *= ss(i,i);
    Q.get_sub_mat(p,m,i,i+1, workv2) *= cc(i,i);
    workv += workv2;
    Q.set_sub_mat(p,i,workv);
  }

  // update Q
  Q.get_sub_mat(p,m,p,m,workM);
  Q.copy_from(workM);
}

/// Performs Gauss-Jordan elimination with partial pivoting on matrix A
MatrixN& LinAlg::gauss_elim(MatrixN& A, Real tol)
{
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  SAFESTATIC FastThreadable<VectorN> rowi, rowmaxi, rowu;
  VectorN& rowix = rowi();
  VectorN& rowmaxix = rowmaxi();
  VectorN& rowux = rowu();

  // look for premature quit
  if (A.rows() == 0 || A.columns() == 0)
    return A;

  // determine the tolerance, if necessary (same way that Octave does)
  if (tol < (Real) 0.0)
    tol = A.norm_inf() * std::max(A.rows(), A.columns()) * std::numeric_limits<Real>::epsilon();

  unsigned i = 0, j = 0;
  while (i < m && j < n)
  {
    // find pivot in column j, starting in row i
    unsigned maxi = i;
    for (unsigned k=i+1; k< m; k++)
      if (std::fabs(A(k,j)) > std::fabs(A(maxi,j)))
        maxi = k;

    if (std::fabs(A(maxi,j)) > tol)
    {
      // swap rows i and maxi
      A.get_row(i, rowix);
      A.get_row(maxi, rowmaxix);
      A.set_row(maxi, rowix);

      // but first, divide each entry in row i by A(i,j)
      rowmaxix /= rowmaxix[j];
      A.set_row(i, rowmaxix);

      for (unsigned u=i+1; u < m; u++)
      {
        A.get_row(u, rowux);
        rowux -= rowmaxix * A(u, j);
        A.set_row(u, rowux);
      }

      i++;
    }
    j++;
  }  

  return A;
}

