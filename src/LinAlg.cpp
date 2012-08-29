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

/// Calls LAPACK function for svd
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

/// Calls LAPACK function for svd
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

/// Calls LAPACK function for svd
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

/// Performs a LDL' factorization of a symmetric, indefinite matrix
/**
 * \param A the matrix A on input; the factorized matrix on output
 */
void LinAlg::factor_LDL(MatrixNN& A, vector<int>& IPIV)
{
  // verify that A is not zero sized
  if (A.size() == 0)
    return;

  // verify that A is symmetric
  #ifndef NDEBUG
  if (!A.is_symmetric(std::sqrt(std::numeric_limits<Real>::epsilon())))
    FILE_LOG(LOG_LINALG) << "LinAlg::factor_LDL() warning: matrix is assymetrical!" << endl;
  #endif

  // resize IPIV
  IPIV.resize(A.size());

  // setup LAPACK args for factorization
  char UPLO = 'L';
  INTEGER N = A.rows();
  INTEGER INFO;

  // get A's data -- we're going to modify it directly
  Real* data = A.begin();

  // alter matrix to put into packed format
  for (unsigned j=0, k=0; j< A.size(); j++)
    for (unsigned i=0; i< A.size(); i++)
      if (i >= j)
         data[k++] = A(i,j);

  // perform the factorization
  sptrf_(&UPLO, &N, A.data(), &IPIV.front(), &INFO);
  assert(INFO >= 0);
}

/// Performs the Cholesky factorization of a matrix
/**
 * \param A the matrix A on input; the factorized (upper triangular) matrix on output
 * \return <b>true</b> if matrix factored successfully, <b>false</b> otherwise
 */
bool LinAlg::factor_chol(MatrixNN& A)
{
  // verify that A is not zero sized
  if (A.size() == 0)
    return true;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.size();
  INTEGER INFO;

  // perform the Cholesky factorization
  potrf_(&UPLO, &N, A.data(), &N, &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    return false;

  // make the matrix upper triangular
  for (unsigned i=1; i< A.size(); i++)
    for (unsigned j=0; j< i; j++)
      A(i,j) = (Real) 0.0;

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

  // call LAPACK to determine work size
  INTEGER LDA = AR.rows();
  INTEGER INFO;
  geqp3_(&M, &N, AR.begin(), &LDA, &PI[0], TAU().begin(), &INFO);
  assert(INFO == 0);

  // correct indices for PI
  for (int i=0; i< N; i++)
    PI[i]--;

  // setup R
  for (unsigned i=0; i< (unsigned) N; i++)
    for (unsigned j=i+1; j< (unsigned) M; j++)
      AR(j,i) = 0.0;

  // resize AR
  AR.resize(std::min(AR.rows(), AR.columns()), AR.columns(), true);
}

/// Performs the QR factorization of a matrix
/**
 * \param AQ the matrix A on input; the matrix R on output 
 * \param Q the matrix Q on output
 */
void LinAlg::factor_QR(MatrixN& AR, MatrixN& Q)
{
  const unsigned m = AR.rows();
  const unsigned n = AR.columns();

  // check for zero sized matrix
  if (m == 0 || n == 0)
    return;

  // determine LAPACK parameters
  INTEGER M = AR.rows();
  INTEGER N = AR.columns();
  INTEGER MINMN = std::min(M, N);
  unsigned min_mn = (unsigned) MINMN; 
 
  // setup tau vector
  SAFESTATIC FastThreadable<VectorN> TAU;
  TAU().resize(min_mn);

  // call LAPACK
  INTEGER LDA = AR.rows();
  INTEGER INFO;
  geqrf_(&M, &N, AR.begin(), &M, TAU().begin(), &INFO);
  assert(INFO == 0);

  // get Q out 
  AR.get_sub_mat(0,AR.rows(),0,std::min(AR.rows(),AR.columns()),Q);
  orgqr_(&M, &M, &MINMN, Q.begin(), &LDA, TAU().begin(), &INFO);
  assert(INFO == 0);
 
  // setup R
  for (unsigned i=0; i< n; i++)
    for (unsigned j=i+1; j< m; j++)
      AR(j,i) = 0.0;

  // resize AR
  AR.resize(std::min(AR.rows(), AR.columns()), AR.columns(), true);
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
MatrixNN& LinAlg::inverse_LU(MatrixNN& M, const vector<int>& IPIV)
{
  if (M.size() == 0)
    return M;

  // call lapack
  INTEGER INFO;
  INTEGER N = M.size();
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
 * \param b the right hand side 
 * \return the solution vector x
 */
VectorN LinAlg::solve_tri(const MatrixNN& A, bool utri, bool transpose_A, const VectorN& b)
{
  VectorN x;
  x = b;
  return solve_tri_fast(A, utri, transpose_A, x);
}

/// Solves a triangular system of linear equations
/**
 * \param A the matrix
 * \param utri if <b>true</b> A is upper triangular (lower triangular otherwise)
 * \param transpose_A if <b>true</b>, solves A'*X = B
 * \param B the right hand side 
 * \return the solution vector X
 */
MatrixN LinAlg::solve_tri(const MatrixNN& A, bool utri, bool transpose_A, const MatrixN& B)
{
  MatrixN X;
  X = B;
  return solve_tri_fast(A, utri, transpose_A, X);
}

/// Solves a triangular system of linear equations
/**
 * \param A the matrix
 * \param utri if <b>true</b> A is upper triangular (lower triangular otherwise)
 * \param transpose_A if <b>true</b>, solves A'*x = b
 * \param xb contains b on entry, x on return
 * \return reference to xb
 */
VectorN& LinAlg::solve_tri_fast(const MatrixNN& A, bool utri, bool transpose_A, VectorN& xb)
{
  if (A.size() != xb.size())
    throw MissizeException();

  if (A.size() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char TRANS = (transpose_A) ? 'T' : 'N';
  char UPLO = (utri) ? 'U' : 'L';
  INTEGER N = A.size();
  INTEGER LDA = A.size();
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
MatrixN& LinAlg::solve_tri_fast(const MatrixNN& A, bool utri, bool transpose_A, MatrixN& XB)
{
  if (A.size() != XB.rows())
    throw MissizeException();

  if (A.size() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char TRANS = (transpose_A) ? 'T' : 'N';
  char UPLO = (utri) ? 'U' : 'L';
  INTEGER N = A.size();
  INTEGER LDA = A.size();
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
VectorN& LinAlg::solve_LDL_fast(const MatrixNN& M, const vector<int>& IPIV, VectorN& xb)
{
  if (M.size() != xb.size())
    throw MissizeException();

  if (M.size() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'L';
  INTEGER N = M.size();
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
MatrixN& LinAlg::solve_LDL_fast(const MatrixNN& M, const vector<int>& IPIV, MatrixN& XB)
{
  if (M.size() != XB.rows())
    throw MissizeException();

  // check for empty matrix
  if (M.size() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.size();
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
VectorN LinAlg::solve_chol(const MatrixNN& M, const VectorN& b)
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
VectorN& LinAlg::solve_chol_fast(const MatrixNN& M, VectorN& xb)
{
  if (M.size() != xb.size())
    throw MissizeException();

  if (M.size() == 0)
    return xb.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.size();
  INTEGER NRHS = 1;
  INTEGER LDA = M.size();
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
MatrixN LinAlg::solve_chol(const MatrixNN& M, const MatrixN& B)
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
MatrixN& LinAlg::solve_chol_fast(const MatrixNN& M, MatrixN& XB)
{
  if (M.size() != XB.rows())
    throw MissizeException();

  // check for empty matrices
  if (M.size() == 0 || XB.columns() == 0)
    return XB.set_zero();

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.size();
  INTEGER NRHS = XB.columns();
  INTEGER LDA = M.size();
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
VectorN LinAlg::solve_LU(const MatrixNN& M, bool transpose, const vector<int>& IPIV, const VectorN& b)
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
VectorN& LinAlg::solve_LU_fast(const MatrixNN& M, bool transpose, const vector<int>& IPIV, VectorN& xb)
{
  if (M.size() == 0)
    return xb.set_zero();

  if (M.size() != xb.size())
    throw MissizeException();

  // setup parameters to LAPACK
  INTEGER N = M.size();
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
MatrixN LinAlg::solve_LU(const MatrixNN& M, bool transpose, const vector<int>& IPIV, const MatrixN& B)
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
MatrixN& LinAlg::solve_LU_fast(const MatrixNN& M, bool transpose, const vector<int>& IPIV, MatrixN& XB)
{
  if (M.size() != XB.rows())
    throw MissizeException();

  // check for empty matrix
  if (M.size() == 0)
    return XB.set_zero();

  // setup parameters to LAPACK
  INTEGER N = M.size();
  INTEGER NRHS = XB.columns();
  char TRANS = (transpose) ? 'T' : 'N';
  
  // call LAPACK
  INTEGER INFO;
  getrs_(&TRANS, &N, &NRHS, (Real*) M.begin(), &N, (INTEGER*) &IPIV.front(), XB.begin(), &N, &INFO);

  return XB;
}

/// Calculates the rank of a matrix
unsigned LinAlg::calc_rank(const MatrixN& A, Real tol)
{
  // look for easy out
  if (A.rows() == 0 || A.columns() == 0)
    return 0;

  // compute the SVD of A
  SAFESTATIC FastThreadable<MatrixNN> U, V;
  SAFESTATIC FastThreadable<VectorN> S;
  MatrixNN& Ux = U();
  MatrixNN& Vx = V();
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
MatrixN& LinAlg::nullspace(const MatrixN& A, MatrixN& nullspace, Real tol)
{
  SAFESTATIC FastThreadable<MatrixNN> U, V;
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
  MatrixNN& Ux = U();
  MatrixNN& Vx = V();
  VectorN& Sx = S();
  svd(A, Ux, Sx, Vx);

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

  return V().get_sub_mat(0,V().rows(),V().columns()-ns,V().columns(), nullspace);
}

/// Computes the condition number of a matrix
Real LinAlg::cond(const MatrixN& A)
{
  MatrixNN U, V;
  VectorN S;
  svd(A, U, S, V);
  return S[0] / S[S.size()-1];
}

/// Iteratively solves a system of equations using the Conjugate-Gradient method
/**
 * \note A must have given form to ensure convergence (e.g., diagonally
 *       dominant, symmetric PD, etc.)
 */
void LinAlg::solve_iterative(const MatrixNN& A, const VectorN& b, VectorN& x, unsigned iter)
{
  VectorN v = b - A*x;
  SAFESTATIC FastThreadable<VectorN> r;
  VectorN& rx = r();
  rx.copy_from(v);
  for (unsigned i=0; i< iter; i++)
  {
    VectorN Av = A*v;
    Real vdot_Av = VectorN::dot(v, Av);
    if (vdot_Av == 0.0)
      return;
    Real alpha = VectorN::dot(v,rx)/vdot_Av;
    x += alpha*v;
    rx -= alpha*Av;
    Real beta = -VectorN::dot(v,A*rx)/vdot_Av;
    v = rx + beta*v;
  }
}  

/// Computes the pseudo-inverse of a square matrix
MatrixNN LinAlg::pseudo_inverse(const MatrixNN& A, Real tol)
{
  MatrixNN result;
  pseudo_inverse(A, result, tol);
  return result;
}

/// Computes the psuedo-inverse of a square matrix
MatrixNN& LinAlg::pseudo_inverse(const MatrixNN& A, MatrixNN& Ainv, Real tol)
{
  // get the dimensionality of A
  unsigned dim = A.rows();

  // init matrices for svd computation 
  SAFESTATIC FastThreadable<MatrixNN> U, V;
  SAFESTATIC FastThreadable<VectorN> S;
  MatrixNN& Ux = U();
  MatrixNN& Vx = V();
  VectorN& Sx = S();

  // compute the svd
  svd(A, Ux, Sx, Vx);
  
  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
    tol = Sx[0] * dim * std::numeric_limits<Real>::epsilon();

  // compute 1/S, according to Numerical Recipes in C, p. 62
  for (unsigned i=0; i< dim; i++)
    Sx[i] = (std::fabs(Sx[i]) > tol) ? (Real) 1.0/Sx[i] : (Real) 0;

  // make a matrix of V * S
  for (unsigned i=0; i< dim; i++)
    for (unsigned j=0; j< dim; j++)
      Vx(i,j) *= Sx[j];

  // compute the pseudo-inverse
  Vx.mult_transpose(Ux, Ainv);

  return Ainv;
}

/// Computes the pseudo-inverse of a non-square matrix
MatrixN LinAlg::pseudo_inverse(const MatrixN& A, Real tol)
{
  MatrixN result;
  pseudo_inverse(A, result, tol);
  return result;
}

/// Computes the psuedo-inverse of a non-square matrix
/**
 * \todo make this faster by optimizing the diagonal matrix multiplication
 */
MatrixN& LinAlg::pseudo_inverse(const MatrixN& A, MatrixN& result, Real tol)
{
  // get the dimensionality of A
  unsigned rows = A.rows();
  unsigned columns = A.columns();

  // check for easy out
  if (rows == 0 || columns == 0)
  {
    result.resize(columns, rows);
    return result;
  }

  // init matrices for svd computation 
  SAFESTATIC FastThreadable<MatrixNN> U, V;
  SAFESTATIC FastThreadable<VectorN> S, row;

  // compute the svd
  MatrixNN& Ux = U();
  MatrixNN& Vx = V();
  VectorN& Sx = S();
  VectorN& rowx = row();
  svd(A, Ux, Sx, Vx);
  
  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
    tol = Sx[0] * std::max(rows,columns) * std::numeric_limits<Real>::epsilon();

  // compute 1/S, according to Numerical Recipes in C, p. 62
  unsigned S_len = Sx.size();

  // create matrix for storing inv(S) * U'
  SAFESTATIC FastThreadable<MatrixN> SM;
  MatrixN& SMx = SM();

  // faster method for computing inv(S) * U'
  for (unsigned i=0; i< S_len; i++)
    Sx[i] = (std::fabs(Sx[i]) > tol) ? (Real) 1.0/Sx[i] : (Real) 0.0;

  // two different methods depending on whether rows > columns
  if (rows < columns)
  {
    // get U' into SM
    SMx.set_zero(columns, rows);
    SMx.set_sub_mat(0, 0, Ux, true);

    // zero remaining rows of SM
    rowx.set_zero(rows);
    for (unsigned i=rows; i< columns; i++)
      SMx.set_row(i, rowx);

    // do the multiplication
    for (unsigned i=0; i< rows; i++)
    {
      SMx.get_row(i, rowx);
      rowx *= Sx[i];
      SMx.set_row(i, rowx);
    }
  }
  else
  {
    // get top rows of U'
    Ux.get_sub_mat(0, rows, 0, columns, SMx, true);
    for (unsigned i=0; i< columns; i++)
    {
      SMx.get_row(i, rowx);
      rowx *= Sx[i];
      SMx.set_row(i, rowx);
    }
  }

// old, proven, slow method
/*
  SM.set_zero(columns,rows);
  for (unsigned i=0; i< S_len; i++)  
    SM(i,i) = (std::fabs(S[i]) > tol) ? (Real) 1.0/S[i] : (Real) 0;
  
  // compute the pseudo-inverse
  return V.mult(SM.mult_transpose(U), result);
*/

  return Vx.mult(SMx, result);
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
bool LinAlg::is_SPSD(const MatrixNN& m, Real tol)
{
  VectorN evals;
  eig_symm(m, evals);

  // make tolerance positive, if it is not already
  tol = std::fabs(tol);

  // check whether all eigenvalues are non-negative to numerical tolerance
  for (unsigned i=0; i< evals.size(); i++)
    if (evals[i] < -tol)
      return false;

  return true;
}

/// Determines whether a matrix is positive-definite
bool LinAlg::is_SPD(const MatrixNN& m, Real tolerance)
{
  // get the eigenvalues of the matrix
  SAFESTATIC FastThreadable<VectorN> evals;
  VectorN& evalsx = evals();
  eig_symm(m, evalsx);

  // make tolerance positive, if it is not already
  tolerance = std::fabs(tolerance);

  // check whether all eigenvalues are positive to numerical tolerance
  for (unsigned i=0; i< evalsx.size(); i++)
    if (evalsx[i] < tolerance)
      return false;

  return true;
}

/// Computes the eigenvalues of the matrix A
/**
 * \param A a matrix
 * \param evals on return, the eigenvalues will be stored here in ascending order
 */
void LinAlg::eig_symm(const MatrixNN& A, VectorN& evals)
{
  // make sure A is not zero sized
  if (A.size() == 0)
  {
    evals.resize(0);
    return;
  }

  // verify that the matrix is symmetric
  #ifndef NDEBUG
  Real max_diff = (Real) 0.0;
  for (unsigned i=0; i< A.size(); i++)
    for (unsigned j=i+1; j< A.size(); j++)
    {
      Real diff = std::fabs(A(i,j) - A(j,i));
      if (diff > max_diff)
        max_diff = diff;
    }
  if (max_diff > NEAR_ZERO)
    std::cerr << "LinAlg::eig_symm() - largest symmetric matrix difference is " << max_diff << std::endl;
  #endif

  // copy A
  SAFESTATIC FastThreadable<MatrixNN> A_copy;
  A_copy() = A;

  // make sure that the eigenvalues array is the proper size
  if (evals.size() != A.rows())
    evals.resize(A.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'N';
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.rows();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, A_copy().data(), &LDA, evals.data(), &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

/// Computes the eigenvalues and eigenvectors of the matrix A
/**
 * \param A a matrix
 * \param evals on return, the eigenvalues will be stored here in ascending order
 * \param evecs on return, the eigenvectors will be stored in order according to the eigenvalues
 */
void LinAlg::eig_symm(const MatrixNN& A, VectorN& evals, MatrixNN& evecs)
{
  // make sure that A is not zero sized
  if (A.size() == 0)
  {
    evals.resize(0);
    evecs.resize(0);
    return;
  }

  // verify that the matrix is symmetric
  #ifndef NDEBUG
  Real max_diff = (Real) 0.0;
  for (unsigned i=0; i< A.size(); i++)
    for (unsigned j=i+1; j< A.size(); j++)
    {
      Real diff = std::fabs(A(i,j) - A(j,i));
      if (diff > max_diff)
        max_diff = diff;
    }
  if (max_diff > NEAR_ZERO)
    std::cerr << "LinAlg::eig_symm() - largest symmetric matrix difference is " << max_diff << std::endl;
  #endif

  // copy A to eigenvectors
  evecs.copy_from(A);

  // make sure that the eigenvalues array is the proper size
  if (evals.size() != A.rows())
    evals.resize(A.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'V';
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.rows();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, evecs.begin(), &LDA, evals.begin(), &INFO);
  assert(INFO == 0);
  
  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

/// Performs singular value decomposition on A
/**
 * The singular value decomposition of A is U*S*V' (' is the transpose 
 * operator); to recompose A, it will be necessary to transpose V before
 * multiplication (i.e., V is returned by the algorithm, not V').
 * Note: passed matrices and vectors U, S, and V are resized as necessary. 
 * \param A the matrix on which the SVD will be performed
 * \param U on output, a A.rows() x A.rows() orthogonal matrix
 * \param S on output, a min(A.rows(), A.columns()) length vector of singular values
 * \param V on output, a A.columns() x A.columns() orthogonal matrix
*/
void LinAlg::svd(const MatrixN& A, MatrixNN& U, VectorN& S, MatrixNN& V)
{
  // make sure that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
  {
    U.resize(A.rows());
    S.resize(0);
    V.resize(A.columns());
    return;
  } 

  // setup U
  if (U.size() != A.rows())
    U.resize(A.rows());

  // setup S
  unsigned minmn = std::min(A.rows(), A.columns());
  if (S.size() != minmn)
    S.resize(minmn);

  // setup V
  if (V.size() != A.columns())
    V.resize(A.columns());

  // setup call to LAPACK
  char JOBZ = 'A';
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  MatrixN A_copy = A;
  INTEGER LDA = A.rows();
  INTEGER LDU = U.rows();
  INTEGER LDVT = V.rows();
  INTEGER INFO;

  // call LAPACK 
  gesdd_(&JOBZ, &M, &N, A_copy.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &INFO);
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
VectorN& LinAlg::solve_symmetric_fast(MatrixNN& A, VectorN& xb)
{
  // make sure that A is not zero sized
  if (A.size() == 0)
    return xb;

  // verify A and b are compatible 
  if (A.columns() != xb.size())
    throw MissizeException();

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = 1;
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK
  sysv_(&UPLO, &N, &NRHS, A.data(), &N, &IPIV().front(), xb.data(), &N, &INFO);

  // mark A as destroyed
  A.resize(0);

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
MatrixN& LinAlg::solve_symmetric_fast(MatrixNN& A, MatrixN& XB)
{
  // make sure that A is not zero sized
  if (A.size() == 0)
    return XB;

  // verify A and b are compatible 
  if (A.columns() != XB.rows())
    throw MissizeException();

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = XB.columns();
  SAFESTATIC FastThreadable<vector<INTEGER> > IPIV;
  IPIV().resize(N);
  INTEGER INFO;

  // call LAPACK
  sysv_(&UPLO, &N, &NRHS, A.data(), &N, &IPIV().front(), XB.data(), &N, &INFO);

  // mark A as destroyed
  A.resize(0);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Solves a symmetric, indefinite square matrix
VectorN LinAlg::solve_symmetric(const MatrixNN& A, const VectorN& b)
{
  // copy A and b
  MatrixNN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call the main solver
  solve_symmetric_fast(A_copy, x);

  return x;
}

/// Solves multiple symmetric, indefinite systems of linear equations
MatrixN LinAlg::solve_symmetric(const MatrixNN& A, const MatrixN& B)
{
  // copy A and B 
  MatrixNN A_copy;
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
MatrixNN& LinAlg::inverse_symmetric(MatrixNN& A)
{
  // verify that A is not zero size
  if (A.size() == 0)
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
  {
    A.resize(0);
    throw SingularException();
  }

  // perform the inversion
  sytri_(&UPLO, &N, A.data(), &N, &IPIV().front(), &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
  {
    A.resize(0);
    throw SingularException();
  }

  // now, make the matrix symmetric
  for (long i=1; i< N; i++)
    for (long j=0; j< i; j++)
      A(i,j) = A(j,i);

  return A;
}

/// Solves a system of equations A*x = b using a symmetric, positive-definite square matrix
/**
 * \param A the matrix coefficients; this matrix will be destroyed on return
 * \param xb on input b, on output x
 */
VectorN& LinAlg::solve_SPD_fast(MatrixNN& A, VectorN& xb)
{
  // verify that A is not zero size
  if (A.size() == 0)
    return xb;

  // verify that A and b are proper size
  if (A.columns() != xb.size())
    throw MissizeException();

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = 1;
  INTEGER INFO;

  // call LAPACK
  posv_(&UPLO, &N, &NRHS, A.data(), &N, xb.data(), &N, &INFO);

  // mark A as destroyed
  A.resize(0);

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
MatrixN& LinAlg::solve_SPD_fast(MatrixNN& A, MatrixN& XB)
{
  // verify that A is not zero size
  if (A.size() == 0)
    return XB;

  // verify that A and b are proper size
  if (A.columns() != XB.rows())
    throw MissizeException();

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;

  // call LAPACK
  posv_(&UPLO, &N, &NRHS, A.data(), &N, XB.data(), &N, &INFO);

  // mark A as destroyed
  A.resize(0);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Solves a system of equations A*X = B using a symmetric, positive-definite square matrix
MatrixN LinAlg::solve_SPD(const MatrixNN& A, const MatrixN& B)
{
  // copy A and B
  MatrixNN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call the main solver
  solve_SPD_fast(A_copy, X);

  return X;
}

/// Solves a system of equations A*x = b using a symmetric, positive-definite square matrix
VectorN LinAlg::solve_SPD(const MatrixNN& A, const VectorN& b)
{
  // copy A and b
  MatrixNN A_copy;
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
MatrixNN& LinAlg::inverse_chol(MatrixNN& A)
{
  // verify that A is not zero sized
  if (A.size() == 0)
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
  for (long i=1; i< N; i++)
    for (long j=0; j< i; j++)
      A(i,j) = A(j,i);

  return A;
}

/// Inverts the symmetric, positive-definite matrix A using Cholesky factorization
/**
 * \param A a square, symmetric positive-definite matrix; contains the inverse
 *        on return
 */
MatrixNN& LinAlg::inverse_SPD(MatrixNN& A)
{
  // verify that A is not zero sized
  if (A.size() == 0)
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
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // now, make the matrix symmetric
  for (long i=1; i< N; i++)
    for (long j=0; j< i; j++)
      A(i,j) = A(j,i);

  return A;
}

/// Inverts the matrix A using LU factorization
/**
 * \param A a square matrix; contains the inverse on return
 */
MatrixNN& LinAlg::inverse(MatrixNN& A)
{
  SAFESTATIC FastThreadable<vector<int> > pivots;
  try
  {
    factor_LU(A, pivots());
    inverse_LU(A, pivots());
  }
  catch (SingularException e)
  {
    A.resize(0);
    throw;
  }
  return A;
}

/// Less robust least squares solver (solves Ax = b)
/**
 * \note this method does not work!
 */
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
 * \param stol the tolerance for determining the rank of A; if stol < 0.0,
 *        stol is computed using machine epsilon
 */
VectorN& LinAlg::solve_LS_fast(MatrixN& A, VectorN& xb, Real stol)
{
  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
    return xb;

  // verify that A and b are appropriate sizes
  if (A.rows() != xb.size())
    throw MissizeException();

  // calculate the pseudo-inverse of A (this is likely not as efficient)
  pseudo_inverse(A, A, stol);
  xb = A * xb;
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
 * \param stol the tolerance for determining the rank of A; if stol < 0.0,
 *        stol is computed using machine epsilon
 */
MatrixN& LinAlg::solve_LS_fast(MatrixN& A, MatrixN& XB, Real stol)
{
  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
    return XB;

  // verify that A and B are appropriate sizes
  if (A.rows() != XB.rows())
    throw MissizeException();

  // calculate the pseudo-inverse of A (this is likely not as efficient)
  pseudo_inverse(A, A, stol);
  XB = A * XB;
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
VectorN LinAlg::solve_LS(const MatrixN& A, const VectorN& b, Real stol)
{
  // copy A and b
  MatrixN A_copy;
  VectorN x;
  A_copy = A;
  x = b;

  // call faster solver
  solve_LS_fast(A_copy, x, stol);

  return x;
}

/// Most robust system of linear equations solver
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param stol the tolerance for determining the rank of A; if stol < 0.0,
 *        stol is computed using machine epsilon
 */
MatrixN LinAlg::solve_LS(const MatrixN& A, const MatrixN& B, Real stol)
{
  MatrixNN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call faster solver
  solve_LS_fast(A_copy, X, stol);

  return X;
}

/// Solves the general system Ax = b
/**
 * \param A a square matrix (destroyed on return)
 * \param xb the vector b on input, the vector x on return
 */
VectorN& LinAlg::solve_fast(MatrixNN& A, VectorN& xb)
{
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
 
  // mark A as destroyed
  A.resize(0);

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
MatrixN& LinAlg::solve_fast(MatrixNN& A, MatrixN& XB)
{  
  // verify that A is not zero size
  if (A.size() == 0)
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

  // mark A as destroyed
  A.resize(0);

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
VectorN LinAlg::solve(const MatrixNN& A, const VectorN& b)
{  
  // copy A and b
  MatrixNN A_copy;
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
MatrixN LinAlg::solve(const MatrixNN& A, const MatrixN& B)
{  
  // copy A and B
  MatrixNN A_copy;
  MatrixN X;
  A_copy = A;
  X = B;

  // call faster solver
  solve_fast(A_copy, X);

  return X;
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

