/*
   ========================================================================
   Definitions and prototypes for LAPACK
   ========================================================================
*/
#ifndef _CLAPACK_H_
#define _CLAPACK_H_

#ifdef BUILD_ARBITRARY_PRECISION
#include <Moby/mpreal.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef int INTEGER;
typedef int LOGICAL;
typedef float SINGLE;
typedef double DOUBLE;
typedef LOGICAL (*__CLPK_L_fp)();
typedef long int __CLPK_ftnlen;

typedef struct { SINGLE r, i; } COMPLEX;
typedef struct { DOUBLE r, i; } DOUBLECOMPLEX;

// *********************************************************************
// arbitrary precision routines start here
// *********************************************************************
 #ifdef BUILD_ARBITRARY_PRECISION
 int aorgqr_(INTEGER *m, INTEGER *n, INTEGER *k, mpfr::mpreal *
	a, INTEGER *lda, mpfr::mpreal* tau, mpfr::mpreal *work, INTEGER *lwork, 
	INTEGER *info);

 int agetrs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	mpfr::mpreal *a, INTEGER *lda, INTEGER *ipiv, mpfr::mpreal *b, INTEGER *
	ldb, INTEGER *info);

 int agetri_(INTEGER *n, mpfr::mpreal *a, INTEGER *lda, INTEGER 
	*ipiv, mpfr::mpreal *work, INTEGER *lwork, INTEGER *info);

 int agetrf_(INTEGER *m, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, INTEGER *ipiv, INTEGER *info);

 int asytrf_(char *uplo, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, INTEGER *ipiv, mpfr::mpreal *work, INTEGER *lwork, INTEGER *info);

 int asytri_(char *uplo, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, INTEGER *ipiv, mpfr::mpreal *work, INTEGER *info);

 int asysv_(char *uplo, INTEGER *n, INTEGER *nrhs, mpfr::mpreal 
	*a, INTEGER *lda, INTEGER *ipiv, mpfr::mpreal *b, INTEGER *ldb, 
	mpfr::mpreal *work, INTEGER *lwork, INTEGER *info);
 
 int aposv_(char *uplo, INTEGER *n, INTEGER *nrhs, mpfr::mpreal 
	*a, INTEGER *lda, mpfr::mpreal *b, INTEGER *ldb, INTEGER *info);

 int apotrf_(char *uplo, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, INTEGER *info);
 
 int apotrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	mpfr::mpreal *a, INTEGER *lda, mpfr::mpreal *b, INTEGER *ldb, INTEGER *
	info);

 int asptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	mpfr::mpreal *ap, INTEGER *ipiv, mpfr::mpreal *b, INTEGER *ldb, INTEGER *
	info);

 int asptrf_(char *uplo, INTEGER *n, mpfr::mpreal *ap, INTEGER *
	ipiv, INTEGER *info);
 
 int aormqr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, mpfr::mpreal *a, INTEGER *lda, mpfr::mpreal *tau, mpfr::mpreal *
	c__, INTEGER *ldc, mpfr::mpreal *work, INTEGER *lwork, INTEGER *info);

 int atrtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, mpfr::mpreal *a, INTEGER *lda, mpfr::mpreal *b, INTEGER *
	ldb, INTEGER *info);

 int agelsd_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	mpfr::mpreal *a, INTEGER *lda, mpfr::mpreal *b, INTEGER *ldb, mpfr::mpreal *
	s, mpfr::mpreal *rcond, INTEGER *rank, mpfr::mpreal *work, INTEGER *lwork,
	 INTEGER *iwork, INTEGER *info);
 

 int agesv_(INTEGER *n, INTEGER *nrhs, mpfr::mpreal *a, INTEGER 
	*lda, INTEGER *ipiv, mpfr::mpreal *b, INTEGER *ldb, INTEGER *info);
 
 int apotri_(char *uplo, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, INTEGER *info);
 
 int asyevd_(char *jobz, char *uplo, INTEGER *n, mpfr::mpreal *
	a, INTEGER *lda, mpfr::mpreal *w, mpfr::mpreal *work, INTEGER *lwork, 
	INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int agesdd_(char *jobz, INTEGER *m, INTEGER *n, mpfr::mpreal *
	a, INTEGER *lda, mpfr::mpreal *s, mpfr::mpreal *u, INTEGER *ldu, 
	mpfr::mpreal *vt, INTEGER *ldvt, mpfr::mpreal *work, INTEGER *lwork, 
	INTEGER *iwork, INTEGER *info);
 

 int ageqrf_(INTEGER *m, INTEGER *n, mpfr::mpreal *a, INTEGER *
	lda, mpfr::mpreal *tau, mpfr::mpreal *work, INTEGER *lwork, INTEGER *info);
 
 #endif
// *********************************************************************
// arbitrary precision routines end here
// *********************************************************************
 
 int cbdsqr_(char *uplo, INTEGER *n, INTEGER *ncvt, INTEGER *
	nru, INTEGER *ncc, SINGLE *d__, SINGLE *e, COMPLEX *vt, INTEGER *ldvt, 
	COMPLEX *u, INTEGER *ldu, COMPLEX *c__, INTEGER *ldc, SINGLE *rwork, 
	INTEGER *info);
 
 int cgbbrd_(char *vect, INTEGER *m, INTEGER *n, INTEGER *ncc,
	 INTEGER *kl, INTEGER *ku, COMPLEX *ab, INTEGER *ldab, SINGLE *d__, 
	SINGLE *e, COMPLEX *q, INTEGER *ldq, COMPLEX *pt, INTEGER *ldpt, 
	COMPLEX *c__, INTEGER *ldc, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgbcon_(char *norm, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 COMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgbequ_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 COMPLEX *ab, INTEGER *ldab, SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE 
	*colcnd, SINGLE *amax, INTEGER *info);
 
 int cgbrfs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *afb, INTEGER *
	ldafb, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *
	ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, INTEGER *
	info);
 
 int cgbsv_(INTEGER *n, INTEGER *kl, INTEGER *ku, INTEGER *
	nrhs, COMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, COMPLEX *b, INTEGER *
	ldb, INTEGER *info);
 
 int cgbsvx_(char *fact, char *trans, INTEGER *n, INTEGER *kl,
	 INTEGER *ku, INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *afb,
	 INTEGER *ldafb, INTEGER *ipiv, char *equed, SINGLE *r__, SINGLE *c__, 
	COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE 
	*ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgbtf2_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 COMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int cgbtrf_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 COMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int cgbtrs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, COMPLEX 
	*b, INTEGER *ldb, INTEGER *info);
 
 int cgebak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, SINGLE *scale, INTEGER *m, COMPLEX *v, INTEGER *ldv, 
	INTEGER *info);
 
 int cgebal_(char *job, INTEGER *n, COMPLEX *a, INTEGER *lda, 
	INTEGER *ilo, INTEGER *ihi, SINGLE *scale, INTEGER *info);
 
 int cgebd2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *d__, SINGLE *e, COMPLEX *tauq, COMPLEX *taup, COMPLEX *work, 
	INTEGER *info);
 
 int cgebrd_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *d__, SINGLE *e, COMPLEX *tauq, COMPLEX *taup, COMPLEX *work, 
	INTEGER *lwork, INTEGER *info);
 
 int cgecon_(char *norm, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *anorm, SINGLE *rcond, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgeequ_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *colcnd, SINGLE *amax, 
	INTEGER *info);
 
 int cgees_(char *jobvs, char *sort, __CLPK_L_fp select, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, INTEGER *sdim, COMPLEX *w, COMPLEX *vs, 
	INTEGER *ldvs, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, LOGICAL *
	bwork, INTEGER *info);
 
 int cgeesx_(char *jobvs, char *sort, __CLPK_L_fp select, char *
	sense, INTEGER *n, COMPLEX *a, INTEGER *lda, INTEGER *sdim, COMPLEX *
	w, COMPLEX *vs, INTEGER *ldvs, SINGLE *rconde, SINGLE *rcondv, COMPLEX *
	work, INTEGER *lwork, SINGLE *rwork, LOGICAL *bwork, INTEGER *info);
 
 int cgeev_(char *jobvl, char *jobvr, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, COMPLEX *w, COMPLEX *vl, INTEGER *ldvl, COMPLEX *vr, 
	INTEGER *ldvr, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *
	info);
 
 int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *w, COMPLEX *vl, 
	INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, INTEGER *ilo, INTEGER *ihi,
	 SINGLE *scale, SINGLE *abnrm, SINGLE *rconde, SINGLE *rcondv, COMPLEX *work, 
	INTEGER *lwork, SINGLE *rwork, INTEGER *info);
 
 int cgegs_(char *jobvsl, char *jobvsr, INTEGER *n, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *alpha, COMPLEX *
	beta, COMPLEX *vsl, INTEGER *ldvsl, COMPLEX *vsr, INTEGER *ldvsr, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *info);
 
 int cgegv_(char *jobvl, char *jobvr, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *alpha, COMPLEX *beta,
	 COMPLEX *vl, INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, COMPLEX *
	work, INTEGER *lwork, SINGLE *rwork, INTEGER *info);
 
 int cgehd2_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, COMPLEX *
	a, INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cgehrd_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, COMPLEX *
	a, INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER 
	*info);
 
 int cgelq2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cgelqf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgels_(char *trans, INTEGER *m, INTEGER *n, INTEGER *
	nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int cgelsx_(INTEGER *m, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *jpvt, SINGLE *rcond,
	 INTEGER *rank, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgelsy_(INTEGER *m, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *jpvt, SINGLE *rcond,
	 INTEGER *rank, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *
	info);
 
 int cgeql2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cgeqlf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgeqp3_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *jpvt, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, SINGLE *
	rwork, INTEGER *info);
 
 int cgeqpf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *jpvt, COMPLEX *tau, COMPLEX *work, SINGLE *rwork, INTEGER *
	info);
 
 int cgeqr2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cgeqrf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgerfs_(char *trans, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *ipiv, COMPLEX *
	b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgerq2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cgerqf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgesc2_(INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *
	rhs, INTEGER *ipiv, INTEGER *jpiv, SINGLE *scale);
 
 int cgesv_(INTEGER *n, INTEGER *nrhs, COMPLEX *a, INTEGER *
	lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cgesvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *
	ipiv, char *equed, SINGLE *r__, SINGLE *c__, COMPLEX *b, INTEGER *ldb, 
	COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgetc2_(INTEGER *n, COMPLEX *a, INTEGER *lda, INTEGER *
	ipiv, INTEGER *jpiv, INTEGER *info);
 
 int cgetf2_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, INTEGER *info);
 
 int cgetrf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, INTEGER *info);
 
 int cgetri_(INTEGER *n, COMPLEX *a, INTEGER *lda, INTEGER *
	ipiv, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgetrs_(char *trans, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int cggbak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, SINGLE *lscale, SINGLE *rscale, INTEGER *m, COMPLEX *v, 
	INTEGER *ldv, INTEGER *info);
 
 int cggbal_(char *job, INTEGER *n, COMPLEX *a, INTEGER *lda, 
	COMPLEX *b, INTEGER *ldb, INTEGER *ilo, INTEGER *ihi, SINGLE *lscale, 
	SINGLE *rscale, SINGLE *work, INTEGER *info);
 
 int cgges_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	selctg, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *
	ldb, INTEGER *sdim, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vsl, 
	INTEGER *ldvsl, COMPLEX *vsr, INTEGER *ldvsr, COMPLEX *work, INTEGER *
	lwork, SINGLE *rwork, LOGICAL *bwork, INTEGER *info);
 
 int cggesx_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	selctg, char *sense, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b,
	 INTEGER *ldb, INTEGER *sdim, COMPLEX *alpha, COMPLEX *beta, COMPLEX *
	vsl, INTEGER *ldvsl, COMPLEX *vsr, INTEGER *ldvsr, SINGLE *rconde, SINGLE 
	*rcondv, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *iwork, 
	INTEGER *liwork, LOGICAL *bwork, INTEGER *info);
 
 int cggev_(char *jobvl, char *jobvr, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *alpha, COMPLEX *beta,
	 COMPLEX *vl, INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, COMPLEX *
	work, INTEGER *lwork, SINGLE *rwork, INTEGER *info);
 
 int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb,
	 COMPLEX *alpha, COMPLEX *beta, COMPLEX *vl, INTEGER *ldvl, COMPLEX *
	vr, INTEGER *ldvr, INTEGER *ilo, INTEGER *ihi, SINGLE *lscale, SINGLE *
	rscale, SINGLE *abnrm, SINGLE *bbnrm, SINGLE *rconde, SINGLE *rcondv, COMPLEX 
	*work, INTEGER *lwork, SINGLE *rwork, INTEGER *iwork, LOGICAL *bwork, 
	INTEGER *info);
 
 int cggglm_(INTEGER *n, INTEGER *m, INTEGER *p, COMPLEX *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *d__, COMPLEX *x, 
	COMPLEX *y, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cgghrd_(char *compq, char *compz, INTEGER *n, INTEGER *
	ilo, INTEGER *ihi, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb,
	 COMPLEX *q, INTEGER *ldq, COMPLEX *z__, INTEGER *ldz, INTEGER *info);
 
 int cgglse_(INTEGER *m, INTEGER *n, INTEGER *p, COMPLEX *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *c__, COMPLEX *d__, 
	COMPLEX *x, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cggqrf_(INTEGER *n, INTEGER *m, INTEGER *p, COMPLEX *a, 
	INTEGER *lda, COMPLEX *taua, COMPLEX *b, INTEGER *ldb, COMPLEX *taub, 
	COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cggrqf_(INTEGER *m, INTEGER *p, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, COMPLEX *taua, COMPLEX *b, INTEGER *ldb, COMPLEX *taub, 
	COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cggsvd_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *n, INTEGER *p, INTEGER *k, INTEGER *l, COMPLEX *a, INTEGER *
	lda, COMPLEX *b, INTEGER *ldb, SINGLE *alpha, SINGLE *beta, COMPLEX *u, 
	INTEGER *ldu, COMPLEX *v, INTEGER *ldv, COMPLEX *q, INTEGER *ldq, 
	COMPLEX *work, SINGLE *rwork, INTEGER *iwork, INTEGER *info);
 
 int cggsvp_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER 
	*ldb, SINGLE *tola, SINGLE *tolb, INTEGER *k, INTEGER *l, COMPLEX *u, 
	INTEGER *ldu, COMPLEX *v, INTEGER *ldv, COMPLEX *q, INTEGER *ldq, 
	INTEGER *iwork, SINGLE *rwork, COMPLEX *tau, COMPLEX *work, INTEGER *
	info);
 
 int cgtcon_(char *norm, INTEGER *n, COMPLEX *dl, COMPLEX *
	d__, COMPLEX *du, COMPLEX *du2, INTEGER *ipiv, SINGLE *anorm, SINGLE *
	rcond, COMPLEX *work, INTEGER *info);
 
 int cgtrfs_(char *trans, INTEGER *n, INTEGER *nrhs, COMPLEX *
	dl, COMPLEX *d__, COMPLEX *du, COMPLEX *dlf, COMPLEX *df, COMPLEX *
	duf, COMPLEX *du2, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *
	x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int cgtsv_(INTEGER *n, INTEGER *nrhs, COMPLEX *dl, COMPLEX *
	d__, COMPLEX *du, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cgtsvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, COMPLEX *dl, COMPLEX *d__, COMPLEX *du, COMPLEX *dlf, COMPLEX *
	df, COMPLEX *duf, COMPLEX *du2, INTEGER *ipiv, COMPLEX *b, INTEGER *
	ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cgttrf_(INTEGER *n, COMPLEX *dl, COMPLEX *d__, COMPLEX *
	du, COMPLEX *du2, INTEGER *ipiv, INTEGER *info);
 
 int cgttrs_(char *trans, INTEGER *n, INTEGER *nrhs, COMPLEX *
	dl, COMPLEX *d__, COMPLEX *du, COMPLEX *du2, INTEGER *ipiv, COMPLEX *
	b, INTEGER *ldb, INTEGER *info);
 
 int cgtts2_(INTEGER *itrans, INTEGER *n, INTEGER *nrhs, 
	COMPLEX *dl, COMPLEX *d__, COMPLEX *du, COMPLEX *du2, INTEGER *ipiv, 
	COMPLEX *b, INTEGER *ldb);
 
 int chbev_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	COMPLEX *ab, INTEGER *ldab, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int chbevd_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	COMPLEX *ab, INTEGER *ldab, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *lrwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int chbevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *kd, COMPLEX *ab, INTEGER *ldab, COMPLEX *q, INTEGER *ldq, 
	SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, INTEGER *
	m, SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *work, SINGLE *rwork, 
	INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int chbgst_(char *vect, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, COMPLEX *ab, INTEGER *ldab, COMPLEX *bb, INTEGER *ldbb, 
	COMPLEX *x, INTEGER *ldx, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int chbgv_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, COMPLEX *ab, INTEGER *ldab, COMPLEX *bb, INTEGER *ldbb, 
	SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int chbgvx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *ka, INTEGER *kb, COMPLEX *ab, INTEGER *ldab, COMPLEX *bb, 
	INTEGER *ldbb, COMPLEX *q, INTEGER *ldq, SINGLE *vl, SINGLE *vu, INTEGER *
	il, INTEGER *iu, SINGLE *abstol, INTEGER *m, SINGLE *w, COMPLEX *z__, 
	INTEGER *ldz, COMPLEX *work, SINGLE *rwork, INTEGER *iwork, INTEGER *
	ifail, INTEGER *info);
 
 int chbtrd_(char *vect, char *uplo, INTEGER *n, INTEGER *kd, 
	COMPLEX *ab, INTEGER *ldab, SINGLE *d__, SINGLE *e, COMPLEX *q, INTEGER *
	ldq, COMPLEX *work, INTEGER *info);
 
 int checon_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, COMPLEX *work, INTEGER *
	info);
 
 int cheev_(char *jobz, char *uplo, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, SINGLE *w, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, 
	INTEGER *info);
 
 int cheevd_(char *jobz, char *uplo, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, SINGLE *w, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, 
	INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int cheevr_(char *jobz, char *range, char *uplo, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *
	iu, SINGLE *abstol, INTEGER *m, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	INTEGER *isuppz, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *
	lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int cheevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *
	iu, SINGLE *abstol, INTEGER *m, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *iwork, INTEGER *
	ifail, INTEGER *info);
 
 int chegs2_(INTEGER *itype, char *uplo, INTEGER *n, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int chegst_(INTEGER *itype, char *uplo, INTEGER *n, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int chegv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, SINGLE *w, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *info);
 
 int chegvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, SINGLE *w, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *lrwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int chegvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, INTEGER *
	m, SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *work, INTEGER *lwork,
	 SINGLE *rwork, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int cherfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *ipiv, COMPLEX *
	b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int chesv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *a,
	 INTEGER *lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *work,
	 INTEGER *lwork, INTEGER *info);
 
 int chesvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *
	ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond,
	 SINGLE *ferr, SINGLE *berr, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, 
	INTEGER *info);
 
 int chetf2_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, INTEGER *info);
 
 int chetrd_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *d__, SINGLE *e, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, 
	INTEGER *info);
 
 int chetrf_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int chetri_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, COMPLEX *work, INTEGER *info);
 
 int chetrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int chgeqz_(char *job, char *compq, char *compz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, COMPLEX *a, INTEGER *lda, COMPLEX *b, 
	INTEGER *ldb, COMPLEX *alpha, COMPLEX *beta, COMPLEX *q, INTEGER *ldq,
	 COMPLEX *z__, INTEGER *ldz, COMPLEX *work, INTEGER *lwork, SINGLE *
	rwork, INTEGER *info);
 
 int chpcon_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, SINGLE *anorm, SINGLE *rcond, COMPLEX *work, INTEGER *info);
 
 int chpev_(char *jobz, char *uplo, INTEGER *n, COMPLEX *ap, 
	SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int chpevd_(char *jobz, char *uplo, INTEGER *n, COMPLEX *ap, 
	SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *work, INTEGER *lwork, 
	SINGLE *rwork, INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, 
	INTEGER *info);
 
 int chpevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	COMPLEX *ap, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *
	abstol, INTEGER *m, SINGLE *w, COMPLEX *z__, INTEGER *ldz, COMPLEX *
	work, SINGLE *rwork, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int chpgst_(INTEGER *itype, char *uplo, INTEGER *n, COMPLEX *
	ap, COMPLEX *bp, INTEGER *info);
 
 int chpgv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, COMPLEX *ap, COMPLEX *bp, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int chpgvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, COMPLEX *ap, COMPLEX *bp, SINGLE *w, COMPLEX *z__, INTEGER *ldz, 
	COMPLEX *work, INTEGER *lwork, SINGLE *rwork, INTEGER *lrwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int chpgvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, COMPLEX *ap, COMPLEX *bp, SINGLE *vl, SINGLE *vu, 
	INTEGER *il, INTEGER *iu, SINGLE *abstol, INTEGER *m, SINGLE *w, COMPLEX *
	z__, INTEGER *ldz, COMPLEX *work, SINGLE *rwork, INTEGER *iwork, 
	INTEGER *ifail, INTEGER *info);
 
 int chprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, COMPLEX *afp, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *x,
	 INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int chpsv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int chpsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *ap, COMPLEX *afp, INTEGER *ipiv, COMPLEX *b, INTEGER *
	ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int chptrd_(char *uplo, INTEGER *n, COMPLEX *ap, SINGLE *d__, 
	SINGLE *e, COMPLEX *tau, INTEGER *info);
 
 int chptrf_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, INTEGER *info);
 
 int chptri_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, COMPLEX *work, INTEGER *info);
 
 int chptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int chsein_(char *side, char *eigsrc, char *initv, LOGICAL *
	select, INTEGER *n, COMPLEX *h__, INTEGER *ldh, COMPLEX *w, COMPLEX *
	vl, INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, INTEGER *mm, INTEGER *
	m, COMPLEX *work, SINGLE *rwork, INTEGER *ifaill, INTEGER *ifailr, 
	INTEGER *info);
 
 int chseqr_(char *job, char *compz, INTEGER *n, INTEGER *ilo,
	 INTEGER *ihi, COMPLEX *h__, INTEGER *ldh, COMPLEX *w, COMPLEX *z__, 
	INTEGER *ldz, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int clabrd_(INTEGER *m, INTEGER *n, INTEGER *nb, COMPLEX *a, 
	INTEGER *lda, SINGLE *d__, SINGLE *e, COMPLEX *tauq, COMPLEX *taup, 
	COMPLEX *x, INTEGER *ldx, COMPLEX *y, INTEGER *ldy);
 
 int clacgv_(INTEGER *n, COMPLEX *x, INTEGER *incx);
 
 int clacon_(INTEGER *n, COMPLEX *v, COMPLEX *x, SINGLE *est, 
	INTEGER *kase);
 
 int clacp2_(char *uplo, INTEGER *m, INTEGER *n, SINGLE *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb);
 
 int clacpy_(char *uplo, INTEGER *m, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, COMPLEX *b, INTEGER *ldb);
 
 int clacrm_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *b, INTEGER *ldb, COMPLEX *c__, INTEGER *ldc, SINGLE *rwork);
 
 int clacrt_(INTEGER *n, COMPLEX *cx, INTEGER *incx, COMPLEX *
	cy, INTEGER *incy, COMPLEX *c__, COMPLEX *s);
 
 int claed0_(INTEGER *qsiz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	COMPLEX *q, INTEGER *ldq, COMPLEX *qstore, INTEGER *ldqs, SINGLE *rwork,
	 INTEGER *iwork, INTEGER *info);
 
 int claed7_(INTEGER *n, INTEGER *cutpnt, INTEGER *qsiz, 
	INTEGER *tlvls, INTEGER *curlvl, INTEGER *curpbm, SINGLE *d__, COMPLEX *
	q, INTEGER *ldq, SINGLE *rho, INTEGER *indxq, SINGLE *qstore, INTEGER *
	qptr, INTEGER *prmptr, INTEGER *perm, INTEGER *givptr, INTEGER *
	givcol, SINGLE *givnum, COMPLEX *work, SINGLE *rwork, INTEGER *iwork, 
	INTEGER *info);
 
 int claed8_(INTEGER *k, INTEGER *n, INTEGER *qsiz, COMPLEX *
	q, INTEGER *ldq, SINGLE *d__, SINGLE *rho, INTEGER *cutpnt, SINGLE *z__, 
	SINGLE *dlamda, COMPLEX *q2, INTEGER *ldq2, SINGLE *w, INTEGER *indxp, 
	INTEGER *indx, INTEGER *indxq, INTEGER *perm, INTEGER *givptr, 
	INTEGER *givcol, SINGLE *givnum, INTEGER *info);
 
 int claein_(LOGICAL *rightv, LOGICAL *noinit, INTEGER *n, 
	COMPLEX *h__, INTEGER *ldh, COMPLEX *w, COMPLEX *v, COMPLEX *b, 
	INTEGER *ldb, SINGLE *rwork, SINGLE *eps3, SINGLE *smlnum, INTEGER *info);
 
 int claesy_(COMPLEX *a, COMPLEX *b, COMPLEX *c__, COMPLEX *
	rt1, COMPLEX *rt2, COMPLEX *evscal, COMPLEX *cs1, COMPLEX *sn1);
 
 int claev2_(COMPLEX *a, COMPLEX *b, COMPLEX *c__, SINGLE *rt1, 
	SINGLE *rt2, SINGLE *cs1, COMPLEX *sn1);
 
 int clags2_(LOGICAL *upper, SINGLE *a1, COMPLEX *a2, SINGLE *a3, 
	SINGLE *b1, COMPLEX *b2, SINGLE *b3, SINGLE *csu, COMPLEX *snu, SINGLE *csv, 
	COMPLEX *snv, SINGLE *csq, COMPLEX *snq);
 
 int clagtm_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *
	alpha, COMPLEX *dl, COMPLEX *d__, COMPLEX *du, COMPLEX *x, INTEGER *
	ldx, SINGLE *beta, COMPLEX *b, INTEGER *ldb);
 
 int clahef_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 COMPLEX *a, INTEGER *lda, INTEGER *ipiv, COMPLEX *w, INTEGER *ldw, 
	INTEGER *info);
 
 int clahqr_(LOGICAL *wantt, LOGICAL *wantz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, COMPLEX *h__, INTEGER *ldh, COMPLEX *w, 
	INTEGER *iloz, INTEGER *ihiz, COMPLEX *z__, INTEGER *ldz, INTEGER *
	info);
 
 int clahrd_(INTEGER *n, INTEGER *k, INTEGER *nb, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *t, INTEGER *ldt, COMPLEX *y, 
	INTEGER *ldy);
 
 int claic1_(INTEGER *job, INTEGER *j, COMPLEX *x, SINGLE *sest,
	 COMPLEX *w, COMPLEX *gamma, SINGLE *sestpr, COMPLEX *s, COMPLEX *c__);
 
 int clals0_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *nrhs, COMPLEX *b, INTEGER *ldb, COMPLEX *bx, 
	INTEGER *ldbx, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, 
	INTEGER *ldgcol, SINGLE *givnum, INTEGER *ldgnum, SINGLE *poles, SINGLE *
	difl, SINGLE *difr, SINGLE *z__, INTEGER *k, SINGLE *c__, SINGLE *s, SINGLE *
	rwork, INTEGER *info);
 
 int clalsa_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *nrhs, COMPLEX *b, INTEGER *ldb, COMPLEX *bx, INTEGER *ldbx, 
	SINGLE *u, INTEGER *ldu, SINGLE *vt, INTEGER *k, SINGLE *difl, SINGLE *difr, 
	SINGLE *z__, SINGLE *poles, INTEGER *givptr, INTEGER *givcol, INTEGER *
	ldgcol, INTEGER *perm, SINGLE *givnum, SINGLE *c__, SINGLE *s, SINGLE *rwork, 
	INTEGER *iwork, INTEGER *info);
 
 int clapll_(INTEGER *n, COMPLEX *x, INTEGER *incx, COMPLEX *
	y, INTEGER *incy, SINGLE *ssmin);
 
 int clapmt_(LOGICAL *forwrd, INTEGER *m, INTEGER *n, COMPLEX 
	*x, INTEGER *ldx, INTEGER *k);
 
 int claqgb_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 COMPLEX *ab, INTEGER *ldab, SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE 
	*colcnd, SINGLE *amax, char *equed);
 
 int claqge_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *colcnd, SINGLE *amax, char *
	equed);
 
 int claqhb_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int claqhe_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int claqhp_(char *uplo, INTEGER *n, COMPLEX *ap, SINGLE *s, 
	SINGLE *scond, SINGLE *amax, char *equed);
 
 int claqp2_(INTEGER *m, INTEGER *n, INTEGER *offset, COMPLEX 
	*a, INTEGER *lda, INTEGER *jpvt, COMPLEX *tau, SINGLE *vn1, SINGLE *vn2, 
	COMPLEX *work);
 
 int claqps_(INTEGER *m, INTEGER *n, INTEGER *offset, INTEGER 
	*nb, INTEGER *kb, COMPLEX *a, INTEGER *lda, INTEGER *jpvt, COMPLEX *
	tau, SINGLE *vn1, SINGLE *vn2, COMPLEX *auxv, COMPLEX *f, INTEGER *ldf);
 
 int claqsb_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int claqsp_(char *uplo, INTEGER *n, COMPLEX *ap, SINGLE *s, 
	SINGLE *scond, SINGLE *amax, char *equed);
 
 int claqsy_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int clar1v_(INTEGER *n, INTEGER *b1, INTEGER *bn, SINGLE *
	sigma, SINGLE *d__, SINGLE *l, SINGLE *ld, SINGLE *lld, SINGLE *gersch, COMPLEX 
	*z__, SINGLE *ztz, SINGLE *mingma, INTEGER *r__, INTEGER *isuppz, SINGLE *
	work);
 
 int clar2v_(INTEGER *n, COMPLEX *x, COMPLEX *y, COMPLEX *z__,
	 INTEGER *incx, SINGLE *c__, COMPLEX *s, INTEGER *incc);
 
 int clarcm_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	COMPLEX *b, INTEGER *ldb, COMPLEX *c__, INTEGER *ldc, SINGLE *rwork);
 
 int clarf_(char *side, INTEGER *m, INTEGER *n, COMPLEX *v, 
	INTEGER *incv, COMPLEX *tau, COMPLEX *c__, INTEGER *ldc, COMPLEX *
	work);
 
 int clarfb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *v, INTEGER *ldv, 
	COMPLEX *t, INTEGER *ldt, COMPLEX *c__, INTEGER *ldc, COMPLEX *work, 
	INTEGER *ldwork);
 
 int clarfg_(INTEGER *n, COMPLEX *alpha, COMPLEX *x, INTEGER *
	incx, COMPLEX *tau);
 
 int clarft_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, COMPLEX *v, INTEGER *ldv, COMPLEX *tau, COMPLEX *t, INTEGER *ldt);
 
 int clarfx_(char *side, INTEGER *m, INTEGER *n, COMPLEX *v, 
	COMPLEX *tau, COMPLEX *c__, INTEGER *ldc, COMPLEX *work);
 
 int clargv_(INTEGER *n, COMPLEX *x, INTEGER *incx, COMPLEX *
	y, INTEGER *incy, SINGLE *c__, INTEGER *incc);
 
 int clarnv_(INTEGER *idist, INTEGER *iseed, INTEGER *n, 
	COMPLEX *x);
 
 int clarrv_(INTEGER *n, SINGLE *d__, SINGLE *l, INTEGER *isplit, 
	INTEGER *m, SINGLE *w, INTEGER *iblock, SINGLE *gersch, SINGLE *tol, 
	COMPLEX *z__, INTEGER *ldz, INTEGER *isuppz, SINGLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int clartg_(COMPLEX *f, COMPLEX *g, SINGLE *cs, COMPLEX *sn, 
	COMPLEX *r__);
 
 int clartv_(INTEGER *n, COMPLEX *x, INTEGER *incx, COMPLEX *
	y, INTEGER *incy, SINGLE *c__, COMPLEX *s, INTEGER *incc);
 
 int clarz_(char *side, INTEGER *m, INTEGER *n, INTEGER *l, 
	COMPLEX *v, INTEGER *incv, COMPLEX *tau, COMPLEX *c__, INTEGER *ldc, 
	COMPLEX *work);
 
 int clarzb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, INTEGER *l, COMPLEX *v, 
	INTEGER *ldv, COMPLEX *t, INTEGER *ldt, COMPLEX *c__, INTEGER *ldc, 
	COMPLEX *work, INTEGER *ldwork);
 
 int clarzt_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, COMPLEX *v, INTEGER *ldv, COMPLEX *tau, COMPLEX *t, INTEGER *ldt);
 
 int clascl_(char *type__, INTEGER *kl, INTEGER *ku, SINGLE *
	cfrom, SINGLE *cto, INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda, 
	INTEGER *info);
 
 int claset_(char *uplo, INTEGER *m, INTEGER *n, COMPLEX *
	alpha, COMPLEX *beta, COMPLEX *a, INTEGER *lda);
 
 int clasr_(char *side, char *pivot, char *direct, INTEGER *m,
	 INTEGER *n, SINGLE *c__, SINGLE *s, COMPLEX *a, INTEGER *lda);
 
 int classq_(INTEGER *n, COMPLEX *x, INTEGER *incx, SINGLE *
	scale, SINGLE *sumsq);
 
 int claswp_(INTEGER *n, COMPLEX *a, INTEGER *lda, INTEGER *
	k1, INTEGER *k2, INTEGER *ipiv, INTEGER *incx);
 
 int clasyf_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 COMPLEX *a, INTEGER *lda, INTEGER *ipiv, COMPLEX *w, INTEGER *ldw, 
	INTEGER *info);
 
 int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, INTEGER *kd, COMPLEX *ab, INTEGER *ldab, COMPLEX *
	x, SINGLE *scale, SINGLE *cnorm, INTEGER *info);
 
 int clatdf_(INTEGER *ijob, INTEGER *n, COMPLEX *z__, INTEGER 
	*ldz, COMPLEX *rhs, SINGLE *rdsum, SINGLE *rdscal, INTEGER *ipiv, INTEGER 
	*jpiv);
 
 int clatps_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, COMPLEX *ap, COMPLEX *x, SINGLE *scale, SINGLE *cnorm,
	 INTEGER *info);
 
 int clatrd_(char *uplo, INTEGER *n, INTEGER *nb, COMPLEX *a, 
	INTEGER *lda, SINGLE *e, COMPLEX *tau, COMPLEX *w, INTEGER *ldw);
 
 int clatrs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *x, SINGLE *scale,
	 SINGLE *cnorm, INTEGER *info);
 
 int clatrz_(INTEGER *m, INTEGER *n, INTEGER *l, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work);
 
 int clatzm_(char *side, INTEGER *m, INTEGER *n, COMPLEX *v, 
	INTEGER *incv, COMPLEX *tau, COMPLEX *c1, COMPLEX *c2, INTEGER *ldc, 
	COMPLEX *work);
 
 int clauu2_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *info);
 
 int clauum_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *info);
 
 int cpbcon_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, SINGLE *anorm, SINGLE *rcond, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int cpbequ_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, SINGLE *s, SINGLE *scond, SINGLE *amax, INTEGER *info);
 
 int cpbrfs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *afb, INTEGER *ldafb, 
	COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *
	berr, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cpbstf_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, INTEGER *info);
 
 int cpbsv_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int cpbsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *kd, 
	INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *afb, INTEGER *
	ldafb, char *equed, SINGLE *s, COMPLEX *b, INTEGER *ldb, COMPLEX *x, 
	INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, COMPLEX *work, 
	SINGLE *rwork, INTEGER *info);
 
 int cpbtf2_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, INTEGER *info);
 
 int cpbtrf_(char *uplo, INTEGER *n, INTEGER *kd, COMPLEX *ab,
	 INTEGER *ldab, INTEGER *info);
 
 int cpbtrs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int cpocon_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 SINGLE *anorm, SINGLE *rcond, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cpoequ_(INTEGER *n, COMPLEX *a, INTEGER *lda, SINGLE *s, 
	SINGLE *scond, SINGLE *amax, INTEGER *info);
 
 int cporfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, COMPLEX *b, INTEGER *ldb,
	 COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, 
	SINGLE *rwork, INTEGER *info);
 
 int cposv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *a,
	 INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cposvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, char *
	equed, SINGLE *s, COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, 
	SINGLE *rcond, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int cpotf2_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *info);
 
 int cpotrf_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *info);
 
 int cpotri_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *info);
 
 int cpotrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cppcon_(char *uplo, INTEGER *n, COMPLEX *ap, SINGLE *anorm,
	 SINGLE *rcond, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cppequ_(char *uplo, INTEGER *n, COMPLEX *ap, SINGLE *s, 
	SINGLE *scond, SINGLE *amax, INTEGER *info);
 
 int cpprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, COMPLEX *afp, COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, 
	SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cppsv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cppsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *ap, COMPLEX *afp, char *equed, SINGLE *s, COMPLEX *b, 
	INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE 
	*berr, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int cpptrf_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	info);
 
 int cpptri_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	info);
 
 int cpptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cptcon_(INTEGER *n, SINGLE *d__, COMPLEX *e, SINGLE *anorm, 
	SINGLE *rcond, SINGLE *rwork, INTEGER *info);
 
 int cptrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *d__,
	 COMPLEX *e, SINGLE *df, COMPLEX *ef, COMPLEX *b, INTEGER *ldb, COMPLEX 
	*x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int cptsv_(INTEGER *n, INTEGER *nrhs, SINGLE *d__, COMPLEX *e, 
	COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cptsvx_(char *fact, INTEGER *n, INTEGER *nrhs, SINGLE *d__,
	 COMPLEX *e, SINGLE *df, COMPLEX *ef, COMPLEX *b, INTEGER *ldb, COMPLEX 
	*x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, COMPLEX *work, 
	SINGLE *rwork, INTEGER *info);
 
 int cpttrf_(INTEGER *n, SINGLE *d__, COMPLEX *e, INTEGER *info);
 
 int cpttrs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *d__,
	 COMPLEX *e, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cptts2_(INTEGER *iuplo, INTEGER *n, INTEGER *nrhs, SINGLE *
	d__, COMPLEX *e, COMPLEX *b, INTEGER *ldb);
 
 int crot_(INTEGER *n, COMPLEX *cx, INTEGER *incx, COMPLEX *
	cy, INTEGER *incy, SINGLE *c__, COMPLEX *s);
 
 int cspcon_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, SINGLE *anorm, SINGLE *rcond, COMPLEX *work, INTEGER *info);
 
 int cspmv_(char *uplo, INTEGER *n, COMPLEX *alpha, COMPLEX *
	ap, COMPLEX *x, INTEGER *incx, COMPLEX *beta, COMPLEX *y, INTEGER *
	incy);
 
 int cspr_(char *uplo, INTEGER *n, COMPLEX *alpha, COMPLEX *x,
	 INTEGER *incx, COMPLEX *ap);
 
 int csprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, COMPLEX *afp, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *x,
	 INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int cspsv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int cspsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *ap, COMPLEX *afp, INTEGER *ipiv, COMPLEX *b, INTEGER *
	ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int csptrf_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, INTEGER *info);
 
 int csptri_(char *uplo, INTEGER *n, COMPLEX *ap, INTEGER *
	ipiv, COMPLEX *work, INTEGER *info);
 
 int csptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	ap, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int csrot_(INTEGER *n, COMPLEX *cx, INTEGER *incx, COMPLEX *
	cy, INTEGER *incy, SINGLE *c__, SINGLE *s);
 
 int csrscl_(INTEGER *n, SINGLE *sa, COMPLEX *sx, INTEGER *incx);
 
 int cstedc_(char *compz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	COMPLEX *z__, INTEGER *ldz, COMPLEX *work, INTEGER *lwork, SINGLE *
	rwork, INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *
	info);
 
 int cstein_(INTEGER *n, SINGLE *d__, SINGLE *e, INTEGER *m, SINGLE 
	*w, INTEGER *iblock, INTEGER *isplit, COMPLEX *z__, INTEGER *ldz, 
	SINGLE *work, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int csteqr_(char *compz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	COMPLEX *z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int csycon_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, COMPLEX *work, INTEGER *
	info);
 
 int csymv_(char *uplo, INTEGER *n, COMPLEX *alpha, COMPLEX *
	a, INTEGER *lda, COMPLEX *x, INTEGER *incx, COMPLEX *beta, COMPLEX *y,
	 INTEGER *incy);
 
 int csyr_(char *uplo, INTEGER *n, COMPLEX *alpha, COMPLEX *x,
	 INTEGER *incx, COMPLEX *a, INTEGER *lda);
 
 int csyrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *ipiv, COMPLEX *
	b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int csysv_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *a,
	 INTEGER *lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *work,
	 INTEGER *lwork, INTEGER *info);
 
 int csysvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *af, INTEGER *ldaf, INTEGER *
	ipiv, COMPLEX *b, INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *rcond,
	 SINGLE *ferr, SINGLE *berr, COMPLEX *work, INTEGER *lwork, SINGLE *rwork, 
	INTEGER *info);
 
 int csytf2_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, INTEGER *info);
 
 int csytrf_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int csytri_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 INTEGER *ipiv, COMPLEX *work, INTEGER *info);
 
 int csytrs_(char *uplo, INTEGER *n, INTEGER *nrhs, COMPLEX *
	a, INTEGER *lda, INTEGER *ipiv, COMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int ctbcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	INTEGER *kd, COMPLEX *ab, INTEGER *ldab, SINGLE *rcond, COMPLEX *work, 
	SINGLE *rwork, INTEGER *info);
 
 int ctbrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *b, 
	INTEGER *ldb, COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, 
	COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int ctbtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, COMPLEX *ab, INTEGER *ldab, COMPLEX *b, 
	INTEGER *ldb, INTEGER *info);
 
 int ctgevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	COMPLEX *vl, INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, INTEGER *mm, 
	INTEGER *m, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int ctgex2_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *q, 
	INTEGER *ldq, COMPLEX *z__, INTEGER *ldz, INTEGER *j1, INTEGER *info);
 
 int ctgexc_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *q, 
	INTEGER *ldq, COMPLEX *z__, INTEGER *ldz, INTEGER *ifst, INTEGER *
	ilst, INTEGER *info);
 
 int ctgsen_(INTEGER *ijob, LOGICAL *wantq, LOGICAL *wantz, 
	LOGICAL *select, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, 
	INTEGER *ldb, COMPLEX *alpha, COMPLEX *beta, COMPLEX *q, INTEGER *ldq,
	 COMPLEX *z__, INTEGER *ldz, INTEGER *m, SINGLE *pl, SINGLE *pr, SINGLE *
	dif, COMPLEX *work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, 
	INTEGER *info);
 
 int ctgsja_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, INTEGER *k, INTEGER *l, COMPLEX *a, INTEGER *
	lda, COMPLEX *b, INTEGER *ldb, SINGLE *tola, SINGLE *tolb, SINGLE *alpha, 
	SINGLE *beta, COMPLEX *u, INTEGER *ldu, COMPLEX *v, INTEGER *ldv, 
	COMPLEX *q, INTEGER *ldq, COMPLEX *work, INTEGER *ncycle, INTEGER *
	info);
 
 int ctgsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	COMPLEX *vl, INTEGER *ldvl, COMPLEX *vr, INTEGER *ldvr, SINGLE *s, SINGLE 
	*dif, INTEGER *mm, INTEGER *m, COMPLEX *work, INTEGER *lwork, INTEGER 
	*iwork, INTEGER *info);
 
 int ctgsy2_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *d__, INTEGER *ldd, COMPLEX *e, INTEGER *lde, 
	COMPLEX *f, INTEGER *ldf, SINGLE *scale, SINGLE *rdsum, SINGLE *rdscal, 
	INTEGER *info);
 
 int ctgsyl_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *d__, INTEGER *ldd, COMPLEX *e, INTEGER *lde, 
	COMPLEX *f, INTEGER *ldf, SINGLE *scale, SINGLE *dif, COMPLEX *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int ctpcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	COMPLEX *ap, SINGLE *rcond, COMPLEX *work, SINGLE *rwork, INTEGER *info);
 
 int ctprfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, COMPLEX *ap, COMPLEX *b, INTEGER *ldb, COMPLEX *x, 
	INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int ctptri_(char *uplo, char *diag, INTEGER *n, COMPLEX *ap, 
	INTEGER *info);
 
 int ctptrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, COMPLEX *ap, COMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int ctrcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	COMPLEX *a, INTEGER *lda, SINGLE *rcond, COMPLEX *work, SINGLE *rwork, 
	INTEGER *info);
 
 int ctrevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, COMPLEX *t, INTEGER *ldt, COMPLEX *vl, INTEGER *ldvl, 
	COMPLEX *vr, INTEGER *ldvr, INTEGER *mm, INTEGER *m, COMPLEX *work, 
	SINGLE *rwork, INTEGER *info);
 
 int ctrexc_(char *compq, INTEGER *n, COMPLEX *t, INTEGER *
	ldt, COMPLEX *q, INTEGER *ldq, INTEGER *ifst, INTEGER *ilst, INTEGER *
	info);
 
 int ctrrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	COMPLEX *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, COMPLEX *work, SINGLE 
	*rwork, INTEGER *info);
 
 int ctrsen_(char *job, char *compq, LOGICAL *select, INTEGER 
	*n, COMPLEX *t, INTEGER *ldt, COMPLEX *q, INTEGER *ldq, COMPLEX *w, 
	INTEGER *m, SINGLE *s, SINGLE *sep, COMPLEX *work, INTEGER *lwork, 
	INTEGER *info);
 
 int ctrsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, COMPLEX *t, INTEGER *ldt, COMPLEX *vl, INTEGER *ldvl, 
	COMPLEX *vr, INTEGER *ldvr, SINGLE *s, SINGLE *sep, INTEGER *mm, INTEGER *
	m, COMPLEX *work, INTEGER *ldwork, SINGLE *rwork, INTEGER *info);
 
 int ctrsyl_(char *trana, char *tranb, INTEGER *isgn, INTEGER 
	*m, INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	COMPLEX *c__, INTEGER *ldc, SINGLE *scale, INTEGER *info);
 
 int ctrti2_(char *uplo, char *diag, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int ctrtri_(char *uplo, char *diag, INTEGER *n, COMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int ctrtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, COMPLEX *a, INTEGER *lda, COMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int ctzrqf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, INTEGER *info);
 
 int ctzrzf_(INTEGER *m, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cung2l_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cung2r_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cungbr_(char *vect, INTEGER *m, INTEGER *n, INTEGER *k, 
	COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int cunghr_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, COMPLEX *
	a, INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER 
	*info);
 
 int cungl2_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cunglq_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cungql_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cungqr_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cungr2_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *info);
 
 int cungrq_(INTEGER *m, INTEGER *n, INTEGER *k, COMPLEX *a, 
	INTEGER *lda, COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cungtr_(char *uplo, INTEGER *n, COMPLEX *a, INTEGER *lda,
	 COMPLEX *tau, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cunm2l_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *info);
 
 int cunm2r_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *info);
 
 int cunmbr_(char *vect, char *side, char *trans, INTEGER *m, 
	INTEGER *n, INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, 
	COMPLEX *c__, INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cunmhr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, COMPLEX *a, INTEGER *lda, COMPLEX *tau, 
	COMPLEX *c__, INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cunml2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *info);
 
 int cunmlq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cunmql_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cunmqr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cunmr2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *info);
 
 int cunmr3_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, COMPLEX *a, INTEGER *lda, COMPLEX *tau, 
	COMPLEX *c__, INTEGER *ldc, COMPLEX *work, INTEGER *info);
 
 int cunmrq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cunmrz_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, COMPLEX *a, INTEGER *lda, COMPLEX *tau, 
	COMPLEX *c__, INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int cunmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, COMPLEX *a, INTEGER *lda, COMPLEX *tau, COMPLEX *c__, 
	INTEGER *ldc, COMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int cupgtr_(char *uplo, INTEGER *n, COMPLEX *ap, COMPLEX *
	tau, COMPLEX *q, INTEGER *ldq, COMPLEX *work, INTEGER *info);
 
 int cupmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, COMPLEX *ap, COMPLEX *tau, COMPLEX *c__, INTEGER *ldc, 
	COMPLEX *work, INTEGER *info);
 
 int dbdsdc_(char *uplo, char *compq, INTEGER *n, DOUBLE *
	d__, DOUBLE *e, DOUBLE *u, INTEGER *ldu, DOUBLE *vt, 
	INTEGER *ldvt, DOUBLE *q, INTEGER *iq, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dbdsqr_(char *uplo, INTEGER *n, INTEGER *ncvt, INTEGER *
	nru, INTEGER *ncc, DOUBLE *d__, DOUBLE *e, DOUBLE *vt, 
	INTEGER *ldvt, DOUBLE *u, INTEGER *ldu, DOUBLE *c__, INTEGER *
	ldc, DOUBLE *work, INTEGER *info);
 
 int ddisna_(char *job, INTEGER *m, INTEGER *n, DOUBLE *
	d__, DOUBLE *sep, INTEGER *info);
 
 int dgbbrd_(char *vect, INTEGER *m, INTEGER *n, INTEGER *ncc,
	 INTEGER *kl, INTEGER *ku, DOUBLE *ab, INTEGER *ldab, DOUBLE *
	d__, DOUBLE *e, DOUBLE *q, INTEGER *ldq, DOUBLE *pt, 
	INTEGER *ldpt, DOUBLE *c__, INTEGER *ldc, DOUBLE *work, 
	INTEGER *info);
 
 int dgbcon_(char *norm, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLE *ab, INTEGER *ldab, INTEGER *ipiv, DOUBLE *anorm, 
	DOUBLE *rcond, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dgbequ_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLE *ab, INTEGER *ldab, DOUBLE *r__, DOUBLE *c__, 
	DOUBLE *rowcnd, DOUBLE *colcnd, DOUBLE *amax, INTEGER *
	info);
 
 int dgbrfs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE *afb, 
	INTEGER *ldafb, INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dgbsv_(INTEGER *n, INTEGER *kl, INTEGER *ku, INTEGER *
	nrhs, DOUBLE *ab, INTEGER *ldab, INTEGER *ipiv, DOUBLE *b, 
	INTEGER *ldb, INTEGER *info);
 
 int dgbsvx_(char *fact, char *trans, INTEGER *n, INTEGER *kl,
	 INTEGER *ku, INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, 
	DOUBLE *afb, INTEGER *ldafb, INTEGER *ipiv, char *equed, 
	DOUBLE *r__, DOUBLE *c__, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dgbtf2_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLE *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int dgbtrf_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLE *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int dgbtrs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, INTEGER *ipiv, 
	DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dgebak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, DOUBLE *scale, INTEGER *m, DOUBLE *v, INTEGER *
	ldv, INTEGER *info);
 
 int dgebal_(char *job, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ilo, INTEGER *ihi, DOUBLE *scale, INTEGER *info);
 
 int dgebd2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *d__, DOUBLE *e, DOUBLE *tauq, DOUBLE *
	taup, DOUBLE *work, INTEGER *info);
 
 int dgebrd_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *d__, DOUBLE *e, DOUBLE *tauq, DOUBLE *
	taup, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgecon_(char *norm, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dgeequ_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *r__, DOUBLE *c__, DOUBLE *rowcnd, DOUBLE 
	*colcnd, DOUBLE *amax, INTEGER *info);
 
 int dgees_(char *jobvs, char *sort, __CLPK_L_fp select, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, INTEGER *sdim, DOUBLE *wr, 
	DOUBLE *wi, DOUBLE *vs, INTEGER *ldvs, DOUBLE *work, 
	INTEGER *lwork, LOGICAL *bwork, INTEGER *info);
 
 int dgeesx_(char *jobvs, char *sort, __CLPK_L_fp select, char *
	sense, INTEGER *n, DOUBLE *a, INTEGER *lda, INTEGER *sdim, 
	DOUBLE *wr, DOUBLE *wi, DOUBLE *vs, INTEGER *ldvs, 
	DOUBLE *rconde, DOUBLE *rcondv, DOUBLE *work, INTEGER *
	lwork, INTEGER *iwork, INTEGER *liwork, LOGICAL *bwork, INTEGER *info);
 
 int dgeev_(char *jobvl, char *jobvr, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *wr, DOUBLE *wi, DOUBLE *vl, 
	INTEGER *ldvl, DOUBLE *vr, INTEGER *ldvr, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *wr, 
	DOUBLE *wi, DOUBLE *vl, INTEGER *ldvl, DOUBLE *vr, 
	INTEGER *ldvr, INTEGER *ilo, INTEGER *ihi, DOUBLE *scale, 
	DOUBLE *abnrm, DOUBLE *rconde, DOUBLE *rcondv, DOUBLE 
	*work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int dgegs_(char *jobvsl, char *jobvsr, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *
	alphar, DOUBLE *alphai, DOUBLE *beta, DOUBLE *vsl, 
	INTEGER *ldvsl, DOUBLE *vsr, INTEGER *ldvsr, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dgegv_(char *jobvl, char *jobvr, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *alphar, 
	DOUBLE *alphai, DOUBLE *beta, DOUBLE *vl, INTEGER *ldvl, 
	DOUBLE *vr, INTEGER *ldvr, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dgehd2_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, 
	INTEGER *info);
 
 int dgehrd_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dgelq2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dgelqf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgels_(char *trans, INTEGER *m, INTEGER *n, INTEGER *
	nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgelsd_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *
	s, DOUBLE *rcond, INTEGER *rank, DOUBLE *work, INTEGER *lwork,
	 INTEGER *iwork, INTEGER *info);
 
 int dgelss_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *
	s, DOUBLE *rcond, INTEGER *rank, DOUBLE *work, INTEGER *lwork,
	 INTEGER *info);
 
 int dgelsx_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *
	jpvt, DOUBLE *rcond, INTEGER *rank, DOUBLE *work, INTEGER *
	info);
 
 int dgelsy_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *
	jpvt, DOUBLE *rcond, INTEGER *rank, DOUBLE *work, INTEGER *
	lwork, INTEGER *info);
 
 int dgeql2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dgeqlf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgeqp3_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *jpvt, DOUBLE *tau, DOUBLE *work, INTEGER *lwork,
	 INTEGER *info);
 
 int dgeqpf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *jpvt, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dgeqr2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dgeqrf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgerfs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, INTEGER *
	ipiv, DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dgerq2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dgerqf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgesc2_(INTEGER *n, DOUBLE *a, INTEGER *lda, 
	DOUBLE *rhs, INTEGER *ipiv, INTEGER *jpiv, DOUBLE *scale);
 
 int dgesdd_(char *jobz, INTEGER *m, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *s, DOUBLE *u, INTEGER *ldu, 
	DOUBLE *vt, INTEGER *ldvt, DOUBLE *work, INTEGER *lwork, 
	INTEGER *iwork, INTEGER *info);
 
 int dgesv_(INTEGER *n, INTEGER *nrhs, DOUBLE *a, INTEGER 
	*lda, INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dgesvd_(char *jobu, char *jobvt, INTEGER *m, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *s, DOUBLE *u, INTEGER *
	ldu, DOUBLE *vt, INTEGER *ldvt, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dgesvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, 
	INTEGER *ipiv, char *equed, DOUBLE *r__, DOUBLE *c__, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dgetc2_(INTEGER *n, DOUBLE *a, INTEGER *lda, INTEGER 
	*ipiv, INTEGER *jpiv, INTEGER *info);
 
 int dgetf2_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, INTEGER *info);
 
 int dgetrf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, INTEGER *info);
 
 int dgetri_(INTEGER *n, DOUBLE *a, INTEGER *lda, INTEGER 
	*ipiv, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dgetrs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, INTEGER *ipiv, DOUBLE *b, INTEGER *
	ldb, INTEGER *info);
 
 int dggbak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, DOUBLE *lscale, DOUBLE *rscale, INTEGER *m, 
	DOUBLE *v, INTEGER *ldv, INTEGER *info);
 
 int dggbal_(char *job, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *b, INTEGER *ldb, INTEGER *ilo, INTEGER *ihi, 
	DOUBLE *lscale, DOUBLE *rscale, DOUBLE *work, INTEGER *
	info);
 
 int dgges_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	delctg, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, INTEGER *sdim, DOUBLE *alphar, DOUBLE *alphai, 
	DOUBLE *beta, DOUBLE *vsl, INTEGER *ldvsl, DOUBLE *vsr, 
	INTEGER *ldvsr, DOUBLE *work, INTEGER *lwork, LOGICAL *bwork, 
	INTEGER *info);
 
 int dggesx_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	delctg, char *sense, INTEGER *n, DOUBLE *a, INTEGER *lda, 
	DOUBLE *b, INTEGER *ldb, INTEGER *sdim, DOUBLE *alphar, 
	DOUBLE *alphai, DOUBLE *beta, DOUBLE *vsl, INTEGER *ldvsl,
	 DOUBLE *vsr, INTEGER *ldvsr, DOUBLE *rconde, DOUBLE *
	rcondv, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *
	liwork, LOGICAL *bwork, INTEGER *info);
 
 int dggev_(char *jobvl, char *jobvr, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *alphar, 
	DOUBLE *alphai, DOUBLE *beta, DOUBLE *vl, INTEGER *ldvl, 
	DOUBLE *vr, INTEGER *ldvr, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *alphar, DOUBLE *alphai, DOUBLE *
	beta, DOUBLE *vl, INTEGER *ldvl, DOUBLE *vr, INTEGER *ldvr, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *lscale, DOUBLE *rscale, 
	DOUBLE *abnrm, DOUBLE *bbnrm, DOUBLE *rconde, DOUBLE *
	rcondv, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, LOGICAL *
	bwork, INTEGER *info);
 
 int dggglm_(INTEGER *n, INTEGER *m, INTEGER *p, DOUBLE *
	a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *d__, 
	DOUBLE *x, DOUBLE *y, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dgghrd_(char *compq, char *compz, INTEGER *n, INTEGER *
	ilo, INTEGER *ihi, DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *q, INTEGER *ldq, DOUBLE *z__, INTEGER *
	ldz, INTEGER *info);
 
 int dgglse_(INTEGER *m, INTEGER *n, INTEGER *p, DOUBLE *
	a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *c__, 
	DOUBLE *d__, DOUBLE *x, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dggqrf_(INTEGER *n, INTEGER *m, INTEGER *p, DOUBLE *
	a, INTEGER *lda, DOUBLE *taua, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *taub, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dggrqf_(INTEGER *m, INTEGER *p, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *taua, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *taub, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dggsvd_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *n, INTEGER *p, INTEGER *k, INTEGER *l, DOUBLE *a, 
	INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *alpha, 
	DOUBLE *beta, DOUBLE *u, INTEGER *ldu, DOUBLE *v, INTEGER 
	*ldv, DOUBLE *q, INTEGER *ldq, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dggsvp_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *tola, DOUBLE *tolb, INTEGER *k, INTEGER 
	*l, DOUBLE *u, INTEGER *ldu, DOUBLE *v, INTEGER *ldv, 
	DOUBLE *q, INTEGER *ldq, INTEGER *iwork, DOUBLE *tau, 
	DOUBLE *work, INTEGER *info);
 
 int dgtcon_(char *norm, INTEGER *n, DOUBLE *dl, 
	DOUBLE *d__, DOUBLE *du, DOUBLE *du2, INTEGER *ipiv, 
	DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dgtrfs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *dl, DOUBLE *d__, DOUBLE *du, DOUBLE *dlf, 
	DOUBLE *df, DOUBLE *duf, DOUBLE *du2, INTEGER *ipiv, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *
	info);
 
 int dgtsv_(INTEGER *n, INTEGER *nrhs, DOUBLE *dl, 
	DOUBLE *d__, DOUBLE *du, DOUBLE *b, INTEGER *ldb, INTEGER 
	*info);
 
 int dgtsvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, DOUBLE *dl, DOUBLE *d__, DOUBLE *du, DOUBLE *
	dlf, DOUBLE *df, DOUBLE *duf, DOUBLE *du2, INTEGER *ipiv, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dgttrf_(INTEGER *n, DOUBLE *dl, DOUBLE *d__, 
	DOUBLE *du, DOUBLE *du2, INTEGER *ipiv, INTEGER *info);
 
 int dgttrs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *dl, DOUBLE *d__, DOUBLE *du, DOUBLE *du2, 
	INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dgtts2_(INTEGER *itrans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *dl, DOUBLE *d__, DOUBLE *du, DOUBLE *du2, 
	INTEGER *ipiv, DOUBLE *b, INTEGER *ldb);
 
 int dhgeqz_(char *job, char *compq, char *compz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *a, INTEGER *lda, DOUBLE *
	b, INTEGER *ldb, DOUBLE *alphar, DOUBLE *alphai, DOUBLE *
	beta, DOUBLE *q, INTEGER *ldq, DOUBLE *z__, INTEGER *ldz, 
	DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dhsein_(char *side, char *eigsrc, char *initv, LOGICAL *
	select, INTEGER *n, DOUBLE *h__, INTEGER *ldh, DOUBLE *wr, 
	DOUBLE *wi, DOUBLE *vl, INTEGER *ldvl, DOUBLE *vr, 
	INTEGER *ldvr, INTEGER *mm, INTEGER *m, DOUBLE *work, INTEGER *
	ifaill, INTEGER *ifailr, INTEGER *info);
 
 int dhseqr_(char *job, char *compz, INTEGER *n, INTEGER *ilo,
	 INTEGER *ihi, DOUBLE *h__, INTEGER *ldh, DOUBLE *wr, 
	DOUBLE *wi, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dlabad_(DOUBLE *small, DOUBLE *large);
 
 int dlabrd_(INTEGER *m, INTEGER *n, INTEGER *nb, DOUBLE *
	a, INTEGER *lda, DOUBLE *d__, DOUBLE *e, DOUBLE *tauq, 
	DOUBLE *taup, DOUBLE *x, INTEGER *ldx, DOUBLE *y, INTEGER 
	*ldy);
 
 int dlacon_(INTEGER *n, DOUBLE *v, DOUBLE *x, 
	INTEGER *isgn, DOUBLE *est, INTEGER *kase);
 
 int dlacpy_(char *uplo, INTEGER *m, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *b, INTEGER *ldb);
 
 int dladiv_(DOUBLE *a, DOUBLE *b, DOUBLE *c__, 
	DOUBLE *d__, DOUBLE *p, DOUBLE *q);
 
 int dlae2_(DOUBLE *a, DOUBLE *b, DOUBLE *c__, 
	DOUBLE *rt1, DOUBLE *rt2);
 
 int dlaebz_(INTEGER *ijob, INTEGER *nitmax, INTEGER *n, 
	INTEGER *mmax, INTEGER *minp, INTEGER *nbmin, DOUBLE *abstol, 
	DOUBLE *reltol, DOUBLE *pivmin, DOUBLE *d__, DOUBLE *
	e, DOUBLE *e2, INTEGER *nval, DOUBLE *ab, DOUBLE *c__, 
	INTEGER *mout, INTEGER *nab, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dlaed0_(INTEGER *icompq, INTEGER *qsiz, INTEGER *n, 
	DOUBLE *d__, DOUBLE *e, DOUBLE *q, INTEGER *ldq, 
	DOUBLE *qstore, INTEGER *ldqs, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dlaed1_(INTEGER *n, DOUBLE *d__, DOUBLE *q, 
	INTEGER *ldq, INTEGER *indxq, DOUBLE *rho, INTEGER *cutpnt, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dlaed2_(INTEGER *k, INTEGER *n, INTEGER *n1, DOUBLE *
	d__, DOUBLE *q, INTEGER *ldq, INTEGER *indxq, DOUBLE *rho, 
	DOUBLE *z__, DOUBLE *dlamda, DOUBLE *w, DOUBLE *q2, 
	INTEGER *indx, INTEGER *indxc, INTEGER *indxp, INTEGER *coltyp, 
	INTEGER *info);
 
 int dlaed3_(INTEGER *k, INTEGER *n, INTEGER *n1, DOUBLE *
	d__, DOUBLE *q, INTEGER *ldq, DOUBLE *rho, DOUBLE *dlamda,
	 DOUBLE *q2, INTEGER *indx, INTEGER *ctot, DOUBLE *w, 
	DOUBLE *s, INTEGER *info);
 
 int dlaed4_(INTEGER *n, INTEGER *i__, DOUBLE *d__, 
	DOUBLE *z__, DOUBLE *delta, DOUBLE *rho, DOUBLE *dlam,
	 INTEGER *info);
 
 int dlaed5_(INTEGER *i__, DOUBLE *d__, DOUBLE *z__, 
	DOUBLE *delta, DOUBLE *rho, DOUBLE *dlam);
 
 int dlaed6_(INTEGER *kniter, LOGICAL *orgati, DOUBLE *
	rho, DOUBLE *d__, DOUBLE *z__, DOUBLE *finit, DOUBLE *
	tau, INTEGER *info);
 
 int dlaed7_(INTEGER *icompq, INTEGER *n, INTEGER *qsiz, 
	INTEGER *tlvls, INTEGER *curlvl, INTEGER *curpbm, DOUBLE *d__, 
	DOUBLE *q, INTEGER *ldq, INTEGER *indxq, DOUBLE *rho, INTEGER 
	*cutpnt, DOUBLE *qstore, INTEGER *qptr, INTEGER *prmptr, INTEGER *
	perm, INTEGER *givptr, INTEGER *givcol, DOUBLE *givnum, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dlaed8_(INTEGER *icompq, INTEGER *k, INTEGER *n, INTEGER 
	*qsiz, DOUBLE *d__, DOUBLE *q, INTEGER *ldq, INTEGER *indxq, 
	DOUBLE *rho, INTEGER *cutpnt, DOUBLE *z__, DOUBLE *dlamda,
	 DOUBLE *q2, INTEGER *ldq2, DOUBLE *w, INTEGER *perm, INTEGER 
	*givptr, INTEGER *givcol, DOUBLE *givnum, INTEGER *indxp, INTEGER 
	*indx, INTEGER *info);
 
 int dlaed9_(INTEGER *k, INTEGER *kstart, INTEGER *kstop, 
	INTEGER *n, DOUBLE *d__, DOUBLE *q, INTEGER *ldq, DOUBLE *
	rho, DOUBLE *dlamda, DOUBLE *w, DOUBLE *s, INTEGER *lds, 
	INTEGER *info);
 
 int dlaeda_(INTEGER *n, INTEGER *tlvls, INTEGER *curlvl, 
	INTEGER *curpbm, INTEGER *prmptr, INTEGER *perm, INTEGER *givptr, 
	INTEGER *givcol, DOUBLE *givnum, DOUBLE *q, INTEGER *qptr, 
	DOUBLE *z__, DOUBLE *ztemp, INTEGER *info);
 
 int dlaein_(LOGICAL *rightv, LOGICAL *noinit, INTEGER *n, 
	DOUBLE *h__, INTEGER *ldh, DOUBLE *wr, DOUBLE *wi, 
	DOUBLE *vr, DOUBLE *vi, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *work, DOUBLE *eps3, DOUBLE *smlnum, DOUBLE *
	bignum, INTEGER *info);
 
 int dlaev2_(DOUBLE *a, DOUBLE *b, DOUBLE *c__, 
	DOUBLE *rt1, DOUBLE *rt2, DOUBLE *cs1, DOUBLE *sn1);
 
 int dlaexc_(LOGICAL *wantq, INTEGER *n, DOUBLE *t, 
	INTEGER *ldt, DOUBLE *q, INTEGER *ldq, INTEGER *j1, INTEGER *n1, 
	INTEGER *n2, DOUBLE *work, INTEGER *info);
 
 int dlag2_(DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *safmin, DOUBLE *scale1, DOUBLE *
	scale2, DOUBLE *wr1, DOUBLE *wr2, DOUBLE *wi);
 
 int dlags2_(LOGICAL *upper, DOUBLE *a1, DOUBLE *a2, 
	DOUBLE *a3, DOUBLE *b1, DOUBLE *b2, DOUBLE *b3, 
	DOUBLE *csu, DOUBLE *snu, DOUBLE *csv, DOUBLE *snv, 
	DOUBLE *csq, DOUBLE *snq);
 
 int dlagtf_(INTEGER *n, DOUBLE *a, DOUBLE *lambda, 
	DOUBLE *b, DOUBLE *c__, DOUBLE *tol, DOUBLE *d__, 
	INTEGER *in, INTEGER *info);
 
 int dlagtm_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *alpha, DOUBLE *dl, DOUBLE *d__, DOUBLE *du, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *beta, DOUBLE *b, INTEGER 
	*ldb);
 
 int dlagts_(INTEGER *job, INTEGER *n, DOUBLE *a, 
	DOUBLE *b, DOUBLE *c__, DOUBLE *d__, INTEGER *in, 
	DOUBLE *y, DOUBLE *tol, INTEGER *info);
 
 int dlagv2_(DOUBLE *a, INTEGER *lda, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *alphar, DOUBLE *alphai, DOUBLE *
	beta, DOUBLE *csl, DOUBLE *snl, DOUBLE *csr, DOUBLE *
	snr);
 
 int dlahqr_(LOGICAL *wantt, LOGICAL *wantz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *h__, INTEGER *ldh, DOUBLE 
	*wr, DOUBLE *wi, INTEGER *iloz, INTEGER *ihiz, DOUBLE *z__, 
	INTEGER *ldz, INTEGER *info);
 
 int dlahrd_(INTEGER *n, INTEGER *k, INTEGER *nb, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *t, INTEGER *ldt, 
	DOUBLE *y, INTEGER *ldy);
 
 int dlaic1_(INTEGER *job, INTEGER *j, DOUBLE *x, 
	DOUBLE *sest, DOUBLE *w, DOUBLE *gamma, DOUBLE *
	sestpr, DOUBLE *s, DOUBLE *c__);
 
 int dlaln2_(LOGICAL *ltrans, INTEGER *na, INTEGER *nw, 
	DOUBLE *smin, DOUBLE *ca, DOUBLE *a, INTEGER *lda, 
	DOUBLE *d1, DOUBLE *d2, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *wr, DOUBLE *wi, DOUBLE *x, INTEGER *ldx, 
	DOUBLE *scale, DOUBLE *xnorm, INTEGER *info);
 
 int dlals0_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *nrhs, DOUBLE *b, INTEGER *ldb, DOUBLE 
	*bx, INTEGER *ldbx, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, 
	INTEGER *ldgcol, DOUBLE *givnum, INTEGER *ldgnum, DOUBLE *
	poles, DOUBLE *difl, DOUBLE *difr, DOUBLE *z__, INTEGER *
	k, DOUBLE *c__, DOUBLE *s, DOUBLE *work, INTEGER *info);
 
 int dlalsa_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *nrhs, DOUBLE *b, INTEGER *ldb, DOUBLE *bx, INTEGER *
	ldbx, DOUBLE *u, INTEGER *ldu, DOUBLE *vt, INTEGER *k, 
	DOUBLE *difl, DOUBLE *difr, DOUBLE *z__, DOUBLE *
	poles, INTEGER *givptr, INTEGER *givcol, INTEGER *ldgcol, INTEGER *
	perm, DOUBLE *givnum, DOUBLE *c__, DOUBLE *s, DOUBLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int dlalsd_(char *uplo, INTEGER *smlsiz, INTEGER *n, INTEGER 
	*nrhs, DOUBLE *d__, DOUBLE *e, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *rcond, INTEGER *rank, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dlamc1_(INTEGER *beta, INTEGER *t, LOGICAL *rnd, LOGICAL 
	*ieee1);
 
 int dlamc2_(INTEGER *beta, INTEGER *t, LOGICAL *rnd, 
	DOUBLE *eps, INTEGER *emin, DOUBLE *rmin, INTEGER *emax, 
	DOUBLE *rmax);
 
 int dlamc4_(INTEGER *emin, DOUBLE *start, INTEGER *base);
 
 int dlamc5_(INTEGER *beta, INTEGER *p, INTEGER *emin, 
	LOGICAL *ieee, INTEGER *emax, DOUBLE *rmax);
 
 int dlamrg_(INTEGER *n1, INTEGER *n2, DOUBLE *a, INTEGER 
	*dtrd1, INTEGER *dtrd2, INTEGER *index);
 
 int dlanv2_(DOUBLE *a, DOUBLE *b, DOUBLE *c__, 
	DOUBLE *d__, DOUBLE *rt1r, DOUBLE *rt1i, DOUBLE *rt2r,
	 DOUBLE *rt2i, DOUBLE *cs, DOUBLE *sn);
 
 int dlapll_(INTEGER *n, DOUBLE *x, INTEGER *incx, 
	DOUBLE *y, INTEGER *incy, DOUBLE *ssmin);
 
 int dlapmt_(LOGICAL *forwrd, INTEGER *m, INTEGER *n, 
	DOUBLE *x, INTEGER *ldx, INTEGER *k);
 
 int dlaqgb_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLE *ab, INTEGER *ldab, DOUBLE *r__, DOUBLE *c__, 
	DOUBLE *rowcnd, DOUBLE *colcnd, DOUBLE *amax, char *equed);
 
 int dlaqge_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *r__, DOUBLE *c__, DOUBLE *rowcnd, DOUBLE 
	*colcnd, DOUBLE *amax, char *equed);
 
 int dlaqp2_(INTEGER *m, INTEGER *n, INTEGER *offset, 
	DOUBLE *a, INTEGER *lda, INTEGER *jpvt, DOUBLE *tau, 
	DOUBLE *vn1, DOUBLE *vn2, DOUBLE *work);
 
 int dlaqps_(INTEGER *m, INTEGER *n, INTEGER *offset, INTEGER 
	*nb, INTEGER *kb, DOUBLE *a, INTEGER *lda, INTEGER *jpvt, 
	DOUBLE *tau, DOUBLE *vn1, DOUBLE *vn2, DOUBLE *auxv, 
	DOUBLE *f, INTEGER *ldf);
 
 int dlaqsb_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, DOUBLE *s, DOUBLE *scond, DOUBLE *amax,
	 char *equed);
 
 int dlaqsp_(char *uplo, INTEGER *n, DOUBLE *ap, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, char *equed);
 
 int dlaqsy_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *s, DOUBLE *scond, DOUBLE *amax, char *equed);
 
 int dlaqtr_(LOGICAL *ltran, LOGICAL *lreal, INTEGER *n, 
	DOUBLE *t, INTEGER *ldt, DOUBLE *b, DOUBLE *w, DOUBLE 
	*scale, DOUBLE *x, DOUBLE *work, INTEGER *info);
 
 int dlar1v_(INTEGER *n, INTEGER *b1, INTEGER *bn, DOUBLE 
	*sigma, DOUBLE *d__, DOUBLE *l, DOUBLE *ld, DOUBLE *
	lld, DOUBLE *gersch, DOUBLE *z__, DOUBLE *ztz, DOUBLE 
	*mingma, INTEGER *r__, INTEGER *isuppz, DOUBLE *work);
 
 int dlar2v_(INTEGER *n, DOUBLE *x, DOUBLE *y, 
	DOUBLE *z__, INTEGER *incx, DOUBLE *c__, DOUBLE *s, 
	INTEGER *incc);
 
 int dlarf_(char *side, INTEGER *m, INTEGER *n, DOUBLE *v,
	 INTEGER *incv, DOUBLE *tau, DOUBLE *c__, INTEGER *ldc, 
	DOUBLE *work);
 
 int dlarfb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *v, INTEGER *
	ldv, DOUBLE *t, INTEGER *ldt, DOUBLE *c__, INTEGER *ldc, 
	DOUBLE *work, INTEGER *ldwork);
 
 int dlarfg_(INTEGER *n, DOUBLE *alpha, DOUBLE *x, 
	INTEGER *incx, DOUBLE *tau);
 
 int dlarft_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, DOUBLE *v, INTEGER *ldv, DOUBLE *tau, DOUBLE *t, 
	INTEGER *ldt);
 
 int dlarfx_(char *side, INTEGER *m, INTEGER *n, DOUBLE *
	v, DOUBLE *tau, DOUBLE *c__, INTEGER *ldc, DOUBLE *work);
 
 int dlargv_(INTEGER *n, DOUBLE *x, INTEGER *incx, 
	DOUBLE *y, INTEGER *incy, DOUBLE *c__, INTEGER *incc);
 
 int dlarnv_(INTEGER *idist, INTEGER *iseed, INTEGER *n, 
	DOUBLE *x);
 
 int dlarrb_(INTEGER *n, DOUBLE *d__, DOUBLE *l, 
	DOUBLE *ld, DOUBLE *lld, INTEGER *ifirst, INTEGER *ilast, 
	DOUBLE *sigma, DOUBLE *reltol, DOUBLE *w, DOUBLE *
	wgap, DOUBLE *werr, DOUBLE *work, INTEGER *iwork, INTEGER *
	info);
 
 int dlarre_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	DOUBLE *tol, INTEGER *nsplit, INTEGER *isplit, INTEGER *m, 
	DOUBLE *w, DOUBLE *woff, DOUBLE *gersch, DOUBLE *work,
	 INTEGER *info);
 
 int dlarrf_(INTEGER *n, DOUBLE *d__, DOUBLE *l, 
	DOUBLE *ld, DOUBLE *lld, INTEGER *ifirst, INTEGER *ilast, 
	DOUBLE *w, DOUBLE *dplus, DOUBLE *lplus, DOUBLE *work,
	 INTEGER *iwork, INTEGER *info);
 
 int dlarrv_(INTEGER *n, DOUBLE *d__, DOUBLE *l, 
	INTEGER *isplit, INTEGER *m, DOUBLE *w, INTEGER *iblock, 
	DOUBLE *gersch, DOUBLE *tol, DOUBLE *z__, INTEGER *ldz, 
	INTEGER *isuppz, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dlartg_(DOUBLE *f, DOUBLE *g, DOUBLE *cs, 
	DOUBLE *sn, DOUBLE *r__);
 
 int dlartv_(INTEGER *n, DOUBLE *x, INTEGER *incx, 
	DOUBLE *y, INTEGER *incy, DOUBLE *c__, DOUBLE *s, INTEGER 
	*incc);
 
 int dlaruv_(INTEGER *iseed, INTEGER *n, DOUBLE *x);
 
 int dlarz_(char *side, INTEGER *m, INTEGER *n, INTEGER *l, 
	DOUBLE *v, INTEGER *incv, DOUBLE *tau, DOUBLE *c__, 
	INTEGER *ldc, DOUBLE *work);
 
 int dlarzb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, INTEGER *l, DOUBLE *v,
	 INTEGER *ldv, DOUBLE *t, INTEGER *ldt, DOUBLE *c__, INTEGER *
	ldc, DOUBLE *work, INTEGER *ldwork);
 
 int dlarzt_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, DOUBLE *v, INTEGER *ldv, DOUBLE *tau, DOUBLE *t, 
	INTEGER *ldt);
 
 int dlas2_(DOUBLE *f, DOUBLE *g, DOUBLE *h__, 
	DOUBLE *ssmin, DOUBLE *ssmax);
 
 int dlascl_(char *type__, INTEGER *kl, INTEGER *ku, 
	DOUBLE *cfrom, DOUBLE *cto, INTEGER *m, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, INTEGER *info);
 
 int dlasd0_(INTEGER *n, INTEGER *sqre, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *u, INTEGER *ldu, DOUBLE *vt, INTEGER *
	ldvt, INTEGER *smlsiz, INTEGER *iwork, DOUBLE *work, INTEGER *
	info);
 
 int dlasd1_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, 
	DOUBLE *d__, DOUBLE *alpha, DOUBLE *beta, DOUBLE *u, 
	INTEGER *ldu, DOUBLE *vt, INTEGER *ldvt, INTEGER *idxq, INTEGER *
	iwork, DOUBLE *work, INTEGER *info);
 
 int dlasd2_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, INTEGER 
	*k, DOUBLE *d__, DOUBLE *z__, DOUBLE *alpha, DOUBLE *
	beta, DOUBLE *u, INTEGER *ldu, DOUBLE *vt, INTEGER *ldvt, 
	DOUBLE *dsigma, DOUBLE *u2, INTEGER *ldu2, DOUBLE *vt2, 
	INTEGER *ldvt2, INTEGER *idxp, INTEGER *idx, INTEGER *idxc, INTEGER *
	idxq, INTEGER *coltyp, INTEGER *info);
 
 int dlasd3_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, INTEGER 
	*k, DOUBLE *d__, DOUBLE *q, INTEGER *ldq, DOUBLE *dsigma, 
	DOUBLE *u, INTEGER *ldu, DOUBLE *u2, INTEGER *ldu2, 
	DOUBLE *vt, INTEGER *ldvt, DOUBLE *vt2, INTEGER *ldvt2, 
	INTEGER *idxc, INTEGER *ctot, DOUBLE *z__, INTEGER *info);
 
 int dlasd4_(INTEGER *n, INTEGER *i__, DOUBLE *d__, 
	DOUBLE *z__, DOUBLE *delta, DOUBLE *rho, DOUBLE *
	sigma, DOUBLE *work, INTEGER *info);
 
 int dlasd5_(INTEGER *i__, DOUBLE *d__, DOUBLE *z__, 
	DOUBLE *delta, DOUBLE *rho, DOUBLE *dsigma, DOUBLE *
	work);
 
 int dlasd6_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, DOUBLE *d__, DOUBLE *vf, DOUBLE *vl, 
	DOUBLE *alpha, DOUBLE *beta, INTEGER *idxq, INTEGER *perm, 
	INTEGER *givptr, INTEGER *givcol, INTEGER *ldgcol, DOUBLE *givnum,
	 INTEGER *ldgnum, DOUBLE *poles, DOUBLE *difl, DOUBLE *
	difr, DOUBLE *z__, INTEGER *k, DOUBLE *c__, DOUBLE *s, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dlasd7_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *k, DOUBLE *d__, DOUBLE *z__, 
	DOUBLE *zw, DOUBLE *vf, DOUBLE *vfw, DOUBLE *vl, 
	DOUBLE *vlw, DOUBLE *alpha, DOUBLE *beta, DOUBLE *
	dsigma, INTEGER *idx, INTEGER *idxp, INTEGER *idxq, INTEGER *perm, 
	INTEGER *givptr, INTEGER *givcol, INTEGER *ldgcol, DOUBLE *givnum,
	 INTEGER *ldgnum, DOUBLE *c__, DOUBLE *s, INTEGER *info);
 
 int dlasd8_(INTEGER *icompq, INTEGER *k, DOUBLE *d__, 
	DOUBLE *z__, DOUBLE *vf, DOUBLE *vl, DOUBLE *difl, 
	DOUBLE *difr, INTEGER *lddifr, DOUBLE *dsigma, DOUBLE *
	work, INTEGER *info);
 
 int dlasd9_(INTEGER *icompq, INTEGER *ldu, INTEGER *k, 
	DOUBLE *d__, DOUBLE *z__, DOUBLE *vf, DOUBLE *vl, 
	DOUBLE *difl, DOUBLE *difr, DOUBLE *dsigma, DOUBLE *
	work, INTEGER *info);
 
 int dlasda_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *sqre, DOUBLE *d__, DOUBLE *e, DOUBLE *u, INTEGER 
	*ldu, DOUBLE *vt, INTEGER *k, DOUBLE *difl, DOUBLE *difr, 
	DOUBLE *z__, DOUBLE *poles, INTEGER *givptr, INTEGER *givcol, 
	INTEGER *ldgcol, INTEGER *perm, DOUBLE *givnum, DOUBLE *c__, 
	DOUBLE *s, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dlasdq_(char *uplo, INTEGER *sqre, INTEGER *n, INTEGER *
	ncvt, INTEGER *nru, INTEGER *ncc, DOUBLE *d__, DOUBLE *e, 
	DOUBLE *vt, INTEGER *ldvt, DOUBLE *u, INTEGER *ldu, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dlasdt_(INTEGER *n, INTEGER *lvl, INTEGER *nd, INTEGER *
	inode, INTEGER *ndiml, INTEGER *ndimr, INTEGER *msub);
 
 int dlaset_(char *uplo, INTEGER *m, INTEGER *n, DOUBLE *
	alpha, DOUBLE *beta, DOUBLE *a, INTEGER *lda);
 
 int dlasq1_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	DOUBLE *work, INTEGER *info);
 
 int dlasq2_(INTEGER *n, DOUBLE *z__, INTEGER *info);
 
 int dlasq3_(INTEGER *i0, INTEGER *n0, DOUBLE *z__, 
	INTEGER *pp, DOUBLE *dmin__, DOUBLE *sigma, DOUBLE *desig,
	 DOUBLE *qmax, INTEGER *nfail, INTEGER *iter, INTEGER *ndiv, 
	LOGICAL *ieee);
 
 int dlasq4_(INTEGER *i0, INTEGER *n0, DOUBLE *z__, 
	INTEGER *pp, INTEGER *n0in, DOUBLE *dmin__, DOUBLE *dmin1, 
	DOUBLE *dmin2, DOUBLE *dn, DOUBLE *dn1, DOUBLE *dn2, 
	DOUBLE *tau, INTEGER *ttype);
 
 int dlasq5_(INTEGER *i0, INTEGER *n0, DOUBLE *z__, 
	INTEGER *pp, DOUBLE *tau, DOUBLE *dmin__, DOUBLE *dmin1, 
	DOUBLE *dmin2, DOUBLE *dn, DOUBLE *dnm1, DOUBLE *dnm2,
	 LOGICAL *ieee);
 
 int dlasq6_(INTEGER *i0, INTEGER *n0, DOUBLE *z__, 
	INTEGER *pp, DOUBLE *dmin__, DOUBLE *dmin1, DOUBLE *dmin2,
	 DOUBLE *dn, DOUBLE *dnm1, DOUBLE *dnm2);
 
 int dlasr_(char *side, char *pivot, char *direct, INTEGER *m,
	 INTEGER *n, DOUBLE *c__, DOUBLE *s, DOUBLE *a, INTEGER *
	lda);
 
 int dlasrt_(char *id, INTEGER *n, DOUBLE *d__, INTEGER *
	info);
 
 int dlassq_(INTEGER *n, DOUBLE *x, INTEGER *incx, 
	DOUBLE *scale, DOUBLE *sumsq);
 
 int dlasv2_(DOUBLE *f, DOUBLE *g, DOUBLE *h__, 
	DOUBLE *ssmin, DOUBLE *ssmax, DOUBLE *snr, DOUBLE *
	csr, DOUBLE *snl, DOUBLE *csl);
 
 int dlaswp_(INTEGER *n, DOUBLE *a, INTEGER *lda, INTEGER 
	*k1, INTEGER *k2, INTEGER *ipiv, INTEGER *incx);
 
 int dlasy2_(LOGICAL *ltranl, LOGICAL *ltranr, INTEGER *isgn, 
	INTEGER *n1, INTEGER *n2, DOUBLE *tl, INTEGER *ldtl, DOUBLE *
	tr, INTEGER *ldtr, DOUBLE *b, INTEGER *ldb, DOUBLE *scale, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *xnorm, INTEGER *info);
 
 int dlasyf_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 DOUBLE *a, INTEGER *lda, INTEGER *ipiv, DOUBLE *w, INTEGER *
	ldw, INTEGER *info);
 
 int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, INTEGER *kd, DOUBLE *ab, INTEGER *ldab, 
	DOUBLE *x, DOUBLE *scale, DOUBLE *cnorm, INTEGER *info);
 
 int dlatdf_(INTEGER *ijob, INTEGER *n, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *rhs, DOUBLE *rdsum, DOUBLE *rdscal, 
	INTEGER *ipiv, INTEGER *jpiv);
 
 int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, DOUBLE *ap, DOUBLE *x, DOUBLE *scale, 
	DOUBLE *cnorm, INTEGER *info);
 
 int dlatrd_(char *uplo, INTEGER *n, INTEGER *nb, DOUBLE *
	a, INTEGER *lda, DOUBLE *e, DOUBLE *tau, DOUBLE *w, 
	INTEGER *ldw);
 
 int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *x, 
	DOUBLE *scale, DOUBLE *cnorm, INTEGER *info);
 
 int dlatrz_(INTEGER *m, INTEGER *n, INTEGER *l, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work);
 
 int dlatzm_(char *side, INTEGER *m, INTEGER *n, DOUBLE *
	v, INTEGER *incv, DOUBLE *tau, DOUBLE *c1, DOUBLE *c2, 
	INTEGER *ldc, DOUBLE *work);
 
 int dlauu2_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *info);
 
 int dlauum_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *info);
 
 int dopgtr_(char *uplo, INTEGER *n, DOUBLE *ap, 
	DOUBLE *tau, DOUBLE *q, INTEGER *ldq, DOUBLE *work, 
	INTEGER *info);
 
 int dopmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, DOUBLE *ap, DOUBLE *tau, DOUBLE *c__, INTEGER 
	*ldc, DOUBLE *work, INTEGER *info);
 
 int dorg2l_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dorg2r_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dorgbr_(char *vect, INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dorghr_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int dorgl2_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dorglq_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dorgql_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dorgqr_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dorgr2_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *info);
 
 int dorgrq_(INTEGER *m, INTEGER *n, INTEGER *k, DOUBLE *
	a, INTEGER *lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dorgtr_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dorm2l_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dorm2r_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dormbr_(char *vect, char *side, char *trans, INTEGER *m, 
	INTEGER *n, INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dormhr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *a, INTEGER *lda, DOUBLE *
	tau, DOUBLE *c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dorml2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dormlq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dormql_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dormqr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dormr2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dormr3_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, DOUBLE *a, INTEGER *lda, DOUBLE *tau, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *work, INTEGER *info);
 
 int dormrq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dormrz_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, DOUBLE *a, INTEGER *lda, DOUBLE *tau, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dormtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *tau, DOUBLE *
	c__, INTEGER *ldc, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dpbcon_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, DOUBLE *anorm, DOUBLE *rcond, DOUBLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int dpbequ_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, DOUBLE *s, DOUBLE *scond, DOUBLE *amax,
	 INTEGER *info);
 
 int dpbrfs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE *afb, INTEGER *ldafb, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *
	info);
 
 int dpbstf_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, INTEGER *info);
 
 int dpbsv_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE *b, INTEGER *ldb, 
	INTEGER *info);
 
 int dpbsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *kd, 
	INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE *afb, 
	INTEGER *ldafb, char *equed, DOUBLE *s, DOUBLE *b, INTEGER *
	ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr,
	 DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dpbtf2_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, INTEGER *info);
 
 int dpbtrf_(char *uplo, INTEGER *n, INTEGER *kd, DOUBLE *
	ab, INTEGER *ldab, INTEGER *info);
 
 int dpbtrs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE *b, INTEGER *ldb, 
	INTEGER *info);
 
 int dpocon_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dpoequ_(INTEGER *n, DOUBLE *a, INTEGER *lda, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, INTEGER *info);
 
 int dporfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *
	info);
 
 int dposv_(char *uplo, INTEGER *n, INTEGER *nrhs, DOUBLE 
	*a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dposvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, 
	char *equed, DOUBLE *s, DOUBLE *b, INTEGER *ldb, DOUBLE *
	x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, DOUBLE *
	berr, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dpotf2_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *info);
 
 int dpotrf_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *info);
 
 int dpotri_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *info);
 
 int dpotrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int dppcon_(char *uplo, INTEGER *n, DOUBLE *ap, 
	DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dppequ_(char *uplo, INTEGER *n, DOUBLE *ap, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, INTEGER *info);
 
 int dpprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *ap, DOUBLE *afp, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dppsv_(char *uplo, INTEGER *n, INTEGER *nrhs, DOUBLE 
	*ap, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dppsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLE *ap, DOUBLE *afp, char *equed, DOUBLE *s, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int dpptrf_(char *uplo, INTEGER *n, DOUBLE *ap, INTEGER *
	info);
 
 int dpptri_(char *uplo, INTEGER *n, DOUBLE *ap, INTEGER *
	info);
 
 int dpptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *ap, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dptcon_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER *info);
 
 int dpteqr_(char *compz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int dptrfs_(INTEGER *n, INTEGER *nrhs, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *df, DOUBLE *ef, DOUBLE *b, INTEGER 
	*ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr,
	 DOUBLE *work, INTEGER *info);
 
 int dptsv_(INTEGER *n, INTEGER *nrhs, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dptsvx_(char *fact, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *d__, DOUBLE *e, DOUBLE *df, DOUBLE *ef, 
	DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *
	rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *
	info);
 
 int dpttrf_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	INTEGER *info);
 
 int dpttrs_(INTEGER *n, INTEGER *nrhs, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dptts2_(INTEGER *n, INTEGER *nrhs, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *b, INTEGER *ldb);
 
 int drscl_(INTEGER *n, DOUBLE *sa, DOUBLE *sx, 
	INTEGER *incx);
 
 int dsbev_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLE *ab, INTEGER *ldab, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *info);
 
 int dsbevd_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLE *ab, INTEGER *ldab, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int dsbevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *kd, DOUBLE *ab, INTEGER *ldab, DOUBLE *q, INTEGER *
	ldq, DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *iu, 
	DOUBLE *abstol, INTEGER *m, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *iwork, INTEGER *ifail, 
	INTEGER *info);
 
 int dsbgst_(char *vect, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, DOUBLE *ab, INTEGER *ldab, DOUBLE *bb, INTEGER *
	ldbb, DOUBLE *x, INTEGER *ldx, DOUBLE *work, INTEGER *info);
 
 int dsbgv_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, DOUBLE *ab, INTEGER *ldab, DOUBLE *bb, INTEGER *
	ldbb, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int dsbgvd_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, DOUBLE *ab, INTEGER *ldab, DOUBLE *bb, INTEGER *
	ldbb, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dsbgvx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *ka, INTEGER *kb, DOUBLE *ab, INTEGER *ldab, DOUBLE *
	bb, INTEGER *ldbb, DOUBLE *q, INTEGER *ldq, DOUBLE *vl, 
	DOUBLE *vu, INTEGER *il, INTEGER *iu, DOUBLE *abstol, INTEGER 
	*m, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int dsbtrd_(char *vect, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLE *ab, INTEGER *ldab, DOUBLE *d__, DOUBLE *e, 
	DOUBLE *q, INTEGER *ldq, DOUBLE *work, INTEGER *info);
 
 int dspcon_(char *uplo, INTEGER *n, DOUBLE *ap, INTEGER *
	ipiv, DOUBLE *anorm, DOUBLE *rcond, DOUBLE *work, INTEGER 
	*iwork, INTEGER *info);
 
 int dspev_(char *jobz, char *uplo, INTEGER *n, DOUBLE *
	ap, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int dspevd_(char *jobz, char *uplo, INTEGER *n, DOUBLE *
	ap, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dspevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLE *ap, DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *
	iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *iwork, INTEGER *ifail, 
	INTEGER *info);
 
 int dspgst_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLE *ap, DOUBLE *bp, INTEGER *info);
 
 int dspgv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLE *ap, DOUBLE *bp, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *info);
 
 int dspgvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLE *ap, DOUBLE *bp, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int dspgvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, DOUBLE *ap, DOUBLE *bp, DOUBLE *vl, 
	DOUBLE *vu, INTEGER *il, INTEGER *iu, DOUBLE *abstol, INTEGER 
	*m, DOUBLE *w, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int dsprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *ap, DOUBLE *afp, INTEGER *ipiv, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dspsv_(char *uplo, INTEGER *n, INTEGER *nrhs, DOUBLE 
	*ap, INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, INTEGER *info);
 
 int dspsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLE *ap, DOUBLE *afp, INTEGER *ipiv, DOUBLE *b, 
	INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *rcond, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dsptrd_(char *uplo, INTEGER *n, DOUBLE *ap, 
	DOUBLE *d__, DOUBLE *e, DOUBLE *tau, INTEGER *info);
 
 int dsptrf_(char *uplo, INTEGER *n, DOUBLE *ap, INTEGER *
	ipiv, INTEGER *info);
 
 int dsptri_(char *uplo, INTEGER *n, DOUBLE *ap, INTEGER *
	ipiv, DOUBLE *work, INTEGER *info);
 
 int dsptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *ap, INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int dstebz_(char *range, char *order, INTEGER *n, DOUBLE 
	*vl, DOUBLE *vu, INTEGER *il, INTEGER *iu, DOUBLE *abstol, 
	DOUBLE *d__, DOUBLE *e, INTEGER *m, INTEGER *nsplit, 
	DOUBLE *w, INTEGER *iblock, INTEGER *isplit, DOUBLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int dstedc_(char *compz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dstegr_(char *jobz, char *range, INTEGER *n, DOUBLE *
	d__, DOUBLE *e, DOUBLE *vl, DOUBLE *vu, INTEGER *il, 
	INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLE *z__, INTEGER *ldz, INTEGER *isuppz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dstein_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	INTEGER *m, DOUBLE *w, INTEGER *iblock, INTEGER *isplit, 
	DOUBLE *z__, INTEGER *ldz, DOUBLE *work, INTEGER *iwork, 
	INTEGER *ifail, INTEGER *info);
 
 int dsteqr_(char *compz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int dsterf_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	INTEGER *info);
 
 int dstev_(char *jobz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int dstevd_(char *jobz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLE *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dstevr_(char *jobz, char *range, INTEGER *n, DOUBLE *
	d__, DOUBLE *e, DOUBLE *vl, DOUBLE *vu, INTEGER *il, 
	INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLE *z__, INTEGER *ldz, INTEGER *isuppz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dstevx_(char *jobz, char *range, INTEGER *n, DOUBLE *
	d__, DOUBLE *e, DOUBLE *vl, DOUBLE *vu, INTEGER *il, 
	INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLE *z__, INTEGER *ldz, DOUBLE *work, INTEGER *iwork, 
	INTEGER *ifail, INTEGER *info);
 
 int dsycon_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, DOUBLE *anorm, DOUBLE *rcond, DOUBLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int dsyev_(char *jobz, char *uplo, INTEGER *n, DOUBLE *a,
	 INTEGER *lda, DOUBLE *w, DOUBLE *work, INTEGER *lwork, 
	INTEGER *info);
 
 int dsyevd_(char *jobz, char *uplo, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLE *w, DOUBLE *work, INTEGER *lwork, 
	INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dsyevr_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *vl, DOUBLE *vu, INTEGER *
	il, INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLE *z__, INTEGER *ldz, INTEGER *isuppz, DOUBLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int dsyevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *vl, DOUBLE *vu, INTEGER *
	il, INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLE *z__, INTEGER *ldz, DOUBLE *work, INTEGER *lwork, 
	INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int dsygs2_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int dsygst_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int dsygv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *w, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dsygvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *w, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int dsygvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER 
	*ldb, DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *iu, 
	DOUBLE *abstol, INTEGER *m, DOUBLE *w, DOUBLE *z__, 
	INTEGER *ldz, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *ifail, INTEGER *info);
 
 int dsyrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, INTEGER *
	ipiv, DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dsysv_(char *uplo, INTEGER *n, INTEGER *nrhs, DOUBLE 
	*a, INTEGER *lda, INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dsysvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *af, INTEGER *ldaf, 
	INTEGER *ipiv, DOUBLE *b, INTEGER *ldb, DOUBLE *x, INTEGER *
	ldx, DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int dsytd2_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *d__, DOUBLE *e, DOUBLE *tau, INTEGER *info);
 
 int dsytf2_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, INTEGER *info);
 
 int dsytrd_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *d__, DOUBLE *e, DOUBLE *tau, DOUBLE *
	work, INTEGER *lwork, INTEGER *info);
 
 int dsytrf_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dsytri_(char *uplo, INTEGER *n, DOUBLE *a, INTEGER *
	lda, INTEGER *ipiv, DOUBLE *work, INTEGER *info);
 
 int dsytrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *a, INTEGER *lda, INTEGER *ipiv, DOUBLE *b, INTEGER *
	ldb, INTEGER *info);
 
 int dtbcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	INTEGER *kd, DOUBLE *ab, INTEGER *ldab, DOUBLE *rcond, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dtbrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE 
	*b, INTEGER *ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dtbtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, DOUBLE *ab, INTEGER *ldab, DOUBLE 
	*b, INTEGER *ldb, INTEGER *info);
 
 int dtgevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *vl, INTEGER *ldvl, DOUBLE *vr, INTEGER *ldvr, INTEGER 
	*mm, INTEGER *m, DOUBLE *work, INTEGER *info);
 
 int dtgex2_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *
	q, INTEGER *ldq, DOUBLE *z__, INTEGER *ldz, INTEGER *j1, INTEGER *
	n1, INTEGER *n2, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dtgexc_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *
	q, INTEGER *ldq, DOUBLE *z__, INTEGER *ldz, INTEGER *ifst, 
	INTEGER *ilst, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
 int dtgsen_(INTEGER *ijob, LOGICAL *wantq, LOGICAL *wantz, 
	LOGICAL *select, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *
	b, INTEGER *ldb, DOUBLE *alphar, DOUBLE *alphai, DOUBLE *
	beta, DOUBLE *q, INTEGER *ldq, DOUBLE *z__, INTEGER *ldz, 
	INTEGER *m, DOUBLE *pl, DOUBLE *pr, DOUBLE *dif, 
	DOUBLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, 
	INTEGER *info);
 
 int dtgsja_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, INTEGER *k, INTEGER *l, DOUBLE *a, 
	INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLE *tola, 
	DOUBLE *tolb, DOUBLE *alpha, DOUBLE *beta, DOUBLE *u, 
	INTEGER *ldu, DOUBLE *v, INTEGER *ldv, DOUBLE *q, INTEGER *
	ldq, DOUBLE *work, INTEGER *ncycle, INTEGER *info);
 
 int dtgsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *vl, INTEGER *ldvl, DOUBLE *vr, INTEGER *ldvr, 
	DOUBLE *s, DOUBLE *dif, INTEGER *mm, INTEGER *m, DOUBLE *
	work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int dtgsy2_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *d__, INTEGER *ldd, 
	DOUBLE *e, INTEGER *lde, DOUBLE *f, INTEGER *ldf, DOUBLE *
	scale, DOUBLE *rdsum, DOUBLE *rdscal, INTEGER *iwork, INTEGER 
	*pq, INTEGER *info);
 
 int dtgsyl_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *c__, INTEGER *ldc, DOUBLE *d__, INTEGER *ldd, 
	DOUBLE *e, INTEGER *lde, DOUBLE *f, INTEGER *ldf, DOUBLE *
	scale, DOUBLE *dif, DOUBLE *work, INTEGER *lwork, INTEGER *
	iwork, INTEGER *info);
 
 int dtpcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	DOUBLE *ap, DOUBLE *rcond, DOUBLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int dtprfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLE *ap, DOUBLE *b, INTEGER *ldb, 
	DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dtptri_(char *uplo, char *diag, INTEGER *n, DOUBLE *
	ap, INTEGER *info);
 
 int dtptrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLE *ap, DOUBLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int dtrcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	DOUBLE *a, INTEGER *lda, DOUBLE *rcond, DOUBLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int dtrevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLE *t, INTEGER *ldt, DOUBLE *vl, INTEGER *
	ldvl, DOUBLE *vr, INTEGER *ldvr, INTEGER *mm, INTEGER *m, 
	DOUBLE *work, INTEGER *info);
 
 int dtrexc_(char *compq, INTEGER *n, DOUBLE *t, INTEGER *
	ldt, DOUBLE *q, INTEGER *ldq, INTEGER *ifst, INTEGER *ilst, 
	DOUBLE *work, INTEGER *info);
 
 int dtrrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *
	ldb, DOUBLE *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int dtrsen_(char *job, char *compq, LOGICAL *select, INTEGER 
	*n, DOUBLE *t, INTEGER *ldt, DOUBLE *q, INTEGER *ldq, 
	DOUBLE *wr, DOUBLE *wi, INTEGER *m, DOUBLE *s, DOUBLE 
	*sep, DOUBLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *
	liwork, INTEGER *info);
 
 int dtrsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLE *t, INTEGER *ldt, DOUBLE *vl, INTEGER *
	ldvl, DOUBLE *vr, INTEGER *ldvr, DOUBLE *s, DOUBLE *sep, 
	INTEGER *mm, INTEGER *m, DOUBLE *work, INTEGER *ldwork, INTEGER *
	iwork, INTEGER *info);
 
 int dtrsyl_(char *trana, char *tranb, INTEGER *isgn, INTEGER 
	*m, INTEGER *n, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *
	ldb, DOUBLE *c__, INTEGER *ldc, DOUBLE *scale, INTEGER *info);
 
 int dtrti2_(char *uplo, char *diag, INTEGER *n, DOUBLE *
	a, INTEGER *lda, INTEGER *info);
 
 int dtrtri_(char *uplo, char *diag, INTEGER *n, DOUBLE *
	a, INTEGER *lda, INTEGER *info);
 
 int dtrtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLE *a, INTEGER *lda, DOUBLE *b, INTEGER *
	ldb, INTEGER *info);
 
 int dtzrqf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, INTEGER *info);
 
 int dtzrzf_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLE *tau, DOUBLE *work, INTEGER *lwork, INTEGER *info);
 
INTEGER icmax1_(INTEGER *n, COMPLEX *cx, INTEGER *incx);
 
INTEGER ieeeck_(INTEGER *ispec, SINGLE *zero, SINGLE *one);
 
INTEGER ilaenv_(INTEGER *ispec, char *name__, char *opts, INTEGER *n1, 
	INTEGER *n2, INTEGER *n3, INTEGER *n4, __CLPK_ftnlen name_len, __CLPK_ftnlen 
	opts_len);
 
INTEGER izmax1_(INTEGER *n, DOUBLECOMPLEX *cx, INTEGER *incx);
 
 int sbdsdc_(char *uplo, char *compq, INTEGER *n, SINGLE *d__, 
	SINGLE *e, SINGLE *u, INTEGER *ldu, SINGLE *vt, INTEGER *ldvt, SINGLE *q, 
	INTEGER *iq, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sbdsqr_(char *uplo, INTEGER *n, INTEGER *ncvt, INTEGER *
	nru, INTEGER *ncc, SINGLE *d__, SINGLE *e, SINGLE *vt, INTEGER *ldvt, SINGLE *
	u, INTEGER *ldu, SINGLE *c__, INTEGER *ldc, SINGLE *work, INTEGER *info);
 
 int sdisna_(char *job, INTEGER *m, INTEGER *n, SINGLE *d__, 
	SINGLE *sep, INTEGER *info);
 
 int sgbbrd_(char *vect, INTEGER *m, INTEGER *n, INTEGER *ncc,
	 INTEGER *kl, INTEGER *ku, SINGLE *ab, INTEGER *ldab, SINGLE *d__, SINGLE *
	e, SINGLE *q, INTEGER *ldq, SINGLE *pt, INTEGER *ldpt, SINGLE *c__, INTEGER 
	*ldc, SINGLE *work, INTEGER *info);
 
 int sgbcon_(char *norm, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 SINGLE *ab, INTEGER *ldab, INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, 
	SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sgbequ_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 SINGLE *ab, INTEGER *ldab, SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *
	colcnd, SINGLE *amax, INTEGER *info);
 
 int sgbrfs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *afb, INTEGER *ldafb,
	 INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *
	ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sgbsv_(INTEGER *n, INTEGER *kl, INTEGER *ku, INTEGER *
	nrhs, SINGLE *ab, INTEGER *ldab, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, 
	INTEGER *info);
 
 int sgbsvx_(char *fact, char *trans, INTEGER *n, INTEGER *kl,
	 INTEGER *ku, INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *afb, 
	INTEGER *ldafb, INTEGER *ipiv, char *equed, SINGLE *r__, SINGLE *c__, 
	SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr,
	 SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sgbtf2_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 SINGLE *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int sgbtrf_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 SINGLE *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int sgbtrs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, INTEGER *ipiv, SINGLE *b, 
	INTEGER *ldb, INTEGER *info);
 
 int sgebak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, SINGLE *scale, INTEGER *m, SINGLE *v, INTEGER *ldv, INTEGER 
	*info);
 
 int sgebal_(char *job, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ilo, INTEGER *ihi, SINGLE *scale, INTEGER *info);
 
 int sgebd2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *d__, SINGLE *e, SINGLE *tauq, SINGLE *taup, SINGLE *work, INTEGER *info);
 
 int sgebrd_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *d__, SINGLE *e, SINGLE *tauq, SINGLE *taup, SINGLE *work, INTEGER *
	lwork, INTEGER *info);
 
 int sgecon_(char *norm, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *anorm, SINGLE *rcond, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sgeequ_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *colcnd, SINGLE *amax, INTEGER 
	*info);
 
 int sgees_(char *jobvs, char *sort, __CLPK_L_fp select, INTEGER *n, 
	SINGLE *a, INTEGER *lda, INTEGER *sdim, SINGLE *wr, SINGLE *wi, SINGLE *vs, 
	INTEGER *ldvs, SINGLE *work, INTEGER *lwork, LOGICAL *bwork, INTEGER *
	info);
 
 int sgeesx_(char *jobvs, char *sort, __CLPK_L_fp select, char *
	sense, INTEGER *n, SINGLE *a, INTEGER *lda, INTEGER *sdim, SINGLE *wr, 
	SINGLE *wi, SINGLE *vs, INTEGER *ldvs, SINGLE *rconde, SINGLE *rcondv, SINGLE *
	work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, LOGICAL *bwork,
	 INTEGER *info);
 
 int sgeev_(char *jobvl, char *jobvr, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *wr, SINGLE *wi, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, 
	INTEGER *ldvr, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *wr, SINGLE *wi, SINGLE *
	vl, INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, INTEGER *ilo, INTEGER *
	ihi, SINGLE *scale, SINGLE *abnrm, SINGLE *rconde, SINGLE *rcondv, SINGLE *work,
	 INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int sgegs_(char *jobvsl, char *jobvsr, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *alphar, SINGLE *alphai, SINGLE 
	*beta, SINGLE *vsl, INTEGER *ldvsl, SINGLE *vsr, INTEGER *ldvsr, SINGLE *
	work, INTEGER *lwork, INTEGER *info);
 
 int sgegv_(char *jobvl, char *jobvr, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *alphar, SINGLE *alphai, SINGLE 
	*beta, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int sgehd2_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgehrd_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgelq2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgelqf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgels_(char *trans, INTEGER *m, INTEGER *n, INTEGER *
	nrhs, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int sgelsd_(INTEGER *m, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *s, SINGLE *rcond, INTEGER *
	rank, SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int sgelss_(INTEGER *m, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *s, SINGLE *rcond, INTEGER *
	rank, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgelsx_(INTEGER *m, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *jpvt, SINGLE *rcond, 
	INTEGER *rank, SINGLE *work, INTEGER *info);
 
 int sgelsy_(INTEGER *m, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *jpvt, SINGLE *rcond, 
	INTEGER *rank, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgeql2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgeqlf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgeqp3_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *jpvt, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgeqpf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *jpvt, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgeqr2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgeqrf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgerfs_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *af, INTEGER *ldaf, INTEGER *ipiv, SINGLE *b, 
	INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int sgerq2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sgerqf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgesc2_(INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *rhs, 
	INTEGER *ipiv, INTEGER *jpiv, SINGLE *scale);
 
 int sgesdd_(char *jobz, INTEGER *m, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *s, SINGLE *u, INTEGER *ldu, SINGLE *vt, INTEGER *ldvt,
	 SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int sgesv_(INTEGER *n, INTEGER *nrhs, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sgesvd_(char *jobu, char *jobvt, INTEGER *m, INTEGER *n, 
	SINGLE *a, INTEGER *lda, SINGLE *s, SINGLE *u, INTEGER *ldu, SINGLE *vt, 
	INTEGER *ldvt, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgesvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, SINGLE *a, INTEGER *lda, SINGLE *af, INTEGER *ldaf, INTEGER *ipiv, 
	char *equed, SINGLE *r__, SINGLE *c__, SINGLE *b, INTEGER *ldb, SINGLE *x, 
	INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int sgetc2_(INTEGER *n, SINGLE *a, INTEGER *lda, INTEGER *ipiv,
	 INTEGER *jpiv, INTEGER *info);
 
 int sgetf2_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, INTEGER *info);
 
 int sgetrf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, INTEGER *info);
 
 int sgetri_(INTEGER *n, SINGLE *a, INTEGER *lda, INTEGER *ipiv,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgetrs_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sggbak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, SINGLE *lscale, SINGLE *rscale, INTEGER *m, SINGLE *v, 
	INTEGER *ldv, INTEGER *info);
 
 int sggbal_(char *job, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *b, INTEGER *ldb, INTEGER *ilo, INTEGER *ihi, SINGLE *lscale, SINGLE 
	*rscale, SINGLE *work, INTEGER *info);
 
 int sgges_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	selctg, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, 
	INTEGER *sdim, SINGLE *alphar, SINGLE *alphai, SINGLE *beta, SINGLE *vsl, 
	INTEGER *ldvsl, SINGLE *vsr, INTEGER *ldvsr, SINGLE *work, INTEGER *lwork,
	 LOGICAL *bwork, INTEGER *info);
 
 int sggesx_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	selctg, char *sense, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, 
	INTEGER *ldb, INTEGER *sdim, SINGLE *alphar, SINGLE *alphai, SINGLE *beta, 
	SINGLE *vsl, INTEGER *ldvsl, SINGLE *vsr, INTEGER *ldvsr, SINGLE *rconde, 
	SINGLE *rcondv, SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *
	liwork, LOGICAL *bwork, INTEGER *info);
 
 int sggev_(char *jobvl, char *jobvr, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *alphar, SINGLE *alphai, SINGLE 
	*beta, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE 
	*alphar, SINGLE *alphai, SINGLE *beta, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, 
	INTEGER *ldvr, INTEGER *ilo, INTEGER *ihi, SINGLE *lscale, SINGLE *rscale,
	 SINGLE *abnrm, SINGLE *bbnrm, SINGLE *rconde, SINGLE *rcondv, SINGLE *work, 
	INTEGER *lwork, INTEGER *iwork, LOGICAL *bwork, INTEGER *info);
 
 int sggglm_(INTEGER *n, INTEGER *m, INTEGER *p, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *d__, SINGLE *x, SINGLE *y, 
	SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sgghrd_(char *compq, char *compz, INTEGER *n, INTEGER *
	ilo, INTEGER *ihi, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE 
	*q, INTEGER *ldq, SINGLE *z__, INTEGER *ldz, INTEGER *info);
 
 int sgglse_(INTEGER *m, INTEGER *n, INTEGER *p, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *c__, SINGLE *d__, SINGLE *x, 
	SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sggqrf_(INTEGER *n, INTEGER *m, INTEGER *p, SINGLE *a, 
	INTEGER *lda, SINGLE *taua, SINGLE *b, INTEGER *ldb, SINGLE *taub, SINGLE *
	work, INTEGER *lwork, INTEGER *info);
 
 int sggrqf_(INTEGER *m, INTEGER *p, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *taua, SINGLE *b, INTEGER *ldb, SINGLE *taub, SINGLE *
	work, INTEGER *lwork, INTEGER *info);
 
 int sggsvd_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *n, INTEGER *p, INTEGER *k, INTEGER *l, SINGLE *a, INTEGER *lda,
	 SINGLE *b, INTEGER *ldb, SINGLE *alpha, SINGLE *beta, SINGLE *u, INTEGER *
	ldu, SINGLE *v, INTEGER *ldv, SINGLE *q, INTEGER *ldq, SINGLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int sggsvp_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, 
	SINGLE *tola, SINGLE *tolb, INTEGER *k, INTEGER *l, SINGLE *u, INTEGER *ldu,
	 SINGLE *v, INTEGER *ldv, SINGLE *q, INTEGER *ldq, INTEGER *iwork, SINGLE *
	tau, SINGLE *work, INTEGER *info);
 
 int sgtcon_(char *norm, INTEGER *n, SINGLE *dl, SINGLE *d__, 
	SINGLE *du, SINGLE *du2, INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, SINGLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int sgtrfs_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *dl,
	 SINGLE *d__, SINGLE *du, SINGLE *dlf, SINGLE *df, SINGLE *duf, SINGLE *du2, 
	INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *
	ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sgtsv_(INTEGER *n, INTEGER *nrhs, SINGLE *dl, SINGLE *d__, 
	SINGLE *du, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sgtsvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, SINGLE *dl, SINGLE *d__, SINGLE *du, SINGLE *dlf, SINGLE *df, SINGLE *duf, 
	SINGLE *du2, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *
	ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int sgttrf_(INTEGER *n, SINGLE *dl, SINGLE *d__, SINGLE *du, SINGLE *
	du2, INTEGER *ipiv, INTEGER *info);
 
 int sgttrs_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *dl,
	 SINGLE *d__, SINGLE *du, SINGLE *du2, INTEGER *ipiv, SINGLE *b, INTEGER *ldb,
	 INTEGER *info);
 
 int sgtts2_(INTEGER *itrans, INTEGER *n, INTEGER *nrhs, SINGLE 
	*dl, SINGLE *d__, SINGLE *du, SINGLE *du2, INTEGER *ipiv, SINGLE *b, INTEGER *
	ldb);
 
 int shgeqz_(char *job, char *compq, char *compz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *
	ldb, SINGLE *alphar, SINGLE *alphai, SINGLE *beta, SINGLE *q, INTEGER *ldq, 
	SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int shsein_(char *side, char *eigsrc, char *initv, LOGICAL *
	select, INTEGER *n, SINGLE *h__, INTEGER *ldh, SINGLE *wr, SINGLE *wi, SINGLE 
	*vl, INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, INTEGER *mm, INTEGER *m, 
	SINGLE *work, INTEGER *ifaill, INTEGER *ifailr, INTEGER *info);
 
 int shseqr_(char *job, char *compz, INTEGER *n, INTEGER *ilo,
	 INTEGER *ihi, SINGLE *h__, INTEGER *ldh, SINGLE *wr, SINGLE *wi, SINGLE *z__,
	 INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int slabad_(SINGLE *small, SINGLE *large);
 
 int slabrd_(INTEGER *m, INTEGER *n, INTEGER *nb, SINGLE *a, 
	INTEGER *lda, SINGLE *d__, SINGLE *e, SINGLE *tauq, SINGLE *taup, SINGLE *x, 
	INTEGER *ldx, SINGLE *y, INTEGER *ldy);
 
 int slacon_(INTEGER *n, SINGLE *v, SINGLE *x, INTEGER *isgn, 
	SINGLE *est, INTEGER *kase);
 
 int slacpy_(char *uplo, INTEGER *m, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb);
 
 int sladiv_(SINGLE *a, SINGLE *b, SINGLE *c__, SINGLE *d__, SINGLE *p, 
	SINGLE *q);
 
 int slae2_(SINGLE *a, SINGLE *b, SINGLE *c__, SINGLE *rt1, SINGLE *rt2);
 
 int slaebz_(INTEGER *ijob, INTEGER *nitmax, INTEGER *n, 
	INTEGER *mmax, INTEGER *minp, INTEGER *nbmin, SINGLE *abstol, SINGLE *
	reltol, SINGLE *pivmin, SINGLE *d__, SINGLE *e, SINGLE *e2, INTEGER *nval, 
	SINGLE *ab, SINGLE *c__, INTEGER *mout, INTEGER *nab, SINGLE *work, INTEGER 
	*iwork, INTEGER *info);
 
 int slaed0_(INTEGER *icompq, INTEGER *qsiz, INTEGER *n, SINGLE 
	*d__, SINGLE *e, SINGLE *q, INTEGER *ldq, SINGLE *qstore, INTEGER *ldqs, 
	SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int slaed1_(INTEGER *n, SINGLE *d__, SINGLE *q, INTEGER *ldq, 
	INTEGER *indxq, SINGLE *rho, INTEGER *cutpnt, SINGLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int slaed2_(INTEGER *k, INTEGER *n, INTEGER *n1, SINGLE *d__, 
	SINGLE *q, INTEGER *ldq, INTEGER *indxq, SINGLE *rho, SINGLE *z__, SINGLE *
	dlamda, SINGLE *w, SINGLE *q2, INTEGER *indx, INTEGER *indxc, INTEGER *
	indxp, INTEGER *coltyp, INTEGER *info);
 
 int slaed3_(INTEGER *k, INTEGER *n, INTEGER *n1, SINGLE *d__, 
	SINGLE *q, INTEGER *ldq, SINGLE *rho, SINGLE *dlamda, SINGLE *q2, INTEGER *
	indx, INTEGER *ctot, SINGLE *w, SINGLE *s, INTEGER *info);
 
 int slaed4_(INTEGER *n, INTEGER *i__, SINGLE *d__, SINGLE *z__, 
	SINGLE *delta, SINGLE *rho, SINGLE *dlam, INTEGER *info);
 
 int slaed5_(INTEGER *i__, SINGLE *d__, SINGLE *z__, SINGLE *delta, 
	SINGLE *rho, SINGLE *dlam);
 
 int slaed6_(INTEGER *kniter, LOGICAL *orgati, SINGLE *rho, 
	SINGLE *d__, SINGLE *z__, SINGLE *finit, SINGLE *tau, INTEGER *info);
 
 int slaed7_(INTEGER *icompq, INTEGER *n, INTEGER *qsiz, 
	INTEGER *tlvls, INTEGER *curlvl, INTEGER *curpbm, SINGLE *d__, SINGLE *q, 
	INTEGER *ldq, INTEGER *indxq, SINGLE *rho, INTEGER *cutpnt, SINGLE *
	qstore, INTEGER *qptr, INTEGER *prmptr, INTEGER *perm, INTEGER *
	givptr, INTEGER *givcol, SINGLE *givnum, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int slaed8_(INTEGER *icompq, INTEGER *k, INTEGER *n, INTEGER 
	*qsiz, SINGLE *d__, SINGLE *q, INTEGER *ldq, INTEGER *indxq, SINGLE *rho, 
	INTEGER *cutpnt, SINGLE *z__, SINGLE *dlamda, SINGLE *q2, INTEGER *ldq2, 
	SINGLE *w, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, SINGLE *
	givnum, INTEGER *indxp, INTEGER *indx, INTEGER *info);
 
 int slaed9_(INTEGER *k, INTEGER *kstart, INTEGER *kstop, 
	INTEGER *n, SINGLE *d__, SINGLE *q, INTEGER *ldq, SINGLE *rho, SINGLE *dlamda,
	 SINGLE *w, SINGLE *s, INTEGER *lds, INTEGER *info);
 
 int slaeda_(INTEGER *n, INTEGER *tlvls, INTEGER *curlvl, 
	INTEGER *curpbm, INTEGER *prmptr, INTEGER *perm, INTEGER *givptr, 
	INTEGER *givcol, SINGLE *givnum, SINGLE *q, INTEGER *qptr, SINGLE *z__, 
	SINGLE *ztemp, INTEGER *info);
 
 int slaein_(LOGICAL *rightv, LOGICAL *noinit, INTEGER *n, 
	SINGLE *h__, INTEGER *ldh, SINGLE *wr, SINGLE *wi, SINGLE *vr, SINGLE *vi, SINGLE 
	*b, INTEGER *ldb, SINGLE *work, SINGLE *eps3, SINGLE *smlnum, SINGLE *bignum, 
	INTEGER *info);
 
 int slaev2_(SINGLE *a, SINGLE *b, SINGLE *c__, SINGLE *rt1, SINGLE *
	rt2, SINGLE *cs1, SINGLE *sn1);
 
 int slaexc_(LOGICAL *wantq, INTEGER *n, SINGLE *t, INTEGER *
	ldt, SINGLE *q, INTEGER *ldq, INTEGER *j1, INTEGER *n1, INTEGER *n2, 
	SINGLE *work, INTEGER *info);
 
 int slag2_(SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, 
	SINGLE *safmin, SINGLE *scale1, SINGLE *scale2, SINGLE *wr1, SINGLE *wr2, SINGLE *
	wi);
 
 int slags2_(LOGICAL *upper, SINGLE *a1, SINGLE *a2, SINGLE *a3, 
	SINGLE *b1, SINGLE *b2, SINGLE *b3, SINGLE *csu, SINGLE *snu, SINGLE *csv, SINGLE *
	snv, SINGLE *csq, SINGLE *snq);
 
 int slagtf_(INTEGER *n, SINGLE *a, SINGLE *lambda, SINGLE *b, SINGLE 
	*c__, SINGLE *tol, SINGLE *d__, INTEGER *in, INTEGER *info);
 
 int slagtm_(char *trans, INTEGER *n, INTEGER *nrhs, SINGLE *
	alpha, SINGLE *dl, SINGLE *d__, SINGLE *du, SINGLE *x, INTEGER *ldx, SINGLE *
	beta, SINGLE *b, INTEGER *ldb);
 
 int slagts_(INTEGER *job, INTEGER *n, SINGLE *a, SINGLE *b, SINGLE 
	*c__, SINGLE *d__, INTEGER *in, SINGLE *y, SINGLE *tol, INTEGER *info);
 
 int slagv2_(SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, 
	SINGLE *alphar, SINGLE *alphai, SINGLE *beta, SINGLE *csl, SINGLE *snl, SINGLE *
	csr, SINGLE *snr);
 
 int slahqr_(LOGICAL *wantt, LOGICAL *wantz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, SINGLE *h__, INTEGER *ldh, SINGLE *wr, SINGLE *
	wi, INTEGER *iloz, INTEGER *ihiz, SINGLE *z__, INTEGER *ldz, INTEGER *
	info);
 
 int slahrd_(INTEGER *n, INTEGER *k, INTEGER *nb, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *t, INTEGER *ldt, SINGLE *y, INTEGER *ldy);
 
 int slaic1_(INTEGER *job, INTEGER *j, SINGLE *x, SINGLE *sest, 
	SINGLE *w, SINGLE *gamma, SINGLE *sestpr, SINGLE *s, SINGLE *c__);
 
 int slaln2_(LOGICAL *ltrans, INTEGER *na, INTEGER *nw, SINGLE *
	smin, SINGLE *ca, SINGLE *a, INTEGER *lda, SINGLE *d1, SINGLE *d2, SINGLE *b, 
	INTEGER *ldb, SINGLE *wr, SINGLE *wi, SINGLE *x, INTEGER *ldx, SINGLE *scale, 
	SINGLE *xnorm, INTEGER *info);
 
 int slals0_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *nrhs, SINGLE *b, INTEGER *ldb, SINGLE *bx, 
	INTEGER *ldbx, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, 
	INTEGER *ldgcol, SINGLE *givnum, INTEGER *ldgnum, SINGLE *poles, SINGLE *
	difl, SINGLE *difr, SINGLE *z__, INTEGER *k, SINGLE *c__, SINGLE *s, SINGLE *
	work, INTEGER *info);
 
 int slalsa_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *nrhs, SINGLE *b, INTEGER *ldb, SINGLE *bx, INTEGER *ldbx, SINGLE *
	u, INTEGER *ldu, SINGLE *vt, INTEGER *k, SINGLE *difl, SINGLE *difr, SINGLE *
	z__, SINGLE *poles, INTEGER *givptr, INTEGER *givcol, INTEGER *ldgcol, 
	INTEGER *perm, SINGLE *givnum, SINGLE *c__, SINGLE *s, SINGLE *work, INTEGER *
	iwork, INTEGER *info);
 
 int slalsd_(char *uplo, INTEGER *smlsiz, INTEGER *n, INTEGER 
	*nrhs, SINGLE *d__, SINGLE *e, SINGLE *b, INTEGER *ldb, SINGLE *rcond, 
	INTEGER *rank, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int slamc1_(INTEGER *beta, INTEGER *t, LOGICAL *rnd, LOGICAL 
	*ieee1);
 
 int slamc2_(INTEGER *beta, INTEGER *t, LOGICAL *rnd, SINGLE *
	eps, INTEGER *emin, SINGLE *rmin, INTEGER *emax, SINGLE *rmax);
 
 int slamc4_(INTEGER *emin, SINGLE *start, INTEGER *base);
 
 int slamc5_(INTEGER *beta, INTEGER *p, INTEGER *emin, 
	LOGICAL *ieee, INTEGER *emax, SINGLE *rmax);
 
 int slamrg_(INTEGER *n1, INTEGER *n2, SINGLE *a, INTEGER *
	strd1, INTEGER *strd2, INTEGER *index);
 
 int slanv2_(SINGLE *a, SINGLE *b, SINGLE *c__, SINGLE *d__, SINGLE *
	rt1r, SINGLE *rt1i, SINGLE *rt2r, SINGLE *rt2i, SINGLE *cs, SINGLE *sn);
 
 int slapll_(INTEGER *n, SINGLE *x, INTEGER *incx, SINGLE *y, 
	INTEGER *incy, SINGLE *ssmin);
 
 int slapmt_(LOGICAL *forwrd, INTEGER *m, INTEGER *n, SINGLE *x,
	 INTEGER *ldx, INTEGER *k);
 
 int slaqgb_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 SINGLE *ab, INTEGER *ldab, SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *
	colcnd, SINGLE *amax, char *equed);
 
 int slaqge_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *r__, SINGLE *c__, SINGLE *rowcnd, SINGLE *colcnd, SINGLE *amax, char *
	equed);
 
 int slaqp2_(INTEGER *m, INTEGER *n, INTEGER *offset, SINGLE *a,
	 INTEGER *lda, INTEGER *jpvt, SINGLE *tau, SINGLE *vn1, SINGLE *vn2, SINGLE *
	work);
 
 int slaqps_(INTEGER *m, INTEGER *n, INTEGER *offset, INTEGER 
	*nb, INTEGER *kb, SINGLE *a, INTEGER *lda, INTEGER *jpvt, SINGLE *tau, 
	SINGLE *vn1, SINGLE *vn2, SINGLE *auxv, SINGLE *f, INTEGER *ldf);
 
 int slaqsb_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int slaqsp_(char *uplo, INTEGER *n, SINGLE *ap, SINGLE *s, SINGLE *
	scond, SINGLE *amax, char *equed);
 
 int slaqsy_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *s, SINGLE *scond, SINGLE *amax, char *equed);
 
 int slaqtr_(LOGICAL *ltran, LOGICAL *lreal, INTEGER *n, SINGLE 
	*t, INTEGER *ldt, SINGLE *b, SINGLE *w, SINGLE *scale, SINGLE *x, SINGLE *work, 
	INTEGER *info);
 
 int slar1v_(INTEGER *n, INTEGER *b1, INTEGER *bn, SINGLE *
	sigma, SINGLE *d__, SINGLE *l, SINGLE *ld, SINGLE *lld, SINGLE *gersch, SINGLE *
	z__, SINGLE *ztz, SINGLE *mingma, INTEGER *r__, INTEGER *isuppz, SINGLE *
	work);
 
 int slar2v_(INTEGER *n, SINGLE *x, SINGLE *y, SINGLE *z__, INTEGER 
	*incx, SINGLE *c__, SINGLE *s, INTEGER *incc);
 
 int slarf_(char *side, INTEGER *m, INTEGER *n, SINGLE *v, 
	INTEGER *incv, SINGLE *tau, SINGLE *c__, INTEGER *ldc, SINGLE *work);
 
 int slarfb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *v, INTEGER *ldv, 
	SINGLE *t, INTEGER *ldt, SINGLE *c__, INTEGER *ldc, SINGLE *work, INTEGER *
	ldwork);
 
 int slarfg_(INTEGER *n, SINGLE *alpha, SINGLE *x, INTEGER *incx, 
	SINGLE *tau);
 
 int slarft_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, SINGLE *v, INTEGER *ldv, SINGLE *tau, SINGLE *t, INTEGER *ldt);
 
 int slarfx_(char *side, INTEGER *m, INTEGER *n, SINGLE *v, 
	SINGLE *tau, SINGLE *c__, INTEGER *ldc, SINGLE *work);
 
 int slargv_(INTEGER *n, SINGLE *x, INTEGER *incx, SINGLE *y, 
	INTEGER *incy, SINGLE *c__, INTEGER *incc);
 
 int slarnv_(INTEGER *idist, INTEGER *iseed, INTEGER *n, SINGLE 
	*x);
 
 int slarrb_(INTEGER *n, SINGLE *d__, SINGLE *l, SINGLE *ld, SINGLE *
	lld, INTEGER *ifirst, INTEGER *ilast, SINGLE *sigma, SINGLE *reltol, SINGLE 
	*w, SINGLE *wgap, SINGLE *werr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int slarre_(INTEGER *n, SINGLE *d__, SINGLE *e, SINGLE *tol, 
	INTEGER *nsplit, INTEGER *isplit, INTEGER *m, SINGLE *w, SINGLE *woff, 
	SINGLE *gersch, SINGLE *work, INTEGER *info);
 
 int slarrf_(INTEGER *n, SINGLE *d__, SINGLE *l, SINGLE *ld, SINGLE *
	lld, INTEGER *ifirst, INTEGER *ilast, SINGLE *w, SINGLE *dplus, SINGLE *
	lplus, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int slarrv_(INTEGER *n, SINGLE *d__, SINGLE *l, INTEGER *isplit, 
	INTEGER *m, SINGLE *w, INTEGER *iblock, SINGLE *gersch, SINGLE *tol, SINGLE *
	z__, INTEGER *ldz, INTEGER *isuppz, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int slartg_(SINGLE *f, SINGLE *g, SINGLE *cs, SINGLE *sn, SINGLE *r__);
 
 int slartv_(INTEGER *n, SINGLE *x, INTEGER *incx, SINGLE *y, 
	INTEGER *incy, SINGLE *c__, SINGLE *s, INTEGER *incc);
 
 int slaruv_(INTEGER *iseed, INTEGER *n, SINGLE *x);
 
 int slarz_(char *side, INTEGER *m, INTEGER *n, INTEGER *l, 
	SINGLE *v, INTEGER *incv, SINGLE *tau, SINGLE *c__, INTEGER *ldc, SINGLE *
	work);
 
 int slarzb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, INTEGER *l, SINGLE *v, 
	INTEGER *ldv, SINGLE *t, INTEGER *ldt, SINGLE *c__, INTEGER *ldc, SINGLE *
	work, INTEGER *ldwork);
 
 int slarzt_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, SINGLE *v, INTEGER *ldv, SINGLE *tau, SINGLE *t, INTEGER *ldt);
 
 int slas2_(SINGLE *f, SINGLE *g, SINGLE *h__, SINGLE *ssmin, SINGLE *
	ssmax);
 
 int slascl_(char *type__, INTEGER *kl, INTEGER *ku, SINGLE *
	cfrom, SINGLE *cto, INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int slasd0_(INTEGER *n, INTEGER *sqre, SINGLE *d__, SINGLE *e, 
	SINGLE *u, INTEGER *ldu, SINGLE *vt, INTEGER *ldvt, INTEGER *smlsiz, 
	INTEGER *iwork, SINGLE *work, INTEGER *info);
 
 int slasd1_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, SINGLE *
	d__, SINGLE *alpha, SINGLE *beta, SINGLE *u, INTEGER *ldu, SINGLE *vt, 
	INTEGER *ldvt, INTEGER *idxq, INTEGER *iwork, SINGLE *work, INTEGER *
	info);
 
 int slasd2_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, INTEGER 
	*k, SINGLE *d__, SINGLE *z__, SINGLE *alpha, SINGLE *beta, SINGLE *u, INTEGER *
	ldu, SINGLE *vt, INTEGER *ldvt, SINGLE *dsigma, SINGLE *u2, INTEGER *ldu2, 
	SINGLE *vt2, INTEGER *ldvt2, INTEGER *idxp, INTEGER *idx, INTEGER *idxc,
	 INTEGER *idxq, INTEGER *coltyp, INTEGER *info);
 
 int slasd3_(INTEGER *nl, INTEGER *nr, INTEGER *sqre, INTEGER 
	*k, SINGLE *d__, SINGLE *q, INTEGER *ldq, SINGLE *dsigma, SINGLE *u, INTEGER *
	ldu, SINGLE *u2, INTEGER *ldu2, SINGLE *vt, INTEGER *ldvt, SINGLE *vt2, 
	INTEGER *ldvt2, INTEGER *idxc, INTEGER *ctot, SINGLE *z__, INTEGER *
	info);
 
 int slasd4_(INTEGER *n, INTEGER *i__, SINGLE *d__, SINGLE *z__, 
	SINGLE *delta, SINGLE *rho, SINGLE *sigma, SINGLE *work, INTEGER *info);
 
 int slasd5_(INTEGER *i__, SINGLE *d__, SINGLE *z__, SINGLE *delta, 
	SINGLE *rho, SINGLE *dsigma, SINGLE *work);
 
 int slasd6_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, SINGLE *d__, SINGLE *vf, SINGLE *vl, SINGLE *alpha, SINGLE *beta,
	 INTEGER *idxq, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, 
	INTEGER *ldgcol, SINGLE *givnum, INTEGER *ldgnum, SINGLE *poles, SINGLE *
	difl, SINGLE *difr, SINGLE *z__, INTEGER *k, SINGLE *c__, SINGLE *s, SINGLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int slasd7_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *k, SINGLE *d__, SINGLE *z__, SINGLE *zw, SINGLE *vf, 
	SINGLE *vfw, SINGLE *vl, SINGLE *vlw, SINGLE *alpha, SINGLE *beta, SINGLE *dsigma,
	 INTEGER *idx, INTEGER *idxp, INTEGER *idxq, INTEGER *perm, INTEGER *
	givptr, INTEGER *givcol, INTEGER *ldgcol, SINGLE *givnum, INTEGER *
	ldgnum, SINGLE *c__, SINGLE *s, INTEGER *info);
 
 int slasd8_(INTEGER *icompq, INTEGER *k, SINGLE *d__, SINGLE *
	z__, SINGLE *vf, SINGLE *vl, SINGLE *difl, SINGLE *difr, INTEGER *lddifr, 
	SINGLE *dsigma, SINGLE *work, INTEGER *info);
 
 int slasd9_(INTEGER *icompq, INTEGER *ldu, INTEGER *k, SINGLE *
	d__, SINGLE *z__, SINGLE *vf, SINGLE *vl, SINGLE *difl, SINGLE *difr, SINGLE *
	dsigma, SINGLE *work, INTEGER *info);
 
 int slasda_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *sqre, SINGLE *d__, SINGLE *e, SINGLE *u, INTEGER *ldu, SINGLE *vt, 
	INTEGER *k, SINGLE *difl, SINGLE *difr, SINGLE *z__, SINGLE *poles, INTEGER *
	givptr, INTEGER *givcol, INTEGER *ldgcol, INTEGER *perm, SINGLE *givnum,
	 SINGLE *c__, SINGLE *s, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int slasdq_(char *uplo, INTEGER *sqre, INTEGER *n, INTEGER *
	ncvt, INTEGER *nru, INTEGER *ncc, SINGLE *d__, SINGLE *e, SINGLE *vt, 
	INTEGER *ldvt, SINGLE *u, INTEGER *ldu, SINGLE *c__, INTEGER *ldc, SINGLE *
	work, INTEGER *info);
 
 int slasdt_(INTEGER *n, INTEGER *lvl, INTEGER *nd, INTEGER *
	inode, INTEGER *ndiml, INTEGER *ndimr, INTEGER *msub);
 
 int slaset_(char *uplo, INTEGER *m, INTEGER *n, SINGLE *alpha, 
	SINGLE *beta, SINGLE *a, INTEGER *lda);
 
 int slasq1_(INTEGER *n, SINGLE *d__, SINGLE *e, SINGLE *work, 
	INTEGER *info);
 
 int slasq2_(INTEGER *n, SINGLE *z__, INTEGER *info);
 
 int slasq3_(INTEGER *i0, INTEGER *n0, SINGLE *z__, INTEGER *pp,
	 SINGLE *dmin__, SINGLE *sigma, SINGLE *desig, SINGLE *qmax, INTEGER *nfail, 
	INTEGER *iter, INTEGER *ndiv, LOGICAL *ieee);
 
 int slasq4_(INTEGER *i0, INTEGER *n0, SINGLE *z__, INTEGER *pp,
	 INTEGER *n0in, SINGLE *dmin__, SINGLE *dmin1, SINGLE *dmin2, SINGLE *dn, 
	SINGLE *dn1, SINGLE *dn2, SINGLE *tau, INTEGER *ttype);
 
 int slasq5_(INTEGER *i0, INTEGER *n0, SINGLE *z__, INTEGER *pp,
	 SINGLE *tau, SINGLE *dmin__, SINGLE *dmin1, SINGLE *dmin2, SINGLE *dn, SINGLE *
	dnm1, SINGLE *dnm2, LOGICAL *ieee);
 
 int slasq6_(INTEGER *i0, INTEGER *n0, SINGLE *z__, INTEGER *pp,
	 SINGLE *dmin__, SINGLE *dmin1, SINGLE *dmin2, SINGLE *dn, SINGLE *dnm1, SINGLE *
	dnm2);
 
 int slasr_(char *side, char *pivot, char *direct, INTEGER *m,
	 INTEGER *n, SINGLE *c__, SINGLE *s, SINGLE *a, INTEGER *lda);
 
 int slasrt_(char *id, INTEGER *n, SINGLE *d__, INTEGER *info);
 
 int slassq_(INTEGER *n, SINGLE *x, INTEGER *incx, SINGLE *scale, 
	SINGLE *sumsq);
 
 int slasv2_(SINGLE *f, SINGLE *g, SINGLE *h__, SINGLE *ssmin, SINGLE *
	ssmax, SINGLE *snr, SINGLE *csr, SINGLE *snl, SINGLE *csl);
 
 int slaswp_(INTEGER *n, SINGLE *a, INTEGER *lda, INTEGER *k1, 
	INTEGER *k2, INTEGER *ipiv, INTEGER *incx);
 
 int slasy2_(LOGICAL *ltranl, LOGICAL *ltranr, INTEGER *isgn, 
	INTEGER *n1, INTEGER *n2, SINGLE *tl, INTEGER *ldtl, SINGLE *tr, INTEGER *
	ldtr, SINGLE *b, INTEGER *ldb, SINGLE *scale, SINGLE *x, INTEGER *ldx, SINGLE 
	*xnorm, INTEGER *info);
 
 int slasyf_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 SINGLE *a, INTEGER *lda, INTEGER *ipiv, SINGLE *w, INTEGER *ldw, INTEGER 
	*info);
 
 int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, INTEGER *kd, SINGLE *ab, INTEGER *ldab, SINGLE *x, 
	SINGLE *scale, SINGLE *cnorm, INTEGER *info);
 
 int slatdf_(INTEGER *ijob, INTEGER *n, SINGLE *z__, INTEGER *
	ldz, SINGLE *rhs, SINGLE *rdsum, SINGLE *rdscal, INTEGER *ipiv, INTEGER *
	jpiv);
 
 int slatps_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, SINGLE *ap, SINGLE *x, SINGLE *scale, SINGLE *cnorm, 
	INTEGER *info);
 
 int slatrd_(char *uplo, INTEGER *n, INTEGER *nb, SINGLE *a, 
	INTEGER *lda, SINGLE *e, SINGLE *tau, SINGLE *w, INTEGER *ldw);
 
 int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *x, SINGLE *scale, SINGLE 
	*cnorm, INTEGER *info);
 
 int slatrz_(INTEGER *m, INTEGER *n, INTEGER *l, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work);
 
 int slatzm_(char *side, INTEGER *m, INTEGER *n, SINGLE *v, 
	INTEGER *incv, SINGLE *tau, SINGLE *c1, SINGLE *c2, INTEGER *ldc, SINGLE *
	work);
 
 int slauu2_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int slauum_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int sopgtr_(char *uplo, INTEGER *n, SINGLE *ap, SINGLE *tau, 
	SINGLE *q, INTEGER *ldq, SINGLE *work, INTEGER *info);
 
 int sopmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, SINGLE *ap, SINGLE *tau, SINGLE *c__, INTEGER *ldc, SINGLE *work, 
	INTEGER *info);
 
 int sorg2l_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sorg2r_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sorgbr_(char *vect, INTEGER *m, INTEGER *n, INTEGER *k, 
	SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER 
	*info);
 
 int sorghr_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorgl2_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sorglq_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorgql_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorgqr_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorgr2_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *info);
 
 int sorgrq_(INTEGER *m, INTEGER *n, INTEGER *k, SINGLE *a, 
	INTEGER *lda, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorgtr_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorm2l_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *info);
 
 int sorm2r_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *info);
 
 int sormbr_(char *vect, char *side, char *trans, INTEGER *m, 
	INTEGER *n, INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, 
	INTEGER *ldc, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormhr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *
	c__, INTEGER *ldc, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sorml2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *info);
 
 int sormlq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormql_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormqr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormr2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *info);
 
 int sormr3_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, 
	INTEGER *ldc, SINGLE *work, INTEGER *info);
 
 int sormrq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormrz_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, 
	INTEGER *ldc, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int sormtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *tau, SINGLE *c__, INTEGER *ldc,
	 SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int spbcon_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, SINGLE *anorm, SINGLE *rcond, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int spbequ_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, SINGLE *s, SINGLE *scond, SINGLE *amax, INTEGER *info);
 
 int spbrfs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *afb, INTEGER *ldafb, SINGLE *b, 
	INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int spbstf_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, INTEGER *info);
 
 int spbsv_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int spbsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *kd, 
	INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *afb, INTEGER *ldafb, 
	char *equed, SINGLE *s, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, 
	SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int spbtf2_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, INTEGER *info);
 
 int spbtrf_(char *uplo, INTEGER *n, INTEGER *kd, SINGLE *ab, 
	INTEGER *ldab, INTEGER *info);
 
 int spbtrs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int spocon_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *anorm, SINGLE *rcond, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int spoequ_(INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *s, SINGLE 
	*scond, SINGLE *amax, INTEGER *info);
 
 int sporfs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *af, INTEGER *ldaf, SINGLE *b, INTEGER *ldb, SINGLE *x,
	 INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int sposv_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sposvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, SINGLE *a, INTEGER *lda, SINGLE *af, INTEGER *ldaf, char *equed, 
	SINGLE *s, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *rcond, 
	SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int spotf2_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int spotrf_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int spotri_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *info);
 
 int spotrs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sppcon_(char *uplo, INTEGER *n, SINGLE *ap, SINGLE *anorm, 
	SINGLE *rcond, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sppequ_(char *uplo, INTEGER *n, SINGLE *ap, SINGLE *s, SINGLE *
	scond, SINGLE *amax, INTEGER *info);
 
 int spprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	SINGLE *afp, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *ferr, 
	SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sppsv_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sppsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, SINGLE *ap, SINGLE *afp, char *equed, SINGLE *s, SINGLE *b, INTEGER *
	ldb, SINGLE *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE 
	*work, INTEGER *iwork, INTEGER *info);
 
 int spptrf_(char *uplo, INTEGER *n, SINGLE *ap, INTEGER *info);
 
 int spptri_(char *uplo, INTEGER *n, SINGLE *ap, INTEGER *info);
 
 int spptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sptcon_(INTEGER *n, SINGLE *d__, SINGLE *e, SINGLE *anorm, 
	SINGLE *rcond, SINGLE *work, INTEGER *info);
 
 int spteqr_(char *compz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int sptrfs_(INTEGER *n, INTEGER *nrhs, SINGLE *d__, SINGLE *e, 
	SINGLE *df, SINGLE *ef, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, 
	SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *info);
 
 int sptsv_(INTEGER *n, INTEGER *nrhs, SINGLE *d__, SINGLE *e, 
	SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sptsvx_(char *fact, INTEGER *n, INTEGER *nrhs, SINGLE *d__,
	 SINGLE *e, SINGLE *df, SINGLE *ef, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER 
	*ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *info);
 
 int spttrf_(INTEGER *n, SINGLE *d__, SINGLE *e, INTEGER *info);
 
 int spttrs_(INTEGER *n, INTEGER *nrhs, SINGLE *d__, SINGLE *e, 
	SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sptts2_(INTEGER *n, INTEGER *nrhs, SINGLE *d__, SINGLE *e, 
	SINGLE *b, INTEGER *ldb);
 
 int srscl_(INTEGER *n, SINGLE *sa, SINGLE *sx, INTEGER *incx);
 
 int ssbev_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	SINGLE *ab, INTEGER *ldab, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work,
	 INTEGER *info);
 
 int ssbevd_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	SINGLE *ab, INTEGER *ldab, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work,
	 INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int ssbevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *kd, SINGLE *ab, INTEGER *ldab, SINGLE *q, INTEGER *ldq, SINGLE *vl,
	 SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, INTEGER *m, SINGLE *
	w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *iwork, INTEGER *
	ifail, INTEGER *info);
 
 int ssbgst_(char *vect, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, SINGLE *ab, INTEGER *ldab, SINGLE *bb, INTEGER *ldbb, SINGLE *
	x, INTEGER *ldx, SINGLE *work, INTEGER *info);
 
 int ssbgv_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, SINGLE *ab, INTEGER *ldab, SINGLE *bb, INTEGER *ldbb, SINGLE *
	w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int ssbgvd_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, SINGLE *ab, INTEGER *ldab, SINGLE *bb, INTEGER *ldbb, SINGLE *
	w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int ssbgvx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *ka, INTEGER *kb, SINGLE *ab, INTEGER *ldab, SINGLE *bb, INTEGER *
	ldbb, SINGLE *q, INTEGER *ldq, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER 
	*iu, SINGLE *abstol, INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE 
	*work, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int ssbtrd_(char *vect, char *uplo, INTEGER *n, INTEGER *kd, 
	SINGLE *ab, INTEGER *ldab, SINGLE *d__, SINGLE *e, SINGLE *q, INTEGER *ldq, 
	SINGLE *work, INTEGER *info);
 
 int sspcon_(char *uplo, INTEGER *n, SINGLE *ap, INTEGER *ipiv, 
	SINGLE *anorm, SINGLE *rcond, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sspev_(char *jobz, char *uplo, INTEGER *n, SINGLE *ap, 
	SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int sspevd_(char *jobz, char *uplo, INTEGER *n, SINGLE *ap, 
	SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER 
	*iwork, INTEGER *liwork, INTEGER *info);
 
 int sspevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	SINGLE *ap, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, 
	INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *
	iwork, INTEGER *ifail, INTEGER *info);
 
 int sspgst_(INTEGER *itype, char *uplo, INTEGER *n, SINGLE *ap,
	 SINGLE *bp, INTEGER *info);
 
 int sspgv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, SINGLE *ap, SINGLE *bp, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, 
	INTEGER *info);
 
 int sspgvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, SINGLE *ap, SINGLE *bp, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int sspgvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, SINGLE *ap, SINGLE *bp, SINGLE *vl, SINGLE *vu, INTEGER *il,
	 INTEGER *iu, SINGLE *abstol, INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *
	ldz, SINGLE *work, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int ssprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	SINGLE *afp, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *
	ldx, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *
	info);
 
 int sspsv_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	INTEGER *ipiv, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sspsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, SINGLE *ap, SINGLE *afp, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE 
	*x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr, SINGLE *berr, SINGLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int ssptrd_(char *uplo, INTEGER *n, SINGLE *ap, SINGLE *d__, 
	SINGLE *e, SINGLE *tau, INTEGER *info);
 
 int ssptrf_(char *uplo, INTEGER *n, SINGLE *ap, INTEGER *ipiv, 
	INTEGER *info);
 
 int ssptri_(char *uplo, INTEGER *n, SINGLE *ap, INTEGER *ipiv, 
	SINGLE *work, INTEGER *info);
 
 int ssptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *ap, 
	INTEGER *ipiv, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int sstebz_(char *range, char *order, INTEGER *n, SINGLE *vl, 
	SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, SINGLE *d__, SINGLE *e, 
	INTEGER *m, INTEGER *nsplit, SINGLE *w, INTEGER *iblock, INTEGER *
	isplit, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int sstedc_(char *compz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int sstegr_(char *jobz, char *range, INTEGER *n, SINGLE *d__, 
	SINGLE *e, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, 
	INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, INTEGER *isuppz, SINGLE *
	work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int sstein_(INTEGER *n, SINGLE *d__, SINGLE *e, INTEGER *m, SINGLE 
	*w, INTEGER *iblock, INTEGER *isplit, SINGLE *z__, INTEGER *ldz, SINGLE *
	work, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int ssteqr_(char *compz, INTEGER *n, SINGLE *d__, SINGLE *e, 
	SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int ssterf_(INTEGER *n, SINGLE *d__, SINGLE *e, INTEGER *info);
 
 int sstev_(char *jobz, INTEGER *n, SINGLE *d__, SINGLE *e, SINGLE *
	z__, INTEGER *ldz, SINGLE *work, INTEGER *info);
 
 int sstevd_(char *jobz, INTEGER *n, SINGLE *d__, SINGLE *e, SINGLE 
	*z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int sstevr_(char *jobz, char *range, INTEGER *n, SINGLE *d__, 
	SINGLE *e, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, 
	INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, INTEGER *isuppz, SINGLE *
	work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int sstevx_(char *jobz, char *range, INTEGER *n, SINGLE *d__, 
	SINGLE *e, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, 
	INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *
	iwork, INTEGER *ifail, INTEGER *info);
 
 int ssycon_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, SINGLE *anorm, SINGLE *rcond, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int ssyev_(char *jobz, char *uplo, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *w, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int ssyevd_(char *jobz, char *uplo, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *w, SINGLE *work, INTEGER *lwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int ssyevr_(char *jobz, char *range, char *uplo, INTEGER *n, 
	SINGLE *a, INTEGER *lda, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, 
	SINGLE *abstol, INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, INTEGER *
	isuppz, SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, 
	INTEGER *info);
 
 int ssyevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	SINGLE *a, INTEGER *lda, SINGLE *vl, SINGLE *vu, INTEGER *il, INTEGER *iu, 
	SINGLE *abstol, INTEGER *m, SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *
	work, INTEGER *lwork, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int ssygs2_(INTEGER *itype, char *uplo, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int ssygst_(INTEGER *itype, char *uplo, INTEGER *n, SINGLE *a, 
	INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int ssygv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *w, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int ssygvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *w, SINGLE *work, 
	INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int ssygvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *
	vl, SINGLE *vu, INTEGER *il, INTEGER *iu, SINGLE *abstol, INTEGER *m, 
	SINGLE *w, SINGLE *z__, INTEGER *ldz, SINGLE *work, INTEGER *lwork, INTEGER 
	*iwork, INTEGER *ifail, INTEGER *info);
 
 int ssyrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, SINGLE *af, INTEGER *ldaf, INTEGER *ipiv, SINGLE *b, 
	INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *
	work, INTEGER *iwork, INTEGER *info);
 
 int ssysv_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int ssysvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, SINGLE *a, INTEGER *lda, SINGLE *af, INTEGER *ldaf, INTEGER *ipiv, 
	SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx, SINGLE *rcond, SINGLE *ferr,
	 SINGLE *berr, SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *
	info);
 
 int ssytd2_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *d__, SINGLE *e, SINGLE *tau, INTEGER *info);
 
 int ssytf2_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, INTEGER *info);
 
 int ssytrd_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *d__, SINGLE *e, SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *
	info);
 
 int ssytrf_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int ssytri_(char *uplo, INTEGER *n, SINGLE *a, INTEGER *lda, 
	INTEGER *ipiv, SINGLE *work, INTEGER *info);
 
 int ssytrs_(char *uplo, INTEGER *n, INTEGER *nrhs, SINGLE *a, 
	INTEGER *lda, INTEGER *ipiv, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int stbcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	INTEGER *kd, SINGLE *ab, INTEGER *ldab, SINGLE *rcond, SINGLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int stbrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *b, INTEGER 
	*ldb, SINGLE *x, INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *work, 
	INTEGER *iwork, INTEGER *info);
 
 int stbtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, SINGLE *ab, INTEGER *ldab, SINGLE *b, INTEGER 
	*ldb, INTEGER *info);
 
 int stgevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *vl, 
	INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, INTEGER *mm, INTEGER *m, SINGLE 
	*work, INTEGER *info);
 
 int stgex2_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, SINGLE 
	*a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *q, INTEGER *ldq, SINGLE *
	z__, INTEGER *ldz, INTEGER *j1, INTEGER *n1, INTEGER *n2, SINGLE *work, 
	INTEGER *lwork, INTEGER *info);
 
 int stgexc_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, SINGLE 
	*a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *q, INTEGER *ldq, SINGLE *
	z__, INTEGER *ldz, INTEGER *ifst, INTEGER *ilst, SINGLE *work, INTEGER *
	lwork, INTEGER *info);
 
 int stgsen_(INTEGER *ijob, LOGICAL *wantq, LOGICAL *wantz, 
	LOGICAL *select, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *
	ldb, SINGLE *alphar, SINGLE *alphai, SINGLE *beta, SINGLE *q, INTEGER *ldq, 
	SINGLE *z__, INTEGER *ldz, INTEGER *m, SINGLE *pl, SINGLE *pr, SINGLE *dif, 
	SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, INTEGER *
	info);
 
 int stgsja_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, INTEGER *k, INTEGER *l, SINGLE *a, INTEGER *lda,
	 SINGLE *b, INTEGER *ldb, SINGLE *tola, SINGLE *tolb, SINGLE *alpha, SINGLE *
	beta, SINGLE *u, INTEGER *ldu, SINGLE *v, INTEGER *ldv, SINGLE *q, INTEGER *
	ldq, SINGLE *work, INTEGER *ncycle, INTEGER *info);
 
 int stgsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *vl, 
	INTEGER *ldvl, SINGLE *vr, INTEGER *ldvr, SINGLE *s, SINGLE *dif, INTEGER *
	mm, INTEGER *m, SINGLE *work, INTEGER *lwork, INTEGER *iwork, INTEGER *
	info);
 
 int stgsy2_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *c__, INTEGER *
	ldc, SINGLE *d__, INTEGER *ldd, SINGLE *e, INTEGER *lde, SINGLE *f, INTEGER 
	*ldf, SINGLE *scale, SINGLE *rdsum, SINGLE *rdscal, INTEGER *iwork, INTEGER 
	*pq, INTEGER *info);
 
 int stgsyl_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *c__, INTEGER *
	ldc, SINGLE *d__, INTEGER *ldd, SINGLE *e, INTEGER *lde, SINGLE *f, INTEGER 
	*ldf, SINGLE *scale, SINGLE *dif, SINGLE *work, INTEGER *lwork, INTEGER *
	iwork, INTEGER *info);
 
 int stpcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	SINGLE *ap, SINGLE *rcond, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int stprfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, SINGLE *ap, SINGLE *b, INTEGER *ldb, SINGLE *x, INTEGER *ldx,
	 SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, INTEGER *info);
 
 int stptri_(char *uplo, char *diag, INTEGER *n, SINGLE *ap, 
	INTEGER *info);
 
 int stptrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, SINGLE *ap, SINGLE *b, INTEGER *ldb, INTEGER *info);
 
 int strcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	SINGLE *a, INTEGER *lda, SINGLE *rcond, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int strevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, SINGLE *t, INTEGER *ldt, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, 
	INTEGER *ldvr, INTEGER *mm, INTEGER *m, SINGLE *work, INTEGER *info);
 
 int strexc_(char *compq, INTEGER *n, SINGLE *t, INTEGER *ldt, 
	SINGLE *q, INTEGER *ldq, INTEGER *ifst, INTEGER *ilst, SINGLE *work, 
	INTEGER *info);
 
 int strrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *x, 
	INTEGER *ldx, SINGLE *ferr, SINGLE *berr, SINGLE *work, INTEGER *iwork, 
	INTEGER *info);
 
 int strsen_(char *job, char *compq, LOGICAL *select, INTEGER 
	*n, SINGLE *t, INTEGER *ldt, SINGLE *q, INTEGER *ldq, SINGLE *wr, SINGLE *wi, 
	INTEGER *m, SINGLE *s, SINGLE *sep, SINGLE *work, INTEGER *lwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int strsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, SINGLE *t, INTEGER *ldt, SINGLE *vl, INTEGER *ldvl, SINGLE *vr, 
	INTEGER *ldvr, SINGLE *s, SINGLE *sep, INTEGER *mm, INTEGER *m, SINGLE *
	work, INTEGER *ldwork, INTEGER *iwork, INTEGER *info);
 
 int strsyl_(char *trana, char *tranb, INTEGER *isgn, INTEGER 
	*m, INTEGER *n, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, SINGLE *
	c__, INTEGER *ldc, SINGLE *scale, INTEGER *info);
 
 int strti2_(char *uplo, char *diag, INTEGER *n, SINGLE *a, 
	INTEGER *lda, INTEGER *info);
 
 int strtri_(char *uplo, char *diag, INTEGER *n, SINGLE *a, 
	INTEGER *lda, INTEGER *info);
 
 int strtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, SINGLE *a, INTEGER *lda, SINGLE *b, INTEGER *ldb, INTEGER *
	info);
 
 int stzrqf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, INTEGER *info);
 
 int stzrzf_(INTEGER *m, INTEGER *n, SINGLE *a, INTEGER *lda, 
	SINGLE *tau, SINGLE *work, INTEGER *lwork, INTEGER *info);
 
 int xerbla_(char *srname, INTEGER *info);
 
 int zbdsqr_(char *uplo, INTEGER *n, INTEGER *ncvt, INTEGER *
	nru, INTEGER *ncc, DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *vt, 
	INTEGER *ldvt, DOUBLECOMPLEX *u, INTEGER *ldu, DOUBLECOMPLEX *c__, 
	INTEGER *ldc, DOUBLE *rwork, INTEGER *info);
 
 int zdrot_(INTEGER *n, DOUBLECOMPLEX *cx, INTEGER *incx, 
	DOUBLECOMPLEX *cy, INTEGER *incy, DOUBLE *c__, DOUBLE *s);
 
 int zdrscl_(INTEGER *n, DOUBLE *sa, DOUBLECOMPLEX *sx, 
	INTEGER *incx);
 
 int zgbbrd_(char *vect, INTEGER *m, INTEGER *n, INTEGER *ncc,
	 INTEGER *kl, INTEGER *ku, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *q, INTEGER *ldq, 
	DOUBLECOMPLEX *pt, INTEGER *ldpt, DOUBLECOMPLEX *c__, INTEGER *ldc, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zgbcon_(char *norm, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, DOUBLE *anorm, 
	DOUBLE *rcond, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zgbequ_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *r__, DOUBLE *c__, 
	DOUBLE *rowcnd, DOUBLE *colcnd, DOUBLE *amax, INTEGER *
	info);
 
 int zgbrfs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *
	afb, INTEGER *ldafb, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zgbsv_(INTEGER *n, INTEGER *kl, INTEGER *ku, INTEGER *
	nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, DOUBLECOMPLEX *
	b, INTEGER *ldb, INTEGER *info);
 
 int zgbsvx_(char *fact, char *trans, INTEGER *n, INTEGER *kl,
	 INTEGER *ku, INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLECOMPLEX *afb, INTEGER *ldafb, INTEGER *ipiv, char *equed, 
	DOUBLE *r__, DOUBLE *c__, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zgbtf2_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int zgbtrf_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, INTEGER *info);
 
 int zgbtrs_(char *trans, INTEGER *n, INTEGER *kl, INTEGER *
	ku, INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *ipiv, 
	DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int zgebak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, DOUBLE *scale, INTEGER *m, DOUBLECOMPLEX *v, 
	INTEGER *ldv, INTEGER *info);
 
 int zgebal_(char *job, INTEGER *n, DOUBLECOMPLEX *a, INTEGER 
	*lda, INTEGER *ilo, INTEGER *ihi, DOUBLE *scale, INTEGER *info);
 
 int zgebd2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *tauq, 
	DOUBLECOMPLEX *taup, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgebrd_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *tauq, 
	DOUBLECOMPLEX *taup, DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *
	info);
 
 int zgecon_(char *norm, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zgeequ_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *r__, DOUBLE *c__, DOUBLE *rowcnd, 
	DOUBLE *colcnd, DOUBLE *amax, INTEGER *info);
 
 int zgees_(char *jobvs, char *sort, __CLPK_L_fp select, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *sdim, DOUBLECOMPLEX *w, 
	DOUBLECOMPLEX *vs, INTEGER *ldvs, DOUBLECOMPLEX *work, INTEGER *lwork,
	 DOUBLE *rwork, LOGICAL *bwork, INTEGER *info);
 
 int zgeesx_(char *jobvs, char *sort, __CLPK_L_fp select, char *
	sense, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *sdim, 
	DOUBLECOMPLEX *w, DOUBLECOMPLEX *vs, INTEGER *ldvs, DOUBLE *
	rconde, DOUBLE *rcondv, DOUBLECOMPLEX *work, INTEGER *lwork, 
	DOUBLE *rwork, LOGICAL *bwork, INTEGER *info);
 
 int zgeev_(char *jobvl, char *jobvr, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *w, DOUBLECOMPLEX *vl, 
	INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, DOUBLECOMPLEX *work, 
	INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *w, 
	DOUBLECOMPLEX *vl, INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *scale, DOUBLE *abnrm, 
	DOUBLE *rconde, DOUBLE *rcondv, DOUBLECOMPLEX *work, INTEGER *
	lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgegs_(char *jobvsl, char *jobvsr, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *beta, DOUBLECOMPLEX *vsl, 
	INTEGER *ldvsl, DOUBLECOMPLEX *vsr, INTEGER *ldvsr, DOUBLECOMPLEX *
	work, INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgegv_(char *jobvl, char *jobvr, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *beta, DOUBLECOMPLEX *vl, INTEGER 
	*ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, DOUBLECOMPLEX *work, INTEGER 
	*lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgehd2_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zgehrd_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zgelq2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgelqf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zgels_(char *trans, INTEGER *m, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zgelsx_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *jpvt, DOUBLE *rcond, INTEGER *rank, DOUBLECOMPLEX *work, 
	DOUBLE *rwork, INTEGER *info);
 
 int zgelsy_(INTEGER *m, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *jpvt, DOUBLE *rcond, INTEGER *rank, DOUBLECOMPLEX *work, 
	INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgeql2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgeqlf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zgeqp3_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *jpvt, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, 
	INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zgeqpf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *jpvt, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, 
	DOUBLE *rwork, INTEGER *info);
 
 int zgeqr2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgeqrf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zgerfs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *ldaf, 
	INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, 
	INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work,
	 DOUBLE *rwork, INTEGER *info);
 
 int zgerq2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgerqf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zgesc2_(INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLECOMPLEX *rhs, INTEGER *ipiv, INTEGER *jpiv, DOUBLE *scale);
 
 int zgesv_(INTEGER *n, INTEGER *nrhs, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *
	info);
 
 int zgesvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *
	ldaf, INTEGER *ipiv, char *equed, DOUBLE *r__, DOUBLE *c__, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zgetc2_(INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	INTEGER *ipiv, INTEGER *jpiv, INTEGER *info);
 
 int zgetf2_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, INTEGER *info);
 
 int zgetrf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, INTEGER *info);
 
 int zgetri_(INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zgetrs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, 
	INTEGER *ldb, INTEGER *info);
 
 int zggbak_(char *job, char *side, INTEGER *n, INTEGER *ilo, 
	INTEGER *ihi, DOUBLE *lscale, DOUBLE *rscale, INTEGER *m, 
	DOUBLECOMPLEX *v, INTEGER *ldv, INTEGER *info);
 
 int zggbal_(char *job, INTEGER *n, DOUBLECOMPLEX *a, INTEGER 
	*lda, DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *ilo, INTEGER *ihi, 
	DOUBLE *lscale, DOUBLE *rscale, DOUBLE *work, INTEGER *
	info);
 
 int zgges_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	delctg, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, INTEGER *sdim, DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *
	beta, DOUBLECOMPLEX *vsl, INTEGER *ldvsl, DOUBLECOMPLEX *vsr, INTEGER 
	*ldvsr, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, 
	LOGICAL *bwork, INTEGER *info);
 
 int zggesx_(char *jobvsl, char *jobvsr, char *sort, __CLPK_L_fp 
	delctg, char *sense, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *sdim, DOUBLECOMPLEX *alpha, 
	DOUBLECOMPLEX *beta, DOUBLECOMPLEX *vsl, INTEGER *ldvsl, 
	DOUBLECOMPLEX *vsr, INTEGER *ldvsr, DOUBLE *rconde, DOUBLE *
	rcondv, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, 
	INTEGER *iwork, INTEGER *liwork, LOGICAL *bwork, INTEGER *info);
 
 int zggev_(char *jobvl, char *jobvr, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *beta, DOUBLECOMPLEX *vl, INTEGER 
	*ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, DOUBLECOMPLEX *work, INTEGER 
	*lwork, DOUBLE *rwork, INTEGER *info);
 
 int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *beta, 
	DOUBLECOMPLEX *vl, INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, 
	INTEGER *ilo, INTEGER *ihi, DOUBLE *lscale, DOUBLE *rscale, 
	DOUBLE *abnrm, DOUBLE *bbnrm, DOUBLE *rconde, DOUBLE *
	rcondv, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, 
	INTEGER *iwork, LOGICAL *bwork, INTEGER *info);
 
 int zggglm_(INTEGER *n, INTEGER *m, INTEGER *p, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *d__, DOUBLECOMPLEX *x, DOUBLECOMPLEX *y, DOUBLECOMPLEX 
	*work, INTEGER *lwork, INTEGER *info);
 
 int zgghrd_(char *compq, char *compz, INTEGER *n, INTEGER *
	ilo, INTEGER *ihi, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *z__, 
	INTEGER *ldz, INTEGER *info);
 
 int zgglse_(INTEGER *m, INTEGER *n, INTEGER *p, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *c__, DOUBLECOMPLEX *d__, DOUBLECOMPLEX *x, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zggqrf_(INTEGER *n, INTEGER *m, INTEGER *p, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *taua, DOUBLECOMPLEX *b,
	 INTEGER *ldb, DOUBLECOMPLEX *taub, DOUBLECOMPLEX *work, INTEGER *
	lwork, INTEGER *info);
 
 int zggrqf_(INTEGER *m, INTEGER *p, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *taua, DOUBLECOMPLEX *b,
	 INTEGER *ldb, DOUBLECOMPLEX *taub, DOUBLECOMPLEX *work, INTEGER *
	lwork, INTEGER *info);
 
 int zggsvd_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *n, INTEGER *p, INTEGER *k, INTEGER *l, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLE *alpha, 
	DOUBLE *beta, DOUBLECOMPLEX *u, INTEGER *ldu, DOUBLECOMPLEX *v, 
	INTEGER *ldv, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *work, 
	DOUBLE *rwork, INTEGER *iwork, INTEGER *info);
 
 int zggsvp_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX 
	*b, INTEGER *ldb, DOUBLE *tola, DOUBLE *tolb, INTEGER *k, 
	INTEGER *l, DOUBLECOMPLEX *u, INTEGER *ldu, DOUBLECOMPLEX *v, INTEGER 
	*ldv, DOUBLECOMPLEX *q, INTEGER *ldq, INTEGER *iwork, DOUBLE *
	rwork, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zgtcon_(char *norm, INTEGER *n, DOUBLECOMPLEX *dl, 
	DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, DOUBLECOMPLEX *du2, INTEGER *
	ipiv, DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *work, 
	INTEGER *info);
 
 int zgtrfs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *dl, DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, 
	DOUBLECOMPLEX *dlf, DOUBLECOMPLEX *df, DOUBLECOMPLEX *duf, 
	DOUBLECOMPLEX *du2, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zgtsv_(INTEGER *n, INTEGER *nrhs, DOUBLECOMPLEX *dl, 
	DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, DOUBLECOMPLEX *b, INTEGER *ldb,
	 INTEGER *info);
 
 int zgtsvx_(char *fact, char *trans, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *dl, DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, 
	DOUBLECOMPLEX *dlf, DOUBLECOMPLEX *df, DOUBLECOMPLEX *duf, 
	DOUBLECOMPLEX *du2, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zgttrf_(INTEGER *n, DOUBLECOMPLEX *dl, DOUBLECOMPLEX *
	d__, DOUBLECOMPLEX *du, DOUBLECOMPLEX *du2, INTEGER *ipiv, INTEGER *
	info);
 
 int zgttrs_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *dl, DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, 
	DOUBLECOMPLEX *du2, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zgtts2_(INTEGER *itrans, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *dl, DOUBLECOMPLEX *d__, DOUBLECOMPLEX *du, 
	DOUBLECOMPLEX *du2, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb);
 
 int zhbev_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *w, DOUBLECOMPLEX *z__, 
	INTEGER *ldz, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zhbevd_(char *jobz, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *w, DOUBLECOMPLEX *z__, 
	INTEGER *ldz, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, 
	INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int zhbevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *kd, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *q, 
	INTEGER *ldq, DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *
	iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, DOUBLECOMPLEX *z__,
	 INTEGER *ldz, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *iwork,
	 INTEGER *ifail, INTEGER *info);
 
 int zhbgst_(char *vect, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *bb, 
	INTEGER *ldbb, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLECOMPLEX *work, 
	DOUBLE *rwork, INTEGER *info);
 
 int zhbgv_(char *jobz, char *uplo, INTEGER *n, INTEGER *ka, 
	INTEGER *kb, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *bb, 
	INTEGER *ldbb, DOUBLE *w, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zhbgvx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	INTEGER *ka, INTEGER *kb, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLECOMPLEX *bb, INTEGER *ldbb, DOUBLECOMPLEX *q, INTEGER *ldq, 
	DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *iu, DOUBLE *
	abstol, INTEGER *m, DOUBLE *w, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *iwork, INTEGER *
	ifail, INTEGER *info);
 
 int zhbtrd_(char *vect, char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *d__, DOUBLE *e, 
	DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zhecon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLE *anorm, DOUBLE *rcond, 
	DOUBLECOMPLEX *work, INTEGER *info);
 
 int zheev_(char *jobz, char *uplo, INTEGER *n, DOUBLECOMPLEX 
	*a, INTEGER *lda, DOUBLE *w, DOUBLECOMPLEX *work, INTEGER *lwork, 
	DOUBLE *rwork, INTEGER *info);
 
 int zheevd_(char *jobz, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *w, DOUBLECOMPLEX *work, 
	INTEGER *lwork, DOUBLE *rwork, INTEGER *lrwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int zheevr_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *vl, DOUBLE *vu, 
	INTEGER *il, INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *
	w, DOUBLECOMPLEX *z__, INTEGER *ldz, INTEGER *isuppz, DOUBLECOMPLEX *
	work, INTEGER *lwork, DOUBLE *rwork, INTEGER *lrwork, INTEGER *
	iwork, INTEGER *liwork, INTEGER *info);
 
 int zheevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *vl, DOUBLE *vu, 
	INTEGER *il, INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *
	w, DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLECOMPLEX *work, INTEGER *
	lwork, DOUBLE *rwork, INTEGER *iwork, INTEGER *ifail, INTEGER *
	info);
 
 int zhegs2_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zhegst_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zhegv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLE *w, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork,
	 INTEGER *info);
 
 int zhegvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLE *w, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork,
	 INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int zhegvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLE *vl, DOUBLE *vu, INTEGER *il, INTEGER *
	iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, DOUBLECOMPLEX *z__,
	 INTEGER *ldz, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork,
	 INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int zherfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *ldaf, 
	INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, 
	INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work,
	 DOUBLE *rwork, INTEGER *info);
 
 int zhesv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zhesvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *
	ldaf, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x,
	 INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zhetf2_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, INTEGER *info);
 
 int zhetrd_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zhetrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *lwork, 
	INTEGER *info);
 
 int zhetri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zhetrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, 
	INTEGER *ldb, INTEGER *info);
 
 int zhgeqz_(char *job, char *compq, char *compz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *
	beta, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *z__, INTEGER *
	ldz, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, INTEGER *
	info);
 
 int zhpcon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zhpev_(char *jobz, char *uplo, INTEGER *n, DOUBLECOMPLEX 
	*ap, DOUBLE *w, DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zhpevd_(char *jobz, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *ap, DOUBLE *w, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, INTEGER *
	lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *info);
 
 int zhpevx_(char *jobz, char *range, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *ap, DOUBLE *vl, DOUBLE *vu, INTEGER *il, 
	INTEGER *iu, DOUBLE *abstol, INTEGER *m, DOUBLE *w, 
	DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLECOMPLEX *work, DOUBLE *
	rwork, INTEGER *iwork, INTEGER *ifail, INTEGER *info);
 
 int zhpgst_(INTEGER *itype, char *uplo, INTEGER *n, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *bp, INTEGER *info);
 
 int zhpgv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *bp, DOUBLE *w, DOUBLECOMPLEX 
	*z__, INTEGER *ldz, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zhpgvd_(INTEGER *itype, char *jobz, char *uplo, INTEGER *
	n, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *bp, DOUBLE *w, DOUBLECOMPLEX 
	*z__, INTEGER *ldz, DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *
	rwork, INTEGER *lrwork, INTEGER *iwork, INTEGER *liwork, INTEGER *
	info);
 
 int zhpgvx_(INTEGER *itype, char *jobz, char *range, char *
	uplo, INTEGER *n, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *bp, DOUBLE *
	vl, DOUBLE *vu, INTEGER *il, INTEGER *iu, DOUBLE *abstol, 
	INTEGER *m, DOUBLE *w, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *iwork, INTEGER *
	ifail, INTEGER *info);
 
 int zhprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, INTEGER *ipiv, DOUBLECOMPLEX *
	b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zhpsv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zhpsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, INTEGER *ipiv, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zhptrd_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLE *d__, DOUBLE *e, DOUBLECOMPLEX *tau, INTEGER *info);
 
 int zhptrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, INTEGER *info);
 
 int zhptri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zhptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zhsein_(char *side, char *eigsrc, char *initv, LOGICAL *
	select, INTEGER *n, DOUBLECOMPLEX *h__, INTEGER *ldh, DOUBLECOMPLEX *
	w, DOUBLECOMPLEX *vl, INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr,
	 INTEGER *mm, INTEGER *m, DOUBLECOMPLEX *work, DOUBLE *rwork, 
	INTEGER *ifaill, INTEGER *ifailr, INTEGER *info);
 
 int zhseqr_(char *job, char *compz, INTEGER *n, INTEGER *ilo,
	 INTEGER *ihi, DOUBLECOMPLEX *h__, INTEGER *ldh, DOUBLECOMPLEX *w, 
	DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zlabrd_(INTEGER *m, INTEGER *n, INTEGER *nb, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *d__, DOUBLE *e, 
	DOUBLECOMPLEX *tauq, DOUBLECOMPLEX *taup, DOUBLECOMPLEX *x, INTEGER *
	ldx, DOUBLECOMPLEX *y, INTEGER *ldy);
 
 int zlacgv_(INTEGER *n, DOUBLECOMPLEX *x, INTEGER *incx);
 
 int zlacon_(INTEGER *n, DOUBLECOMPLEX *v, DOUBLECOMPLEX *x, 
	DOUBLE *est, INTEGER *kase);
 
 int zlacp2_(char *uplo, INTEGER *m, INTEGER *n, DOUBLE *
	a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb);
 
 int zlacpy_(char *uplo, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb);
 
 int zlacrm_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *b, INTEGER *ldb, DOUBLECOMPLEX *c__, 
	INTEGER *ldc, DOUBLE *rwork);
 
 int zlacrt_(INTEGER *n, DOUBLECOMPLEX *cx, INTEGER *incx, 
	DOUBLECOMPLEX *cy, INTEGER *incy, DOUBLECOMPLEX *c__, DOUBLECOMPLEX *
	s);
 
 int zlaed0_(INTEGER *qsiz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *qstore, 
	INTEGER *ldqs, DOUBLE *rwork, INTEGER *iwork, INTEGER *info);
 
 int zlaed7_(INTEGER *n, INTEGER *cutpnt, INTEGER *qsiz, 
	INTEGER *tlvls, INTEGER *curlvl, INTEGER *curpbm, DOUBLE *d__, 
	DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLE *rho, INTEGER *indxq, 
	DOUBLE *qstore, INTEGER *qptr, INTEGER *prmptr, INTEGER *perm, 
	INTEGER *givptr, INTEGER *givcol, DOUBLE *givnum, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *iwork, INTEGER *info);
 
 int zlaed8_(INTEGER *k, INTEGER *n, INTEGER *qsiz, 
	DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLE *d__, DOUBLE *rho, 
	INTEGER *cutpnt, DOUBLE *z__, DOUBLE *dlamda, DOUBLECOMPLEX *
	q2, INTEGER *ldq2, DOUBLE *w, INTEGER *indxp, INTEGER *indx, 
	INTEGER *indxq, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, 
	DOUBLE *givnum, INTEGER *info);
 
 int zlaein_(LOGICAL *rightv, LOGICAL *noinit, INTEGER *n, 
	DOUBLECOMPLEX *h__, INTEGER *ldh, DOUBLECOMPLEX *w, DOUBLECOMPLEX *v, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLE *rwork, DOUBLE *eps3, 
	DOUBLE *smlnum, INTEGER *info);
 
 int zlaesy_(DOUBLECOMPLEX *a, DOUBLECOMPLEX *b, 
	DOUBLECOMPLEX *c__, DOUBLECOMPLEX *rt1, DOUBLECOMPLEX *rt2, 
	DOUBLECOMPLEX *evscal, DOUBLECOMPLEX *cs1, DOUBLECOMPLEX *sn1);
 
 int zlaev2_(DOUBLECOMPLEX *a, DOUBLECOMPLEX *b, 
	DOUBLECOMPLEX *c__, DOUBLE *rt1, DOUBLE *rt2, DOUBLE *cs1,
	 DOUBLECOMPLEX *sn1);
 
 int zlags2_(LOGICAL *upper, DOUBLE *a1, DOUBLECOMPLEX *
	a2, DOUBLE *a3, DOUBLE *b1, DOUBLECOMPLEX *b2, DOUBLE *b3,
	 DOUBLE *csu, DOUBLECOMPLEX *snu, DOUBLE *csv, DOUBLECOMPLEX *
	snv, DOUBLE *csq, DOUBLECOMPLEX *snq);
 
 int zlagtm_(char *trans, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *alpha, DOUBLECOMPLEX *dl, DOUBLECOMPLEX *d__, 
	DOUBLECOMPLEX *du, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *beta, 
	DOUBLECOMPLEX *b, INTEGER *ldb);
 
 int zlahef_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *w, 
	INTEGER *ldw, INTEGER *info);
 
 int zlahqr_(LOGICAL *wantt, LOGICAL *wantz, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLECOMPLEX *h__, INTEGER *ldh, 
	DOUBLECOMPLEX *w, INTEGER *iloz, INTEGER *ihiz, DOUBLECOMPLEX *z__, 
	INTEGER *ldz, INTEGER *info);
 
 int zlahrd_(INTEGER *n, INTEGER *k, INTEGER *nb, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *t, 
	INTEGER *ldt, DOUBLECOMPLEX *y, INTEGER *ldy);
 
 int zlaic1_(INTEGER *job, INTEGER *j, DOUBLECOMPLEX *x, 
	DOUBLE *sest, DOUBLECOMPLEX *w, DOUBLECOMPLEX *gamma, DOUBLE *
	sestpr, DOUBLECOMPLEX *s, DOUBLECOMPLEX *c__);
 
 int zlals0_(INTEGER *icompq, INTEGER *nl, INTEGER *nr, 
	INTEGER *sqre, INTEGER *nrhs, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *bx, INTEGER *ldbx, INTEGER *perm, INTEGER *givptr, 
	INTEGER *givcol, INTEGER *ldgcol, DOUBLE *givnum, INTEGER *ldgnum,
	 DOUBLE *poles, DOUBLE *difl, DOUBLE *difr, DOUBLE *
	z__, INTEGER *k, DOUBLE *c__, DOUBLE *s, DOUBLE *rwork, 
	INTEGER *info);
 
 int zlalsa_(INTEGER *icompq, INTEGER *smlsiz, INTEGER *n, 
	INTEGER *nrhs, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *bx, 
	INTEGER *ldbx, DOUBLE *u, INTEGER *ldu, DOUBLE *vt, INTEGER *
	k, DOUBLE *difl, DOUBLE *difr, DOUBLE *z__, DOUBLE *
	poles, INTEGER *givptr, INTEGER *givcol, INTEGER *ldgcol, INTEGER *
	perm, DOUBLE *givnum, DOUBLE *c__, DOUBLE *s, DOUBLE *
	rwork, INTEGER *iwork, INTEGER *info);
 
 int zlapll_(INTEGER *n, DOUBLECOMPLEX *x, INTEGER *incx, 
	DOUBLECOMPLEX *y, INTEGER *incy, DOUBLE *ssmin);
 
 int zlapmt_(LOGICAL *forwrd, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *x, INTEGER *ldx, INTEGER *k);
 
 int zlaqgb_(INTEGER *m, INTEGER *n, INTEGER *kl, INTEGER *ku,
	 DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *r__, DOUBLE *c__, 
	DOUBLE *rowcnd, DOUBLE *colcnd, DOUBLE *amax, char *equed);
 
 int zlaqge_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *r__, DOUBLE *c__, DOUBLE *rowcnd, 
	DOUBLE *colcnd, DOUBLE *amax, char *equed);
 
 int zlaqhb_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *s, DOUBLE *scond, 
	DOUBLE *amax, char *equed);
 
 int zlaqhe_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *s, DOUBLE *scond, DOUBLE *amax, 
	char *equed);
 
 int zlaqhp_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, char *equed);
 
 int zlaqp2_(INTEGER *m, INTEGER *n, INTEGER *offset, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *jpvt, DOUBLECOMPLEX *tau, 
	DOUBLE *vn1, DOUBLE *vn2, DOUBLECOMPLEX *work);
 
 int zlaqps_(INTEGER *m, INTEGER *n, INTEGER *offset, INTEGER 
	*nb, INTEGER *kb, DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *jpvt, 
	DOUBLECOMPLEX *tau, DOUBLE *vn1, DOUBLE *vn2, DOUBLECOMPLEX *
	auxv, DOUBLECOMPLEX *f, INTEGER *ldf);
 
 int zlaqsb_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *s, DOUBLE *scond, 
	DOUBLE *amax, char *equed);
 
 int zlaqsp_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, char *equed);
 
 int zlaqsy_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *s, DOUBLE *scond, DOUBLE *amax, 
	char *equed);
 
 int zlar1v_(INTEGER *n, INTEGER *b1, INTEGER *bn, DOUBLE 
	*sigma, DOUBLE *d__, DOUBLE *l, DOUBLE *ld, DOUBLE *
	lld, DOUBLE *gersch, DOUBLECOMPLEX *z__, DOUBLE *ztz, 
	DOUBLE *mingma, INTEGER *r__, INTEGER *isuppz, DOUBLE *work);
 
 int zlar2v_(INTEGER *n, DOUBLECOMPLEX *x, DOUBLECOMPLEX *y, 
	DOUBLECOMPLEX *z__, INTEGER *incx, DOUBLE *c__, DOUBLECOMPLEX *s, 
	INTEGER *incc);
 
 int zlarcm_(INTEGER *m, INTEGER *n, DOUBLE *a, INTEGER *
	lda, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *c__, INTEGER *ldc,
	 DOUBLE *rwork);
 
 int zlarf_(char *side, INTEGER *m, INTEGER *n, DOUBLECOMPLEX 
	*v, INTEGER *incv, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *c__, INTEGER *
	ldc, DOUBLECOMPLEX *work);
 
 int zlarfb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, DOUBLECOMPLEX *v, INTEGER 
	*ldv, DOUBLECOMPLEX *t, INTEGER *ldt, DOUBLECOMPLEX *c__, INTEGER *
	ldc, DOUBLECOMPLEX *work, INTEGER *ldwork);
 
 int zlarfg_(INTEGER *n, DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *
	x, INTEGER *incx, DOUBLECOMPLEX *tau);
 
 int zlarft_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, DOUBLECOMPLEX *v, INTEGER *ldv, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	t, INTEGER *ldt);
 
 int zlarfx_(char *side, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *v, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *c__, INTEGER *
	ldc, DOUBLECOMPLEX *work);
 
 int zlargv_(INTEGER *n, DOUBLECOMPLEX *x, INTEGER *incx, 
	DOUBLECOMPLEX *y, INTEGER *incy, DOUBLE *c__, INTEGER *incc);
 
 int zlarnv_(INTEGER *idist, INTEGER *iseed, INTEGER *n, 
	DOUBLECOMPLEX *x);
 
 int zlarrv_(INTEGER *n, DOUBLE *d__, DOUBLE *l, 
	INTEGER *isplit, INTEGER *m, DOUBLE *w, INTEGER *iblock, 
	DOUBLE *gersch, DOUBLE *tol, DOUBLECOMPLEX *z__, INTEGER *ldz,
	 INTEGER *isuppz, DOUBLE *work, INTEGER *iwork, INTEGER *info);
 
 int zlartg_(DOUBLECOMPLEX *f, DOUBLECOMPLEX *g, DOUBLE *
	cs, DOUBLECOMPLEX *sn, DOUBLECOMPLEX *r__);
 
 int zlartv_(INTEGER *n, DOUBLECOMPLEX *x, INTEGER *incx, 
	DOUBLECOMPLEX *y, INTEGER *incy, DOUBLE *c__, DOUBLECOMPLEX *s, 
	INTEGER *incc);
 
 int zlarz_(char *side, INTEGER *m, INTEGER *n, INTEGER *l, 
	DOUBLECOMPLEX *v, INTEGER *incv, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	c__, INTEGER *ldc, DOUBLECOMPLEX *work);
 
 int zlarzb_(char *side, char *trans, char *direct, char *
	storev, INTEGER *m, INTEGER *n, INTEGER *k, INTEGER *l, DOUBLECOMPLEX 
	*v, INTEGER *ldv, DOUBLECOMPLEX *t, INTEGER *ldt, DOUBLECOMPLEX *c__, 
	INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *ldwork);
 
 int zlarzt_(char *direct, char *storev, INTEGER *n, INTEGER *
	k, DOUBLECOMPLEX *v, INTEGER *ldv, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	t, INTEGER *ldt);
 
 int zlascl_(char *type__, INTEGER *kl, INTEGER *ku, 
	DOUBLE *cfrom, DOUBLE *cto, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *info);
 
 int zlaset_(char *uplo, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *beta, DOUBLECOMPLEX *a, INTEGER *
	lda);
 
 int zlasr_(char *side, char *pivot, char *direct, INTEGER *m,
	 INTEGER *n, DOUBLE *c__, DOUBLE *s, DOUBLECOMPLEX *a, 
	INTEGER *lda);
 
 int zlassq_(INTEGER *n, DOUBLECOMPLEX *x, INTEGER *incx, 
	DOUBLE *scale, DOUBLE *sumsq);
 
 int zlaswp_(INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	INTEGER *k1, INTEGER *k2, INTEGER *ipiv, INTEGER *incx);
 
 int zlasyf_(char *uplo, INTEGER *n, INTEGER *nb, INTEGER *kb,
	 DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *w, 
	INTEGER *ldw, INTEGER *info);
 
 int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, INTEGER *kd, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLECOMPLEX *x, DOUBLE *scale, DOUBLE *cnorm, INTEGER *info);
 
 int zlatdf_(INTEGER *ijob, INTEGER *n, DOUBLECOMPLEX *z__, 
	INTEGER *ldz, DOUBLECOMPLEX *rhs, DOUBLE *rdsum, DOUBLE *
	rdscal, INTEGER *ipiv, INTEGER *jpiv);
 
 int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *x, DOUBLE *
	scale, DOUBLE *cnorm, INTEGER *info);
 
 int zlatrd_(char *uplo, INTEGER *n, INTEGER *nb, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *e, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *w, INTEGER *ldw);
 
 int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *x, 
	DOUBLE *scale, DOUBLE *cnorm, INTEGER *info);
 
 int zlatrz_(INTEGER *m, INTEGER *n, INTEGER *l, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work);
 
 int zlatzm_(char *side, INTEGER *m, INTEGER *n, 
	DOUBLECOMPLEX *v, INTEGER *incv, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	c1, DOUBLECOMPLEX *c2, INTEGER *ldc, DOUBLECOMPLEX *work);
 
 int zlauu2_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int zlauum_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int zpbcon_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *anorm, DOUBLE *
	rcond, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zpbequ_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *s, DOUBLE *scond, 
	DOUBLE *amax, INTEGER *info);
 
 int zpbrfs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *afb, INTEGER *
	ldafb, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx,
	 DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *
	rwork, INTEGER *info);
 
 int zpbstf_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *info);
 
 int zpbsv_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *b, INTEGER *
	ldb, INTEGER *info);
 
 int zpbsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *kd, 
	INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *afb, 
	INTEGER *ldafb, char *equed, DOUBLE *s, DOUBLECOMPLEX *b, INTEGER 
	*ldb, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *
	ferr, DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, 
	INTEGER *info);
 
 int zpbtf2_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *info);
 
 int zpbtrf_(char *uplo, INTEGER *n, INTEGER *kd, 
	DOUBLECOMPLEX *ab, INTEGER *ldab, INTEGER *info);
 
 int zpbtrs_(char *uplo, INTEGER *n, INTEGER *kd, INTEGER *
	nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLECOMPLEX *b, INTEGER *
	ldb, INTEGER *info);
 
 int zpocon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zpoequ_(INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, INTEGER *info);
 
 int zporfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *ldaf, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *
	rwork, INTEGER *info);
 
 int zposv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zposvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *
	ldaf, char *equed, DOUBLE *s, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zpotf2_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int zpotrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int zpotri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *info);
 
 int zpotrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zppcon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *work, DOUBLE 
	*rwork, INTEGER *info);
 
 int zppequ_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLE *s, DOUBLE *scond, DOUBLE *amax, INTEGER *info);
 
 int zpprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, DOUBLECOMPLEX *b, INTEGER *ldb,
	 DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int zppsv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int zppsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, char *equed, DOUBLE *
	s, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zpptrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *info);
 
 int zpptri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *info);
 
 int zpptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int zptcon_(INTEGER *n, DOUBLE *d__, DOUBLECOMPLEX *e, 
	DOUBLE *anorm, DOUBLE *rcond, DOUBLE *rwork, INTEGER *
	info);
 
 int zptrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *d__, DOUBLECOMPLEX *e, DOUBLE *df, DOUBLECOMPLEX *ef, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *
	rwork, INTEGER *info);
 
 int zptsv_(INTEGER *n, INTEGER *nrhs, DOUBLE *d__, 
	DOUBLECOMPLEX *e, DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int zptsvx_(char *fact, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *d__, DOUBLECOMPLEX *e, DOUBLE *df, DOUBLECOMPLEX *ef, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zpttrf_(INTEGER *n, DOUBLE *d__, DOUBLECOMPLEX *e, 
	INTEGER *info);
 
 int zpttrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *d__, DOUBLECOMPLEX *e, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zptts2_(INTEGER *iuplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLE *d__, DOUBLECOMPLEX *e, DOUBLECOMPLEX *b, INTEGER *ldb);
 
 int zrot_(INTEGER *n, DOUBLECOMPLEX *cx, INTEGER *incx, 
	DOUBLECOMPLEX *cy, INTEGER *incy, DOUBLE *c__, DOUBLECOMPLEX *s);
 
 int zspcon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, DOUBLE *anorm, DOUBLE *rcond, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zspmv_(char *uplo, INTEGER *n, DOUBLECOMPLEX *alpha, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *x, INTEGER *incx, DOUBLECOMPLEX *
	beta, DOUBLECOMPLEX *y, INTEGER *incy);
 
 int zspr_(char *uplo, INTEGER *n, DOUBLECOMPLEX *alpha, 
	DOUBLECOMPLEX *x, INTEGER *incx, DOUBLECOMPLEX *ap);
 
 int zsprfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, INTEGER *ipiv, DOUBLECOMPLEX *
	b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int zspsv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zspsvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *afp, INTEGER *ipiv, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int zsptrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, INTEGER *info);
 
 int zsptri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zsptrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *ap, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int zstedc_(char *compz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLECOMPLEX *work, 
	INTEGER *lwork, DOUBLE *rwork, INTEGER *lrwork, INTEGER *iwork, 
	INTEGER *liwork, INTEGER *info);
 
 int zstein_(INTEGER *n, DOUBLE *d__, DOUBLE *e, 
	INTEGER *m, DOUBLE *w, INTEGER *iblock, INTEGER *isplit, 
	DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLE *work, INTEGER *iwork, 
	INTEGER *ifail, INTEGER *info);
 
 int zsteqr_(char *compz, INTEGER *n, DOUBLE *d__, 
	DOUBLE *e, DOUBLECOMPLEX *z__, INTEGER *ldz, DOUBLE *work, 
	INTEGER *info);
 
 int zsycon_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLE *anorm, DOUBLE *rcond, 
	DOUBLECOMPLEX *work, INTEGER *info);
 
 int zsymv_(char *uplo, INTEGER *n, DOUBLECOMPLEX *alpha, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *x, INTEGER *incx, 
	DOUBLECOMPLEX *beta, DOUBLECOMPLEX *y, INTEGER *incy);
 
 int zsyr_(char *uplo, INTEGER *n, DOUBLECOMPLEX *alpha, 
	DOUBLECOMPLEX *x, INTEGER *incx, DOUBLECOMPLEX *a, INTEGER *lda);
 
 int zsyrfs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *ldaf, 
	INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, 
	INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work,
	 DOUBLE *rwork, INTEGER *info);
 
 int zsysv_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int zsysvx_(char *fact, char *uplo, INTEGER *n, INTEGER *
	nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *af, INTEGER *
	ldaf, INTEGER *ipiv, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x,
	 INTEGER *ldx, DOUBLE *rcond, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, INTEGER *lwork, DOUBLE *rwork, INTEGER *info);
 
 int zsytf2_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, INTEGER *info);
 
 int zsytrf_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *lwork, 
	INTEGER *info);
 
 int zsytri_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zsytrs_(char *uplo, INTEGER *n, INTEGER *nrhs, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *ipiv, DOUBLECOMPLEX *b, 
	INTEGER *ldb, INTEGER *info);
 
 int ztbcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	INTEGER *kd, DOUBLECOMPLEX *ab, INTEGER *ldab, DOUBLE *rcond, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int ztbrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, 
	DOUBLE *ferr, DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *
	rwork, INTEGER *info);
 
 int ztbtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *kd, INTEGER *nrhs, DOUBLECOMPLEX *ab, INTEGER *ldab, 
	DOUBLECOMPLEX *b, INTEGER *ldb, INTEGER *info);
 
 int ztgevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER 
	*ldb, DOUBLECOMPLEX *vl, INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *
	ldvr, INTEGER *mm, INTEGER *m, DOUBLECOMPLEX *work, DOUBLE *rwork,
	 INTEGER *info);
 
 int ztgex2_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	INTEGER *j1, INTEGER *info);
 
 int ztgexc_(LOGICAL *wantq, LOGICAL *wantz, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *z__, INTEGER *ldz, 
	INTEGER *ifst, INTEGER *ilst, INTEGER *info);
 
 int ztgsen_(INTEGER *ijob, LOGICAL *wantq, LOGICAL *wantz, 
	LOGICAL *select, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLECOMPLEX *alpha, DOUBLECOMPLEX *
	beta, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *z__, INTEGER *
	ldz, INTEGER *m, DOUBLE *pl, DOUBLE *pr, DOUBLE *dif, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *iwork, INTEGER *liwork, 
	INTEGER *info);
 
 int ztgsja_(char *jobu, char *jobv, char *jobq, INTEGER *m, 
	INTEGER *p, INTEGER *n, INTEGER *k, INTEGER *l, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, DOUBLE *tola, 
	DOUBLE *tolb, DOUBLE *alpha, DOUBLE *beta, DOUBLECOMPLEX *
	u, INTEGER *ldu, DOUBLECOMPLEX *v, INTEGER *ldv, DOUBLECOMPLEX *q, 
	INTEGER *ldq, DOUBLECOMPLEX *work, INTEGER *ncycle, INTEGER *info);
 
 int ztgsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER 
	*ldb, DOUBLECOMPLEX *vl, INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *
	ldvr, DOUBLE *s, DOUBLE *dif, INTEGER *mm, INTEGER *m, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *iwork, INTEGER *info);
 
 int ztgsy2_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *d__, INTEGER *ldd, 
	DOUBLECOMPLEX *e, INTEGER *lde, DOUBLECOMPLEX *f, INTEGER *ldf, 
	DOUBLE *scale, DOUBLE *rdsum, DOUBLE *rdscal, INTEGER *
	info);
 
 int ztgsyl_(char *trans, INTEGER *ijob, INTEGER *m, INTEGER *
	n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *d__, INTEGER *ldd, 
	DOUBLECOMPLEX *e, INTEGER *lde, DOUBLECOMPLEX *f, INTEGER *ldf, 
	DOUBLE *scale, DOUBLE *dif, DOUBLECOMPLEX *work, INTEGER *
	lwork, INTEGER *iwork, INTEGER *info);
 
 int ztpcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	DOUBLECOMPLEX *ap, DOUBLE *rcond, DOUBLECOMPLEX *work, DOUBLE 
	*rwork, INTEGER *info);
 
 int ztprfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *b, INTEGER *ldb, 
	DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, DOUBLE *berr, 
	DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int ztptri_(char *uplo, char *diag, INTEGER *n, 
	DOUBLECOMPLEX *ap, INTEGER *info);
 
 int ztptrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *b, INTEGER *ldb, 
	INTEGER *info);
 
 int ztrcon_(char *norm, char *uplo, char *diag, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLE *rcond, DOUBLECOMPLEX *
	work, DOUBLE *rwork, INTEGER *info);
 
 int ztrevc_(char *side, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLECOMPLEX *t, INTEGER *ldt, DOUBLECOMPLEX *vl, 
	INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, INTEGER *mm, INTEGER 
	*m, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *info);
 
 int ztrexc_(char *compq, INTEGER *n, DOUBLECOMPLEX *t, 
	INTEGER *ldt, DOUBLECOMPLEX *q, INTEGER *ldq, INTEGER *ifst, INTEGER *
	ilst, INTEGER *info);
 
 int ztrrfs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *x, INTEGER *ldx, DOUBLE *ferr, 
	DOUBLE *berr, DOUBLECOMPLEX *work, DOUBLE *rwork, INTEGER *
	info);
 
 int ztrsen_(char *job, char *compq, LOGICAL *select, INTEGER 
	*n, DOUBLECOMPLEX *t, INTEGER *ldt, DOUBLECOMPLEX *q, INTEGER *ldq, 
	DOUBLECOMPLEX *w, INTEGER *m, DOUBLE *s, DOUBLE *sep, 
	DOUBLECOMPLEX *work, INTEGER *lwork, INTEGER *info);
 
 int ztrsna_(char *job, char *howmny, LOGICAL *select, 
	INTEGER *n, DOUBLECOMPLEX *t, INTEGER *ldt, DOUBLECOMPLEX *vl, 
	INTEGER *ldvl, DOUBLECOMPLEX *vr, INTEGER *ldvr, DOUBLE *s, 
	DOUBLE *sep, INTEGER *mm, INTEGER *m, DOUBLECOMPLEX *work, 
	INTEGER *ldwork, DOUBLE *rwork, INTEGER *info);
 
 int ztrsyl_(char *trana, char *tranb, INTEGER *isgn, INTEGER 
	*m, INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLE *scale, 
	INTEGER *info);
 
 int ztrti2_(char *uplo, char *diag, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *info);
 
 int ztrtri_(char *uplo, char *diag, INTEGER *n, 
	DOUBLECOMPLEX *a, INTEGER *lda, INTEGER *info);
 
 int ztrtrs_(char *uplo, char *trans, char *diag, INTEGER *n, 
	INTEGER *nrhs, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *b, 
	INTEGER *ldb, INTEGER *info);
 
 int ztzrqf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, INTEGER *info);
 
 int ztzrzf_(INTEGER *m, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zung2l_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zung2r_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zungbr_(char *vect, INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zunghr_(INTEGER *n, INTEGER *ilo, INTEGER *ihi, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zungl2_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zunglq_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zungql_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zungqr_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zungr2_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zungrq_(INTEGER *m, INTEGER *n, INTEGER *k, 
	DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zungtr_(char *uplo, INTEGER *n, DOUBLECOMPLEX *a, 
	INTEGER *lda, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zunm2l_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zunm2r_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zunmbr_(char *vect, char *side, char *trans, INTEGER *m, 
	INTEGER *n, INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX 
	*tau, DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *
	lwork, INTEGER *info);
 
 int zunmhr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *ilo, INTEGER *ihi, DOUBLECOMPLEX *a, INTEGER *lda, 
	DOUBLECOMPLEX *tau, DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *
	work, INTEGER *lwork, INTEGER *info);
 
 int zunml2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zunmlq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zunmql_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zunmqr_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zunmr2_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *info);
 
 int zunmr3_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX 
	*tau, DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *
	info);
 
 int zunmrq_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zunmrz_(char *side, char *trans, INTEGER *m, INTEGER *n, 
	INTEGER *k, INTEGER *l, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX 
	*tau, DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *
	lwork, INTEGER *info);
 
 int zunmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, DOUBLECOMPLEX *a, INTEGER *lda, DOUBLECOMPLEX *tau, 
	DOUBLECOMPLEX *c__, INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *lwork,
	 INTEGER *info);
 
 int zupgtr_(char *uplo, INTEGER *n, DOUBLECOMPLEX *ap, 
	DOUBLECOMPLEX *tau, DOUBLECOMPLEX *q, INTEGER *ldq, DOUBLECOMPLEX *
	work, INTEGER *info);
 
 int zupmtr_(char *side, char *uplo, char *trans, INTEGER *m, 
	INTEGER *n, DOUBLECOMPLEX *ap, DOUBLECOMPLEX *tau, DOUBLECOMPLEX *c__,
	 INTEGER *ldc, DOUBLECOMPLEX *work, INTEGER *info);

#ifdef __cplusplus
}
#endif

#endif /* _CLAPACK_H_ */
