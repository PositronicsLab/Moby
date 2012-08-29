/* sasum.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <cmath>
#include "f2c-ap.h"

extern "C" {

// necessary for f2c
void MAIN__() {}

doublereal aasum_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static integer i__, m, mp1, nincx;
    static real stemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     takes the sum of the absolute values. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	stemp += (r__1 = sx[i__], dabs(r__1));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	stemp += (r__1 = sx[i__], dabs(r__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
	stemp = stemp + (r__1 = sx[i__], dabs(r__1)) + (r__2 = sx[i__ + 1], 
		dabs(r__2)) + (r__3 = sx[i__ + 2], dabs(r__3)) + (r__4 = sx[
		i__ + 3], dabs(r__4)) + (r__5 = sx[i__ + 4], dabs(r__5)) + (
		r__6 = sx[i__ + 5], dabs(r__6));
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* aasum_ */

/* saxpy.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int aaxpy_(integer *n, real *sa, real *sx, integer *incx, 
	real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     SAXPY constant times a vector plus a vector. */
/*     uses unrolled loop for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*sa == (real) 0.f) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] += *sa * sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] += *sa * sx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	sy[i__] += *sa * sx[i__];
	sy[i__ + 1] += *sa * sx[i__ + 1];
	sy[i__ + 2] += *sa * sx[i__ + 2];
	sy[i__ + 3] += *sa * sx[i__ + 3];
/* L50: */
    }
    return 0;
} /* aaxpy_ */

/* scabs1.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal acabs1_(complex *z__)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *);

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SCABS1 computes absolute value of a complex number */

/*     .. Intrinsic Functions .. */
/*     .. */
    ret_val = (r__1 = z__->r, dabs(r__1)) + (r__2 = r_imag(z__), dabs(r__2));
    return ret_val;
} /* scabs1_ */

/* scasum.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal acasum_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    static integer i__, nincx;
    static real stemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     takes the sum of the absolute values of a complex vector and */
/*     returns a single precision result. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__;
	stemp = stemp + (r__1 = cx[i__3].r, dabs(r__1)) + (r__2 = r_imag(&cx[
		i__]), dabs(r__2));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */

L20:
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	stemp = stemp + (r__1 = cx[i__1].r, dabs(r__1)) + (r__2 = r_imag(&cx[
		i__]), dabs(r__2));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
} /* acasum_ */

/* scnrm2.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal acnrm2_(integer *n, complex *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1;

    /* Builtin functions */
    double r_imag(complex *), sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static real ssq, temp, norm, scale;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SCNRM2 returns the euclidean norm of a vector via the function */
/*  name, so that */

/*     SCNRM2 := sqrt( conjg( x' )*x ) */



/*  -- This version written on 25-October-1982. */
/*     Modified on 14-October-1993 to inline the call to CLASSQ. */
/*     Sven Hammarling, Nag Ltd. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = (real) 0.f;
    } else {
	scale = (real) 0.f;
	ssq = (real) 1.f;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL CLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    i__3 = ix;
	    if (x[i__3].r != (real) 0.f) {
		i__3 = ix;
		temp = (r__1 = x[i__3].r, dabs(r__1));
		if (scale < temp) {
/* Computing 2nd power */
		    r__1 = scale / temp;
		    ssq = ssq * (r__1 * r__1) + (real) 1.f;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    r__1 = temp / scale;
		    ssq += r__1 * r__1;
		}
	    }
	    if (r_imag(&x[ix]) != (real) 0.f) {
		temp = (r__1 = r_imag(&x[ix]), dabs(r__1));
		if (scale < temp) {
/* Computing 2nd power */
		    r__1 = scale / temp;
		    ssq = ssq * (r__1 * r__1) + (real) 1.f;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    r__1 = temp / scale;
		    ssq += r__1 * r__1;
		}
	    }
/* L10: */
	}
	norm = scale * std::sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of SCNRM2. */

} /* acnrm2_ */

/* scopy.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int acopy_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] = sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] = sx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	sy[i__] = sx[i__];
	sy[i__ + 1] = sx[i__ + 1];
	sy[i__ + 2] = sx[i__ + 2];
	sy[i__ + 3] = sx[i__ + 3];
	sy[i__ + 4] = sx[i__ + 4];
	sy[i__ + 5] = sx[i__ + 5];
	sy[i__ + 6] = sx[i__ + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */

/* sdot.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal adot_(integer *n, real *sx, integer *incx, real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static real stemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[ix] * sy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[i__] * sy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] + sx[
		i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] + sx[i__ + 
		4] * sy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* adot_ */

/* sdsdot.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal adsdot_(integer *n, real *sb, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static integer i__, ns, kx, ky;
    static doublereal dsdot;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  PURPOSE */
/*  ======= */

/*  Compute the inner product of two vectors with extended */
/*  precision accumulation. */

/*  Returns S.P. result with dot product accumulated in D.P. */
/*  SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY), */
/*  where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*  defined in a similar way using INCY. */

/*  AUTHOR */
/*  ====== */
/*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA), */
/*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL) */

/*  ARGUMENTS */
/*  ========= */

/*  N      (input) INTEGER */
/*         number of elements in input vector(s) */

/*  SB     (input) REAL */
/*         single precision scalar to be added to inner product */

/*  SX     (input) REAL array, dimension (N) */
/*         single precision vector with N elements */

/*  INCX   (input) INTEGER */
/*         storage spacing between elements of SX */

/*  SY     (input) REAL array, dimension (N) */
/*         single precision vector with N elements */

/*  INCY   (input) INTEGER */
/*         storage spacing between elements of SY */

/*  SDSDOT (output) REAL */
/*         single precision dot product (SB if N .LE. 0) */

/*  REFERENCES */
/*  ========== */

/*  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*  Krogh, Basic linear algebra subprograms for Fortran */
/*  usage, Algorithm No. 539, Transactions on Mathematical */
/*  Software 5, 3 (September 1979), pp. 308-323. */

/*  REVISION HISTORY  (YYMMDD) */
/*  ========================== */

/*  791001  DATE WRITTEN */
/*  890531  Changed all specific intrinsics to generic.  (WRB) */
/*  890831  Modified array declarations.  (WRB) */
/*  890831  REVISION DATE from Version 3.2 */
/*  891214  Prologue converted to Version 4.0 format.  (BAB) */
/*  920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*  920501  Reformatted the REFERENCES section.  (WRB) */
/*  070118  Reformat to LAPACK coding style */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    dsdot = *sb;
    if (*n <= 0) {
	goto L30;
    }
    if (*incx == *incy && *incx > 0) {
	goto L40;
    }

/*     Code for unequal or nonpositive increments. */

    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dsdot += (doublereal) sx[kx] * (doublereal) sy[ky];
	kx += *incx;
	ky += *incy;
/* L10: */
    }
L30:
    ret_val = dsdot;
    return ret_val;

/*     Code for equal and positive increments. */

L40:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dsdot += (doublereal) sx[i__] * (doublereal) sy[i__];
/* L50: */
    }
    ret_val = dsdot;
    return ret_val;
} /* adsdot_ */

/* sgbmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int agbmv_(char *trans, integer *m, integer *n, integer *kl, 
	integer *ku, real *alpha, real *a, integer *lda, real *x, integer *
	incx, real *beta, real *y, integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    static real temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGBMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n band matrix, with kl sub-diagonals and ku super-diagonals. */

/*  Arguments */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  KL     - INTEGER. */
/*           On entry, KL specifies the number of sub-diagonals of the */
/*           matrix A. KL must satisfy  0 .le. KL. */
/*           Unchanged on exit. */

/*  KU     - INTEGER. */
/*           On entry, KU specifies the number of super-diagonals of the */
/*           matrix A. KU must satisfy  0 .le. KU. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading ( kl + ku + 1 ) by n part of the */
/*           array A must contain the matrix of coefficients, supplied */
/*           column by column, with the leading diagonal of the matrix in */
/*           row ( ku + 1 ) of the array, the first super-diagonal */
/*           starting at position 2 in row ku, the first sub-diagonal */
/*           starting at position 1 in row ( ku + 2 ), and so on. */
/*           Elements in the array A that do not correspond to elements */
/*           in the band matrix (such as the top left ku by ku triangle) */
/*           are not referenced. */
/*           The following program segment will transfer a band matrix */
/*           from conventional full matrix storage to band storage: */

/*                 DO 20, J = 1, N */
/*                    K = KU + 1 - J */
/*                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) */
/*                       A( K + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           ( kl + ku + 1 ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry, the incremented array Y must contain the */
/*           vector y. On exit, Y is overwritten by the updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*kl < 0) {
	info = 4;
    } else if (*ku < 0) {
	info = 5;
    } else if (*lda < *kl + *ku + 1) {
	info = 8;
    } else if (*incx == 0) {
	info = 10;
    } else if (*incy == 0) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("SGBMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == (real) 0.f && *beta == (real) 1.f) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the band part of A. */

/*     First form  y := beta*y. */

    if (*beta != (real) 1.f) {
	if (*incy == 1) {
	    if (*beta == (real) 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = (real) 0.f;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == (real) 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == (real) 0.f) {
	return 0;
    }
    kup1 = *ku + 1;
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.f) {
		    temp = *alpha * x[jx];
		    k = kup1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__4 = min(i__5,i__6);
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			y[i__] += temp * a[k + i__ + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.f) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    k = kup1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__3 = min(i__5,i__6);
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			y[iy] += temp * a[k + i__ + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		if (j > *ku) {
		    ky += *incy;
		}
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		k = kup1 - j;
/* Computing MAX */
		i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
		    temp += a[k + i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		ix = kx;
		k = kup1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__4 = min(i__5,i__6);
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    temp += a[k + i__ + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
		if (j > *ku) {
		    kx += *incx;
		}
/* L120: */
	    }
	}
    }

    return 0;

/*     End of SGBMV . */

} /* sgbmv_ */

/* sgemm.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int agemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *
	ldb, real *beta, real *c__, integer *ldc, ftnlen transa_len, ftnlen 
	transb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static logical nota, notb;
    static real temp;
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*  where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X', */

/*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

/*  Arguments */
/*  ========== */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n',  op( A ) = A. */

/*              TRANSA = 'T' or 't',  op( A ) = A'. */

/*              TRANSA = 'C' or 'c',  op( A ) = A'. */

/*           Unchanged on exit. */

/*  TRANSB - CHARACTER*1. */
/*           On entry, TRANSB specifies the form of op( B ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSB = 'N' or 'n',  op( B ) = B. */

/*              TRANSB = 'T' or 't',  op( B ) = B'. */

/*              TRANSB = 'C' or 'c',  op( B ) = B'. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies  the number  of rows  of the  matrix */
/*           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N  specifies the number  of columns of the matrix */
/*           op( B ) and the number of columns of the matrix C. N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry,  K  specifies  the number of columns of the matrix */
/*           op( A ) and the number of rows of the matrix op( B ). K must */
/*           be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is */
/*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/*           part of the array  A  must contain the matrix  A,  otherwise */
/*           the leading  k by m  part of the array  A  must contain  the */
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/*           least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is */
/*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/*           part of the array  B  must contain the matrix  B,  otherwise */
/*           the leading  n by k  part of the array  B  must contain  the */
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/*           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - REAL             array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must */
/*           contain the matrix  C,  except when  beta  is zero, in which */
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n  matrix */
/*           ( alpha*op( A )*op( B ) + beta*C ). */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
/*     and  columns of  A  and the  number of  rows  of  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! lsame_(transa, "C", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    transa, "T", (ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! notb && ! lsame_(transb, "C", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(transb, "T", (ftnlen)1, (ftnlen)1)) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("SGEMM ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == (real) 0.f || *k == 0) && *beta == (real) 1.f) {
	return 0;
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == (real) 0.f) {
	if (*beta == (real) 0.f) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = (real) 0.f;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.f;
/* L50: */
		    }
		} else if (*beta != (real) 1.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[l + j * b_dim1] != (real) 0.f) {
			temp = *alpha * b[l + j * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = (real) 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
/* L100: */
		    }
		    if (*beta == (real) 0.f) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := alpha*A*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.f;
/* L130: */
		    }
		} else if (*beta != (real) 1.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[j + l * b_dim1] != (real) 0.f) {
			temp = *alpha * b[j + l * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = (real) 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
/* L180: */
		    }
		    if (*beta == (real) 0.f) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }

    return 0;

/*     End of SGEMM . */

} /* agemm_ */

/* sgemv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int agemv_(char *trans, integer *m, integer *n, real *alpha, 
	real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static real temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Arguments */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the */
/*           updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("SGEMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == (real) 0.0 && *beta == (real) 1.0) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != (real) 1.0) {
	if (*incy == 1) {
	    if (*beta == 0.0) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = (real) 0.0;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == (real) 0.0) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = (real) 0.0;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == (real) 0.0) {
	return 0;
    }
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a[i__ + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = (real) 0.0;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = (real) 0.0;
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of SGEMV . */

} /* agemv_ */

/* sger.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int ager_(integer *m, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jy, kx, info;
    static real temp;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGER   performs the rank 1 operation */

/*     A := alpha*x*y' + A, */

/*  where alpha is a scalar, x is an m element vector, y is an n element */
/*  vector and A is an m by n matrix. */

/*  Arguments */
/*  ========== */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the m */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. On exit, A is */
/*           overwritten by the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SGER  ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == (real) 0.0) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != (real) 0.0) {
		temp = *alpha * y[jy];
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[i__] * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != (real) 0.0) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[ix] * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

/*     End of SGER  . */

} /* ager_ */

/* snrm2.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



doublereal anrm2_(integer *n, real *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static real ssq, norm, scale, absxi;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SNRM2 returns the euclidean norm of a vector via the function */
/*  name, so that */

/*     SNRM2 := sqrt( x'*x ). */

/*  Further Details */
/*  =============== */

/*  -- This version written on 25-October-1982. */
/*     Modified on 14-October-1993 to inline the call to SLASSQ. */
/*     Sven Hammarling, Nag Ltd. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = (real) 0.0;
    } else if (*n == 1) {
	norm = dabs(x[1]);
    } else {
	scale = (real) 0.0;
	ssq = (real) 1.0;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != (real) 0.0) {
		absxi = (r__1 = x[ix], dabs(r__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    r__1 = scale / absxi;
		    ssq = ssq * (r__1 * r__1) + (real) 1.0;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    r__1 = absxi / scale;
		    ssq += r__1 * r__1;
		}
	    }
/* L10: */
	}
	norm = scale * std::sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of SNRM2. */

} /* anrm2_ */

/* srot.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int arot_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy, real *c__, real *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
    static real stemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     applies a plane rotation. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = *c__ * sx[ix] + *s * sy[iy];
	sy[iy] = *c__ * sy[iy] - *s * sx[ix];
	sx[ix] = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = *c__ * sx[i__] + *s * sy[i__];
	sy[i__] = *c__ * sy[i__] - *s * sx[i__];
	sx[i__] = stemp;
/* L30: */
    }
    return 0;
} /* arot_ */

/* srotg.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Table of constant values */

static real c_b4 = (real) 1.0;

/* Subroutine */ int arotg_(real *sa, real *sb, real *c__, real *s)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static real r__, z__, roe, scale;

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     construct givens plane rotation. */
/*     jack dongarra, linpack, 3/11/78. */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    roe = *sb;
    if (dabs(*sa) > dabs(*sb)) {
	roe = *sa;
    }
    scale = dabs(*sa) + dabs(*sb);
    if (scale != (real) 0.0) {
	goto L10;
    }
    *c__ = (real) 1.0;
    *s = (real) 0.0;
    r__ = (real) 0.0;
    z__ = (real) 0.0;
    goto L20;
L10:
/* Computing 2nd power */
    r__1 = *sa / scale;
/* Computing 2nd power */
    r__2 = *sb / scale;
    r__ = scale * std::sqrt(r__1 * r__1 + r__2 * r__2);
    r__ = r_sign(&c_b4, &roe) * r__;
    *c__ = *sa / r__;
    *s = *sb / r__;
    z__ = (real) 1.0;
    if (dabs(*sa) > dabs(*sb)) {
	z__ = *s;
    }
    if (dabs(*sb) >= dabs(*sa) && *c__ != (real) 0.0) {
	z__ = (real) 1.0 / *c__;
    }
L20:
    *sa = r__;
    *sb = z__;
    return 0;
} /* arotg_ */

/* srotm.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int arotm_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy, real *sparam)
{
    /* Initialized data */

    static real zero = (real) 0.0;
    static real two = (real) 2.0;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static real w, z__;
    static integer kx, ky;
    static real sh11, sh12, sh21, sh22, sflag;
    static integer nsteps;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */

/*     (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN */
/*     (DX**T) */

/*     SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/*     LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY. */
/*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */

/*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/*     SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM. */


/*  Arguments */
/*  ========= */

/*  N      (input) INTEGER */
/*         number of elements in input vector(s) */

/*  SX     (input/output) REAL array, dimension N */
/*         double precision vector with 5 elements */

/*  INCX   (input) INTEGER */
/*         storage spacing between elements of SX */

/*  SY     (input/output) REAL array, dimension N */
/*         double precision vector with N elements */

/*  INCY   (input) INTEGER */
/*         storage spacing between elements of SY */

/*  SPARAM (input/output)  REAL array, dimension 5 */
/*     SPARAM(1)=SFLAG */
/*     SPARAM(2)=SH11 */
/*     SPARAM(3)=SH21 */
/*     SPARAM(4)=SH12 */
/*     SPARAM(5)=SH22 */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Data statements .. */
    /* Parameter adjustments */
    --sparam;
    --sy;
    --sx;

    /* Function Body */
/*     .. */

    sflag = sparam[1];
    if (*n <= 0 || sflag + two == zero) {
	goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
	goto L70;
    }

    nsteps = *n * *incx;
    if (sflag < (real) 0.0) {
	goto L50;
    } else if (sflag == 0) {
	goto L10;
    } else {
	goto L30;
    }
L10:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w + z__ * sh12;
	sy[i__] = w * sh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w * sh11 + z__;
	sy[i__] = -w + sh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w * sh11 + z__ * sh12;
	sy[i__] = w * sh21 + z__ * sh22;
/* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }

    if (sflag < (real) 0.0) {
	goto L120;
    } else if (sflag == 0) {
	goto L80;
    } else {
	goto L100;
    }
L80:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w + z__ * sh12;
	sy[ky] = w * sh21 + z__;
	kx += *incx;
	ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w * sh11 + z__;
	sy[ky] = -w + sh22 * z__;
	kx += *incx;
	ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w * sh11 + z__ * sh12;
	sy[ky] = w * sh21 + z__ * sh22;
	kx += *incx;
	ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* srotm_ */

/* srotmg.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int arotmg_(real *sd1, real *sd2, real *sx1, real *sy1, real 
	*sparam)
{
    /* Initialized data */

    static real zero = (real) 0.0;
    static real one = (real) 1.0;
    static real two = (real) 2.0;
    static real gam = (real) 4096.0;
    static real gamsq = (real) 16777200.0;
    static real rgamsq = (real) 5.96046e-8;

    /* Format strings */
    static char fmt_120[] = "";
    static char fmt_150[] = "";
    static char fmt_180[] = "";
    static char fmt_210[] = "";

    /* System generated locals */
    real r__1;

    /* Local variables */
    static real su, sp1, sp2, sq1, sq2, sh11, sh12, sh21, sh22;
    static integer igo;
    static real sflag, stemp;

    /* Assigned format variables */
    static char *igo_fmt;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/*     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)* */
/*     SY2)**T. */
/*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */

/*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/*     LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22 */
/*     RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE */
/*     VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.) */

/*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE */
/*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE */
/*     OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM. */


/*  Arguments */
/*  ========= */


/*  SD1    (input/output) REAL */

/*  SD2    (input/output) REAL */

/*  SX1    (input/output) REAL */

/*  SY1    (input) REAL */


/*  SPARAM (input/output)  REAL array, dimension 5 */
/*     SPARAM(1)=SFLAG */
/*     SPARAM(2)=SH11 */
/*     SPARAM(3)=SH21 */
/*     SPARAM(4)=SH12 */
/*     SPARAM(5)=SH22 */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */

    /* Parameter adjustments */
    --sparam;

    /* Function Body */
/*     .. */
    if (! (*sd1 < zero)) {
	goto L10;
    }
/*       GO ZERO-H-D-AND-SX1.. */
    goto L60;
L10:
/*     CASE-SD1-NONNEGATIVE */
    sp2 = *sd2 * *sy1;
    if (! (sp2 == zero)) {
	goto L20;
    }
    sflag = -two;
    goto L260;
/*     REGULAR-CASE.. */
L20:
    sp1 = *sd1 * *sx1;
    sq2 = sp2 * *sy1;
    sq1 = sp1 * *sx1;

    if (! (dabs(sq1) > dabs(sq2))) {
	goto L40;
    }
    sh21 = -(*sy1) / *sx1;
    sh12 = sp2 / sp1;

    su = one - sh12 * sh21;

    if (! (su <= zero)) {
	goto L30;
    }
/*         GO ZERO-H-D-AND-SX1.. */
    goto L60;
L30:
    sflag = zero;
    *sd1 /= su;
    *sd2 /= su;
    *sx1 *= su;
/*         GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (sq2 < zero)) {
	goto L50;
    }
/*         GO ZERO-H-D-AND-SX1.. */
    goto L60;
L50:
    sflag = one;
    sh11 = sp1 / sp2;
    sh22 = *sx1 / *sy1;
    su = one + sh11 * sh22;
    stemp = *sd2 / su;
    *sd2 = *sd1 / su;
    *sd1 = stemp;
    *sx1 = *sy1 * su;
/*         GO SCALE-CHECK */
    goto L100;
/*     PROCEDURE..ZERO-H-D-AND-SX1.. */
L60:
    sflag = -one;
    sh11 = zero;
    sh12 = zero;
    sh21 = zero;
    sh22 = zero;

    *sd1 = zero;
    *sd2 = zero;
    *sx1 = zero;
/*         RETURN.. */
    goto L220;
/*     PROCEDURE..FIX-H.. */
L70:
    if (! (sflag >= zero)) {
	goto L90;
    }

    if (! (sflag == zero)) {
	goto L80;
    }
    sh11 = one;
    sh22 = one;
    sflag = -one;
    goto L90;
L80:
    sh21 = -one;
    sh12 = one;
    sflag = -one;
L90:
    switch (igo) {
	case 0: goto L120;
	case 1: goto L150;
	case 2: goto L180;
	case 3: goto L210;
    }
/*     PROCEDURE..SCALE-CHECK */
L100:
L110:
    if (! (*sd1 <= rgamsq)) {
	goto L130;
    }
    if (*sd1 == zero) {
	goto L160;
    }
    igo = 0;
    igo_fmt = fmt_120;
/*              FIX-H.. */
    goto L70;
L120:
/* Computing 2nd power */
    r__1 = gam;
    *sd1 *= r__1 * r__1;
    *sx1 /= gam;
    sh11 /= gam;
    sh12 /= gam;
    goto L110;
L130:
L140:
    if (! (*sd1 >= gamsq)) {
	goto L160;
    }
    igo = 1;
    igo_fmt = fmt_150;
/*              FIX-H.. */
    goto L70;
L150:
/* Computing 2nd power */
    r__1 = gam;
    *sd1 /= r__1 * r__1;
    *sx1 *= gam;
    sh11 *= gam;
    sh12 *= gam;
    goto L140;
L160:
L170:
    if (! (dabs(*sd2) <= rgamsq)) {
	goto L190;
    }
    if (*sd2 == zero) {
	goto L220;
    }
    igo = 2;
    igo_fmt = fmt_180;
/*              FIX-H.. */
    goto L70;
L180:
/* Computing 2nd power */
    r__1 = gam;
    *sd2 *= r__1 * r__1;
    sh21 /= gam;
    sh22 /= gam;
    goto L170;
L190:
L200:
    if (! (dabs(*sd2) >= gamsq)) {
	goto L220;
    }
    igo = 3;
    igo_fmt = fmt_210;
/*              FIX-H.. */
    goto L70;
L210:
/* Computing 2nd power */
    r__1 = gam;
    *sd2 /= r__1 * r__1;
    sh21 *= gam;
    sh22 *= gam;
    goto L200;
L220:
    if (sflag < (real) 0.0) {
	goto L250;
    } else if (sflag == 0) {
	goto L230;
    } else {
	goto L240;
    }
L230:
    sparam[3] = sh21;
    sparam[4] = sh12;
    goto L260;
L240:
    sparam[2] = sh11;
    sparam[5] = sh22;
    goto L260;
L250:
    sparam[2] = sh11;
    sparam[3] = sh21;
    sparam[4] = sh12;
    sparam[5] = sh22;
L260:
    sparam[1] = sflag;
    return 0;
} /* srotmg_ */

/* ssbmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asbmv_(char *uplo, integer *n, integer *k, real *alpha, 
	real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSBMV  performs the matrix-vector  operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric band matrix, with k super-diagonals. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the band matrix A is being supplied as */
/*           follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  being supplied. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  being supplied. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry, K specifies the number of super-diagonals of the */
/*           matrix A. K must satisfy  0 .le. K. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/*           by n part of the array A must contain the upper triangular */
/*           band part of the symmetric matrix, supplied column by */
/*           column, with the leading diagonal of the matrix in row */
/*           ( k + 1 ) of the array, the first super-diagonal starting at */
/*           position 2 in row k, and so on. The top left k by k triangle */
/*           of the array A is not referenced. */
/*           The following program segment will transfer the upper */
/*           triangular part of a symmetric band matrix from conventional */
/*           full matrix storage to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = K + 1 - J */
/*                    DO 10, I = MAX( 1, J - K ), J */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/*           by n part of the array A must contain the lower triangular */
/*           band part of the symmetric matrix, supplied column by */
/*           column, with the leading diagonal of the matrix in row 1 of */
/*           the array, the first sub-diagonal starting at position 1 in */
/*           row 2, and so on. The bottom right k by k triangle of the */
/*           array A is not referenced. */
/*           The following program segment will transfer the lower */
/*           triangular part of a symmetric band matrix from conventional */
/*           full matrix storage to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = 1 - J */
/*                    DO 10, I = J, MIN( N, J + K ) */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           ( k + 1 ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the */
/*           vector y. On exit, Y is overwritten by the updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k < 0) {
	info = 3;
    } else if (*lda < *k + 1) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("SSBMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0 && *beta == (real) 1.0) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of the array A */
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != (real) 1.0) {
	if (*incy == 1) {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = (real) 0.0;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = (real) 0.0;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == (real) 0.0) {
	return 0;
    }
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when upper triangle of A is stored. */

	kplus1 = *k + 1;
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		l = kplus1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *k;
		i__4 = j - 1;
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		ix = kx;
		iy = ky;
		l = kplus1 - j;
/* Computing MAX */
		i__4 = 1, i__2 = j - *k;
		i__3 = j - 1;
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a[kplus1 + j * a_dim1] + *alpha * 
			temp2;
		jx += *incx;
		jy += *incy;
		if (j > *k) {
		    kx += *incx;
		    ky += *incy;
		}
/* L80: */
	    }
	}
    } else {

/*        Form  y  when lower triangle of A is stored. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		y[j] += temp1 * a[j * a_dim1 + 1];
		l = 1 - j;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		y[jy] += temp1 * a[j * a_dim1 + 1];
		l = 1 - j;
		ix = jx;
		iy = jy;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of SSBMV . */

} /* asbmv_ */

/* sscal.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int ascal_(integer *n, real *sa, real *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, mp1, nincx;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	sx[i__] = *sa * sx[i__];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	sx[i__] = *sa * sx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	sx[i__] = *sa * sx[i__];
	sx[i__ + 1] = *sa * sx[i__ + 1];
	sx[i__ + 2] = *sa * sx[i__ + 2];
	sx[i__ + 3] = *sa * sx[i__ + 3];
	sx[i__ + 4] = *sa * sx[i__ + 4];
/* L50: */
    }
    return 0;
} /* ascal_ */

/* sspmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int aspmv_(char *uplo, integer *n, real *alpha, real *ap, 
	real *x, integer *incx, real *beta, real *y, integer *incy, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSPMV  performs the matrix-vector operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric matrix, supplied in packed form. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  AP     - REAL             array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. On exit, Y is overwritten by the updated */
/*           vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --y;
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 6;
    } else if (*incy == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SSPMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0 && *beta == (real) 1.0) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

    if (*beta != (real) 1.0) {
	if (*incy == 1) {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = (real) 0.0;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = (real) 0.0;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == (real) 0.0) {
	return 0;
    }
    kk = 1;
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		k = kk;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * ap[k];
		    temp2 += ap[k] * x[i__];
		    ++k;
/* L50: */
		}
		y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
		kk += j;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		ix = kx;
		iy = ky;
		i__2 = kk + j - 2;
		for (k = kk; k <= i__2; ++k) {
		    y[iy] += temp1 * ap[k];
		    temp2 += ap[k] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += j;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when AP contains the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		y[j] += temp1 * ap[kk];
		k = kk + 1;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * ap[k];
		    temp2 += ap[k] * x[i__];
		    ++k;
/* L90: */
		}
		y[j] += *alpha * temp2;
		kk += *n - j + 1;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		y[jy] += temp1 * ap[kk];
		ix = jx;
		iy = jy;
		i__2 = kk + *n - j;
		for (k = kk + 1; k <= i__2; ++k) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * ap[k];
		    temp2 += ap[k] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += *n - j + 1;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of SSPMV . */

} /* aspmv_ */

/* sspr2.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int aspr2_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSPR2  performs the symmetric rank 2 operation */

/*     A := alpha*x*y' + alpha*y*x' + A, */

/*  where alpha is a scalar, x and y are n element vectors and A is an */
/*  n by n symmetric matrix, supplied in packed form. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  AP     - REAL             array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the upper triangular part of the */
/*           updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the lower triangular part of the */
/*           updated matrix. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;
    --y;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("SSPR2 ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0 || y[j] != (real) 0.0) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    k = kk;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0 || y[jy] != (real) 0.0) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0 || y[j] != (real) 0.0) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    k = kk;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0 || y[jy] != (real) 0.0) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of SSPR2 . */

} /* aspr2_ */

/* sspr.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int aspr_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSPR    performs the symmetric rank 1 operation */

/*     A := alpha*x*x' + A, */

/*  where alpha is a real scalar, x is an n element vector and A is an */
/*  n by n symmetric matrix, supplied in packed form. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  AP     - REAL             array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the upper triangular part of the */
/*           updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the lower triangular part of the */
/*           updated matrix. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    }
    if (info != 0) {
	xerbla_("SSPR  ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0) {
	return 0;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    ix = kx;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    ix = jx;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of SSPR  . */

} /* sspr_ */

/* sswap.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int aswap_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static real stemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = sx[ix];
	sx[ix] = sy[iy];
	sy[iy] = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = sx[i__];
	sx[i__] = sy[i__];
	sy[i__] = stemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	stemp = sx[i__];
	sx[i__] = sy[i__];
	sy[i__] = stemp;
	stemp = sx[i__ + 1];
	sx[i__ + 1] = sy[i__ + 1];
	sy[i__ + 1] = stemp;
	stemp = sx[i__ + 2];
	sx[i__ + 2] = sy[i__ + 2];
	sy[i__ + 2] = stemp;
/* L50: */
    }
    return 0;
} /* sswap_ */

/* ssymm.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asymm_(char *side, char *uplo, integer *m, integer *n, 
	real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta,
	 real *c__, integer *ldc, ftnlen side_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYMM  performs one of the matrix-matrix operations */

/*     C := alpha*A*B + beta*C, */

/*  or */

/*     C := alpha*B*A + beta*C, */

/*  where alpha and beta are scalars,  A is a symmetric matrix and  B and */
/*  C are  m by n matrices. */

/*  Arguments */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry,  SIDE  specifies whether  the  symmetric matrix  A */
/*           appears on the  left or right  in the  operation as follows: */

/*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C, */

/*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C, */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/*           triangular  part  of  the  symmetric  matrix   A  is  to  be */
/*           referenced as follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of the */
/*                                  symmetric matrix is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of the */
/*                                  symmetric matrix is to be referenced. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies the number of rows of the matrix  C. */
/*           M  must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix C. */
/*           N  must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is */
/*           m  when  SIDE = 'L' or 'l'  and is  n otherwise. */
/*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of */
/*           the array  A  must contain the  symmetric matrix,  such that */
/*           when  UPLO = 'U' or 'u', the leading m by m upper triangular */
/*           part of the array  A  must contain the upper triangular part */
/*           of the  symmetric matrix and the  strictly  lower triangular */
/*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/*           the leading  m by m  lower triangular part  of the  array  A */
/*           must  contain  the  lower triangular part  of the  symmetric */
/*           matrix and the  strictly upper triangular part of  A  is not */
/*           referenced. */
/*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of */
/*           the array  A  must contain the  symmetric matrix,  such that */
/*           when  UPLO = 'U' or 'u', the leading n by n upper triangular */
/*           part of the array  A  must contain the upper triangular part */
/*           of the  symmetric matrix and the  strictly  lower triangular */
/*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/*           the leading  n by n  lower triangular part  of the  array  A */
/*           must  contain  the  lower triangular part  of the  symmetric */
/*           matrix and the  strictly upper triangular part of  A  is not */
/*           referenced. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - REAL             array of DIMENSION ( LDB, n ). */
/*           Before entry, the leading  m by n part of the array  B  must */
/*           contain the matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - REAL             array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must */
/*           contain the matrix  C,  except when  beta  is zero, in which */
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n updated */
/*           matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Set NROWA as the number of rows of A. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

    info = 0;
    if (! lsame_(side, "L", (ftnlen)1, (ftnlen)1) && ! lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,*m)) {
	info = 9;
    } else if (*ldc < max(1,*m)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("SSYMM ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == (real) 0.0 && *beta == (real) 1.0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == (real) 0.0) {
	if (*beta == (real) 0.0) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = (real) 0.0;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = *alpha * b[i__ + j * b_dim1];
		    temp2 = (real) 0.0;
		    i__3 = i__ - 1;
		    for (k = 1; k <= i__3; ++k) {
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
/* L50: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = *m; i__ >= 1; --i__) {
		    temp1 = *alpha * b[i__ + j * b_dim1];
		    temp2 = (real) 0.0;
		    i__2 = *m;
		    for (k = i__ + 1; k <= i__2; ++k) {
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
/* L80: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
    } else {

/*        Form  C := alpha*B*A + beta*C. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    temp1 = *alpha * a[j + j * a_dim1];
	    if (*beta == (real) 0.0) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
/* L110: */
		}
	    } else {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] + 
			    temp1 * b[i__ + j * b_dim1];
/* L120: */
		}
	    }
	    i__2 = j - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a[k + j * a_dim1];
		} else {
		    temp1 = *alpha * a[j + k * a_dim1];
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
/* L130: */
		}
/* L140: */
	    }
	    i__2 = *n;
	    for (k = j + 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a[j + k * a_dim1];
		} else {
		    temp1 = *alpha * a[k + j * a_dim1];
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
    }

    return 0;

/*     End of SSYMM . */

} /* ssymm_ */

/* ssymv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asymv_(char *uplo, integer *n, real *alpha, real *a, 
	integer *lda, real *x, integer *incx, real *beta, real *y, integer *
	incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYMV  performs the matrix-vector  operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric matrix. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. On exit, Y is overwritten by the updated */
/*           vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("SSYMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0 && *beta == (real) 1.0) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

    if (*beta != (real) 1.0) {
	if (*incy == 1) {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = (real) 0.0;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = (real) 0.0;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == (real) 0.0) {
	return 0;
    }
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = (real) 0.0;
		y[j] += temp1 * a[j + j * a_dim1];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = (real) 0.0;
		y[jy] += temp1 * a[j + j * a_dim1];
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of SSYMV . */

} /* ssymv_ */

/* ssyr2.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asyr2_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYR2  performs the symmetric rank 2 operation */

/*     A := alpha*x*y' + alpha*y*x' + A, */

/*  where alpha is a scalar, x and y are n element vectors and A is an n */
/*  by n symmetric matrix. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. On exit, the */
/*           upper triangular part of the array A is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. On exit, the */
/*           lower triangular part of the array A is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SSYR2 ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0 || y[j] != (real) 0.0) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0 || y[jy] != (real) 0.0) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0 || y[j] != (real) 0.0) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0 || y[jy] != (real) 0.0) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of SSYR2 . */

} /* ssyr2_ */

/* ssyr2k.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta,
	 real *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYR2K  performs one of the symmetric rank 2k operations */

/*     C := alpha*A*B' + alpha*B*A' + beta*C, */

/*  or */

/*     C := alpha*A'*B + alpha*B'*A + beta*C, */

/*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/*  and  A and B  are  n by k  matrices  in the  first  case  and  k by n */
/*  matrices in the second case. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/*           triangular  part  of the  array  C  is to be  referenced  as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry,  TRANS  specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' + */
/*                                        beta*C. */

/*              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A + */
/*                                        beta*C. */

/*              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A + */
/*                                        beta*C. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N specifies the order of the matrix C.  N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/*           of  columns  of the  matrices  A and B,  and on  entry  with */
/*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/*           of rows of the matrices  A and B.  K must be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is */
/*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/*           part of the array  A  must contain the matrix  A,  otherwise */
/*           the leading  k by n  part of the array  A  must contain  the */
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/*           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/*           be at least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is */
/*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/*           part of the array  B  must contain the matrix  B,  otherwise */
/*           the leading  k by n  part of the array  B  must contain  the */
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/*           then  LDB must be at least  max( 1, n ), otherwise  LDB must */
/*           be at least  max( 1, k ). */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. */
/*           Unchanged on exit. */

/*  C      - REAL             array of DIMENSION ( LDC, n ). */
/*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/*           upper triangular part of the array C must contain the upper */
/*           triangular part  of the  symmetric matrix  and the strictly */
/*           lower triangular part of C is not referenced.  On exit, the */
/*           upper triangular part of the array  C is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/*           lower triangular part of the array C must contain the lower */
/*           triangular part  of the  symmetric matrix  and the strictly */
/*           upper triangular part of C is not referenced.  On exit, the */
/*           lower triangular part of the array  C is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */


/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    info = 0;
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,nrowa)) {
	info = 9;
    } else if (*ldc < max(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("SSYR2K", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == (real) 0.0 || *k == 0) && *beta == (real) 1.0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == (real) 0.0) {
	if (upper) {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B' + alpha*B*A' + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.0) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L90: */
		    }
		} else if (*beta != (real) 1.0) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != (real) 0.0 || b[j + l * b_dim1] != (real) 0.0) 
			    {
			temp1 = *alpha * b[j + l * b_dim1];
			temp2 = *alpha * a[j + l * a_dim1];
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.0) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L140: */
		    }
		} else if (*beta != (real) 1.0) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != (real) 0.0 || b[j + l * b_dim1] != (real) 0.0) 
			    {
			temp1 = *alpha * b[j + l * b_dim1];
			temp2 = *alpha * a[j + l * a_dim1];
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*B + alpha*B'*A + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = (real) 0.0;
		    temp2 = (real) 0.0;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L190: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp1 = (real) 0.0;
		    temp2 = (real) 0.0;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L220: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of SSYR2K. */

} /* ssyr2k_ */

/* ssyr.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asyr_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYR   performs the symmetric rank 1 operation */

/*     A := alpha*x*x' + A, */

/*  where alpha is a real scalar, x is an n element vector and A is an */
/*  n by n symmetric matrix. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. On exit, the */
/*           upper triangular part of the array A is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. On exit, the */
/*           lower triangular part of the array A is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*lda < max(1,*n)) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("SSYR  ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == (real) 0.0) {
	return 0;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in upper triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0) {
		    temp = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[i__] * temp;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    ix = kx;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[ix] * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in lower triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != (real) 0.0) {
		    temp = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[i__] * temp;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != (real) 0.0) {
		    temp = *alpha * x[jx];
		    ix = jx;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[ix] * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of SSYR  . */

} /* ssyr_ */

/* ssyrk.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int asyrk_(char *uplo, char *trans, integer *n, integer *k, 
	real *alpha, real *a, integer *lda, real *beta, real *c__, integer *
	ldc, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYRK  performs one of the symmetric rank k operations */

/*     C := alpha*A*A' + beta*C, */

/*  or */

/*     C := alpha*A'*A + beta*C, */

/*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix */
/*  in the second case. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/*           triangular  part  of the  array  C  is to be  referenced  as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry,  TRANS  specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C. */

/*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C. */

/*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N specifies the order of the matrix C.  N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/*           of  columns   of  the   matrix   A,   and  on   entry   with */
/*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/*           of rows of the matrix  A.  K must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is */
/*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/*           part of the array  A  must contain the matrix  A,  otherwise */
/*           the leading  k by n  part of the array  A  must contain  the */
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/*           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/*           be at least  max( 1, k ). */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. */
/*           Unchanged on exit. */

/*  C      - REAL             array of DIMENSION ( LDC, n ). */
/*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/*           upper triangular part of the array C must contain the upper */
/*           triangular part  of the  symmetric matrix  and the strictly */
/*           lower triangular part of C is not referenced.  On exit, the */
/*           upper triangular part of the array  C is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/*           lower triangular part of the array C must contain the lower */
/*           triangular part  of the  symmetric matrix  and the strictly */
/*           upper triangular part of C is not referenced.  On exit, the */
/*           lower triangular part of the array  C is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    info = 0;
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldc < max(1,*n)) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("SSYRK ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == (real) 0.0 || *k == 0) && *beta == (real) 1.0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == (real) 0.0) {
	if (upper) {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == (real) 0.0) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A' + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.0) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L90: */
		    }
		} else if (*beta != (real) 1.0) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != (real) 0.0) {
			temp = *alpha * a[j + l * a_dim1];
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == (real) 0.0) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = (real) 0.0;
/* L140: */
		    }
		} else if (*beta != (real) 1.0) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != (real) 0.0) {
			temp = *alpha * a[j + l * a_dim1];
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*A + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = (real) 0.0;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
/* L190: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp = (real) 0.0;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
/* L220: */
		    }
		    if (*beta == (real) 0.0) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of SSYRK . */

} /* ssyrk_ */

/* stbmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atbmv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, real *a, integer *lda, real *x, integer *incx, ftnlen 
	uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STBMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, */
/*  upper or lower triangular band matrix, with ( k + 1 ) diagonals. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry with UPLO = 'U' or 'u', K specifies the number of */
/*           super-diagonals of the matrix A. */
/*           On entry with UPLO = 'L' or 'l', K specifies the number of */
/*           sub-diagonals of the matrix A. */
/*           K must satisfy  0 .le. K. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/*           by n part of the array A must contain the upper triangular */
/*           band part of the matrix of coefficients, supplied column by */
/*           column, with the leading diagonal of the matrix in row */
/*           ( k + 1 ) of the array, the first super-diagonal starting at */
/*           position 2 in row k, and so on. The top left k by k triangle */
/*           of the array A is not referenced. */
/*           The following program segment will transfer an upper */
/*           triangular band matrix from conventional full matrix storage */
/*           to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = K + 1 - J */
/*                    DO 10, I = MAX( 1, J - K ), J */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/*           by n part of the array A must contain the lower triangular */
/*           band part of the matrix of coefficients, supplied column by */
/*           column, with the leading diagonal of the matrix in row 1 of */
/*           the array, the first sub-diagonal starting at position 1 in */
/*           row 2, and so on. The bottom right k by k triangle of the */
/*           array A is not referenced. */
/*           The following program segment will transfer a lower */
/*           triangular band matrix from conventional full matrix storage */
/*           to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = 1 - J */
/*                    DO 10, I = J, MIN( N, J + K ) */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Note that when DIAG = 'U' or 'u' the elements of the array A */
/*           corresponding to the diagonal elements of the matrix are not */
/*           referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           ( k + 1 ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < *k + 1) {
	info = 7;
    } else if (*incx == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("STBMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX   too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*         Form  x := A*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			l = kplus1 - j;
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__4 = j - 1;
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			    x[i__] += temp * a[l + i__ + j * a_dim1];
/* L10: */
			}
			if (nounit) {
			    x[j] *= a[kplus1 + j * a_dim1];
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			l = kplus1 - j;
/* Computing MAX */
			i__4 = 1, i__2 = j - *k;
			i__3 = j - 1;
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			    x[ix] += temp * a[l + i__ + j * a_dim1];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a[kplus1 + j * a_dim1];
			}
		    }
		    jx += *incx;
		    if (j > *k) {
			kx += *incx;
		    }
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			l = 1 - j;
/* Computing MIN */
			i__1 = *n, i__3 = j + *k;
			i__4 = j + 1;
			for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
			    x[i__] += temp * a[l + i__ + j * a_dim1];
/* L50: */
			}
			if (nounit) {
			    x[j] *= a[j * a_dim1 + 1];
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			l = 1 - j;
/* Computing MIN */
			i__4 = *n, i__1 = j + *k;
			i__3 = j + 1;
			for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
			    x[ix] += temp * a[l + i__ + j * a_dim1];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a[j * a_dim1 + 1];
			}
		    }
		    jx -= *incx;
		    if (*n - j >= *k) {
			kx -= *incx;
		    }
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    l = kplus1 - j;
		    if (nounit) {
			temp *= a[kplus1 + j * a_dim1];
		    }
/* Computing MAX */
		    i__4 = 1, i__1 = j - *k;
		    i__3 = max(i__4,i__1);
		    for (i__ = j - 1; i__ >= i__3; --i__) {
			temp += a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    kx -= *incx;
		    ix = kx;
		    l = kplus1 - j;
		    if (nounit) {
			temp *= a[kplus1 + j * a_dim1];
		    }
/* Computing MAX */
		    i__4 = 1, i__1 = j - *k;
		    i__3 = max(i__4,i__1);
		    for (i__ = j - 1; i__ >= i__3; --i__) {
			temp += a[l + i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    temp = x[j];
		    l = 1 - j;
		    if (nounit) {
			temp *= a[j * a_dim1 + 1];
		    }
/* Computing MIN */
		    i__1 = *n, i__2 = j + *k;
		    i__4 = min(i__1,i__2);
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
			temp += a[l + i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    temp = x[jx];
		    kx += *incx;
		    ix = kx;
		    l = 1 - j;
		    if (nounit) {
			temp *= a[j * a_dim1 + 1];
		    }
/* Computing MIN */
		    i__1 = *n, i__2 = j + *k;
		    i__4 = min(i__1,i__2);
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
			temp += a[l + i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STBMV . */

} /* stbmv_ */

/* stbsv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int stbsv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, real *a, integer *lda, real *x, integer *incx, ftnlen 
	uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STBSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular band matrix, with ( k + 1 ) */
/*  diagonals. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry with UPLO = 'U' or 'u', K specifies the number of */
/*           super-diagonals of the matrix A. */
/*           On entry with UPLO = 'L' or 'l', K specifies the number of */
/*           sub-diagonals of the matrix A. */
/*           K must satisfy  0 .le. K. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/*           by n part of the array A must contain the upper triangular */
/*           band part of the matrix of coefficients, supplied column by */
/*           column, with the leading diagonal of the matrix in row */
/*           ( k + 1 ) of the array, the first super-diagonal starting at */
/*           position 2 in row k, and so on. The top left k by k triangle */
/*           of the array A is not referenced. */
/*           The following program segment will transfer an upper */
/*           triangular band matrix from conventional full matrix storage */
/*           to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = K + 1 - J */
/*                    DO 10, I = MAX( 1, J - K ), J */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/*           by n part of the array A must contain the lower triangular */
/*           band part of the matrix of coefficients, supplied column by */
/*           column, with the leading diagonal of the matrix in row 1 of */
/*           the array, the first sub-diagonal starting at position 1 in */
/*           row 2, and so on. The bottom right k by k triangle of the */
/*           array A is not referenced. */
/*           The following program segment will transfer a lower */
/*           triangular band matrix from conventional full matrix storage */
/*           to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = 1 - J */
/*                    DO 10, I = J, MIN( N, J + K ) */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Note that when DIAG = 'U' or 'u' the elements of the array A */
/*           corresponding to the diagonal elements of the matrix are not */
/*           referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           ( k + 1 ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < *k + 1) {
	info = 7;
    } else if (*incx == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("STBSV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed by sequentially with one pass through A. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.0) {
			l = kplus1 - j;
			if (nounit) {
			    x[j] /= a[kplus1 + j * a_dim1];
			}
			temp = x[j];
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i__ = j - 1; i__ >= i__1; --i__) {
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    kx -= *incx;
		    if (x[jx] != (real) 0.0) {
			ix = kx;
			l = kplus1 - j;
			if (nounit) {
			    x[jx] /= a[kplus1 + j * a_dim1];
			}
			temp = x[jx];
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i__ = j - 1; i__ >= i__1; --i__) {
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
			    ix -= *incx;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.0) {
			l = 1 - j;
			if (nounit) {
			    x[j] /= a[j * a_dim1 + 1];
			}
			temp = x[j];
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    kx += *incx;
		    if (x[jx] != (real) 0.0) {
			ix = kx;
			l = 1 - j;
			if (nounit) {
			    x[jx] /= a[j * a_dim1 + 1];
			}
			temp = x[jx];
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
			    ix += *incx;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A')*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    l = kplus1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *k;
		    i__4 = j - 1;
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			temp -= a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    if (nounit) {
			temp /= a[kplus1 + j * a_dim1];
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    l = kplus1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *k;
		    i__3 = j - 1;
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			temp -= a[l + i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= a[kplus1 + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx += *incx;
		    if (j > *k) {
			kx += *incx;
		    }
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    l = 1 - j;
/* Computing MIN */
		    i__1 = *n, i__3 = j + *k;
		    i__4 = j + 1;
		    for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
			temp -= a[l + i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    if (nounit) {
			temp /= a[j * a_dim1 + 1];
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__1 = j + *k;
		    i__3 = j + 1;
		    for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
			temp -= a[l + i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= a[j * a_dim1 + 1];
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    if (*n - j >= *k) {
			kx -= *incx;
		    }
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STBSV . */

} /* stbsv_ */

/* stpmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atpmv_(char *uplo, char *trans, char *diag, integer *n, 
	real *ap, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STPMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, */
/*  upper or lower triangular matrix, supplied in packed form. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - REAL             array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("STPMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x:= A*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			k = kk;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * ap[k];
			    ++k;
/* L10: */
			}
			if (nounit) {
			    x[j] *= ap[kk + j - 1];
			}
		    }
		    kk += j;
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			i__2 = kk + j - 2;
			for (k = kk; k <= i__2; ++k) {
			    x[ix] += temp * ap[k];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= ap[kk + j - 1];
			}
		    }
		    jx += *incx;
		    kk += j;
/* L40: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			k = kk;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * ap[k];
			    --k;
/* L50: */
			}
			if (nounit) {
			    x[j] *= ap[kk - *n + j];
			}
		    }
		    kk -= *n - j + 1;
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			i__1 = kk - (*n - (j + 1));
			for (k = kk; k >= i__1; --k) {
			    x[ix] += temp * ap[k];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= ap[kk - *n + j];
			}
		    }
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= ap[kk];
		    }
		    k = kk - 1;
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += ap[k] * x[i__];
			--k;
/* L90: */
		    }
		    x[j] = temp;
		    kk -= j;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= ap[kk];
		    }
		    i__1 = kk - j + 1;
		    for (k = kk - 1; k >= i__1; --k) {
			ix -= *incx;
			temp += ap[k] * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    kk -= j;
/* L120: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= ap[kk];
		    }
		    k = kk + 1;
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += ap[k] * x[i__];
			++k;
/* L130: */
		    }
		    x[j] = temp;
		    kk += *n - j + 1;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= ap[kk];
		    }
		    i__2 = kk + *n - j;
		    for (k = kk + 1; k <= i__2; ++k) {
			ix += *incx;
			temp += ap[k] * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
		    kk += *n - j + 1;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STPMV . */

} /* stpmv_ */

/* stpsv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atpsv_(char *uplo, char *trans, char *diag, integer *n, 
	real *ap, real *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STPSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix, supplied in packed form. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - REAL             array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("STPSV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.0) {
			if (nounit) {
			    x[j] /= ap[kk];
			}
			temp = x[j];
			k = kk - 1;
			for (i__ = j - 1; i__ >= 1; --i__) {
			    x[i__] -= temp * ap[k];
			    --k;
/* L10: */
			}
		    }
		    kk -= j;
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != (real) 0.0) {
			if (nounit) {
			    x[jx] /= ap[kk];
			}
			temp = x[jx];
			ix = jx;
			i__1 = kk - j + 1;
			for (k = kk - 1; k >= i__1; --k) {
			    ix -= *incx;
			    x[ix] -= temp * ap[k];
/* L30: */
			}
		    }
		    jx -= *incx;
		    kk -= j;
/* L40: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.0) {
			if (nounit) {
			    x[j] /= ap[kk];
			}
			temp = x[j];
			k = kk + 1;
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * ap[k];
			    ++k;
/* L50: */
			}
		    }
		    kk += *n - j + 1;
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != (real) 0.0) {
			if (nounit) {
			    x[jx] /= ap[kk];
			}
			temp = x[jx];
			ix = jx;
			i__2 = kk + *n - j;
			for (k = kk + 1; k <= i__2; ++k) {
			    ix += *incx;
			    x[ix] -= temp * ap[k];
/* L70: */
			}
		    }
		    jx += *incx;
		    kk += *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    k = kk;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= ap[k] * x[i__];
			++k;
/* L90: */
		    }
		    if (nounit) {
			temp /= ap[kk + j - 1];
		    }
		    x[j] = temp;
		    kk += j;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    i__2 = kk + j - 2;
		    for (k = kk; k <= i__2; ++k) {
			temp -= ap[k] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= ap[kk + j - 1];
		    }
		    x[jx] = temp;
		    jx += *incx;
		    kk += j;
/* L120: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    k = kk;
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= ap[k] * x[i__];
			--k;
/* L130: */
		    }
		    if (nounit) {
			temp /= ap[kk - *n + j];
		    }
		    x[j] = temp;
		    kk -= *n - j + 1;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    i__1 = kk - (*n - (j + 1));
		    for (k = kk; k >= i__1; --k) {
			temp -= ap[k] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= ap[kk - *n + j];
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STPSV . */

} /* stpsv_ */

/* strmm.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atrmm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, 
	integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, info;
    static real temp;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRMM  performs one of the matrix-matrix operations */

/*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */

/*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or */
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */

/*     op( A ) = A   or   op( A ) = A'. */

/*  Arguments */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry,  SIDE specifies whether  op( A ) multiplies B from */
/*           the left or right as follows: */

/*              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */

/*              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular */
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at */
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/*           zero then  A is not referenced and  B need not be set before */
/*           entry. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m */
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/*           upper triangular part of the array  A must contain the upper */
/*           triangular matrix  and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/*           lower triangular part of the array  A must contain the lower */
/*           triangular matrix  and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/*           A  are not referenced either,  but are assumed to be  unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - REAL             array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must */
/*           contain the matrix  B,  and  on exit  is overwritten  by the */
/*           transformed matrix. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    info = 0;
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	info = 2;
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 3;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("STRMM ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == (real) 0.0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = (real) 0.0;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*A*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != (real) 0.0) {
			    temp = *alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
/* L30: */
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != (real) 0.0) {
			    temp = *alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = *m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {

/*           Form  B := alpha*A'*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L90: */
			}
			b[i__ + j * b_dim1] = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = *m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L120: */
			}
			b[i__ + j * b_dim1] = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*A. */

	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != (real) 0.0) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != (real) 0.0) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != (real) 0.0) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != (real) 1.0) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != (real) 0.0) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != (real) 1.0) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }

    return 0;

/*     End of STRMM . */

} /* strmm_ */

/* strmv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atrmv_(char *uplo, char *trans, char *diag, integer *n, 
	real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, */
/*  upper or lower triangular matrix. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular matrix and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular matrix and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("STRMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
/* L10: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.0) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != (real) 0.0) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= *incx;
			temp += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += *incx;
			temp += a[i__ + j * a_dim1] * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STRMV . */

} /* strmv_ */

/* strsm.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atrsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, 
	integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, info;
    static real temp;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRSM  solves one of the matrix equations */

/*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */

/*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */

/*     op( A ) = A   or   op( A ) = A'. */

/*  The matrix X is overwritten on B. */

/*  Arguments */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry, SIDE specifies whether op( A ) appears on the left */
/*           or right of X as follows: */

/*              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */

/*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular */
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at */
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/*           zero then  A is not referenced and  B need not be set before */
/*           entry. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m */
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/*           upper triangular part of the array  A must contain the upper */
/*           triangular matrix  and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/*           lower triangular part of the array  A must contain the lower */
/*           triangular matrix  and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/*           A  are not referenced either,  but are assumed to be  unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - REAL             array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must */
/*           contain  the  right-hand  side  matrix  B,  and  on exit  is */
/*           overwritten by the solution matrix  X. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */


/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    info = 0;
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	info = 2;
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 3;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("STRSM ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == (real) 0.0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = (real) 0.0;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*inv( A )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != (real) 1.0) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L30: */
			}
		    }
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != (real) 0.0) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__2 = k - 1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L40: */
			    }
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != (real) 1.0) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L70: */
			}
		    }
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != (real) 0.0) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__3 = *m;
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L80: */
			    }
			}
/* L90: */
		    }
/* L100: */
		}
	    }
	} else {

/*           Form  B := alpha*inv( A' )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__3 = i__ - 1;
			for (k = 1; k <= i__3; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L110: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L120: */
		    }
/* L130: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__2 = *m;
			for (k = i__ + 1; k <= i__2; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L140: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L150: */
		    }
/* L160: */
		}
	    }
	}
    } else {
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*inv( A ). */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != (real) 1.0) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L170: */
			}
		    }
		    i__2 = j - 1;
		    for (k = 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != (real) 0.0) {
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L180: */
			    }
			}
/* L190: */
		    }
		    if (nounit) {
			temp = (real) 1.0 / a[j + j * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L200: */
			}
		    }
/* L210: */
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != (real) 1.0) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L220: */
			}
		    }
		    i__1 = *n;
		    for (k = j + 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != (real) 0.0) {
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    if (nounit) {
			temp = (real) 1.0 / a[j + j * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    }
	} else {

/*           Form  B := alpha*B*inv( A' ). */

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = (real) 1.0 / a[k + k * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L270: */
			}
		    }
		    i__1 = k - 1;
		    for (j = 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != (real) 0.0) {
			    temp = a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L280: */
			    }
			}
/* L290: */
		    }
		    if (*alpha != (real) 1.0) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L300: */
			}
		    }
/* L310: */
		}
	    } else {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    if (nounit) {
			temp = (real) 1.0 / a[k + k * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L320: */
			}
		    }
		    i__2 = *n;
		    for (j = k + 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != (real) 0.0) {
			    temp = a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L330: */
			    }
			}
/* L340: */
		    }
		    if (*alpha != (real) 1.0) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L350: */
			}
		    }
/* L360: */
		}
	    }
	}
    }

    return 0;

/*     End of STRSM . */

} /* strsm_ */

/* strsv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Subroutine */ int atrsv_(char *uplo, char *trans, char *diag, integer *n, 
	real *a, integer *lda, real *x, integer *incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ATRSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Arguments */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular matrix and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular matrix and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
	info = 2;
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("STRSV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != (real) 0.f) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			for (i__ = j - 1; i__ >= 1; --i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != (real) 0.f) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			for (i__ = j - 1; i__ >= 1; --i__) {
			    ix -= *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != (real) 0.f) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != (real) 0.f) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    ix += *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of STRSV . */

} /* strsv_ */

/* lsame.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */

/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */

/*     Test if the characters are equal */

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */

/* xerbla.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int xerbla_(char *srname, integer *info, ftnlen srname_len)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu"
	    "mber \002,i2,\002 had \002,\002an illegal value\002)";

    /* Builtin functions */
//    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */


//    s_wsfe(&io___1);
//    do_fio(&c__1, srname, (ftnlen)6);
//    do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
//    e_wsfe();

//    s_stop("", (ftnlen)0);


/*     End of XERBLA */

    return 0;
} /* xerbla_ */

} // end extern "C"

