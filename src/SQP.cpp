/* SLSQP: Sequentional Least Squares Programming (aka sequential quadratic programming SQP)
   method for nonlinearly constrained nonlinear optimization, by Dieter Kraft (1991).
   Fortran released under a free (BSD) license by ACM to the SciPy project and used there.
   C translation via f2c + hand-cleanup and incorporation into NLopt by S. G. Johnson (2009). 
*/

// Table of constant values 
#include <iostream>
#include <limits>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <Moby/Types.h>
#include <Moby/cblas.h>

using boost::shared_ptr;
using namespace Moby;

//      ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM. 
//      TRANSACTIONS ON MATHEMATICAL SOFTWARE, 
//      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281. 
//      http://doi.acm.org/10.1145/192115.192124 


//      http://permalink.gmane.org/gmane.comp.python.scientific.devel/6725 
//      ------ 
//      From: Deborah Cotton <cotton@hq.acm.org> 
//      Date: Fri, 14 Sep 2007 12:35:55 -0500 
//      Subject: RE: Algorithm License requested 
//      To: Alan Isaac 

//      Prof. Issac, 

//      In that case, then because the author consents to [the ACM] releasing 
//      the code currently archived at http://www.netlib.org/toms/733 under the 
//      BSD license, the ACM hereby releases this code under the BSD license. 

//      Regards, 

//      Deborah Cotton, Copyright & Permissions 
//      ACM Publications 
//      2 Penn Plaza, Suite 701** 
//      New York, NY 10121-0701 
//      permissions@acm.org 
//      212.869.7440 ext. 652 
//      Fax. 212.869.0481 
//      ------ 

static const int c__0 = 0;
static const int c__1 = 1;
static const int c__2 = 2;

static void h12_(const int *mode, int *lpivot, int *l1, 
		 int *m, Real *u, const int *iue, Real *up, 
		 Real *c__, const int *ice, const int *icv, const int *ncv)
{
    // Initialized data 

    const Real one = (Real) 1.;

    // System generated locals 
    int u_dim1, u_offset, i__1, i__2;
    Real d__1;

    // Local variables 
    Real b;
    int i__, j, i2, i3, i4;
    Real cl, sm;
    int incr;
    Real clinv;

//     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 
//     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 
//     CONSTRUCTION AND/OR APPLICATION OF A SINGLE 
//     HOUSEHOLDER TRANSFORMATION  Q = I + U*(U**T)/B 
//     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 . 
//     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
//     L1,M   IF L1 <= M   THE TRANSFORMATION WILL BE CONSTRUCTED TO 
//            ZERO ELEMENTS INDEXED FROM L1 THROUGH M. 
//            IF L1 > M THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. 
//     U(),IUE,UP 
//            ON ENTRY TO H1 U() STORES THE PIVOT VECTOR. 
//            IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS. 
//            ON EXIT FROM H1 U() AND UP STORE QUANTITIES DEFINING 
//            THE VECTOR U OF THE HOUSEHOLDER TRANSFORMATION. 
//            ON ENTRY TO H2 U() AND UP 
//            SHOULD STORE QUANTITIES PREVIOUSLY COMPUTED BY H1. 
//            THESE WILL NOT BE MODIFIED BY H2. 
//     C()    ON ENTRY TO H1 OR H2 C() STORES A MATRIX WHICH WILL BE 
//            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER 
//            TRANSFORMATION IS TO BE APPLIED. 
//            ON EXIT C() STORES THE SET OF TRANSFORMED VECTORS. 
//     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). 
//     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). 
//     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. 
//            IF NCV <= 0 NO OPERATIONS WILL BE DONE ON C(). 
    // Parameter adjustments 
    u_dim1 = *iue;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --c__;

    // Function Body 
    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
	goto L80;
    }
    cl = (d__1 = u[*lpivot * u_dim1 + 1], std::fabs(d__1));
    if (*mode == 2) {
	goto L30;
    }
//     ****** CONSTRUCT THE TRANSFORMATION ****** 
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
	sm = (d__1 = u[j * u_dim1 + 1], std::fabs(d__1));
// L10: 
	cl = std::max(sm,cl);
    }
    if (cl <= 0.0) {
	goto L80;
    }
    clinv = one / cl;
// Computing 2nd power 
    d__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = d__1 * d__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
// L20: 
// Computing 2nd power 
	d__1 = u[j * u_dim1 + 1] * clinv;
	sm += d__1 * d__1;
    }
    cl *= std::sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] > 0.0) {
	cl = -cl;
    }
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L40;
//     ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C ****** 
L30:
    if (cl <= 0.0) {
	goto L80;
    }
L40:
    if (*ncv <= 0) {
	goto L80;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
    if (b >= 0.0) {
	goto L80;
    }
    b = one / b;
    i2 = 1 - *icv + *ice * (*lpivot - 1);
    incr = *ice * (*l1 - *lpivot);
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	i2 += *icv;
	i3 = i2 + incr;
	i4 = i3;
	sm = c__[i2] * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    sm += c__[i3] * u[i__ * u_dim1 + 1];
// L50: 
	    i3 += *ice;
	}
	if (sm == 0.0) {
	    goto L70;
	}
	sm *= b;
	c__[i2] += sm * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    c__[i4] += sm * u[i__ * u_dim1 + 1];
// L60: 
	    i4 += *ice;
	}
L70:
	;
    }
L80:
    return;
} // h12_ 

static void nnls_(Real *a, int *mda, int *m, int *
	n, Real *b, Real *x, Real *rnorm, Real *w, 
	Real *z__, int *indx, int *mode)
{
    // Initialized data 

    const Real one = (Real) 1.;
    const Real factor = (Real) .01;

    // System generated locals 
    int a_dim1, a_offset, i__1, i__2;
    Real d__1;

    // Local variables 
    Real c__;
    int i__, j, k, l;
    Real s, t;
    int ii, jj, ip, iz, jz;
    Real up;
    int iz1, iz2, npp1, iter;
    Real wmax, alpha, asave;
    int itmax, izmax, nsetp;
    Real unorm;

//     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY: 
//     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974 
//      **********   NONNEGATIVE LEAST SQUARES   ********** 
//     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN 
//     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM 
//                  A*X = B  SUBJECT TO  X >= 0 
//     A(),MDA,M,N 
//            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A(). 
//            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A. 
//            ON EXIT A() CONTAINS THE PRODUCT Q*A, 
//            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED 
//            IMPLICITLY BY THIS SUBROUTINE. 
//            EITHER M>=N OR M<N IS PERMISSIBLE. 
//            THERE IS NO RESTRICTION ON THE RANK OF A. 
//     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B. 
//            ON EXIT B() CONTAINS Q*B. 
//     X()    ON ENTRY X() NEED NOT BE INITIALIZED. 
//            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR. 
//     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE 
//            RESIDUAL VECTOR. 
//     W()    AN N-ARRAY OF WORKING SPACE. 
//            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR. 
//            W WILL SATISFY W(I)=0 FOR ALL I IN SET P 
//            AND W(I)<=0 FOR ALL I IN SET Z 
//     Z()    AN M-ARRAY OF WORKING SPACE. 
//     INDX()AN INT WORKING ARRAY OF LENGTH AT LEAST N. 
//            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS 
//            P AND Z AS FOLLOWS: 
//            INDX(1)    THRU INDX(NSETP) = SET P. 
//            INDX(IZ1)  THRU INDX (IZ2)  = SET Z. 
//            IZ1=NSETP + 1 = NPP1, IZ2=N. 
//     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING: 
//            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. 
//            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG, 
//                 EITHER M <= 0 OR N <= 0. 
//            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS. 
    // Parameter adjustments 
    --z__;
    --b;
    --indx;
    --w;
    --x;
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    // Function Body 
//     revised          Dieter Kraft, March 1983 
    *mode = 2;
    if (*m <= 0 || *n <= 0) {
	goto L290;
    }
    *mode = 1;
    iter = 0;
    itmax = *n * 3;
// STEP ONE (INITIALIZE) 
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
// L100: 
	indx[i__] = i__;
    }
    iz1 = 1;
    iz2 = *n;
    nsetp = 0;
    npp1 = 1;
    x[1] = (Real) 0.0;
    CBLAS::copy(*n, &x[1], 0, &x[1], 1);
// STEP TWO (COMPUTE DUAL VARIABLES) 
// .....ENTRY LOOP A 
L110:
    if (iz1 > iz2 || nsetp >= *m) {
	goto L280;
    }
    i__1 = iz2;
    for (iz = iz1; iz <= i__1; ++iz) {
	j = indx[iz];
// L120: 
	i__2 = *m - nsetp;
	w[j] = CBLAS::dot(i__2, &a[npp1 + j * a_dim1], 1, &b[npp1], 1)
		;
    }
// STEP THREE (TEST DUAL VARIABLES) 
L130:
    wmax = (Real) 0.0;
    i__2 = iz2;
    for (iz = iz1; iz <= i__2; ++iz) {
	j = indx[iz];
	if (w[j] <= wmax) {
	    goto L140;
	}
	wmax = w[j];
	izmax = iz;
L140:
	;
    }
// .....EXIT LOOP A 
    if (wmax <= (Real) 0.0) {
	goto L280;
    }
    iz = izmax;
    j = indx[iz];
// STEP FOUR (TEST INDX J FOR LINEAR DEPENDENCY) 
    asave = a[npp1 + j * a_dim1];
    i__2 = npp1 + 1;
    h12_(&c__1, &npp1, &i__2, m, &a[j * a_dim1 + 1], &c__1, &up, &z__[1], &
	    c__1, &c__1, &c__0);
    unorm = CBLAS::nrm2(nsetp, &a[j * a_dim1 + 1], 1);
    t = factor * (d__1 = a[npp1 + j * a_dim1], std::fabs(d__1));
    d__1 = unorm + t;
    if (d__1 - unorm <= (Real) 0.0) {
	goto L150;
    }
    CBLAS::copy(*m, &b[1], 1, &z__[1], 1);
    i__2 = npp1 + 1;
    h12_(&c__2, &npp1, &i__2, m, &a[j * a_dim1 + 1], &c__1, &up, &z__[1], &
	    c__1, &c__1, &c__1);
    if (z__[npp1] / a[npp1 + j * a_dim1] > (Real) 0.0) {
	goto L160;
    }
L150:
    a[npp1 + j * a_dim1] = asave;
    w[j] = (Real) 0.0;
    goto L130;
// STEP FIVE (ADD COLUMN) 
L160:
    CBLAS::copy(*m, &z__[1], 1, &b[1], 1);
    indx[iz] = indx[iz1];
    indx[iz1] = j;
    ++iz1;
    nsetp = npp1;
    ++npp1;
    i__2 = iz2;
    for (jz = iz1; jz <= i__2; ++jz) {
	jj = indx[jz];
// L170: 
	h12_(&c__2, &nsetp, &npp1, m, &a[j * a_dim1 + 1], &c__1, &up, &a[jj * 
		a_dim1 + 1], &c__1, mda, &c__1);
    }
    k = std::min(npp1,*mda);
    w[j] = (Real) 0.0;
    i__2 = *m - nsetp;
    CBLAS::copy(i__2, &w[j], 0, &a[k + j * a_dim1], 1);
// STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM) 
// .....ENTRY LOOP B 
L180:
    for (ip = nsetp; ip >= 1; --ip) {
	if (ip == nsetp) {
	    goto L190;
	}
	d__1 = -z__[ip + 1];
	CBLAS::axpy(ip, d__1, &a[jj * a_dim1 + 1], 1, &z__[1], 1);
L190:
	jj = indx[ip];
// L200: 
	z__[ip] /= a[ip + jj * a_dim1];
    }
    ++iter;
    if (iter <= itmax) {
	goto L220;
    }
L210:
    *mode = 3;
    goto L280;
// STEP SEVEN TO TEN (STEP LENGTH ALGORITHM) 
L220:
    alpha = one;
    jj = 0;
    i__2 = nsetp;
    for (ip = 1; ip <= i__2; ++ip) {
	if (z__[ip] > (Real) 0.0) {
	    goto L230;
	}
	l = indx[ip];
	t = -x[l] / (z__[ip] - x[l]);
	if (alpha < t) {
	    goto L230;
	}
	alpha = t;
	jj = ip;
L230:
	;
    }
    i__2 = nsetp;
    for (ip = 1; ip <= i__2; ++ip) {
	l = indx[ip];
// L240: 
	x[l] = (one - alpha) * x[l] + alpha * z__[ip];
    }
// .....EXIT LOOP B 
    if (jj == 0) {
	goto L110;
    }
// STEP ELEVEN (DELETE COLUMN) 
    i__ = indx[jj];
L250:
    x[i__] = (Real) 0.0;
    ++jj;
    i__2 = nsetp;
    for (j = jj; j <= i__2; ++j) {
	ii = indx[j];
	indx[j - 1] = ii;
	CBLAS::rotg(a[j - 1 + ii * a_dim1], a[j + ii * a_dim1], c__, s);
	t = a[j - 1 + ii * a_dim1];
	CBLAS::rot(*n, &a[j - 1 + a_dim1], *mda, &a[j + a_dim1], *mda, c__, s);
	a[j - 1 + ii * a_dim1] = t;
	a[j + ii * a_dim1] = (Real) 0.0;
// L260: 
	CBLAS::rot(1, &b[j - 1], 1, &b[j], 1, c__, s);
    }
    npp1 = nsetp;
    --nsetp;
    --iz1;
    indx[iz1] = i__;
    if (nsetp <= 0) {
	goto L210;
    }
    i__2 = nsetp;
    for (jj = 1; jj <= i__2; ++jj) {
	i__ = indx[jj];
	if (x[i__] <= (Real) 0.0) {
	    goto L250;
	}
// L270: 
    }
    CBLAS::copy(*m, &b[1], 1, &z__[1], 1);
    goto L180;
// STEP TWELVE (SOLUTION) 
L280:
    k = std::min(npp1,*m);
    i__2 = *m - nsetp;
    *rnorm = CBLAS::nrm2(i__2, &b[k], 1);
    if (npp1 > *m) {
	w[1] = (Real) 0.0;
	CBLAS::copy(*n, &w[1], 0, &w[1], 1);
    }
// END OF SUBROUTINE NNLS 
L290:
    return;
} // nnls_ 

static void ldp_(Real *g, int *mg, int *m, int *n, 
	Real *h__, Real *x, Real *xnorm, Real *w, 
	int *indx, int *mode)
{
    // Initialized data 

    const Real one = (Real) 1.;

    // System generated locals 
    int g_dim1, g_offset, i__1, i__2;
    Real d__1;

    // Local variables 
    int i__, j, n1, if__, iw, iy, iz;
    Real fac;
    Real rnorm;
    int iwdual;

//                     T 
//     MINIMIZE   1/2 X X    SUBJECT TO   G * X >= H. 
//       C.L. LAWSON, R.J. HANSON: 'SOLVING LEAST SQUARES PROBLEMS' 
//       PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY, 1974. 
//     PARAMETER DESCRIPTION: 
//     G(),MG,M,N   ON ENTRY G() STORES THE M BY N MATRIX OF 
//                  LINEAR INEQUALITY CONSTRAINTS. G() HAS FIRST 
//                  DIMENSIONING PARAMETER MG 
//     H()          ON ENTRY H() STORES THE M VECTOR H REPRESENTING 
//                  THE RIGHT SIDE OF THE INEQUALITY SYSTEM 
//     REMARK: G(),H() WILL NOT BE CHANGED DURING CALCULATIONS BY LDP 
//     X()          ON ENTRY X() NEED NOT BE INITIALIZED. 
//                  ON EXIT X() STORES THE SOLUTION VECTOR X IF MODE=1. 
//     XNORM        ON EXIT XNORM STORES THE EUCLIDIAN NORM OF THE 
//                  SOLUTION VECTOR IF COMPUTATION IS SUCCESSFUL 
//     W()          W IS A ONE DIMENSIONAL WORKING SPACE, THE LENGTH 
//                  OF WHICH SHOULD BE AT LEAST (M+2)*(N+1) + 2*M 
//                  ON EXIT W() STORES THE LAGRANGE MULTIPLIERS 
//                  ASSOCIATED WITH THE CONSTRAINTS 
//                  AT THE SOLUTION OF PROBLEM LDP 
//     INDX()      INDX() IS A ONE DIMENSIONAL INT WORKING SPACE 
//                  OF LENGTH AT LEAST M 
//     MODE         MODE IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
//                  MEANINGS: 
//          MODE=1: SUCCESSFUL COMPUTATION 
//               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N.LE.0) 
//               3: ITERATION COUNT EXCEEDED BY NNLS 
//               4: INEQUALITY CONSTRAINTS INCOMPATIBLE 
    // Parameter adjustments 
    --indx;
    --h__;
    --x;
    g_dim1 = *mg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --w;

    // Function Body 
    *mode = 2;
    if (*n <= 0) {
	goto L50;
    }
//  STATE DUAL PROBLEM 
    *mode = 1;
    x[1] = (Real) 0.0;
    CBLAS::copy(*n, &x[1], 0, &x[1], 1);
    *xnorm = (Real) 0.0;
    if (*m == 0) {
	goto L50;
    }
    iw = 0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++iw;
// L10: 
	    w[iw] = g[j + i__ * g_dim1];
	}
	++iw;
// L20: 
	w[iw] = h__[j];
    }
    if__ = iw + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++iw;
// L30: 
	w[iw] = (Real) 0.0;
    }
    w[iw + 1] = one;
    n1 = *n + 1;
    iz = iw + 2;
    iy = iz + n1;
    iwdual = iy + *m;
//  SOLVE DUAL PROBLEM 
    nnls_(&w[1], &n1, &n1, m, &w[if__], &w[iy], &rnorm, &w[iwdual], &w[iz], &
	    indx[1], mode);
    if (*mode != 1) {
	goto L50;
    }
    *mode = 4;
    if (rnorm <= (Real) 0.0) {
	goto L50;
    }
//  COMPUTE SOLUTION OF PRIMAL PROBLEM 
    fac = one - CBLAS::dot(*m, &h__[1], 1, &w[iy], 1);
    d__1 = one + fac;
    if (d__1 - one <= (Real) 0.0) {
	goto L50;
    }
    *mode = 1;
    fac = one / fac;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
// L40: 
	x[j] = fac * CBLAS::dot(*m, &g[j * g_dim1 + 1], 1, &w[iy], 1);
    }
    *xnorm = CBLAS::nrm2(*n, &x[1], 1);
//  COMPUTE LAGRANGE MULTIPLIERS FOR PRIMAL PROBLEM 
    w[1] = (Real) 0.0;
    CBLAS::copy(*m, &w[1], 0, &w[1], 1);
    CBLAS::axpy(*m, fac, &w[iy], 1, &w[1], 1);
//  END OF SUBROUTINE LDP 
L50:
    return;
} // ldp_ 

static void lsi_(Real *e, Real *f, Real *g, 
	Real *h__, int *le, int *me, int *lg, int *mg, 
	int *n, Real *x, Real *xnorm, Real *w, int *
	jw, int *mode)
{
    // Initialized data 

    const Real epmach = std::numeric_limits<Real>::epsilon();
    const Real one = (Real) 1.;

    // System generated locals 
    int e_dim1, e_offset, g_dim1, g_offset, i__1, i__2, i__3;
    Real d__1;

    // Local variables 
    int i__, j;
    Real t;

//     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF 
//     INEQUALITY CONSTRAINED LINEAR LEAST SQUARES PROBLEM: 
//                    MIN ||E*X-F|| 
//                     X 
//                    S.T.  G*X >= H 
//     THE ALGORITHM IS BASED ON QR DECOMPOSITION AS DESCRIBED IN 
//     CHAPTER 23.5 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS 
//     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM 
//     ARE NECESSARY 
//     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N) 
//     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  ) 
//     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N) 
//     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  ) 
//     DIM(X) :   N 
//     DIM(W) :   (N+1)*(MG+2) + 2*MG 
//     DIM(JW):   LG 
//     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS E, F, G, AND H. 
//     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE. 
//     X     STORES THE SOLUTION VECTOR 
//     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM 
//     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST 
//           MG ELEMENTS 
//     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: 
//          MODE=1: SUCCESSFUL COMPUTATION 
//               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) 
//               3: ITERATION COUNT EXCEEDED BY NNLS 
//               4: INEQUALITY CONSTRAINTS INCOMPATIBLE 
//               5: MATRIX E IS NOT OF FULL RANK 
//     03.01.1980, DIETER KRAFT: CODED 
//     20.03.1987, DIETER KRAFT: REVISED TO FORTRAN 77 
    // Parameter adjustments 
    --f;
    --jw;
    --h__;
    --x;
    g_dim1 = *lg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    e_dim1 = *le;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    --w;

    // Function Body 
//  QR-FACTORS OF E AND APPLICATION TO F 
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
// Computing MIN 
	i__2 = i__ + 1;
	j = std::min(i__2,*n);
	i__2 = i__ + 1;
	i__3 = *n - i__;
	h12_(&c__1, &i__, &i__2, me, &e[i__ * e_dim1 + 1], &c__1, &t, &e[j * 
		e_dim1 + 1], &c__1, le, &i__3);
// L10: 
	i__2 = i__ + 1;
	h12_(&c__2, &i__, &i__2, me, &e[i__ * e_dim1 + 1], &c__1, &t, &f[1], &
		c__1, &c__1, &c__1);
    }
//  TRANSFORM G AND H TO GET LEAST DISTANCE PROBLEM 
    *mode = 5;
    i__2 = *mg;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if ((d__1 = e[j + j * e_dim1], std::fabs(d__1)) < epmach) {
		goto L50;
	    }
// L20: 
	    i__3 = j - 1;
	    g[i__ + j * g_dim1] = (g[i__ + j * g_dim1] - CBLAS::dot(i__3, &g[
		    i__ + g_dim1], *lg, &e[j * e_dim1 + 1], 1)) / e[j + j *
		     e_dim1];
	}
// L30: 
	h__[i__] -= CBLAS::dot(*n, &g[i__ + g_dim1], *lg, &f[1], 1);
    }
//  SOLVE LEAST DISTANCE PROBLEM 
    ldp_(&g[g_offset], lg, mg, n, &h__[1], &x[1], xnorm, &w[1], &jw[1], mode);
    if (*mode != 1) {
	goto L50;
    }
//  SOLUTION OF ORIGINAL PROBLEM 
    CBLAS::axpy(*n, one, &f[1], 1, &x[1], 1);
    for (i__ = *n; i__ >= 1; --i__) {
// Computing MIN 
	i__2 = i__ + 1;
	j = std::min(i__2,*n);
// L40: 
	i__2 = *n - i__;
	x[i__] = (x[i__] - CBLAS::dot(i__2, &e[i__ + j * e_dim1], *le, &x[j], 1))
	     / e[i__ + i__ * e_dim1];
    }
// Computing MIN 
    i__2 = *n + 1;
    j = std::min(i__2,*me);
    i__2 = *me - *n;
    t = CBLAS::nrm2(i__2, &f[j], 1);
    *xnorm = std::sqrt(*xnorm * *xnorm + t * t);
//  END OF SUBROUTINE LSI 
L50:
    return;
} // lsi_ 

static void hfti_(Real *a, int *mda, int *m, int *
	n, Real *b, int *mdb, const int *nb, Real *tau, int 
	*krank, Real *rnorm, Real *h__, Real *g, int *
	ip)
{
    // Initialized data 

    const Real factor = (Real) .001;

    // System generated locals 
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    Real d__1;

    // Local variables 
    int i__, j, k, l;
    int jb, kp1;
    Real tmp, hmax;
    int lmax, ldiag;

//     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN: 
//     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 
//     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 
//     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A 
//                      OF THE LEAST SQUARES PROBLEM AX = B. 
//                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY 
//                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED. 
//                      THERE IS NO RESTRICTION ON THE RANK OF A. 
//                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE. 
//     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE 
//                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST 
//                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE 
//                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN 
//                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X. 
//                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED 
//                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N), 
//                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR 
//                      DOUBLE SUBSCRIPTED. 
//     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK 
//                      DETERMINATION, PROVIDED BY THE USER. 
//     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE. 
//     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDIAN 
//                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM 
//                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B. 
//     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N. 
//     IP()             INT ARRAY OF WORKING SPACE OF LENGTH >= N 
//                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS 
    // Parameter adjustments 
    --ip;
    --g;
    --h__;
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --rnorm;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    // Function Body 
    k = 0;
    ldiag = std::min(*m,*n);
    if (ldiag <= 0) {
	goto L270;
    }
//   COMPUTE LMAX 
    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
	if (j == 1) {
	    goto L20;
	}
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
// Computing 2nd power 
	    d__1 = a[j - 1 + l * a_dim1];
	    h__[l] -= d__1 * d__1;
// L10: 
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
	}
	d__1 = hmax + factor * h__[lmax];
	if (d__1 - hmax > (Real) 0.0) {
	    goto L50;
	}
L20:
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
	    h__[l] = (Real) 0.0;
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
// L30: 
// Computing 2nd power 
		d__1 = a[i__ + l * a_dim1];
		h__[l] += d__1 * d__1;
	    }
// L40: 
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
	}
	hmax = h__[lmax];
//   COLUMN INTERCHANGES IF NEEDED 
L50:
	ip[j] = lmax;
	if (ip[j] == j) {
	    goto L70;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tmp = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
// L60: 
	    a[i__ + lmax * a_dim1] = tmp;
	}
	h__[lmax] = h__[j];
//   J-TH TRANSFORMATION AND APPLICATION TO A AND B 
L70:
// Computing MIN 
	i__2 = j + 1;
	i__ = std::min(i__2,*n);
	i__2 = j + 1;
	i__3 = *n - j;
	h12_(&c__1, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &a[i__ *
		 a_dim1 + 1], &c__1, mda, &i__3);
// L80: 
	i__2 = j + 1;
	h12_(&c__2, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &b[
		b_offset], &c__1, mdb, nb);
    }
//   DETERMINE PSEUDORANK 
    i__2 = ldiag;
    for (j = 1; j <= i__2; ++j) {
// L90: 
	if ((d__1 = a[j + j * a_dim1], std::fabs(d__1)) <= *tau) {
	    goto L100;
	}
    }
    k = ldiag;
    goto L110;
L100:
    k = j - 1;
L110:
    kp1 = k + 1;
//   NORM OF RESIDUALS 
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
// L130: 
	i__1 = *m - k;
	rnorm[jb] = CBLAS::nrm2(i__1, &b[kp1 + jb * b_dim1], 1);
    }
    if (k > 0) {
	goto L160;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
// L150: 
	    b[i__ + jb * b_dim1] = (Real) 0.0;
	}
    }
    goto L270;
L160:
    if (k == *n) {
	goto L180;
    }
//   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS 
    for (i__ = k; i__ >= 1; --i__) {
// L170: 
	i__2 = i__ - 1;
	h12_(&c__1, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &a[
		a_offset], mda, &c__1, &i__2);
    }
L180:
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
//   SOLVE K*K TRIANGULAR SYSTEM 
	for (i__ = k; i__ >= 1; --i__) {
// Computing MIN 
	    i__1 = i__ + 1;
	    j = std::min(i__1,*n);
// L210: 
	    i__1 = k - i__;
	    b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1] - CBLAS::dot(i__1, &
		    a[i__ + j * a_dim1], *mda, &b[j + jb * b_dim1], 1)) / 
		    a[i__ + i__ * a_dim1];
	}
//   COMPLETE SOLUTION VECTOR 
	if (k == *n) {
	    goto L240;
	}
	i__1 = *n;
	for (j = kp1; j <= i__1; ++j) {
// L220: 
	    b[j + jb * b_dim1] = (Real) 0.0;
	}
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
// L230: 
	    h12_(&c__2, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &b[jb *
		     b_dim1 + 1], &c__1, mdb, &c__1);
	}
//   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES 
L240:
	for (j = ldiag; j >= 1; --j) {
	    if (ip[j] == j) {
		goto L250;
	    }
	    l = ip[j];
	    tmp = b[l + jb * b_dim1];
	    b[l + jb * b_dim1] = b[j + jb * b_dim1];
	    b[j + jb * b_dim1] = tmp;
L250:
	    ;
	}
    }
L270:
    *krank = k;
} // hfti_ 

static void lsei_(Real *c__, Real *d__, Real *e, 
	Real *f, Real *g, Real *h__, int *lc, int *
	mc, int *le, int *me, int *lg, int *mg, int *n, 
	Real *x, Real *xnrm, Real *w, int *jw, int *
	mode)
{
    // Initialized data 

    const Real epmach = std::numeric_limits<Real>::epsilon();

    // System generated locals 
    int c_dim1, c_offset, e_dim1, e_offset, g_dim1, g_offset, i__1, i__2, 
	    i__3;
    Real d__1;

    // Local variables 
    int i__, j, k, l;
    Real t;
    int ie, if__, ig, iw, mc1;
    int krank;

//     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF 
//     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI : 
//                MIN ||E*X - F|| 
//                 X 
//                S.T.  C*X  = D, 
//                      G*X >= H. 
//     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C 
//     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS. 
//     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM 
//     ARE NECESSARY 
//     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N) 
//     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  ) 
//     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N) 
//     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  ) 
//     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N) 
//     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  ) 
//     DIM(X) :   FORMAL (N   ),    ACTUAL (N   ) 
//     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI 
//              +(N-MC+1)*(MG+2)+2*MG     for LSI 
//     DIM(JW):   MAX(MG,L) 
//     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H. 
//     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE. 
//     X     STORES THE SOLUTION VECTOR 
//     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM 
//     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST 
//           MC+MG ELEMENTS 
//     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: 
//          MODE=1: SUCCESSFUL COMPUTATION 
//               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) 
//               3: ITERATION COUNT EXCEEDED BY NNLS 
//               4: INEQUALITY CONSTRAINTS INCOMPATIBLE 
//               5: MATRIX E IS NOT OF FULL RANK 
//               6: MATRIX C IS NOT OF FULL RANK 
//               7: RANK DEFECT IN HFTI 
//     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN 
//     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN 
    // Parameter adjustments 
    --d__;
    --f;
    --h__;
    --x;
    g_dim1 = *lg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    e_dim1 = *le;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    c_dim1 = *lc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --w;
    --jw;

    // Function Body 
    *mode = 2;
    if (*mc > *n) {
	goto L75;
    }
    l = *n - *mc;
    mc1 = *mc + 1;
    iw = (l + 1) * (*mg + 2) + (*mg << 1) + *mc;
    ie = iw + *mc + 1;
    if__ = ie + *me * l;
    ig = if__ + *me;
//  TRIANGULARIZE C AND APPLY FACTORS TO E AND G 
    i__1 = *mc;
    for (i__ = 1; i__ <= i__1; ++i__) {
// Computing MIN 
	i__2 = i__ + 1;
	j = std::min(i__2,*lc);
	i__2 = i__ + 1;
	i__3 = *mc - i__;
	h12_(&c__1, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &
		c__[j + c_dim1], lc, &c__1, &i__3);
	i__2 = i__ + 1;
	h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &e[
		e_offset], le, &c__1, me);
// L10: 
	i__2 = i__ + 1;
	h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &g[
		g_offset], lg, &c__1, mg);
    }
//  SOLVE C*X=D AND MODIFY F 
    *mode = 6;
    i__2 = *mc;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if ((d__1 = c__[i__ + i__ * c_dim1], std::fabs(d__1)) < epmach) {
	    goto L75;
	}
	i__1 = i__ - 1;
	x[i__] = (d__[i__] - CBLAS::dot(i__1, &c__[i__ + c_dim1], *lc, &x[1], 1)) 
	     / c__[i__ + i__ * c_dim1];
// L15: 
    }
    *mode = 1;
    w[mc1] = 0.0;
    i__2 = *mg; // BUGFIX for *mc == *n: changed from *mg - *mc, SGJ 2010
    CBLAS::copy(i__2, &w[mc1], 0, &w[mc1], 1);
    if (*mc == *n) {
	goto L50;
    }
    i__2 = *me;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L20: 
	w[if__ - 1 + i__] = f[i__] - CBLAS::dot(*mc, &e[i__ + e_dim1], *le, &x[1], 1);
    }
//  STORE TRANSFORMED E & G 
    i__2 = *me;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L25: 
	CBLAS::copy(l, &e[i__ + mc1 * e_dim1], *le, &w[ie - 1 + i__], *me);
    }
    i__2 = *mg;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L30: 
	CBLAS::copy(l, &g[i__ + mc1 * g_dim1], *lg, &w[ig - 1 + i__], *mg);
    }
    if (*mg > 0) {
	goto L40;
    }
//  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS 
    *mode = 7;
    k = std::max(*le,*n);
    t = std::sqrt(epmach);
    hfti_(&w[ie], me, me, &l, &w[if__], &k, &c__1, &t, &krank, xnrm, &w[1], &
	    w[l + 1], &jw[1]);
    CBLAS::copy(l, &w[if__], 1, &x[mc1], 1);
    if (krank != l) {
	goto L75;
    }
    *mode = 1;
    goto L50;
//  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM 
L40:
    i__2 = *mg;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L45: 
	h__[i__] -= CBLAS::dot(*mc, &g[i__ + g_dim1], *lg, &x[1], 1);
    }
    lsi_(&w[ie], &w[if__], &w[ig], &h__[1], me, me, mg, mg, &l, &x[mc1], xnrm,
	     &w[mc1], &jw[1], mode);
    if (*mc == 0) {
	goto L75;
    }
    t = CBLAS::nrm2(*mc, &x[1], 1);
    *xnrm = std::sqrt(*xnrm * *xnrm + t * t);
    if (*mode != 1) {
	goto L75;
    }
//  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS 
L50:
    i__2 = *me;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L55: 
	f[i__] = CBLAS::dot(*n, &e[i__ + e_dim1], *le, &x[1], 1) - f[i__];
    }
    i__2 = *mc;
    for (i__ = 1; i__ <= i__2; ++i__) {
// L60: 
	d__[i__] = CBLAS::dot(*me, &e[i__ * e_dim1 + 1], 1, &f[1], 1) - 
		CBLAS::dot(*mg, &g[i__ * g_dim1 + 1], 1, &w[mc1], 1);
    }
    for (i__ = *mc; i__ >= 1; --i__) {
// L65: 
	i__2 = i__ + 1;
	h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &x[
		1], &c__1, &c__1, &c__1);
    }
    for (i__ = *mc; i__ >= 1; --i__) {
// Computing MIN 
	i__2 = i__ + 1;
	j = std::min(i__2,*lc);
	i__2 = *mc - i__;
	w[i__] = (d__[i__] - CBLAS::dot(i__2, &c__[j + i__ * c_dim1], 1, &
		w[j], 1)) / c__[i__ + i__ * c_dim1];
// L70: 
    }
//  END OF SUBROUTINE LSEI 
L75:
    return;
} // lsei_ 

static void lsq_(int *m, int *meq, int *n, int *nl, 
	int *la, Real *l, Real *g, Real *a, Real *
	b, const Real *xl, const Real *xu, Real *x, Real *y, 
	Real *w, int *jw, int *mode)
{
    // Initialized data 

    const Real one = (Real) 1.;

    // System generated locals 
    int a_dim1, a_offset, i__1, i__2;
    Real d__1;

    // Local variables 
    int i__, i1, i2, i3, i4, m1, n1, n2, n3, ic, id, ie, if__, ig, ih, il,
	     im, ip, iu, iw;
    Real diag;
    int mineq;
    Real xnorm;

//   MINIMIZE with respect to X 
//             ||E*X - F|| 
//                                      1/2  T 
//   WITH UPPER TRIANGULAR MATRIX E = +D   *L , 
//                                      -1/2  -1 
//                     AND VECTOR F = -D    *L  *G, 
//  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE 
//  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS 
// 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L 
//   SUBJECT TO 
//             A(J)*X - B(J) = 0 ,         J=1,...,MEQ, 
//             A(J)*X - B(J) >=0,          J=MEQ+1,...,M, 
//             XL(I) <= X(I) <= XU(I),     I=1,...,N, 
//     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU. 
//     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N) 
//     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION: 
//     DIM(W) =        (3*N+M)*(N+1)                        for LSQ 
//                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI 
//                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI 
//                      with MINEQ = M - MEQ + 2*N 
//     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE. 
//     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR 
//     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION 
//           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS) 
//     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: 
//          MODE=1: SUCCESSFUL COMPUTATION 
//               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) 
//               3: ITERATION COUNT EXCEEDED BY NNLS 
//               4: INEQUALITY CONSTRAINTS INCOMPATIBLE 
//               5: MATRIX E IS NOT OF FULL RANK 
//               6: MATRIX C IS NOT OF FULL RANK 
//               7: RANK DEFECT IN HFTI 
//     coded            Dieter Kraft, april 1987 
//     revised                        march 1989 
    // Parameter adjustments 
    --y;
    --x;
    --xu;
    --xl;
    --g;
    --l;
    --b;
    a_dim1 = *la;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --jw;

    // Function Body 
    n1 = *n + 1;
    mineq = *m - *meq;
    m1 = mineq + *n + *n;
//  determine whether to solve problem 
//  with inconsistent linerarization (n2=1) 
//  or not (n2=0) 
    n2 = n1 * *n / 2 + 1;
    if (n2 == *nl) {
	n2 = 0;
    } else {
	n2 = 1;
    }
    n3 = *n - n2;
//  RECOVER MATRIX E AND VECTOR F FROM L AND G 
    i2 = 1;
    i3 = 1;
    i4 = 1;
    ie = 1;
    if__ = *n * *n + 1;
    i__1 = n3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = n1 - i__;
	diag = std::sqrt(l[i2]);
	w[i3] = (Real) 0.0;
	CBLAS::copy(i1, &w[i3], 0, &w[i3], 1);
	i__2 = i1 - n2;
	CBLAS::copy(i__2, &l[i2], 1, &w[i3], *n);
	i__2 = i1 - n2;
	CBLAS::scal(i__2, diag, &w[i3], *n);
	w[i3] = diag;
	i__2 = i__ - 1;
	w[if__ - 1 + i__] = (g[i__] - CBLAS::dot(i__2, &w[i4], 1, &w[if__]
		, 1)) / diag;
	i2 = i2 + i1 - n2;
	i3 += n1;
	i4 += *n;
// L10: 
    }
    if (n2 == 1) {
	w[i3] = l[*nl];
	w[i4] = (Real) 0.0;
	CBLAS::copy(n3, &w[i4], 0, &w[i4], 1);
	w[if__ - 1 + *n] = 0.0;
    }
    d__1 = -one;
    CBLAS::scal(*n, d__1, &w[if__], 1);
    ic = if__ + *n;
    id = ic + *meq * *n;
    if (*meq > 0) {
//  RECOVER MATRIX C FROM UPPER PART OF A 
	i__1 = *meq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    CBLAS::copy(*n, &a[i__ + a_dim1], *la, &w[ic - 1 + i__], *meq);
// L20: 
	}
//  RECOVER VECTOR D FROM UPPER PART OF B 
	CBLAS::copy(*meq, &b[1], 1, &w[id], 1);
	d__1 = -one;
	CBLAS::scal(*meq, d__1, &w[id], 1);
    }
    ig = id + *meq;
    if (mineq > 0) {
//  RECOVER MATRIX G FROM LOWER PART OF A 
	i__1 = mineq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    CBLAS::copy(*n, &a[*meq + i__ + a_dim1], *la, &w[ig - 1 + i__], m1);
// L30: 
	}
    }
//  AUGMENT MATRIX G BY +I AND -I 
    ip = ig + mineq;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[ip - 1 + i__] = (Real) 0.0;
	CBLAS::copy(*n, &w[ip - 1 + i__], 0, &w[ip - 1 + i__], m1);
// L40: 
    }
    i__1 = m1 + 1;
    // SGJ, 2010: skip constraints for infinite bounds 
    for (i__ = 1; i__ <= *n; ++i__)
	 if (!std::isinf(xl[i__])) w[(ip - i__1) + i__ * i__1] = +1.0;
    // Old code: w[ip] = one; CBLAS::copy(n, &w[ip], 0, &w[ip], i__1); 
    im = ip + *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[im - 1 + i__] = (Real) 0.0;
	CBLAS::copy(*n, &w[im - 1 + i__], 0, &w[im - 1 + i__], m1);
// L50: 
    }
    i__1 = m1 + 1;
    // SGJ, 2010: skip constraints for infinite bounds 
    for (i__ = 1; i__ <= *n; ++i__)
	 if (!std::isinf(xu[i__])) w[(im - i__1) + i__ * i__1] = (Real) -1.0;
    // Old code: w[im] = -one;  CBLAS::copy(n, &w[im], 0, &w[im], i__1); 
    ih = ig + m1 * *n;
    if (mineq > 0) {
//  RECOVER H FROM LOWER PART OF B 
	CBLAS::copy(mineq, &b[*meq + 1], 1, &w[ih], 1);
	d__1 = -one;
	CBLAS::scal(mineq, d__1, &w[ih], 1);
    }
//  AUGMENT VECTOR H BY XL AND XU 
    il = ih + mineq;
    iu = il + *n;
    // SGJ, 2010: skip constraints for infinite bounds 
    for (i__ = 1; i__ <= *n; ++i__) {
	 w[(il-1) + i__] = std::isinf(xl[i__]) ? 0 : xl[i__];
	 w[(iu-1) + i__] = std::isinf(xu[i__]) ? 0 : -xu[i__];
    }
    // Old code: CBLAS::copy(n, &xl[1], 1, &w[il], 1);
                 CBLAS::copy(*n, &xu[1], 1, &w[iu], 1);
		 d__1 = -one; CBLAS::scal(*n, d__1, &w[iu], 1); 
    iw = iu + *n;
    i__1 = std::max(1,*meq);
    lsei_(&w[ic], &w[id], &w[ie], &w[if__], &w[ig], &w[ih], &i__1, meq, n, n, 
	    &m1, &m1, n, &x[1], &xnorm, &w[iw], &jw[1], mode);
    if (*mode == 1) {
//   restore Lagrange multipliers 
	CBLAS::copy(*m, &w[iw], 1, &y[1], 1);
	CBLAS::copy(n3, &w[iw + *m], 1, &y[*m + 1], 1);
	CBLAS::copy(n3, &w[iw + *m + *n], 1, &y[*m + n3 + 1], 1);

	// SGJ, 2010: make sure bound constraints are satisfied, since
	//   roundoff error sometimes causes slight violations and
	//   NLopt guarantees that bounds are strictly obeyed 
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	     if (x[i__] < xl[i__]) x[i__] = xl[i__];
	     else if (x[i__] > xu[i__]) x[i__] = xu[i__];
	}
    }
//   END OF SUBROUTINE LSQ 
} // lsq_ 

static void ldl_(int *n, Real *a, Real *z__, 
	Real *sigma, Real *w)
{
    // Initialized data 

    const Real one = (Real) 1.;
    const Real four = (Real) 4.;
    const Real epmach = std::numeric_limits<Real>::epsilon();

    // System generated locals 
    int i__1, i__2;

    // Local variables 
    int i__, j;
    Real t, u, v;
    int ij;
    Real tp, beta, gamma_, alpha, delta;

//   LDL     LDL' - RANK-ONE - UPDATE 
//   PURPOSE: 
//           UPDATES THE LDL' FACTORS OF MATRIX A BY RANK-ONE MATRIX 
//           SIGMA*Z*Z' 
//   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION) 
//     N     : ORDER OF THE COEFFICIENT MATRIX A 
//   * A     : POSITIVE DEFINITE MATRIX OF DIMENSION N; 
//             ONLY THE LOWER TRIANGLE IS USED AND IS STORED COLUMN BY 
//             COLUMN AS ONE DIMENSIONAL ARRAY OF DIMENSION N*(N+1)/2. 
//   * Z     : VECTOR OF DIMENSION N OF UPDATING ELEMENTS 
//     SIGMA : SCALAR FACTOR BY WHICH THE MODIFYING DYADE Z*Z' IS 
//             MULTIPLIED 
//   OUTPUT ARGUMENTS: 
//     A     : UPDATED LDL' FACTORS 
//   WORKING ARRAY: 
//     W     : VECTOR OP DIMENSION N (USED ONLY IF SIGMA .LT. ZERO) 
//   METHOD: 
//     THAT OF FLETCHER AND POWELL AS DESCRIBED IN : 
//     FLETCHER,R.,(1974) ON THE MODIFICATION OF LDL' FACTORIZATION. 
//     POWELL,M.J.D.      MATH.COMPUTATION 28, 1067-1078. 
//   IMPLEMENTED BY: 
//     KRAFT,D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME 
//               D-8031  OBERPFAFFENHOFEN 
//   STATUS: 15. JANUARY 1980 
//   SUBROUTINES REQUIRED: NONE 
    // Parameter adjustments 
    --w;
    --z__;
    --a;

    // Function Body 
    if (*sigma == (Real) 0.0) {
	goto L280;
    }
    ij = 1;
    t = one / *sigma;
    if (*sigma > (Real) 0.0) {
	goto L220;
    }
// PREPARE NEGATIVE UPDATE 
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
// L150: 
	w[i__] = z__[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v = w[i__];
	t += v * v / a[ij];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++ij;
// L160: 
	    w[j] -= v * a[ij];
	}
// L170: 
	++ij;
    }
    if (t >= 0.0) {
	t = epmach / *sigma;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + 1 - i__;
	ij -= i__;
	u = w[j];
	w[j] = t;
// L210: 
	t -= u * u / a[ij];
    }
L220:
// HERE UPDATING BEGINS 
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v = z__[i__];
	delta = v / a[ij];
	if (*sigma < (Real) 0.0) {
	    tp = w[i__];
	}
	else // if (*sigma > 0.0), since *sigma != 0 from above  
        {
	    tp = t + delta * v;
	}
	alpha = tp / t;
	a[ij] = alpha * a[ij];
	if (i__ == *n) {
	    goto L280;
	}
	beta = delta / tp;
	if (alpha > four) {
	    goto L240;
	}
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++ij;
	    z__[j] -= v * a[ij];
// L230: 
	    a[ij] += beta * z__[j];
	}
	goto L260;
L240:
	gamma_ = t / tp;
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++ij;
	    u = a[ij];
	    a[ij] = gamma_ * u + beta * z__[j];
// L250: 
	    z__[j] -= v * u;
	}
L260:
	++ij;
// L270: 
	t = tp;
    }
L280:
    return;
// END OF LDL 
} // ldl_ 

struct slsqpb_state {
    Real t, f0, h1, h2, h3, h4;
    int n1, n2, n3;
    Real t0, gs;
    Real tol;
    int line;
    Real alpha;
    int iexact;
    int incons, ireset, itermx;
    Real *x0;

    slsqpb_state()
    {
      t = f0 = h1 = h2 = h3 = h4 = t0 = gs = tol = alpha = (Real) 0.0;
      n1 = n2 = n3 = line = iexact = incons = ireset = itermx = 0;
      x0 = NULL;
    }
};

#define SS(var) state->var = var
#define SAVE_STATE \
     SS(t); SS(f0); SS(h1); SS(h2); SS(h3); SS(h4);	\
     SS(n1); SS(n2); SS(n3); \
     SS(t0); SS(gs); \
     SS(tol); \
     SS(line); \
     SS(alpha); \
     SS(iexact); \
     SS(incons); SS(ireset); SS(itermx)

#define RS(var) var = state->var
#define RESTORE_STATE \
     RS(t); RS(f0); RS(h1); RS(h2); RS(h3); RS(h4);	\
     RS(n1); RS(n2); RS(n3); \
     RS(t0); RS(gs); \
     RS(tol); \
     RS(line); \
     RS(alpha); \
     RS(iexact); \
     RS(incons); RS(ireset); RS(itermx)

static void slsqpb_(int *m, int *meq, int *la, int *
		    n, Real *x, const Real *xl, const Real *xu, Real *f, 
		    Real *c__, Real *g, Real *a, Real *acc, 
		    int *iter, int *mode, Real *r__, Real *l, 
		    Real *x0, Real *mu, Real *s, Real *u, 
		    Real *v, Real *w, int *iw, 
		    slsqpb_state *state)
{
    // Initialized data 

    const Real one = (Real) 1.;
    const Real alfmin = (Real) .1;
    const Real hun = (Real) 100.;
    const Real ten = (Real) 10.;
    const Real two = (Real) 2.;

    // System generated locals 
    int a_dim1, a_offset, i__1, i__2;
    Real d__1, d__2;

    // Local variables 
    int i__, j, k;

    // saved state from one call to the next;
    //   SGJ 2010: save/restore via state parameter, to make re-entrant. 
    Real t, f0, h1, h2, h3, h4;
    int n1, n2, n3;
    Real t0, gs;
    Real tol;
    int line;
    Real alpha;
    int iexact;
    int incons, ireset, itermx;
    RESTORE_STATE;

//   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS 
//        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  - 
//                      BODY SUBROUTINE FOR SLSQP 
//     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ 
//                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ 
//                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI 
//                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1 
    // Parameter adjustments 
    --mu;
    --c__;
    --v;
    --u;
    --s;
    --x0;
    --l;
    --r__;
    a_dim1 = *la;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --g;
    --xu;
    --xl;
    --x;
    --w;
    --iw;

    // Function Body 
    if (*mode == -1) {
	goto L260;
    } else if (*mode == 0) {
	goto L100;
    } else {
	goto L220;
    }
L100:
    itermx = *iter;
    if (*acc >= (Real) 0.0) {
	iexact = 0;
    } else {
	iexact = 1;
    }
    *acc = std::fabs(*acc);
    tol = ten * *acc;
    *iter = 0;
    ireset = 0;
    n1 = *n + 1;
    n2 = n1 * *n / 2;
    n3 = n2 + 1;
    s[1] = (Real) 0.0;
    mu[1] = (Real) 0.0;
    CBLAS::copy(*n, &s[1], 0, &s[1], 1);
    CBLAS::copy(*m, &mu[1], 0, &mu[1], 1);
//   RESET BFGS MATRIX 
L110:
    ++ireset;
    if (ireset > 5) {
	goto L255;
    }
    l[1] = (Real) 0.0;
    CBLAS::copy(n2, &l[1], 0, &l[1], 1);
    j = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[j] = one;
	j = j + n1 - i__;
// L120: 
    }
//   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE 
L130:
    ++(*iter);
    *mode = 9;
    if (*iter > itermx && itermx > 0) { // SGJ 2010: ignore if itermx <= 0 
	goto L330;
    }
//   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM 
    CBLAS::copy(*n, &xl[1], 1, &u[1], 1);
    CBLAS::copy(*n, &xu[1], 1, &v[1], 1);
    d__1 = -one;
    CBLAS::axpy(*n, d__1, &x[1], 1, &u[1], 1);
    d__1 = -one;
    CBLAS::axpy(*n, d__1, &x[1], 1, &v[1], 1);
    h4 = one;
    lsq_(m, meq, n, &n3, la, &l[1], &g[1], &a[a_offset], &c__[1], &u[1], &v[1]
	    , &s[1], &r__[1], &w[1], &iw[1], mode);

//   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION 
    if (*mode == 6) {
	if (*n == *meq) {
	    *mode = 4;
	}
    }
    if (*mode == 4) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    if (j <= *meq) {
		a[j + n1 * a_dim1] = -c__[j];
	    } else {
// Computing MAX 
		d__1 = -c__[j];
		a[j + n1 * a_dim1] = std::max(d__1,0.0);
	    }
// L140: 
	}
	s[1] = (Real) 0.0;
	CBLAS::copy(*n, &s[1], 0, &s[1], 1);
	h3 = (Real) 0.0;
	g[n1] = (Real) 0.0;
	l[n3] = hun;
	s[n1] = one;
	u[n1] = (Real) 0.0;
	v[n1] = one;
	incons = 0;
L150:
	lsq_(m, meq, &n1, &n3, la, &l[1], &g[1], &a[a_offset], &c__[1], &u[1],
		 &v[1], &s[1], &r__[1], &w[1], &iw[1], mode);
	h4 = one - s[n1];
	if (*mode == 4) {
	    l[n3] = ten * l[n3];
	    ++incons;
	    if (incons > 5) {
		goto L330;
	    }
	    goto L150;
	} else if (*mode != 1) {
	    goto L330;
	}
    } else if (*mode != 1) {
	goto L330;
    }
//   UPDATE MULTIPLIERS FOR L1-TEST 
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__] = g[i__] - CBLAS::dot(*m, &a[i__ * a_dim1 + 1], 1, &r__[1], 1);
// L160: 
    }
    f0 = *f;
    CBLAS::copy(*n, &x[1], 1, &x0[1], 1);
    gs = CBLAS::dot(*n, &g[1], 1, &s[1], 1);
    h1 = std::fabs(gs);
    h2 = (Real) 0.0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (j <= *meq) {
	    h3 = c__[j];
	} else {
	    h3 = (Real) 0.0;
	}
// Computing MAX 
	d__1 = -c__[j];
	h2 += std::max(d__1,h3);
	h3 = (d__1 = r__[j], std::fabs(d__1));
// Computing MAX 
	d__1 = h3, d__2 = (mu[j] + h3) / two;
	mu[j] = std::max(d__1,d__2);
	h1 += h3 * (d__1 = c__[j], std::fabs(d__1));
// L170: 
    }
//   CHECK CONVERGENCE 
    *mode = 0;
    if (h1 < *acc && h2 < *acc) {
	goto L330;
    }
    h1 = (Real) 0.0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (j <= *meq) {
	    h3 = c__[j];
	} else {
	    h3 = (Real) 0.0;
	}
// Computing MAX 
	d__1 = -c__[j];
	h1 += mu[j] * std::max(d__1,h3);
// L180: 
    }
    t0 = *f + h1;
    h3 = gs - h1 * h4;
    *mode = 8;
    if (h3 >= (Real) 0.0) {
	goto L110;
    }
//   LINE SEARCH WITH AN L1-TESTFUNCTION 
    line = 0;
    alpha = one;
    if (iexact == 1) {
	goto L210;
    }
//   INEXACT LINESEARCH 
L190:
    ++line;
    h3 = alpha * h3;
    CBLAS::scal(*n, alpha, &s[1], 1);
    CBLAS::copy(*n, &x0[1], 1, &x[1], 1);
    CBLAS::axpy(*n, one, &s[1], 1, &x[1], 1);
    
    // SGJ 2010: ensure roundoff doesn't push us past bound constraints 
    i__1 = *n; for (i__ = 1; i__ <= i__1; ++i__) {
	 if (x[i__] < xl[i__]) x[i__] = xl[i__];
	 else if (x[i__] > xu[i__]) x[i__] = xu[i__];
    }

    // SGJ 2010: optimizing for the common case where the inexact line
    //   search succeeds in one step, use special mode = -2 here to
    //   eliminate a a subsequent unnecessary mode = -1 call, at the 
    //   expense of extra gradient evaluations when more than one inexact
    //   line-search step is required 
    *mode = line == 1 ? -2 : 1;
    goto L330;
L200:
    if (h1 <= h3 / ten || line > 10) {
	goto L240;
    }
// Computing MAX 
    d__1 = h3 / (two * (h3 - h1));
    alpha = std::max(d__1,alfmin);
    goto L190;
//   EXACT LINESEARCH 
L210:
    *mode = 9; // will yield nlopt_failure 
    return;
    CBLAS::scal(*n, alpha, &s[1], 1);
    goto L240;
//   CALL FUNCTIONS AT CURRENT X 
L220:
    t = *f;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (j <= *meq) {
	    h1 = c__[j];
	} else {
	    h1 = (Real) 0.0;
	}
// Computing MAX 
	d__1 = -c__[j];
	t += mu[j] * std::max(d__1,h1);
// L230: 
    }
    h1 = t - t0;
    switch (iexact + 1) {
	case 1:  goto L200;
	case 2:  goto L210;
    }
//   CHECK CONVERGENCE 
L240:
    h3 = 0.0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (j <= *meq) {
	    h1 = c__[j];
	} else {
	    h1 = (Real) 0.0;
	}
// Computing MAX 
	d__1 = -c__[j];
	h3 += std::max(d__1,h1);
// L250: 
    }
    if (((d__1 = *f - f0, std::fabs(d__1)) < *acc || CBLAS::nrm2(*n, &s[1], 1) < *
	    acc) && h3 < *acc) {
	*mode = 0;
    } else {
	*mode = -1;
    }
    goto L330;
//   CHECK relaxed CONVERGENCE in case of positive directional derivative 
L255:
    if (((d__1 = *f - f0, std::fabs(d__1)) < tol || CBLAS::nrm2(*n, &s[1], 1) < tol)
	     && h3 < tol) {
	*mode = 0;
    } else {
	*mode = 8;
    }
    goto L330;
//   CALL JACOBIAN AT CURRENT X 
//   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA 
L260:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u[i__] = g[i__] - CBLAS::dot(*m, &a[i__ * a_dim1 + 1], 1, &r__[1], 1) - v[i__];
// L270: 
    }
//   L'*S 
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h1 = (Real) 0.0;
	++k;
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++k;
	    h1 += l[k] * s[j];
// L280: 
	}
	v[i__] = s[i__] + h1;
// L290: 
    }
//   D*L'*S 
    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__] = l[k] * v[i__];
	k = k + n1 - i__;
// L300: 
    }
//   L*D*L'*S 
    for (i__ = *n; i__ >= 1; --i__) {
	h1 = (Real) 0.0;
	k = i__;
	i__1 = i__ - 1;
	for (j = 1; j <= i__1; ++j) {
	    h1 += l[k] * v[j];
	    k = k + *n - j;
// L310: 
	}
	v[i__] += h1;
// L320: 
    }
    h1 = CBLAS::dot(*n, &s[1], 1, &u[1], 1);
    h2 = CBLAS::dot(*n, &s[1], 1, &v[1], 1);
    h3 = h2 * (Real) .2;
    if (h1 < h3) {
	h4 = (h2 - h3) / (h2 - h1);
	h1 = h3;
	CBLAS::scal(*n, h4, &u[1], 1);
	d__1 = one - h4;
	CBLAS::axpy(*n, d__1, &v[1], 1, &u[1], 1);
    }
    d__1 = one / h1;
    ldl_(n, &l[1], &u[1], &d__1, &v[1]);
    d__1 = -one / h2;
    ldl_(n, &l[1], &v[1], &d__1, &u[1]);
//   END OF MAIN ITERATION 
    goto L130;
//   END OF SLSQPB 
L330:
    SAVE_STATE;
} // slsqpb_ 

// *********************************************************************** 
//                              optimizer                               * 
// *********************************************************************** 
static void slsqp(int *m, int *meq, int *la, int *n,
		  Real *x, const Real *xl, const Real *xu, Real *f, 
		  Real *c__, Real *g, Real *a, Real *acc, 
		  int *iter, int *mode, Real *w, int *l_w__, int *
		  jw, int *l_jw__, slsqpb_state *state)
{
    // System generated locals 
    int a_dim1, a_offset, i__1, i__2;

    // Local variables 
    int n1, il, im, ir, is, iu, iv, iw, ix, mineq;

//   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING 
//            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS 
// *********************************************************************** 
// *                                                                     * 
// *                                                                     * 
// *            A NONLINEAR PROGRAMMING METHOD WITH                      * 
// *            QUADRATIC  PROGRAMMING  SUBPROBLEMS                      * 
// *                                                                     * 
// *                                                                     * 
// *  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   * 
// *                                                                     * 
// *            MINIMIZE    F(X)                                         * 
// *                                                                     * 
// *            SUBJECT TO  C (X) .EQ. 0  ,  J = 1,...,MEQ               * 
// *                         J                                           * 
// *                                                                     * 
// *                        C (X) .GE. 0  ,  J = MEQ+1,...,M             * 
// *                         J                                           * 
// *                                                                     * 
// *                        XL .LE. X .LE. XU , I = 1,...,N.             * 
// *                          I      I       I                           * 
// *                                                                     * 
// *  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              * 
// *  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              * 
// *  WITHIN THE STEPLENGTH ALGORITHM.                                   * 
// *                                                                     * 
// *    PARAMETER DESCRIPTION:                                           * 
// *    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    * 
// *                                                                     * 
// *    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .GE. 0      * 
// *    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .GE. 0 * 
// *    LA             SEE A, LA .GE. MAX(M,1)                           * 
// *    N              IS THE NUMBER OF VARIBLES, N .GE. 1               * 
// *  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  * 
// *                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     * 
// *                   STORES THE SOLUTION VECTOR X IF MODE = 0.         * 
// *    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  * 
// *    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  * 
// *    F              IS THE VALUE OF THE OBJECTIVE FUNCTION.           * 
// *    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         * 
// *                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              * 
// *                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       * 
// *                   which must be GREATER OR EQUAL MAX(1,M).          * 
// *    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      * 
// *                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        * 
// *                   GREATER OR EQUAL N+1.                             * 
// *    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  * 
// *                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        * 
// *                   A() HAS FIRST DIMENSIONING PARAMETER LA,          * 
// *                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          * 
// *    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     * 
// *  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             * 
// *                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,* 
// *                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      * 
// *  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      * 
// *                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  * 
// *  * MODE           MODE CONTROLS CALCULATION:                        * 
// *                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   * 
// *                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS* 
// *                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN* 
// *                   WITH MODE .NE. IABS(1) TAKES PLACE.               * 
// *                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     * 
// *                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATED 
// *                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS * 
// *                   OF SQP.                                           * 
// *                   EVALUATION MODES:                                 * 
// *        MODE = -2,-1: GRADIENT EVALUATION, (G&A)                     * 
// *                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               * 
// *                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED * 
// *                1: FUNCTION EVALUATION, (F&C)                        * 
// *                                                                     * 
// *                   FAILURE MODES:                                    * 
// *                2: NUMBER OF EQUALITY CONTRAINTS LARGER THAN N       * 
// *                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        * 
// *                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               * 
// *                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               * 
// *                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               * 
// *                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI* 
// *                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    * 
// *                9: MORE THAN ITER ITERATIONS IN SQP                  * 
// *             >=10: WORKING SPACE W OR JW TOO SMALL,                  * 
// *                   W SHOULD BE ENLARGED TO L_W=MODE/1000             * 
// *                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       * 
// *  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           * 
// *                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        * 
// *                   (3*N1+M)*(N1+1)                        for LSQ    * 
// *                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    * 
// *                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   * 
// *                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB * 
// *                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          * 
// *        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO * 
// *                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    * 
// *                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    * 
// ####################################################################### 
//     INT LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ 
//     PARAMETER (M=... , MEQ=... , N=...  ) 
//     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1) 
//     PARAMETER (LEN_W= 
//    $           (3*N1+M)*(N1+1) 
//    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ 
//    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1 
//    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1, 
//    $           LEN_JW=MINEQ) 
//     DOUBLE PRECISION W(LEN_W) 
//     INT          JW(LEN_JW) 
// ####################################################################### 
// *                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    * 
// *                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        * 
// *                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   * 
// *                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    * 
// *                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR* 
// *                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        * 
// *                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   * 
// *                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        * 
// *                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       * 
// *                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       * 
// *                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  * 
// *                   THE SEARCH DIRECTION TO THE SOLUTION X*           * 
// *  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INT WORKING SPACE   * 
// *                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       * 
// *                   MINEQ                                             * 
// *                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          * 
// *                                                                     * 
// *  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 * 
// *     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          * 
// *     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          * 
// *     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              * 
// *                                                                     * 
// *        SOLUTION OF THE QUADRATIC PROGRAM                            * 
// *                QPSOL IS RECOMMENDED:                                * 
// *     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      * 
// *     USER'S GUIDE FOR SOL/QPSOL:                                     * 
// *     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    * 
// *     TECHNICAL REPORT SOL 83-7, JULY 1983                            * 
// *     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          * 
// *     STANFORD, CA 94305                                              * 
// *     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                * 
// *     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               * 
// *                                                                     * 
// *     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      * 
// *     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   * 
// *     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   * 
// *     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          * 
// *     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. * 
// *                                                                     * 
// *     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         * 
// *                                                                     * 
// *     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               * 
// *     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    * 
// *                                                                     * 
// *  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               * 
// *  as described in Dieter Kraft: A Software Package for               * 
// *                                Sequential Quadratic Programming     * 
// *                                DFVLR-FB 88-28, 1988                 * 
// *  which should be referenced if the user publishes results of SLSQP  * 
// *                                                                     * 
// *  DATE:           APRIL - OCTOBER, 1981.                             * 
// *  STATUS:         DECEMBER, 31-ST, 1984.                             * 
// *  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTAN 77        * 
// *  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       * 
// *  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       * 
// *  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      * 
// *                                         accepts Statement Functions * 
// *  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         * 
// *                                         FTN77/386 COMPILER VERS 2.40* 
// *                                         in protected mode           * 
// *                                                                     * 
// *********************************************************************** 
// *                                                                     * 
// *  Copyright 1991: Dieter Kraft, FHM                                  * 
// *                                                                     * 
// *********************************************************************** 
//     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ 
//                    +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ          for LSI 
//                    +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1        for LSEI 
//                    + N1*N/2 + 2*M + 3*N +3*N1 + 1           for SLSQPB 
//                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1 
//   CHECK LENGTH OF WORKING ARRAYS 
    // Parameter adjustments 
    --c__;
    a_dim1 = *la;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --g;
    --xu;
    --xl;
    --x;
    --w;
    --jw;

    // Function Body 
    n1 = *n + 1;
    mineq = *m - *meq + n1 + n1;
    il = (n1 * 3 + *m) * (n1 + 1) + (n1 - *meq + 1) * (mineq + 2) + (mineq << 
	    1) + (n1 + mineq) * (n1 - *meq) + (*meq << 1) + n1 * *n / 2 + (*m 
	    << 1) + *n * 3 + (n1 << 2) + 1;
// Computing MAX 
    i__1 = mineq, i__2 = n1 - *meq;
    im = std::max(i__1,i__2);
    if (*l_w__ < il || *l_jw__ < im) {
	*mode = std::max(10,il) * 1000;
	*mode += std::max(10,im);
	return;
    }
//   PREPARE DATA FOR CALLING SQPBDY  -  INITIAL ADDRESSES IN W 
    im = 1;
    il = im + std::max(1,*m);
    il = im + *la;
    ix = il + n1 * *n / 2 + 1;
    ir = ix + *n;
    is = ir + *n + *n + std::max(1,*m);
    is = ir + *n + *n + *la;
    iu = is + n1;
    iv = iu + n1;
    iw = iv + n1;
    slsqpb_(m, meq, la, n, &x[1], &xl[1], &xu[1], f, &c__[1], &g[1], &a[
	    a_offset], acc, iter, mode, &w[ir], &w[il], &w[ix], &w[im], &w[is]
	    , &w[iu], &w[iv], &w[iw], &jw[1], state);
    state->x0 = &w[ix];
    return;
} // slsqp_ 

static void length_work(int *LEN_W, int *LEN_JW, int M, int MEQ, int N)
{
     int N1 = N+1, MINEQ = M-MEQ+N1+N1;
     *LEN_W = (3*N1+M)*(N1+1) 
	  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1
          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1;
     *LEN_JW = MINEQ;
}

void SLSQP(unsigned M, unsigned MEQ, unsigned LA, unsigned N, Real* X, const Real* XL, const Real* XU, Real F, Real* C, Real* G, Real* A, Real ACC, unsigned& ITER, int& MODE, Real* W, const unsigned L_W, int* JW, const unsigned L_JW, shared_ptr<void>& state)
{
  shared_ptr<slsqpb_state> sstate;

  #ifdef THREADED 
  std::cerr << "SLSQP() was not written to be thread-safe!" << std::endl;
  #endif

  // setup the state of the optimizer
  if (!state)
  {
    sstate = shared_ptr<slsqpb_state>(new slsqpb_state);
    state = sstate;
  }
  else
    sstate = boost::static_pointer_cast<slsqpb_state>(state);

  slsqp((int*) &M, (int*) &MEQ, (int*) &LA, (int*) &N, X, XL, XU, &F, C, G, A, &ACC, (int*) &ITER, &MODE, W, (int*) &L_W, JW, (int*) &L_JW, sstate.get());
}


