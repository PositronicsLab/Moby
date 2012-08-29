#include <stdlib.h>
#include <cmath>
#include "f2c-ap.h"
#include <cstring>
#include <malloc.h>

extern "C" {

doublereal d_lg10(doublereal *x)
{
  const doublereal log10e = "0.43429448190325182765";
  return log10e * std::log(*x);
}

doublereal pow_dd(doublereal* ap, doublereal* bp)
{
  return std::pow(*ap, *bp);
}

doublereal d_sign(doublereal *a, doublereal *b)
{
  doublereal x;
  x = (*a >= (doublereal) 0 ? *a : - *a);
  return( *b >= (doublereal) 0 ? x : -x);
}

doublereal pow_di(doublereal *ap, integer *bp)
{
  doublereal pw, x;
  integer n;
  unsigned long u;

  pw = (doublereal) 1;
  x = *ap;
  n = *bp;

  if(n != 0)
  {
    if(n < 0)
    {
      n = -n;
      x = (doublereal) 1/x;
    }
    for(u = n; ; )
    {
      if(u & 01)
        pw *= x;
      if(u >>= 1)
        x *= x;
      else
        break;
    }
  }
  return(pw);
}

integer pow_ii(integer *ap, integer *bp)
{
  integer pw, x, n;
  unsigned long u;

  x = *ap;
  n = *bp;

  if (n <= 0) 
  {
    if (n == 0 || x == 1)
      return 1;
    if (x != -1)
      return x == 0 ? 1/x : 0;
    n = -n;
  }
  u = n;
  for (pw = 1; ; )
  {
    if(u & 01)
      pw *= x;
    if(u >>= 1)
      x *= x;
    else
      break;
  }
  return(pw);
}

void s_cat(char *lp, char *rpp[], ftnint rnp[], ftnint *np, ftnlen ll)
{
  ftnlen i, nc;
  char *rp;
  ftnlen n = *np;
#ifndef NO_OVERWRITE
  ftnlen L, m;
  char *lp0, *lp1;

  lp0 = 0;
  lp1 = lp;
  L = ll;
  i = 0;
  while(i < n) 
  {
    rp = rpp[i];
    m = rnp[i++];
    if (rp >= lp1 || rp + m <= lp) 
    {
      if ((L -= m) <= 0) 
      {
        n = i;
        break;
      }
      lp1 += m;
      continue;
    }
    lp0 = lp;
    lp = lp1 = (char*) malloc(L = ll);
    break;
  }
  lp1 = lp;
#endif /* NO_OVERWRITE */
  for(i = 0 ; i < n ; ++i) 
  {
    nc = ll;
    if(rnp[i] < nc)
      nc = rnp[i];
    ll -= nc;
    rp = rpp[i];
    while(--nc >= 0)
      *lp++ = *rp++;
  }
  while(--ll >= 0)
    *lp++ = ' ';
#ifndef NO_OVERWRITE
  if (lp0) 
  {
    memcpy(lp0, lp1, L);
    free(lp1);
  }
#endif
}

doublereal r_sign(real *a, real *b)
{
  doublereal x;
  x = (*a >= (doublereal) 0 ? *a : - *a);
  return( *b >= (doublereal) 0 ? x : -x);
}

doublereal r_imag(complex* z)
{
  return z->i;
}

integer i_nint(real* x)
{
  return (integer) (*x >= (real) 0 ? std::floor(*x + (real) 0.5) : -std::floor((real) 0.5 - *x));
}

integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
  register unsigned char *a, *aend, *b, *bend;
  a = (unsigned char *)a0;
  b = (unsigned char *)b0;
  aend = a + la;
  bend = b + lb;

  if(la <= lb)
  {
    while(a < aend)
      if(*a != *b)
        return( *a - *b );
      else
      { 
        ++a; 
        ++b; 
      }

    while(b < bend)
      if(*b != ' ')
        return( ' ' - *b );
      else ++b;
  }
  else
  {
    while(b < bend)
      if(*a == *b)
      { 
        ++a; 
        ++b; 
      }
      else
        return( *a - *b );

    while(a < aend)
      if(*a != ' ')
        return(*a - ' ');
      else  ++a;
  }

  return(0);
}

int s_stop(char *s, ftnlen n)
{
  int i;

  if(n > 0)
  {
    fprintf(stderr, "STOP ");
    for(i = 0; i<n ; ++i)
      putc(*s++, stderr);
    fprintf(stderr, " statement executed\n");
  }
  exit(0);

/* We cannot avoid (useless) compiler diagnostics here:    */
/* some compilers complain if there is no return statement,  */
/* and others complain that this one cannot be reached.    */

  return 0; /* NOT REACHED */
}

void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
  register char *aend, *bend;

  aend = a + la;

  if(la <= lb)
#ifndef NO_OVERWRITE
    if (a <= b || a >= b + la)
#endif
      while(a < aend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else
      for(b += la; a < aend; )
        *--aend = *--b;
#endif
  else 
  {
    bend = b + lb;
#ifndef NO_OVERWRITE
    if (a <= b || a >= bend)
#endif
      while(b < bend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else 
    {
      a += lb;
      while(b < bend)
        *--a = *--bend;
      a += lb;
    }
#endif
    while(a < aend)
      *a++ = ' ';
  }
}

} // extern "C"

/*
real f__cabs(real r, real imag)
{
  if (r < (real) 0)
    r = -r;
  if (imag < (real) 0)
    imag = -imag;
  if (imag > r)
    std::swap(r, imag);
  if ((r+imag) == r)
    return r;
  return r*std::sqrt(1.0 + imag*imag/(r*r));
}

real c_abs(complex* z)
{
  return f__cabs(z->r, z->i);
}   

void c_sin(complex* r, complex* z)
{
  real zi = z->i, zr = z->r;
  r->r = std::sin(zr) * std::cosh(zi);
  r->i = std::cos(zr) * std::sinh(zi);
}

void c_cos(complex* r, complex* z)
{
  real zi = z->i, zr = z->r;
  r->r = std::cos(zr) * std::cosh(zi);
  r->i = -std::sin(zr) * std::sinh(zi);
}

void c_sqrt(complex* r, complex* z)
{
  real mag, t;
  real zi = z->i, zr = z->r;

  if( (mag = f__cabs(zr, zi)) == (real) 0.)
    r->r = r->i = 0.;
  else if(zr > 0)
    {
    r->r = t = std::sqrt((real) 0.5 * (mag + zr) );
    t = zi / t;
    r->i = (real) 0.5 * t;
    }
  else
    {
    t = sqrt((real) 0.5 * (mag - zr) );
    if(zi < (real) 0)
      t = -t;
    r->i = t;
    t = zi / t;
    r->r = (real) 0.5 * t;
    }
  }

void c_div(complex* c, complex* a, complex* b)
{
  real ratio, den;
  real abr, abi, cr;
  if( (abr = b->r) < (real) 0.)
    abr = - abr;
  if( (abi = b->i) < (real) 0.)
    abi = - abi;
  if( abr <= abi )
  {
    if(abi == (real) 0) 
    {
#ifdef IEEE_COMPLEX_DIVIDE
      real af, bf;
      af = bf = abr;
      if (a->i != (real) 0 || a->r != (real) 0)
        af = 1.;
      c->i = c->r = af / bf;
      return;
#else
      throw "complex division by zero";
#endif
      }
    ratio = (real) b->r / b->i ;
    den = b->i * (1 + ratio*ratio);
    cr = (a->r*ratio + a->i) / den;
    c->i = (a->i*ratio - a->r) / den;
  }
  else
  {
    ratio = (real) b->i / b->r ;
    den = b->r * (1 + ratio*ratio);
    cr = (a->r + a->i*ratio) / den;
    c->i = (a->i - a->r*ratio) / den;
  }
  c->r = cr;
}

void c_exp(complex* r, complex* z)
{
  real expx, zi = z->i;
  expx = std::exp(z->r);
  r->r = expx * std::cos(zi);
  r->i = expx * std::sin(zi);
}

void c_log(complex* r, complex* z)
{
  real zi, zr;
  r->i = std::atan2(zi = z->i, zr = z->r);
  r->r = std::log(f__cabs(zr, zi));
}

doublereal d_abs(doublereal *x)
{
if(*x >= 0)
  return(*x);
return(- *x);
}

doublereal d_acos(doublereal *x)
{
return( acos(*x) );
}

doublereal d_asin(doublereal *x)
{
return( asin(*x) );
}

doublereal d_atan(doublereal *x)
{
return( atan(*x) );
}

doublereal d_atn2(doublereal *x, doublereal *y)
{
return( atan2(*x,*y) );
}

void d_cnjg(doublecomplex *r, doublecomplex *z)
{
  doublereal zi = z->i;
  r->r = z->r;
  r->i = -zi;
}

doublereal d_cos(doublereal *x)
{
return( cos(*x) );
}
doublereal d_cosh(doublereal *x)
{
return( cosh(*x) );
}
doublereal d_dim(doublereal *a, doublereal *b)
{
return( *a > *b ? *a - *b : (doublereal) 0);
}
doublereal derf_(doublereal *x)
{
return( erf(*x) );
}

doublereal derfc_(doublereal *x)
{
return( erfc(*x) );
}

doublereal d_exp(doublereal *x)
{
return( exp(*x) );
}

 int
y_rsk(void)
{
  if(f__curunit->uend || f__curunit->url <= f__recpos
    || f__curunit->url == 1) return 0;
  do {
    getc(f__cf);
  } while(++f__recpos < f__curunit->url);
  return 0;
}

 int
y_getc(void)
{
  int ch;
  if(f__curunit->uend) return(-1);
  if((ch=getc(f__cf))!=EOF)
  {
    f__recpos++;
    if(f__curunit->url>=f__recpos ||
      f__curunit->url==1)
      return(ch);
    else  return(' ');
  }
  if(feof(f__cf))
  {
    f__curunit->uend=1;
    errno=0;
    return(-1);
  }
  err(f__elist->cierr,errno,"readingd");
}

 static int
y_rev(void)
{
  if (f__recpos < f__hiwater)
    f__recpos = f__hiwater;
  if (f__curunit->url > 1)
    while(f__recpos < f__curunit->url)
      (*f__putn)(' ');
  if (f__recpos)
    f__putbuf(0);
  f__recpos = 0;
  return(0);
}

 static int
y_err(void)
{
  err(f__elist->cierr, 110, "dfe");
}

 static int
y_newrec(void)
{
  y_rev();
  f__hiwater = f__cursor = 0;
  return(1);
}

c_dfe(cilist *a)
{
  f__sequential=0;
  f__formatted=f__external=1;
  f__elist=a;
  f__cursor=f__scale=f__recpos=0;
  f__curunit = &f__units[a->ciunit];
  if(a->ciunit>MXUNIT || a->ciunit<0)
    err(a->cierr,101,"startchk");
  if(f__curunit->ufd==NULL && fk_open(DIR,FMT,a->ciunit))
    err(a->cierr,104,"dfe");
  f__cf=f__curunit->ufd;
  if(!f__curunit->ufmt) err(a->cierr,102,"dfe")
  if(!f__curunit->useek) err(a->cierr,104,"dfe")
  f__fmtbuf=a->cifmt;
  if(a->cirec <= 0)
    err(a->cierr,130,"dfe")
  FSEEK(f__cf,(OFF_T)f__curunit->url * (a->cirec-1),SEEK_SET);
  f__curunit->uend = 0;
  return(0);
}
integer s_rdfe(cilist *a)
{
  int n;
  if(!f__init) f_init();
  f__reading=1;
  if(n=c_dfe(a))return(n);
  if(f__curunit->uwrt && f__nowreading(f__curunit))
    err(a->cierr,errno,"read start");
  f__getn = y_getc;
  f__doed = rd_ed;
  f__doned = rd_ned;
  f__dorevert = f__donewrec = y_err;
  f__doend = y_rsk;
  if(pars_f(f__fmtbuf)<0)
    err(a->cierr,100,"read start");
  fmt_bg();
  return(0);
}

integer s_wdfe(cilist *a)
{
  int n;
  if(!f__init) f_init();
  f__reading=0;
  if(n=c_dfe(a)) return(n);
  if(f__curunit->uwrt != 1 && f__nowwriting(f__curunit))
    err(a->cierr,errno,"startwrt");
  f__putn = x_putc;
  f__doed = w_ed;
  f__doned= w_ned;
  f__dorevert = y_err;
  f__donewrec = y_newrec;
  f__doend = y_rev;
  if(pars_f(f__fmtbuf)<0)
    err(a->cierr,100,"startwrt");
  fmt_bg();
  return(0);
}
integer e_rdfe(void)
{
  en_fio();
  return 0;
}
integer e_wdfe(void)
{
  return en_fio();
}

doublereal d_imag(doublecomplex *z)
{
return(z->i);
}

doublereal d_int(doublereal *x)
{
return( (*x>0) ? floor(*x) : -floor(- *x) );
}

#define log10e 0.43429448190325182765

doublereal d_lg10(doublereal *x)
{
return( log10e * log(*x) );
}

doublereal d_log(doublereal *x)
{
return( log(*x) );
}

doublereal d_mod(doublereal *x, doublereal *y)
{
#ifdef IEEE_drem
  doublereal xa, ya, z;
  if ((ya = *y) < 0.)
    ya = -ya;
  z = drem(xa = *x, ya);
  if (xa > 0) {
    if (z < 0)
      z += ya;
    }
  else if (z > 0)
    z -= ya;
  return z;
#else
  doublereal quotient;
  if( (quotient = *x / *y) >= 0)
    quotient = floor(quotient);
  else
    quotient = -floor(-quotient);
  return(*x - (*y) * quotient );
#endif
}

doublereal d_nint(doublereal *x)
{
return( (*x)>=0 ?
  floor(*x + .5) : -floor(.5 - *x) );
}

integer do_lio(ftnint *type, ftnint *number, char *ptr, ftnlen len)
{
  return((*f__lioproc)(number,ptr,len,*type));
}

doublereal d_prod(real *x, real *y)
{
return( (*x) * (*y) );
}
*/

/*
doublereal d_sin(doublereal *x)
{
return( sin(*x) );
}

doublereal d_sinh(doublereal *x)
{
return( sinh(*x) );
}

doublereal d_sqrt(doublereal *x)
{
return( sqrt(*x) );
}

doublereal d_tan(doublereal *x)
{
return( tan(*x) );
}

doublereal d_tanh(doublereal *x)
{
return( tanh(*x) );
}

#ifndef REAL
#define REAL double
#endif

#ifndef USE_CLOCK
#define _INCLUDE_POSIX_SOURCE  // for HP-UX 
#define _INCLUDE_XOPEN_SOURCE  // for HP-UX 
#include "sys/types.h"
#include "sys/times.h"
#ifdef __cplusplus
extern "C" {
#endif
#endif

#undef Hz
#ifdef CLK_TCK
#define Hz CLK_TCK
#else
#ifdef HZ
#define Hz HZ
#else
#define Hz 60
#endif
#endif

 REAL
#ifdef KR_headers
dtime_(tarray) real *tarray;
#else
dtime_(real *tarray)
#endif
{
#ifdef USE_CLOCK
#ifndef CLOCKS_PER_SECOND
#define CLOCKS_PER_SECOND Hz
#endif
  static doublereal t0;
  doublereal t = clock();
  tarray[1] = 0;
  tarray[0] = (t - t0) / CLOCKS_PER_SECOND;
  t0 = t;
  return tarray[0];
#else
  struct tms t;
  static struct tms t0;

  times(&t);
  tarray[0] = (double)(t.tms_utime - t0.tms_utime) / Hz;
  tarray[1] = (double)(t.tms_stime - t0.tms_stime) / Hz;
  t0 = t;
  return tarray[0] + tarray[1];
#endif
  }

c_due(cilist *a)
{
  if(!f__init) f_init();
  f__sequential=f__formatted=f__recpos=0;
  f__external=1;
  f__curunit = &f__units[a->ciunit];
  if(a->ciunit>=MXUNIT || a->ciunit<0)
    err(a->cierr,101,"startio");
  f__elist=a;
  if(f__curunit->ufd==NULL && fk_open(DIR,UNF,a->ciunit) ) err(a->cierr,104,"due");
  f__cf=f__curunit->ufd;
  if(f__curunit->ufmt) err(a->cierr,102,"cdue")
  if(!f__curunit->useek) err(a->cierr,104,"cdue")
  if(f__curunit->ufd==NULL) err(a->cierr,114,"cdue")
  if(a->cirec <= 0)
    err(a->cierr,130,"due")
  FSEEK(f__cf,(OFF_T)(a->cirec-1)*f__curunit->url,SEEK_SET);
  f__curunit->uend = 0;
  return(0);
}

integer s_rdue(cilist *a)
{
  int n;
  f__reading=1;
  if(n=c_due(a)) return(n);
  if(f__curunit->uwrt && f__nowreading(f__curunit))
    err(a->cierr,errno,"read start");
  return(0);
}

integer s_wdue(cilist *a)
{
  int n;
  f__reading=0;
  if(n=c_due(a)) return(n);
  if(f__curunit->uwrt != 1 && f__nowwriting(f__curunit))
    err(a->cierr,errno,"write start");
  return(0);
}
integer e_rdue(void)
{
  if(f__curunit->url==1 || f__recpos==f__curunit->url)
    return(0);
  FSEEK(f__cf,(OFF_T)(f__curunit->url-f__recpos),SEEK_CUR);
  if(FTELL(f__cf)%f__curunit->url)
    err(f__elist->cierr,200,"syserr");
  return(0);
}
integer e_wdue(void)
{
#ifdef ALWAYS_FLUSH
  if (fflush(f__cf))
    err(f__elist->cierr,errno,"write end");
#endif
  return(e_rdue());
}
#ifdef __cplusplus
}
#endif
*/

