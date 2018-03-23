/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *   Enrico Fucile  - 19.06.2007                                           *
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"


static void init_ibm_table();

typedef struct ibm_table_t ibm_table_t;

struct ibm_table_t {
  int inited;
  double e[128];
  double v[128];
  double vmin;
  double vmax;
};

static ibm_table_t ibm_table={ 0,{0,},{0,} };

static void init_ibm_table() {
  if (!ibm_table.inited) {
    unsigned long i;
    unsigned long mmin=0x100000;
    unsigned long mmax=0xffffff;
    double e=1;
    for (i=1; i<=57;i++) { e*=16; ibm_table.e[i+70]=e;ibm_table.v[i+70]=e*mmin;}
    ibm_table.e[70]=1;
    ibm_table.v[70]=mmin;
    e=1;
    for (i=1; i<=70;i++) { e/=16; ibm_table.e[-i+70]=e;ibm_table.v[-i+70]=e*mmin;}
    ibm_table.vmin=ibm_table.v[0];
    ibm_table.vmax=ibm_table.e[127]*mmax;
    ibm_table.inited=1;
    /*for (i=0;i<128;i++) printf("++++ ibm_table.v[%d]=%g\n",i,ibm_table.v[i]);*/
  }
}

static void binary_search(double xx[], const unsigned long n, double x, unsigned long *j)
{
  /*These routine works only on ascending ordered arrays*/
  unsigned long ju,jm,jl;
  jl=0;
  ju=n;
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    /* printf("jl=%lu jm=%lu ju=%lu\n",jl,jm,ju); */
    /* printf("xx[jl]=%.10e xx[jm]=%.10e xx[ju]=%.10e\n",xx[jl],xx[jm],xx[ju]); */
    if (x >= xx[jm]) jl=jm;
    else ju=jm;
  }
  *j=jl;
}

unsigned long grib_ibm_to_long(double x) {
  unsigned long s = 0;
  unsigned long mmax = 0xffffff;
  unsigned long mmin = 0x800000;
  unsigned long m = mmax;
  unsigned long e=0;
  double rmmax=mmax+0.5;

  if (!ibm_table.inited) init_ibm_table();

  /* printf("\ngrib_ibm_to_long: x=%.20e\n",x); */
  if (x < 0)  {  s  = 1; x = -x; }

  /* Underflow */
  if (x < ibm_table.vmin) {
    /*printf("grib_ibm_to_long: (x < ibm_table.vmin) x=%.20e vmin=%.20e v=0x%lX\n",x,ibm_table.vmin,(s<<31));*/
    return (s << 31);
  }

  /* Overflow */
  Assert(x <= ibm_table.vmax);
/*
  if (x > ibm_table.vmax) {
    printf("grib_ibm_to_long: x=%.20e > xmax=%.20e\n",x,ibm_table.vmax);
  return 0;
  }
*/

  binary_search(ibm_table.v, 127, x, &e);

  /* printf("grib_ibm_to_long: e=%ld\n",e); */

  x/=ibm_table.e[e];

  /* printf("grib_ibm_to_long: x=%.20e\n",x); */

  while(x < mmin ) { x *= 16; e--;
    /* printf("grib_ibm_to_long (e--): x=%.20e e=%ld \n",x,e); */
  }

  while(x > rmmax ) { x /= 16; e++;
    /* printf("grib_ibm_to_long (e++): x=%.20e e=%ld \n",x,e); */
  }

  m=x+0.5;
  /* printf("grib_ibm_to_long: m=0x%lX (%lu) x=%.10e \n",m,m,x ); */
  if ( m > mmax ) { e++; m=0x800000;
     /* printf("grib_ibm_to_long: ( m > mmax ) m=0x%lX (%lu) x=%.10e \n",m,m,x ); */
  }

  /* printf("grib_ibm_to_long: s=%lu c=%lu (0x%lX) m=%lu (0x%lX)\n",s,e,e,m,m ); */

  return (s << 31) | ( e << 24 ) | m;
}

double grib_ibmfloat_error(double x) {
  unsigned long e=0;
  
  if (!ibm_table.inited) init_ibm_table();

  if (x < 0)  x = -x;

  /* Underflow */
  if (x <= ibm_table.vmin) return ibm_table.vmin;

  /* Overflow */
  Assert(x <= ibm_table.vmax);

  binary_search(ibm_table.v, 127, x, &e);
  
  return ibm_table.e[e] ;
}

double grib_long_to_ibm(unsigned long x){
  unsigned long s = x  & 0x80000000;
  unsigned long c = (x & 0x7f000000) >> 24;
  unsigned long m = (x & 0x00ffffff);

  double val = m;

  if (!ibm_table.inited) init_ibm_table();

  /*if(x == 0) return 0;*/
  if (c==0 && m <= 1 ) return 0;

  val*=ibm_table.e[c];

  if(s) val = -val;

  return val;
}

double grib_ibm_table_e(unsigned long e) {
  if (!ibm_table.inited) init_ibm_table();
  return ibm_table.e[e];
}

double grib_ibm_table_v(unsigned long e) {
  if (!ibm_table.inited) init_ibm_table();
  return ibm_table.v[e];
}

unsigned long grib_ibm_nearest_smaller_to_long(double x)
{
   unsigned long l;
   unsigned long e;
   unsigned long m ;
   unsigned long s;
   unsigned long mmin = 0x100000;
   double y,eps;

   if(x == 0) return 0;

   if (!ibm_table.inited) init_ibm_table();

   l=grib_ibm_to_long(x);
   y=grib_long_to_ibm(l);

   if ( x < y ) {
     if ( x < 0 && -x < ibm_table.vmin ) {
       l=0x80100000;
     } else {
       e = (l & 0x7f000000) >> 24;
       m = (l & 0x00ffffff);
       s  = l  & 0x80000000;

       if ( m == mmin ) {
         /* printf("grib_ibm_nearest_smaller_to_long: m == mmin (0x%lX) e=%lu\n",m,e);  */
          e = s ? e : e-1;
          if (e<0) e=0;
          if (e>127) e=127;
          /* printf("grib_ibm_nearest_smaller_to_long: e=%lu \n",e);  */
       }

       eps=ibm_table.e[e];

       /* printf("grib_ibm_nearest_smaller_to_long: x<y\n"); */
       l=grib_ibm_to_long(y-eps);
       /* printf("grib_ibm_nearest_smaller_to_long: grib_ibm_to_long(y-eps)=0x%lX y=%.10e eps=%.10e x=%.10e\n",l,y,eps,x); */
     }
   }

   if (x<grib_long_to_ibm(l)) {
     l=grib_ibm_to_long(x-eps);
	 if (x<grib_long_to_ibm(l)) {
		 printf("grib_ibm_nearest_smaller_to_long: x=%.20e grib_long_to_ibm(0x%lX)=%.20e\n",x,l,grib_long_to_ibm(l));
		 Assert(x >= grib_long_to_ibm(l));
	 }
   }

   return l;
}

int grib_nearest_smaller_ibm_float(double a,double* ret) {
  unsigned long l=0;

  if (!ibm_table.inited) init_ibm_table();

  if (a>ibm_table.vmax) return GRIB_INTERNAL_ERROR;

  l=grib_ibm_nearest_smaller_to_long(a);
  *ret=grib_long_to_ibm(l);
  return GRIB_SUCCESS;
}
