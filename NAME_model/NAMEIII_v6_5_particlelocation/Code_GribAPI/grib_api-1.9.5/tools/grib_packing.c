/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_debug
 *
 * Description: 
 *
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "grib_api_internal.h"
void usage(const char* p)
{
	printf("Usage: %s [-n minimum] [-x maximum] [-b number of bits] [-d decimal scaling] [-2] [--] [values ..]\n",p);
	exit(1);
}

unsigned long pack(double v,double ref,long bscale,long dscale)
{
	double b = grib_power(-bscale,2);
	double d = grib_power( dscale,10);
	return (v*d - ref)*b + 0.5;
}

double unpack(unsigned long v,double ref,long bscale,long dscale)
{
	double b = grib_power( bscale,2);
	double d = grib_power(-dscale,10);
	return (v*b+ref)*d;
}

void try(const char* txt,double v,double ref,long bscale,long dscale)
{
	double x;
	if(txt) printf("%-18s  : ",txt);
	x = unpack(pack(v,ref,bscale,dscale),ref,bscale,dscale);
	printf("%g (err=%g)\n",x,v-x); 
}

int main(int argc, char *argv[])
{
	extern char *optarg;
	extern int optind, opterr, optopt;
	double min = 0,max = 100;
	long bits = 16;
	int g2 = 0;
	long dscale = 0;
	double ref;
	double decimal_scale;
	long bscale;
	double dmin,dmax;
	int err=0;

	int c;

	while((c = getopt(argc, argv,"2n:x:b:d:")) != EOF) {
		switch(c) {

			case 'n':
				min = atof(optarg);
				break;

			case 'x':
				max = atof(optarg);
				break;

			case 'b':
				bits = atol(optarg);
				break;

			case 'd':
				dscale = atol(optarg);
				break;

			case '2':
				g2=1;
				break;

			default:
				usage(argv[0]);
				break;
		}
	}

	decimal_scale = grib_power(dscale,10) ;

	dmin = min*decimal_scale;
	dmax = max*decimal_scale;

	if(g2) {
    double* pref=&ref;
    grib_nearest_smaller_ieee_float(dmin,pref);
  }
	else
	  grib_nearest_smaller_ibm_float(dmin,&ref);

	bscale = grib_get_binary_scale_fact(dmax,ref,bits,&err);

	printf("\n");

	printf("Minimum value       : %g\n",min);
	printf("Maximum value       : %g\n",max);
	printf("Bits per values     : %ld\n",bits);
	printf("Binary scale factor : %ld\n",bscale);
	printf("Decimal scale factor: %ld\n",dscale);
	printf("Encoding            : %s\n",g2 ? "IEEE" : "IBM");
	printf("\n");

	if(ref != dmin)
	{
		printf("WARNING: reference value != scaled minimum value\n");
		printf("         reference value = %g\n",ref);
		printf("         scaled minimum  = %g\n",dmin);
		printf("         difference      = %g\n",ref-dmin);
		printf("\n");
	}




	printf("Packing range       : [%0*lx - %0*lx]\n",
		(int)(bits+7)/8*2,pack(min,ref,bscale,dscale),
		(int)(bits+7)/8*2,pack(max,ref,bscale,dscale));

	printf("Bucket size         : %g\n",unpack(1,ref,bscale,dscale) - unpack(0,ref,bscale,dscale));

	printf("\n");
	printf("Packing/unpacking   :\n");
	printf("\n");
	try("Minimum",min,ref,bscale,dscale);
	try("Maximum",max,ref,bscale,dscale);

	for(;optind < argc; optind++)
		try(argv[optind],atof(argv[optind]),ref,bscale,dscale);


	return 0;
}


