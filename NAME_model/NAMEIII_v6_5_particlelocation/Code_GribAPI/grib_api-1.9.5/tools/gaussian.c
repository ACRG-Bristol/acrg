/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "grib_api.h"

void usage(const char *prog)
{
	fprintf(stderr,"%s: N\n",prog);
	exit(1);
}
int main(int argc, const char *argv[])
{
	double *latitudes;
	int count = 0;
	int n,i;


	if(argc != 2) usage(argv[0]);

	n = atol(argv[1]);
	count = n*2;

	latitudes = (double*)malloc(count*sizeof(double));
	assert(latitudes);

	assert(grib_get_gaussian_latitudes(n,latitudes) == 0);

	for(i = 0; i < count ; i++)
		printf("%g\n",latitudes[i]);

	return 0;
}
