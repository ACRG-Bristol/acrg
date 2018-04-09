/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: print_data
 *
 * Description: prints all the data contained in a grib file
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "grib_api.h"

void usage(char* prog) {
	printf("usage: %s out.grib\n",prog);
	exit(1);
}

int main(int argc, char** argv) {
  int i;
  double *values = NULL;
  size_t values_len= 0;

  char* filename ;
  grib_handle *h = NULL;
  double d,e;
  long count;

  if (argc!=2) usage(argv[0]);
  /* create new handle from a message in a file*/
  h = grib_handle_new_from_template(0,"regular_ll_pl_grib1");
  if (h == NULL) {
    printf("Error: unable to create handle from file %s\n",filename);
	exit(1);
  }

  values_len=10000;
  values = malloc(values_len*sizeof(double));
  d=10e-8;
  e=d;
  count=1;
  for (i=0;i<values_len;i++) {
 	if (count>1000) {e*=10; count=1;}
	values[i]=d;
	printf("%g \n",values[i]);
	d+=e;
	count++;
  }

  GRIB_CHECK(grib_set_long(h,"bitsPerValue",16),0);

  /* get data values*/
  GRIB_CHECK(grib_set_double_array(h,"values",values,values_len),0);

  grib_write_message(h,argv[1],"w");

  free(values);
  grib_handle_delete(h);

  return 0;
}
