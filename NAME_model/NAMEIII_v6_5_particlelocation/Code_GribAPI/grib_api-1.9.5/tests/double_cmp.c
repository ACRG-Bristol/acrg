/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: double_cmp
 *
 * Description: compare values retrieved with grib_get_double_array and
 *               grib_get_double_element
 *
 * Author: Cristian D. Codorean
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "grib_api.h"

void usage(char* prog) {
  printf("usage: %s filename\n",prog);
  exit(1);
}

int main(int argc, char** argv) {
  int err = 0,i;
  double *values = NULL;
  size_t values_len= 0;
  double element_value=0;

  FILE* in = NULL;
  char* filename ;
  grib_handle *h = NULL;

  if (argc<2) usage(argv[0]);
  filename=argv[1];

  in = fopen(filename,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",filename);
    return 1;
  }

  /* create new handle from a message in a file*/
  h = grib_handle_new_from_file(0,in,&err);
  if (h == NULL) {
    printf("Error: unable to create handle from file %s\n",filename);
  }


  /* get the size of the values array*/
  GRIB_CHECK(grib_get_size(h,"values",&values_len),0);

  values = malloc(values_len*sizeof(double));

  /* get data values*/
  GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0);

  if (values_len < 100) exit(1);

  for (i = 0; i < 100; i++) {
	  /* get double element */
     GRIB_CHECK(grib_get_double_element(h,"values",i,&element_value),0);

     if (element_value != values[i]) {
	    exit(1);
	 }
  }

  free(values);

  grib_handle_delete(h);

  fclose(in);
  return 0;
}
