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
  printf("usage: %s filename\n",prog);
  exit(1);
}

int main(int argc, char** argv) {
  int err = 0,i;
  double *values = NULL;
  double max,min,average;
  size_t values_len= 0;

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

  for(i = 0; i < values_len; i++)
    printf("%d %.10e\n",i+1,values[i]);

  free(values);


  GRIB_CHECK(grib_get_double(h,"max",&max),0);
  GRIB_CHECK(grib_get_double(h,"min",&min),0);
  GRIB_CHECK(grib_get_double(h,"average",&average),0);

  printf("%d values found in %s\n",(int)values_len,filename);
  printf("max=%.10e\nmin=%.10e\naverage=%.10e\n",max,min,average);

  grib_handle_delete(h);

  fclose(in);
  return 0;
}
