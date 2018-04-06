/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: set_bitmap
 *
 * Description: how to set a bitmap in a grib message
 *
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "grib_api.h"

int main(int argc, char** argv) {
  int err = 0;
  size_t size=0;

  FILE* in = NULL;
  char* infile = "../data/regular_latlon_surface.grib1";
  FILE* out = NULL;
  char* outfile = "out.grib1";
  grib_handle *h = NULL;
  const void* buffer = NULL;
  size_t values_len;
  double* values;
  double missing=9999;
  int i=0;

  in = fopen(infile,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",infile);
    return 1;
  }

  out = fopen(outfile,"w");
  if(!in) {
    printf("ERROR: unable to open file %s\n",outfile);
    return 1;
  }

  h = grib_handle_new_from_file(0,in,&err);
  if (h == NULL) {
    printf("Error: unable to create handle from file %s\n",infile);
  }

  GRIB_CHECK(grib_set_double(h,"missingValue",missing),0);

  /* get the size of the values array*/
  GRIB_CHECK(grib_get_size(h,"values",&values_len),0);

  values = malloc(values_len*sizeof(double));

  /* get data values*/
  GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0);

  GRIB_CHECK(grib_set_long(h,"bitmapPresent",1),0);

  for(i = 0; i < 10; i++)
    values[i]=missing;

  GRIB_CHECK(grib_set_double_array(h,"values",values,values_len),0);

  /* get the coded message in a buffer */
  GRIB_CHECK(grib_get_message(h,&buffer,&size),0);

  /* write the buffer in a file*/
  if(fwrite(buffer,1,size,out) != size)
  {
     perror(outfile);
     exit(1);
  }

  /* delete handle */
  grib_handle_delete(h);

  fclose(in);
  fclose(out);

  return 0;
}
