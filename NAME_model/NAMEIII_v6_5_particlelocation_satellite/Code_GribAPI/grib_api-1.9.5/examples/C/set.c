/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: set
 *
 * Description: how to set key values.
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
  long centre=80;
  long long_value=0;
  char string_value[100];
  size_t len = sizeof(string_value)/sizeof(char);
  size_t size=0;

  FILE* in = NULL;
  char* infile = "../../data/regular_latlon_surface.grib1";
  FILE* out = NULL;
  char* outfile = "out.grib1";
  grib_handle *h = NULL;
  const void* buffer = NULL;

  in = fopen(infile,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",infile);
    return 1;
  }

  out = fopen(outfile,"w");
  if(!out) {
    printf("ERROR: unable to open file %s\n",outfile);
    return 1;
  }

  /* create a new handle from a message in a file */
  h = grib_handle_new_from_file(0,in,&err);
  if (h == NULL) {
    printf("Error: unable to create handle from file %s\n",infile);
  }

  /* set centre as a long */
  GRIB_CHECK(grib_set_long(h,"centre",centre),0);

  /* get centre as a long */
  GRIB_CHECK(grib_get_long(h,"centre",&long_value),0);
  printf("centre long value=%ld\n",long_value);

  /* get centre as a string */
  GRIB_CHECK(grib_get_string(h,"centre",string_value,&len),0);
  printf("centre string value=%s\n",string_value);

  /* get the coded message in a buffer */
  GRIB_CHECK(grib_get_message(h,&buffer,&size),0);

  /* write the buffer in a file*/
  if(fwrite(buffer,1,size,out) != size) 
  {
     perror(argv[1]);
     exit(1);
  }

  /* delete handle */
  grib_handle_delete(h);

  fclose(in);
  fclose(out);

  return 0;
}
