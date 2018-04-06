/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: set_missing
 *
 * Description: how to set a key as missing.
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

  FILE* in = NULL;
  char* infile = "../../data/reduced_gaussian_pressure_level.grib2";
  FILE* out = NULL;
  char* outfile = "out_surface_level.grib2";
  grib_handle *h = NULL;
  const void* buffer = NULL;
  size_t size=0;
  char str[]="sfc";
  size_t str_len=3;

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

  /* create a new handle from a message in a file */
  h = grib_handle_new_from_file(0,in,&err);
  if (h == NULL) {
    printf("Error: unable to create handle from file %s\n",infile);
  }

  GRIB_CHECK(grib_set_string(h,"typeOfFirstFixedSurface",str,&str_len),0);
  GRIB_CHECK(grib_set_missing(h,"scaleFactorOfFirstFixedSurface"),0);
  GRIB_CHECK(grib_set_missing(h,"scaledValueOfFirstFixedSurface"),0);

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
