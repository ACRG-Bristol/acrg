/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: big2gribex
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

#include "grib_api.h"
void usage(char*);

int main(int argc, char *argv[])
{
  grib_handle* h = NULL;
  FILE* inf = NULL;
  FILE* ouf = NULL;
  char* infile=0;
  char* oufile=0;
  int err = 0;
  size_t size;
  const void* buffer=NULL;

  if (argc < 3) usage(argv[0]);
  infile=argv[1];
  oufile=argv[2];

  inf = fopen(infile,"r");
  if(!inf) {
     perror(infile);
     exit(1);
  }

  ouf = fopen(oufile,"w");
  if(!ouf) {
     perror(oufile);
     exit(1);
  }

  while((h = grib_handle_new_from_file(0,inf,&err)) != NULL) {
      grib_update_sections_lengths(h);
	  GRIB_CHECK(grib_get_message(h,&buffer,&size),0);
	  if(fwrite(buffer,1,size,ouf) != size) {
        perror(oufile);
		exit(1);
	  }
  }

  grib_handle_delete(h);

  if(fclose(inf)) {
     perror(infile);
     exit(1);
  }

  if(fclose(ouf)) {
     perror(oufile);
     exit(1);
  }


  return 0;
}

void usage(char* progname) {
  printf("\nUsage: %s grib_in grib_out\n",progname);
  exit(1);
}
