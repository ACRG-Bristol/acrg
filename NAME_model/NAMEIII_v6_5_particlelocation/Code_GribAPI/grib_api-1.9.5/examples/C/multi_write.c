/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: multi  
 *
 * Description: How to encode grib messages containing multiple
 *              fields. 
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "grib_api.h"

void usage(char* prog) {
  printf("usage: %s in.grib out.grib\n",prog);
  exit(1);
}

int main(int argc, char** argv) {
  int err = 0;
  FILE* in = NULL;
  FILE* of = NULL;
  long step;
  char* filename=NULL;
  char* ofilename=NULL;
  grib_handle *h = NULL;
  grib_multi_handle *mh=NULL;

  if (argc < 3) usage(argv[0]);
  filename=argv[1];
  ofilename=argv[2];

  /* open input file */
  in = fopen(filename,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",filename);
    return 1;
  }

  /* new grib handle from input file */
  h = grib_handle_new_from_file(0,in,&err);  
  GRIB_CHECK(err,0);
  /* create a new empty multi field handle */
  mh=grib_multi_handle_new(0);
  if (!mh) {
    printf("unable to create multi field handle\n");
	exit(1);
  }

  for (step=12;step<=120;step+=12) {
    /* set step */
    grib_set_long(h,"step",step);
	/* append h to mh repeating the sections preceding section 4 */
	grib_multi_handle_append(h,4,mh);
  }

  /* open output file */
  of=fopen(ofilename,"w");
  if(!of) {
	printf("ERROR: unable to open file %s\n",ofilename);
    exit(1);
  }

  /* write multi fields handle to output file */
  grib_multi_handle_write(mh,of);
  fclose(of);

  /* clean memory */
  grib_handle_delete(h);
  grib_multi_handle_delete(mh);

  fclose(in);
  return 0;
}
