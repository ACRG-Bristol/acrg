/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: iterator
 *
 * Description: how to use an iterator on lat/lon/values.
 *
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grib_api.h"  

void usage(char* prog) {
  printf("Usage: %s grib_file\n",prog);
  exit(1);
}

int main(int argc, char** argv) {
  FILE* in = NULL;
  int err = 0;
  double lat,lon,value,missingValue=0;
  int n=0;
  char* filename = NULL;

  /* Message handle. Required in all the grib_api calls acting on a message.*/
  grib_handle *h = NULL;   
  /* Iterator on lat/lon/values.*/
  grib_iterator* iter=NULL;  
  
  if (argc != 2) usage(argv[0]);

  filename=strdup(argv[1]);

  in = fopen(filename,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",filename);
    return 1;
  }

  /* Loop on all the messages in a file.*/
  while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) { 
	/* Check of errors after reading a message. */
    if (err != GRIB_SUCCESS) GRIB_CHECK(err,0);                       

	/* Get the double representing the missing value in the field. */
	GRIB_CHECK(grib_get_double(h,"missingValue",&missingValue),0);  

	/* A new iterator on lat/lon/values is created from the message handle h. */
	iter=grib_iterator_new(h,0,&err);                                     
    if (err != GRIB_SUCCESS) GRIB_CHECK(err,0);                       

    n = 0;
	/* Loop on all the lat/lon/values. */
    while(grib_iterator_next(iter,&lat,&lon,&value)) {   
	  /* You can now print lat and lon,  */
      printf("- %d - lat=%f lon=%f value=",n,lat,lon);   
	  /* decide what to print if a missing value is found. */
      if (value == missingValue ) printf("missing\n");   
	  /* and print the value if is not missing. */
	  else printf("%f\n",value);
      n++;
    }

	/* At the end the iterator is deleted to free memory. */
    grib_iterator_delete(iter);               

    /* At the end the grib_handle is deleted to free memory. */
    grib_handle_delete(h);            
  }


  fclose(in);

  return 0;
}
