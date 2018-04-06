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
 * Description: How to decode grib messages containing multiple
 *              fields. Try to turn on and off multi support to
 *              see the difference. Default is OFF.
 *  			For all the tools defalut is multi support ON.
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
  long parameterCategory=0,parameterNumber=0,discipline=0;
  FILE* in = NULL;
  char* filename = "../../data/multi.grib2";
  grib_handle *h = NULL;

  /* turn on support for multi fields messages */
  grib_multi_support_on(0);

  /* turn off support for multi fields messages */
  /* grib_multi_support_off(0); */

  in = fopen(filename,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",filename);
    return 1;
  }


  while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) {

    GRIB_CHECK(err,0);

    GRIB_CHECK(grib_get_long(h,"discipline",&discipline),0);
    printf("discipline=%ld\n",discipline);

    GRIB_CHECK(grib_get_long(h,"parameterCategory",&parameterCategory),0);
    printf("parameterCategory=%ld\n",parameterCategory);

    GRIB_CHECK(grib_get_long(h,"parameterNumber",&parameterNumber),0);
    printf("parameterNumber=%ld\n",parameterNumber);

    if ( discipline == 0 && parameterCategory==2) {
    if (parameterNumber == 2) printf("-------- u -------\n");
    if (parameterNumber == 3) printf("-------- v -------\n");
    }
  }

  grib_handle_delete(h);

  fclose(in);
  return 0;
}
