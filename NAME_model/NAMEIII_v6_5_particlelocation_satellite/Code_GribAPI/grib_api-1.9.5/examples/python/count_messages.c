/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: get
 *
 * Description: how to get values using keys.
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
  int err = 0;

  size_t i = 0;
  size_t nmsg = 0;

  double latitudeOfFirstGridPointInDegrees;
  double longitudeOfFirstGridPointInDegrees;
  double latitudeOfLastGridPointInDegrees;
  double longitudeOfLastGridPointInDegrees;

  double jDirectionIncrementInDegrees;
  double iDirectionIncrementInDegrees;

  long numberOfPointsAlongAParallel;
  long numberOfPointsAlongAMeridian;

  size_t values_len= 0;
  double average;
  double min;
  double max;

  FILE* in = NULL;
  char* filename;
  grib_handle *h = NULL;

  if (argc<2) usage(argv[0]);
  filename=argv[1];

  in = fopen(filename,"r");
  if(!in) {
    printf("ERROR: unable to open file %s\n",filename);
    return 1;
  }

  while((h = grib_handle_new_from_file(0,in,&err)) != NULL) {

    GRIB_CHECK(err,0);

    nmsg++;
    printf("processing message number %d\n",nmsg);

    /* get as a long*/
    GRIB_CHECK(grib_get_long(h,"numberOfPointsAlongAParallel",&numberOfPointsAlongAParallel),0);
    printf("numberOfPointsAlongAParallel=%ld\n",numberOfPointsAlongAParallel);

    /* get as a long*/
    GRIB_CHECK(grib_get_long(h,"numberOfPointsAlongAMeridian",&numberOfPointsAlongAMeridian),0);
    printf("numberOfPointsAlongAMeridian=%ld\n",numberOfPointsAlongAMeridian);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"latitudeOfFirstGridPointInDegrees",&latitudeOfFirstGridPointInDegrees),0);
    printf("latitudeOfFirstGridPointInDegrees=%g\n",latitudeOfFirstGridPointInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"longitudeOfFirstGridPointInDegrees",&longitudeOfFirstGridPointInDegrees),0);
    printf("longitudeOfFirstGridPointInDegrees=%g\n",longitudeOfFirstGridPointInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"latitudeOfLastGridPointInDegrees",&latitudeOfLastGridPointInDegrees),0);
    printf("latitudeOfLastGridPointInDegrees=%g\n",latitudeOfLastGridPointInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"longitudeOfLastGridPointInDegrees",&longitudeOfLastGridPointInDegrees),0);
    printf("longitudeOfLastGridPointInDegrees=%g\n",longitudeOfLastGridPointInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"jDirectionIncrementInDegrees",&jDirectionIncrementInDegrees),0);
    printf("jDirectionIncrementInDegrees=%g\n",jDirectionIncrementInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"iDirectionIncrementInDegrees",&iDirectionIncrementInDegrees),0);
    printf("iDirectionIncrementInDegrees=%g\n",iDirectionIncrementInDegrees);

    /* get as a double*/
    GRIB_CHECK(grib_get_double(h,"average",&average),0);
    GRIB_CHECK(grib_get_double(h,"min",&min),0);
    GRIB_CHECK(grib_get_double(h,"max",&max),0);

    /* get the size of the values array*/
    GRIB_CHECK(grib_get_size(h,"values",&values_len),0);

    printf("There are %d, average is %g, min is %g, max is %g\n",(int)values_len,average,min,max);

    for (i=0;i<100;i++) printf("-");
    printf("\n");

    grib_handle_delete(h);

  }

  fclose(in);
  return 0;
}
