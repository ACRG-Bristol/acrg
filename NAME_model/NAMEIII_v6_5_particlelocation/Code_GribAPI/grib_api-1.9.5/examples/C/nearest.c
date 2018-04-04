/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: fieldset
 *
 * Description: how to use a fieldset.
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
  printf("Usage: %s grib_file grib_file ...\n",prog);
  exit(1);
}

int main(int argc, char** argv) {
  int err = 0;
  long step=0;
  size_t nfiles;
  int i=0;
  grib_fieldset* set=NULL;
  grib_handle* h=NULL;
  char param[20]={0,};
  size_t len=20;
  double lats[4]={0,};
  double lons[4]={0,};
  double values[4]={0,};
  double distances[4]={0,};
  int indexes[4]={0,};
  char* order_by="param,step";

  size_t size=4;
  double lat=-40,lon=15;
  int mode=0;
  int count;
  char** filenames;
  grib_nearest* nearest=NULL;

  if (argc < 2) usage(argv[0]);

  nfiles=argc-1;
  filenames=(char**)malloc(sizeof(char*)*nfiles);
  for (i=0;i<nfiles;i++)
    filenames[i]=(char*)strdup(argv[i+1]);

  set=grib_fieldset_new_from_files(0,filenames,nfiles,0,0,0,order_by,&err);
  GRIB_CHECK(err,0);

  printf("\nordering by %s\n",order_by);
  printf("\n%d fields in the fieldset\n",grib_fieldset_count(set));
  printf("n,step,param\n");

  mode=GRIB_NEAREST_SAME_GRID |  GRIB_NEAREST_SAME_POINT;
  count=1;
  while ((h=grib_fieldset_next_handle(set,&err))!=NULL) {
    GRIB_CHECK(grib_get_long(h,"step",&step),0);
	len=20;
    GRIB_CHECK(grib_get_string(h,"param",param,&len),0);

    printf("%d %ld %s  ",count,step,param);
    if (!nearest) nearest=grib_nearest_new(h,&err);
    GRIB_CHECK(err,0);
    GRIB_CHECK(grib_nearest_find(nearest,h,lat,lon,mode,lats,lons,values,distances,indexes,&size),0);
    for (i=0;i<4;i++) printf("%d %.2f %.2f %g %g - ",
         (int)indexes[i],lats[i],lons[i],distances[i],values[i]);
    printf("\n");

    grib_handle_delete(h);
    count++;
  }

  if (nearest) grib_nearest_delete(nearest);

  if (set) grib_fieldset_delete(set);

  return 0;
}
