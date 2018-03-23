/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: keys_iterator
 *
 * Description:
 * Example on how to use keys_iterator functions and the
 * grib_keys_iterator structure to get all the available
 * keys in a message.
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "grib_api.h"

#define MAX_KEY_LEN  255
#define MAX_VAL_LEN  1024

static void usage(char* progname);

int main(int argc, char *argv[])
{
  /* To skip read only and not coded keys
     unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_SKIP_READ_ONLY ||
     GRIB_KEYS_ITERATOR_SKIP_COMPUTED;
  */
  unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;

  /* valid name_spaces are ls and mars */
  char* name_space="ls";

  /* name_space=NULL to get all the keys */
  /* char* name_space=0; */

  FILE* f;
  grib_handle* h=NULL;
  grib_keys_iterator* kiter=NULL;
  int err=0;
  int grib_count=0;

  char value[MAX_VAL_LEN];
  size_t vlen=MAX_VAL_LEN;

  if (argc != 2) usage(argv[0]);

  f = fopen(argv[1],"r");
  if(!f) {
    perror(argv[1]);
    exit(1);
  }

  while((h = grib_handle_new_from_file(0,f,&err)) != NULL) {

    grib_count++;
    printf("-- GRIB N. %d --\n",grib_count);
    if(!h) {
      printf("ERROR: Unable to create grib handle\n");
      exit(1);
    }

    kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
    if (!kiter) {
      printf("ERROR: Unable to create keys iterator\n");
      exit(1);
    }

    while(grib_keys_iterator_next(kiter))
    {
      const char* name = grib_keys_iterator_get_name(kiter);
      vlen=MAX_VAL_LEN;
      bzero(value,vlen);
      GRIB_CHECK(grib_get_string(h,name,value,&vlen),name);
      printf("%s = %s\n",name,value);
    }

    grib_keys_iterator_delete(kiter);

  }

  return 0;

}

static void usage(char* progname) {
  printf("\nUsage: %s grib_file\n",progname);
  exit(1);
}

