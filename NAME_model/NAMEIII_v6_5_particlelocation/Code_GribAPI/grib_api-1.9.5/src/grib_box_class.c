/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 */


#include "grib_api_internal.h"

#define NUMBER(x) (sizeof(x)/sizeof(x[0]))

/* This file is generated my ./make_class.pl */
#include "grib_box_class.h"

struct table_entry
{
  char              *type;
  grib_box_class   **cclass;
};

static struct table_entry table[] =
{
  /* This file is generated my ./make_class.pl */
#include "grib_box_factory.h"
};

grib_box* grib_box_factory( grib_handle* h, grib_arguments* args)
{
  int i;
  int ret = GRIB_SUCCESS;
  char* type = (char*)grib_arguments_get_name(h,args,0);

  for(i = 0; i < NUMBER(table) ; i++)
    if(strcmp(type,table[i].type) == 0)
    {
      grib_box_class* c = *(table[i].cclass);
      grib_box*       it = (grib_box*) grib_context_malloc_clear(h->context,c->size);
      it->cclass              = c;
      ret = GRIB_SUCCESS;
      ret = grib_box_init(it,h,args);
      if(ret ==  GRIB_SUCCESS) return it;
      grib_context_log(h->context,GRIB_LOG_ERROR,"grib_box_factory: error %d instantiating box %s", ret, table[i].type);
      grib_box_delete(it);
      return NULL;
    }

  grib_context_log(h->context,GRIB_LOG_ERROR,"grib_box_factory : Unknown type : %s for box", type);

  return NULL;
}
