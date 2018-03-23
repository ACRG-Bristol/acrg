/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"

/* This file is generated my ./make_class.pl */
#include "grib_dumper_class.h"


struct table_entry
{
  char                  *type;
  grib_dumper_class   **cclass;
};


static struct table_entry table[] =
{
  /* This file is generated my ./make_class.pl */
#include "grib_dumper_factory.h"
};

#define NUMBER(x) (sizeof(x)/sizeof(x[0]))

grib_dumper* grib_dumper_factory(const char* op, grib_handle* h, FILE* out, unsigned long option_flags,void* arg)
{
  int i;
  for(i = 0; i < NUMBER(table) ; i++)
    if(strcmp(op,table[i].type) == 0)
    {
      grib_dumper_class* c = *(table[i].cclass);
      grib_dumper*       d = (grib_dumper*) grib_context_malloc_clear(h->context,c->size);
      d->depth             = 0;
      d->handle            = h;
      d->cclass            = c;
      d->option_flags           = option_flags;
      d->arg               = arg;
      d->out               = out;
      grib_init_dumper(d);
      grib_context_log(h->context,GRIB_LOG_DEBUG,"Creating dumper of type : %s ", op);
      return d;
    }
  grib_context_log(h->context,GRIB_LOG_ERROR,"Unknown type : %s for dumper", op);
  return NULL;
}

void grib_dump_accessors_block(grib_dumper* dumper,grib_block_of_accessors* block)
{
  grib_accessor* a = block->first;
  while(a)
  {
    grib_print_accessor(a,dumper);
    a = a->next;
  }
}

int grib_print       (grib_handle* h, const char* name, grib_dumper *d ){

  grib_accessor* act = grib_find_accessor(h, name);

  if(act){
    grib_print_accessor(act, d );
    return  GRIB_SUCCESS;
  }
  return GRIB_NOT_FOUND;
}

void grib_dump_content(grib_handle* h, FILE* f,const char* mode,unsigned long option_flags,void *data)
{
  grib_dumper *dumper;
  dumper =  grib_dumper_factory(mode?mode:"serialize",h,f,option_flags,data);
  grib_dump_header(dumper,h);
  grib_dump_accessors_block(dumper,h->root->block);
  grib_dump_footer(dumper,h);
  grib_dumper_delete(dumper);

}

