/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/**************************************
 *  Enrico Fucile
 **************************************/


#include "grib_api_internal.h"


static void init_dumpers(grib_dumper_class* c,grib_dumper* d)
{
  if(c) {
    grib_dumper_class *s = c->super ? *(c->super) : NULL;
    if(!c->inited)
    {
      if(c->init_class) c->init_class(c);
      c->inited = 1;
    }
    init_dumpers(s,d);
    if(c->init) c->init(d);
  }
}

void grib_init_dumper(grib_dumper* d)
{
  init_dumpers(d->cclass,d);
}

void grib_dumper_delete(grib_dumper* d)
{
  grib_dumper_class *c = d->cclass;
  grib_context *ctx = d->handle->context;
  while(c)
  {
    grib_dumper_class *s = c->super ? *(c->super) : NULL;
    if(c->destroy) c->destroy(d);
    c = s;
  }
  grib_context_free(ctx,d);
}


void  grib_dump_long(grib_dumper* d, grib_accessor* a, const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_long)
    {
      c->dump_long(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_double(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_double)
    {
      c->dump_double(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_string(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_string)
    {
      c->dump_string(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_label(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_label)
    {
      c->dump_label(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_bytes(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_bytes)
    {
      c->dump_bytes(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_bits(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_bits)
    {
      c->dump_bits(d, a, comment);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}


void  grib_dump_section(grib_dumper* d,grib_accessor* a,grib_block_of_accessors* block)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_section)
    {
      c->dump_section(d, a, block);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_values(grib_dumper* d,grib_accessor* a)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->dump_values)
    {
      c->dump_values(d, a);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
  Assert(0);
}

void  grib_dump_header(grib_dumper* d,grib_handle* h)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->header)
    {
      c->header(d,h);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
}

void  grib_dump_footer(grib_dumper* d,grib_handle* h)
{
  grib_dumper_class *c = d->cclass;
  while(c)
  {
    if(c->footer)
    {
      c->footer(d,h);
      return;
    }
    c = c->super ? *(c->super) : NULL;
  }
}
