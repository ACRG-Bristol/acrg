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
#include <ctype.h>
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = dumper
   IMPLEMENTS = dump_long;dump_bits
   IMPLEMENTS = dump_double;dump_string
   IMPLEMENTS = dump_bytes;dump_values
   IMPLEMENTS = dump_label;dump_section
   IMPLEMENTS = init;destroy
   MEMBERS = long section_offset
   MEMBERS = long begin
   MEMBERS = long end
   END_CLASS_DEF

 */


/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "dumper.class" and rerun ./make_class.pl

*/

static void init_class      (grib_dumper_class*);
static int init            (grib_dumper* d);
static int destroy         (grib_dumper*);
static void dump_long       (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_bits       (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_double     (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_string     (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_bytes      (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_values     (grib_dumper* d, grib_accessor* a);
static void dump_label      (grib_dumper* d, grib_accessor* a,const char* comment);
static void dump_section    (grib_dumper* d, grib_accessor* a,grib_block_of_accessors* block);

typedef struct grib_dumper_wmo {
    grib_dumper          dumper;  
/* Members defined in wmo */
	long section_offset;
	long begin;
	long end;
} grib_dumper_wmo;


static grib_dumper_class _grib_dumper_class_wmo = {
    0,                              /* super                     */
    "wmo",                              /* name                      */
    sizeof(grib_dumper_wmo),     /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    &init,                               /* init                      */
    &destroy,                            /* free mem                       */
    &dump_long,                          /* dump long         */
    &dump_double,                        /* dump double    */
    &dump_string,                        /* dump string    */
    &dump_label,                         /* dump labels  */
    &dump_bytes,                         /* dump bytes  */
    &dump_bits,                          /* dump bits   */
    &dump_section,                       /* dump section      */
    &dump_values,                        /* dump values   */
    0,                             /* header   */
    0,                             /* footer   */
};

grib_dumper_class* grib_dumper_class_wmo = &_grib_dumper_class_wmo;

/* END_CLASS_IMP */
static void set_begin_end(grib_dumper* d,grib_accessor* a);
static void print_offset(FILE* out,long begin,long end);
static void print_hexadecimal(FILE* out,unsigned long flags,grib_accessor* a);

static void init_class      (grib_dumper_class* c){}

static int  init(grib_dumper* d)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  self->section_offset=0;

  return GRIB_SUCCESS;
}

static int  destroy  (grib_dumper* d){
  return GRIB_SUCCESS;
}


static void aliases(grib_dumper* d,grib_accessor* a)
{
int i;
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;

  if( (d->option_flags & GRIB_DUMP_FLAG_ALIASES) == 0)
    return;

  if(a->all_names[1])
  {
    char *sep = "";
    fprintf(self->dumper.out," [");

    for(i = 1; i < MAX_ACCESSOR_NAMES; i++)
    {
      if(a->all_names[i])
      {
        if(a->all_name_spaces[i])
          fprintf(self->dumper.out,"%s%s.%s", sep,a->all_name_spaces[i],a->all_names[i]);
        else
          fprintf(self->dumper.out,"%s%s", sep,a->all_names[i]);
      }
    sep = ", ";
    }
    fprintf(self->dumper.out,"]");
  }
}

static void dump_long(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  long value; size_t size = 1;
  int err = grib_unpack_long(a,&value,&size);

  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  if( (a->flags & GRIB_ACCESSOR_FLAG_READ_ONLY) != 0 &&
      (d->option_flags & GRIB_DUMP_FLAG_READ_ONLY) == 0)
    return;

  set_begin_end(d,a);

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/

  print_offset(self->dumper.out,self->begin,self->end);

  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  if( ((a->flags & GRIB_ACCESSOR_FLAG_CAN_BE_MISSING) != 0) && grib_is_missing_internal(a) )
    fprintf(self->dumper.out,"%s = MISSING",a->name);
  else
    fprintf(self->dumper.out,"%s = %ld",a->name,value);

  print_hexadecimal(self->dumper.out,d->option_flags,a);

  if(comment) fprintf(self->dumper.out," [%s]",comment);
  if(err)
    fprintf(self->dumper.out," *** ERR=%d (%s)",err,grib_get_error_message(err));

  aliases(d,a);


  fprintf(self->dumper.out,"\n");
}

static int test_bit(long a, long b) {return a&(1<<b);}

static void dump_bits(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  int i;
  long value; size_t size = 1;
  int err = grib_unpack_long(a,&value,&size);


  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  set_begin_end(d,a);

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  print_offset(self->dumper.out,self->begin,self->end);
  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  fprintf(self->dumper.out,"%s = %ld [",a->name,value);

  for(i=0;i<(a->length*8);i++) {
    if(test_bit(value,a->length*8-i-1))
      fprintf(self->dumper.out,"1");
    else
      fprintf(self->dumper.out,"0");
  }
/*
  if(comment)
    fprintf(self->dumper.out,":%s]",comment);
  else
*/
  fprintf(self->dumper.out,"]");

  if (err==0) print_hexadecimal(self->dumper.out,d->option_flags,a);

  if(err)
    fprintf(self->dumper.out," *** ERR=%d (%s)",err,grib_get_error_message(err));

  aliases(d,a);
  fprintf(self->dumper.out,"\n");
}

static void dump_double(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  double value; size_t size = 1;
  int err = grib_unpack_double(a,&value,&size);


  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  set_begin_end(d,a);

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/

  print_offset(self->dumper.out,self->begin,self->end);
  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  if( ((a->flags & GRIB_ACCESSOR_FLAG_CAN_BE_MISSING) != 0) && grib_is_missing_internal(a))
    fprintf(self->dumper.out,"%s = MISSING",a->name);
  else
    fprintf(self->dumper.out,"%s = %g",a->name,value);
  /*if(comment) fprintf(self->dumper.out," [%s]",comment);*/

  if (err==0) print_hexadecimal(self->dumper.out,d->option_flags,a);

  if(err)
    fprintf(self->dumper.out," *** ERR=%d (%s)",err,grib_get_error_message(err));
  aliases(d,a);
  fprintf(self->dumper.out,"\n");
}

static void dump_string(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  char value[1024]; size_t size = sizeof(value);
  int err = grib_unpack_string(a,value,&size);
  char *p = value;

  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  set_begin_end(d,a);

  while(*p) { if(!isprint(*p)) *p = '.'; p++; }

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  print_offset(self->dumper.out,self->begin,self->end);
  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  fprintf(self->dumper.out,"%s = %s",a->name,value);

  if (err==0) print_hexadecimal(self->dumper.out,d->option_flags,a);

  /*if(comment) fprintf(self->dumper.out," [%s]",comment);*/
  if(err)
    fprintf(self->dumper.out," *** ERR=%d (%s)",err,grib_get_error_message(err));
  aliases(d,a);
  fprintf(self->dumper.out,"\n");
}

static void dump_bytes(grib_dumper* d,grib_accessor* a,const char* comment)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  int i,k,err =0;
  int more = 0;
  size_t size = a->length;
  unsigned char* buf = grib_context_malloc(d->handle->context,size);

  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  set_begin_end(d,a);

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  print_offset(self->dumper.out,self->begin,self->end);
  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  fprintf(self->dumper.out,"%s = %ld",a->name,a->length);
  aliases(d,a);
  fprintf(self->dumper.out," {");

  if(!buf)
  {
    if(size == 0)
      fprintf(self->dumper.out,"}\n");
    else
      fprintf(self->dumper.out," *** ERR cannot malloc(%ld) }\n",(long)size);
    return;
  }

  print_hexadecimal(self->dumper.out,d->option_flags,a);

  fprintf(self->dumper.out,"\n");

  err = grib_unpack_bytes(a,buf,&size);
  if(err){
    grib_context_free(d->handle->context,buf);
    fprintf(self->dumper.out," *** ERR=%d (%s) \n}",err,grib_get_error_message(err));
    return ;
  }

  if(size > 100) {
    more = size - 100;
    size = 100;
  }

  k = 0;
  /* if(size > 100) size = 100;  */
  while(k < size)
  {
    int j;
    for(i = 0; i < d->depth + 3 ; i++) fprintf(self->dumper.out," ");
    for(j = 0; j < 16 && k < size; j++, k++)
    {
      fprintf(self->dumper.out,"%02x",buf[k]);
      if(k != size-1)
        fprintf(self->dumper.out,", ");
    }
    fprintf(self->dumper.out,"\n");
  }

  if(more)
  {
    for(i = 0; i < d->depth + 3 ; i++) fprintf(self->dumper.out," ");
    fprintf(self->dumper.out,"... %d more values\n",more);
  }

  for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");
  fprintf(self->dumper.out,"} # %s %s \n",a->creator->op, a->name);
  grib_context_free(d->handle->context,buf);
}

static void dump_values(grib_dumper* d,grib_accessor* a)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  int k,err =0;
  int more = 0;
  double*  buf = NULL;
  size_t size=0;

  if( a->length == 0  &&
      (d->option_flags & GRIB_DUMP_FLAG_CODED) != 0)
    return;

  size=grib_value_count(a);

  if(size == 1){
    dump_double(d,a,NULL);
    return ;
  }
  buf = grib_context_malloc(d->handle->context,size * sizeof(double));

  set_begin_end(d,a);

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  print_offset(self->dumper.out,self->begin,self->end);
  if ((d->option_flags & GRIB_DUMP_FLAG_TYPE) != 0)
        fprintf(self->dumper.out,"%s ",a->creator->op);

  fprintf(self->dumper.out,"%s = (%ld,%ld)",a->name,(long)size,a->length);
  aliases(d,a);
  fprintf(self->dumper.out," {");

  if(!buf)
  {
    if(size == 0)
      fprintf(self->dumper.out,"}\n");
    else
      fprintf(self->dumper.out," *** ERR cannot malloc(%ld) }\n",(long)size);
    return;
  }

  fprintf(self->dumper.out,"\n");

  err =  grib_unpack_double(a,buf,&size);

  if(err){
    grib_context_free(d->handle->context,buf);
    fprintf(self->dumper.out," *** ERR=%d (%s) \n}",err,grib_get_error_message(err));
    return ;
  }

  if(size > 100) {
    more = size - 100;
    size = 100;
  }


  k = 0;
  while(k < size)
  {
#if 1
    int j;
    /*for(i = 0; i < d->depth + 3 ; i++) fprintf(self->dumper.out," ");*/
    for(j = 0; j < 8 && k < size; j++, k++)
    {
      fprintf(self->dumper.out,"%.10e",buf[k]);
      if(k != size-1)
        fprintf(self->dumper.out,", ");
    }
    fprintf(self->dumper.out,"\n");
#else

    fprintf(self->dumper.out,"%d %g\n",k,buf[k]);

#endif

  }
  if(more)
  {
    /*for(i = 0; i < d->depth + 3 ; i++) fprintf(self->dumper.out," ");*/
    fprintf(self->dumper.out,"... %d more values\n",more);
  }

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  fprintf(self->dumper.out,"} # %s %s \n",a->creator->op, a->name);
  grib_context_free(d->handle->context,buf);
}

static void dump_label(grib_dumper* d,grib_accessor* a,const char* comment)
{
  /*grib_dumper_wmo *self = (grib_dumper_wmo*)d;

  for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");
  fprintf(self->dumper.out,"----> %s %s %s\n",a->creator->op, a->name,comment?comment:"");*/
}

static void dump_section(grib_dumper* d,grib_accessor* a,grib_block_of_accessors* block)
{
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  grib_section* s = a->sub_section;
  int is_wmo_section=0;
  char* upper=NULL;
  char tmp[512];
  char *p=NULL,*q=NULL;
  if (!strncmp(a->name,"section",7)) is_wmo_section=1;

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  if (is_wmo_section) {
    upper=(char*)malloc(strlen(a->name)+1);
    p=(char*)a->name;
    q=upper;
    while (*p != '\0') {
      *q=toupper(*p);
      q++;
      p++;
    }
    *q='\0';
    sprintf(tmp,"%s ( length=%ld, padding=%ld )",upper,(long)s->length,(long)s->padding);
    fprintf(self->dumper.out,"======================   %-35s   ======================\n",tmp);
    free(upper);
    self->section_offset=a->offset;
  } else {

  }

  /*printf("------------- section_offset = %ld\n",self->section_offset);*/
  d->depth += 3;
  grib_dump_accessors_block(d,block);
  d->depth -= 3;

  /*for(i = 0; i < d->depth ; i++) fprintf(self->dumper.out," ");*/
  /*fprintf(self->dumper.out,"<===== %s %s\n",a->creator->op, a->name);*/
}

static void set_begin_end(grib_dumper* d,grib_accessor* a) {
  grib_dumper_wmo *self = (grib_dumper_wmo*)d;
  if ((d->option_flags & GRIB_DUMP_FLAG_OCTECT) != 0) {

    self->begin=a->offset-self->section_offset+1;
    self->end=grib_get_next_position_offset(a)-self->section_offset;
  } else {
    self->begin=a->offset;
    self->end=grib_get_next_position_offset(a);
  }
}

static void print_offset(FILE* out,long begin,long end) {
  char tmp[50];
  if (begin == end)
    fprintf(out,"%-10ld" ,begin);
  else {
    sprintf(tmp,"%ld-%ld" ,begin,end);
    fprintf(out,"%-10s",tmp);
  }
}

static void print_hexadecimal(FILE* out,unsigned long flags,grib_accessor* a) {
  int i=0;
  unsigned long offset=0;
  if ((flags & GRIB_DUMP_FLAG_HEXADECIMAL) != 0 && a->length != 0) {
    fprintf(out," (");
    offset=a->offset;
    for (i=0;i<a->length;i++) {
      fprintf(out," 0x%.2X",a->parent->h->buffer->data[offset]);
      offset++;
    }
    fprintf(out," )");
  }
}
