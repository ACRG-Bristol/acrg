/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/
/***********************************************************
 *  Enrico Fucile
 ***********************************************************/

#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_gen
   IMPLEMENTS = unpack_string;pack_string
   IMPLEMENTS = unpack_long;pack_long
   IMPLEMENTS = unpack_double;pack_double
   IMPLEMENTS = init;dump
   IMPLEMENTS = value_count
   IMPLEMENTS = get_native_type
   IMPLEMENTS = compare
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int  get_native_type(grib_accessor*);
static int pack_double(grib_accessor*, const double* val,size_t *len);
static int pack_long(grib_accessor*, const long* val,size_t *len);
static int pack_string(grib_accessor*, const char*, size_t *len);
static int unpack_double(grib_accessor*, double* val,size_t *len);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long value_count(grib_accessor*);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static int compare(grib_accessor*, grib_accessor*);

typedef struct grib_accessor_ascii {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in ascii */
} grib_accessor_ascii;

extern grib_accessor_class* grib_accessor_class_gen;

static grib_accessor_class _grib_accessor_class_ascii = {
    &grib_accessor_class_gen,                      /* super                     */
    "ascii",                      /* name                      */
    sizeof(grib_accessor_ascii),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    &dump,                       /* describes himself         */
    0,                /* get length of section     */
    &value_count,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    &get_native_type,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    &pack_long,                  /* grib_pack procedures long      */
    &unpack_long,                /* grib_unpack procedures long    */
    &pack_double,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
    &pack_string,                /* grib_pack procedures string    */
    &unpack_string,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    0,            /* pack_expression */
    0,              /* notify_change   */
    0,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    &compare,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_ascii = &_grib_accessor_class_ascii;


static void init_class(grib_accessor_class* c)
{
	c->next_offset	=	(*(c->super))->next_offset;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
	c->update_size	=	(*(c->super))->update_size;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a, const long len , grib_arguments* arg )
{
  a->length = len;
  Assert(a->length>=0);
}

static long value_count(grib_accessor* a){
  return a->length +1;
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
  grib_dump_string(dumper,a,NULL);
}

static int  get_native_type(grib_accessor* a){
  return GRIB_TYPE_STRING;
}

static int unpack_string(grib_accessor* a, char* val, size_t *len)
{

  int i = 0;

  if(len[0] < (a->length+1))
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "unpack_string: Wrong size (%d) for %s it contains %d values ", len[0], a->name , a->length+1 );
    len[0] = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  for ( i = 0; i < a->length; i++)
    val[i] = a->parent->h->buffer->data[a->offset+i];
  val[i] = 0;
  len[0] = i;
  return GRIB_SUCCESS;
}

static int pack_string(grib_accessor* a, const char* val, size_t *len)
{

  int i = 0;
  if(len[0] > (a->length)+1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "pack_string: Wrong size (%d) for %s it contains %d values ", len[0], a->name , a->length+1 );
    len[0] = 0;
    return GRIB_BUFFER_TOO_SMALL;
  }

  for ( i = 0; i < a->length; i++)
  {
    if( i < len[0] )
      a->parent->h->buffer->data[a->offset+i] = val[i];
    else
      a->parent->h->buffer->data[a->offset+i] = 0;
  }

  return GRIB_SUCCESS;
}

static int pack_long(grib_accessor* a, const long*  v, size_t *len){
  grib_context_log(a->parent->h->context,GRIB_LOG_ERROR, " Should not pack %s as long", a->name);
  return GRIB_NOT_IMPLEMENTED;
}

static int pack_double(grib_accessor* a, const double*v, size_t *len){
  grib_context_log(a->parent->h->context,GRIB_LOG_ERROR, " Should not pack %s  as double", a->name);
  return GRIB_NOT_IMPLEMENTED;
}


static int  unpack_long   (grib_accessor* a, long*  v, size_t *len){

  char val[1024];
  size_t l = sizeof(val);
  char  *last = NULL;
  grib_unpack_string (a , val, &l);

  *v = strtol(val,&last,10);

  if(*last == 0)
  {
    grib_context_log(a->parent->h->context,GRIB_LOG_DEBUG, " Casting string %s to long", a->name);
    return GRIB_SUCCESS;
  }

  return GRIB_INVALID_TYPE;
}

static int unpack_double (grib_accessor* a, double*v, size_t *len){
  char val[1024];
  size_t l = sizeof(val);
  char  *last = NULL;
  grib_unpack_string (a , val, &l);

  *v = strtod(val,&last);

  if(*last == 0)
  {
    grib_context_log(a->parent->h->context,GRIB_LOG_DEBUG, " Casting string %s to long", a->name);
    return GRIB_SUCCESS;
  }

  return GRIB_NOT_IMPLEMENTED;
}


static int compare(grib_accessor* a,grib_accessor* b) {
  int retval=0;
  char *aval=0;
  char *bval=0;

  size_t alen = (size_t)grib_value_count(a);
  size_t blen = (size_t)grib_value_count(b);

  if (alen != blen) return GRIB_COUNT_MISMATCH;

  aval=grib_context_malloc(a->parent->h->context,alen*sizeof(char));
  bval=grib_context_malloc(b->parent->h->context,blen*sizeof(char));

  grib_unpack_string(a,aval,&alen);
  grib_unpack_string(b,bval,&blen);

  retval = GRIB_SUCCESS;
  if (strcmp(aval,bval)) retval = GRIB_STRING_VALUE_MISMATCH;

  grib_context_free(a->parent->h->context,aval);
  grib_context_free(b->parent->h->context,bval);

  return retval;
}

