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


#include <assert.h>
#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_gen
   IMPLEMENTS = unpack_long;get_native_type
   IMPLEMENTS = pack_long
   IMPLEMENTS = init
   MEMBERS=const char*    argument
   MEMBERS=long    start
   MEMBERS=long    len
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int  get_native_type(grib_accessor*);
static int pack_long(grib_accessor*, const long* val,size_t *len);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_bits {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in bits */
	const char*    argument;
	long    start;
	long    len;
} grib_accessor_bits;

extern grib_accessor_class* grib_accessor_class_gen;

static grib_accessor_class _grib_accessor_class_bits = {
    &grib_accessor_class_gen,                      /* super                     */
    "bits",                      /* name                      */
    sizeof(grib_accessor_bits),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    0,                       /* describes himself         */
    0,                /* get length of section     */
    0,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    &get_native_type,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    &pack_long,                  /* grib_pack procedures long      */
    &unpack_long,                /* grib_unpack procedures long    */
    0,                /* grib_pack procedures double    */
    0,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    0,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    0,            /* pack_expression */
    0,              /* notify_change   */
    0,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_bits = &_grib_accessor_class_bits;


static void init_class(grib_accessor_class* c)
{
	c->dump	=	(*(c->super))->dump;
	c->next_offset	=	(*(c->super))->next_offset;
	c->value_count	=	(*(c->super))->value_count;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_double	=	(*(c->super))->pack_double;
	c->unpack_double	=	(*(c->super))->unpack_double;
	c->pack_string	=	(*(c->super))->pack_string;
	c->unpack_string	=	(*(c->super))->unpack_string;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
	c->update_size	=	(*(c->super))->update_size;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->compare	=	(*(c->super))->compare;
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long l, grib_arguments* c)
{
  grib_accessor_bits* self = (grib_accessor_bits*)a;
  int n = 0;

  self->argument = grib_arguments_get_name(a->parent->h,c,n++);
  self->start    = grib_arguments_get_long(a->parent->h,c,n++);
  self->len    = grib_arguments_get_long(a->parent->h,c,n++);

  assert(self->len <= sizeof(long)*8);

  a->length=0;

}

static int    unpack_long   (grib_accessor* a, long* val, size_t *len)
{
  grib_accessor_bits* self = (grib_accessor_bits*)a;
  grib_accessor* x=NULL;
  unsigned char* p=NULL;
  grib_handle* h=a->parent->h;
  long start,length;
  int ret=0;

  if(*len < 1) return GRIB_WRONG_ARRAY_SIZE;

  start=self->start;
  length=self->len;
  
  x=grib_find_accessor(a->parent->h,self->argument);
  if (!x) return GRIB_NOT_FOUND;

  p  = h->buffer->data + grib_byte_offset(x);
  *val=grib_decode_unsigned_long(p,&start,length);

  *len=1;

  return ret;
}


static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_bits* self = (grib_accessor_bits*)a;
  grib_accessor* x=NULL;
  grib_handle* h=a->parent->h;
  unsigned char* p=NULL;
  long start,length;

  if(*len != 1) return GRIB_WRONG_ARRAY_SIZE;

  start  = self->start;
  length = self->len;
  
  x=grib_find_accessor(a->parent->h,self->argument);
  if (!x) return GRIB_NOT_FOUND;

  p=h->buffer->data + grib_byte_offset(x);
  return grib_encode_unsigned_longb(p,*val,&start,length);

}

static int  get_native_type(grib_accessor* a){
  return GRIB_TYPE_BYTES;
}




