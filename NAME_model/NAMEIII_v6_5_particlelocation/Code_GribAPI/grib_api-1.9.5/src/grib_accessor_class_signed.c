/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_long
   IMPLEMENTS = unpack_long;pack_long
   IMPLEMENTS = init;dump
   IMPLEMENTS = next_offset
   IMPLEMENTS = byte_count
   IMPLEMENTS = value_count
   IMPLEMENTS = byte_offset
   IMPLEMENTS = update_size; is_missing
   MEMBERS    = grib_arguments* arg
   MEMBERS    = int nbytes;
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int is_missing(grib_accessor*);
static int pack_long(grib_accessor*, const long* val,size_t *len);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static long byte_count(grib_accessor*);
static long byte_offset(grib_accessor*);
static long next_offset(grib_accessor*);
static long value_count(grib_accessor*);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static void update_size(grib_accessor*,size_t);

typedef struct grib_accessor_signed {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in signed */
	grib_arguments* arg;
	int nbytes;
} grib_accessor_signed;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_signed = {
    &grib_accessor_class_long,                      /* super                     */
    "signed",                      /* name                      */
    sizeof(grib_accessor_signed),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    &dump,                       /* describes himself         */
    &next_offset,                /* get length of section     */
    &value_count,                /* get number of values      */
    &byte_count,                 /* get number of bytes      */
    &byte_offset,                /* get offset to bytes           */
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    &is_missing,               /* grib_pack procedures long      */
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
    &update_size,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_signed = &_grib_accessor_class_signed;


static void init_class(grib_accessor_class* c)
{
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->pack_double	=	(*(c->super))->pack_double;
	c->unpack_double	=	(*(c->super))->unpack_double;
	c->pack_string	=	(*(c->super))->pack_string;
	c->unpack_string	=	(*(c->super))->unpack_string;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
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

static void init(grib_accessor* a, const long len , grib_arguments* arg )
{
  grib_accessor_signed* self = (grib_accessor_signed*)a;
  self->arg = arg;
  a->length = len * grib_value_count(a);
  self->nbytes = len;
  Assert(a->length>=0);

}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
  long rlen = grib_value_count(a);
  if(rlen == 1)
    grib_dump_long(dumper,a,NULL);
  else
    grib_dump_values(dumper,a);
}

static long ones[] = {
  0,
  -0x7f,
  -0x7fff,
  -0x7fffff,
  -0x7fffffff,
};


static int unpack_long(grib_accessor* a, long* val, size_t *len)
{

  grib_accessor_signed* self = (grib_accessor_signed*)a;
  unsigned long rlen = grib_value_count(a);
  unsigned long i = 0;

  long pos = a->offset;
  long missing = 0;


  if(*len < rlen)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, " wrong size for %s it contains %d values ", a->name , rlen);
    *len = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  if(a->flags & GRIB_ACCESSOR_FLAG_CAN_BE_MISSING)
  {
    Assert(self->nbytes <= 4);
    missing = ones[self->nbytes];
  }


  for(i=0; i< rlen;i++){
    val[i] = (long)grib_decode_signed_long(a->parent->h->buffer->data , pos, self->nbytes);
    if(missing)
      if(val[i] == missing)
        val[i] = GRIB_MISSING_LONG;
    pos += self->nbytes;
  }

  *len = rlen;
  return GRIB_SUCCESS;
}

static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_signed* self = (grib_accessor_signed*)a;
  int ret = 0;
  long off = 0;
  unsigned long rlen = grib_value_count(a);
  size_t buflen  = 0;
  unsigned char *buf = NULL;
  unsigned long i = 0;
  long missing = 0;


  if(*len < 1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , 1 );
    len[0] = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  if(a->flags & GRIB_ACCESSOR_FLAG_CAN_BE_MISSING)
      {
            Assert(self->nbytes <= 4);
                missing = ones[self->nbytes];
                }

  if (rlen == 1){
    long v = val[0];
    if(missing)
         if(v == GRIB_MISSING_LONG)
                     v = missing;


    off = a->offset;
    ret = grib_encode_signed_long(a->parent->h->buffer->data, v ,  off,  a->length);
    if (ret == GRIB_SUCCESS) len[0] = 1;
    if (*len > 1)  grib_context_log(a->parent->h->context, GRIB_LOG_WARNING, "grib_accessor_signed : Trying to pack %d values in a scalar %s, packing first value",  *len, a->name  );
    len[0] = 1;
    return ret;
  }

   /* TODO: We assume that there are no missing values if there are more that 1 value */



  buflen = *len*a->length;

  buf = grib_context_malloc(a->parent->h->context,buflen);

  for(i=0; i < *len;i++){
    grib_encode_signed_long(buf, val[i] ,  off,  a->length);
    off+=  a->length;
  }
  ret = grib_set_long_internal(a->parent->h,grib_arguments_get_name(a->parent->h,self->arg,0),*len);

  if(ret == GRIB_SUCCESS)
    grib_buffer_replace(a, buf, buflen,1,1);
  else
    *len = 0;

  grib_context_free(a->parent->h->context,buf);
  return ret;
}

static long byte_count(grib_accessor* a){
  return a->length;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_signed* self = (grib_accessor_signed*)a;
  long len = 0;
  int ret =0;
  if(!self->arg) return 1;
  ret = grib_get_long_internal(a->parent->h,grib_arguments_get_name(a->parent->h,self->arg,0),&len);
  if(ret == GRIB_SUCCESS)  return len;
  return 1;
}

static long byte_offset(grib_accessor* a){
  return a->offset;
}

static void update_size(grib_accessor* a,size_t s)
{
  a->length = s;
  Assert(a->length>=0);
}

static long next_offset(grib_accessor* a){
  return grib_byte_offset(a)+grib_byte_count(a);
}

static int is_missing(grib_accessor* a){
  int i=0;
  unsigned char ff=0xff;
  unsigned long offset=a->offset;

  if (a->length==0) {
    Assert(a->vvalue!=NULL);
	return a->vvalue->missing;
  }

  for (i=0;i<a->length;i++) {
     if (a->parent->h->buffer->data[offset] != ff) return 0;
     offset++;
  }

  return 1;
}
