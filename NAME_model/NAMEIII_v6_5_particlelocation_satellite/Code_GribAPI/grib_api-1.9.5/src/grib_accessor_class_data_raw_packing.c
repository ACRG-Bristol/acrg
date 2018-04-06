/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/
/*****************************
 *  Enrico Fucile
 ****************************/

#include "grib_api_internal.h"

#define PRE_PROCESSING_NONE       0
#define PRE_PROCESSING_DIFFERENCE 1
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_values
   IMPLEMENTS = init
   IMPLEMENTS = unpack_double
   IMPLEMENTS = unpack_double_element
   IMPLEMENTS = pack_double
   IMPLEMENTS = value_count
   MEMBERS=const char*   number_of_values
   MEMBERS=const char*   precision
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int pack_double(grib_accessor*, const double* val,size_t *len);
static int unpack_double(grib_accessor*, double* val,size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static int unpack_double_element(grib_accessor*,size_t i, double* val);

typedef struct grib_accessor_data_raw_packing {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in values */
	int  carg;
	const char* seclen;
	const char* offsetdata;
	const char* offsetsection;
	int dirty;
/* Members defined in data_raw_packing */
	const char*   number_of_values;
	const char*   precision;
} grib_accessor_data_raw_packing;

extern grib_accessor_class* grib_accessor_class_values;

static grib_accessor_class _grib_accessor_class_data_raw_packing = {
    &grib_accessor_class_values,                      /* super                     */
    "data_raw_packing",                      /* name                      */
    sizeof(grib_accessor_data_raw_packing),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    0,                       /* describes himself         */
    0,                /* get length of section     */
    &value_count,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    0,                /* grib_unpack procedures long    */
    &pack_double,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
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
    &unpack_double_element,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_data_raw_packing = &_grib_accessor_class_data_raw_packing;


static void init_class(grib_accessor_class* c)
{
	c->dump	=	(*(c->super))->dump;
	c->next_offset	=	(*(c->super))->next_offset;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->unpack_long	=	(*(c->super))->unpack_long;
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
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long v, grib_arguments* args)
{
  grib_accessor_data_raw_packing *self =(grib_accessor_data_raw_packing*)a;

  self->number_of_values      = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->precision       = grib_arguments_get_name(a->parent->h,args,self->carg++);
  a->flags |= GRIB_ACCESSOR_FLAG_DATA;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_data_raw_packing *self =(grib_accessor_data_raw_packing*)a;
  long n_vals= 0;
  if(grib_get_long_internal(a->parent->h,self->number_of_values,&n_vals) != GRIB_SUCCESS)
    return 0;

  return n_vals;
}

static int  unpack_double(grib_accessor* a, double* val, size_t *len)
{
  grib_accessor_data_raw_packing *self =(grib_accessor_data_raw_packing*)a;
  unsigned char* buf = NULL;
  int bytes = 0;
  size_t nvals = 0;
  long inlen = grib_byte_count(a);

  long precision = 0;

  int code = GRIB_SUCCESS;

  if((code = grib_get_long_internal(a->parent->h,self->precision,&precision))
      != GRIB_SUCCESS)
    return code;

  self->dirty=0;

  buf =  (unsigned char*)a->parent->h->buffer->data;
  buf += grib_byte_offset(a);

  switch(precision)
  {
    case 1:
      bytes = 4;
      break;

    case 2:
      bytes = 8;
      break;

    default:
      return GRIB_NOT_IMPLEMENTED;
      break;
  }

  nvals = inlen / bytes;

  if(*len < nvals)
    return GRIB_ARRAY_TOO_SMALL;

  code=grib_ieee_decode_array(a->parent->h->context,buf,nvals,bytes,val);

  *len = nvals;

  return code;
}

static int pack_double(grib_accessor* a, const double* val, size_t *len)
{
  grib_accessor_data_raw_packing *self =(grib_accessor_data_raw_packing*)a;

  int bytes = 0;
  unsigned char* buffer = NULL;

  long precision = 0;

  double*  values = (double*)val;
  size_t inlen = *len;

  int free_buffer = 0;
  int free_values = 0;

  int code = GRIB_SUCCESS;


  size_t bufsize = 0;

  if (*len ==0) return GRIB_NO_VALUES;

  if((code = grib_get_long_internal(a->parent->h,self->precision,&precision))
      != GRIB_SUCCESS)
    return code;

  self->dirty=1;

  switch(precision)
  {
    case 1:
      bytes = 4;
      break;

    case 2:
      bytes = 8;
      break;

    default:
      code = GRIB_NOT_IMPLEMENTED;
      goto clean_up;
      break;
  }

  bufsize = bytes*inlen;

  buffer = grib_context_malloc(a->parent->h->context, bufsize);

  if(!buffer)
  {
    code = GRIB_OUT_OF_MEMORY;
    goto clean_up;
  }

  code=grib_ieee_encode_array(a->parent->h->context,values,inlen,bytes,buffer);

clean_up:
  if(free_buffer) free(buffer);
  if(free_values) free(values);

  grib_buffer_replace(a, buffer, bufsize,1,1);

  grib_context_buffer_free(a->parent->h->context,buffer);

  code = grib_set_long(a->parent->h,self->number_of_values, inlen);
  if(code==GRIB_READ_ONLY) code=0;

  return code;

}

static int  unpack_double_element(grib_accessor* a, size_t idx, double* val) {
  int ret=0;
  grib_accessor_data_raw_packing *self =(grib_accessor_data_raw_packing*)a;
  unsigned char* buf = NULL;
  int bytes = 0;
  size_t nvals = 0;
  long inlen = grib_byte_count(a);
  long pos;

  long precision = 0;

  if((ret = grib_get_long_internal(a->parent->h,self->precision,&precision))
      != GRIB_SUCCESS)
    return ret;

  self->dirty=0;

  buf =  (unsigned char*)a->parent->h->buffer->data;
  buf += grib_byte_offset(a);

  switch(precision)
  {
    case 1:
      bytes = 4;
      break;

    case 2:
      bytes = 8;
      break;

    default:
      return GRIB_NOT_IMPLEMENTED;
      break;
  }

  pos=bytes*idx;
  
  Assert(pos<=inlen);
  
  nvals = 1;
  buf+=pos;
  
  ret=grib_ieee_decode_array(a->parent->h->context,buf,nvals,bytes,val);

  return ret;
}
