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
   SUPER      = grib_accessor_class_data_simple_packing
   IMPLEMENTS = init
   IMPLEMENTS = pack_double
   IMPLEMENTS = pack_bytes 
   IMPLEMENTS = value_count
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int pack_bytes(grib_accessor*,const unsigned char*, size_t *len);
static int pack_double(grib_accessor*, const double* val,size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_data_g2simple_packing {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in values */
	int  carg;
	const char* seclen;
	const char* offsetdata;
	const char* offsetsection;
	int dirty;
/* Members defined in data_simple_packing */
	const char*  units_factor;
	const char*  units_bias;
	const char*  changing_precision;
	const char*  number_of_values;
	const char*  bits_per_value;
	const char*  reference_value;
	const char*  binary_scale_factor;
	const char*  decimal_scale_factor;
/* Members defined in data_g2simple_packing */
} grib_accessor_data_g2simple_packing;

extern grib_accessor_class* grib_accessor_class_data_simple_packing;

static grib_accessor_class _grib_accessor_class_data_g2simple_packing = {
    &grib_accessor_class_data_simple_packing,                      /* super                     */
    "data_g2simple_packing",                      /* name                      */
    sizeof(grib_accessor_data_g2simple_packing),  /* size                      */
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
    0,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    0,              /* grib_unpack procedures string  */
    &pack_bytes,                 /* grib_pack procedures bytes     */
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


grib_accessor_class* grib_accessor_class_data_g2simple_packing = &_grib_accessor_class_data_g2simple_packing;


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
	c->unpack_double	=	(*(c->super))->unpack_double;
	c->pack_string	=	(*(c->super))->pack_string;
	c->unpack_string	=	(*(c->super))->unpack_string;
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

static void init(grib_accessor* a,const long v, grib_arguments* args)
{
  a->flags |= GRIB_ACCESSOR_FLAG_DATA;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_data_g2simple_packing *self =(grib_accessor_data_g2simple_packing*)a;
  long n_vals= 0;
  if(grib_get_long_internal(a->parent->h,self->number_of_values,&n_vals) != GRIB_SUCCESS)
    return 0;
  return n_vals;
}



static int pack_double(grib_accessor* a, const double* cval, size_t *len)
{
  grib_accessor_data_g2simple_packing* self =  (grib_accessor_data_g2simple_packing*)a;
  grib_accessor_class* super = *(a->cclass->super);
  size_t n_vals = *len;
  double reference_value = 0;
  long   binary_scale_factor = 0;
  long   bits_per_value = 0;
  long   decimal_scale_factor = 0;
  double decimal = 1;
  size_t buflen = 0;
  unsigned char*  buf = NULL; 
  unsigned char*  encoded = NULL; 
  double divisor = 1;
  long off = 0;
  int ret =0;
  double units_factor=1.0;
  double units_bias=0.0;
  double* val=(double*)cval;
  int i;

  if(*len == 0) {
    grib_buffer_replace(a, NULL, 0,1,1);
    return GRIB_SUCCESS;
  }

  if(ret == GRIB_SUCCESS)
    ret = grib_set_long_internal(a->parent->h,self->number_of_values, *len);

  if (ret!=GRIB_SUCCESS) return ret;

	 if(self->units_factor &&
	 (grib_get_double_internal(a->parent->h,self->units_factor,&units_factor)== GRIB_SUCCESS)) {
	  grib_set_double_internal(a->parent->h,self->units_factor,1.0);
	 }

	 if(self->units_bias &&
		   (grib_get_double_internal(a->parent->h,self->units_bias,&units_bias)== GRIB_SUCCESS)) {
		 grib_set_double_internal(a->parent->h,self->units_bias,0.0);
	}
		   
    if (units_factor != 1.0) {
        if (units_bias != 0.0)
            for (i=0;i<n_vals;i++) val[i]=val[i]*units_factor+units_bias;
        else
            for (i=0;i<n_vals;i++) val[i]*=units_factor;
    } else if (units_bias != 0.0)
        for (i=0;i<n_vals;i++) val[i]+=units_bias;

  ret = super->pack_double(a,val,len);
  switch (ret) {
    case GRIB_CONSTANT_FIELD:
      grib_buffer_replace(a, NULL, 0,1,1);
      return GRIB_SUCCESS;
    case GRIB_SUCCESS:
      break;
    default:
      grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,"unable to compute packing parameters\n");
      return ret;
  }

  if((ret = grib_get_double_internal(a->parent->h,self->reference_value, &reference_value))
      != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h,self->binary_scale_factor, &binary_scale_factor))
      != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h,self->bits_per_value,&bits_per_value)) !=
      GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h,self->decimal_scale_factor, &decimal_scale_factor))
     != GRIB_SUCCESS)
    return ret;

  decimal = grib_power(decimal_scale_factor,10) ;
  divisor = grib_power(-binary_scale_factor,2);

  buflen = (((bits_per_value*n_vals)+7)/8)*sizeof(unsigned char);
  buf = grib_context_buffer_malloc_clear(a->parent->h->context,buflen);
  encoded = buf;

  grib_encode_double_array(n_vals,val,bits_per_value,reference_value,decimal,divisor,encoded,&off);

  grib_context_log(a->parent->h->context, GRIB_LOG_DEBUG,
      "grib_accessor_data_g2simple_packing : pack_double : packing %s, %d values", a->name, n_vals);

  grib_buffer_replace(a, buf, buflen,1,1);

  grib_context_buffer_free(a->parent->h->context,buf);

  return ret;
}

static int pack_bytes(grib_accessor* a, const unsigned char* val, size_t *len)
{
  size_t length = *len;
  grib_buffer_replace(a, val, length,1,1);
  return GRIB_SUCCESS;
}

