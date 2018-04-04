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
   SUPER      = grib_accessor_class_data_shsimple_packing
   IMPLEMENTS = init
   IMPLEMENTS = unpack_double
   IMPLEMENTS = pack_double
   IMPLEMENTS = value_count
   MEMBERS=const char*  numberOfValues
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

typedef struct grib_accessor_data_g2shsimple_packing {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in data_shsimple_packing */
	const char*  coded_values;
	const char*  real_part;
	int dirty;
/* Members defined in data_g2shsimple_packing */
	const char*  numberOfValues;
} grib_accessor_data_g2shsimple_packing;

extern grib_accessor_class* grib_accessor_class_data_shsimple_packing;

static grib_accessor_class _grib_accessor_class_data_g2shsimple_packing = {
    &grib_accessor_class_data_shsimple_packing,                      /* super                     */
    "data_g2shsimple_packing",                      /* name                      */
    sizeof(grib_accessor_data_g2shsimple_packing),  /* size                      */
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
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_data_g2shsimple_packing = &_grib_accessor_class_data_g2shsimple_packing;


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
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long v, grib_arguments* args)
{
  grib_accessor_data_g2shsimple_packing *self =(grib_accessor_data_g2shsimple_packing*)a;

  self->numberOfValues  = grib_arguments_get_name(a->parent->h,args,2);
  a->flags |= GRIB_ACCESSOR_FLAG_DATA;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_data_g2shsimple_packing *self =(grib_accessor_data_g2shsimple_packing*)a;
  long len = 0;
  int ret =   grib_get_long(a->parent->h,self->numberOfValues,&len);

  if( ret == GRIB_SUCCESS)
      return len;
   return 0;
}


static int  unpack_double(grib_accessor* a, double* val, size_t *len)
{
  grib_accessor_data_g2shsimple_packing* self =  (grib_accessor_data_g2shsimple_packing*)a;
  int err =  GRIB_SUCCESS;

  size_t n_vals = 0;

  if((err = grib_get_size(a->parent->h,self->coded_values,&n_vals)) != GRIB_SUCCESS)
    return err;

  self->dirty=0;

    /* n_vals = coded_n_vals+1; */

  if(*len < n_vals)
  {
    *len = n_vals;
    return GRIB_ARRAY_TOO_SMALL;
  }

  if((err = grib_get_double_internal(a->parent->h,self->real_part,val)) != GRIB_SUCCESS)
    return err;

  val++;

  if((err = grib_get_double_array_internal(a->parent->h,self->coded_values,val,&n_vals)) != GRIB_SUCCESS)
    return err;

  *len =  n_vals;

  return err;
}

static int pack_double(grib_accessor* a, const double* val, size_t *len)
{
  grib_accessor_data_g2shsimple_packing* self =  (grib_accessor_data_g2shsimple_packing*)a;
  int err =  GRIB_SUCCESS;

  size_t coded_n_vals = *len-1;
  size_t n_vals = *len;

  if (*len ==0) return GRIB_NO_VALUES;

  self->dirty=1;

  if((err = grib_set_double_internal(a->parent->h,self->real_part,*val)) != GRIB_SUCCESS)
    return err;
  {
    /* Make sure we can decode it again */
    double ref = 1e-100;
    grib_get_double_internal(a->parent->h,self->real_part,&ref);
    Assert(ref == *val);
  }

  val++;

  if((err = grib_set_double_array_internal(a->parent->h,self->coded_values,val,coded_n_vals)) != GRIB_SUCCESS)
    return err;

  *len =  n_vals;

    if((err = grib_set_long_internal(a->parent->h,self->numberOfValues,(long)n_vals)) != GRIB_SUCCESS)
    return err;

  return err;
}

