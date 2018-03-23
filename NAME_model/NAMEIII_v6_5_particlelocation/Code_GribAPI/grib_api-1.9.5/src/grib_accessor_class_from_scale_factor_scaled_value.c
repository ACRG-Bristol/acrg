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
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_double
   IMPLEMENTS = unpack_double;is_missing
   IMPLEMENTS = init
   MEMBERS=const char*    scaleFactor
   MEMBERS=const char*    scaledValue
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int is_missing(grib_accessor*);
static int unpack_double(grib_accessor*, double* val,size_t *len);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_from_scale_factor_scaled_value {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in double */
/* Members defined in from_scale_factor_scaled_value */
	const char*    scaleFactor;
	const char*    scaledValue;
} grib_accessor_from_scale_factor_scaled_value;

extern grib_accessor_class* grib_accessor_class_double;

static grib_accessor_class _grib_accessor_class_from_scale_factor_scaled_value = {
    &grib_accessor_class_double,                      /* super                     */
    "from_scale_factor_scaled_value",                      /* name                      */
    sizeof(grib_accessor_from_scale_factor_scaled_value),  /* size                      */
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
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    &is_missing,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    0,                /* grib_unpack procedures long    */
    0,                /* grib_pack procedures double    */
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


grib_accessor_class* grib_accessor_class_from_scale_factor_scaled_value = &_grib_accessor_class_from_scale_factor_scaled_value;


static void init_class(grib_accessor_class* c)
{
	c->dump	=	(*(c->super))->dump;
	c->next_offset	=	(*(c->super))->next_offset;
	c->value_count	=	(*(c->super))->value_count;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->unpack_long	=	(*(c->super))->unpack_long;
	c->pack_double	=	(*(c->super))->pack_double;
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
  grib_accessor_from_scale_factor_scaled_value* self = (grib_accessor_from_scale_factor_scaled_value*)a;
  int n = 0;

  self->scaleFactor = grib_arguments_get_name(a->parent->h,c,n++);
  self->scaledValue = grib_arguments_get_name(a->parent->h,c,n++);

  a->flags |= GRIB_ACCESSOR_FLAG_READ_ONLY;

}

static int    unpack_double   (grib_accessor* a, double* val, size_t *len) {
  grib_accessor_from_scale_factor_scaled_value* self = (grib_accessor_from_scale_factor_scaled_value*)a;
  int ret = 0;
  long scaleFactor=0;
  long scaledValue=0;
 
  if((ret = grib_get_long_internal(a->parent->h, self->scaleFactor,&scaleFactor)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->scaledValue,&scaledValue)) != GRIB_SUCCESS)
    return ret;

  *val=scaledValue;
  while (scaleFactor <0) {*val*=10;scaleFactor++;}
  while (scaleFactor >0) {*val/=10;scaleFactor--;}

  if (ret == GRIB_SUCCESS) *len = 1;

  return ret;
}


static int is_missing(grib_accessor* a){
  grib_accessor_from_scale_factor_scaled_value* self = (grib_accessor_from_scale_factor_scaled_value*)a;
  int ret = 0;
  long scaleFactor=0;
  long scaledValue=0;
 
  if((ret = grib_get_long_internal(a->parent->h, self->scaleFactor,&scaleFactor))
      != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->scaledValue,&scaledValue))
      != GRIB_SUCCESS)
    return ret;

  return ((scaleFactor == GRIB_MISSING_LONG) | (scaledValue == GRIB_MISSING_LONG));
}

