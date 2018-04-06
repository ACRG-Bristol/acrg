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
   MEMBERS=const char* verification_yearmonth
   MEMBERS=const char* base_date
   MEMBERS=const char* day
   MEMBERS=const char* hour
   MEMBERS=const char* fcmonth
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int pack_long(grib_accessor*, const long* val,size_t *len);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_g1forecastmonth {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in g1forecastmonth */
	const char* verification_yearmonth;
	const char* base_date;
	const char* day;
	const char* hour;
	const char* fcmonth;
} grib_accessor_g1forecastmonth;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_g1forecastmonth = {
    &grib_accessor_class_long,                      /* super                     */
    "g1forecastmonth",                      /* name                      */
    sizeof(grib_accessor_g1forecastmonth),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    &dump,                       /* describes himself         */
    0,                /* get length of section     */
    0,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    0,            /* get native type               */
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


grib_accessor_class* grib_accessor_class_g1forecastmonth = &_grib_accessor_class_g1forecastmonth;


static void init_class(grib_accessor_class* c)
{
	c->next_offset	=	(*(c->super))->next_offset;
	c->value_count	=	(*(c->super))->value_count;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
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
  grib_accessor_g1forecastmonth* self = (grib_accessor_g1forecastmonth*)a;
  int n = 0;

  self->verification_yearmonth   = grib_arguments_get_name(a->parent->h,c,n++);
  self->base_date                = grib_arguments_get_name(a->parent->h,c,n++);
  self->day                      = grib_arguments_get_name(a->parent->h,c,n++);
  self->hour                     = grib_arguments_get_name(a->parent->h,c,n++);
  self->fcmonth                  = grib_arguments_get_name(a->parent->h,c,n++);
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
  grib_dump_long(dumper,a,NULL);
}

static int unpack_long(grib_accessor* a, long* val, size_t *len)
{
  int ret=0;
  grib_accessor_g1forecastmonth* self = (grib_accessor_g1forecastmonth*)a;

  long verification_yearmonth = 0;
  long base_yearmonth         = 0;
  long base_date              = 0;
  long day                    = 0;
  long hour                   = 0;

  long vyear  = 0;
  long vmonth = 0;
  long byear  = 0;
  long bmonth = 0;

  long fcmonth = 0;
  long gribForecastMonth = 0;

  if ((ret=grib_get_long_internal(a->parent->h,
       self->verification_yearmonth,&verification_yearmonth))!=GRIB_SUCCESS)
    return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->base_date,&base_date))!=GRIB_SUCCESS)
    return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->day,&day))!=GRIB_SUCCESS)
    return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->hour,&hour))!=GRIB_SUCCESS)
    return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->fcmonth,&gribForecastMonth))!=GRIB_SUCCESS)
    return ret;

  base_yearmonth = base_date / 100;

  vyear  = verification_yearmonth / 100;
  vmonth = verification_yearmonth % 100;
  byear  = base_yearmonth / 100;
  bmonth = base_yearmonth % 100;

  fcmonth = (vyear - byear) * 12 + (vmonth - bmonth);
  if(day == 1 && hour == 0)
    fcmonth++;

  if(gribForecastMonth != 0 && gribForecastMonth!=fcmonth) {
	  grib_context_log(a->parent->h->context,GRIB_LOG_FATAL,"%s=%ld (%s-%s)=%ld",self->fcmonth,
					   gribForecastMonth,self->base_date,self->verification_yearmonth,fcmonth);
	  Assert(gribForecastMonth == fcmonth);
  }

  *val = fcmonth;

  return GRIB_SUCCESS;
}

/* TODO: Check for a valid date */

static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_g1forecastmonth* self = (grib_accessor_g1forecastmonth*)a;
  return grib_set_long_internal(a->parent->h,self->fcmonth, *val);

}
