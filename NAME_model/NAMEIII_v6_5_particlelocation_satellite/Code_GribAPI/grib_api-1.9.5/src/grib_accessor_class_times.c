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
   SUPER      = grib_accessor_class_unsigned
   IMPLEMENTS = unpack_long;pack_long; value_count
   IMPLEMENTS = init
   MEMBERS= const char*      value
   MEMBERS= const char*      factor
   MEMBERS= const char*      divisor
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
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_times {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in unsigned */
	long nbytes;
	grib_arguments* arg;
/* Members defined in times */
	const char*      value;
	const char*      factor;
	const char*      divisor;
} grib_accessor_times;

extern grib_accessor_class* grib_accessor_class_unsigned;

static grib_accessor_class _grib_accessor_class_times = {
    &grib_accessor_class_unsigned,                      /* super                     */
    "times",                      /* name                      */
    sizeof(grib_accessor_times),  /* size                      */
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


grib_accessor_class* grib_accessor_class_times = &_grib_accessor_class_times;


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
  grib_accessor_times* self = (grib_accessor_times*)a;
  int n = 0;

  self->value = grib_arguments_get_name(a->parent->h,c,n++);
  self->factor = grib_arguments_get_name(a->parent->h,c,n++);
  self->divisor = grib_arguments_get_name(a->parent->h,c,n++);
  a->length=0;
}

static int    unpack_long   (grib_accessor* a, long* val, size_t *len)
{
  grib_accessor_times* self = (grib_accessor_times*)a;
  int ret = 0;
  long factor=0;
  long divisor=1;
  int *err=&ret;
  long value = 0;

  if(*len < 1)
    return GRIB_ARRAY_TOO_SMALL;

  if (grib_is_missing(a->parent->h,self->value,err)!=0) {
    *val=GRIB_MISSING_LONG;
    return GRIB_SUCCESS;
  }
  if(ret ) return ret;

  ret = grib_get_long_internal(a->parent->h, self->factor,&factor);
  if(ret ) return ret;
  if (self->divisor)
  	ret = grib_get_long_internal(a->parent->h, self->divisor,&divisor);
  if(ret ) return ret;
  ret = grib_get_long_internal(a->parent->h, self->value,&value);
  if(ret ) return ret;
  /* printf("factor=%ld divisor=%ld value=%ld\n",factor,divisor,value); */

  *val =( (double)value * (double)factor) / (double)divisor;

  *len = 1;

  return ret;
}


static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_times* self = (grib_accessor_times*)a;
  int ret = 0;
  long value = 0;
  long factor,v,divisor=1;

  if (*val==GRIB_MISSING_LONG) 
    return grib_set_missing(a->parent->h,self->value);
  
  ret = grib_get_long_internal(a->parent->h, self->factor,&factor);
  if(ret ) return ret;
  if (self->divisor)
  	ret = grib_get_long_internal(a->parent->h, self->divisor,&divisor);
  if(ret ) return ret;

  /*Assert((*val%self->factor)==0);*/
  v=*val*divisor;
  if ((v%factor)==0) {
      value = v/factor;
  } else {
    value = v > 0 ? ((double)v)/factor+0.5 :
                       ((double)v)/factor-0.5;
    /* grib_context_log(a->parent->h->context,GRIB_LOG_WARNING,"%s/%ld = %ld/%ld = %ld. Rounding to convert key.",a->name,self->factor,*val,self->factor,value); */
  } 

  ret = grib_set_long_internal(a->parent->h, self->value,value);
  if(ret ) return ret;

  *len = 1;

  return ret;
}

static long value_count(grib_accessor* a) { return 1;}

