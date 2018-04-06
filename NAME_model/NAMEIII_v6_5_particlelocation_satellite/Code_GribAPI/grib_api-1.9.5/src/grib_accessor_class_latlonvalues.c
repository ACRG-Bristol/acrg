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
   IMPLEMENTS = unpack_double;
   IMPLEMENTS = value_count
   IMPLEMENTS = init
   MEMBERS =const char* values
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int unpack_double(grib_accessor*, double* val,size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_latlonvalues {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in double */
/* Members defined in latlonvalues */
	const char* values;
} grib_accessor_latlonvalues;

extern grib_accessor_class* grib_accessor_class_double;

static grib_accessor_class _grib_accessor_class_latlonvalues = {
    &grib_accessor_class_double,                      /* super                     */
    "latlonvalues",                      /* name                      */
    sizeof(grib_accessor_latlonvalues),  /* size                      */
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


grib_accessor_class* grib_accessor_class_latlonvalues = &_grib_accessor_class_latlonvalues;


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
  grib_accessor_latlonvalues* self = (grib_accessor_latlonvalues*)a;
  int n = 0;

  self->values = grib_arguments_get_name(a->parent->h,c,n++);

  a->flags |= GRIB_ACCESSOR_FLAG_READ_ONLY;

}

static int    unpack_double   (grib_accessor* a, double* val, size_t *len)
{
  grib_context* c=a->parent->h->context;
  int ret = 0;
  double* v=val;
  double lat,lon,value;
  size_t size=0;
  grib_iterator* iter=grib_iterator_new(a->parent->h,0,&ret);

  size=value_count(a);

  if (*len<size) return GRIB_ARRAY_TOO_SMALL;

  if (ret!=GRIB_SUCCESS) {
    if (iter) grib_iterator_delete(iter);
    grib_context_log(c,GRIB_LOG_ERROR,"unable to create iterator");
    return ret;
  }

  while(grib_iterator_next(iter,&lat,&lon,&value)) {
     *(v++)=lat;*(v++)=lon;*(v++)=value;
  }

  grib_iterator_delete(iter);

  *len=size;

  return ret;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_latlonvalues* self = (grib_accessor_latlonvalues*)a;
  grib_handle* h=a->parent->h;
  int ret;
  size_t size;
  if ((ret=grib_get_size(h,self->values,&size))!=GRIB_SUCCESS) {
    grib_context_log(h->context,GRIB_LOG_ERROR,"unable to get size of %s",self->values);
    return ret;
  }

  return (long)(3*size);

}


