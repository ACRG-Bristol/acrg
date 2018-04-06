/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*************************************************
 * Enrico Fucile
 ***********************************************/

#include "grib_api_internal.h"
#include <math.h>
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_long
   IMPLEMENTS = unpack_long;pack_long
   IMPLEMENTS = init
   MEMBERS=const char*                  N
   MEMBERS=const char*                  Ni
   MEMBERS=const char*                  di
   MEMBERS=const char*                  latfirst
   MEMBERS=const char*                  lonfirst
   MEMBERS=const char*                  latlast
   MEMBERS=const char*                  lonlast
   MEMBERS=const char*                  basic_angle
   MEMBERS=const char*                  subdivision
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
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_global_gaussian {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in global_gaussian */
	const char*                  N;
	const char*                  Ni;
	const char*                  di;
	const char*                  latfirst;
	const char*                  lonfirst;
	const char*                  latlast;
	const char*                  lonlast;
	const char*                  basic_angle;
	const char*                  subdivision;
} grib_accessor_global_gaussian;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_global_gaussian = {
    &grib_accessor_class_long,                      /* super                     */
    "global_gaussian",                      /* name                      */
    sizeof(grib_accessor_global_gaussian),  /* size                      */
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


grib_accessor_class* grib_accessor_class_global_gaussian = &_grib_accessor_class_global_gaussian;


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
  grib_accessor_global_gaussian* self = (grib_accessor_global_gaussian*)a;
  int n = 0;

  self->N        = grib_arguments_get_name(a->parent->h,c,n++);
  self->Ni        = grib_arguments_get_name(a->parent->h,c,n++);
  self->di         = grib_arguments_get_name(a->parent->h,c,n++);
  self->latfirst      = grib_arguments_get_name(a->parent->h,c,n++);
  self->lonfirst     = grib_arguments_get_name(a->parent->h,c,n++);
  self->latlast       = grib_arguments_get_name(a->parent->h,c,n++);
  self->lonlast      = grib_arguments_get_name(a->parent->h,c,n++);
  self->basic_angle         = grib_arguments_get_name(a->parent->h,c,n++);
  self->subdivision         = grib_arguments_get_name(a->parent->h,c,n++);
}


static int unpack_long(grib_accessor* a, long* val, size_t *len)
{
  grib_accessor_global_gaussian* self = (grib_accessor_global_gaussian*)a;
  int ret = GRIB_SUCCESS;
  long latfirst,latlast,lonfirst,lonlast,basic_angle,subdivision,N,Ni;
  double dlatfirst,dlatlast,dlonfirst,dlonlast,d;
  double* lats;
  long factor;
  grib_context* c=a->parent->h->context;

  if (self->basic_angle && self->subdivision) {

	  factor=1000000;
	  if((ret = grib_get_long_internal(a->parent->h, self->basic_angle,&basic_angle)) != GRIB_SUCCESS)
		return ret;

	  if((ret = grib_get_long_internal(a->parent->h, self->subdivision,&subdivision)) != GRIB_SUCCESS)
		return ret;

	  if (basic_angle !=0 || subdivision !=0 || subdivision == GRIB_MISSING_LONG) {
		*val=0;
		return ret;
	  }
  } else {
  	factor=1000;
  }

  if((ret = grib_get_long_internal(a->parent->h, self->N,&N)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->Ni,&Ni)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->latfirst,&latfirst)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->lonfirst,&lonfirst)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->latlast,&latlast)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->lonlast,&lonlast)) != GRIB_SUCCESS)
    return ret;

  dlatfirst=((double)latfirst)/factor;
  dlatlast=((double)latlast)/factor;
  dlonfirst=((double)lonfirst)/factor;
  dlonlast=((double)lonlast)/factor;

  lats=(double*)grib_context_malloc(c,sizeof(double)*N*2);
  if (!lats) {
  	grib_context_log(c,GRIB_LOG_FATAL,
		"global_gaussian: unable to allocate %d bytes",sizeof(double)*N*2);
  }
  if((ret = grib_get_gaussian_latitudes(N, lats)) != GRIB_SUCCESS)
      return ret;

  if (Ni == GRIB_MISSING_LONG ) Ni=N*4;
  d=fabs(lats[0]-lats[1]);
  if ( (fabs(dlatfirst-lats[0]) >= d ) ||
       (fabs(dlatlast+lats[0]) >= d )  ||
        dlonfirst != 0                 ||
        fabs(dlonlast  - (360.0-360.0/Ni)) > 360.0/Ni
     ) 
 {
	/* not global */
	*val=0;
  } else {
  	/* global */
	*val=1;
  }

  grib_context_free(c,lats);

  return ret;
}

static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_global_gaussian* self = (grib_accessor_global_gaussian*)a;
  int ret=GRIB_SUCCESS;
  long latfirst,latlast,lonfirst,lonlast,di,diold,basic_angle=0,subdivision=0,N,Ni;
  long factor;
  double* lats;
  grib_context* c=a->parent->h->context;

  if (*val == 0) return ret;

  if (self->basic_angle)  {
	  factor=1000000;
	  if((ret = grib_set_long_internal(a->parent->h, self->basic_angle,basic_angle)) != GRIB_SUCCESS)
		return ret;

	  if((ret = grib_set_long_internal(a->parent->h, self->subdivision,subdivision)) != GRIB_SUCCESS)
		return ret;
  } else factor=1000;

  if((ret = grib_get_long_internal(a->parent->h, self->N,&N)) != GRIB_SUCCESS)
    return ret;
  if (N==0) return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->Ni,&Ni)) != GRIB_SUCCESS)
    return ret;
  if (Ni == GRIB_MISSING_LONG ) Ni=N*4;
  if (Ni==0) return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->di,&diold)) != GRIB_SUCCESS)
    return ret;

  lats=(double*)grib_context_malloc(c,sizeof(double)*N*2);
  if (!lats) {
  	grib_context_log(c,GRIB_LOG_FATAL,
		"global_gaussian: unable to allocate %d bytes",sizeof(double)*N*2);
  }
  if((ret = grib_get_gaussian_latitudes(N, lats)) != GRIB_SUCCESS)
      return ret;

  latfirst=lats[0]*factor;
  latlast=-latfirst;
  lonfirst=0;
  di=(360*factor)/Ni;
  lonlast=(360*factor)-di;

  grib_context_free(c,lats);

  if((ret = grib_set_long_internal(a->parent->h, self->latfirst,latfirst)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_set_long_internal(a->parent->h, self->lonfirst,lonfirst)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_set_long_internal(a->parent->h, self->latlast,latlast)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_set_long_internal(a->parent->h, self->lonlast,lonlast)) != GRIB_SUCCESS)
    return ret;

  if (diold != GRIB_MISSING_LONG)
  	if((ret = grib_set_long_internal(a->parent->h, self->di,di)) != GRIB_SUCCESS)
		return ret;

  return GRIB_SUCCESS;
}



