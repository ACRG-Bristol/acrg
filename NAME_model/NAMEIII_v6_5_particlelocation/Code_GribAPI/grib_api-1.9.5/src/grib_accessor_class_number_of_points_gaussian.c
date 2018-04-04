/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/**************************************
 * Enrico Fucile
 ************************************/

#include "grib_api_internal.h"


/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_long
   IMPLEMENTS = unpack_long
   IMPLEMENTS = init
   MEMBERS = const char* ni
   MEMBERS = const char* nj
   MEMBERS = const char* plpresent
   MEMBERS = const char* pl
   MEMBERS = const char* order
   MEMBERS = const char* lat_first
   MEMBERS = const char* lon_first
   MEMBERS = const char* lat_last
   MEMBERS = const char* lon_last
   END_CLASS_DEF
 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int unpack_long(grib_accessor*, long* val,size_t *len);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_number_of_points_gaussian {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in number_of_points_gaussian */
	const char* ni;
	const char* nj;
	const char* plpresent;
	const char* pl;
	const char* order;
	const char* lat_first;
	const char* lon_first;
	const char* lat_last;
	const char* lon_last;
} grib_accessor_number_of_points_gaussian;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_number_of_points_gaussian = {
    &grib_accessor_class_long,                      /* super                     */
    "number_of_points_gaussian",                      /* name                      */
    sizeof(grib_accessor_number_of_points_gaussian),  /* size                      */
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
    0,                  /* grib_pack procedures long      */
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


grib_accessor_class* grib_accessor_class_number_of_points_gaussian = &_grib_accessor_class_number_of_points_gaussian;


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
	c->pack_long	=	(*(c->super))->pack_long;
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

#define EFDEBUG 0

static void init(grib_accessor* a,const long l, grib_arguments* c)
{
  int n=0;
  grib_accessor_number_of_points_gaussian* self = (grib_accessor_number_of_points_gaussian*)a;
  self->ni = grib_arguments_get_name(a->parent->h,c,n++);
  self->nj = grib_arguments_get_name(a->parent->h,c,n++);
  self->plpresent = grib_arguments_get_name(a->parent->h,c,n++);
  self->pl = grib_arguments_get_name(a->parent->h,c,n++);
  self->order = grib_arguments_get_name(a->parent->h,c,n++);
  self->lat_first = grib_arguments_get_name(a->parent->h,c,n++);
  self->lon_first = grib_arguments_get_name(a->parent->h,c,n++);
  self->lat_last = grib_arguments_get_name(a->parent->h,c,n++);
  self->lon_last = grib_arguments_get_name(a->parent->h,c,n++);
  a->flags  |= GRIB_ACCESSOR_FLAG_READ_ONLY;
  a->flags |= GRIB_ACCESSOR_FLAG_FUNCTION;
  a->length=0;
}

static int  unpack_long(grib_accessor* a, long* val, size_t *len){
  int ret=GRIB_SUCCESS;
  long ni=0,nj=0,plpresent=0,order=0;
  size_t plsize=0;
  double* lats={0,};
  double lat_first,lat_last,lon_first,lon_last;
  float d;
  long* pl=NULL;
  long* plsave=NULL;
  int j=0,i=0;
  long row_count;
  long ilon_first=0,ilon_last=0;
  double lon_first_row=0,lon_last_row=0;


  grib_accessor_number_of_points_gaussian* self = (grib_accessor_number_of_points_gaussian*)a;
  grib_context* c=a->parent->h->context;

  if((ret = grib_get_long_internal(a->parent->h, self->ni,&ni)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->nj,&nj)) != GRIB_SUCCESS)
    return ret;

  if((ret = grib_get_long_internal(a->parent->h, self->plpresent,&plpresent)) != GRIB_SUCCESS)
    return ret;

  if (nj == 0) return GRIB_GEOCALCULUS_PROBLEM;

  if (plpresent) {
    /*reduced*/
    if((ret = grib_get_long_internal(a->parent->h, self->order,&order)) != GRIB_SUCCESS)
      return ret;
    if((ret = grib_get_double_internal(a->parent->h, self->lat_first,&lat_first)) != GRIB_SUCCESS)
      return ret;
    if((ret = grib_get_double_internal(a->parent->h, self->lon_first,&lon_first)) != GRIB_SUCCESS)
      return ret;
    if((ret = grib_get_double_internal(a->parent->h, self->lat_last,&lat_last)) != GRIB_SUCCESS)
      return ret;
    if((ret = grib_get_double_internal(a->parent->h, self->lon_last,&lon_last)) != GRIB_SUCCESS)
      return ret;

    lats=(double*)grib_context_malloc(a->parent->h->context,sizeof(double)*order*2);
    if((ret = grib_get_gaussian_latitudes(order, lats)) != GRIB_SUCCESS)
      return ret;

    if((ret = grib_get_size(a->parent->h,self->pl,&plsize)) != GRIB_SUCCESS)
      return ret;

    pl=(long*)grib_context_malloc_clear(c,sizeof(long)*plsize);
    plsave=pl;
    grib_get_long_array_internal(a->parent->h,self->pl,pl, &plsize);

    if (lon_last<0) lon_last+=360;
    if (lon_first<0) lon_last+=360;

    d=fabs(lats[0]-lats[1]);
    if ( (fabs(lat_first-lats[0]) >= d ) ||
         (fabs(lat_last+lats[0]) >= d )  ||
         lon_first != 0                 ||
         fabs(lon_last  - (360.0-90.0/order)) > 90.0/order
         ) {
      /*sub area*/
#if EFDEBUG
      printf("-------- subarea fabs(lat_first-lats[0])=%g d=%g\n",fabs(lat_first-lats[0]),d);
      printf("-------- subarea fabs(lat_last+lats[0])=%g d=%g\n",fabs(lat_last+lats[0]),d);
      printf("-------- subarea lon_last=%g order=%ld 360.0-90.0/order=%g\n",
      lon_last,order,360.0-90.0/order);
      printf("-------- subarea lon_first=%g fabs(lon_last  -( 360.0-90.0/order))=%g 90.0/order=%g\n",
      lon_first,fabs(lon_last  - (360.0-90.0/order)),90.0/order);
#endif
      *val=0;
      for (j=0;j<nj;j++) {
        row_count=0;
#if EFDEBUG
        printf("--  %d ",j);
#endif
		grib_get_reduced_row(pl[j],lon_first,lon_last,&row_count,&ilon_first,&ilon_last);
        lon_first_row=((ilon_first)*360.0)/pl[j];
        lon_last_row=((ilon_last)*360.0)/pl[j];
        *val+=row_count;
#if EFDEBUG
        printf("        ilon_first=%ld lon_first=%.10e ilon_last=%ld lon_last=%.10e count=%ld row_count=%ld\n",
            ilon_first,lon_first_row,ilon_last,lon_last_row,*val,row_count);
#endif
      }

    } else {
      *val=0;
      for (i=0;i<plsize;i++) *val+=pl[i];
    }
  } else {
    /*regular*/
    *val=ni*nj;
  }
#if EFDEBUG
  printf("DEBUG:     number_of_points_gaussian=%ld plpresent=%ld plsize=%d\n",*val,plpresent,plsize);
  for (i=0;i<plsize;i++) printf(" DEBUG: pl[%d]=%ld\n",i,pl[i]);
#endif
  if (lats) grib_context_free(c,lats);
  if (plsave) grib_context_free(c,plsave);
  return ret;
}


