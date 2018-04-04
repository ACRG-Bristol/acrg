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
   SUPER      = grib_accessor_class_double
   IMPLEMENTS = pack_double;unpack_string
   IMPLEMENTS = init;dump
   IMPLEMENTS = unpack_double
   MEMBERS=const char* laf
   MEMBERS=const char* lof
   MEMBERS=const char* lal
   MEMBERS=const char* lol
   MEMBERS=const char* div
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
static int unpack_string (grib_accessor*, char*, size_t *len);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_g1area {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in double */
/* Members defined in g1area */
	const char* laf;
	const char* lof;
	const char* lal;
	const char* lol;
	const char* div;
} grib_accessor_g1area;

extern grib_accessor_class* grib_accessor_class_double;

static grib_accessor_class _grib_accessor_class_g1area = {
    &grib_accessor_class_double,                      /* super                     */
    "g1area",                      /* name                      */
    sizeof(grib_accessor_g1area),  /* size                      */
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
    0,                  /* grib_pack procedures long      */
    0,                /* grib_unpack procedures long    */
    &pack_double,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    &unpack_string,              /* grib_unpack procedures string  */
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


grib_accessor_class* grib_accessor_class_g1area = &_grib_accessor_class_g1area;


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
	c->pack_long	=	(*(c->super))->pack_long;
	c->unpack_long	=	(*(c->super))->unpack_long;
	c->pack_string	=	(*(c->super))->pack_string;
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
	grib_accessor_g1area* self = (grib_accessor_g1area*)a; 
	int n = 0;

	self->laf = grib_arguments_get_name(a->parent->h,c,n++);
	self->lof = grib_arguments_get_name(a->parent->h,c,n++);
	self->lal = grib_arguments_get_name(a->parent->h,c,n++);
	self->lol = grib_arguments_get_name(a->parent->h,c,n++);

}


static int    pack_double   (grib_accessor* a, const double* val, size_t *len)
{
	grib_accessor_g1area* self = (grib_accessor_g1area*)a;
	int ret = 0;

	ret = grib_set_double_internal(a->parent->h, self->laf,val[0]);
	if(ret ) return ret;
	ret = grib_set_double_internal(a->parent->h, self->lof,val[1]);
  if(ret ) return ret;
	ret = grib_set_double_internal(a->parent->h, self->lal,val[2]);
  if(ret) return ret;

	ret = grib_set_double_internal(a->parent->h, self->lol,val[3]);
  if(ret ) return ret;

	if (ret == GRIB_SUCCESS) *len = 4;

	return ret;
}
static int    unpack_double   (grib_accessor* a,  double* val, size_t *len)
{
	grib_accessor_g1area* self = (grib_accessor_g1area*)a;
	int ret = 0;

	if(*len < 4){
		*len = 4;
		return GRIB_BUFFER_TOO_SMALL;
	}

	ret = grib_get_double_internal(a->parent->h, self->laf,val++);
  if(ret) return ret;

	ret = grib_get_double_internal(a->parent->h, self->lof,val++);
  if(ret) return ret;

	ret = grib_get_double_internal(a->parent->h, self->lal,val++);
  if(ret) return ret;

	ret = grib_get_double_internal(a->parent->h, self->lol,val);
  if(ret ) return ret;

	if (ret == GRIB_SUCCESS) *len = 4;

	return ret;
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
	grib_dump_string(dumper,a,NULL);
}

static int    unpack_string(grib_accessor* a, char* val, size_t *len)
{   
	grib_accessor_g1area* self = (grib_accessor_g1area*)a;
	int ret = 0;

	double laf,lof,lal,lol;

	ret = grib_get_double_internal(a->parent->h, self->laf,&laf);
  if(ret) return ret;
	ret = grib_get_double_internal(a->parent->h, self->lof,&lof);
  if(ret) return ret;
	ret = grib_get_double_internal(a->parent->h, self->lal,&lal);
  if(ret) return ret;
	ret = grib_get_double_internal(a->parent->h, self->lol,&lol);
  if(ret) return ret;
	if(*len < 60)
	{
		grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, " Buffer too smalle for %s (%d) ", a->name ,*len);
		len = 0;
		return GRIB_BUFFER_TOO_SMALL;
	}


	sprintf(val,"N:%3.5f W:%3.5f S:%3.5f E:%3.5f",((float)laf),((float)lof),((float)lal),((float)lol));

	len[0] = strlen(val);
	return GRIB_SUCCESS;
}

