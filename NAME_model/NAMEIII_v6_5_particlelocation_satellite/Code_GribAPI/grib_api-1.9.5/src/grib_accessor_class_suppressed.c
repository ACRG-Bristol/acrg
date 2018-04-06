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
   IMPLEMENTS = value_count;get_native_type
   IMPLEMENTS = pack_string;unpack_string
   IMPLEMENTS = pack_long;unpack_long
   IMPLEMENTS = pack_double;unpack_double
   IMPLEMENTS = init
   MEMBERS= grib_arguments* args  
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int  get_native_type(grib_accessor*);
static int pack_double(grib_accessor*, const double* val,size_t *len);
static int pack_long(grib_accessor*, const long* val,size_t *len);
static int pack_string(grib_accessor*, const char*, size_t *len);
static int unpack_double(grib_accessor*, double* val,size_t *len);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_suppressed {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in suppressed */
	grib_arguments* args;
} grib_accessor_suppressed;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_suppressed = {
    &grib_accessor_class_long,                      /* super                     */
    "suppressed",                      /* name                      */
    sizeof(grib_accessor_suppressed),  /* size                      */
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
    &get_native_type,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    &pack_long,                  /* grib_pack procedures long      */
    &unpack_long,                /* grib_unpack procedures long    */
    &pack_double,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
    &pack_string,                /* grib_pack procedures string    */
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


grib_accessor_class* grib_accessor_class_suppressed = &_grib_accessor_class_suppressed;


static void init_class(grib_accessor_class* c)
{
	c->dump	=	(*(c->super))->dump;
	c->next_offset	=	(*(c->super))->next_offset;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
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
	grib_accessor_suppressed* self = (grib_accessor_suppressed*)a; 
	self->args = c;
  a->flags |= GRIB_ACCESSOR_FLAG_READ_ONLY;
  a->flags |= GRIB_ACCESSOR_FLAG_FUNCTION;
  a->length=0;
}

static void log_message(grib_accessor* a) {
  grib_accessor_suppressed* self = (grib_accessor_suppressed*)a;
  int i=0;
  
  grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
                   "key %s is unvailable in this version.",a->name);
  grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
                   "Please use the following keys:");
  while (grib_arguments_get_name(a->parent->h,self->args,i)) {
    grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
                     "\t- %s",
                     grib_arguments_get_name(a->parent->h,self->args,i));
    i++;
  }
}

static int pack_string(grib_accessor* a, const char* val, size_t *len){

  log_message(a);
	return GRIB_NOT_IMPLEMENTED;

}

static int pack_long(grib_accessor* a, const long* val, size_t *len){

  log_message(a);
  return GRIB_NOT_IMPLEMENTED;

}

static int pack_double(grib_accessor* a, const double* val, size_t *len){

  log_message(a);
  return GRIB_NOT_IMPLEMENTED;

}


static int    unpack_string(grib_accessor* a, char* val, size_t *len)
{   
  /*log_message(a);*/
	return GRIB_NOT_IMPLEMENTED;
}

static int unpack_long(grib_accessor* a, long* val, size_t *len)
{   
  log_message(a);
  return GRIB_NOT_IMPLEMENTED;
}

static int unpack_double(grib_accessor* a, double* val, size_t *len)
{
  log_message(a);
  return GRIB_NOT_IMPLEMENTED;
}


static long value_count(grib_accessor* a) { return 1; }

static int  get_native_type(grib_accessor* a){
  return GRIB_TYPE_STRING;
}

