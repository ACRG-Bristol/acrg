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
   SUPER      = grib_accessor_class_gen
   IMPLEMENTS = unpack_double;pack_double
   IMPLEMENTS = unpack_string;pack_string
   IMPLEMENTS = unpack_long;pack_long;destroy
   IMPLEMENTS = init;dump;value_count;get_native_type
   IMPLEMENTS = compare
   MEMBERS=double dval
   MEMBERS=char*  cval
   MEMBERS=int    type
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
static void destroy(grib_context*,grib_accessor*);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static int compare(grib_accessor*, grib_accessor*);

typedef struct grib_accessor_variable {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in variable */
	double dval;
	char*  cval;
	int    type;
} grib_accessor_variable;

extern grib_accessor_class* grib_accessor_class_gen;

static grib_accessor_class _grib_accessor_class_variable = {
    &grib_accessor_class_gen,                      /* super                     */
    "variable",                      /* name                      */
    sizeof(grib_accessor_variable),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    &destroy,                    /* free mem                       */
    &dump,                       /* describes himself         */
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
    &compare,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_variable = &_grib_accessor_class_variable;


static void init_class(grib_accessor_class* c)
{
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
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a, const long length , grib_arguments* args )
{
  grib_accessor_variable* self = (grib_accessor_variable*)a;
  grib_expression *expression = grib_arguments_get_expression(a->parent->h,args,0);
  const char* p = 0;
  size_t len = 1;
  long l;
  int ret=0;
  double d;
  char tmp[1024];

  a->length = 0;
  self->type = grib_expression_native_type(a->parent->h,expression);

  switch(self->type)
  {
    case GRIB_TYPE_DOUBLE:
      grib_expression_evaluate_double(a->parent->h,expression,&d);
      pack_double(a,&d,&len);
      break;

    case GRIB_TYPE_LONG:
      grib_expression_evaluate_long(a->parent->h,expression,&l);
      pack_long(a,&l,&len);
      break;

    default:
      len = sizeof(tmp);
      p = grib_expression_evaluate_string(a->parent->h,expression,tmp,&len,&ret);
      if (ret != GRIB_SUCCESS) {
        grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,"unable to evaluate %s as string",a->name);
        Assert(0);
      }
      len = strlen(p)+1;
      pack_string(a,p,&len);
      break;
  }
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;
  switch(self->type)
  {
    case GRIB_TYPE_DOUBLE:
      grib_dump_double(dumper,a,NULL);
      break;

    case GRIB_TYPE_LONG:
      grib_dump_long(dumper,a,NULL);
      break;

    default:
      grib_dump_string(dumper,a,NULL);
      break;
  }
}

static int pack_double(grib_accessor* a, const double* val, size_t *len)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;

  if(*len != 1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , 1 );
    *len = 1;
    return GRIB_ARRAY_TOO_SMALL;
  }

  self->dval = *val;
  self->type = ((long)*val == *val) ? GRIB_TYPE_LONG : GRIB_TYPE_DOUBLE;

  return GRIB_SUCCESS;
}

static int pack_long(grib_accessor* a, const long* val, size_t *len)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;

  if(*len != 1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , 1 );
    *len = 1;
    return GRIB_ARRAY_TOO_SMALL;
  }

  self->dval = *val;
  self->type = GRIB_TYPE_LONG;

  return GRIB_SUCCESS;
}

static int unpack_double(grib_accessor* a, double* val, size_t *len)
{
  grib_accessor_variable *ac = (grib_accessor_variable*)a;

  if(*len < 1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , 1 );
    *len = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }
  *val = ac->dval;
  *len = 1;
  return GRIB_SUCCESS;

}
static int unpack_long(grib_accessor* a, long* val, size_t *len)
{
  grib_accessor_variable *ac = (grib_accessor_variable*)a;

  if(*len < 1)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , 1 );
    *len = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }
  *val = (long)ac->dval;
  *len = 1;
  return GRIB_SUCCESS;

}

static int get_native_type(grib_accessor* a)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;
  return self->type;
}

static void destroy(grib_context* c,grib_accessor* a)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;
  grib_context_free(c,self->cval);
}

static int unpack_string (grib_accessor* a, char* val, size_t *len){
  grib_accessor_variable *self = (grib_accessor_variable*)a;

  char buf[80];
  char *p = buf;
  size_t slen ;

  if(self->type == GRIB_TYPE_STRING) {
    p = self->cval;
  } else {
    sprintf(p,"%g",self->dval);
  }

  slen = strlen(p) +1;
  if(*len < slen)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Variable unpack_string Wrong size for %s it is %d bytes big (len=%d)", a->name , slen ,*len);
    *len = slen;
    return GRIB_BUFFER_TOO_SMALL;
  }
  strcpy(val,p);
  *len = slen;

  return GRIB_SUCCESS;
}


static int pack_string(grib_accessor* a, const char* val, size_t *len)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;
  grib_context *c = a->parent->h->context;

  grib_context_free(c,self->cval);
  self->cval = grib_context_strdup(c,val);
  self->dval = atof(self->cval);
  self->type = GRIB_TYPE_STRING;
  return GRIB_SUCCESS;
}

static long value_count(grib_accessor* a)
{
  grib_accessor_variable *self = (grib_accessor_variable*)a;
  if(self->type == GRIB_TYPE_STRING)
    return strlen(self->cval) +1;
  else
    return 1;
}

static int compare(grib_accessor* a, grib_accessor* b) {
  int retval=0;
  double *aval=0;
  double *bval=0;

  size_t alen = (size_t)grib_value_count(a);
  size_t blen = (size_t)grib_value_count(b);

  if (alen != blen) return GRIB_COUNT_MISMATCH;

  aval=grib_context_malloc(a->parent->h->context,alen*sizeof(double));
  bval=grib_context_malloc(b->parent->h->context,blen*sizeof(double));

  grib_unpack_double(a,aval,&alen);
  grib_unpack_double(b,bval,&blen);

  retval = GRIB_SUCCESS;
  while (alen != 0) {
    if (*bval != *aval) retval = GRIB_DOUBLE_VALUE_MISMATCH;
    alen--;
  }

  grib_context_free(a->parent->h->context,aval);
  grib_context_free(b->parent->h->context,bval);

  return retval;

}

