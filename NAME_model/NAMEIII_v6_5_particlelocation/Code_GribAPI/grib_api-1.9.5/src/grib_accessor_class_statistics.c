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
   SUPER      = grib_accessor_class_abstract_vector
   IMPLEMENTS = unpack_double; destroy
   IMPLEMENTS = unpack_string
   IMPLEMENTS = value_count;compare
   IMPLEMENTS = init
   MEMBERS = const char* values
   MEMBERS = const char* missing_value
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int unpack_double(grib_accessor*, double* val,size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long value_count(grib_accessor*);
static void destroy(grib_context*,grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static int compare(grib_accessor*, grib_accessor*);

typedef struct grib_accessor_statistics {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in double */
/* Members defined in abstract_vector */
	double* v;
	int number_of_elements;
/* Members defined in statistics */
	const char* values;
	const char* missing_value;
} grib_accessor_statistics;

extern grib_accessor_class* grib_accessor_class_abstract_vector;

static grib_accessor_class _grib_accessor_class_statistics = {
    &grib_accessor_class_abstract_vector,                      /* super                     */
    "statistics",                      /* name                      */
    sizeof(grib_accessor_statistics),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    &destroy,                    /* free mem                       */
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


grib_accessor_class* grib_accessor_class_statistics = &_grib_accessor_class_statistics;


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

static void init(grib_accessor* a,const long l, grib_arguments* c)
{
  grib_accessor_statistics* self = (grib_accessor_statistics*)a;
  int n = 0;

  self->missing_value = grib_arguments_get_name(a->parent->h,c,n++);
  self->values = grib_arguments_get_name(a->parent->h,c,n++);
  a->flags  |= GRIB_ACCESSOR_FLAG_READ_ONLY;
  a->flags |= GRIB_ACCESSOR_FLAG_FUNCTION;
  a->flags |= GRIB_ACCESSOR_FLAG_HIDDEN;

  self->number_of_elements=8;
  self->v=(double*)grib_context_malloc(a->parent->h->context,
                sizeof(double)*self->number_of_elements);

  a->length=0;
  a->dirty=1;
}

static int    unpack_double   (grib_accessor* a, double* val, size_t *len)
{
  grib_accessor_statistics* self = (grib_accessor_statistics*)a;
  int ret = 0,i=0;
  double* values;
  size_t size=0;
  double max,min,avg,sd,value,skew,kurt,x;
  double missing=0;
  long number_of_missing=0;
  grib_context* c=a->parent->h->context;
  grib_handle* h=a->parent->h;

  if (!a->dirty) return GRIB_SUCCESS;

  if ( (ret=grib_get_size(h,self->values,&size)) != GRIB_SUCCESS) return ret;

  grib_context_log(a->parent->h->context,GRIB_LOG_DEBUG,
    "grib_accessor_statistics: computing statistics for %d values",size);

  if((ret=grib_get_double(a->parent->h,self->missing_value,&missing))
       != GRIB_SUCCESS) return ret;

  values=(double*)grib_context_malloc_clear(c,size*sizeof(double));
  if (!values) return GRIB_OUT_OF_MEMORY;

  if((ret = grib_get_double_array_internal(h,self->values,values,&size))
       != GRIB_SUCCESS) {
        grib_context_free(c,values);
        return ret;
  }

  number_of_missing=0;
  i=0;
  while (i<size && values[i] == missing ) {i++;number_of_missing++;}
  max=values[i];
  min=values[i];
  avg=values[i];
  for (i=number_of_missing+1;i<size;i++) {
    value=values[i];
    if (value > max && value != missing) max=value; 
    if (value < min && value != missing) min=value; 
    if (value != missing) avg+=value;                      
    else number_of_missing++;
  }

  avg/=(size-number_of_missing);

  sd=0; skew=0; kurt=0;
  for (i=0;i<size;i++) {
    if (values[i] != missing) {
      x=(avg-values[i])*(avg-values[i]);
      sd+=x;
      x*=(avg-values[i]);
      skew+=x;
      kurt+=x*(avg-values[i]);
    }
  }

  kurt=kurt;
  if (size-number_of_missing!=0)
	  sd=sqrt(sd/(size-number_of_missing));
  skew=skew;
  
  a->dirty=0;

  grib_context_free(c,values);

  self->v[0]=max;
  self->v[1]=min;
  self->v[2]=avg;
  self->v[3]=number_of_missing;
  self->v[4]=sd;
  self->v[5]=skew;
  self->v[6]=kurt;
  self->v[7]= sd == 0 ? 1 : 0;

  for (i=0;i<self->number_of_elements;i++)
    val[i]=self->v[i];
  
  return ret;
}

static long value_count(grib_accessor* a) {
  grib_accessor_statistics* self = (grib_accessor_statistics*)a;
  return self->number_of_elements;
}

static void destroy(grib_context* c,grib_accessor* a)
{
  grib_accessor_statistics* self = (grib_accessor_statistics*)a;
  grib_context_free(c,self->v);
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

  b->dirty=1;
  a->dirty=1;

  grib_unpack_double(a,aval,&alen);
  grib_unpack_double(b,bval,&blen);

  retval = GRIB_SUCCESS;
  while (alen != 0) {
    if (*bval != *aval) retval = GRIB_DOUBLE_VALUE_MISMATCH;
    alen--;
  }

  grib_context_free(a->parent->h->context,aval);
  grib_context_free(b->parent->h->context,bval);

  return GRIB_SUCCESS;
}

static int unpack_string(grib_accessor*a , char*  v, size_t *len){
	return GRIB_NOT_IMPLEMENTED;
}


