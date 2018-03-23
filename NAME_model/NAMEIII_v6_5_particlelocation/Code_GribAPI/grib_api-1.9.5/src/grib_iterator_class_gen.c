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
   CLASS      = iterator
   IMPLEMENTS = destroy
   IMPLEMENTS = has_next
   IMPLEMENTS = init
   IMPLEMENTS = reset
   MEMBERS     =   long carg
   MEMBERS     =  const char* missingValue;
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "iterator.class" and rerun ./make_class.pl

*/


static void init_class              (grib_iterator_class*);

static int init               (grib_iterator* i,grib_handle*,grib_arguments*);
static int destroy            (grib_iterator* i);
static int reset              (grib_iterator* i);
static long has_next          (grib_iterator* i);


typedef struct grib_iterator_gen{
  grib_iterator it;
/* Members defined in gen */
	long carg;
	const char* missingValue;
} grib_iterator_gen;


static grib_iterator_class _grib_iterator_class_gen = {
    0,                    /* super                     */
    "gen",                    /* name                      */
    sizeof(grib_iterator_gen),/* size of instance          */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                     /* constructor               */
    &destroy,                  /* destructor                */
    0,                     /* Next Value                */
    0,                 /*  Previous Value           */
    &reset,                    /* Reset the counter         */
    &has_next,                 /* has next values           */
};

grib_iterator_class* grib_iterator_class_gen = &_grib_iterator_class_gen;


static void init_class(grib_iterator_class* c)
{
}
/* END_CLASS_IMP */

static int init(grib_iterator* i,grib_handle *h, grib_arguments* args)
{
  grib_iterator_gen* self = (grib_iterator_gen*) i;
  size_t dli=0;
  int ret = GRIB_SUCCESS;
  const char* rawdat  = NULL;
  const char* snumberOfPoints=NULL;
  long numberOfPoints=0;
  self->carg = 1;
  
  snumberOfPoints = grib_arguments_get_name(h,args,self->carg++);
  self->missingValue  = grib_arguments_get_name(h,args,self->carg++);
  rawdat      = grib_arguments_get_name(h,args,self->carg++);

  i->h    = h; /* We may not need to keep them */
  i->args = args;
  if( (ret =  grib_get_size(h,rawdat,&dli))!= GRIB_SUCCESS) return ret;

  if( (ret =  grib_get_long_internal(h,snumberOfPoints,&numberOfPoints))
       != GRIB_SUCCESS)
    return ret;

  if (numberOfPoints!=dli) {
    grib_context_log(h->context,GRIB_LOG_ERROR,"%s != size(%s) (%ld!=%ld)",
                     snumberOfPoints,rawdat,numberOfPoints,dli);
    return GRIB_WRONG_GRID;
  }
  i->nv = dli;
  i->data = grib_context_malloc(h->context,(i->nv)*sizeof(double));

  if( (ret = grib_get_double_array_internal(h,rawdat,i->data ,&(i->nv))))
    return ret;

  i->e = -1;

  return ret;
}

static int reset(grib_iterator* i){
  i->e = -1;
  return 0;
}

static int destroy(grib_iterator* ei){
  const grib_context *c = ei->h->context;
  grib_context_free(c,ei->data);
  return 1;
}


static long has_next(grib_iterator* i){
  if(i->data == NULL) return 0;
  return   i->nv - i->e;
}
