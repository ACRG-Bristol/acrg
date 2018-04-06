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
#include <math.h>

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = iterator
   SUPER      = grib_iterator_class_gen
   IMPLEMENTS = previous;next
   IMPLEMENTS = init;destroy
   MEMBERS     =  double   *las
   MEMBERS     =  double   *los
   MEMBERS     =  long      nap
   MEMBERS     =  long      nam
   MEMBERS     =  long iScansNegatively
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
static int next               (grib_iterator* i, double *lat, double *lon, double *val);
static int previous           (grib_iterator* ei, double *lat, double *lon, double *val);
static int destroy            (grib_iterator* i);


typedef struct grib_iterator_regular{
  grib_iterator it;
/* Members defined in gen */
	long carg;
	const char* missingValue;
/* Members defined in regular */
	double   *las;
	double   *los;
	long      nap;
	long      nam;
	long iScansNegatively;
} grib_iterator_regular;

extern grib_iterator_class* grib_iterator_class_gen;

static grib_iterator_class _grib_iterator_class_regular = {
    &grib_iterator_class_gen,                    /* super                     */
    "regular",                    /* name                      */
    sizeof(grib_iterator_regular),/* size of instance          */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                     /* constructor               */
    &destroy,                  /* destructor                */
    &next,                     /* Next Value                */
    &previous,                 /*  Previous Value           */
    0,                    /* Reset the counter         */
    0,                 /* has next values           */
};

grib_iterator_class* grib_iterator_class_regular = &_grib_iterator_class_regular;


static void init_class(grib_iterator_class* c)
{
	c->reset	=	(*(c->super))->reset;
	c->has_next	=	(*(c->super))->has_next;
}
/* END_CLASS_IMP */


static int next(grib_iterator* i, double *lat, double *lon, double *val){
  grib_iterator_regular* self = (grib_iterator_regular*)i;

  if((long)i->e >= (long)(i->nv-1))  return 0;

  i->e++;

  *lat = self->las[(long)floor(i->e/self->nap)];
  *lon = self->los[(long)i->e%self->nap];
  *val = i->data[i->e];

  return 1;
}


static int previous(grib_iterator* i, double *lat, double *lon, double *val){
  grib_iterator_regular* self = (grib_iterator_regular*)i;

  if(i->e < 0)      return 0;
  *lat = self->las[(long)floor(i->e/self->nap)];
  *lon = self->los[i->e%self->nap];
  *val = i->data[i->e];
  i->e--;

  return 1;
}

static int destroy(grib_iterator* i){
  grib_iterator_regular* self = (grib_iterator_regular*)i;
  const grib_context *c = i->h->context;
  grib_context_free(c,self->las);
  grib_context_free(c,self->los);
  return GRIB_SUCCESS;
}

static int init(grib_iterator* i,grib_handle* h,grib_arguments* args)
{
  grib_iterator_regular* self = (grib_iterator_regular*)i;
  int ret = GRIB_SUCCESS;

  long nap;
  long nam;
  double idir;

  double lof;
  long loi;
  
  const char* longoffirst = grib_arguments_get_name(h,args,self->carg++);
  const char* idirec      = grib_arguments_get_name(h,args,self->carg++);
  const char* nalpar      = grib_arguments_get_name(h,args,self->carg++);
  const char* nalmer      = grib_arguments_get_name(h,args,self->carg++);
  const char* iScansNegatively  = grib_arguments_get_name(h,args,self->carg++);

  
  if((ret = grib_get_double_internal(h,longoffirst,   &lof))) return ret;

  if((ret = grib_get_double_internal(h,idirec,        &idir))) return ret;

  if((ret = grib_get_long_internal(h,nalpar,          &nap))) return ret;
  if((ret = grib_get_long_internal(h,nalmer,          &nam))) return ret;
  if((ret = grib_get_long_internal(h,iScansNegatively,&self->iScansNegatively)))
     return ret;

  if (self->iScansNegatively) {
    idir=-idir;
  } else {
	if (lof+(nap-2)*idir>360) lof-=360;
    else if (lof+nap*idir>360) idir=360.0/(float)nap;
  }

  self->nap = nap;
  self->nam = nam;

  self->las = grib_context_malloc(h->context,nam*sizeof(double));
  self->los = grib_context_malloc(h->context,nap*sizeof(double));

  for( loi = 0; loi < nap; loi++ )  {
    self->los[loi] = lof;
    lof += idir ;
  }

  return ret;
}
